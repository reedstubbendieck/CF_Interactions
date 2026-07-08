import argparse, gzip, os, shutil, subprocess
import pandas as pd

CF = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

# load the metadata file
def load_tsv(file_path):
    return pd.read_csv(file_path, sep='\t')

# builds Bowtie2 index, skipped if already built
def build_bowtie2_index(reference_fasta, index_prefix):
    index_files = [index_prefix + ext for ext in [".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2", ".rev.1.bt2", ".rev.2.bt2"]]
    if not all([os.path.exists(file) for file in index_files]):
        subprocess.run(["bowtie2-build", reference_fasta, index_prefix], check=True)

# download the SRA runs in the metadata file
def download_sra(run_id, output_dir):
    dump_call = f"fasterq-dump -O {output_dir} -e 8 {run_id}"
    subprocess.call(dump_call, shell=True)

# run fastp for quality control of reads
def run_fastp(output_dir, run_id):
    fastq_files = [os.path.join(output_dir, f) for f in os.listdir(output_dir) if f.startswith(run_id) and f.endswith('.fastq')]
    if len(fastq_files) == 1:
        subprocess.run(["fastp", "-i", fastq_files[0], "-o", fastq_files[0].replace(".fastq", "_trimmed.fastq.gz")], check=True)
    elif len(fastq_files) == 2:
        subprocess.run(["fastp", "-i", fastq_files[0], "-I", fastq_files[1],
                        "-o", fastq_files[0].replace(".fastq", "_trimmed.fastq.gz"),
                        "-O", fastq_files[1].replace(".fastq", "_trimmed.fastq.gz")], check=True)

# remove potential contaminating human reads using Bowtie2
def remove_human_reads(output_dir, run_id, human_ref_index):
    fastq_files = [os.path.join(output_dir, f) for f in os.listdir(output_dir)
                   if f.startswith(run_id) and f.endswith('_trimmed.fastq.gz')]

    sam_output = os.path.join(output_dir, run_id + "_aligned_to_human.sam")
    bam_output = sam_output.replace(".sam", ".bam")

    if len(fastq_files) == 1:
        unmapped_fastq_filename = run_id + "_non_human.fastq.gz"
        unmapped_fastq_path = os.path.join(output_dir, unmapped_fastq_filename)

        subprocess.run(["bowtie2", "-x", human_ref_index, "-U", fastq_files[0], "-S", sam_output], check=True)
        subprocess.run(["samtools", "view", "-b", "-f", "4", sam_output, "-o", bam_output], check=True)

        # Safely write uncompressed FASTQ and then compress it
        uncompressed_fastq = unmapped_fastq_path.replace(".gz", "")
        with open(uncompressed_fastq, 'wb') as f_out:
            subprocess.run(["samtools", "fastq", bam_output], stdout=f_out, check=True)

        with open(uncompressed_fastq, 'rb') as f_in, gzip.open(unmapped_fastq_path, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)

        #os.remove(uncompressed_fastq)

    elif len(fastq_files) == 2:
        unmapped_fastq_1 = os.path.join(output_dir, run_id + "_non_human_1.fastq.gz")
        unmapped_fastq_2 = os.path.join(output_dir, run_id + "_non_human_2.fastq.gz")

        subprocess.run(["bowtie2", "-x", human_ref_index, "-1", fastq_files[0], "-2", fastq_files[1], "-S", sam_output], check=True)
        subprocess.run(["samtools", "view", "-b", "-f", "12", sam_output, "-o", bam_output], check=True)
        subprocess.run(["samtools", "fastq", "-f", "12", bam_output, "-1", unmapped_fastq_1, "-2", unmapped_fastq_2], check=True)

    else:
        raise RuntimeError(f"[ERROR] Expected 1 or 2 trimmed fastq files, found {len(fastq_files)} for {run_id}")

# dynamically select and run Mothur with correct batch file based on read type
def run_mothur(run_id, output_dir, mothur_batch_file_pe, mothur_batch_file_se):
    files_in_dir = os.listdir(output_dir)
    se_mode = any(f.startswith(run_id) and f.endswith("_non_human.fastq.gz") for f in files_in_dir)

    mothur_output_dir = os.path.join(os.path.dirname(output_dir.rstrip("/")), "mothur_output")
    mothur_output_reads_dir = os.path.join(mothur_output_dir, "reads")
    mothur_output_taxonomy_dir = os.path.join(mothur_output_dir, "taxonomy")
    if not os.path.exists(mothur_output_reads_dir):
        os.makedirs(mothur_output_reads_dir)
    if not os.path.exists(mothur_output_taxonomy_dir):
        os.makedirs(mothur_output_taxonomy_dir)

    if se_mode:
        print(f"Running Mothur for {run_id} using single-end batch file.")
        updated_batch = os.path.join(output_dir, f"{run_id}_mothur_batch.txt")
        update_mothur_batch_file(mothur_batch_file_se, updated_batch, run_id, output_dir)
        try:
            subprocess.run(["mothur", updated_batch], check=True)
        except subprocess.CalledProcessError:
            pass
    else:
        print(f"Running Mothur for {run_id} using paired-end batch file.")
        try:
            subprocess.run(["mothur", mothur_batch_file_pe], check=True)
        except subprocess.CalledProcessError:
            pass

    if se_mode:
        fasta_file = os.path.join(output_dir, f"{run_id}_non_human.unique.good.filter.precluster.denovo.vsearch.pick.fasta")
        taxonomy_file = os.path.join(output_dir, f"{run_id}_non_human.unique.good.filter.precluster.denovo.vsearch.pick.mothur.wang.pick.taxonomy")
        taxonomy_summary_file = os.path.join(output_dir, f"{run_id}_non_human.unique.good.filter.precluster.denovo.vsearch.pick.mothur.wang.pick.tax.summary")
    else:
        fasta_file = os.path.join(output_dir, "stability.trim.contigs.good.filter.precluster.pick.fasta")
        taxonomy_file = os.path.join(output_dir, "stability.trim.contigs.good.filter.precluster.pick.mothur.wang.pick.taxonomy")
        taxonomy_summary_file = os.path.join(output_dir, "stability.trim.contigs.good.filter.precluster.pick.mothur.wang.tax.summary")

    if os.path.exists(fasta_file):
        compressed_fasta = os.path.join(mothur_output_reads_dir, f"{run_id}_reads.fna.tar.gz")
        subprocess.run(["tar", "-czf", compressed_fasta, "-C", os.path.dirname(fasta_file), os.path.basename(fasta_file)], check=True)

    if os.path.exists(taxonomy_file):
        compressed_taxonomy = os.path.join(mothur_output_taxonomy_dir, f"{run_id}.taxonomy.tar.gz")
        subprocess.run(["tar", "-czf", compressed_taxonomy, "-C", os.path.dirname(taxonomy_file), os.path.basename(taxonomy_file)], check=True)

    new_taxonomy_summary = os.path.join(mothur_output_dir, f"{run_id}.tax.summary")
    if os.path.exists(taxonomy_summary_file):
        os.rename(taxonomy_summary_file, new_taxonomy_summary)


def update_mothur_batch_file(template_path, output_path, run_id, output_dir):
    with open(template_path, 'r') as f:
        lines = f.readlines()

    updated_lines = []
    for line in lines:
        if line.strip().startswith("fastq.info("):
            fastq_path = os.path.join(output_dir, f"{run_id}_non_human.fastq")
            new_line = f"fastq.info(fastq={fastq_path})\n"
            updated_lines.append(new_line)
        else:
            updated_lines.append(line)

    with open(output_path, 'w') as f:
        f.writelines(updated_lines)

# cleanup files after each iteration (part 1)
def cleanup_prefiles(output_dir, run_id):
    files_to_remove = [os.path.join(output_dir, f) for f in os.listdir(output_dir)
                       if f.startswith(run_id) and (
                           #f.endswith(".fastq") or
                           f.endswith("_trimmed.fastq.gz") or
                           f.endswith("_aligned_to_human.sam") or
                           f.endswith("_aligned_to_human.bam")
                       )]
    for file in files_to_remove:
        os.remove(file)

# cleanup files after each iteration (part 2)
def cleanup_postfiles(output_dir, run_id):
    files_to_remove = [os.path.join(output_dir, f) for f in os.listdir(output_dir)
                       if f.startswith(run_id)]
    for file in files_to_remove:
        os.remove(file)

    for file in os.listdir("."):
        if file in ["fastp.json", "fastp.html"] or file.startswith("mothur."):
            os.remove(file)

# main function
def main(tsv_file, output_dir, reference_fasta, human_ref_index, mothur_batch_file_pe, mothur_batch_file_se, local_reads=None):
    os.makedirs(output_dir, exist_ok=True)
    df = load_tsv(tsv_file)

    # mothur output lives beside the reads output directory
    mothur_output_dir = os.path.join(os.path.dirname(output_dir.rstrip("/\\")), "mothur_output")

    build_bowtie2_index(reference_fasta, human_ref_index)

    for run_id in df['Run']:
        if os.path.exists(output_dir):
            for filename in os.listdir(output_dir):
                file_path = os.path.join(output_dir, filename)
                try:
                    if os.path.isfile(file_path):
                        os.unlink(file_path)
                except Exception as e:
                    print(f"Failed to delete {file_path}. Reason: {e}")

        try:
            if os.path.exists(os.path.join(mothur_output_dir, f"{run_id}.tax.summary")):
                print(f"Run {run_id} already processed, skipping.")
                continue

            # use a local copy of the reads if one is provided, otherwise download from the SRA
            def local(name):
                return os.path.join(local_reads, name) if local_reads else None
            if local_reads and os.path.exists(local(f"{run_id}_1.fastq")) and os.path.exists(local(f"{run_id}_2.fastq")):
                shutil.copyfile(local(f"{run_id}_1.fastq"), os.path.join(output_dir, f"{run_id}_1.fastq"))
                shutil.copyfile(local(f"{run_id}_2.fastq"), os.path.join(output_dir, f"{run_id}_2.fastq"))
            elif local_reads and os.path.exists(local(f"{run_id}_R1.fastq.gz")) and os.path.exists(local(f"{run_id}_R2.fastq.gz")):
                for suffix in ["R1", "R2"]:
                    gz_path = local(f"{run_id}_{suffix}.fastq.gz")
                    out_path = os.path.join(output_dir, f"{run_id}_{suffix}.fastq.gz")
                    shutil.copyfile(gz_path, out_path)
                    with gzip.open(out_path, 'rb') as f_in, open(out_path.replace(".gz", ""), 'wb') as f_out:
                        shutil.copyfileobj(f_in, f_out)
                    os.remove(out_path)
            elif local_reads and os.path.exists(local(f"{run_id}.fastq")):
                shutil.copyfile(local(f"{run_id}.fastq"), os.path.join(output_dir, f"{run_id}.fastq"))
            elif local_reads and os.path.exists(local(f"{run_id}.fastq.gz")):
                gz_path = local(f"{run_id}.fastq.gz")
                out_path = os.path.join(output_dir, f"{run_id}.fastq.gz")
                shutil.copyfile(gz_path, out_path)
                with gzip.open(out_path, 'rb') as f_in, open(out_path.replace(".gz", ""), 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)
                os.remove(out_path)
            else:
                download_sra(run_id, output_dir)

            if not os.listdir(output_dir):
                print(f"Download failed for {run_id}. Output directory is empty. Skipping to next run.")
                continue

            run_fastp(output_dir, run_id)
            remove_human_reads(output_dir, run_id, human_ref_index)
            cleanup_prefiles(output_dir, run_id)
            run_mothur(run_id, output_dir, mothur_batch_file_pe, mothur_batch_file_se)
            cleanup_postfiles(output_dir, run_id)

        except Exception as e:
            print(f"Error processing run {run_id}: {e}. Skipping to next run.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process sequencing data pipeline.")
    parser.add_argument("-i", "--input", required=True, help="Path to the TSV file with metadata.")
    parser.add_argument("-o", "--output", required=True, help="Output directory for processed files.")
    parser.add_argument("--mothur_pe", required=True, help="Path to the paired-end Mothur batch file.")
    parser.add_argument("--mothur_se", required=True, help="Path to the single-end Mothur batch file.")
    parser.add_argument("--local-reads", default=os.environ.get("LOCAL_READS"),
                        help="Optional directory of pre-downloaded FASTQ files to use instead of downloading from the SRA "
                             "(falls back to the LOCAL_READS environment variable).")
    args = parser.parse_args()

    tsv_file = args.input
    output_dir = args.output
    mothur_batch_file_pe = args.mothur_pe
    mothur_batch_file_se = args.mothur_se
    reference_fasta = os.path.join(CF, "rawData", "references", "GCA_009914755.4_T2T-CHM13v2.0_genomic.fna")
    human_ref_index = os.path.join(CF, "rawData", "references", "human_genome_index")

    main(tsv_file, output_dir, reference_fasta, human_ref_index, mothur_batch_file_pe, mothur_batch_file_se, args.local_reads)