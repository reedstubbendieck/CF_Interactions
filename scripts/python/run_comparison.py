#!/usr/bin/env python3
# usage:
#   python3 run_comparison.py               # skips finished heavy steps
#   python3 run_comparison.py --force        # re-run every stage
#   python3 run_comparison.py --skip-tools   # only rebuild tables + plotdata

import argparse, contextlib, csv, glob, io, json, math, os, subprocess, sys
from collections import defaultdict, deque

@contextlib.contextmanager
def safe_open(path):
    buf = io.StringIO()
    yield buf
    try:
        with open(path, "w", newline="") as fh:
            fh.write(buf.getvalue())
    except PermissionError:
        print(f"    WARN: {os.path.basename(path)} is locked; kept existing copy")

CF = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
GC = os.path.join(CF, "derivedData", "genome_comparison")   # tool outputs, tables, plotdata
G = os.path.join(CF, "rawData", "genome_sequences")
REF = os.path.join(G, "ATCC10145.fna")
BAKTA_DB = os.environ.get("BAKTA_DB", "")          # required only for the Bakta stage
THREADS = int(os.environ.get("THREADS", 0)) or (os.cpu_count() or 4)
# conda env names (override per system); tool stages are skipped if outputs exist
ENV_BAKTA = os.environ.get("BAKTA_ENV", "bakta")
ENV_SNIPPY = os.environ.get("SNIPPY_ENV", "snippy")
ENV_ANTISMASH = os.environ.get("ANTISMASH_ENV", "antismash")
ENV_BLAST = os.environ.get("BLAST_ENV", "ncbi_tools")

# isolate -> (donor, phenotype); reference isolate per donor is the first inhibitory one
ISOLATES = {"NB0017": ("HID07", "inhibitory"), "NB0078": ("HID07", "inhibitory"),
            "NB1046": ("HID07", "non-inhibitory"), "NB0219": ("HID08", "inhibitory"),
            "RSD0018": ("HID08", "inhibitory"), "NB1386": ("HID08", "non-inhibitory")}
DONORS = {  # figure / plotdata (reference = non-inhibitory, ref first)
    "HID07": dict(panel="A", ref="NB1046", isolates=["NB1046", "NB0017", "NB0078"], noninh="NB1046"),
    "HID08": dict(panel="B", ref="NB1386", isolates=["NB1386", "NB0219", "RSD0018"], noninh="NB1386"),
}
DONORS_TAB = {  # tables (reference = inhibitory, ref first)
    "HID07": dict(panel="A", ref="NB0017", isolates=["NB0017", "NB0078", "NB1046"], noninh="NB1046"),
    "HID08": dict(panel="B", ref="NB0219", isolates=["NB0219", "RSD0018", "NB1386"], noninh="NB1386"),
}

def distinguishes_noninh(v, cfg):
    ref, noninh = cfg["ref"], cfg["noninh"]
    carriers = v["carriers"]
    if noninh == ref:
        return set(i for i in cfg["isolates"] if i != ref) <= carriers
    other_inh = set(i for i in cfg["isolates"] if i not in (ref, noninh))
    return (noninh in carriers) and not (other_inh & carriers)

# antiSMASH standard colours
CAT_COLOR = {"NRPS": "#2E8B57", "PKS": "#F4A460", "RiPP": "#4169E1",
             "terpene": "#800080", "saccharide": "#F2F2F2", "other": "#191970"}
PROD_COLOR = {"NRP-metallophore": "#DC143C", "opine-like-metallophore": "#DC143C",
              "phenazine": "#DDA0DD", "hserlactone": "#D2B48C", "NAGGN": "#9ACD32",
              "betalactone": "#800080", "terpene-precursor": "#800080"}
PROD_CAT = {"CDPS": "NRPS", "NAGGN": "other", "NRP-metallophore": "NRPS", "NRPS": "NRPS",
            "NRPS-like": "NRPS", "RiPP-like": "RiPP", "betalactone": "other",
            "hserlactone": "other", "hydrogen-cyanide": "other",
            "opine-like-metallophore": "other", "phenazine": "other",
            "redox-cofactor": "RiPP", "terpene-precursor": "terpene"}

def sh(cmd, log=None):
    print("  $", cmd if isinstance(cmd, str) else " ".join(cmd))
    out = open(log, "w") if log else None
    subprocess.run(cmd, shell=isinstance(cmd, str), check=True, stdout=out,
                   stderr=subprocess.STDOUT if out else None)
    if out:
        out.close()

def conda(env, args, log=None):
    sh(["conda", "run", "-n", env] + args, log=log)

# stage 1: fastANI
def run_fastani(force):
    out = f"{GC}/fastani"; os.makedirs(out, exist_ok=True)
    if not force and os.path.exists(f"{out}/all_vs_all.tsv"):
        print("[1] fastANI: cached"); return
    print("[1] fastANI")
    isos = [i for i in ISOLATES]
    with open(f"{out}/queries.txt", "w") as fh:
        fh.write("\n".join(f"{G}/{i}.fna" for i in isos) + "\n")
    with open(f"{out}/all_genomes.txt", "w") as fh:
        fh.write("\n".join([REF] + [f"{G}/{i}.fna" for i in isos]) + "\n")
    sh(f"fastANI --ql {out}/queries.txt --ref {REF} -o {out}/isolates_vs_ATCC10145.tsv", f"{GC}/logs/fastani_ref.log")
    sh(f"fastANI --ql {out}/all_genomes.txt --rl {out}/all_genomes.txt -o {out}/all_vs_all.tsv", f"{GC}/logs/fastani_ava.log")

# stage 2: Bakta
def run_bakta(force):
    print("[2] Bakta")
    for iso in ISOLATES:
        od = f"{GC}/bakta/{iso}"
        if not force and os.path.exists(f"{od}/{iso}.gbff"):
            print(f"    {iso}: cached"); continue
        if not BAKTA_DB:
            sys.exit("ERROR: set the BAKTA_DB environment variable to the Bakta database path")
        conda(ENV_BAKTA, ["bakta", "--db", BAKTA_DB, "--genus", "Pseudomonas",
                        "--species", "aeruginosa", "--strain", iso, "--complete",
                        "--prefix", iso, "--output", od, "--threads", str(THREADS),
                        "--force", f"{G}/{iso}.fna"], log=f"{GC}/logs/bakta_{iso}.log")

# stage 3: snippy
def run_snippy(force):
    print("[3] snippy")
    for cfg, base in ((DONORS, "snippy"), (DONORS_TAB, "snippy_tab")):
        for donor, d in cfg.items():
            refgbff = f"{GC}/bakta/{d['ref']}/{d['ref']}.gbff"
            ddir = f"{GC}/{base}/{donor}"; os.makedirs(ddir, exist_ok=True)
            for q in [i for i in d["isolates"] if i != d["ref"]]:
                if not force and os.path.exists(f"{ddir}/{q}/snps.tab"):
                    print(f"    {base}/{donor}/{q}: cached"); continue
                conda(ENV_SNIPPY, ["snippy", "--outdir", f"{ddir}/{q}", "--ref", refgbff,
                                 "--ctgs", f"{G}/{q}.fna", "--cpus", str(THREADS), "--force"],
                      log=f"{GC}/logs/snippy_{base}_{donor}_{q}.log")

# stage 4: antiSMASH
def run_antismash(force):
    print("[4] antiSMASH 8")
    for iso in ISOLATES:
        od = f"{GC}/antismash/{iso}"
        if not force and os.path.exists(f"{od}/{iso}.json"):
            print(f"    {iso}: cached"); continue
        conda(ENV_ANTISMASH, ["antismash", "--taxon", "bacteria", "--genefinding-tool", "none",
                            "--cb-knownclusters", "--cc-mibig", "--asf", "--pfam2go",
                            "--cpus", str(THREADS), "--output-dir", od, "--output-basename", iso,
                            f"{GC}/bakta/{iso}/{iso}.gbff"], log=f"{GC}/logs/antismash_{iso}.log")

# stage 5: origins
def genome_len(iso):
    n = 0
    for line in open(f"{G}/{iso}.fna"):
        if not line.startswith(">"):
            n += len(line.strip())
    return n

# position in each donor reference that aligns to reference (ATCC) nt 1
def run_origins(force):
    out = f"{GC}/metadata/origin_offsets.tsv"
    if not force and os.path.exists(out):
        print("[5] origins: cached")
    else:
        print("[5] origins (blastn ref head vs donor references)")
        head = f"{GC}/logs/_refhead.fa"
        seq = "".join(l.strip() for l in open(REF) if not l.startswith(">"))
        open(head, "w").write(">ATCC_head\n" + seq[:5000] + "\n")
        rows = []
        for donor, d in DONORS.items():
            iso = d["ref"]
            res = subprocess.run(["conda", "run", "-n", ENV_BLAST, "blastn", "-query", head,
                                  "-subject", f"{G}/{iso}.fna", "-outfmt", "6 qstart qend sstart send sstrand"],
                                 capture_output=True, text=True, check=True).stdout.splitlines()
            hits = [r.split("\t") for r in res if r.strip() and r.split("\t")[0].strip().isdigit()]
            best = min(hits, key=lambda x: int(x[0]))       # HSP covering query pos 1
            rows.append([iso, genome_len(iso), best[2], best[4]])
        with open(out, "w", newline="") as fh:
            w = csv.writer(fh, delimiter="\t")
            w.writerow(["isolate", "genome_length", "ref_nt1_maps_to_pos", "strand"])
            w.writerows(rows)
    origins = {}
    for r in csv.DictReader(open(out), delimiter="\t"):
        origins[r["isolate"]] = (int(r["genome_length"]), int(r["ref_nt1_maps_to_pos"]))
    return origins

# stage 6: tables
ANN = ["EFFECT", "LOCUS_TAG", "GENE", "PRODUCT"]

def load_snps(path):
    return list(csv.DictReader(open(path), delimiter="\t")) if os.path.isfile(path) else []

def vclass(effect):
    e = effect or ""
    if not e:
        return "Intergenic"
    if "frameshift" in e:
        return "Frameshift"
    if "synonymous" in e:
        return "Silent"
    if any(k in e for k in ("missense", "stop_gained", "stop_lost",
                            "initiator_codon", "start_lost", "inframe")):
        return "Missense"
    return "Non-Coding"

# merge per-query snps.tab into {donor: [variant dicts]} for a config
def collect_variants(cfg, snippy_dir):
    out = {}
    for donor, d in cfg.items():
        ref = d["ref"]; variants = {}
        for iso in d["isolates"]:
            if iso == ref:
                continue
            for r in load_snps(f"{GC}/{snippy_dir}/{donor}/{iso}/snps.tab"):
                key = (r["CHROM"], int(r["POS"]), r["TYPE"], r["REF"], r["ALT"])
                v = variants.setdefault(key, dict(CHROM=r["CHROM"], POS=int(r["POS"]), TYPE=r["TYPE"],
                                                  REF=r["REF"], ALT=r["ALT"], carriers=set(),
                                                  **{c: r.get(c, "") for c in ANN}))
                v["carriers"].add(iso)
                for c in ANN:
                    if not v[c] and r.get(c):
                        v[c] = r[c]
        out[donor] = sorted(variants.values(), key=lambda v: (v["CHROM"], v["POS"]))
    return out


def named_regions(ref, canon):
    pools = defaultdict(deque)
    for c in sorted(canon, key=lambda c: c["num"]):
        pools[c["sig"]].append(c)
    regs = []
    for s, e, prods in sorted(load_areas(ref), key=lambda x: x[0]):
        pool = pools[frozenset(prods)]
        c = pool.popleft() if pool else None
        regs.append(dict(start=s, end=e, num=(c["num"] if c else ""),
                         label=(c["label"] if c else fallback_label(prods))))
    return regs

def build_variant_tables(cfg, snippy_dir, canon):
    print("[6] variant tables (reference = inhibitory isolate)")
    tab = f"{GC}/tables"; os.makedirs(tab, exist_ok=True)
    summary = []
    dv = collect_variants(cfg, snippy_dir)
    for donor, d in cfg.items():
        ref = d["ref"]; rows = dv[donor]
        cols = [f"{i}[{ISOLATES[i][1][:3]}]" for i in d["isolates"]]

        # differences.tsv (genotype matrix)
        with safe_open(f"{tab}/{donor}_differences.tsv") as fh:
            w = csv.writer(fh, delimiter="\t")
            w.writerow(["CHROM", "POS", "TYPE", "REF", "ALT"] + cols + ANN)
            for v in rows:
                geno = [v["REF"] if i == ref else (v["ALT"] if i in v["carriers"] else v["REF"])
                        for i in d["isolates"]]
                w.writerow([v["CHROM"], v["POS"], v["TYPE"], v["REF"], v["ALT"]] + geno + [v[c] for c in ANN])

        # mutations that distinguish the non-inhibitory isolate from its inhibitory relatives
        spec = [v for v in rows if distinguishes_noninh(v, d)]
        with safe_open(f"{tab}/{donor}_nonInhibitory_specific.tsv") as fh:
            w = csv.writer(fh, delimiter="\t")
            w.writerow(["CHROM", "POS", "TYPE", "REF", "ALT", "EFFECT", "GENE", "PRODUCT", "LOCUS_TAG"])
            for v in spec:
                w.writerow([v["CHROM"], v["POS"], v["TYPE"], v["REF"], v["ALT"],
                            v["EFFECT"], v["GENE"], v["PRODUCT"], v["LOCUS_TAG"]])

        # genes affected
        genes = {}
        for v in rows:
            k = (v["GENE"] or "-", v["PRODUCT"] or "intergenic/unannotated", v["LOCUS_TAG"] or "-")
            g = genes.setdefault(k, dict(n=0, types=[], eff=set(), pos=[]))
            g["n"] += 1; g["types"].append(v["TYPE"]); g["pos"].append(str(v["POS"]))
            g["eff"].add((v["EFFECT"] or "-").split(" ")[0])
        with safe_open(f"{tab}/{donor}_genes_affected.tsv") as fh:
            w = csv.writer(fh, delimiter="\t")
            w.writerow(["GENE", "PRODUCT", "LOCUS_TAG", "n_variants", "types", "effects", "positions"])
            for (gene, prod, lt), g in genes.items():
                w.writerow([gene, prod, lt, g["n"], ",".join(g["types"]),
                            ",".join(sorted(g["eff"])), ",".join(g["pos"])])

        # pairwise variant-distance matrix among the donor's isolates
        def geno(v, iso):
            return v["REF"] if iso == ref else (v["ALT"] if iso in v["carriers"] else v["REF"])
        with safe_open(f"{tab}/{donor}_pairwise_variant_distance.tsv") as fh:
            w = csv.writer(fh, delimiter="\t")
            w.writerow([""] + d["isolates"])
            for a in d["isolates"]:
                w.writerow([a] + [sum(1 for v in rows if geno(v, a) != geno(v, b))
                                  for b in d["isolates"]])

        # variants inside a BGC, named from this (inhibitory) reference's antiSMASH
        regs = named_regions(ref, canon)
        def bgc_at(pos):
            for r in regs:
                if r["start"] <= pos <= r["end"]:
                    return r["num"], r["label"]
            return "", ""
        in_bgc = [(v, *bgc_at(v["POS"]), distinguishes_noninh(v, d)) for v in rows]
        in_bgc = [x for x in in_bgc if x[2]]     # keep only those inside a BGC
        with safe_open(f"{tab}/{donor}_variants_in_BGCs.tsv") as fh:
            w = csv.writer(fh, delimiter="\t")
            w.writerow(["donor", "isolate", "gene", "position", "type", "class",
                        "effect", "BGC_number", "BGC", "non_inhibitory_specific"])
            for v, bn, bl, s in sorted(in_bgc, key=lambda x: not x[3]):
                w.writerow([donor, d["noninh"] if s else "inhibitory/shared", v["GENE"] or "(intergenic)",
                            v["POS"], v["TYPE"], vclass(v["EFFECT"]), v["EFFECT"] or "",
                            bn, bl, "yes" if s else "no"])
        print(f"    {donor}: {len(spec)} variants distinguish {d['noninh']}, "
              f"{sum(1 for x in in_bgc if x[3])} of them inside a BGC")

        n_snp = sum(1 for v in rows if v["TYPE"] == "snp")
        n_indel = sum(1 for v in rows if v["TYPE"] in ("ins", "del"))
        summary.append([donor, ref, len(rows), n_snp, n_indel, len(rows) - n_snp - n_indel])

    with safe_open(f"{tab}/summary_counts.tsv") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["donor", "reference_isolate", "total_variants", "snps", "indels", "other"])
        w.writerows(summary)

def region_color(products):
    for p in products:
        if p in PROD_COLOR:
            return PROD_COLOR[p]
    return CAT_COLOR.get(PROD_CAT.get(products[0], "other"), "#191970")

BGC_NAMES = {1: "Pyochelin", 6: "Azetidomonamide", 8: "Pyoverdine",
             9: "AMB", 10: "HCN", 13: "Pyocyanin", 16: "Pseudopaline"}


ACRONYMS = {"nrps": "NRPS", "ripp": "RiPP", "hcn": "HCN", "hsl": "HSL",
            "naggn": "NAGGN", "cdps": "CDPS", "nrp": "NRP", "pks": "PKS",
            "amb": "AMB", "hserlactone": "HSL", "betalactone": "β-lactone"}


def _tc(tok):
    return ACRONYMS.get(tok.lower(), tok.capitalize())


def fallback_label(products):
    return "-".join("-".join(_tc(t) for t in p.split("-")) for p in products)


def load_areas(iso):
    d = json.load(open(f"{GC}/antismash/{iso}/{iso}.json"))
    return [(a["start"] + 1, a["end"], a["products"]) for a in d["records"][0]["areas"]]


def frac_fn(L, origin):
    return lambda pos: ((pos - origin) % L) / L


def build_canonical(origins):
    ref = DONORS["HID07"]["ref"]; L, origin = origins[ref]; f = frac_fn(L, origin)
    regs = []
    for s, e, prods in load_areas(ref):
        fs, fe = f(s), f(e)
        regs.append(dict(f_mid=(fs + ((fe - fs) % 1.0) / 2) % 1.0,
                         sig=frozenset(prods), prods=prods, color=region_color(prods)))
    regs.sort(key=lambda r: r["f_mid"])
    canon = []
    for i, r in enumerate(regs, 1):
        canon.append(dict(num=i, sig=r["sig"], color=r["color"], f_mid=r["f_mid"],
                          label=BGC_NAMES.get(i, fallback_label(r["prods"]))))
    return canon

def align_transform(canon, regs):
    from collections import Counter
    csig = Counter(c["sig"] for c in canon)
    uniq_f = {c["sig"]: c["f_mid"] for c in canon if csig[c["sig"]] == 1}
    pairs = [(uniq_f[r["sig"]], r["f_mid"]) for r in regs if r["sig"] in uniq_f]
    best = (1, 0.0, 1e9)
    for o in (1, -1):
        sx = sum(math.sin(2 * math.pi * (fa - o * fd)) for fa, fd in pairs)
        cx = sum(math.cos(2 * math.pi * (fa - o * fd)) for fa, fd in pairs)
        off = (math.atan2(sx, cx) / (2 * math.pi)) % 1.0
        score = sum(1 - math.cos(2 * math.pi * (fa - (o * fd + off))) for fa, fd in pairs)
        if score < best[2]:
            best = (o, off, score)
    return best[0], best[1]

def build_bgc_table():
    print("[6] BGC table (all isolates)")
    with safe_open(f"{GC}/tables/bgc_regions_all_isolates.tsv") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["isolate", "region", "start", "end", "length_bp", "products"])
        for iso in ISOLATES:
            for i, (s, e, prods) in enumerate(load_areas(iso), 1):
                w.writerow([iso, i, s, e, e - s, ",".join(prods)])

# stage 7: plotdata
def export_plotdata(origins, cfg, snippy_dir, canon):
    print("[7] plotdata for figure (reference = non-inhibitory isolate)")
    pd = f"{GC}/plotdata"; os.makedirs(pd, exist_ok=True)
    dv = collect_variants(cfg, snippy_dir)

    # single shared BGC legend
    with safe_open(f"{pd}/bgc_legend.csv") as fh:
        w = csv.writer(fh)
        w.writerow(["num", "label", "color"])
        for c in canon:
            w.writerow([c["num"], c["label"], c["color"]])

    meta = []
    for donor, d in cfg.items():
        ref = d["ref"]; L, origin = origins[ref]; frac = frac_fn(L, origin)

        areas = load_areas(ref)
        regs = []
        for s, e, prods in areas:
            fs, fe = frac(s), frac(e)
            regs.append(dict(start=s, end=e, f_start=fs, f_end=fe,
                             f_mid=(fs + ((fe - fs) % 1.0) / 2) % 1.0,
                             sig=frozenset(prods), products=",".join(prods)))
        orient, offset = align_transform(canon, regs)   # reorient panel B onto A

        canon_by_sig = defaultdict(list)
        for c in canon:
            canon_by_sig[c["sig"]].append(c)
        regs_by_sig = defaultdict(list)
        for r in sorted(regs, key=lambda r: ((orient * r["f_mid"]) + offset) % 1.0):
            regs_by_sig[r["sig"]].append(r)
        for sig, rlist in regs_by_sig.items():
            for r, c in zip(rlist, canon_by_sig[sig]):
                r["num"], r["label"], r["color"] = c["num"], c["label"], c["color"]
        with safe_open(f"{pd}/bgc_{donor}.csv") as fh:
            w = csv.writer(fh)
            w.writerow(["num", "start", "end", "f_start", "f_end", "f_mid", "label", "color", "products"])
            for r in sorted(regs, key=lambda r: r["num"]):
                w.writerow([r["num"], r["start"], r["end"], f"{r['f_start']:.6f}",
                            f"{r['f_end']:.6f}", f"{r['f_mid']:.6f}", r["label"], r["color"], r["products"]])

        # variant ticks
        def bgc_at(pos):
            for r in regs:
                if r["start"] <= pos <= r["end"]:
                    return r["num"], r["label"]
            return "", ""
        with safe_open(f"{pd}/mut_{donor}.csv") as fh:
            w = csv.writer(fh)
            w.writerow(["pos", "f", "class", "inbgc", "noninh_specific", "gene", "type", "effect"])
            for v in dv[donor]:
                bn, bl = bgc_at(v["POS"])
                w.writerow([v["POS"], f"{frac(v['POS']):.6f}", vclass(v["EFFECT"]),
                            int(bool(bl)), int(distinguishes_noninh(v, d)),
                            v["GENE"] or "", v["TYPE"], v["EFFECT"] or ""])
        meta.append([donor, d["panel"], ref, L, origin, len(regs),
                     len(dv[donor]), orient, f"{offset:.6f}"])

    with safe_open(f"{pd}/meta.csv") as fh:
        w = csv.writer(fh)
        w.writerow(["donor", "panel", "reference", "length", "origin", "n_bgc",
                    "n_variants", "orient", "offset"])
        w.writerows(meta)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--force", action="store_true", help="re-run finished heavy steps")
    ap.add_argument("--skip-tools", action="store_true", help="only rebuild tables + plotdata")
    a = ap.parse_args()
    os.makedirs(f"{GC}/logs", exist_ok=True)
    if not a.skip_tools:
        run_fastani(a.force)
        run_bakta(a.force)
        run_snippy(a.force)
        run_antismash(a.force)
    origins = run_origins(a.force)
    canon = build_canonical(origins)
    build_variant_tables(DONORS_TAB, "snippy_tab", canon)
    build_bgc_table()
    export_plotdata(origins, DONORS, "snippy", canon)

if __name__ == "__main__":
    main()
