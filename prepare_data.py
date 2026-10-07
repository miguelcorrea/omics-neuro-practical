# /// script
# requires-python = ">=3.10"
# dependencies = [
#     "numpy",
#     "pandas",
#     "pydeseq2==0.5.4",
# ]
# ///
"""Instructor-only: build the small, clean data files the student notebook reads.

Run once with:  uv run prepare_data.py

Downloads (cached in data_raw/):
  * GEO GSE113754 count matrix + series matrix (Ingiosi et al. 2019, eLife 8:e42819)
  * NCBI mouse gene_info (Entrez ID -> gene symbol)
  * g:Profiler mouse GMT (GO gene sets, annotations propagated up the GO graph)
  * go-basic.obo (to label each GO term as BP / MF / CC)

Writes (to public/, which is bundled with the notebook):
  samples.csv, genes.csv, counts.csv.gz, de_results.csv.gz, go_terms.csv.gz,
  marker_genes.csv, data_versions.txt
"""

import datetime as dt
import gzip
import json
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

HERE = Path(__file__).parent
RAW = HERE / "data_raw"
OUT = HERE / "public"
RAW.mkdir(exist_ok=True)
OUT.mkdir(exist_ok=True)

URLS = {
    "counts.txt.gz": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE113nnn/GSE113754/suppl/GSE113754_RNASeq_S3_WT_SD5_HC5_Counts.txt.gz",
    "series_matrix.txt.gz": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE113nnn/GSE113754/matrix/GSE113754_series_matrix.txt.gz",
    "Mus_musculus.gene_info.gz": "https://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/Mus_musculus.gene_info.gz",
    "gprofiler_mmusculus.gmt": "https://biit.cs.ut.ee/gprofiler/static/gprofiler_full_mmusculus.name.gmt",
    "go-basic.obo": "https://current.geneontology.org/ontology/go-basic.obo",
}

# Expression filter: keep genes with CPM > 1 in at least MIN_SAMPLES samples
# (5 = the size of one experimental group).
MIN_CPM, MIN_SAMPLES = 1.0, 5

# GO term size limits (counted within all genes in the count table). The
# notebook lets students narrow this further.
GO_MIN, GO_MAX = 5, 2000

# Curated genes for dot plots, with a one-line note shown to students.
MARKERS = {
    "Arc": "Activity-regulated gene; rises with wakefulness",
    "Fos": "Classic immediate-early gene (neuronal activity)",
    "Egr1": "Immediate-early transcription factor",
    "Nr4a1": "Immediate-early nuclear receptor",
    "Homer1": "Synaptic scaffold; short isoform Homer1a is a sleep-loss marker",
    "Sgk1": "Glucocorticoid-induced kinase (stress hormone response)",
    "Fkbp5": "Glucocorticoid-responsive gene",
    "Hspa5": "ER chaperone (BiP); unfolded-protein response",
    "Per1": "Clock gene",
    "Per2": "Clock gene",
    "Dbp": "Clock-controlled transcription factor",
    "Cirbp": "Cold-inducible RNA-binding protein; falls with sleep loss",
    "Rbm3": "Cold-inducible RNA-binding protein",
    "Shank3": "The mutated gene (exon 21 deleted in the mutants)",
    "Gapdh": "Housekeeping gene (should not change)",
}


def fetch(name: str) -> Path:
    path = RAW / name
    if not path.exists():
        print(f"downloading {name} ...")
        # Some servers (geneontology.org) reject urllib's default User-Agent.
        req = urllib.request.Request(URLS[name], headers={"User-Agent": "omics-teaching/1.0"})
        with urllib.request.urlopen(req, timeout=120) as r:
            path.write_bytes(r.read())
    return path


def read_sample_metadata(path: Path) -> pd.DataFrame:
    """Parse the GEO series matrix into one row per sample."""
    fields: dict[str, list[str]] = {}
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if not line.startswith("!Sample_"):
                continue
            key, *vals = line.rstrip("\n").split("\t")
            vals = [v.strip('"') for v in vals]
            if key == "!Sample_characteristics_ch1":
                k = vals[0].split(":", 1)[0].strip().lower()
                fields[k] = [v.split(":", 1)[1].strip() for v in vals]
            elif key in ("!Sample_geo_accession", "!Sample_title", "!Sample_description"):
                fields[key.removeprefix("!Sample_")] = vals
    return pd.DataFrame(fields)


def build_samples(counts: pd.DataFrame, series_matrix: Path) -> pd.DataFrame:
    meta = read_sample_metadata(series_matrix)
    # The count columns are named like the GEO 'description' field (e.g. WTHC5_1);
    # match by name, never by position (the orders differ).
    meta = meta.set_index("description").loc[counts.columns]
    m = meta.index.str.extract(r"^(WT|S3)(HC|SD)5_(\d)$")
    samples = pd.DataFrame(
        {
            "sample": meta.index,
            "gsm": meta["geo_accession"].values,
            "genotype": m[0].map({"WT": "WT", "S3": "Shank3"}).values,
            "condition": m[1].map({"HC": "Sleep", "SD": "Sleep-deprived"}).values,
            # GEO: "The protocol was repeated over the course of 5 days ... so that
            # 5 independent SD experiments were represented". We ASSUME replicate k
            # = experiment day k (supported by PCA; not stated explicitly by GEO).
            "day": ("Day " + m[2]).values,
            "geo_genotype": meta["genotype"].values,
            "geo_title": meta["title"].values,
        }
    )
    # Cross-check our parsed labels against GEO's own annotations.
    geo_geno = samples["geo_genotype"].str.contains("Shank3").map({True: "Shank3", False: "WT"})
    assert (geo_geno == samples["genotype"]).all(), "genotype mismatch vs GEO"
    geo_sd = samples["geo_title"].str.contains("sleep deprivation")
    assert ((samples["condition"] == "Sleep-deprived") == geo_sd).all(), "condition mismatch vs GEO"
    samples["group"] = samples["genotype"] + " | " + samples["condition"]
    samples["library_size"] = counts.sum().values
    return samples.drop(columns=["geo_genotype", "geo_title"])


def build_genes(counts: pd.DataFrame) -> pd.DataFrame:
    info = pd.read_csv(
        fetch("Mus_musculus.gene_info.gz"),
        sep="\t",
        usecols=["GeneID", "Symbol", "description", "type_of_gene"],
        index_col="GeneID",
    )
    genes = info.reindex(counts.index)
    missing = genes["Symbol"].isna()
    print(f"{missing.sum()} Entrez IDs no longer in gene_info; labelled by ID")
    genes.loc[missing, "Symbol"] = "Entrez:" + genes.index[missing].astype(str)
    genes = genes.rename(columns={"Symbol": "symbol", "type_of_gene": "gene_type"})
    genes.index.name = "entrez"
    assert genes["symbol"].is_unique, "duplicate symbols"
    return genes


def run_deseq2(counts: pd.DataFrame, samples: pd.DataFrame) -> tuple[pd.DataFrame, pd.Series]:
    """Two models x four contrasts. Returns a tidy table and DESeq2 baseMean."""
    meta = samples.set_index("sample")[["group", "day"]].copy()
    meta["group"] = meta["group"].str.replace(" | ", "_", regex=False).str.replace("-", "", regex=False)
    contrasts = {
        "Sleep loss in WT": ("WT_Sleepdeprived", "WT_Sleep"),
        "Sleep loss in Shank3": ("Shank3_Sleepdeprived", "Shank3_Sleep"),
        "Shank3 vs WT (sleep)": ("Shank3_Sleep", "WT_Sleep"),
        "Shank3 vs WT (sleep-deprived)": ("Shank3_Sleepdeprived", "WT_Sleepdeprived"),
    }
    models = {"ignore day": "~group", "account for day": "~day + group"}
    rows, base_mean = [], None
    for model, design in models.items():
        print(f"PyDESeq2: {design}")
        dds = DeseqDataSet(counts=counts.T, metadata=meta, design=design, quiet=True)
        dds.deseq2()
        for name, (num, ref) in contrasts.items():
            ds = DeseqStats(dds, contrast=["group", num, ref], quiet=True)
            ds.summary()
            res = ds.results_df
            res.index = res.index.astype(int)  # PyDESeq2 turns gene IDs into strings
            if base_mean is None:
                base_mean = res["baseMean"]
            n_na = res["padj"].isna().sum()
            rows.append(
                pd.DataFrame(
                    {
                        "entrez": res.index,
                        "contrast": name,
                        "model": model,
                        "log2FC": res["log2FoldChange"].values,
                        "lfcSE": res["lfcSE"].values,
                        # NA p-values (Cook's outliers) / NA padj (independent
                        # filtering) are shown to students as "not significant".
                        "pvalue": res["pvalue"].fillna(1.0).values,
                        "padj": res["padj"].fillna(1.0).values,
                    }
                )
            )
            print(f"  {name:32s} padj<0.05: {(res['padj'] < 0.05).sum():5d}  (padj NA: {n_na})")
    return pd.concat(rows, ignore_index=True), base_mean


def build_go(gene_symbols: set[str]) -> pd.DataFrame:
    ns_map = {"biological_process": "BP", "molecular_function": "MF", "cellular_component": "CC"}
    namespace, term_id = {}, None
    with open(fetch("go-basic.obo")) as fh:
        for line in fh:
            if line.startswith("id: GO:"):
                term_id = line[4:].strip()
            elif line.startswith("namespace:") and term_id:
                ns = ns_map.get(line.split(":", 1)[1].strip())
                if ns:
                    namespace[term_id] = ns
    rows = []
    with open(fetch("gprofiler_mmusculus.gmt")) as fh:
        for line in fh:
            term, name, *members = line.rstrip("\n").split("\t")
            if not term.startswith("GO:") or term not in namespace:
                continue
            members = sorted(set(members) & gene_symbols)
            if GO_MIN <= len(members) <= GO_MAX:
                rows.append((term, name, namespace[term], len(members), " ".join(members)))
    go = pd.DataFrame(rows, columns=["term_id", "name", "ontology", "n_genes", "genes"])
    print(f"GO terms kept: {len(go)}  {go['ontology'].value_counts().to_dict()}")
    return go


def gprofiler_version() -> str:
    try:
        req = urllib.request.Request(
            "https://biit.cs.ut.ee/gprofiler/api/util/data_versions",
            headers={"User-Agent": "omics-teaching/1.0"},
        )
        with urllib.request.urlopen(req, timeout=20) as r:
            info = json.load(r)
        go = info["sources"]["GO:BP"]["version"].replace("\n", "; ")
        return f"g:Profiler {info['gprofiler_version']} (Ensembl {info['biomart_version']}; GO {go})"
    except Exception as exc:  # version info is nice-to-have only
        return f"g:Profiler version unavailable ({exc})"


def main() -> None:
    raw = pd.read_csv(fetch("counts.txt.gz"), sep="\t", index_col=0)
    raw = raw.drop(columns="Length")
    raw.columns = raw.columns.str.strip()  # GEO column names carry trailing spaces
    raw.index.name = "entrez"

    samples = build_samples(raw, fetch("series_matrix.txt.gz"))
    genes = build_genes(raw)

    cpm = raw / raw.sum() * 1e6
    genes["expressed"] = ((cpm > MIN_CPM).sum(axis=1) >= MIN_SAMPLES).values
    expressed = raw.loc[genes["expressed"]]
    print(f"genes: {len(raw)} total, {len(expressed)} expressed (CPM>{MIN_CPM} in >={MIN_SAMPLES} samples)")

    de, base_mean = run_deseq2(expressed, samples)
    genes["baseMean"] = base_mean.reindex(genes.index)

    go = build_go(set(genes["symbol"]))

    markers = pd.DataFrame({"symbol": list(MARKERS), "note": list(MARKERS.values())})
    sym2id = genes.reset_index().set_index("symbol")["entrez"]
    markers["entrez"] = markers["symbol"].map(sym2id)
    markers = markers[markers["entrez"].isin(expressed.index)]

    # ---- self-checks -------------------------------------------------------
    assert len(samples) == 20 and samples.groupby("group").size().eq(5).all()
    assert samples.notna().all().all()
    wt = de[(de.contrast == "Sleep loss in WT") & (de.model == "ignore day")].set_index("entrez")
    for sym, sign in [("Arc", 1), ("Fos", 1), ("Cirbp", -1)]:
        assert np.sign(wt.loc[sym2id[sym], "log2FC"]) == sign, f"{sym} direction unexpected"

    # ---- write -------------------------------------------------------------
    samples.to_csv(OUT / "samples.csv", index=False)
    genes.reset_index().to_csv(OUT / "genes.csv.gz", index=False, float_format="%.4g")
    expressed.to_csv(OUT / "counts.csv.gz", float_format="%d")
    de.to_csv(OUT / "de_results.csv.gz", index=False, float_format="%.4g")
    go.to_csv(OUT / "go_terms.csv.gz", index=False)
    markers.to_csv(OUT / "marker_genes.csv", index=False)
    (OUT / "data_versions.txt").write_text(
        f"Built {dt.date.today().isoformat()} by prepare_data.py\n"
        f"Counts: GEO GSE113754 (Ingiosi et al. 2019 eLife 8:e42819)\n"
        f"Gene symbols: NCBI Mus_musculus.gene_info (downloaded {dt.date.today().isoformat()})\n"
        f"GO sets: {gprofiler_version()}; ontology namespaces from go-basic.obo\n"
        f"Expression filter: CPM>{MIN_CPM} in >={MIN_SAMPLES} samples; GO term size {GO_MIN}-{GO_MAX}\n"
        f"DE: PyDESeq2 0.5.4, designs ~group and ~day + group\n"
    )
    for f in sorted(OUT.iterdir()):
        print(f"  {f.name:24s} {f.stat().st_size / 1e6:6.2f} MB")


if __name__ == "__main__":
    main()
