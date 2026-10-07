# Sleep deprivation RNA-seq practical

An interactive [marimo](https://marimo.io) notebook for the MSc course *Omics for
Neuroscience*. Students with no programming background explore real bulk RNA-seq data
using only sliders, dropdowns, switches and buttons. The notebook teaches:

1. **Multiple hypothesis testing**: false positives and the false discovery rate (FDR).
2. **Differential expression**: effect size vs significance vs biological relevance.
3. **GO enrichment**: over-representation analysis, backgrounds and thresholds.
4. **PCA**: what PCs capture, preprocessing, batch effects and circular analysis.

**Data:** mouse prefrontal cortex after 5 h sleep deprivation vs sleep, in wild-type and
*Shank3* mutant mice, with 5 mice per group. Source: GEO
[GSE113754](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE113754), from Ingiosi *et al.*
(2019) *eLife* 8:e42819.

---

## For students

Open the link your instructor gives you (Chrome or Edge recommended). The first load takes
about 20–30 s while Python starts inside your browser. Nothing needs installing. Download
your answers at the end (Section 6), because they are lost if you reload the page.

## For instructors

All commands need [uv](https://docs.astral.sh/uv/). Each notebook declares its own
dependencies (PEP 723), so no environment setup is needed.

| Task | Command |
|---|---|
| Preview the student view (code hidden) | `uvx marimo@0.25.1 run --sandbox sleep_omics.py` |
| Edit the notebook | `uvx marimo@0.25.1 edit --sandbox sleep_omics.py` |
| Rebuild the data in `public/` | `uv run prepare_data.py` (≈1 min; downloads to `data_raw/`) |
| Check the notebook | `uvx marimo@0.25.1 check sleep_omics.py` and `uv run sleep_omics.py` |
| Build the in-browser site | `uvx marimo@0.25.1 export html-wasm sleep_omics.py -o site --mode run --no-show-code` |
| Test the site locally | `python -m http.server --directory site` → open http://localhost:8000 |

**Publishing.** Push this folder to a GitHub repository and enable *Settings → Pages →
Source: GitHub Actions*. The workflow in `.github/workflows/pages.yml` builds the site and
publishes it on every push. Share the Pages URL with students. Alternatives: any static web
server (upload the contents of `site/`), or an *Open in molab* link
(`https://molab.marimo.io/github/<owner>/<repo>/blob/main/sleep_omics.py`, add `/wasm`
to run it without a login). The site must be served over HTTP; opening `index.html`
directly from disk does not work.

### Files

```
prepare_data.py   instructor-only data build (GEO → clean tables, PyDESeq2, GO sets)
sleep_omics.py    the marimo notebook (students see the app view, code hidden)
theme.css         page fonts (IBM Plex, embedded) and heading/link colours; chart colours
                  are set in the palette cell near the top of sleep_omics.py
public/           small data files the notebook loads (~5 MB, committed)
  samples.csv  genes.csv.gz  counts.csv.gz  de_results.csv.gz
  go_terms.csv.gz  marker_genes.csv  data_versions.txt
  fonts/OFL.txt   licence for the IBM Plex fonts embedded in theme.css
worksheet.md      printable student question sheet
data_raw/         download cache (git-ignored)
site/             built web version (git-ignored; built by CI)
```

### What runs live and what is pre-computed

| Live in the browser (fast) | Pre-computed by `prepare_data.py` |
|---|---|
| Simulated t-tests, BH correction (Section 2a) | DESeq2 for 4 comparisons × 2 models (Section 3) |
| Welch t-tests on log2 CPM, label shuffles (2b), replicate subsampling (3b) | Gene symbols, expression filter |
| GO hypergeometric tests + BH (Section 4) | GO gene sets (g:Profiler, propagated) |
| PCA by SVD, the circularity demo (Section 5) | |

In each section the controls are grouped under "Analysis settings" (they change the calculation:
comparison, model, normalisation, test) and "Display settings" (they only filter or display
existing results).

## Answer key

The answer key (expected values for every exercise) is kept in `instructor_notes.md`,
which is git-ignored so that it is not published with the repository.

## Assumptions and caveats

- **Experiment day.** GEO lists the mice as replicates 1–5 of an experiment "repeated over
  the course of 5 days". The notebook assumes replicate *k* = day *k*. The PCA supports this
  (PC3 is 76% explained by it), but GEO doesn't state it. Confirm with the authors before
  presenting it as fact.
- **Methods differ from the paper.** The paper used edgeR + RUVSeq (k = 5). The notebook uses
  PyDESeq2 (Section 3) and Welch t-tests (Sections 2b and 3b). The numbers differ from the
  publication; the notebook says so in its "About the data" box.
- **GO annotation drifts between releases.** The GO file in `public/` is frozen (g:Profiler
  e114_eg62_p19, GO release 2026-01-23). Rebuilding the data may shift term names and
  counts slightly.
- **Browser delivery** needs internet on first load (Pyodide comes from a CDN). Check that the
  campus network allows `cdn.jsdelivr.net`. Offline fallback: `uvx marimo run --sandbox
  sleep_omics.py` on lab machines.
- **Answers are not stored** anywhere. Students download a `.txt` file and hand it in (e.g. via
  the LMS).
- About 1,300 Entrez IDs in the 2018 count table no longer exist in NCBI Gene. They are
  labelled `Entrez:<id>`.

### Design notes

- Data are stored as `.csv.gz` rather than Parquet, so the browser version does not need to
  load `pyarrow`.
- PCA uses log2(CPM + 1) rather than a DESeq2 variance-stabilising transformation; it is easier
  to explain and gives the same picture for these data.
- There are no per-section reset buttons; reloading the page resets all controls.
