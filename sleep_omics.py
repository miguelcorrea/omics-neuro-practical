# /// script
# requires-python = ">=3.10"
# dependencies = [
#     "marimo==0.25.1",
#     "altair",
#     "numpy",
#     "pandas",
#     "scipy",
# ]
# ///

import marimo

__generated_with = "0.25.1"
app = marimo.App(
    width="medium",
    app_title="Sleep deprivation RNA-seq practical",
    css_file="theme.css",
)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    # Sleep deprivation and gene expression in mouse prefrontal cortex

    *Omics for Neuroscience, computer practical*

    In this study, mice were kept awake for 5 hours at the start of their normal sleep period
    by gentle handling, while control mice were left to sleep. The prefrontal cortex was then
    collected and gene expression was measured by bulk RNA-sequencing. Half of the mice carried
    a mutation in *Shank3*, a synaptic gene associated with autism.

    We will ask which genes respond to sleep loss, what those genes have in common, and whether
    the *Shank3* mutation changes anything. Along the way we will look at how to avoid false
    conclusions when testing about 15,000 genes at once.

    Data: GEO [GSE113754](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE113754);
    Ingiosi et al. (2019), [eLife 8:e42819](https://doi.org/10.7554/eLife.42819).
    """)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## How to use this notebook

    You do not need to write any code. Each section has controls (sliders, menus, switches and
    buttons), and the results update as soon as you change them.

    The controls are grouped under two headings. **Analysis settings** change the calculation
    itself, for example which groups are compared, the statistical model or the normalisation.
    **Display settings** only filter or display results that have already been calculated.

    Each section starts with a prediction. Make it before you change anything, then work through
    the questions and write your answer in the box at the end of the section. In Section 6 you
    can download all your answers as a text file. Answers are not stored, so download them before
    you close or reload the page; reloading also resets all the controls.
    """)
    return


@app.cell(hide_code=True)
def _():
    import io

    import altair as alt
    import marimo as mo
    import numpy as np
    import pandas as pd
    from scipy import sparse, stats

    return alt, io, mo, np, pd, sparse, stats


@app.cell(hide_code=True)
def _(alt):
    # Colour scheme: sleep = indigo (night), sleep-deprived = amber (day). Colours follow
    # meaning, never rank. Each categorical pair was checked for colour-blind separation;
    # page fonts and heading colours are set in theme.css.
    C_SLEEP, C_SD = "#4a4fb0", "#b87800"
    C_UP, C_DOWN, C_NS = "#c9302c", "#2b6cb0", "#cbd5e1"
    C_TRUE, C_NULL = "#178f80", "#94a3b8"
    C_WT, C_MUT = "#178f80", "#9c4fb0"
    C_HIGHLIGHT, C_ACCENT = "#178f80", "#4a4fb0"
    DAY_COLORS = ["#0072B2", "#E69F00", "#009E73", "#CC79A7", "#56B4E9"]  # Okabe-Ito; days are also labelled
    SEQ_RANGE = ["#e0e3ff", "#2e3192"]  # sequential (light to dark indigo)
    INK, INK_MUTED, INK_STRONG = "#334155", "#64748b", "#1e293b"
    GRID, AXIS = "#e2e8f0", "#cbd5e1"

    COND_SCALE = alt.Scale(domain=["Sleep", "Sleep-deprived"], range=[C_SLEEP, C_SD])
    GROUP_ORDER = ["WT | Sleep", "WT | Sleep-deprived", "Shank3 | Sleep", "Shank3 | Sleep-deprived"]

    def style(chart):
        """Recessive axes and grid, readable labels, legend on top."""
        return (
            chart.configure(font="IBM Plex Sans")  # same typeface as the page (theme.css)
            .configure_axis(
                gridColor=GRID,
                domainColor=AXIS,
                tickColor=AXIS,
                labelColor=INK,
                titleColor=INK,
                labelFontSize=12,
                titleFontSize=12,
                titleFontWeight="normal",
            )
            .configure_view(strokeWidth=0)
            .configure_legend(orient="top", labelColor=INK, titleColor=INK, labelFontSize=12, titleFontSize=12)
            .configure_title(color=INK, fontSize=13, anchor="start", subtitleColor=INK_MUTED)
        )

    return (
        COND_SCALE,
        C_ACCENT,
        C_DOWN,
        C_HIGHLIGHT,
        C_MUT,
        C_NS,
        C_NULL,
        C_SD,
        C_TRUE,
        C_UP,
        C_WT,
        DAY_COLORS,
        GROUP_ORDER,
        INK_MUTED,
        INK_STRONG,
        SEQ_RANGE,
        style,
    )


@app.cell(hide_code=True)
def _(mo):
    def control_groups(analysis=(), display=()):
        """Lay out controls under 'Analysis settings' and 'Display settings'.

        Each argument is a list of rows; each row is a list of controls.
        """
        _parts = []
        for _title, _rows in (("Analysis settings", analysis), ("Display settings", display)):
            if _rows:
                _parts.append(mo.md(f"**{_title}**"))
                _parts.extend(mo.hstack(_row, justify="start", gap=2, wrap=True) for _row in _rows)
        return mo.vstack(_parts)

    return (control_groups,)


@app.cell(hide_code=True)
async def _(io, mo, pd):
    async def _read_bytes(name):
        _path = mo.notebook_location() / "public" / name
        _url = str(_path)
        if _url.startswith(("http://", "https://")):  # running in the browser (WebAssembly)
            from pyodide.http import pyfetch

            _resp = await pyfetch(_url)
            return await _resp.bytes()
        with open(_path, "rb") as _fh:
            return _fh.read()

    async def _csv(name, **kw):
        _data = await _read_bytes(name)
        return pd.read_csv(io.BytesIO(_data), compression="gzip" if name.endswith(".gz") else None, **kw)

    samples = await _csv("samples.csv")
    genes = await _csv("genes.csv.gz", index_col="entrez")
    counts = (await _csv("counts.csv.gz", index_col="entrez"))[samples["sample"]]
    de_all = await _csv("de_results.csv.gz")
    go_terms = await _csv("go_terms.csv.gz")
    markers = await _csv("marker_genes.csv")
    data_versions = (await _read_bytes("data_versions.txt")).decode()
    return counts, data_versions, de_all, genes, go_terms, markers, samples


@app.cell(hide_code=True)
def _(counts, genes, np, samples):
    # Normalisation: counts per million (CPM) using each mouse's total reads,
    # then log2(CPM + 1). Only genes detected in the experiment are kept.
    lib_size = samples.set_index("sample")["library_size"]
    cpm = counts / lib_size * 1e6
    logcpm = np.log2(cpm + 1)
    symbol = genes["symbol"]
    entrez_of = dict(zip(genes["symbol"], genes.index))
    return cpm, entrez_of, logcpm, symbol


@app.cell(hide_code=True)
def _(np, stats):
    def bh(p):
        """Benjamini-Hochberg adjusted p-values (q-values)."""
        p = np.asarray(p, dtype=float)
        n = p.size
        if n == 0:
            return p
        order = np.argsort(p)
        ranked = p[order] * n / np.arange(1, n + 1)
        q = np.minimum.accumulate(ranked[::-1])[::-1]
        out = np.empty(n)
        out[order] = np.minimum(q, 1.0)
        return out

    def ttest_rows(a, b):
        """Welch t-test for every row (gene) at once. Returns t and p."""
        t, p = stats.ttest_ind(a, b, axis=1, equal_var=False)
        return np.nan_to_num(t), np.nan_to_num(p, nan=1.0)

    def fmt(n):
        return f"{n:,}"

    return bh, fmt, ttest_rows


@app.cell(hide_code=True)
def _(
    alt,
    COND_SCALE,
    counts,
    cpm,
    GROUP_ORDER,
    INK_STRONG,
    logcpm,
    np,
    samples,
    symbol,
):
    def gene_dotplot(entrez, scale="log2 CPM", title_extra=""):
        """One dot per mouse, grouped by genotype x condition, with group means."""
        _source = {"Raw counts": counts, "CPM": cpm, "log2 CPM": logcpm}[scale]
        _d = samples.copy()
        _d["value"] = _source.loc[entrez, _d["sample"]].values
        _rng = np.random.default_rng(int(entrez) % 9973)
        _d["jitter"] = _rng.uniform(-0.25, 0.25, len(_d))
        _ylab = {"Raw counts": "Raw read count", "CPM": "Counts per million (CPM)", "log2 CPM": "log2(CPM + 1)"}[scale]
        _base = alt.Chart(_d).encode(x=alt.X("group:N", sort=GROUP_ORDER, title=None, axis=alt.Axis(labelAngle=0, labelLimit=180)))
        _pts = _base.mark_circle(size=90, opacity=0.9, stroke="white", strokeWidth=1).encode(
            xOffset=alt.XOffset("jitter:Q", scale=alt.Scale(domain=[-1, 1])),
            y=alt.Y("value:Q", title=_ylab, scale=alt.Scale(zero=scale != "log2 CPM")),
            color=alt.Color("condition:N", scale=COND_SCALE, title="Condition"),
            tooltip=["sample", "genotype", "condition", "day", alt.Tooltip("value:Q", format=",.2f", title=scale)],
        )
        _mean = _base.mark_tick(color=INK_STRONG, thickness=2, size=36).encode(y="mean(value):Q", tooltip=[alt.Tooltip("mean(value):Q", format=",.2f", title="group mean")])
        return (_pts + _mean).properties(
            width=520, height=260, title=alt.TitleParams(f"{symbol.get(entrez, entrez)}{title_extra}", subtitle="Each dot is one mouse; the black bar is the group mean")
        )

    return (gene_dotplot,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ---
    ## 1. The experiment and the data

    The design has two genotypes (wild type and *Shank3* mutant) and two conditions (sleep and
    sleep deprivation), with 5 male mice aged 8–12 weeks in each of the four groups. Each mouse
    is a **biological replicate**. The experiment was repeated on 5 separate days; we assume one
    mouse from each group per day (see "Which mouse was processed on which day?" in Section 3).

    | | Sleep (home cage) | Sleep-deprived (5 h) |
    |---|---|---|
    | Wild type (WT) | 5 mice | 5 mice |
    | *Shank3* mutant | 5 mice | 5 mice |

    RNA-seq gives a **count** for each gene in each mouse: the number of sequencing reads assigned
    to that gene. The mice were not all sequenced to the same depth, however. The total number of
    reads per mouse (the library size) varies, and a mouse with twice as many reads will have
    roughly twice the count for every gene, whatever its biology.

    To compare mice fairly, each count is divided by that mouse's library size and multiplied by
    one million. The result is **counts per million (CPM)**: the number of reads a gene would have
    if every mouse had been sequenced to exactly one million reads.

    Expression values span a huge range, from a handful of reads to hundreds of thousands, and
    changes in expression are multiplicative: a gene doubles or halves. We therefore often use
    **log2 CPM**, more precisely log2(CPM + 1); the +1 avoids taking the log of zero. On a log2
    scale a doubling is always a step of +1, whether the gene is rare or abundant, so lowly and
    highly expressed genes can be compared on the same axis. The gene plot below lets you switch
    between the three scales.
    """)
    return


@app.cell(hide_code=True)
def _(COND_SCALE, GROUP_ORDER, alt, mo, samples, style):
    _d = samples.assign(million_reads=samples["library_size"] / 1e6)
    _d["order"] = _d["group"].map({g: i for i, g in enumerate(GROUP_ORDER)})
    _bars = (
        alt.Chart(_d)
        .mark_bar(cornerRadiusTopLeft=4, cornerRadiusTopRight=4, width={"band": 0.8})
        .encode(
            x=alt.X("sample:N", sort=alt.EncodingSortField("order"), title="Mouse (sample)", axis=alt.Axis(labelAngle=-60)),
            y=alt.Y("million_reads:Q", title="Total reads (millions)"),
            color=alt.Color("condition:N", scale=COND_SCALE, title="Condition"),
            tooltip=["sample", "group", "day", alt.Tooltip("million_reads:Q", format=".1f", title="million reads")],
        )
        .properties(
            width=620,
            height=220,
            title=alt.TitleParams(
                "Library size per mouse",
                subtitle="Sample names give genotype (WT or S3), condition (HC = home-cage sleep, SD = sleep-deprived) and replicate 1–5",
            ),
        )
    )
    mo.vstack(
        [
            style(_bars),
            mo.accordion({"Sample table": mo.ui.table(samples, selection=None, show_download=False)}),
        ]
    )
    return


@app.cell(hide_code=True)
def _(counts, markers, mo, symbol):
    _all = {symbol[e]: int(e) for e in counts.index}  # detected genes only
    gene_pick = mo.ui.dropdown(
        options=dict(sorted(_all.items(), key=lambda kv: kv[0].lower())),
        value="Arc",
        searchable=True,
        label="Gene (type to search)",
    )
    scale_pick = mo.ui.radio(["Raw counts", "CPM", "log2 CPM"], value="Raw counts", inline=True, label="Expression scale")
    pred_1 = mo.ui.radio(
        ["Higher", "Lower", "Unchanged"],
        label="Before you change anything: *Arc* is induced by neuronal activity. After 5 hours awake, will its expression be higher, lower or unchanged?",
    )
    mo.vstack([pred_1, mo.md("Genes to try: " + ", ".join(f"*{s}*" for s in markers["symbol"]))])
    return gene_pick, pred_1, scale_pick


@app.cell(hide_code=True)
def _(control_groups, gene_dotplot, gene_pick, markers, mo, scale_pick, style):
    _note = markers.set_index("symbol")["note"].to_dict()
    _sym = gene_pick.selected_key
    mo.vstack(
        [
            control_groups(analysis=[[scale_pick]], display=[[gene_pick]]),
            mo.md(f"*{_sym}*: {_note[_sym]}") if _sym in _note else mo.md(""),
            style(gene_dotplot(gene_pick.value, scale_pick.value)),
        ]
    )
    return


@app.cell(hide_code=True)
def _(mo, pred_1):
    _fb = (
        "Correct." if pred_1.value == "Higher" else ("Not quite. Look at the plot." if pred_1.value else "Make your prediction above first.")
    )
    mo.vstack(
        [
            mo.md(
                f"""
    #### Questions
    1. Was your prediction for *Arc* right? {_fb}
    2. Select *Shank3* and keep the scale on raw counts. The sleep-deprived mice appear to have
       less *Shank3*. Compare with the library sizes above, then switch to CPM. What changed,
       and why?
    3. Select *Fos* on the log2 CPM scale. How much do the five mice within a group differ? Why
       does the experiment need several mice per group?
    4. Compare *Cirbp* and *Gapdh*. Which one changes with sleep loss?
    """
            ),
            mo.accordion(
                {
                    "Show explanation": mo.md(
                        """
    - *Arc* increases about 4-fold after sleep loss (a log2 fold change of about 2). It is a
      well-characterised activity-regulated gene.
    - The sleep-deprived WT mice have smaller libraries (23–31 million reads, compared with
      36–65 million for the sleeping WT mice), so their raw counts are lower for almost every
      gene. CPM divides each count by the library size, and the apparent difference in *Shank3*
      disappears. Without normalisation, thousands of genes would wrongly appear down-regulated.
    - Mice within a group differ because of biological variability. A difference between groups
      is convincing only when it is large compared with this spread, and the spread can only be
      estimated from several biological replicates.
    - *Cirbp* decreases after sleep loss; *Gapdh* does not change.
    """
                    )
                }
            ),
        ]
    )
    return


@app.cell(hide_code=True)
def _(mo):
    ans_1 = mo.ui.text_area(placeholder="Section 1: normalisation and replicates", rows=3, full_width=True, label="Your answer")
    ans_1
    return (ans_1,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ---
    ## 2. Testing thousands of genes

    For each gene we can run a statistical test, such as a t-test, of whether average expression
    differs between two groups. The test returns a **p-value**: the probability, assuming there is
    no true difference and the test's assumptions hold, of obtaining a test statistic at least as
    extreme as the one observed. Because thousands of genes are tested at once, p-values are
    adjusted for multiple testing (typically by controlling the **false discovery rate**, FDR),
    and a gene is called significant if its adjusted p-value falls below a chosen threshold, such
    as 0.05.

    Why adjust? Even when nothing has changed, a raw p-value falls below 0.05 in about 5% of tests
    by chance. For a single test this risk is small; across 15,000 genes it produces hundreds of
    false positives. The xkcd comic [Significant](https://xkcd.com/882/) makes the same point with
    jelly beans. This section shows the problem first in simulated data, then in the real mice.

    ### 2a. Simulation

    In this simulation the computer generates expression data for two groups of mice, and you
    choose how many genes truly change. Because the data are simulated, we know which significant
    genes are real effects and which are **false positives**.
    """)
    return


@app.cell(hide_code=True)
def _(mo):
    sim_genes = mo.ui.slider(100, 20000, step=100, value=10000, debounce=True, show_value=True, label="Number of genes tested")
    sim_true = mo.ui.slider(0, 30, step=1, value=0, debounce=True, show_value=True, label="Genes that truly change (%)")
    sim_effect = mo.ui.slider(0.25, 2.0, step=0.25, value=1.0, debounce=True, show_value=True, label="Size of the true change (log2 fold change)")
    sim_n = mo.ui.slider(2, 10, step=1, value=5, debounce=True, show_value=True, label="Mice per group")
    sim_method = mo.ui.radio(["Raw p-value", "BH-adjusted p-value (FDR)"], value="Raw p-value", label="Call a gene significant using")
    sim_alpha = mo.ui.slider(steps=[0.001, 0.01, 0.05, 0.1, 0.2], value=0.05, show_value=True, label="Threshold")
    sim_rerun = mo.ui.button(value=0, on_click=lambda v: v + 1, label="New random data", kind="neutral")
    pred_2 = mo.ui.radio(
        ["About 0", "About 50", "About 500", "About 5,000"],
        label="Before you change anything: 10,000 genes are tested, none of them truly changes, and genes are called significant at raw *p* < 0.05. How many will be significant?",
    )
    pred_2
    return (
        pred_2,
        sim_alpha,
        sim_effect,
        sim_genes,
        sim_method,
        sim_n,
        sim_rerun,
        sim_true,
    )


@app.cell(hide_code=True)
def _(bh, np, sim_effect, sim_genes, sim_n, sim_rerun, sim_true, ttest_rows):
    # Each gene: log2 expression ~ Normal(mean, sd = 0.5) in each mouse.
    # "True" genes are shifted by the chosen log2 fold change in group B.
    _rng = np.random.default_rng(2024 + sim_rerun.value)
    _G, _n = sim_genes.value, sim_n.value
    sim_is_true = np.zeros(_G, dtype=bool)
    sim_is_true[: int(round(_G * sim_true.value / 100))] = True
    _a = _rng.normal(0, 0.5, size=(_G, _n))
    _b = _rng.normal(0, 0.5, size=(_G, _n)) + np.where(sim_is_true, sim_effect.value, 0.0)[:, None]
    _, sim_p = ttest_rows(_b, _a)
    sim_q = bh(sim_p)
    return sim_is_true, sim_p, sim_q


@app.cell(hide_code=True)
def _(
    alt,
    C_NULL,
    C_TRUE,
    control_groups,
    fmt,
    INK_STRONG,
    mo,
    np,
    pd,
    sim_alpha,
    sim_effect,
    sim_genes,
    sim_is_true,
    sim_method,
    sim_n,
    sim_p,
    sim_q,
    sim_rerun,
    sim_true,
    style,
):
    _use_q = sim_method.value.startswith("BH")
    _called = (sim_q if _use_q else sim_p) < sim_alpha.value
    _tp = int((_called & sim_is_true).sum())
    _fp = int((_called & ~sim_is_true).sum())
    _fn = int((~_called & sim_is_true).sum())
    _n_called = _tp + _fp
    _fdp = _fp / _n_called if _n_called else 0.0

    _hist = pd.DataFrame({"p": sim_p, "truth": np.where(sim_is_true, "Truly changed", "No real change")})
    _hist["bin"] = np.minimum((_hist["p"] * 20).astype(int), 19) / 20
    _hist = _hist.groupby(["bin", "truth"]).size().reset_index(name="genes")
    _hist["bin_end"] = _hist["bin"] + 0.05
    _null_level = (~sim_is_true).sum() / 20
    _chart = alt.Chart(_hist).mark_bar(stroke="white", strokeWidth=1).encode(
        x=alt.X("bin:Q", bin="binned", title="p-value", scale=alt.Scale(domain=[0, 1])),
        x2="bin_end",
        y=alt.Y("genes:Q", stack=True, title="Number of genes"),
        color=alt.Color("truth:N", scale=alt.Scale(domain=["No real change", "Truly changed"], range=[C_NULL, C_TRUE]), title=None),
        order=alt.Order("truth:N", sort="descending"),
        tooltip=["truth", "genes", alt.Tooltip("bin:Q", title="p from")],
    )
    _flat = alt.Chart(pd.DataFrame({"y": [_null_level]})).mark_rule(strokeDash=[4, 4], color=INK_STRONG).encode(y="y:Q")
    _hist_chart = (_chart + _flat).properties(
        width=430, height=230, title=alt.TitleParams("Distribution of p-values", subtitle="Dashed line: the level expected from genes with no real change")
    )

    _rule = f"{'BH-adjusted p' if _use_q else 'raw p'} < {sim_alpha.value}"
    _tiles = mo.hstack(
        [
            mo.stat(fmt(_n_called), label="Called significant", caption=_rule, bordered=True),
            mo.stat(fmt(_fp), label="False positives", caption="no real change, but called", bordered=True),
            mo.stat(fmt(_tp), label="True positives", caption="real changes called", bordered=True),
            mo.stat(fmt(_fn), label="Missed", caption="real changes not called", bordered=True),
            mo.stat(f"{_fdp:.0%}", label="False discovery proportion", caption="false positives / called", bordered=True),
        ],
        widths="equal",
    )
    _sentence = (
        f"{fmt(sim_genes.value)} genes tested with {sim_n.value} mice per group, of which "
        f"{fmt(int(sim_is_true.sum()))} truly change (log2 fold change {sim_effect.value}). "
        f"Using {_rule}, {fmt(_n_called)} genes were called significant, and {fmt(_fp)} of them "
        f"({_fdp:.0%}) are false positives."
    )
    if sim_true.value == 0 and not _use_q:
        _sentence += f" With no real changes, about {sim_alpha.value:.0%} of genes pass the threshold by chance."
    mo.vstack(
        [
            control_groups(
                analysis=[[sim_genes, sim_true], [sim_effect, sim_n], [sim_method, sim_rerun]],
                display=[[sim_alpha]],
            ),
            mo.md(f"<small>Random seed: {2024 + sim_rerun.value}</small>"),
            _tiles,
            mo.md(_sentence),
            style(_hist_chart),
        ]
    )
    return


@app.cell(hide_code=True)
def _(mo, pred_2):
    _fb = {"About 500": "Correct: 5% of 10,000 is 500.", None: "Make your prediction above first."}.get(
        pred_2.value, "Not quite. Run the simulation with 0% true changes and count."
    )
    mo.vstack(
        [
            mo.md(
                f"""
    #### Questions
    1. With 0% true changes and raw p-values, how many genes are called significant? Press
       "New random data" a few times. {_fb}
    2. Keep 0% true changes and switch to BH-adjusted p-values. What happens to the false
       positives?
    3. Set 10% of genes to truly change. Compare raw and BH-adjusted p-values at a threshold of
       0.05, then try BH thresholds of 0.01, 0.05 and 0.2. Does the false discovery proportion
       stay close to the threshold?
    4. With BH at 0.05 and 10% true changes, vary the number of mice per group from 2 to 10.
       What happens to the number of missed genes?
    5. How do the truly changed genes alter the shape of the p-value histogram?
    """
            ),
            mo.accordion(
                {
                    "Show explanation": mo.md(
                        """
    - When nothing changes, p-values are spread evenly between 0 and 1, so 5% fall below 0.05:
      about 500 of 10,000 genes, all of them false positives.
    - The Benjamini–Hochberg (BH) procedure adjusts the p-values so that, among the genes called
      significant, the expected proportion of false positives stays below the chosen threshold.
      This proportion is the false discovery rate (FDR). With no real effects, BH usually calls no
      genes at all.
    - With real effects, BH keeps the proportion of false positives near the target. The FDR is a
      property of the whole list, not a guarantee about any single gene.
    - More replicates give more power, so fewer real changes are missed at the same FDR.
    - Real effects add a peak of small p-values near 0 on top of the flat background.
    """
                    )
                }
            ),
        ]
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ### 2b. Real data

    Here we run a t-test on log2 CPM values for each of the ~14,700 genes detected in the
    experiment. Two comparisons are available: sleep deprivation, where we expect a large effect,
    and genotype, where any effect is likely to be small. The shuffle switch assigns the mice to
    the two groups at random. This removes any real link between group and expression, so the
    result shows what the real data look like when nothing is going on.

    A t-test on log2 CPM is a shortcut used here because it is fast and easy to follow. A real
    analysis would use a method designed for RNA-seq counts, such as DESeq2 (Section 3), edgeR or
    limma-voom.
    """)
    return


@app.cell(hide_code=True)
def _(control_groups, mo):
    real_contrast = mo.ui.dropdown(
        {
            "Shank3 vs WT (sleeping mice)": ("genotype", "Sleep"),
            "Sleep-deprived vs sleeping (WT mice)": ("condition", "WT"),
        },
        value="Shank3 vs WT (sleeping mice)",
        label="Comparison",
    )
    real_shuffle = mo.ui.switch(label="Shuffle group labels")
    real_reshuffle = mo.ui.button(value=0, on_click=lambda v: v + 1, label="New shuffle")
    real_alpha = mo.ui.slider(steps=[0.001, 0.01, 0.05, 0.1], value=0.05, show_value=True, label="Threshold")
    control_groups(analysis=[[real_contrast], [real_shuffle, real_reshuffle]], display=[[real_alpha]])
    return real_alpha, real_contrast, real_reshuffle, real_shuffle


@app.cell(hide_code=True)
def _(
    bh,
    logcpm,
    np,
    real_contrast,
    real_reshuffle,
    real_shuffle,
    samples,
    ttest_rows,
):
    _factor, _subset = real_contrast.value
    if _factor == "genotype":
        _s = samples[samples["condition"] == _subset]
        _is_b = (_s["genotype"] == "Shank3").to_numpy()
    else:
        _s = samples[samples["genotype"] == _subset]
        _is_b = (_s["condition"] == "Sleep-deprived").to_numpy()
    if real_shuffle.value:
        _is_b = np.random.default_rng(7 + real_reshuffle.value).permutation(_is_b)
    _x = logcpm[_s["sample"]].to_numpy()
    _, real_p = ttest_rows(_x[:, _is_b], _x[:, ~_is_b])
    real_q = bh(real_p)
    real_groups = list(zip(_s["sample"], np.where(_is_b, "group B", "group A")))
    return real_groups, real_p, real_q


@app.cell(hide_code=True)
def _(
    alt,
    C_NULL,
    fmt,
    INK_STRONG,
    mo,
    np,
    pd,
    real_alpha,
    real_contrast,
    real_groups,
    real_p,
    real_q,
    real_shuffle,
    style,
):
    _raw = int((real_p < real_alpha.value).sum())
    _adj = int((real_q < real_alpha.value).sum())
    _expected = real_alpha.value * len(real_p)
    _h = pd.DataFrame({"bin": np.minimum((real_p * 20).astype(int), 19) / 20}).value_counts().reset_index(name="genes")
    _h["bin_end"] = _h["bin"] + 0.05
    _chart = (
        alt.Chart(_h)
        .mark_bar(color=C_NULL, stroke="white", strokeWidth=1)
        .encode(
            x=alt.X("bin:Q", bin="binned", title="p-value", scale=alt.Scale(domain=[0, 1])),
            x2="bin_end",
            y=alt.Y("genes:Q", title="Number of genes"),
            tooltip=["genes", alt.Tooltip("bin:Q", title="p from")],
        )
    )
    _flat = alt.Chart(pd.DataFrame({"y": [len(real_p) / 20]})).mark_rule(strokeDash=[4, 4], color=INK_STRONG).encode(y="y:Q")
    _title = real_contrast.selected_key + (", shuffled labels" if real_shuffle.value else "")
    _plot = (_chart + _flat).properties(width=430, height=230, title=alt.TitleParams(_title, subtitle="Dashed line: the level expected if nothing changes"))
    _groups = pd.DataFrame(real_groups, columns=["mouse", "assigned to"])
    mo.vstack(
        [
            mo.hstack(
                [
                    mo.stat(fmt(len(real_p)), label="Genes tested", bordered=True),
                    mo.stat(fmt(_raw), label=f"Raw p < {real_alpha.value}", caption=f"about {fmt(int(_expected))} expected by chance", bordered=True),
                    mo.stat(fmt(_adj), label=f"BH-adjusted p < {real_alpha.value}", caption="after correction for multiple testing", bordered=True),
                ],
                widths="equal",
            ),
            style(_plot),
            mo.accordion({"Group assignment for each mouse": mo.ui.table(_groups, selection=None, show_download=False, page_size=10)}),
        ]
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.vstack(
        [
            mo.md(
                """
    #### Questions
    1. For Shank3 vs WT with the real labels, how many genes have raw *p* < 0.05? How many would
       you expect by chance, and how many remain after BH correction?
    2. Turn on the shuffle and draw a few new shuffles. Compare the counts and the histogram with
       question 1. What does this suggest about the genotype genes with *p* < 0.05?
    3. Switch to sleep-deprived vs sleeping. Describe the histogram. Is the effect of sleep loss
       real?
    4. Shuffle the sleep-deprivation labels as well. What happens to the peak near zero?
    """
            ),
            mo.accordion(
                {
                    "Show explanation": mo.md(
                        """
    - The genotype comparison gives about 400–450 genes at raw *p* < 0.05, fewer than the roughly
      730 expected by chance, and none remain after BH correction. Shuffled labels give similar
      numbers. The histogram has no peak near 0; it even rises towards 1, meaning there are fewer
      small p-values than chance alone would give, a sign that this simple test is conservative
      with only 5 mice per group. Any genotype effect is too small for it to detect. Section 3
      shows that a better model can find a few genotype genes.
    - Sleep deprivation gives a large peak near p = 0, and thousands of genes remain after BH
      correction, so the effect is real and strong. Shuffling removes the peak.
    - A list of genes with *p* < 0.05 is not informative on its own. It has to be corrected for
      multiple testing and compared with what chance alone would produce.
    """
                    )
                }
            ),
        ]
    )
    return


@app.cell(hide_code=True)
def _(mo):
    ans_2 = mo.ui.text_area(placeholder="Section 2: why is correction for multiple testing needed, and what did the shuffle show?", rows=3, full_width=True, label="Your answer")
    ans_2
    return (ans_2,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ---
    ## 3. Differential expression

    For differential expression (DE) analysis we use DESeq2, a standard method for RNA-seq. Its
    input is the table of raw read counts (whole numbers for the ~14,700 detected genes), not CPM
    or log values. DESeq2 does its own normalisation for library size, by estimating a "size
    factor" for each mouse, and models the counts with a negative binomial distribution, which
    describes how much read counts vary between replicate mice. It also shares information across
    genes, so that variability can be estimated reliably with only 5 mice per group. The two
    models (with and without the experiment day) were fitted in advance. For each gene the
    results give:

    - the **log2 fold change (log2FC)**, the size and direction of the effect: +1 means doubled,
      −1 halved and 0 unchanged;
    - the **adjusted p-value (padj)**: the DESeq2 p-value after BH correction (Section 2);
    - the **baseMean**: the average of the normalised counts across all 20 mice, that is, how
      strongly the gene is expressed overall.

    The volcano plot shows all three. Click a point to plot that gene's expression in each mouse.
    """)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.accordion(
        {
            "What is a log2 fold change?": mo.md(
                """
    The fold change compares the two group averages. Taking its log2 makes increases and
    decreases symmetric:

    | change | fold change | log2FC |
    |---|---|---|
    | 4 times higher | 4 | +2 |
    | 2 times higher | 2 | +1 |
    | 20% higher | 1.2 | +0.26 |
    | no change | 1 | 0 |
    | halved | 0.5 | −1 |
    """
            ),
            "What is baseMean used for?": mo.md(
                """
    baseMean is not part of the test itself; it describes how much data the test had to work
    with. A gene with a baseMean of 5 has only a handful of reads per mouse, so a few reads more
    or less can produce a large fold change by chance, and its estimates are noisy. Analysts
    therefore often set aside very lowly expressed genes.

    In this notebook baseMean is used in two places: the "Minimum average expression" slider
    removes genes below the chosen baseMean from the selected list, and the MA plot shows it on
    the x-axis. DESeq2 also uses it internally: its "independent filtering" step sets aside genes
    with a very low baseMean before the BH correction, which reduces the number of tests.
    """
            ),
            "What is an MA plot?": mo.md(
                """
    An MA plot shows each gene's log2 fold change (M, y-axis) against its average expression (A,
    x-axis). The names come from microarray analysis, where M was the difference between two log
    intensities ("minus") and A their average ("add").

    The volcano plot shows how significant each change is; the MA plot shows how the size of the
    changes depends on how strongly a gene is expressed. Lowly expressed genes, on the left,
    scatter widely in fold change, mostly because of counting noise. Highly expressed genes, on
    the right, show smaller but more reliable changes.
    """
            ),
            "Why account for the day of the experiment?": mo.md(
                """
    The experiment was run on 5 days, with one mouse from each group per day. Mice processed on
    the same day share small technical and environmental differences, known as a **batch
    effect**. If the model is told which day each mouse came from, it can remove this variation
    and compare the groups more precisely. Because every day includes every group, day is not
    confounded with the comparisons of interest.
    """
            ),
            "Which mouse was processed on which day?": mo.md(
                """
    The GEO record describes the design as follows: "The protocol was repeated over the course of
    5 days to obtain 5 animals per group, so that 5 independent SD experiments were represented in
    each time point." The individual samples are labelled only as replicates 1–5, for example
    "Wildtype homecage control harvested at ZT hour 5, replicate 1". Neither GEO nor the paper
    states which replicate was collected on which day.

    This notebook assumes that replicate *k* of every group was collected on day *k*. This is an
    inference, not a recorded fact. The data support it: in the PCA (Section 5) the replicate
    number explains about 76% of the variation along PC3 and 35% along PC2, which is what a batch
    effect would look like. Read "day" in this notebook as "probable experiment day".
    """
            ),
        }
    )
    return


@app.cell(hide_code=True)
def _(control_groups, counts, de_all, mo, symbol):
    de_contrast = mo.ui.dropdown(list(de_all["contrast"].unique()), value="Sleep loss in WT", label="Comparison")
    de_model = mo.ui.radio(["ignore day", "account for day"], value="ignore day", inline=True, label="DESeq2 model")
    de_fdr = mo.ui.slider(steps=[0.001, 0.005, 0.01, 0.05, 0.1, 0.25], value=0.05, show_value=True, label="FDR threshold (padj below)")
    de_lfc = mo.ui.slider(0, 3, step=0.1, value=0, show_value=True, label="Minimum |log2 fold change|")
    de_minexpr = mo.ui.slider(steps=[0, 10, 50, 100, 500, 1000], value=0, show_value=True, label="Minimum average expression (baseMean)")
    de_labels = mo.ui.switch(value=True, label="Label marker genes that pass the thresholds")
    de_find = mo.ui.dropdown(sorted(symbol[counts.index], key=str.lower), value=None, searchable=True, label="Find a gene")
    control_groups(
        analysis=[[de_contrast, de_model]],
        display=[[de_fdr, de_lfc], [de_minexpr, de_labels], [de_find]],
    )
    return de_contrast, de_fdr, de_find, de_labels, de_lfc, de_minexpr, de_model


@app.cell(hide_code=True)
def _(de_all, de_contrast, de_fdr, de_lfc, de_minexpr, de_model, genes, np):
    de = de_all[(de_all["contrast"] == de_contrast.value) & (de_all["model"] == de_model.value)].copy()
    de = de.join(genes[["symbol", "baseMean", "description"]], on="entrez")
    de["neglog10_padj"] = -np.log10(de["padj"].clip(lower=1e-300))
    _pass = (de["padj"] < de_fdr.value) & (de["log2FC"].abs() >= de_lfc.value) & (de["baseMean"] >= de_minexpr.value)
    de["call"] = np.where(_pass & (de["log2FC"] > 0), "Up", np.where(_pass & (de["log2FC"] < 0), "Down", "Not selected"))
    de_settings = f"{de_contrast.value}; model: {de_model.value}; padj < {de_fdr.value}; |log2FC| ≥ {de_lfc.value}; baseMean ≥ {de_minexpr.value}"
    return de, de_settings


@app.cell(hide_code=True)
def _(
    alt,
    C_DOWN,
    C_HIGHLIGHT,
    C_NS,
    C_UP,
    de,
    de_fdr,
    de_find,
    de_labels,
    de_lfc,
    INK_STRONG,
    markers,
    mo,
    np,
    pd,
    style,
):
    _cap = 60  # very small padj values are drawn at the top edge
    # In the browser (WebAssembly) chart data is embedded as JSON, so keep only the
    # columns each chart needs, rounded, to stay well under marimo's output size limit.
    _d = pd.DataFrame(
        {
            "entrez": de["entrez"],
            "symbol": de["symbol"],
            "log2FC": de["log2FC"].round(3),
            "y": de["neglog10_padj"].clip(upper=_cap).round(2),
            "padj": de["padj"].map(lambda v: float(f"{v:.2e}")),
            "baseMean": de["baseMean"].round(0),
            "call": de["call"],
        }
    )
    _color = alt.Color("call:N", scale=alt.Scale(domain=["Up", "Down", "Not selected"], range=[C_UP, C_DOWN, C_NS]), title="Selected genes")
    _tip = ["symbol", alt.Tooltip("log2FC:Q", format=".2f"), alt.Tooltip("padj:Q", format=".2e"), alt.Tooltip("baseMean:Q", format=",.0f")]
    _pick = alt.selection_point(name="pick", fields=["entrez"], on="click", empty=False)
    _pts = (
        alt.Chart(_d)
        .mark_circle(size=22, opacity=0.65)
        .encode(
            x=alt.X("log2FC:Q", title="log2 fold change"),
            y=alt.Y("y:Q", title="−log10 adjusted p-value"),
            color=_color,
            size=alt.condition(_pick, alt.value(120), alt.value(22)),
            tooltip=_tip,
        )
        .add_params(_pick)
    )
    _hline = alt.Chart(pd.DataFrame({"y": [-np.log10(de_fdr.value)]})).mark_rule(strokeDash=[4, 4], color=INK_STRONG).encode(y="y:Q")
    _layers = [_pts, _hline]
    if de_lfc.value > 0:
        _layers.append(alt.Chart(pd.DataFrame({"x": [-de_lfc.value, de_lfc.value]})).mark_rule(strokeDash=[4, 4], color=INK_STRONG).encode(x="x:Q"))
    if de_labels.value:  # only markers that made it into the list, to limit clutter
        _m = _d[_d["symbol"].isin(markers["symbol"]) & (_d["call"] != "Not selected")]
        _layers.append(alt.Chart(_m).mark_point(shape="circle", size=60, filled=False, color=INK_STRONG, strokeWidth=1.2).encode(x="log2FC:Q", y="y:Q", tooltip=_tip))
        _layers.append(alt.Chart(_m).mark_text(align="left", dx=6, dy=-6, fontSize=11, color=INK_STRONG).encode(x="log2FC:Q", y="y:Q", text="symbol:N"))
    if de_find.value:
        _f = _d[_d["symbol"] == de_find.value]
        _layers.append(alt.Chart(_f).mark_point(shape="circle", size=220, filled=False, color=C_HIGHLIGHT, strokeWidth=2.5).encode(x="log2FC:Q", y="y:Q", tooltip=_tip))
        _layers.append(alt.Chart(_f).mark_text(align="left", dx=10, dy=-10, fontSize=13, fontWeight="bold", color=C_HIGHLIGHT).encode(x="log2FC:Q", y="y:Q", text="symbol:N"))
    volcano = mo.ui.altair_chart(
        style(
            alt.layer(*_layers).properties(
                width=560,
                height=380,
                title=alt.TitleParams("Volcano plot", subtitle=f"Dashed lines show the current thresholds; points above {_cap} are drawn at the top edge"),
            )
        ),
        chart_selection="point",
        legend_selection=False,
    )
    _ma_d = pd.DataFrame({"symbol": _d["symbol"], "log10_baseMean": np.log10(_d["baseMean"] + 1).round(3), "log2FC": _d["log2FC"], "call": _d["call"]})
    _ma = (
        alt.Chart(_ma_d)
        .mark_circle(size=18, opacity=0.6)
        .encode(x=alt.X("log10_baseMean:Q", title="log10(average expression + 1)"), y=alt.Y("log2FC:Q", title="log2 fold change"), color=_color, tooltip=["symbol", alt.Tooltip("log2FC:Q", format=".2f")])
        .properties(width=560, height=340, title=alt.TitleParams("MA plot", subtitle="Fold change against average expression; note the wide spread at low expression"))
    )
    # Large charts go through mo.ui.altair_chart, which lifts Altair's 5,000-row limit.
    ma_plot = mo.ui.altair_chart(style(_ma), chart_selection=False, legend_selection=False)
    return ma_plot, volcano


@app.cell(hide_code=True)
def _(de, de_settings, fmt, ma_plot, markers, mo, volcano):
    _up = int((de["call"] == "Up").sum())
    _down = int((de["call"] == "Down").sum())
    _known = de[de["symbol"].isin(markers["symbol"])]
    _known_in = ", ".join(f"*{s}*" for s in sorted(_known.loc[_known["call"] != "Not selected", "symbol"])) or "none"
    mo.vstack(
        [
            mo.hstack(
                [
                    mo.stat(fmt(_up + _down), label="Genes selected", caption="pass all thresholds", bordered=True),
                    mo.stat(fmt(_up), label="Up-regulated", caption="higher in the first group", bordered=True),
                    mo.stat(fmt(_down), label="Down-regulated", caption="lower in the first group", bordered=True),
                ],
                widths="equal",
            ),
            mo.md(f"<small>Current settings: {de_settings}<br>Marker genes in the list: {_known_in}</small>"),
            mo.ui.tabs({"Volcano plot": volcano, "MA plot": ma_plot}, lazy=True),
        ]
    )
    return


@app.cell(hide_code=True)
def _(de, mo):
    _sel = de[de["call"] != "Not selected"].sort_values("padj")
    _tbl = _sel[["symbol", "log2FC", "padj", "baseMean", "description"]].assign(fold_change=lambda t: 2 ** t["log2FC"])
    de_table = mo.ui.table(
        _tbl[["symbol", "log2FC", "fold_change", "padj", "baseMean", "description"]].reset_index(drop=True),
        selection="single",
        page_size=8,
        format_mapping={"log2FC": "{:.2f}", "fold_change": "{:.2f}", "padj": "{:.2e}", "baseMean": "{:,.0f}"},
        label="Selected genes, most significant first (tick a row to plot it)",
        show_download=False,
    )
    de_table
    return (de_table,)


@app.cell(hide_code=True)
def _(de, de_find, de_table, entrez_of, gene_dotplot, mo, style, volcano):
    _clicked = (volcano.selections or {}).get("pick", {}).get("entrez") or []
    if _clicked:
        _e, _how = int(_clicked[0]), "selected in the volcano plot"
    elif len(de_table.value):
        _e, _how = int(entrez_of[de_table.value.iloc[0]["symbol"]]), "selected in the table"
    elif de_find.value:
        _e, _how = int(entrez_of[de_find.value]), "from the gene search"
    else:
        _e, _how = int(de.sort_values("padj").iloc[0]["entrez"]), "the most significant gene; click a point to choose another"
    _r = de.set_index("entrez").loc[_e]
    mo.vstack(
        [
            mo.md(f"Gene shown: *{_r['symbol']}* ({_how}). log2FC = {_r['log2FC']:.2f}, padj = {_r['padj']:.2e}, baseMean = {_r['baseMean']:,.0f}"),
            style(gene_dotplot(_e, "log2 CPM")),
        ]
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.vstack(
        [
            mo.md(
                """
    #### Questions
    1. For sleep loss in WT, raise the minimum |log2FC| from 0 to 1. How many genes remain, and
       which marker genes drop out? Should small but very consistent changes be ignored?
    2. Find a gene far to the left or right of the volcano plot but low down, and click it. Why
       is it not significant? Look at its baseMean and at the individual mice.
    3. Find a gene high on the plot but close to log2FC = 0. Is a 10% change biologically
       important? What else would you want to know?
    4. Open the MA plot. Where are the most extreme fold changes? Set the minimum baseMean
       to 100. Which genes disappear?
    5. Switch the model to "account for day". How many more sleep-loss genes are found in WT?
       Then choose "Shank3 vs WT (sleep)" and compare the two models. Do you recognise any of the
       genotype genes? *Tef*, *Per3* and *Nfil3* are worth looking up.
    6. Find *Shank3* in the "Shank3 vs WT" comparison. Is it differentially expressed? Why might
       the mutated gene itself show no change?
    """
            ),
            mo.accordion(
                {
                    "Show explanation": mo.md(
                        """
    - Statistical significance and biological importance are different things. A large fold
      change based on few reads is uncertain, while a small fold change can be highly
      significant when variability is low. Fold-change cut-offs shorten the list, but the
      cut-off is an arbitrary choice and should be reported.
    - Fold changes of lowly expressed genes are noisy, which is why such genes are filtered and
      why the MA plot fans out on the left.
    - Accounting for day removes day-to-day variation and adds about 800 sleep-loss genes. For
      genotype, the list grows from about 10 to about 20 genes and includes the circadian
      transcription factors *Tef*, *Per3* and *Nfil3*. This matches the main finding of the
      original paper ("Shank3 modulates sleep and expression of circadian transcription
      factors"), which used a more thorough batch-correction method (RUVSeq).
    - The mutation deletes one exon of *Shank3*. Gene-level counts sum the reads over all exons,
      so the total barely changes and the analysis misses the engineered change.
    - Two of the strongest genotype genes, *Acr* and *Rabl2*, lie next to *Shank3* on mouse
      chromosome 15. Genes near an engineered mutation can change because of the modification
      itself, for example DNA from another mouse strain carried along with the mutation, rather
      than because of Shank3 function. The genomic location of a DE gene is worth checking.
    """
                    )
                }
            ),
        ]
    )
    return


@app.cell(hide_code=True)
def _(bh, control_groups, fmt, logcpm, mo, np, samples, ttest_rows):
    # For each n, draw n of the 5 WT mice per condition (12 random draws), run the
    # quick t-test and count genes at BH q < 0.05.
    _wt = samples[samples["genotype"] == "WT"]
    _sd = _wt.loc[_wt["condition"] == "Sleep-deprived", "sample"].to_numpy()
    _hc = _wt.loc[_wt["condition"] == "Sleep", "sample"].to_numpy()
    _x = logcpm
    _rng = np.random.default_rng(11)
    _rows = []
    for _n in [2, 3, 4, 5]:
        for _draw in range(1 if _n == 5 else 12):
            _a = _rng.choice(_sd, _n, replace=False)
            _b = _rng.choice(_hc, _n, replace=False)
            _, _p = ttest_rows(_x[_a].to_numpy(), _x[_b].to_numpy())
            _rows.append((_n, int((bh(_p) < 0.05).sum())))
    power_curve = _rows
    reps_n = mo.ui.slider(2, 5, step=1, value=5, show_value=True, label="Mice per group")
    mo.vstack(
        [
            mo.md(
                f"""
    ### 3b. Number of replicates

    What if the experiment had used fewer mice? For each group size *n* from 2 to 5, the computer
    picks *n* of the 5 sleep-deprived WT mice and *n* of the 5 sleeping WT mice at random, runs a
    t-test on log2 CPM for every gene, and counts the genes at FDR < 0.05. A different random
    choice of mice gives a different answer, so for n = 2, 3 and 4 this is repeated 12 times, and
    each repeat is one dot in the plot. For n = 5 there is only one possible choice (all the
    mice), so there is a single dot. This t-test is simpler than DESeq2, so it finds fewer genes
    than above: {fmt(power_curve[-1][1])} sleep-loss genes with all 5 mice.
    """
            ),
            control_groups(display=[[reps_n]]),
        ]
    )
    return power_curve, reps_n


@app.cell(hide_code=True)
def _(
    alt,
    C_ACCENT,
    C_NULL,
    INK_STRONG,
    mo,
    np,
    pd,
    power_curve,
    reps_n,
    style,
):
    _d = pd.DataFrame(power_curve, columns=["mice per group", "genes found"])
    _d["chosen"] = _d["mice per group"] == reps_n.value
    # Spread the repeats sideways so dots with the same value (e.g. 0 genes) stay visible.
    _d["jitter"] = np.where(_d["mice per group"] == 5, 0.0, np.random.default_rng(3).uniform(-0.6, 0.6, len(_d)))
    _pts = (
        alt.Chart(_d)
        .mark_circle(size=70, stroke="white", strokeWidth=1)
        .encode(
            x=alt.X("mice per group:O", axis=alt.Axis(labelAngle=0)),
            xOffset=alt.XOffset("jitter:Q", scale=alt.Scale(domain=[-1, 1])),
            y=alt.Y("genes found:Q", title="Genes found (FDR < 0.05)"),
            color=alt.condition("datum.chosen", alt.value(C_ACCENT), alt.value(C_NULL)),
            tooltip=["mice per group", "genes found"],
        )
    )
    _mean = alt.Chart(_d).mark_line(color=INK_STRONG, strokeWidth=2).encode(x="mice per group:O", y="mean(genes found):Q")
    _sub = _d[_d["chosen"]]["genes found"]
    mo.vstack(
        [
            style(
                (_pts + _mean).properties(
                    width=420,
                    height=230,
                    title=alt.TitleParams("Genes found against number of mice per group", subtitle="Each dot is one random choice of mice, spread sideways so overlapping dots stay visible; the line joins the averages"),
                )
            ),
            mo.md(
                (
                    f"With all 5 mice per group there is only one possible selection: {_sub.max():,} genes. "
                    if reps_n.value == 5
                    else f"With {reps_n.value} mice per group, the 12 random choices found between {_sub.min():,} and {_sub.max():,} genes. "
                )
                + "Why do repeated selections with the same *n* give different results, and why does the number rise so steeply with *n*?"
            ),
        ]
    )
    return


@app.cell(hide_code=True)
def _(mo):
    ans_3 = mo.ui.text_area(placeholder="Section 3: significance, effect size and biological relevance; accounting for day; replicates", rows=3, full_width=True, label="Your answer")
    ans_3
    return (ans_3,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ---
    ## 4. Gene Ontology enrichment

    A list of hundreds of genes is hard to interpret directly. The **Gene Ontology (GO)** is a
    controlled vocabulary that annotates genes with the biological processes they take part in
    (for example "protein folding"), their molecular functions and their cellular locations.

    **Over-representation analysis** asks, for each GO term, whether the gene list contains more
    genes from that term than a random list of the same size would (a hypergeometric test).
    Thousands of terms are tested, so the results again need correction for multiple testing.

    The term size slider sets which GO terms are tested, according to how many background genes
    each term contains. Very small terms can look strongly enriched on the strength of two or
    three genes, and very large terms (such as "cellular process", with thousands of genes) are
    too general to be informative. A range of 15–500 genes is a common default. Narrowing the
    range also reduces the number of terms tested, and with it the correction for multiple
    testing.

    The gene list is taken from the settings in Section 3 (comparison, model, FDR, fold-change
    and expression thresholds). To change the list, change those settings and return here.
    """)
    return


@app.cell(hide_code=True)
def _(genes, go_terms, np, sparse):
    # Sparse membership matrix: GO terms (rows) x gene symbols (columns), built once.
    go_universe = genes["symbol"].to_numpy()
    _col = {s: i for i, s in enumerate(go_universe)}
    _rows, _cols = [], []
    for _i, _members in enumerate(go_terms["genes"].str.split(" ")):
        _idx = [_col[g] for g in _members if g in _col]
        _rows.extend([_i] * len(_idx))
        _cols.extend(_idx)
    go_matrix = sparse.csr_matrix((np.ones(len(_rows), dtype=np.int32), (_rows, _cols)), shape=(len(go_terms), len(go_universe)))
    go_detected = genes["expressed"].to_numpy()
    return go_detected, go_matrix, go_universe


@app.cell(hide_code=True)
def _(control_groups, mo):
    go_direction = mo.ui.radio(["Up", "Down", "Both"], value="Up", inline=True, label="Genes to test")
    go_ont = mo.ui.dropdown({"Biological process": "BP", "Molecular function": "MF", "Cellular component": "CC"}, value="Biological process", label="GO category")
    go_bg = mo.ui.radio(
        ["Genes detected in this experiment (recommended)", "All genes in the mouse annotation"],
        value="Genes detected in this experiment (recommended)",
        label="Background (the genes a random list is drawn from)",
    )
    go_size = mo.ui.range_slider(steps=[5, 10, 15, 25, 50, 100, 200, 500, 1000, 2000], value=[15, 500], show_value=True, label="GO term size (number of genes)")
    go_random = mo.ui.switch(label="Negative control: replace the list with a random list of the same size")
    go_reroll = mo.ui.button(value=0, on_click=lambda v: v + 1, label="New random list")
    go_fdr = mo.ui.slider(steps=[0.01, 0.05, 0.1, 0.25], value=0.05, show_value=True, label="Term FDR threshold")
    go_top = mo.ui.slider(5, 30, step=5, value=15, show_value=True, label="Number of terms shown")
    control_groups(
        analysis=[[go_direction, go_ont], [go_size], [go_bg], [go_random, go_reroll]],
        display=[[go_fdr, go_top]],
    )
    return (
        go_bg,
        go_direction,
        go_fdr,
        go_ont,
        go_random,
        go_reroll,
        go_size,
        go_top,
    )


@app.cell(hide_code=True)
def _(
    bh,
    de,
    go_bg,
    go_detected,
    go_direction,
    go_matrix,
    go_ont,
    go_random,
    go_reroll,
    go_size,
    go_terms,
    go_universe,
    np,
    stats,
):
    _bg = go_detected if go_bg.value.startswith("Genes detected") else np.ones(len(go_universe), dtype=bool)
    _wanted = {"Up": ["Up"], "Down": ["Down"], "Both": ["Up", "Down"]}[go_direction.value]
    _hit_symbols = set(de.loc[de["call"].isin(_wanted), "symbol"])
    _hits = np.isin(go_universe, list(_hit_symbols)) & _bg
    go_n_list = int(_hits.sum())
    if go_random.value and go_n_list:
        _rng = np.random.default_rng(99 + go_reroll.value)
        _hits = np.zeros(len(go_universe), dtype=bool)
        _hits[_rng.choice(np.flatnonzero(_bg), go_n_list, replace=False)] = True
    _K = go_matrix @ _bg.astype(np.int32)  # term size within background
    _k = go_matrix @ _hits.astype(np.int32)  # term genes in my list
    _N, _n = int(_bg.sum()), int(_hits.sum())
    _keep = (go_terms["ontology"].to_numpy() == go_ont.value) & (_K >= go_size.value[0]) & (_K <= go_size.value[1])
    _p = stats.hypergeom.sf(_k[_keep] - 1, _N, _K[_keep], _n) if _n else np.ones(_keep.sum())
    go_res = go_terms.loc[_keep, ["term_id", "name"]].copy()
    go_res["in_list"] = _k[_keep]
    go_res["term_size"] = _K[_keep]
    go_res["fold_enrichment"] = np.where(_K[_keep] > 0, (_k[_keep] / max(_n, 1)) / (_K[_keep] / _N), 0.0)
    go_res["p"] = _p
    go_res["q"] = bh(_p)
    go_res = go_res.sort_values(["q", "p"]).reset_index(drop=True)
    go_hits = _hits
    go_N = _N
    return go_N, go_hits, go_n_list, go_res


@app.cell(hide_code=True)
def _(
    alt,
    de_settings,
    fmt,
    go_direction,
    go_fdr,
    go_N,
    go_n_list,
    go_random,
    go_res,
    go_top,
    mo,
    np,
    SEQ_RANGE,
    style,
):
    _sig = go_res[go_res["q"] < go_fdr.value]
    _show = go_res.head(go_top.value).copy()
    _show["neglog10_q"] = -np.log10(_show["q"].clip(lower=1e-300))
    _show["label"] = _show["name"].str.slice(0, 55)
    _which = {"Up": "up-regulated", "Down": "down-regulated", "Both": "up- and down-regulated"}[go_direction.value]
    _random = " (replaced by a random list of the same size)" if go_random.value else ""
    _summary = mo.md(
        f"Testing {fmt(go_n_list)} {_which} genes{_random} against a background of {fmt(go_N)} genes. "
        f"{fmt(len(go_res))} GO terms tested, {fmt(len(_sig))} significant at term FDR < {go_fdr.value}.<br>"
        f"<small>Gene list settings: {de_settings}</small>"
    )
    if go_n_list < 5:
        _out = mo.callout(mo.md("The gene list has fewer than 5 genes. Loosen the thresholds in Section 3 or choose a comparison with more changes."), kind="warn")
    else:
        _chart = (
            alt.Chart(_show)
            .mark_circle(stroke="white", strokeWidth=1)
            .encode(
                y=alt.Y("label:N", sort=alt.EncodingSortField("q", order="ascending"), title=None, axis=alt.Axis(labelLimit=380)),
                x=alt.X("fold_enrichment:Q", title="Fold enrichment (observed / expected)", scale=alt.Scale(zero=True)),
                size=alt.Size("in_list:Q", title="Genes in list", scale=alt.Scale(range=[40, 400])),
                color=alt.Color("neglog10_q:Q", title="−log10(q)", scale=alt.Scale(range=SEQ_RANGE)),
                opacity=alt.condition(f"datum.q < {go_fdr.value}", alt.value(1.0), alt.value(0.35)),
                tooltip=["term_id", "name", "in_list", "term_size", alt.Tooltip("fold_enrichment:Q", format=".2f"), alt.Tooltip("q:Q", format=".2e")],
            )
            .properties(width=380, height=max(160, 22 * len(_show)), title=alt.TitleParams("Top GO terms", subtitle="Faded points are not significant at the chosen term FDR"))
        )
        _out = style(_chart)
    mo.vstack([_summary, _out])
    return


@app.cell(hide_code=True)
def _(go_res, go_top, mo):
    go_table = mo.ui.table(
        go_res.head(max(go_top.value, 30))[["term_id", "name", "in_list", "term_size", "fold_enrichment", "q"]],
        selection="single",
        page_size=8,
        format_mapping={"fold_enrichment": "{:.2f}", "q": "{:.2e}"},
        label="Tick a GO term to show its genes on the volcano plot",
        show_download=False,
    )
    go_table
    return (go_table,)


@app.cell(hide_code=True)
def _(
    C_HIGHLIGHT,
    C_NS,
    alt,
    de,
    go_hits,
    go_table,
    go_terms,
    go_universe,
    mo,
    pd,
    style,
):
    if len(go_table.value) == 0:
        _out = mo.md("*Tick a term in the table to show where its genes fall on the volcano plot.*")
    else:
        _tid = go_table.value.iloc[0]["term_id"]
        _members = set(go_terms.loc[go_terms["term_id"] == _tid, "genes"].iloc[0].split(" "))
        _list = set(go_universe[go_hits])
        _d = pd.DataFrame({"symbol": de["symbol"], "log2FC": de["log2FC"].round(3), "y": de["neglog10_padj"].clip(upper=60).round(2)})
        _d["group"] = pd.Series("other genes", index=_d.index).where(~_d["symbol"].isin(_members), "in this GO term")
        _d = _d.sort_values("group", ascending=False)
        _overlap = sorted(_members & _list)
        _chart = (
            alt.Chart(_d)
            .mark_circle()
            .encode(
                x=alt.X("log2FC:Q", title="log2 fold change"),
                y=alt.Y("y:Q", title="−log10 adjusted p-value"),
                color=alt.Color("group:N", scale=alt.Scale(domain=["in this GO term", "other genes"], range=[C_HIGHLIGHT, C_NS]), title=None),
                size=alt.condition("datum.group == 'in this GO term'", alt.value(45), alt.value(12)),
                opacity=alt.condition("datum.group == 'in this GO term'", alt.value(0.95), alt.value(0.35)),
                tooltip=["symbol", alt.Tooltip("log2FC:Q", format=".2f")],
            )
            .properties(width=480, height=300, title=alt.TitleParams(go_table.value.iloc[0]["name"], subtitle="Do the genes in this term change in the same direction?"))
        )
        _out = mo.vstack(
            [
                mo.ui.altair_chart(style(_chart), chart_selection=False, legend_selection=False),
                mo.md(f"{len(_overlap)} genes from the list are in this term: " + (", ".join(f"*{g}*" for g in _overlap[:80]) + (" …" if len(_overlap) > 80 else ""))),
            ]
        )
    _out
    return


@app.cell(hide_code=True)
def _(mo):
    mo.vstack(
        [
            mo.md(
                """
    #### Questions
    1. Before looking, which biological processes do you expect sleep loss to affect? Compare
       your prediction with the default result (sleep loss in WT, up-regulated genes, no
       fold-change cut-off).
    2. In Section 3, set the minimum |log2FC| to 1 (then 0.5) and return. How does the top of the
       list change? Which result would you report: the one from the full list or the one from the
       strongly changed genes?
    3. The term "response to glucocorticoid" (stress hormones) may appear. Is this caused by sleep
       loss or by the stress of being kept awake by handling? What control experiment could
       distinguish the two?
    4. Switch the background to all genes in the mouse annotation. How many terms are significant
       now, compared with the default? Why is this background inappropriate here?
    5. At |log2FC| ≥ 0.5 the term "female pregnancy" appears, although all the mice are male.
       What does this suggest about how GO term names should be read?
    6. Turn on the random-list control and draw several random lists. How many terms are
       significant?
    7. Similar terms often appear together, such as "protein folding" and "chaperone-mediated
       protein folding". Are these independent findings?
    8. Try the down-regulated genes and the Shank3 vs WT comparisons. What do you find?
    """
            ),
            mo.accordion(
                {
                    "Show explanation": mo.md(
                        """
    - The full list of up-regulated genes points to protein folding and the unfolded-protein
      response (chaperones such as *Hspa5*), a known effect of sleep loss in the cortex. With
      |log2FC| ≥ 1, the strongly changed genes point instead to glucocorticoid and steroid-hormone
      responses (*Sgk1*, *Fkbp5*). The thresholds used to build the list change the conclusion,
      so they should be reported with the result.
    - Gentle handling keeps the mice awake but may also raise stress-hormone levels. Enrichment
      analysis cannot separate the two; only the experimental design can, for example by using
      other methods of keeping mice awake or by measuring corticosterone.
    - The background should contain only genes that could have been detected. With all genes as
      the background, any term containing brain-expressed genes looks enriched, simply because
      those genes are detectable in brain tissue. For the default list, the number of significant
      terms rises from about 40 to about 1,100.
    - GO term names reflect the context in which genes were first described. "Female pregnancy"
      contains hormone-responsive genes, so the genes in a term are more informative than its
      name.
    - Random lists give almost no significant terms, as expected when the background is
      appropriate and the results are corrected for multiple testing.
    - GO terms are nested and overlap, so related terms share genes. Several similar terms
      usually reflect one finding.
    - Enrichment does not show whether a pathway is activated or inhibited. It is also biased
      towards well-studied genes and depends on the annotation release.
    """
                    )
                }
            ),
        ]
    )
    return


@app.cell(hide_code=True)
def _(mo):
    ans_4 = mo.ui.text_area(placeholder="Section 4: what did the enrichment analysis show, and how did the thresholds and background change it?", rows=3, full_width=True, label="Your answer")
    ans_4
    return (ans_4,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ---
    ## 5. Principal component analysis

    Each mouse is described by about 14,700 expression values, too many to inspect directly.
    **Principal component analysis (PCA)** finds new axes, the principal components (PCs), that
    capture as much of the variation between mice as possible. PC1 captures the most variation,
    PC2 the next most, and so on. Each mouse becomes a point, and mice with similar expression
    profiles lie close together.

    PCA does not use the group labels; the colours are added afterwards. If the colours separate
    along a PC, that variable is associated with a major source of variation in the data.
    """)
    return


@app.cell(hide_code=True)
def _(control_groups, mo):
    pca_transform = mo.ui.radio(["CPM (no log)", "log2 CPM"], value="log2 CPM", inline=True, label="Values")
    pca_scale = mo.ui.radio(["Centre only", "Standardise each gene"], value="Centre only", inline=True, label="Scaling")
    pca_top = mo.ui.dropdown({"100": 100, "500": 500, "1,000": 1000, "5,000": 5000, "all (~14,700)": 0}, value="500", label="Genes used (most variable)")
    pca_samples = mo.ui.dropdown(["All 20 mice", "WT mice only", "Sleeping mice only", "Sleep-deprived mice only"], value="All 20 mice", label="Mice")
    pca_color = mo.ui.dropdown(["condition", "genotype", "day", "library size"], value="condition", label="Colour by")
    pca_x = mo.ui.dropdown(["PC1", "PC2", "PC3", "PC4"], value="PC1", label="x-axis")
    pca_y = mo.ui.dropdown(["PC1", "PC2", "PC3", "PC4"], value="PC2", label="y-axis")
    pred_5 = mo.ui.radio(
        ["Genotype (WT vs Shank3)", "Condition (sleep vs sleep-deprived)", "The day of the experiment", "None of these"],
        label="Before you change anything: what will PC1 separate?",
    )
    mo.vstack(
        [
            pred_5,
            control_groups(
                analysis=[[pca_transform, pca_scale], [pca_top, pca_samples]],
                display=[[pca_color, pca_x, pca_y]],
            ),
        ]
    )
    return (
        pca_color,
        pca_samples,
        pca_scale,
        pca_top,
        pca_transform,
        pca_x,
        pca_y,
        pred_5,
    )


@app.cell(hide_code=True)
def _(np, pd):
    def run_pca(values, n_top=500, standardise=False):
        """values: genes x samples. Returns scores (samples x PCs), % variance, loadings, genes used."""
        _v = values.var(axis=1)
        _use = _v.sort_values(ascending=False).index if not n_top else _v.nlargest(n_top).index
        _use = _use[_v.loc[_use] > 0]
        _x = values.loc[_use].T.to_numpy()
        _x = _x - _x.mean(axis=0)
        if standardise:
            _x = _x / _x.std(axis=0, ddof=1)
        _u, _s, _vt = np.linalg.svd(_x, full_matrices=False)
        _k = min(10, len(_s))
        scores = pd.DataFrame(_u[:, :_k] * _s[:_k], index=values.columns, columns=[f"PC{i + 1}" for i in range(_k)])
        pct = 100 * _s**2 / (_s**2).sum()
        loadings = pd.DataFrame(_vt[:_k].T, index=_use, columns=scores.columns)
        return scores, pct[:_k], loadings

    def r2_by(scores, labels):
        """% of each PC's variance explained by a grouping (one-way ANOVA R²) or numeric variable (r²)."""
        out = {}
        for pc in scores.columns:
            y = scores[pc].to_numpy()
            if labels.dtype.kind in "if":
                out[pc] = 100 * np.corrcoef(y, labels.to_numpy())[0, 1] ** 2
            else:
                fitted = pd.Series(y).groupby(labels.to_numpy()).transform("mean").to_numpy()
                out[pc] = 100 * (1 - ((y - fitted) ** 2).sum() / ((y - y.mean()) ** 2).sum())
        return out

    return r2_by, run_pca


@app.cell(hide_code=True)
def _(
    cpm,
    logcpm,
    pca_samples,
    pca_scale,
    pca_top,
    pca_transform,
    run_pca,
    samples,
):
    _filter = {
        "All 20 mice": samples["sample"].notna(),
        "WT mice only": samples["genotype"] == "WT",
        "Sleeping mice only": samples["condition"] == "Sleep",
        "Sleep-deprived mice only": samples["condition"] == "Sleep-deprived",
    }[pca_samples.value]
    pca_meta = samples[_filter].set_index("sample")
    _vals = (logcpm if pca_transform.value == "log2 CPM" else cpm)[pca_meta.index]
    pca_scores, pca_pct, pca_loadings = run_pca(_vals, pca_top.value, pca_scale.value.startswith("Standardise"))
    return pca_loadings, pca_meta, pca_pct, pca_scores


@app.cell(hide_code=True)
def _(
    alt,
    C_ACCENT,
    C_MUT,
    C_WT,
    COND_SCALE,
    DAY_COLORS,
    INK_MUTED,
    INK_STRONG,
    mo,
    pca_color,
    pca_loadings,
    pca_meta,
    pca_pct,
    pca_scores,
    pca_x,
    pca_y,
    pd,
    r2_by,
    SEQ_RANGE,
    style,
    symbol,
):
    _d = pca_scores.join(pca_meta[["genotype", "condition", "day", "library_size"]])
    _d["library size"] = _d["library_size"] / 1e6
    _d["sample"] = _d.index
    _d["day_label"] = _d["day"].str.replace("Day ", "D")
    _xl = f"{pca_x.value} ({pca_pct[int(pca_x.value[2:]) - 1]:.0f}% of variance)"
    _yl = f"{pca_y.value} ({pca_pct[int(pca_y.value[2:]) - 1]:.0f}% of variance)"
    _color = {
        "condition": alt.Color("condition:N", scale=COND_SCALE, title="Condition"),
        "genotype": alt.Color("genotype:N", scale=alt.Scale(domain=["WT", "Shank3"], range=[C_WT, C_MUT]), title="Genotype"),
        "day": alt.Color("day:N", scale=alt.Scale(domain=[f"Day {i}" for i in range(1, 6)], range=DAY_COLORS), title="Experiment day"),
        "library size": alt.Color("library size:Q", scale=alt.Scale(range=SEQ_RANGE), title="Million reads"),
    }[pca_color.value]
    _shape = alt.Shape("genotype:N", scale=alt.Scale(domain=["WT", "Shank3"], range=["circle", "triangle-up"]), title="Genotype")
    _base = alt.Chart(_d).encode(x=alt.X(f"{pca_x.value}:Q", title=_xl), y=alt.Y(f"{pca_y.value}:Q", title=_yl))
    _pts = _base.mark_point(size=130, filled=True, opacity=0.9, stroke="white", strokeWidth=1).encode(
        color=_color, shape=_shape, tooltip=["sample", "genotype", "condition", "day", alt.Tooltip("library size:Q", format=".1f")]
    )
    _layers = [_pts]
    if pca_color.value == "day":  # secondary encoding: day written next to each dot
        _layers.append(_base.mark_text(dx=9, dy=-8, fontSize=10, color=INK_STRONG).encode(text="day_label:N"))
    _scatter = alt.layer(*_layers).properties(width=440, height=360, title=alt.TitleParams("PCA of the mice", subtitle="Each point is one mouse; triangles are Shank3 mutants"))

    _scree = (
        alt.Chart(pd.DataFrame({"PC": [f"PC{i + 1}" for i in range(len(pca_pct))], "pct": pca_pct}))
        .mark_bar(color=C_ACCENT, cornerRadiusTopLeft=4, cornerRadiusTopRight=4)
        .encode(x=alt.X("PC:N", sort=None, axis=alt.Axis(labelAngle=0), title=None), y=alt.Y("pct:Q", title="% of variance"), tooltip=["PC", alt.Tooltip("pct:Q", format=".1f")])
        .properties(width=300, height=150, title="Variance captured by each PC")
    )
    _ld = pca_loadings[pca_x.value]
    _top = _ld.abs().nlargest(10).index
    _ldf = pd.DataFrame({"gene": [symbol.get(g, g) for g in _top], "loading": _ld.loc[_top].to_numpy()})
    _load = (
        alt.Chart(_ldf)
        .mark_bar(color=INK_MUTED, cornerRadius=3)
        .encode(y=alt.Y("gene:N", sort=alt.EncodingSortField("loading", op="sum", order="descending"), title=None), x=alt.X("loading:Q", title=f"Weight on {pca_x.value}"), tooltip=["gene", alt.Tooltip("loading:Q", format=".3f")])
        .properties(width=300, height=170, title=f"Genes with the largest weights on {pca_x.value}")
    )
    # How much of each PC is explained by each variable?
    _assoc = pd.DataFrame(
        {
            "condition": r2_by(pca_scores, pca_meta["condition"]),
            "genotype": r2_by(pca_scores, pca_meta["genotype"]),
            "day": r2_by(pca_scores, pca_meta["day"]),
            "library size": r2_by(pca_scores, pca_meta["library_size"].astype(float)),
        }
    ).T[["PC1", "PC2", "PC3", "PC4"]].round(0)
    mo.vstack(
        [
            mo.hstack([style(_scatter), mo.vstack([style(_scree), style(_load)])], justify="start", gap=1, wrap=True),
            mo.md("Percentage of each PC's variation explained by each variable (100 = completely, 0 = not at all):"),
            mo.ui.table(_assoc.reset_index(names="variable"), selection=None, show_download=False, show_column_summaries=False),
        ]
    )
    return


@app.cell(hide_code=True)
def _(mo, pred_5):
    _fb = {"Condition (sleep vs sleep-deprived)": "Correct.", None: "Make your prediction above first."}.get(
        pred_5.value, "Not quite. Colour by each variable and see."
    )
    mo.vstack(
        [
            mo.md(
                f"""
    #### Questions
    1. Colour by condition, then by genotype. Which one does PC1 separate? {_fb}
    2. Plot PC2 against PC3 and colour by day. What drives these components? Use the table below
       the plot. Is this biology or a property of the experiment?
    3. Show only the sleeping mice and colour by genotype. Do WT and Shank3 mice separate on any
       PC? Would that show that the mutation has no effect? (Compare with Section 3.)
    4. Switch the values to CPM without the log. What happens to the variance explained and to
       the genes with the largest weights on PC1? Why?
    5. Compare the 100, 500 and all most variable genes, and try standardising each gene. Does
       the picture change?
    """
            ),
            mo.accordion(
                {
                    "Show explanation": mo.md(
                        """
    - With the default settings, PC1 (about 50% of the variance) separates sleep-deprived from
      sleeping mice. Sleep loss is the main source of variation in these samples.
    - PC2 and PC3 largely follow the experiment day, a batch effect. The second-largest source of
      variation is technical rather than biological. This is why the design included every group
      on every day, and why the DE model can account for day.
    - Genotype does not separate on PC1–PC3. With all 20 mice it is associated only with PC4,
      which carries about 4% of the variance, and among the sleeping mice PC1 follows the day
      rather than the genotype. PCA shows the largest sources of variation; small effects on a
      limited number of genes (Section 3) may appear only on minor components or not at all. The
      absence of a cluster does not mean the absence of an effect.
    - Without the log transformation, PC1 explains about 87% of the variance and is driven by a
      few very abundant mitochondrial transcripts (*COX1*, *CYTB*, *ND1*). The PCA then mainly
      describes those genes. The log transformation puts genes on a comparable scale.
    - The number of genes and the scaling change the picture. A PCA plot is a summary that
      depends on preprocessing, not a statistical test.
    """
                    )
                }
            ),
        ]
    )
    return


@app.cell(hide_code=True)
def _(control_groups, mo):
    trap_on = mo.ui.switch(label="Run PCA on genes selected with random groups")
    trap_reroll = mo.ui.button(value=0, on_click=lambda v: v + 1, label="New random groups")
    mo.vstack(
        [
            mo.md(
                """
    ### 5b. Selecting genes with the group labels

    A common error is to select the genes that differ most between two groups, then show a PCA
    or heatmap of those genes as evidence that the groups separate. In this demonstration we take
    the 10 sleeping mice, split them at random into groups A and B, select the 200 genes with the
    largest A-vs-B differences, and run a PCA on those genes.
    """
            ),
            control_groups(analysis=[[trap_on, trap_reroll]]),
        ]
    )
    return trap_on, trap_reroll


@app.cell(hide_code=True)
def _(
    alt,
    C_MUT,
    C_WT,
    logcpm,
    mo,
    np,
    run_pca,
    samples,
    style,
    trap_on,
    trap_reroll,
    ttest_rows,
):
    if not trap_on.value:
        _out = mo.md("*Turn on the switch to run the demonstration.*")
    else:
        _s = samples[samples["condition"] == "Sleep"]["sample"].to_numpy()
        _fake = np.random.default_rng(5 + trap_reroll.value).permutation(np.array(["A"] * 5 + ["B"] * 5))
        _x = logcpm[_s]
        _x = _x[_x.var(axis=1) > 0]
        _t, _ = ttest_rows(_x.loc[:, _fake == "A"].to_numpy(), _x.loc[:, _fake == "B"].to_numpy())
        _chosen = _x.index[np.argsort(-np.abs(_t))[:200]]
        _scores, _pct, _ = run_pca(_x.loc[_chosen], 0)
        _scores_all, _pct_all, _ = run_pca(_x, 500)

        def _plot(sc, pct, title):
            _d = sc[["PC1", "PC2"]].assign(random_group=_fake, sample=sc.index)
            return (
                alt.Chart(_d)
                .mark_point(size=130, filled=True, stroke="white", strokeWidth=1)
                .encode(
                    x=alt.X("PC1:Q", title=f"PC1 ({pct[0]:.0f}%)"),
                    y=alt.Y("PC2:Q", title=f"PC2 ({pct[1]:.0f}%)"),
                    color=alt.Color("random_group:N", scale=alt.Scale(domain=["A", "B"], range=[C_WT, C_MUT]), title="Random group"),
                    shape=alt.Shape("random_group:N", scale=alt.Scale(domain=["A", "B"], range=["circle", "square"]), title="Random group"),
                    tooltip=["sample", "random_group"],
                )
                .properties(width=300, height=260, title=title)
            )

        _out = mo.vstack(
            [
                mo.hstack(
                    [
                        style(_plot(_scores, _pct, "Genes selected using groups A and B")),
                        style(_plot(_scores_all, _pct_all, "PCA on the 500 most variable genes")),
                    ],
                    justify="start",
                    wrap=True,
                ),
                mo.callout(
                    mo.md(
                        "Groups A and B are random, yet they separate clearly in the left plot. Selecting genes using the labels and then "
                        "showing that the labels separate is circular: the result is guaranteed. The PCA on the 500 most variable genes "
                        "(right) shows no such separation. Genes used for a PCA must not be selected with the group labels that the plot "
                        "is meant to show."
                    ),
                    kind="warn",
                ),
            ]
        )
    _out
    return


@app.cell(hide_code=True)
def _(mo):
    ans_5 = mo.ui.text_area(placeholder="Section 5: what do PC1 and PC2 represent here, and why must clusters be interpreted with care?", rows=3, full_width=True, label="Your answer")
    ans_5
    return (ans_5,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ---
    ## 6. Summary and report

    Mark the safeguards that you think this analysis included, and discuss them briefly in your
    report.
    """)
    return


@app.cell(hide_code=True)
def _(mo):
    checklist = mo.ui.array(
        [
            mo.ui.switch(label="Biological replicates (several mice per group)"),
            mo.ui.switch(label="Normalisation for library size"),
            mo.ui.switch(label="Filtering of lowly expressed genes"),
            mo.ui.switch(label="Batch (experiment day) considered"),
            mo.ui.switch(label="Multiple-testing correction (FDR)"),
            mo.ui.switch(label="Effect sizes reported, not only p-values"),
            mo.ui.switch(label="Enrichment background limited to detected genes"),
            mo.ui.switch(label="PCA clusters interpreted with caution"),
        ]
    )
    student_name = mo.ui.text(placeholder="Name or student ID", label="Name / ID")
    final_report = mo.ui.text_area(
        placeholder="One paragraph: what does 5 hours of sleep loss do to gene expression in the prefrontal cortex, how confident are you, and what are the main caveats? Does Shank3 matter?",
        rows=6,
        full_width=True,
        label="Report",
    )
    mo.vstack([mo.md("Safeguards in this analysis:"), *[checklist[_i] for _i in range(8)], student_name, final_report])
    return checklist, final_report, student_name


@app.cell(hide_code=True)
def _(
    ans_1,
    ans_2,
    ans_3,
    ans_4,
    ans_5,
    checklist,
    final_report,
    mo,
    pred_1,
    pred_2,
    pred_5,
    student_name,
):
    _items = [
        "replicates",
        "normalisation",
        "low-expression filter",
        "batch/day",
        "FDR",
        "effect sizes",
        "enrichment background",
        "PCA caution",
    ]
    _ticked = [name for name, on in zip(_items, checklist.value) if on]
    _text = "\n".join(
        [
            "Omics for Neuroscience: sleep deprivation RNA-seq practical",
            f"Name / ID: {student_name.value}",
            "",
            f"Prediction (Arc after sleep loss): {pred_1.value}",
            f"Prediction (false positives): {pred_2.value}",
            f"Prediction (PC1): {pred_5.value}",
            "",
            "1. The experiment and the data:", ans_1.value, "",
            "2. Testing thousands of genes:", ans_2.value, "",
            "3. Differential expression:", ans_3.value, "",
            "4. Gene Ontology enrichment:", ans_4.value, "",
            "5. Principal component analysis:", ans_5.value, "",
            "Safeguards marked: " + (", ".join(_ticked) or "none"),
            "",
            "Report:", final_report.value,
        ]
    )
    _done = sum(bool(a.value.strip()) for a in [ans_1, ans_2, ans_3, ans_4, ans_5, final_report])
    mo.vstack(
        [
            mo.md(f"{_done} of 6 answers written."),
            mo.download(data=_text.encode("utf-8"), filename="omics_practical_answers.txt", mimetype="text/plain", label="Download answers"),
        ]
    )
    return


@app.cell(hide_code=True)
def _(data_versions, mo):
    mo.accordion(
        {
            "About the data and methods": mo.md(
                f"""
    - Data: prefrontal-cortex bulk RNA-seq from GEO
      [GSE113754](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE113754); Ingiosi AM
      et al. (2019) Shank3 modulates sleep and expression of circadian transcription factors.
      [eLife 8:e42819](https://doi.org/10.7554/eLife.42819). Reads were counted per gene by the
      authors (featureCounts, mm10).
    - Filtering: genes with more than 1 count per million in at least 5 mice (about 14,700 of
      27,000).
    - Sections 2 and 3b: Welch t-tests on log2(CPM + 1) with Benjamini–Hochberg correction,
      chosen because they are fast and simple to explain.
    - Section 3: DESeq2 (PyDESeq2 0.5.4) with the designs `~ group` and `~ day + group`,
      computed in advance. The original paper used edgeR with RUVSeq normalisation, so the
      numbers differ somewhat from the publication.
    - Section 4: hypergeometric over-representation test on GO gene sets from g:Profiler
      (annotations propagated up the GO hierarchy), with BH correction across terms.
    - Experiment day: GEO lists the mice as replicates 1–5 of an experiment repeated over 5 days.
      We assume that the replicate number is the day.

    ```
    {data_versions}
    ```
    """
            )
        }
    )
    return


if __name__ == "__main__":
    app.run()
