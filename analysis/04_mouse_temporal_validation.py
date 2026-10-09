
# ============================================================
# FINAL ANALYSIS: INDEPENDENT LIVER FIBROSIS TIME-COURSE
# GEO: GSE222576 | Mouse CCl4 fibrosis | Automatic download
# ============================================================

!pip -q install GEOparse mygene statsmodels scipy seaborn

import os, re, warnings
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import GEOparse
import mygene

from pathlib import Path
from scipy.stats import linregress, ttest_ind
from statsmodels.stats.multitest import multipletests

warnings.filterwarnings("ignore", category=FutureWarning)
sns.set_theme(style="whitegrid")

OUT = Path("/content/liver_fibrosis/temporal_validation")
OUT.mkdir(parents=True, exist_ok=True)

GENES = [
    "Lox", "Loxl1", "Loxl2", "Loxl3",
    "Acta2", "Ccn2", "Thbs1", "Nedd9",
    "Col1a1", "Col3a1", "Tgfb1", "Itgb1",
    "Adamts1", "Timp1", "Mmp2", "Mmp9",
    "Fn1", "Pdgfrb", "Fos", "Egr1"
]

MODULES = {
    "Fibrogenic activation":
        ["Acta2", "Ccn2", "Pdgfrb", "Tgfb1"],
    "Collagen / ECM":
        ["Col1a1", "Col3a1", "Fn1"],
    "Crosslinking":
        ["Lox", "Loxl1", "Loxl2", "Loxl3"],
    "Matrix remodeling":
        ["Thbs1", "Adamts1", "Timp1", "Mmp2", "Mmp9"],
    "Adhesion / response":
        ["Itgb1", "Nedd9", "Fos", "Egr1"]
}

print("=" * 65)
print("STEP 1: AUTOMATICALLY DOWNLOAD GEO")
print("=" * 65)

geo = GEOparse.get_GEO(
    geo="GSE222576",
    destdir=str(OUT),
    annotate_gpl=False,
    silent=True
)

print("Samples downloaded:", len(geo.gsms))

# ------------------------------------------------------------
# 2. SAMPLE METADATA + NORMALIZED EXPRESSION
# ------------------------------------------------------------

records = []
columns = {}

for gsm_id, gsm in geo.gsms.items():

    title = gsm.metadata.get("title", [""])[0]
    m = re.search(r"(\d+)\s*weeks?", title, re.I)

    if not m:
        raise ValueError(f"Cannot identify time point: {title}")

    week = int(m.group(1))
    control = "control" in title.lower()

    condition = "Control" if control else "CCl4"
    label = "Control" if control else f"W{week}"

    table = gsm.table.copy()

    if not {"ID_REF", "VALUE"}.issubset(table.columns):
        raise ValueError(f"Missing expression data for {gsm_id}")

    values = pd.to_numeric(
        table["VALUE"], errors="coerce"
    )

    s = pd.Series(
        values.to_numpy(),
        index=table["ID_REF"].astype(str).to_numpy()
    ).groupby(level=0).mean()

    columns[gsm_id] = s

    records.append({
        "sample": gsm_id,
        "title": title,
        "week": week,
        "condition": condition,
        "group": label
    })

meta = pd.DataFrame(records).set_index("sample")
expression = pd.DataFrame(columns)

print("\nSample groups:")
print(meta.groupby(["condition", "week"]).size())

print("\nExpression matrix:", expression.shape)

# These values were RMA processed by the original study.
# DO NOT apply log2 again.

if len(meta) != 18:
    raise ValueError("Expected 18 samples: inspect GEO metadata.")

if set(meta.loc[meta.condition == "CCl4", "week"]) != {2,4,6,8,10}:
    raise ValueError("Unexpected treatment time points.")

# ------------------------------------------------------------
# 3. AUTOMATIC MOUSE GENE ID MAPPING
# ------------------------------------------------------------

print("\nSTEP 2: Mapping mouse gene symbols to Entrez IDs")

mg = mygene.MyGeneInfo()

mapped = mg.querymany(
    GENES,
    scopes="symbol",
    fields="symbol,entrezgene",
    species="mouse",
    as_dataframe=False,
    verbose=False
)

gene_to_probe = {}

for item in mapped:
    if item.get("notfound"):
        continue

    query = str(item.get("query", ""))
    symbol = str(item.get("symbol", ""))

    if query.lower() != symbol.lower():
        continue

    entrez = item.get("entrezgene", item.get("_id"))
    if entrez is None:
        continue

    # GEO uses custom Entrez-based CDF probe-set IDs.
    probe = f"{int(entrez)}_at"

    if probe in expression.index:
        gene_to_probe[query] = probe

print("Genes successfully mapped:", len(gene_to_probe))
print("Missing:", sorted(set(GENES) - set(gene_to_probe)))

if len(gene_to_probe) < 5:
    raise ValueError(
        "Too few probes mapped. Check GPL22598 CDF annotation."
    )

gene_expr = pd.DataFrame({
    gene: expression.loc[probe]
    for gene, probe in gene_to_probe.items()
}).T

gene_expr.columns.name = "sample"

gene_expr.to_csv(OUT / "mouse_gene_expression.csv")
meta.to_csv(OUT / "mouse_sample_metadata.csv")

# ------------------------------------------------------------
# 4. TIME-DEPENDENT GENE ANALYSIS
# ------------------------------------------------------------

print("\nSTEP 3: Testing temporal associations")

treated = meta[meta.condition == "CCl4"].copy()
control = meta[meta.condition == "Control"].copy()

results = []

for gene in gene_expr.index:

    x = treated["week"].to_numpy(float)
    y = gene_expr.loc[gene, treated.index].to_numpy(float)

    valid = np.isfinite(y)

    if valid.sum() < 6:
        continue

    fit = linregress(x[valid], y[valid])

    # Independent comparison:
    # Controls collected at week 10 versus treated at week 10.
    y_control = gene_expr.loc[
        gene, control.index
    ].dropna().to_numpy(float)

    w10 = treated.index[treated.week == 10]

    y_10 = gene_expr.loc[
        gene, w10
    ].dropna().to_numpy(float)

    if len(y_control) >= 2 and len(y_10) >= 2:
        contrast_p = ttest_ind(
            y_10, y_control,
            equal_var=False
        ).pvalue

        contrast_effect = y_10.mean() - y_control.mean()
    else:
        contrast_p = np.nan
        contrast_effect = np.nan

    results.append({
        "gene": gene,
        "slope_per_week": fit.slope,
        "trend_p": fit.pvalue,
        "trend_r": fit.rvalue,
        "W10_minus_control": contrast_effect,
        "W10_vs_control_p": contrast_p
    })

results = pd.DataFrame(results)

for pcol, qcol in [
    ("trend_p", "trend_FDR"),
    ("W10_vs_control_p", "W10_vs_control_FDR")
]:
    results[qcol] = np.nan
    valid = results[pcol].notna()

    if valid.any():
        results.loc[valid, qcol] = multipletests(
            results.loc[valid, pcol],
            method="fdr_bh"
        )[1]

results = results.sort_values("trend_FDR")
results.to_csv(OUT / "mouse_temporal_statistics.csv", index=False)

print("\nTEMPORAL ANALYSIS:")
print(results.round(5).to_string(index=False))

# ------------------------------------------------------------
# 5. PLOT EXPRESSION ACROSS ACTUAL EXPERIMENTAL WEEKS
# ------------------------------------------------------------

print("\nSTEP 4: Generating time-course plots")

focus = [
    g for g in [
        "Lox", "Acta2", "Ccn2",
        "Col1a1", "Col3a1", "Thbs1",
        "Nedd9", "Timp1"
    ]
    if g in gene_expr.index
]

fig, axes = plt.subplots(
    2, 4, figsize=(17, 9), sharex=True
)

for ax, gene in zip(axes.flat, focus):

    values = gene_expr.loc[gene]

    # Treated sample means + sample-level observations
    for week in [2, 4, 6, 8, 10]:
        ids = treated.index[treated.week == week]
        yy = values.loc[ids].to_numpy(float)
        ax.scatter(
            np.repeat(week, len(yy)),
            yy, alpha=0.50, s=30, color="tab:blue"
        )

    means = [
        values.loc[
            treated.index[treated.week == week]
        ].mean()
        for week in [2, 4, 6, 8, 10]
    ]

    ax.plot(
        [2, 4, 6, 8, 10],
        means,
        "-o",
        color="tab:blue",
        linewidth=2,
        label="CCl4"
    )

    # The control group was also collected at week 10.
    control_mean = values.loc[control.index].mean()

    ax.axhline(
        control_mean,
        linestyle="--",
        color="gray",
        label="10w control mean"
    )

    ax.set_title(gene, fontweight="bold")
    ax.set_xlabel("Weeks of CCl4 treatment")
    ax.set_ylabel("RMA expression")

for ax in axes.flat[len(focus):]:
    ax.axis("off")

axes.flat[0].legend(fontsize=8)
fig.suptitle(
    "GSE222576 | Liver fibrosis transcriptional time course",
    fontsize=15
)
plt.tight_layout()
plt.savefig(
    OUT / "gene_time_courses.png",
    dpi=200, bbox_inches="tight"
)
plt.show()

# ------------------------------------------------------------
# 6. GENE MODULES: EXPLORATORY TRANSCRIPTIONAL SCORES
# ------------------------------------------------------------

print("\nSTEP 5: Computing gene module scores")

module_rows = []

for module, members in MODULES.items():

    present = [g for g in members if g in gene_expr.index]

    if len(present) < 2:
        continue

    X = gene_expr.loc[present].astype(float)

    # Standardize each gene across all 18 samples.
    mu = X.mean(axis=1)
    sd = X.std(axis=1).replace(0, np.nan)

    Z = X.sub(mu, axis=0).div(sd, axis=0)
    scores = Z.mean(axis=0)

    for sample, score in scores.items():
        module_rows.append({
            "sample": sample,
            "module": module,
            "score": score,
            "week": int(meta.loc[sample, "week"]),
            "condition": meta.loc[sample, "condition"],
            "group": meta.loc[sample, "group"]
        })

module_df = pd.DataFrame(module_rows)
module_df.to_csv(OUT / "module_scores.csv", index=False)

if not module_df.empty:

    plt.figure(figsize=(11, 6))

    order = ["Control", "W2", "W4", "W6", "W8", "W10"]

    sns.lineplot(
        data=module_df[module_df.condition == "CCl4"],
        x="week",
        y="score",
        hue="module",
        estimator="mean",
        errorbar=None,
        marker="o",
        linewidth=2
    )

    plt.title("Exploratory ECM / activation module trajectories")
    plt.xlabel("Weeks of CCl4 exposure")
    plt.ylabel("Mean standardized expression score")
    plt.tight_layout()
    plt.savefig(
        OUT / "module_trajectories.png",
        dpi=200, bbox_inches="tight"
    )
    plt.show()

# ------------------------------------------------------------
# 7. OPTIONAL COMPARISON WITH OUR HUMAN DATASET
# ------------------------------------------------------------

print("\nSTEP 6: Human-mouse directional comparison")

human_path = Path(
    "/content/liver_fibrosis/fibrosis_gene_associations.csv"
)

if human_path.exists():

    human = pd.read_csv(human_path)
    human["gene_key"] = human["gene"].str.upper()
    mouse = results.copy()
    mouse["gene_key"] = mouse["gene"].str.upper()

    comparison = mouse.merge(
        human[["gene_key", "stage_beta", "stage_FDR"]],
        on="gene_key",
        how="inner"
    )

    comparison["same_direction"] = (
        np.sign(comparison["slope_per_week"])
        == np.sign(comparison["stage_beta"])
    )

    comparison.to_csv(
        OUT / "human_mouse_comparison.csv",
        index=False
    )

    print(comparison[
        ["gene", "stage_beta", "slope_per_week",
         "stage_FDR", "trend_FDR", "same_direction"]
    ].to_string(index=False))

else:
    print("Human file not found: comparison skipped.")

# ------------------------------------------------------------
# 8. FINAL OUTPUT
# ------------------------------------------------------------

print("\n" + "=" * 65)
print("ANALYSIS COMPLETE")
print("=" * 65)

print("Saved results in:", OUT)
print("\nOutput files:")

for f in sorted(OUT.iterdir()):
    if f.is_file() and f.suffix in [".csv", ".png"]:
        print(" -", f.name)

print("""
INTERPRETATION RULES:
1. A positive slope means expression increases with
   treatment duration in the sampled mouse groups.
2. FDR is corrected across the selected genes, not
   the entire mouse transcriptome.
3. Mouse weeks are exposure durations, not repeated
   measurements of the same individual.
4. Controls were collected at week 10, NOT week 0.
5. Transcript expression is NOT a measurement of
   LOX enzyme activity, crosslinks or ECM stiffness.
6. Temporal association alone does NOT prove causality.
""")
