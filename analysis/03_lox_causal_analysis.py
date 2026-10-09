
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
from scipy.stats import spearmanr
from statsmodels.stats.multitest import multipletests
import statsmodels.api as sm
from IPython.display import display

BASE = Path("/content/liver_fibrosis")

expr = pd.read_csv(
    BASE / "log_expression.csv", index_col=0
)
pca = pd.read_csv(
    BASE / "pca_results.csv", index_col=0
)

# Match sample titles
common = [s for s in pca.index if s in expr.columns]
expr = expr[common].astype(float)
pca = pca.loc[common].copy()

# Ordered fibrosis stages
stage_map = {"F0-2": 0, "F3": 1, "F4": 2}
pca["stage_order"] = pca["disease group"].map(stage_map)

if pca["stage_order"].isna().any():
    raise ValueError("Unexpected fibrosis-stage labels.")

# 1. PC2-associated genes
pc2 = pca["PC2"].astype(float)

pc2_results = []

for gene, values in expr.iterrows():
    rho, pval = spearmanr(values.values, pc2.values)
    pc2_results.append((gene, rho, pval))

pc2_results = pd.DataFrame(
    pc2_results,
    columns=["gene", "PC2_rho", "PC2_p"]
).dropna()

pc2_results["PC2_FDR"] = multipletests(
    pc2_results["PC2_p"],
    method="fdr_bh"
)[1]

pc2_results = pc2_results.sort_values("PC2_FDR")

# 2. Etiology-adjusted fibrosis association
# Linear model: expression ~ stage + disease etiology
# Stage is coded as an ordered exploratory trend.
design = pd.get_dummies(
    pca[["stage_order", "disease"]],
    columns=["disease"],
    drop_first=True,
    dtype=float
)

design = sm.add_constant(design).astype(float)

Y = expr.T.values
X = design.values

# Matrix-based ordinary least squares
beta = np.linalg.lstsq(X, Y, rcond=None)[0]
residual = Y - X @ beta

n, k = X.shape
df = n - k

if df <= 0:
    raise ValueError("Insufficient residual degrees of freedom.")

sigma2 = (residual**2).sum(axis=0) / df
cov_diag = np.diag(np.linalg.pinv(X.T @ X))

stage_idx = list(design.columns).index("stage_order")
stage_beta = beta[stage_idx]

stage_se = np.sqrt(
    sigma2 * cov_diag[stage_idx]
)

from scipy.stats import t as student_t

tstat = np.divide(
    stage_beta,
    stage_se,
    out=np.full_like(stage_beta, np.nan),
    where=stage_se > 0
)

pvals = 2 * student_t.sf(np.abs(tstat), df)

association = pd.DataFrame({
    "gene": expr.index,
    "stage_beta": stage_beta,
    "stage_p": pvals
}).dropna()

association["stage_FDR"] = multipletests(
    association["stage_p"],
    method="fdr_bh"
)[1]

# 3. Combine results
combined = association.merge(
    pc2_results,
    on="gene",
    how="inner"
)

combined = combined.sort_values("stage_FDR")

combined.to_csv(
    BASE / "fibrosis_gene_associations.csv",
    index=False
)

print("TOP FIBROSIS-ASSOCIATED GENES")
display(combined.head(20))

print("\nSignificant stage-associated genes (FDR < 0.05):")
print((combined["stage_FDR"] < 0.05).sum())

# 4. Plot key ECM-related genes
targets = [
    "COL1A1", "COL3A1", "ACTA2",
    "TGFB1", "LOX", "ITGB1"
]

present = [g for g in targets if g in expr.index]

fig, axes = plt.subplots(
    2, 3, figsize=(15, 9)
)

for ax, gene in zip(axes.flat, present):
    values = expr.loc[gene]

    groups = [
        values[pca["disease group"] == stage].values
        for stage in ["F0-2", "F3", "F4"]
    ]

    ax.boxplot(groups, tick_labels=["F0-2", "F3", "F4"])
    ax.set_title(gene)
    ax.set_ylabel("log2(expression + 1)")

for ax in axes.flat[len(present):]:
    ax.axis("off")

plt.suptitle("ECM and Mechanosignaling Candidate Genes")
plt.tight_layout()
plt.show()

print("\nECM candidate statistics:")
display(
    combined[
        combined["gene"].isin(targets)
    ].sort_values("stage_FDR")
)

print("\nSaved:", BASE / "fibrosis_gene_associations.csv")
