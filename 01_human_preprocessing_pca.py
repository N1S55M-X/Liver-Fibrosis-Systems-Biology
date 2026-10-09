
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from IPython.display import display

# =============================================
# 1. LOAD DATA
# =============================================

BASE = Path("/content/liver_fibrosis")

expression = pd.read_csv(
    "/content/GSE276114_raw.count.txt",
    sep="\t",
    index_col=0
)

metadata = pd.read_csv(
    BASE / "sample_metadata.csv"
)

print("Expression matrix:", expression.shape)
print("Metadata:", metadata.shape)

# =============================================
# 2. MATCH BY SAMPLE TITLE
# =============================================

expression.columns = expression.columns.astype(str).str.strip()
metadata["title"] = metadata["title"].astype(str).str.strip()

if metadata["title"].duplicated().any():
    raise ValueError("Duplicate sample titles in metadata.")

if expression.columns.duplicated().any():
    raise ValueError("Duplicate sample names in expression matrix.")

meta = metadata.set_index("title")

common = [
    name for name in expression.columns
    if name in meta.index
]

print("Matched samples:", len(common))

if len(common) < 3:
    raise ValueError("Sample titles still do not match.")

expression = expression[common].apply(
    pd.to_numeric, errors="raise"
)

meta = meta.loc[common].copy()

assert list(expression.columns) == list(meta.index)

print("\nFibrosis distribution:")
print(meta["disease group"].value_counts())

print("\nDisease distribution:")
print(meta["disease"].value_counts())

# =============================================
# 3. QUALITY CONTROL
# =============================================

expression = expression.loc[
    ~expression.index.duplicated()
]

if expression.isna().any().any():
    raise ValueError("Missing expression values detected.")

if (expression < 0).any().any():
    raise ValueError("Negative expression values detected.")

print("\nExpression summary:")
print(expression.stack().describe())

# =============================================
# 4. FILTER + EXPLORATORY NORMALIZATION
# =============================================

keep = (expression > 1).sum(axis=1) >= max(
    3, int(0.1 * expression.shape[1])
)

filtered = expression.loc[keep]

# Exploratory transformation only:
# Does not assume integer counts or re-normalize
# potentially preprocessed values using library sizes.
log_expression = np.log2(filtered + 1)

print("\nGenes retained:", len(filtered))

# =============================================
# 5. PCA
# =============================================

X = log_expression.T

top_genes = X.var(axis=0).nlargest(
    min(2000, X.shape[1])
).index

X_scaled = StandardScaler().fit_transform(
    X[top_genes]
)

pca = PCA(n_components=2)
pcs = pca.fit_transform(X_scaled)

results = pd.DataFrame(
    pcs,
    columns=["PC1", "PC2"],
    index=common
).join(meta)

# =============================================
# 6. VISUALIZE
# =============================================

fig, axes = plt.subplots(1, 2, figsize=(15, 6))

for stage in ["F0-2", "F3", "F4"]:
    subset = results[
        results["disease group"] == stage
    ]

    axes[0].scatter(
        subset["PC1"],
        subset["PC2"],
        label=f"{stage} (n={len(subset)})",
        s=45,
        alpha=0.75
    )

axes[0].set_title("Fibrosis Severity")

for disease in sorted(results["disease"].dropna().unique()):
    subset = results[
        results["disease"] == disease
    ]

    axes[1].scatter(
        subset["PC1"],
        subset["PC2"],
        label=f"{disease} (n={len(subset)})",
        s=45,
        alpha=0.75
    )

axes[1].set_title("Disease Etiology")

for ax in axes:
    ax.set_xlabel(
        f"PC1 ({pca.explained_variance_ratio_[0]*100:.1f}%)"
    )
    ax.set_ylabel(
        f"PC2 ({pca.explained_variance_ratio_[1]*100:.1f}%)"
    )
    ax.legend()
    ax.grid(alpha=0.2)

plt.tight_layout()
plt.show()

# =============================================
# 7. SAVE RESULTS
# =============================================

results.to_csv(BASE / "pca_results.csv")
log_expression.to_csv(
    BASE / "log_expression.csv"
)

print("\nANALYSIS COMPLETE")
print("Samples analyzed:", len(common))
print("Genes analyzed:", len(filtered))

print("\nVariance explained:")
print(pca.explained_variance_ratio_)

print("\nSaved results in:", BASE)
