#!/bin/bash
# Generates minimal test inputs for interpret/pathway_enrichment.
# Run from the repository root:
#   bash src/interpret/pathway_enrichment/test_data/test_data_script.sh

set -eo pipefail

REPO_ROOT=$(git rev-parse --show-toplevel)
cd "$REPO_ROOT"

OUT="src/interpret/pathway_enrichment/test_data"

python3 - <<'PYEOF'
import os, sys
import numpy as np
import pandas as pd
import anndata as ad
import mudata as mu

out = "src/interpret/pathway_enrichment/test_data"
os.makedirs(out, exist_ok=True)
rng = np.random.default_rng(42)

# ── synthetic h5mu (100 cells × 200 genes) ──────────────────────────────────
n_cells, n_genes = 100, 200
gene_names = [f"GENE{i:04d}" for i in range(n_genes)]
counts = rng.negative_binomial(5, 0.3, size=(n_cells, n_genes)).astype("float32")
adata = ad.AnnData(X=counts)
adata.var_names = gene_names
adata.obs_names = [f"cell_{i}" for i in range(n_cells)]
adata.obs["subpopulation"] = pd.Categorical(
    ["ExN1"] * 50 + ["ExN2"] * 50
)
mdata = mu.MuData({"rna": adata})
mdata.write_h5mu(f"{out}/input.h5mu")
print(f"Wrote {out}/input.h5mu  ({n_cells} cells × {n_genes} genes)")

# ── synthetic DESeq2-like CSV ────────────────────────────────────────────────
lfc = rng.normal(0, 2, size=n_genes)
pval_raw = rng.uniform(0, 1, size=n_genes)
padj = np.clip(pval_raw * n_genes / np.arange(1, n_genes + 1)[np.argsort(pval_raw).argsort()], 0, 1)
de = pd.DataFrame({
    "log2FoldChange": lfc,
    "pvalue": pval_raw,
    "padj": padj,
    "baseMean": rng.uniform(10, 1000, size=n_genes),
    "lfcSE": rng.uniform(0.1, 0.5, size=n_genes),
    "stat": lfc / rng.uniform(0.1, 0.5, size=n_genes),
}, index=gene_names)
de.index.name = "gene"
de.to_csv(f"{out}/deseq2_results.csv")
print(f"Wrote {out}/deseq2_results.csv  ({n_genes} genes, "
      f"{(padj < 0.05).sum()} significant at padj<0.05)")

# ── tiny GMT file (3 gene sets, entirely from the synthetic gene list) ───────
gmt_path = f"{out}/test_gene_sets.gmt"
with open(gmt_path, "w") as fh:
    fh.write("PATHWAY_A\tna\t" + "\t".join(gene_names[:30]) + "\n")
    fh.write("PATHWAY_B\tna\t" + "\t".join(gene_names[30:70]) + "\n")
    fh.write("PATHWAY_C\tna\t" + "\t".join(gene_names[70:120]) + "\n")
print(f"Wrote {gmt_path}  (3 gene sets)")
PYEOF

echo ""
echo "Test data generated. Example viash run commands:"
echo ""
echo "# prerank GSEA with local GMT:"
echo "viash run src/interpret/pathway_enrichment/config.vsh.yaml -- \\"
echo "  --input       $OUT/input.h5mu \\"
echo "  --input_degenes $OUT/deseq2_results.csv \\"
echo "  --method      prerank \\"
echo "  --gene_sets   $OUT/test_gene_sets.gmt \\"
echo "  --permutation_num 100 \\"
echo "  --output      $OUT/output_prerank.h5mu \\"
echo "  --output_csv_dir $OUT/results_prerank/"
echo ""
echo "# ORA with local GMT:"
echo "viash run src/interpret/pathway_enrichment/config.vsh.yaml -- \\"
echo "  --input       $OUT/input.h5mu \\"
echo "  --input_degenes $OUT/deseq2_results.csv \\"
echo "  --method      ora \\"
echo "  --gene_sets   $OUT/test_gene_sets.gmt \\"
echo "  --output      $OUT/output_ora.h5mu \\"
echo "  --output_csv_dir $OUT/results_ora/"
