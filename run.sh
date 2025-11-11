#!/usr/bin/env bash
set -euo pipefail

CFG=${1:-config/config.yaml}

echo "[1/6] QC & preprocessing"
Rscript scripts/01_qc_preprocess.R --config $CFG

echo "[2/6] Integration & clustering"
Rscript scripts/02_integration_clustering.R --config $CFG

echo "[3/6] Annotation"
Rscript scripts/03_annotation.R --config $CFG

echo "[4/6] DEG & enrichment"
Rscript scripts/04_deg_enrichment.R --config $CFG

echo "[5/6] Cell–cell communication"
Rscript scripts/05_cellchat.R --config $CFG

echo "[6/6] GSVA/GSEA"
Rscript scripts/06_gsva_gsea.R --config $CFG

echo "✅ Pipeline finished. Results in ./results"
