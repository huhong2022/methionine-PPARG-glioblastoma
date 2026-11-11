# Complete execution checklist

Run the analyses from a fresh R session and from the repository root. Do not change the input column names documented in `README.md`.

## 1. Prepare input files

Place cohort tables and Seurat objects directly in `data/`. Place the six GSE233304 H5 files in `data/GSE233304/`. Place the MSigDB file `c5.go.v7.4.symbols.gmt` in `gene_lists/`.

Before running, confirm that these required files exist:

```text
data/TCGA_methionine_scores.tsv
data/TCGA_clinical.tsv
data/TCGA_all_glioma_scores.tsv
data/TCGA_GBM_IDHWT_scores.tsv
data/TCGA_GBM_IDHWT_immunosuppression.tsv
data/TCGA_GBM_IDHWT_TPM.tsv
data/CGGA_methionine_scores.tsv
data/CGGA_GBM_PPARG_CD274.tsv
data/Gravendeel_PPARG_CD274.tsv
data/GSE84465_annotated.rds
data/U251_TRANSPORTERS_TPM.tsv
data/U251_all_gene_DEG.tsv
data/U251_transporter_DEG.tsv
gene_lists/c5.go.v7.4.symbols.gmt
```

The six required GSE233304 files are listed in the main `README.md`.

## 2. Run the R analyses

From the repository root:

```bash
Rscript 01_bulk_RNA_analysis.R .
Rscript 02_ssGSEA_GSEA.R all .
Rscript 03_scRNA_GSE84465_plots.R .
Rscript 04_scRNA_GSE233304.R . preprocess
Rscript 05_U251_DEG_visualization.R .
```

The first GSE233304 run must use `preprocess`. After it completes successfully, annotation and plots can be regenerated without repeating preprocessing:

```bash
Rscript 04_scRNA_GSE233304.R . reuse
```

Each R script writes to its corresponding folder under `results/`. Retain all generated `sessionInfo.txt` files.

## 3. Run the Fiji macro

Open Fiji and select **Plugins → Macros → Run**, then choose `06_PPARg_ImageJ_macro.ijm`. Select the parent folder containing the paired DAPI and PPARγ TIFF files.

The macro writes to `PPARg_Nuclear_Analysis/`. Review every file ending in:

```text
_DAPI_MASK_CHECK.tif
_DAPI_SELECTED_ROIS_MASK.tif
_NUCLEAR_ROI_CHECK.tif
```

Do not accept fluorescence measurements from images with visibly merged nuclei, extensive missed nuclei, or background incorrectly classified as nuclei.

## 4. Check the expected results

| Workflow | Minimum result to inspect |
|---|---|
| `01_bulk_RNA_analysis.R` | Cohort violin plots, immune heatmap, and correlation plots |
| `02_ssGSEA_GSEA.R` | ssGSEA scores, ESTIMATE scores, both target GSEA plots, and GSEA result tables |
| `03_scRNA_GSE84465_plots.R` | Cell-type tSNE, three FeaturePlots, and PPARG DotPlot |
| `04_scRNA_GSE233304.R` | Split UMAP, final annotation table, PPARG patient summary, paired-test table, and paired plot |
| `05_U251_DEG_visualization.R` | Whole-transcriptome volcano plot, transporter volcano plot, and transporter heatmap |
| `06_PPARg_ImageJ_macro.ijm` | Pairing report, processing status, per-cell table, image summary, masks, and ROI overlays |

## 5. Report a failed run

For any failure, record:

1. Script name and command used.
2. Complete error message, starting from the first warning or error.
3. `sessionInfo()` output for R errors, or Fiji version for macro errors.
4. Input filename and its column names; do not send the full large matrix unless requested.
5. For plotting problems, include the generated plot or a screenshot.
