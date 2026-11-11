# Reproducible glioma analysis code

This repository contains reproducible analysis code supporting the reported bulk RNA-seq, single-cell RNA-seq, and U251 RNA-seq figures:

- `01_bulk_RNA_analysis.R`: bulk-cohort violin plots, immune-gene heatmap, and correlation plots.
- `02_ssGSEA_GSEA.R`: methionine ssGSEA, ESTIMATE, and PPARG-stratified GSEA.
- `03_scRNA_GSE84465_plots.R`: GSE84465 cell-type, FeaturePlot, and DotPlot visualizations using the published annotations and tSNE embedding.
- `04_scRNA_GSE233304.R`: GSE233304 preprocessing, global annotation, split UMAP, and patient-level paired PPARG analysis.
- `05_U251_DEG_visualization.R`: U251 NC-versus-MR transporter heatmap and volcano plots.
- `06_PPARg_ImageJ_macro.ijm`: recursive DAPI/PPARγ image pairing, nuclear segmentation, quality-control images, and nuclear fluorescence quantification.

## Directory layout

```text
.
├── 01_bulk_RNA_analysis.R
├── 02_ssGSEA_GSEA.R
├── 03_scRNA_GSE84465_plots.R
├── 04_scRNA_GSE233304.R
├── 05_U251_DEG_visualization.R
├── 06_PPARg_ImageJ_macro.ijm
├── CITATION.cff
├── FILE_NAMING_CHANGES.md
├── LICENSE
├── RUN_ALL_CHECKLIST.md
├── gene_lists/
│   ├── methionine_metabolism_activity_gene_set.gmt
│   ├── human_transporters_detected_in_U251.txt
│   └── c5.go.v7.4.symbols.gmt             # user downloads from MSigDB
├── data/                                  # not tracked by Git by default
├── results/                               # generated; not tracked by Git
├── software_versions.txt
└── README.md
```

## Required input files

All tab-separated files must retain their original sample IDs.

| File | Structure | Used by |
|---|---|---|
| `data/TCGA_GBM_IDHWT_TPM.tsv` | First column: gene symbol; remaining columns: untransformed TPM for the 149 IDH-wildtype GBM samples | GSEA; optional ssGSEA/ESTIMATE rerun |
| `data/TCGA_methionine_scores.tsv` | Columns: `id`, `Methionine metabolism Score`; 701 TCGA tumors, 5 TCGA normal tissues, and 205 GTEx controls | TCGA/GTEx panel of S Fig. 1A |
| `data/TCGA_clinical.tsv` | Patient-level table with columns `Id`, `Grade`, and `IDH status`; `Id` is the 12-character TCGA patient ID | Automatic TCGA grade assignment |
| `data/TCGA_all_glioma_scores.tsv` | Rows: all 701 retained TCGA tumor samples; includes `id`, `Methionine metabolism Score`, `ImmuneScore`, `PPARG`, and `CD274` | S Fig. 1C, S Fig. 6C left, S Fig. 11A TCGA |
| `data/TCGA_GBM_IDHWT_scores.tsv` | Rows: the curated 149 IDH-wildtype GBM samples; includes `id`, `Methionine metabolism Score`, `PPARG`, and `CD274` | S Fig. 6C right |
| `data/TCGA_GBM_IDHWT_immunosuppression.tsv` | Same 149 sample IDs; columns `id`, `Methionine metabolism Score`, and the 40 immune-suppression genes | S Fig. 1B |
| `data/CGGA_methionine_scores.tsv` | Columns: `id`, `Methionine metabolism Score`, `Grade`, `IDH`; includes grade-known tumors and non-glioma controls | CGGA panel of S Fig. 1A |
| `data/CGGA_GBM_PPARG_CD274.tsv` | Columns `Sample`, `Histology`, `Subtype`, `PPARG`, `CD274`; 139 CGGA GBM samples | S Fig. 11A CGGA |
| `data/Gravendeel_PPARG_CD274.tsv` | Columns `Sample`, `Histology`, `Subtype`, `PPARG`, `CD274` | S Fig. 11A Gravendeel |
| `data/GSE84465_annotated.rds` | Seurat object containing the published `cell_type` metadata and original `tsne` reduction | Fig. 2E-F |
| `data/GSE233304/*.h5` | Six 10x Genomics HDF5 matrices listed below; three matched patients, each with tumor CD45+ and PBMC CD45+ libraries | GSE233304 preprocessing and annotation |
| `data/U251_TRANSPORTERS_TPM.tsv` | First column: `Gene_Name`; columns N1–N4 are NC and T1–T4 are MR, with four biological replicates per group | S Fig. 13D |
| `data/U251_all_gene_DEG.tsv` | Differential-expression results with gene symbol, log2 fold change, and raw P value | S Fig. 13B |
| `data/U251_transporter_DEG.tsv` | Rename from `volcano.txt`; transporter candidates represented in the U251 RNA-seq dataset | Fig. 7A |

Optional inputs extend S Fig. 11A to additional subsets:

| File | Structure | Used by |
|---|---|---|
| `data/CGGA_all_glioma_PPARG_CD274.tsv` | Columns `Sample`, `PPARG`, `CD274`; `Grade` or `Histology` recommended | S Fig. 11A CGGA all-glioma subset |
| `data/Gravendeel_clinical.tsv` | Sample identifiers with available `Grade`, `Histology`, and/or `IDH` fields | S Fig. 11A Gravendeel GBM subset |

`c5.go.v7.4.symbols.gmt` is distributed by MSigDB and is not bundled here. Download the exact v7.4 symbol GMT and place it in `gene_lists/`.

The GSE233304 files must be placed in `data/GSE233304/` with these names:

| GEO file | Patient | Source |
|---|---|---|
| `GSM7421861_S2_filtered_feature_bc_matrix.h5` | Patient 1 | Tumor CD45+ |
| `GSM7421862_S3_filtered_feature_bc_matrix.h5` | Patient 1 | PBMC CD45+ |
| `GSM7421864_S5_filtered_feature_bc_matrix.h5` | Patient 2 | Tumor CD45+ |
| `GSM7421865_S6_filtered_feature_bc_matrix.h5` | Patient 2 | PBMC CD45+ |
| `GSM7421869_S10_filtered_feature_bc_matrix.h5` | Patient 9 | Tumor CD45+ |
| `GSM7421870_S11_filtered_feature_bc_matrix.h5` | Patient 9 | PBMC CD45+ |

## Analysis definitions

- TPM matrices are transformed as `log2(TPM + 1)` inside the scripts.
- The methionine metabolism activity score is computed by ssGSEA using the supplied 38-gene set.
- TCGA tumor grades are assigned by matching the first 12 characters of each sample barcode to `TCGA_clinical.tsv`. Primary (`01`) and recurrent (`02`) tumor samples are retained. Samples without a resolvable grade are excluded only from grade-stratified plots and are recorded in `TCGA_grade_assignment_audit.tsv`.
- The TCGA/GTEx normal group contains the 5 TCGA sample-type-11 tissues and 205 GTEx controls. The CGGA normal group contains the 20 samples from CGGA dataset `mRNA sequencing (non-glioma as control)`.
- The immune-suppression heatmap uses the curated 149 IDH-wildtype GBM samples and includes every immune gene present in its input table. The left strip displays Pearson R only when nominal `P < 0.05`; non-significant cells are left blank. The heatmap itself uses per-gene row z-scores.
- PPARG GSEA divides samples at the median of log2-transformed PPARG expression. The ranking metric is limma's moderated t statistic for `PPARG-high minus PPARG-low`. PPARG itself is removed from the ranked vector.
- GSEA positive NES indicates enrichment in PPARG-high samples; negative NES indicates enrichment in PPARG-low samples.
- GSE84465 plots use the cell annotations and original tSNE coordinates stored in the supplied Seurat object. The script does not alter cell-type labels or recalculate dimensionality reduction, clustering, marker genes, or annotations.
- Fig. 2E combines the cell-type DimPlot with FeaturePlots for `PPARG`, `FARP1`, and `S100A16`. Fig. 2F shows a PPARG DotPlot grouped by `cell_type`.
- GSE233304 contains three matched patients and six libraries. The script uses `SCTransform`, PCA, Harmony integration, shared-nearest-neighbor clustering, and UMAP before assigning cluster labels from predefined immune marker modules. Related macrophage modules are merged into `TAMs`, while classical monocytes are reported as `Mono/Mac`.
- PPARG comparisons use each patient as the biological replicate. Expression is summarized within patient and cell group before paired two-sided t-tests. The cell-level violin plot is descriptive and carries no cell-level inferential test.
- U251 volcano thresholds follow the figure legend: raw `P < 0.05` and `|log2FC| >= 1`. The plotting script derives Down/Not significant/Up from these thresholds instead of trusting the source `group` column, which only records fold-change direction.
- The transporter heatmap uses `log2(TPM + 1)` followed by per-gene row z-scoring.
- The human-transporter reference list contains 576 unique candidate genes. Of these, 253 have corresponding quantified entries in the U251 RNA-seq dataset and are therefore represented in the heatmap and transporter-focused volcano plot. The number 253 describes the dataset intersection, not the size of the human-transporter reference list.
- The ImageJ macro treats CH1 as DAPI and CH2 as PPARγ for legacy paired TIFF files. For exported channel files, `c2x` is DAPI and `c1x` is PPARγ; merged `c0-3x` files are ignored.
- Nuclear ROIs are obtained from an 8-bit smoothed DAPI copy using an adjusted Otsu threshold, morphological closing, hole filling, and predefined area/circularity filters. Raw DAPI and PPARγ intensities are not altered by segmentation processing.
- PPARγ background is estimated from the 10th percentile of the image histogram. The macro reports background-corrected nuclear mean fluorescence, nuclear integrated fluorescence, perinuclear mean fluorescence, and the nuclear/perinuclear ratio. Edge nuclei without a complete perinuclear ring are excluded from fluorescence summaries.

## Data preparation

- All tab-separated inputs use UTF-8 encoding and retain the original sample identifiers.
- Expression matrices contain gene symbols in the first column and samples in subsequent columns.
- TPM inputs remain untransformed; scripts apply `log2(TPM + 1)` where required.
- Cohort-level score tables use the exact column names shown above.

## Running the analyses

For a complete input, execution, output, and error-reporting checklist, see `RUN_ALL_CHECKLIST.md`.

From the repository root:

```bash
Rscript 01_bulk_RNA_analysis.R
Rscript 02_ssGSEA_GSEA.R ssgsea
Rscript 02_ssGSEA_GSEA.R estimate
Rscript 02_ssGSEA_GSEA.R gsea
Rscript 03_scRNA_GSE84465_plots.R
Rscript 04_scRNA_GSE233304.R . preprocess
Rscript 05_U251_DEG_visualization.R
```

For GSE233304, the first argument is the repository root and the second is either `preprocess` or `reuse`. `preprocess` rebuilds the Seurat object from the six H5 files. After a successful first run, `reuse` starts from `results/04_scRNA_GSE233304/objects/GSE233304_preprocessed.rds` and reruns annotation, plots, and paired statistics.

Run `06_PPARg_ImageJ_macro.ijm` from Fiji using **Plugins → Macros → Run**. Select the parent image directory when prompted. The macro searches that directory and up to four subdirectory levels, then writes results to a new `PPARg_Nuclear_Analysis` folder. Inspect every `DAPI_MASK_CHECK` and `NUCLEAR_ROI_CHECK` image before accepting the quantitative output.

Each script writes to its own subdirectory under `results/` and records `sessionInfo.txt` there.

## Validation status

All six workflows were successfully executed with the study inputs. The R analyses were validated under R 4.5.0 on Windows 10, and the image-analysis macro was validated with ImageJ 1.54f and Bio-Formats 7.0.0. Exact package versions are recorded in `software_versions.txt`. The scripts include input validation, sample-assignment audits, analysis summaries, and `sessionInfo()` output. The complete human-transporter reference list contains 576 candidates; the visualized subset consists of the 253 candidates with quantified entries in the U251 RNA-seq dataset.

## Citation

This repository supports the manuscript *Methionine Availability Drives m6A-Dependent Regulation of PPARG in Macrophages to Promote Immunosuppression in Glioblastoma*, currently under revision and not yet accepted for publication. Hong Hu (胡鸿) is the repository submitter. Citation metadata are provided in `CITATION.cff`; the formal journal citation and DOI will be added after acceptance.

## License

The analysis code is released under the MIT License. See `LICENSE` for details.

