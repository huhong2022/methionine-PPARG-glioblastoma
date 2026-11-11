# File naming changes

This list records source names that were replaced by stable, descriptive repository names. Numeric suffixes added by repeated downloads, local project paths, and ambiguous names such as `volcano.txt` are not used in the repository.

| Source or local name | Repository name | Note |
|---|---|---|
| `GSE233304_6paired_core_final_annotation(1).R` | `04_scRNA_GSE233304.R` | Preprocessing and annotation section merged into the public script |
| `单细胞患者样本配对检验(1).R` | `04_scRNA_GSE233304.R` | Patient-level PPARG paired-analysis section merged into the same script |
| `Tight_DAPI_Nuclear_Analysis_CH2_Red_NewExport_Recursive_SizeFix_v4(4).ijm` | `06_PPARg_ImageJ_macro.ijm` | Version history and local-project wording removed; CH2 is identified as PPARγ |
| `TCGA-GBM-IDH-WT.txt` (gene-by-sample TPM matrix) | `data/TCGA_GBM_IDHWT_TPM.tsv` | Untransformed TPM input |
| `TCGA-GBM-IDH-WT.txt` (sample-level score table) | `data/TCGA_GBM_IDHWT_scores.tsv` | Curated 149-sample IDH-wildtype GBM table |
| `IGC cor heatmap.txt` | `data/TCGA_GBM_IDHWT_immunosuppression.tsv` | Methionine score and immunosuppression-gene table |
| `TCGA+GTEx.txt` | `data/TCGA_methionine_scores.tsv` | TCGA/GTEx score table used for the cohort plot |
| `TPM_symbol(N_T)707.txt` | Input represented by the documented TCGA TPM tables | The title was inaccurate; the matrix contains 706 samples |
| `TRANSPORTERS.txt` | `data/U251_TRANSPORTERS_TPM.tsv` | N1-N4 are NC; T1-T4 are MR |
| `total_diff.txt` | `data/U251_all_gene_DEG.tsv` | Whole-transcriptome differential-expression results |
| `volcano.txt` | `data/U251_transporter_DEG.tsv` | Transporter-focused differential-expression results |
| Annotated GSE84465 Seurat object | `data/GSE84465_annotated.rds` | Contains the original annotation and tSNE coordinates |

The six GSE233304 GEO H5 names are retained unchanged because their GSM and sample identifiers are required for traceability.
