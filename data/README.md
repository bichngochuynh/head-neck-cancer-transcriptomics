# Data

No raw data is committed to this repository. Both datasets are public and are
either downloaded automatically by the scripts or fetched once by hand.

| Dataset | Accession / source | How to obtain | Used by |
|---|---|---|---|
| TCGA-HNSC bulk RNA-seq (STAR counts, primary tumours) | [GDC TCGA-HNSC](https://portal.gdc.cancer.gov/projects/TCGA-HNSC) | Downloaded automatically by `bulk_seq_analysis/01_download_HNSC_TCGAbiolinks.R` (`TCGAbiolinks`) into `GDC_HNSC/` | `bulk_seq_analysis/` |
| TCGA-HNSC clinical (HPV status, OS, PFS) | [cBioPortal – hnsc_tcga_pan_can_atlas_2018](https://www.cbioportal.org/study/summary?id=hnsc_tcga_pan_can_atlas_2018) | Download `hnsc_tcga_pan_can_atlas_2018_clinical_data.tsv` and place it in `GDC_HNSC/cBioPortal/` | `bulk_seq_analysis/02_...R` |
| HNSCC scRNA-seq | [GEO GSE181919](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE181919) | Download `GSE181919_UMI_counts.txt.gz` and `GSE181919_Barcode_metadata.txt.gz` from the GEO supplementary files into the project root | `scRNA_analysis/01_...R` |

Both expression datasets are already processed count matrices (gene x
sample / gene x cell), so no FASTQ alignment step is required to reproduce
the analysis.
