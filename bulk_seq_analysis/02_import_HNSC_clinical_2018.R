# Import & inspect TCGA-HNSC clinical data (cBioPortal, PanCancer Atlas 2018)

clinical <- read.delim(
  "GDC_HNSC/cBioPortal/hnsc_tcga_pan_can_atlas_2018_clinical_data.tsv",
  check.names = FALSE
)

## Dimensions
dim(clinical)               # patients x clinical fields
nrow(clinical); ncol(clinical)

## Column names
colnames(clinical)

## First look at the data
head(clinical, 5)
str(clinical[1:10])

## Key identifiers
head(clinical$`Patient ID`, 10)
head(clinical$`Sample ID`, 10)

## Quick distributions of a few key clinical fields
table(clinical$Sex, useNA = "always")
table(clinical$`Sample Type`, useNA = "always")
table(clinical$Subtype, useNA = "always")           # includes HPV+/HPV- status
table(clinical$`American Joint Committee on Cancer Tumor Stage Code`, useNA = "always")
summary(clinical$`Diagnosis Age`)

## Missingness overview across all columns
na_counts <- sapply(clinical, function(x) sum(is.na(x) | x == "" | x == "[Not Available]"))
sort(na_counts, decreasing = TRUE)

## ---------------------------------------------------------------------
## Anatomic subtype breakdown (oral cavity / oropharynx / larynx / etc.)
## cBioPortal's "Tumor Disease Anatomic Site" is uninformative here (always
## "Head and Neck"), so use the GDC clinical file's tissue_or_organ_of_origin,
## which records the specific primary site per patient.
clinical_gdc <- read.delim("GDC_HNSC/HNSC_clinical_GDC.tsv", check.names = FALSE)

subsite_map <- c(
  "Tongue, NOS" = "Oral Cavity",
  "Floor of mouth, NOS" = "Oral Cavity",
  "Mouth, NOS" = "Oral Cavity",
  "Cheek mucosa" = "Oral Cavity",
  "Gum, NOS" = "Oral Cavity",
  "Hard palate" = "Oral Cavity",
  "Lip, NOS" = "Oral Cavity",
  "Anterior floor of mouth" = "Oral Cavity",
  "Lower gum" = "Oral Cavity",
  "Retromolar area" = "Oral Cavity",
  "Upper Gum" = "Oral Cavity",
  "Ventral surface of tongue, NOS" = "Oral Cavity",
  "Palate, NOS" = "Oral Cavity",
  "Border of tongue" = "Oral Cavity",
  "Tonsil, NOS" = "Oropharynx",
  "Base of tongue, NOS" = "Oropharynx",
  "Oropharynx, NOS" = "Oropharynx",
  "Posterior wall of oropharynx" = "Oropharynx",
  "Larynx, NOS" = "Larynx",
  "Supraglottis" = "Larynx",
  "Hypopharynx, NOS" = "Hypopharynx",
  "Overlapping lesion of lip, oral cavity and pharynx" = "Overlapping/Other",
  "Mandible" = "Overlapping/Other"
)

clinical_gdc$hnc_subtype <- subsite_map[clinical_gdc$tissue_or_organ_of_origin]
table(clinical_gdc$hnc_subtype, useNA = "always")
nrow(clinical_gdc)          # total patients in current GDC clinical pull

## ---------------------------------------------------------------------
## Restrict to the 520 primary-tumor samples actually used in the bulk
## RNA-seq analysis (HNSC_RNAseq_counts.tsv also carries 44 matched normal
## and 2 metastatic samples, which inflate the clinical table above).
counts_header <- scan(
  "GDC_HNSC/HNSC_RNAseq_counts.tsv",
  what = character(), nlines = 1, sep = "\t", quiet = TRUE
)
sample_barcodes <- counts_header[-1]                       # drop "gene_id"
sample_type_code <- substr(sample_barcodes, 14, 15)
table(sample_type_code)                                    # 01 tumor, 06 metastatic, 11 normal

tumor_barcodes <- sample_barcodes[sample_type_code == "01"]
tumor_patients <- substr(tumor_barcodes, 1, 12)             # e.g. "TCGA-CQ-6221"
length(unique(tumor_patients))                              # should be 520

clinical_tumor <- clinical_gdc[clinical_gdc$submitter_id %in% tumor_patients, ]
nrow(clinical_tumor)                                        # 520 patients

table(clinical_tumor$hnc_subtype, useNA = "always")
