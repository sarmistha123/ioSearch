% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Omic1_rmNoOmic2.R
\name{Omic1_rmNoOmic2}
\alias{Omic1_rmNoOmic2}
\title{Create omic1 data matrix for both disease groups by collecting information from all subject-specific omic1 files.}
\usage{
Omic1_rmNoOmic2(
  samplesheet,
  data.type1,
  data.type2,
  sample.type,
  sampledata_folder,
  clinical_samplesheet,
  pattern,
  omic1.expr.column,
  samplesheet.sample.type.column,
  samplesheet.data.category.column,
  samplesheet.case.id.column,
  samplesheet.data.type.column,
  clin.samplesheet.case.id.column,
  clin.samplesheet.stage.column,
  clin.samplesheet.age.column
)
}
\arguments{
\item{samplesheet}{Samplesheet matrix downloaded from TCGA or constructed with TCGAbiolinks (as shown in the example).}

\item{data.type1}{Omic1 data type name (as in TCGA).}

\item{data.type2}{Omic2 data type name (as in TCGA).}

\item{sample.type}{Sample type name (Tumor type such as "Primary Tumor", "Solid Tumor" etc.).}

\item{sampledata_folder}{Omics data folder downloaded from TCGA or using TCGAbiolinks (as shown in the example).}

\item{clinical_samplesheet}{Clinical data matrix downloaded from TCGA or constructed with TCGAbiolinks (as shown in the example).}

\item{pattern}{String to extract omic1 data files.}

\item{omic1.expr.column}{Column number containing the omic1 expression values in the omic-specific datafile}

\item{samplesheet.sample.type.column}{Column number containing the sample type in the samplesheet file.}

\item{samplesheet.data.category.column}{Column number containing the data category in the samplesheet file.}

\item{samplesheet.case.id.column}{Column number containing the case IDs in the samplesheet file.}

\item{samplesheet.data.type.column}{Column number containing the data type in the samplesheet file.}

\item{clin.samplesheet.case.id.column}{Column number containing the case IDs in the clinical file.}

\item{clin.samplesheet.stage.column}{Column number containing the tumor stage in the clinical file.}

\item{clin.samplesheet.age.column}{Column number containing the age values in the clinical file.}
}
\value{
Individual matrices for phenotype, subject IDs in each disease group, omic1 expression data for each disease group.
}
\description{
Create omic1 data matrix for both disease groups by collecting information from all subject-specific omic1 files.
}
