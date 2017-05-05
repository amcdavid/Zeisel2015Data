
 
#' Zeisel2015
#' Single cell RNAseq data from Zeisel2015
#' @docType package
#' @aliases Zeisel2015Data-package
#' @name Zeisel2015Data
#' @description  UMI-deduplicated mouse cortex data, and ERCC-spike ins from Zeisel, et al 2015 (DOI: 10.1126/science.aaa1934).
#' @details Wild type (CD-1) and transgenic mice (5HT3EGFP on a CD-1 background)
#' between postnatal 21 and 31 days of both sexes were used. Somatosensory cortex or CA1 hippocampal
#' brain regions were dissociated into a single cell suspension, ERCC and UMI added and cDNA libraries generated on Fluidigm C1. Cells were selected to be valid for analysis if they passed manual inspection of images and had more than 2500 total detected RNA molecules,  This resulted in 3315 valid cells. About 310 cells were excluded from further analysis as they clustered separately from everything else and did not express any distinct markers. Based on their expression profiles, these cells were probably neurons of low quality, resulting in 3005 remaining cells.
#' Data was downloaded from http://linnarssonlab.org/cortex/ on 01-Feb-2017.
#' @import data.table
#' @import Biostrings
#' @import stringr
#' @seealso \link{zeisel2015}
NULL
NULL
