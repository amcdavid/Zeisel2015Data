library(data.table)
library(Biostrings)
library(stringr)
library(SummarizedExperiment)
library(PenultimateGEXContainer)
data_dir = "../inst/extdata/"
endo = fread(file.path(data_dir, "expression_mRNA_17-Aug-2014.txt"))
ercc = fread(file.path(data_dir, "expression_spikes_17-Aug-2014.txt"))

parse_linnearson <- function(dt, nc, nr) {
    #one blank line separator and starting index is inclusive so offset by 2
  dt_expr = as.matrix(dt[(2+nr):nrow(dt),(1+nc):ncol(dt), with=F])
  dt_expr_num = as.numeric(dt_expr)
  stopifnot(!any(is.na(dt_expr_num)))
  dim(dt_expr_num) = dim(dt_expr)
  dt_cdata = as.data.table(t(dt[seq_len(nr), (nc+1):ncol(dt), with=F]))
  setnames(dt_cdata, unlist(dt[seq_len(nr), nc, with=F]))
  dt_fdata = dt[(2+nr):nrow(dt), seq_len(nc), with=F]
  list(expr=dt_expr_num, cdat=dt_cdata, fdat=dt_fdata)
}

endo_list = parse_linnearson(endo, 2, 9)
ercc_list = parse_linnearson(ercc, 2, 9)
## Mix 1 used.
ercc_conc = fread(file.path(data_dir, "ercc_concentrations.txt"))
ercc_conc[,':='(`Re-sort ID` = NULL,
                subgroup = NULL,
                `concentration in Mix 2 (attomoles/ul)` = NULL,
                `expected fold-change ratio` = NULL,
                `log2(Mix 1/Mix 2)` = NULL)]
ercc_seq = DataFrame(fread(file.path(data_dir, 'ercc_seq.txt')))
ercc_seq$Sequence = DNAStringSet(ercc_seq$Sequence)
ercc_seq = within(ercc_seq, {
        transcript_len_bp = nchar(Sequence)
        gc = alphabetFrequency(Sequence)
        AfreqCenter = gc[,1]/transcript_len_bp - .25
        CfreqCenter = gc[,2]/transcript_len_bp - .25
        GfreqCenter = gc[,3]/transcript_len_bp - .25
        TfreqCenter = gc[,4]/transcript_len_bp - .25
        gc = GenBank =X5prime_assay = X3prime_assay = NULL
        `ERCC ID` = ERCC_ID
        ERCC_ID = NULL
        Sequence=NULL
})
ercc_conc = merge(ercc_conc, ercc_seq)
ercc_conc$V1 = ercc_conc$`ERCC ID`
ercc_conc = merge(ercc_conc, ercc_list$fdat)
fdat_all = rbind(endo_list$fdat, as.data.frame(ercc_conc), fill=TRUE)
fdat_all$primerid = make.unique(fdat_all$V1)
setnames(fdat_all, c('concentration.in.Mix.1..attomoles.ul.'), c('attomoles_per_ul'))
expr_all = log2(rbind(endo_list$expr, ercc_list$expr)+1)
endo_list$cdat$cdr_endo = colSums(endo_list$expr>0)
endo_list$cdat$cdr_ercc = colSums(ercc_list$expr>0)
endo_list$cdat$`mRNA molecules` = as.numeric(endo_list$cdat$`total mRNA mol`)
endo_list$cdat$pgex_cell_size_um = as.numeric(endo_list$cdat$`diameter`)
endo_list$cdat$pgex_filter_study_pass = TRUE
endo_list$cdat$pgex_biosample = endo_list$cdat$pgex_batch = str_split_fixed(endo_list$cdat$cell_id, "_", 2)[,1]
rownames(expr_all) = fdat_all$primerid
se = SummarizedExperiment(assay=expr_all, rowData=fdat_all, colData=endo_list$cdat)

pgex_exp = list(pgex_platform = 'fluidigm_c1',
                pgex_chemistry='SMARTer',
                pgex_has_umi = TRUE,
                pgex_ercc_version='ERCC RNA Spike-In Mix 1',
                pgex_molecule = 'polyA RNA',
                pgex_geo_accession='GSE60361',
                pgex_pmid_accession=25700174L,
                pgex_is_public=TRUE,
                pgex_transformation='shifted log2')
data = pgex_to_list(PGEX_promote(PGEXContainerProto(se, pgex_experiment = pgex_exp)))
batch_tab <- with(endo_list$cdat, table(pgex_batch, interaction(age, sex)))
rowSums(batch_tab>0)


keepDataObjects('data')
 
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

#' Data from Zeisel 2015
#'@name zeisel2015
#'@docType data
#'@format a \code{list} containing the following fields
#'\describe{
#'\item{exprs}{expression data, rows are genes, columns cells}
#'\item{cdata}{cell covariates}
#'\item{fdata}{feature covariates}
#'}
#' The attributes
#' \describe{
#'\item{has_umi}{Are these data UMI-deduplicated?}
#'\item{has_ercc}{Does the genes contain ERCC spike ins?}
#'\item{trans}{What transformation was applied to the (expected) counts?}
#'}
#' are set accordingly.
#' The 
#'@source Zeisel 2015  (DOI: 10.1126/science.aaa1934)

