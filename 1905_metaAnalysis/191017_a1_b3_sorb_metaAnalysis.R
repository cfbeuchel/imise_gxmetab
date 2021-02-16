#' ---
#' title: "Meta-Analysis of Gx-Metabolite Association - Data Preparation"
#' author: "Carl Beuchel"
#' date: "`r format(Sys.time(), '%d %B, %Y')`"
#' output:
#'     html_document:
#'      fig_captions: yes
#'      theme: spacelab #sandstone #flatfly #spacelab
#'      highlight: pygments
#'      toc: TRUE
#'      toc_depth: 3
#'      code_folding: hide
#'      number_sections: TRUE
#'      toc_float:
#'          smooth_scroll: FALSE
#' ---
#+ include=F
#==========# 
# Initiate #
#==========#

#+ set.LibPath, include=F
# define alternative package directory
r_on_server <- TRUE
if (r_on_server == TRUE) {
  bp <- "/net/ifs1/san_projekte/projekte/genstat/"
  computer <- "forostar"
  .libPaths(
    paste0(
      bp,
      "07_programme/rpackages/",
      computer
    )
  )
  # set number of cores
  nc <- 10
} else {
  if(grepl(x= getwd(),pattern =  "mnt")){
    bp <- "/mnt/mount1/genstat"
  } else {
    bp <- "/net/ifs1/san_projekte/projekte/genstat/"
  }
  # set number of cores
  nc <- parallel::detectCores()-1
}

#+ set.proj.root, include=F
# set project root
root.dir <- paste0(bp, "02_projekte/1703_ge_metab_a1_b3_sorbs/1805_gxMetaboliteAssociation/")
setwd(root.dir)

#+ load.packages, include=F
for (i in c(
  "data.table",
  "CarlHelpR",
  "toolboxH",
  "R.utils",
  "here"
)){
  suppressPackageStartupMessages(
    library(i, character.only=TRUE
    ))}

#+ script.setup, include=F
library("knitr", lib.loc = .libPaths()[1])
library("evaluate", lib.loc = .libPaths()[1])
opts_chunk$set(
  cache = F,
  results = "markup",
  echo = T,
  include = T,
  message = T,
  warning = T
)
# set root directory not to file path
opts_knit$set(root.dir=root.dir)
source(here("../functions/option_setup.R"))
source(here("../functions/meta_analysis.R"))
option_setup()

#+ include=F
#==========#
#   DONE   #
#==========#

#' # Data

#' ## Adjusted Limma results with covariates in Model

# LIMMA
nf <- newest_file("a1_b3_sorb_limmaAssocWithCovarsResults", subfolder = "res", print_full = T)
res.covar.limma <- fread(nf)

# limma
setnames(res.covar.limma, 
         old = c("ilmn", "n_individuals", "logFC", "beta.se", "P.Value"),
         new = c("markerID", "n", "beta", "se", "p")
           )

# create dummy place holders, since I do not need them
res.covar.limma[, `:=`(
  codedAll = "A",
  noncodedAll = "B",
  eaf = 0.2,
  imputed = FALSE,
  infoscore = 0.9
)]

#=================================================================#

#' # Meta Analysis

#' ## Pre-adjusted metabolites with covariates in LIMMA

#' Loop through each metabolite

# calculate the meta analysis
res.list.meta <- lapply(unique(res.covar.limma$var), function(i){
  res.list.meta <- meta_analysis(
    phenotype = i,
    assocData = res.covar.limma)
  return(res.list.meta)
})

# consolidate
res.meta <- rbindlist(res.list.meta)

# remove columns that are of no use
colfilt <- grepl(pattern = "nWeightedMAF|minMAF|nWeightedInfoScore|minInfoScore|eaf.|maf.|imputed.|infoscore.",
     x = names(res.meta))
res.meta <- res.meta[,.SD,.SDcols=names(res.meta)[!colfilt]]

# loop through each cohort and the fixed and random effects result
cohort.cols <- paste0(".", unique(res.covar.limma$cohort))
for(i in c(cohort.cols, "FEM", "REM")){
  
  pcol <- paste0("p", i)
  bcol <- paste0("beta", i)
  secol <- paste0("se", i)
  
  # calculate the confidence Intervals/hierarchical fdr
  set(x = res.meta, j = paste0("ci.l", i), value = res.meta[[bcol]] - 1.96 * res.meta[[secol]])
  set(x = res.meta, j = paste0("ci.u", i), value = res.meta[[bcol]] + 1.96 * res.meta[[secol]])
  
  # add hierarchical fdr
  hfdr <- addHierarchFDR(pvalues = res.meta[!is.na(get(pcol))][[pcol]],
                         categs = res.meta[!is.na(get(pcol)), metabolite],
                         fdrmethod_level1 = "BH",
                         fdrmethod_level2 = "BH",
                         correctionLevel1 = "BB")
  
  # filter for NA in pcol T/F
  na.filt <- which(!is.na(res.meta[[pcol]]))
  set(x = res.meta, i = na.filt, j = paste0("hfdr.1", i), value = hfdr$fdr_level1)
  set(x = res.meta, i = na.filt, j = paste0("hfdr.2", i), value = hfdr$fdr_level2)
  set(x = res.meta, i = na.filt, j = paste0("hfdr.sig", i), value = hfdr$hierarch_fdr5proz)
  
}

# save results
save_csv_carl(file = res.meta,
              file_name = "a1_b3_sorb_metaAnalysisWithCovarsLimmaResults",
              subfolder = "res")

#' # Session Overview
devtools::session_info()
