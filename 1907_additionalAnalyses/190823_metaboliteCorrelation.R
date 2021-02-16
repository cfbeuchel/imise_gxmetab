#' ---
#' title: "Correlation of metabolites (and subsequent clustering of gene association)"
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
#'         smooth_scroll: FALSE
#' ---
#+ include=F
#==========# 
# Initiate #
#==========#

#+ set.LibPath, include=F
# define alternative package directory
r_on_server <- FALSE
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
  if(grepl(x = getwd(), pattern =  "mnt")){
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
  "here",
  "ggplot2",
  "corrplot",
  "magrittr",
  "ComplexHeatmap"
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
source(here("../functions/beta_plot.R"))
source(here("../functions/all_annotation_tables.R"))
an <- all_annotation_tables(mount.point = "/net/ifs2/san_projekte/projekte/genstat/")
option_setup()

#+ include=F
#==========#
#   DONE   #
#==========#

#' # Load metabolite data

#' The data was cleaned in STEP1
load("res/dat_for_STEP2.RData", verbose = TRUE)

# all metabolites
am <- an$metab.annot$metabolite
an$metab.annot[,pheno.class.q:=pheno.class]
an$metab.annot[grepl("^Q",metabolite),pheno.class.q:="quotient"]

#' # Show total correlation of metabolities  in each study

#' ## single study metabolite correlation

cor.a1 <- metab.a1[,..am] %>% na.omit() %>% cor(method = "pearson")
cor.a1[upper.tri(cor.a1, diag = FALSE)] %>% hist(breaks=50, main="Pearson correlation of LIFE Adult Metabolites")
diag(cor.a1) <- 0
cor.a1[lower.tri(cor.a1)] <- 0

#' ## Joint study correlation

# cam = correlation all metabolites
cam <- rbindlist(
  list(
    metab.a1[,..am],
    metab.b3.heart[,..am],
    metab.b3.ami[,..am],
    metab.sorb[,..am]
  )
) %>% cor(method = "pearson", use = "pairwise")