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
  "ComplexHeatmap",
  "circlize"
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
source(here("../functions/plot_metabolite_correlation.R"))
source(here("../functions/all_annotation_tables.R"))
an <- all_annotation_tables(mount.point = "/net/ifs2/san_projekte/projekte/genstat/")
option_setup()

#+ include=F
#==========#
#   DONE   #
#==========#

#' # Load data

#' The data was cleaned in STEP1
load("res/dat_for_STEP2.RData", verbose = TRUE)

# all metabolites
am <- an$metab.annot$metabolite
an$metab.annot[,pheno.class.q:=pheno.class]
an$metab.annot[grepl("^Q",metabolite),pheno.class.q:="quotient"]

# load supplementary table 1
st1 <- fread("obj/191209_metabOverview.csv")

#' # Show total correlation of metabolities  in each study

#' Ideas for correct displaying: 
#' - seperate clusters visually
#' - 

#' ## single study metabolite correlation

#' ### Adult

# calculate correlation of all metabolites
cor.a1 <- metab.nc.a1[,..am] %>% 
  cor(method = "pearson", use = "pairwise")
cor.a1[upper.tri(cor.a1, diag = FALSE)] %>%
  hist(
    breaks=50, 
    main="Pearsons correlation of LIFE Adult Metabolites"
  )

# use plotting function
plot_metabolite_correlation(
  mat = cor.a1, 
  st1 = st1
  )

#' ### HEART AMI

# calculate correlation of all metabolites
cor.b3.ami <- metab.nc.b3.ami[,..am] %>% 
  cor(method = "pearson", use = "pairwise")
cor.b3.ami[upper.tri(cor.a1, diag = FALSE)] %>%
  hist(
    breaks=50, 
    main="Pearsons correlation of LIFE HEART AMI Metabolites"
  )

# use plotting function
plot_metabolite_correlation(
  mat = cor.b3.ami, 
  st1 = st1
)

#' ### HEART

# calculate correlation of all metabolites
cor.b3.heart <- metab.nc.b3.heart[,..am] %>% 
  cor(method = "pearson", use = "pairwise")
cor.b3.heart[upper.tri(cor.b3.heart, diag = FALSE)] %>%
  hist(
    breaks=50, 
    main="Pearsons correlation of LIFE Heart Metabolites"
  )

# use plotting function
plot_metabolite_correlation(
  mat = cor.b3.heart, 
  st1 = st1
)

#' ### Sorb

# calculate correlation of all metabolites
cor.sorb <- metab.nc.sorb[,..am] %>% 
  cor(method = "pearson", use = "pairwise")
cor.sorb[upper.tri(cor.sorb, diag = FALSE)] %>%
  hist(
    breaks=50, 
    main="Pearsons correlation of Sorb study Metabolites"
  )

# use plotting function
plot_metabolite_correlation(
  mat = cor.sorb, 
  st1 = st1
)

#' ## Joint study correlation

# cam = correlation all metabolites
cam <- rbindlist(
  list(
    metab.nc.a1[,..am],
    metab.nc.b3.heart[,..am],
    metab.nc.b3.ami[,..am],
    metab.nc.sorb[,..am]
  )
) %>% cor(method = "pearson", use = "pairwise")

#save for use in other script
save(list = "cam", file = "obj/metabCorrelation.RData")

cam[upper.tri(cam, diag = FALSE)] %>%
  hist(
    breaks=50, 
    main="Pearsons correlation of joint study Metabolites"
  )

#' This will be dominated by the larger studies:

tiff(dpx("all_metaboliteCorrelation.tiff", "plt/"), 
     width = 10, 
     height = 10, 
     res=300, 
     unit="in")

# correlation matrix
plot_metabolite_correlation(
  mat = cam, 
  st1 = st1
)

dev.off()

#' # Session Info

devtools::session_info()
