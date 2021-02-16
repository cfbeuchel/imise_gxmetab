#' ---
#' title: "Association of Gene Expression with Metabolite Data"
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
  computer <- "amanMRO"
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
  "magrittr",
  "here",
  "readxl",
  "ggplot2",
  "extrafont",
  "lumi",
  "limma"
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
source(here("../functions/beta_plot_cohort.R"))
source(here("../functions/option_setup.R"))
option_setup()
source(paste0(bp, "07_programme/rtools/1807_gx_tools/AllgFunktionenExpressionspipeline_1906.R"))
source(paste0(bp, "07_programme/rtools/1409_pathway_enrichment/hypergeometricEnrichment_KEGGneu_GO_Reactome_DOSE_viaEntrez_180608.R"))
source("../functions/each_cohort_limma.R")

#+ include=F
#==========#
#   DONE   #
#==========#

#' # Load

#' ## Gx probe and sample Annotations

#' 190617: I re-did the preprocessing of all data, removing diabetes status and log-BMI from the adjusted covariates. Therefore I need to repeat the whole analysis. Additionally, I want to seperate data preparation and analysis. 

#' ### LIFE Adult

sample.adult <- fread(paste0(bp, "02_projekte/1705_ge_lifea1v2/180719_a1_gePreProForMetab/tosend_NoBMI_NoDiabetes/s401_1_sampleannot_HT12v4.txt"))
probe.adult  <- fread(paste0(bp, "02_projekte/1705_ge_lifea1v2/180719_a1_gePreProForMetab/tosend_NoBMI_NoDiabetes/s401_1_probeannot_HT12v4.txt"))

#' ### LIFE Heart

sample.heart.old <- fread(paste0(bp, "/02_projekte/1703_ge_metab_a1_b3_sorbs/1805_gxMetaboliteAssociation/dat/s12_sampleannot_HT12v4.txt"))
sample.heart <- fread(paste0(bp, "02_projekte/1102_ge_lifeb3/07_prepro_ami_allb3_for_metab/tosend_NoBMI_NoDiabetes/sampleannot_HT12v4.txt"))
probe.heart  <- fread(paste0(bp, "02_projekte/1102_ge_lifeb3/07_prepro_ami_allb3_for_metab/tosend_NoBMI_NoDiabetes/probeannot_HT12v4.txt"))

#' ### Sorb

sample.sorb <- fread(paste0(bp, "02_projekte/1706_ge_sorben_v2/181009_sorb_gePreProForMetab/tosend_NoBMI_NoDiabetes/s401_1_sampleannot_HT12v4.txt"))
probe.sorb  <- fread(paste0(bp, "02_projekte/1706_ge_sorben_v2/181009_sorb_gePreProForMetab/tosend_NoBMI_NoDiabetes/s401_1_probeannot_HT12v4.txt"))

#' ## No Covariate adjusted gx probe annotation

#' ### LIFE Adult

sample.adult.nc <- fread(paste0(bp, "02_projekte/1705_ge_lifea1v2/180719_a1_gePreProForMetab/tosend_noCovariates/s401_1_sampleannot_HT12v4.txt"))
probe.adult.nc  <- fread(paste0(bp, "02_projekte/1705_ge_lifea1v2/180719_a1_gePreProForMetab/tosend_noCovariates/s401_1_probeannot_HT12v4.txt"))

#' ### LIFE Heart

sample.heart.nc <- fread(paste0(bp, "02_projekte/1102_ge_lifeb3/07_prepro_ami_allb3_for_metab/tosend_noCovariates/sampleannot_HT12v4.txt"))
probe.heart.nc  <- fread(paste0(bp, "02_projekte/1102_ge_lifeb3/07_prepro_ami_allb3_for_metab/tosend_noCovariates/probeannot_HT12v4.txt"))

#' ### Sorb

sample.sorb.nc <- fread(paste0(bp, "02_projekte/1706_ge_sorben_v2/181009_sorb_gePreProForMetab/tosend_noCovariates/s401_1_sampleannot_HT12v4.txt"))
probe.sorb.nc  <- fread(paste0(bp, "02_projekte/1706_ge_sorben_v2/181009_sorb_gePreProForMetab/tosend_noCovariates/s401_1_probeannot_HT12v4.txt"))

#' ## Metabolite Annotation

# load the functino that loads all annotation tables
source(paste0(bp, "02_projekte/1703_ge_metab_a1_b3_sorbs/functions/all_annotation_tables.R"))
annot <- all_annotation_tables(mount.point = bp)

#' ## Gx Data and Metabolite Data

#' this object was saved in step 1 and contains all formatted objects for LIMMA analysis

what.dat <- load("res/dat_for_STEP2.RData")
what.dat
load("res/dat_for_STEP2.RData")

#' ## Load Covars

#'  load the covariate names to be adjusted to
directory <- "/net/ifs1/san_projekte/projekte/genstat/02_projekte/1703_ge_metab_a1_b3_sorbs/170829_ConfounderAnalysis/190207_ExtraAnalysis/results"
nf <- newest_file(look_for = "FullModelWithoutQuotientsInclHighMissings", directory = directory, print_full = T)
adjust <- fread(nf)
cm <- adjust[covariate %nin% c("diabetes.status.tri", "log.bmi"), ]
cm <- cm$covariate

#' load covariate data for easier merging to mega-eset
covar.directory <- paste0(bp, "/02_projekte/1703_ge_metab_a1_b3_sorbs/171113_CovarAnnotTab/results/RData_objects/")
covar.file.name <- "a1_b3_sorb_allCovariatesFinal"
load.this.covar <- newest_file(look_for = covar.file.name, directory = covar.directory, print_full = T)
what <- load(load.this.covar)
what
load(load.this.covar)
dat.covar <- all.covars

# correct cohort identifier
dat.covar[id %in% covar.b3.ami$id, cohort := "ami"]

#' # Association Analyses

#' I only need the approach that uses combat considered and limma-adjusted covariates!

#' First merge the metabolite data to the ESet and then compute association. All steps integrated in this function. Function output will be omitted due to length. This will compute a linear regression model Gx-probe~Metabolite for each possible combination of the 97 metabolites and the QC-OK probes in each of the cohorts.

#' ## Single cohorts with covariates

#' This is for the preadjusted data to be analysed with covariates in the LIMMA model

# merge to eset
for(i in cm) {
# for(i in adjust$covariate) {
  pData(eset.a1.f)[,i] <- dat.covar[match_hk(rownames(pData(eset.a1.f)), dat.covar$id), get(i)]
  pData(eset.b3.ami.f)[,i] <- dat.covar[match_hk(rownames(pData(eset.b3.ami.f)), dat.covar$id), get(i)]
  pData(eset.b3.heart.f)[,i] <- dat.covar[match_hk(rownames(pData(eset.b3.heart.f)), dat.covar$id), get(i)]
  pData(eset.sorb.f)[,i] <- dat.covar[match_hk(rownames(pData(eset.sorb.f)), dat.covar$id), get(i)]
}

# is a factor!
pData(eset.a1.f)$sex <- as.character(pData(eset.a1.f)$sex)
pData(eset.b3.ami.f)$sex <- as.character(pData(eset.b3.ami.f)$sex)
pData(eset.b3.heart.f)$sex <- as.character(pData(eset.b3.heart.f)$sex)
pData(eset.sorb.f)$sex <- as.character(pData(eset.sorb.f)$sex)

#' # Table 1
# Table 1 ----

# TODO: EDIT ----

# create a summary of the data
source(here("../functions/cov_char.R"))

dat.covar[, bmi := exp(log.bmi)]

# filter for correct ids
dat.covar[id %in% covar.b3.ami$id, cohort := "ami"]
table(dat.covar$cohort)

# filter to total overlap
table1dat <- dat.covar[id %in% c(
  metab.a1$id, 
  metab.b3.ami$id,
  metab.b3.heart$id,
  metab.sorb$id
)]
table(table1dat$cohort)

# all covars and continuous covars
covars <- c("age", "fasting.hours", "sex",
            "bmi", "hematocrit", "monocytes.percent",
             "neutrophils.percent", "diabetes.status.tri"
             )

cont.covars <- c(
  "age", "fasting.hours", "bmi",
  "hematocrit", "monocytes.percent",
  "neutrophils.percent")

# run the function
char <- cov_char(dat = table1dat,
                 covars = covars,
                 continuousCovars = cont.covars)

# format continuous studies
tab1cont <- char$continuous[,.(
  "Covariate" = term,
  "LIFE Adult" = combined.iqr_a1,
  "LIFE Heart" = combined.iqr_b3,
  "LIFE AMI" = combined.iqr_ami,
  "Sorb Study" = combined.iqr_sorb
)]
tab1cont[Covariate == "fasting.hours", `Sorb Study` := ">8"]

# format continuous studies
tab1bin <- char$discrete[,.(
  "Covariate" = term,
  "LIFE Adult" = all.n_a1,
  "LIFE Heart" = all.n_b3,
  "LIFE AMI" = all.n_ami,
  "Sorb Study" = all.n_sorb
)]

# create header with N and location
tab1head <- data.table(
  "Covariate" = c("Area of collection", "N"),
  "LIFE Adult" = c("Leipzig, Germany", as.character(nrow(metab.a1))),
  "LIFE Heart" = c("Leipzig, Germany", as.character(nrow(metab.b3.heart))),
  "LIFE AMI" = c("Leipzig, Germany", as.character(nrow(metab.b3.ami))),
  "Sorb Study" = c("Upper Lusatia", as.character(nrow(metab.sorb)))
)

# bind everything together
table1 <- rbindlist(
  list(tab1head, tab1cont,tab1bin)
)

# proper names
table1$Covariate <- c(
  "Area of collection", 
  "N",
  "Age (years)", 
  "BMI (kg/m²)", 
  "Fasting hours (hours)", 
  "Hematocrit (%)",
  "Monocytes (%)", 
  "Neutrophils (%)", 
  "Diabetes status (yes)", 
  "Sex (female)"
  )

# reorder rows
table1 <- table1[c(1,2,3,10,4,9,5,6,7,8),]

toolboxH::WriteXLS_hk(
  x = table1,
  SheetNames = "Study Characteristics",
  ExcelFileName = dpx("study_characteristics.xls", "ppr/supp_tables/")
)

#+ limma.covar.a1, results="hide", message=FALSE, include=FALSE
# Adult
limma.covar.a1 <- each_cohort_limma(
  gx.dat = eset.a1.f,
  gx.annot = probe.adult,
  metab.dat = metab.a1,
  metabolites = annot$metab.annot$metabolite,
  cohortID = "Adult w/ covariates",
  covariates = cm
)

#+ limma.covar.ami, results="hide", message=FALSE, include=FALSE
# Heart-AMI
limma.covar.b3.ami <- each_cohort_limma(
  gx.dat = eset.b3.ami.f,
  gx.annot = probe.heart,
  metab.dat = metab.b3.ami,
  metabolites = annot$metab.annot$metabolite,
  cohortID = "Heart-AMI w/ covariates",
  covariates = cm
)

#+ limma.covar.heart, results="hide", message=FALSE, include=FALSE
# Heart
limma.covar.b3 <- each_cohort_limma(
  gx.dat = eset.b3.heart.f,
  gx.annot = probe.heart,
  metab.dat = metab.b3.heart,
  metabolites = annot$metab.annot$metabolite,
  cohortID = "Heart w/ covariates",
  covariates = cm
  
)

#+ limma.covar.sorb, results="hide", message=FALSE, include=FALSE
# Sorb
limma.covar.sorb <- each_cohort_limma(
  gx.dat = eset.sorb.f,
  gx.annot = probe.sorb,
  metab.dat = metab.sorb,
  metabolites = annot$metab.annot$metabolite,
  cohortID = "Sorb w/ covariates",
  covariates = setdiff(cm, "fasting.hours")
  
)

# take the output from all cohorts and rbindlist them first for each metab and then for the cohorts
# save the association results in one table
limma.covar.res <- rbindlist(
  list(
    limma.covar.a1$res,
    limma.covar.b3$res,
    limma.covar.b3.ami$res,
    limma.covar.sorb$res
  )
)

# get beta mse from CI
limma.covar.res[, beta.se := (CI.R - logFC) / 1.96 ]

# save the summary with the eta1 statistic
limma.covar.summary <- rbindlist(
  list(
    limma.covar.a1$summary,
    limma.covar.b3$summary,
    limma.covar.b3.ami$summary,
    limma.covar.sorb$summary
  )
)

# save both the summary and the association statistics
save_csv_carl(file = limma.covar.summary, 
              file_name = "a1_b3_sorb_limmaAssocWithCovarsSummary", 
              subfolder = "res")
# limma.covar.res <- fread("res/191017_a1_b3_sorb_limmaAssocWithCovarsResults.csv")
# limma.covar.summary <- fread("res/191017_a1_b3_sorb_limmaAssocWithCovarsSummary.csv")
save_csv_carl(file = limma.covar.res, 
              file_name = "a1_b3_sorb_limmaAssocWithCovarsResults",
              subfolder = "res")

# prepare table with summary
mytab1 <- copy(limma.covar.summary)
oldnames <- names(limma.covar.summary)
mytab1[, `:=`(
  `Study` = cohort,
  `Metabolite` = metabolite,
  `Eta1` = eta1,
  `Eta1 SE` = eta.se,
  `Number of Probes` = n_probes,
  `Number of Genes` = n_genes,
  `Minimal p-value` = min_p,
  `Bonferroni-significant genes` = n_genes_bonf,
  `Minimum q-value` = min_qval,
  `Number of genes at FDR=5%` = n_genes_q5,
  `Number of genes at FDR=20%` = n_genes_q20,
  `Minimal p-value of upregulated genes` =  min_p_up,
  `Bonferroni-significant upregulated genes` =  n_genes_bonf_up,
  `Minimum q-value of upregulated genes` = min_qval_up,
  `Number of upregulated genes at FDR=5%` = n_genes_q5_up,
  `Number of upregulated genes at FDR=20%` = n_genes_q20_up,
  `Minimal p-value of downregulated genes` = min_p_down,
  `Bonferroni-significant downregulated genes` = n_genes_bonf_down,
  `Minimum q-value of downregulated genes` = min_qval_down,
  `Number of downregulated genes at FDR=5%` = n_genes_q5_down,
  `Number of downregulated genes at FDR=20%` = n_genes_q20_down,
  `Number of genes up- and downregulated genes at FDR=5%` = n_genes_q5_upANDdown,
  `Number of genes up- and downregulated genes at FDR=20%` = n_genes_q20_upANDdown
)]

# remove redundant names
set(
  x = mytab1, 
  j = oldnames, 
  value = NULL
)

# rename cohorts
mytab1[, Study := gsub(x = Study, pattern = "\\s.+", replacement = "")]
mytab1$Study %>% unique

# rename metabolites
load("/net/ifs1/san_projekte/projekte/genstat/02_projekte/1703_ge_metab_a1_b3_sorbs/1805_gxMetaboliteAssociation/obj/200306_metaboliteNames.RData")

# replace R id metabolite names with the correct ones
m <- match(mytab1$Metabolite,agn.all$rid)
mytab1[, `Metabolite abbreviation` := agn.all[m, abbr]]
mytab1[, Metabolite := agn.all[m, full.name]]
mynames <- c("Study", "Metabolite abbreviation", "Metabolite")
mytab1 <- mytab1[,.SD,.SDcols=c(mynames, setdiff(names(mytab1), mynames))]
WriteXLS_hk(x = mytab1, 
            ExcelFileName = "ppr/supp_tables/201104_singleStudySummary.xls", 
            SheetNames = "Single Study Association")


# prepare table with all association results
mytab2 <- copy(limma.covar.res)
oldnames <- names(limma.covar.res)

# replace R id metabolite names with the correct ones
m <- match(mytab2$var,agn.all$rid)
mytab2[, metabr := agn.all[m, abbr]]
mytab2[, met := agn.all[m, full.name]]

# order and rename columns
mytab2[, `:=`(
  Study = cohort,
  `Metabolite abbreviation` = metabr,
  Metabolite = met,
  ILMN = ilmn,
  `Entrez ID` = EntrezReannotated,
  `HGNC Symbol` = SymbolReannotated_orgHsEg,
  `T-statistic` = t,
  `log-Fold Change` = logFC,
  `Fold Change` = FC,
  `Standard Error` = beta.se,
  `Lower 95% CI` = CI.L,
  `Upper 95% CI` = CI.R,
  `Average Expression` = AveExpr,
  `P-value` = P.Value,
  `Level 1 adjusted P-value` = hfdr.1,
  `Level 2 adjusted P-value` = hfdr.2,
  `Significance at hFDR=5%` = sig.hfdr,
  N = n_individuals
)]

# remove old names
set(x = mytab2, j = oldnames, value = NULL)
mytab2$metabr <- NULL
mytab2$met <- NULL

# rename cohorts
mytab2[, Study := gsub(x = Study, pattern = "\\s.+", replacement = "")]
mytab2$Study %>% unique

fwrite(x = mytab2,
       file = "ppr/supp_tables/201106_allSingleStudySumStats.tsv",
       sep = "\t",
       quote = "auto"
)


devtools::session_info()
