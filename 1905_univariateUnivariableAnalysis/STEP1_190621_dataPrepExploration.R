#' ---
#' title: "Variable Selection in Gene Expression Data"
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
  if(grepl(x= getwd(),pattern =  "/home/carl/")){
    bp <- "/home/carl/imise/projekte/genstat/"
  } else {
    bp <- "/net/ifs2/san_projekte/projekte/genstat/"
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
  "readxl",
  "ggplot2",
  "extrafont",
  "lumi",
  "limma",
  "MatrixEQTL"
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
source(here("../functions/lm_via_MatrixEQTL.R"))
option_setup()

#+ include=F
#==========#
#   DONE   #
#==========#

#' # Data

#' ## Annotation and ID Data

# load the functino that loads all annotation tables
source(paste0(bp, "02_projekte/1703_ge_metab_a1_b3_sorbs/functions/all_annotation_tables.R"))
annot <- all_annotation_tables(mount.point = bp)

#' ## IDs

# load adult IDs
ids.a1 <- fread(paste0(bp, "/02_projekte/1703_ge_metab_a1_b3_sorbs/results/s601_1_CONFIDENTIAL_schluessel.txt"))

# load Heart IDs
ids.b3 <- readxl::read_excel("/net/ifs1/san_projekte/projekte/genstat/02_projekte/1511_lifeb3_gwas_gefaesse/pheno/170898_Rekrutierungsmodus.xlsx")
setDT(ids.b3)
ids.b3 <- ids.b3[!is.na(ID),]

#' ## Gx probe and sample annotation

#' 190617: I re-did the preprocessing of all data, removing diabetes status and log-BMI from the adjusted covariates. Therefore I need to repeat the whole analysis. Additionally, I want to seperate data preparation and analysis. 

#' ### LIFE Adult

sample.adult <- fread(paste0(bp, "02_projekte/1705_ge_lifea1v2/180719_a1_gePreProForMetab/tosend_NoBMI_NoDiabetes/s401_1_sampleannot_HT12v4.txt"))
probe.adult  <- fread(paste0(bp, "02_projekte/1705_ge_lifea1v2/180719_a1_gePreProForMetab/tosend_NoBMI_NoDiabetes/s401_1_probeannot_HT12v4.txt"))
sample.adult.nc <- fread(paste0(bp, "02_projekte/1705_ge_lifea1v2/180719_a1_gePreProForMetab/tosend_noCovariates/s401_1_sampleannot_HT12v4.txt"))
probe.adult.nc  <- fread(paste0(bp, "02_projekte/1705_ge_lifea1v2/180719_a1_gePreProForMetab/tosend_noCovariates/s401_1_probeannot_HT12v4.txt"))

#' ### LIFE Heart
sample.heart.old <- fread(paste0(bp, "/02_projekte/1703_ge_metab_a1_b3_sorbs/1805_gxMetaboliteAssociation/dat/s12_sampleannot_HT12v4.txt"))
sample.heart <- fread(paste0(bp, "02_projekte/1102_ge_lifeb3/07_prepro_ami_allb3_for_metab/tosend_NoBMI_NoDiabetes/sampleannot_HT12v4.txt"))
probe.heart  <- fread(paste0(bp, "02_projekte/1102_ge_lifeb3/07_prepro_ami_allb3_for_metab/tosend_NoBMI_NoDiabetes/probeannot_HT12v4.txt"))
sample.heart.nc <- fread(paste0(bp, "02_projekte/1102_ge_lifeb3/07_prepro_ami_allb3_for_metab/tosend_noCovariates/sampleannot_HT12v4.txt"))
probe.heart.nc  <- fread(paste0(bp, "02_projekte/1102_ge_lifeb3/07_prepro_ami_allb3_for_metab/tosend_noCovariates/probeannot_HT12v4.txt"))

#' 190621: Consider the subgroup in the Heart cohort - AMI and the rest as
#' seperate cohort. Annotation is found here

# should be identical to the gx annotation
group.heart <- fread("../171113_CovarAnnotTab/results/190620_amiGroupMatching.csv")

#' ### Sorb

sample.sorb <- fread(paste0(bp, "02_projekte/1706_ge_sorben_v2/181009_sorb_gePreProForMetab/tosend_NoBMI_NoDiabetes/s401_1_sampleannot_HT12v4.txt"))
probe.sorb  <- fread(paste0(bp, "02_projekte/1706_ge_sorben_v2/181009_sorb_gePreProForMetab/tosend_NoBMI_NoDiabetes/s401_1_probeannot_HT12v4.txt"))
sample.sorb.nc <- fread(paste0(bp, "02_projekte/1706_ge_sorben_v2/181009_sorb_gePreProForMetab/tosend_noCovariates/s401_1_sampleannot_HT12v4.txt"))
probe.sorb.nc  <- fread(paste0(bp, "02_projekte/1706_ge_sorben_v2/181009_sorb_gePreProForMetab/tosend_noCovariates/s401_1_probeannot_HT12v4.txt"))

#' ## Gx Data

#' ### Gx Remapping QC data

#' This file is a list of all QC-OK probes. Filtering for those is the current
#' strategy.

gx.qc.ok <- fread(paste0(bp, "07_programme/rtools/1807_gx_tools/mappingHT12/results/102_8_remappingHT12_INGENUITY_goodilmns.txt"))

#' ### LIFE Adult

load(paste0(bp, "02_projekte/1705_ge_lifea1v2/180719_a1_gePreProForMetab/tosend_NoBMI_NoDiabetes/s401_1_expressionset_preprocessed.Rdata"))

# rename and delete the original
eset.a1 <- eset_preproc
rm(eset_preproc)

# unadjusted data
load(paste0(bp, "02_projekte/1705_ge_lifea1v2/180719_a1_gePreProForMetab/tosend_noCovariates/s401_1_expressionset_preprocessed.Rdata"))

# rename and delete the original
eset.nc.a1 <- eset_preproc
rm(eset_preproc)

#' ### LIFE Heart

load(paste0(bp, "/02_projekte/1102_ge_lifeb3/07_prepro_ami_allb3_for_metab/tosend_NoBMI_NoDiabetes/expressionset_preprocessed.Rdata"))

# rename and delete the original
eset.b3 <- eset_preproc
rm(eset_preproc)

load(paste0(bp, "/02_projekte/1102_ge_lifeb3/07_prepro_ami_allb3_for_metab/tosend_noCovariates/expressionset_preprocessed.Rdata"))

# rename and delete the original
eset.nc.b3 <- eset_preproc
rm(eset_preproc)

#' ### Sorb cohort

#'  This is the unrelated sorb gx expression matrix
nf <- newest_file(
  look_for = "sorb_GxRelationshipAdjust",
  directory = paste0(bp, "/02_projekte/1703_ge_metab_a1_b3_sorbs/170623_sorben_UnrelateGxMetab/Results/"),
  print_full = T
)

# there are still -1000s in this file
gx.sorb <- fread(nf)

# remove the -10000
invisible(
  change_in_dt(
    dat = gx.sorb,
    from = -10000,
    to = NA,
    change_in_dat = T
  )
)

# check for remaining -10000s
negs <- gx.sorb[, lapply(.SD, function(i){sum(i<= -100, na.rm=TRUE)})]
any(negs >0)

#'  This is the unadjusted and unrelated sorb gx expression matrix 
nf <- newest_file(
  look_for = "sorb_GxRelationshipAdjust",
  directory = paste0(bp, "/02_projekte/1703_ge_metab_a1_b3_sorbs/170623_sorben_UnrelateGxMetab/Results_noCovariates/"),
  print_full = T
)

# there are still -1000s in this file
gx.nc.sorb <- fread(nf)

# remove the -10000
invisible(
  change_in_dt(
    dat = gx.nc.sorb,
    from = -10000,
    to = NA,
    change_in_dat = T
  )
)

# check for remaining -10000s
negs <- gx.nc.sorb[, lapply(.SD, function(i){sum(i<= -100, na.rm=TRUE)})]
any(negs >0)

# load(paste0(bp, "02_projekte/1706_ge_sorben_v2/181009_sorb_gePreProForMetab/tosend/s401_1_expressionset_preprocessed.Rdata"))

#' ## Metabolite Data

# these are the metabolites adjusted for the gx covariates
metab.dir <- paste0(bp, "/02_projekte/1703_ge_metab_a1_b3_sorbs/170824_MetabAnnotTab/1906_adjustMetabsForGx_NoBMI_NoDiabetes/results")
# metab.dir <- paste0(bp, "/02_projekte/1703_ge_metab_a1_b3_sorbs/170829_ConfounderAnalysis/190207_ExtraAnalysis/tosend")
metab.loc <- newest_file(look_for = "a1_b3_sorb_metabolites.csv",
                         directory = metab.dir,
                         print_full = TRUE)
dat.metab <- fread(metab.loc)

#' ## Metabolite Data Not adjusted for covariates

# these are the metabolites not adjusted for the gx covariates
metab.dir <- paste0(bp, "/02_projekte/1703_ge_metab_a1_b3_sorbs/170824_MetabAnnotTab/1907_adjustMetabsForGx_noCovariates/results")
# metab.dir <- paste0(bp, "/02_projekte/1703_ge_metab_a1_b3_sorbs/170829_ConfounderAnalysis/190207_ExtraAnalysis/tosend")
metab.loc <- newest_file(look_for = "a1_b3_sorb_metabolites.csv",
                         directory = metab.dir,
                         print_full = TRUE)
dat.metab.nc <- fread(metab.loc)

#' ## Covariate Data 

#' For a sensitivity analysis I want to add the adjusted covariates to the limma model. I need to load them from wherever.

# load
covar.directory <- paste0(bp, "/02_projekte/1703_ge_metab_a1_b3_sorbs/171113_CovarAnnotTab/results/RData_objects/")
covar.file.name <- "a1_b3_sorb_allCovariatesFinal"
load.this.covar <- newest_file(look_for = covar.file.name, directory = covar.directory, print_full = T)
what <- load(load.this.covar)
what
load(load.this.covar)

#' Also load the covariate names to be adjusted for
directory <- "/net/ifs1/san_projekte/projekte/genstat/02_projekte/1703_ge_metab_a1_b3_sorbs/170829_ConfounderAnalysis/190207_ExtraAnalysis/results"
nf <- newest_file(look_for = "FullModelWithoutQuotientsInclHighMissings", directory = directory, print_full = T)
adjust <- fread(nf)
cm <- adjust[covariate %nin% c("diabetes.status.tri", "log.bmi"), ]
cm <- cm$covariate

#' # Explore data

#' ## ID overlaps

#' ### Adult

# get the preprocessed IDs, find overlap
metab.ids.a1 <- dat.metab[cohort=="a2",id]
gx.ids.a1 <- sample.adult[in_study==T, finalID]
gx.metab.overlap.a1 <- venn2(metab.ids.a1, gx.ids.a1)

#' ### Heart

# 190513 I renamed heart gx samples in the gx prepro, therefore no matching is
# necessary in this step. I just rename the finalID column in the sample
# annotation to 'aliquot' for downstream compatibility

sample.heart[, aliquot := finalID]

# get the overlap
metab.ids.b3 <- dat.metab[cohort=="b3",id]
gx.ids.b3 <- sample.heart[in_study==T, finalID]
gx.metab.overlap.b3 <- venn2(metab.ids.b3, gx.ids.b3)

#' Done. Now I have the Heart data properly matched with all necessary IDs for
#' my Gx-Metab association.

#' ### Sorb cohort

#' First I need to create an expression set from my Sorb gx data: For the
#' relationship adjustment, the data was processed as a DT, but for processing
#' in LIMMA, I'd need the data in the lumi ExpressionSet format.

# get sample and feature names
gx.sorb.probes <- names(gx.sorb)[-1]
gx.sorb.samples <- gx.sorb$barcode
gx.sorb.pop <- gx.sorb$pop
gx.sorb$pop <- NULL

# gx data needs to be transposed
gx.sorb.t <- dcast.data.table(melt(gx.sorb,
                                   id.vars="barcode",
                                   variable.name="probes"),
                              formula=probes~barcode)

# remove untransposed gx.sorb
rm(gx.sorb)

# create matrix from DT
gx.sorb.t.m <- as.matrix(gx.sorb.t, rownames=1)
rm(gx.sorb.t)

# initiate ESet
eset.sorb <- ExpressionSet(assayData = gx.sorb.t.m)
rm(gx.sorb.t.m)

#' Now that I have consolidated the relationship adjusted GX data from the sorb
#' cohort, I will figure out the gx-metabolite overlap

# get the metab IDs of the sorb
metab.ids <- dat.metab[cohort=="sorb",id]

# get the sorb gx OK samples
gx.ids.qc.ok <- sample.sorb[in_study==T, old_ID]

# get the overlap
gx.metab.overlap.sorb <- venn2(as.character(metab.ids), as.character(gx.ids.qc.ok))

#' ## Unify ID for joint modelling in the metabolites

#' ### Adult

#' I should just overwrite the IDs in the metabolite file so I don't have to
#' mess around with the ExpressionSet

head(gx.metab.overlap.a1$q1)

# take only the metab data I have qc ok gx data for
metab.a1 <- dat.metab[cohort=="a1" & id %in% gx.metab.overlap.a1$q1, ]

#' ### Heart

head(gx.metab.overlap.b3$q1)
metab.b3 <- dat.metab[cohort=="b3" & id %in% gx.metab.overlap.b3$q1,]

#' I want to add the subcohort identifier to the heart metabolites to seperate
#' the two as distinct cohorts

# create a binary identifier ami T/F
filt.ami <- group.heart$group == "ami"
group.heart$ami <- ifelse(filt.ami, "ami", "heart")

# match the ami identifier to the metabolites
m1 <- match(metab.b3$id, group.heart$aliquot)
metab.b3$subgroup <- group.heart[m1, ami]

# seperate heart into ami and non-ami
metab.b3.all <- metab.b3
metab.b3.ami <- metab.b3[subgroup=="ami", ]
metab.b3.heart <- metab.b3[subgroup!="ami", ]

#' ### Sorb

head(gx.metab.overlap.sorb$q1)
metab.sorb <- dat.metab[cohort=="sorb" & id %in% gx.metab.overlap.sorb$q1, ]

#' ## Unify ID and QC IDs in the ESets

#' I have filtered the metabolite data to match the gx data, including a unique
#' identifier column (gx.id), which allows me to map each metabolite sample to
#' their respective gx sample. Now I need to filter the respective ESet down
#' and allign the data.

#' ### Adult

# these are the metab ids
gx.order <- match(metab.a1$id, sampleNames(eset.a1))

# these are the probes with good QC criteria (see probe annotation xlsx sheet
# readme) for explanation of column differences
table(probe.adult$goodexpressedprobe_all)
table(probe.adult$perfectprobe_all) # probes including old remapping

# get all QC-OK probes from the annotation
good.probes.a1 <- probe.adult[goodexpressedprobe_all==T, ilmn]
good.probes.ipa.remap.a1 <- good.probes.a1[good.probes.a1 %in% gx.qc.ok$mappingOK_INGENUITY]
table(good.probes.a1 %in% gx.qc.ok$mappingOK_INGENUITY)
length(good.probes.ipa.remap.a1)

# these are the filtered gx probe names
probe.filter.a1 <- featureNames(eset.a1) %in% good.probes.ipa.remap.a1
table(probe.filter.a1)

# filter probes and filter and re-arrange samples
eset.a1.f <- eset.a1[probe.filter.a1,gx.order]

#' ### Heart

#' Set the correct filters for the Heart cohort. Create two sepereate Esets for
#' ami/non-ami Heart samples

#' #### AMI

# reorder metabolites
metab.order.b3ami <- match(metab.b3.ami$id, sampleNames(eset.b3))

# get the good probes from the annotation
table(probe.heart$goodexpressedprobe_allsubgroups)
table(probe.heart$goodexpressedprobe_ami)
good.probes.b3 <- probe.heart[goodexpressedprobe_ami==T,ilmn]
good.probes.ipa.remap.b3 <- good.probes.b3[good.probes.b3 %in% gx.qc.ok$mappingOK_INGENUITY]
table(good.probes.b3 %in% gx.qc.ok$mappingOK_INGENUITY)
length(good.probes.ipa.remap.b3)
probe.filter.b3 <- featureNames(eset.b3) %in% good.probes.ipa.remap.b3
table(probe.filter.b3)

# filter probes and filter and re-arrange samples
eset.b3.ami.f <- eset.b3[probe.filter.b3,metab.order.b3ami]

#' #### Non-AMI

# reorder metabolites
metab.order.b3.heart <- match(metab.b3.heart$id, sampleNames(eset.b3))

# use the heart subgroup qc
# get the good probes from the annotation
table(probe.heart$goodexpressedprobe_allsubgroups)
table(probe.heart$goodexpressedprobe_heart)
good.probes.b3 <- probe.heart[goodexpressedprobe_heart==T,ilmn]
good.probes.ipa.remap.b3 <- good.probes.b3[good.probes.b3 %in% gx.qc.ok$mappingOK_INGENUITY]
table(good.probes.b3 %in% gx.qc.ok$mappingOK_INGENUITY)
length(good.probes.ipa.remap.b3)
probe.filter.b3 <- featureNames(eset.b3) %in% good.probes.ipa.remap.b3
table(probe.filter.b3)

# filter probes and filter and re-arrange samples
eset.b3.heart.f <- eset.b3[probe.filter.b3,metab.order.b3.heart]

#' ### Sorb

# reorder metabolites
metab.order.sorb <- match(metab.sorb$id, sampleNames(eset.sorb))

# get the good probes from the annotation
table(probe.sorb$goodexpressedprobe_allsubgroups)
good.probes.sorb <- probe.sorb[goodexpressedprobe_allsubgroups==T,ilmn]
good.probes.ipa.remap.sorb <- good.probes.sorb[good.probes.sorb %in% gx.qc.ok$mappingOK_INGENUITY]
table(good.probes.sorb %in% gx.qc.ok$mappingOK_INGENUITY)
length(good.probes.ipa.remap.sorb)
probe.filter.sorb <- featureNames(eset.sorb) %in% good.probes.ipa.remap.sorb
table(probe.filter.sorb)

# filter probes and filter and re-arrange samples
eset.sorb.f <- eset.sorb[probe.filter.sorb,metab.order.sorb ]

#' ## Unify covariate Info

#' ### Adult

# get the respective covariates for each cohort
m1 <- match(metab.a1$id, all.covars[cohort=="a1", id])
metab.a1 %>% dim
m1 %>% length
covar.a1 <- all.covars[cohort=="a1", ][m1, .SD, .SDcols = c("id", cm)]
dim(covar.a1)
dim(metab.a1)
dim(eset.a1.f)

#' ### Heart

m2 <- match(metab.b3.heart$id, all.covars[cohort=="b3", id])
covar.b3.heart <- all.covars[cohort=="b3", ][m2, .SD, .SDcols = c("id", cm)]
dim(covar.b3.heart)
dim(metab.b3.heart)
dim(eset.b3.heart.f)

#' ### AMI

m3 <- match(metab.b3.ami$id, all.covars[cohort=="b3", id])
covar.b3.ami <- all.covars[cohort=="b3", ][m3, .SD, .SDcols = c("id", cm)]
dim(covar.b3.ami)
dim(metab.b3.ami)
dim(eset.b3.ami.f)

#' ### Sorb

m4 <- match(metab.sorb$id, all.covars[cohort=="sorb", id])
covar.sorb <- all.covars[cohort=="sorb", ][m4, .SD, .SDcols = c("id", setdiff(cm, "fasting.hours"))]
dim(covar.sorb)
dim(metab.sorb)
dim(eset.sorb.f)

#=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=#
### NOW EVERYTHING FOR THE UNADJUSTED DATA  ###
#=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=#

#' # Explore data

#' ## ID overlaps

#' ### Adult

# get the preprocessed IDs, find overlap
metab.ids.a1 <- dat.metab.nc[cohort=="a1",id]
gx.ids.a1 <- sample.adult.nc[in_study==T, finalID]
gx.metab.overlap.a1 <- venn2(metab.ids.a1, gx.ids.a1)

#' ### Heart

#' set the aliquot ID for heart metabolites

v1 <- venn2(ids.b3$LEH_ID, dat.metab.nc[cohort=="b3", id])

# remove the three ids not found
dat.metab.nc <- dat.metab.nc[id %nin% v1$q3, ]

# match the aliquot instead LEH id
matched1 <- match(dat.metab.nc[cohort=="b3", id], ids.b3$LEH_ID)
dat.metab.nc[cohort=="b3", id := ids.b3[matched1, `BE-NR`]]

sample.heart.nc[, aliquot := finalID]

# get the overlap
metab.ids.b3 <- dat.metab.nc[cohort=="b3",id]
gx.ids.b3 <- sample.heart.nc[in_study==T, finalID]
gx.metab.overlap.b3 <- venn2(metab.ids.b3, gx.ids.b3)

#' Done. Now I have the Heart data properly matched with all necessary IDs for
#' my Gx-Metab association.

#' ### Sorb cohort

#' First I need to create an expression set from my Sorb gx data: For the
#' relationship adjustment, the data was processed as a DT, but for processing
#' in LIMMA, I'd need the data in the lumi ExpressionSet format.

# get sample and feature names
gx.sorb.probes <- names(gx.nc.sorb)[-1]
gx.sorb.samples <- gx.nc.sorb$barcode
gx.sorb.pop <- gx.nc.sorb$pop
gx.nc.sorb$pop <- NULL

# gx data needs to be transposed
gx.nc.sorb.t <- dcast.data.table(melt(gx.nc.sorb,
                                   id.vars="barcode",
                                   variable.name="probes"),
                              formula=probes~barcode)

# remove untransposed gx.sorb
rm(gx.nc.sorb)

# create matrix from DT
gx.nc.sorb.t.m <- as.matrix(gx.nc.sorb.t, rownames=1)
rm(gx.nc.sorb.t)

# initiate ESet
eset.nc.sorb <- ExpressionSet(assayData = gx.nc.sorb.t.m)
rm(gx.nc.sorb.t.m)

#' Now that I have consolidated the relationship adjusted GX data from the sorb
#' cohort, I will figure out the gx-metabolite overlap

# get the metab IDs of the sorb
metab.ids <- dat.metab.nc[cohort=="sorb",id]

# get the sorb gx OK samples
gx.ids.qc.ok <- sample.sorb.nc[in_study==T, old_ID]

# get the overlap
gx.metab.overlap.sorb <- venn2(as.character(metab.ids), as.character(gx.ids.qc.ok))

#' ## Unify ID for joint modelling in the metabolites

#' ### Adult

#' I should just overwrite the IDs in the metabolite file so I don't have to
#' mess around with the ExpressionSet

head(gx.metab.overlap.a1$q1)

# take only the metab data I have qc ok gx data for
metab.nc.a1 <- dat.metab.nc[cohort=="a1" & id %in% gx.metab.overlap.a1$q1, ]

#' ### Heart

head(gx.metab.overlap.b3$q1)
metab.nc.b3 <- dat.metab.nc[cohort=="b3" & id %in% gx.metab.overlap.b3$q1,]

#' I want to add the subcohort identifier to the heart metabolites to seperate
#' the two as distinct cohorts

# create a binary identifier ami T/F
filt.ami <- group.heart$group == "ami"
group.heart$ami <- ifelse(filt.ami, "ami", "heart")

# match the ami identifier to the metabolites
m1 <- match(metab.nc.b3$id, group.heart$aliquot)
metab.nc.b3$subgroup <- group.heart[m1, ami]

# seperate heart into ami and non-ami
metab.nc.b3.all <- metab.nc.b3
metab.nc.b3.ami <- metab.nc.b3[subgroup=="ami", ]
metab.nc.b3.heart <- metab.nc.b3[subgroup!="ami", ]

#' ### Sorb

head(gx.metab.overlap.sorb$q1)
metab.nc.sorb <- dat.metab.nc[cohort=="sorb" & id %in% gx.metab.overlap.sorb$q1, ]

#' ## Unify ID and QC IDs in the ESets

#' I have filtered the metabolite data to match the gx data, including a unique
#' identifier column (gx.id), which allows me to map each metabolite sample to
#' their respective gx sample. Now I need to filter the respective ESet down
#' and allign the data.

#' ### Adult

# these are the metab ids
gx.order <- match(metab.nc.a1$id, sampleNames(eset.nc.a1))

# these are the probes with good QC criteria (see probe annotation xlsx sheet
# readme) for explanation of column differences
table(probe.adult.nc$goodexpressedprobe_allsubgroups)
table(probe.adult.nc$perfectprobe_allsubgroups) # probes including old remapping

# get all QC-OK probes from the annotation
good.probes.a1 <- probe.adult.nc[goodexpressedprobe_allsubgroups==T, ilmn]
good.probes.ipa.remap.a1 <- good.probes.a1[good.probes.a1 %in% gx.qc.ok$mappingOK_INGENUITY]
table(good.probes.a1 %in% gx.qc.ok$mappingOK_INGENUITY)
length(good.probes.ipa.remap.a1)

# these are the filtered gx probe names
probe.filter.a1 <- featureNames(eset.nc.a1) %in% good.probes.ipa.remap.a1
table(probe.filter.a1)

# filter probes and filter and re-arrange samples
eset.nc.a1.f <- eset.nc.a1[probe.filter.a1, gx.order]

#' ### Heart

#' Set the correct filters for the Heart cohort. Create two sepereate Esets for
#' ami/non-ami Heart samples

#' #### AMI

# reorder metabolites
metab.order.b3ami <- match(metab.nc.b3.ami$id, sampleNames(eset.nc.b3))

# get the good probes from the annotation
table(probe.heart$goodexpressedprobe_allsubgroups)
table(probe.heart$goodexpressedprobe_ami)
good.probes.b3 <- probe.heart.nc[goodexpressedprobe_ami==T,ilmn]
good.probes.ipa.remap.b3 <- good.probes.b3[good.probes.b3 %in% gx.qc.ok$mappingOK_INGENUITY]
table(good.probes.b3 %in% gx.qc.ok$mappingOK_INGENUITY)
length(good.probes.ipa.remap.b3)
probe.filter.b3 <- featureNames(eset.nc.b3) %in% good.probes.ipa.remap.b3
table(probe.filter.b3)

# filter probes and filter and re-arrange samples
eset.nc.b3.ami.f <- eset.nc.b3[probe.filter.b3,metab.order.b3ami]

#' #### Non-AMI

# reorder metabolites
metab.order.b3.heart <- match(metab.nc.b3.heart$id, sampleNames(eset.nc.b3))

# use the heart subgroup qc
# get the good probes from the annotation
table(probe.heart.nc$goodexpressedprobe_allsubgroups)
table(probe.heart.nc$goodexpressedprobe_heart)
good.probes.b3 <- probe.heart.nc[goodexpressedprobe_heart==T,ilmn]
good.probes.ipa.remap.b3 <- good.probes.b3[good.probes.b3 %in% gx.qc.ok$mappingOK_INGENUITY]
table(good.probes.b3 %in% gx.qc.ok$mappingOK_INGENUITY)
length(good.probes.ipa.remap.b3)
probe.filter.b3 <- featureNames(eset.nc.b3) %in% good.probes.ipa.remap.b3
table(probe.filter.b3)

# filter probes and filter and re-arrange samples
eset.nc.b3.heart.f <- eset.nc.b3[probe.filter.b3, metab.order.b3.heart]

#' ### Sorb

# reorder metabolites
metab.order.sorb <- match(metab.nc.sorb$id, sampleNames(eset.nc.sorb))

# get the good probes from the annotation
table(probe.sorb.nc$goodexpressedprobe_allsubgroups)
good.probes.sorb <- probe.sorb.nc[goodexpressedprobe_allsubgroups==T,ilmn]
good.probes.ipa.remap.sorb <- good.probes.sorb[good.probes.sorb %in% gx.qc.ok$mappingOK_INGENUITY]
table(good.probes.sorb %in% gx.qc.ok$mappingOK_INGENUITY)
length(good.probes.ipa.remap.sorb)
probe.filter.sorb <- featureNames(eset.nc.sorb) %in% good.probes.ipa.remap.sorb
table(probe.filter.sorb)

# filter probes and filter and re-arrange samples
eset.nc.sorb.f <- eset.nc.sorb[probe.filter.sorb,metab.order.sorb ]

#' ## Unify covariate Info

#' ### Adult

# get the respective covariates for each cohort
m1 <- match(metab.nc.a1$id, all.covars[cohort=="a1", id])
metab.nc.a1 %>% dim
m1 %>% length
covar.nc.a1 <- all.covars[cohort=="a1", ][m1, .SD, .SDcols = c("id", cm)]
dim(covar.nc.a1)
dim(metab.nc.a1)
dim(eset.nc.a1.f)

#' ### Heart

m2 <- match(metab.nc.b3.heart$id, all.covars[cohort=="b3", id])
covar.nc.b3.heart <- all.covars[cohort=="b3", ][m2, .SD, .SDcols = c("id", cm)]
dim(covar.nc.b3.heart)
dim(metab.nc.b3.heart)
dim(eset.nc.b3.heart.f)

#' ### AMI

m3 <- match(metab.nc.b3.ami$id, all.covars[cohort=="b3", id])
covar.nc.b3.ami <- all.covars[cohort=="b3", ][m3, .SD, .SDcols = c("id", cm)]
dim(covar.nc.b3.ami)
dim(metab.nc.b3.ami)
dim(eset.nc.b3.ami.f)

#' ### Sorb

m4 <- match(metab.nc.sorb$id, all.covars[cohort=="sorb", id])
covar.nc.sorb <- all.covars[cohort=="sorb", ][m4, .SD, .SDcols = c("id", setdiff(cm, "fasting.hours"))]
dim(covar.nc.sorb)
dim(metab.nc.sorb)
dim(eset.nc.sorb.f)


#' # Save

#' Save all the data objects for analysis in the next step

save(file = "res/dat_for_STEP2.RData", 
     list = c(
       "eset.a1.f",
       "eset.b3.ami.f",
       "eset.b3.heart.f",
       "eset.sorb.f",
       "metab.a1",
       "metab.b3.ami",
       "metab.b3.heart",
       "metab.sorb",
       "covar.a1",
       "covar.b3.ami",
       "covar.b3.heart",
       "covar.sorb",
       "eset.nc.a1.f", # sample size without previous covariate adjustment!
       "eset.nc.b3.ami.f",
       "eset.nc.b3.heart.f",
       "eset.nc.sorb.f",
       "metab.nc.a1",
       "metab.nc.b3.ami",
       "metab.nc.b3.heart",
       "metab.nc.sorb",
       "covar.nc.a1",
       "covar.nc.b3.ami",
       "covar.nc.b3.heart",
       "covar.nc.sorb"
     )
)

#' # Session Info

devtools::session_info()

