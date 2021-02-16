#' ---
#' title: "Statistics necessary for mediation analysis"
#' author: "Carl Beuchel"
#' date: "`r format(Sys.time(), '%d %B, %Y')`"
#' output:
#'     html_document:
#'      code_download: true
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
  "here",
  "magrittr",
  "limma",
  "ggplot2",
  "lumi"
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
source(here("../functions/meta_analysis.R"))
source(here("../functions/meta_mediation_analysis.R"))
source(here("../functions/all_annotation_tables.R"))
source(here("../functions/option_setup.R"))
source(here("../functions/glm_loop.R"))
source(here("../functions/glm_mediation_loop.R"))
an <- all_annotation_tables(mount.point = "/net/ifs1/san_projekte/projekte/genstat/")
theme_set(theme_light(base_size = 15, base_family = "Helvetica"))
setDTthreads(1)
option_setup()

#+ include=F
#==========#
#   DONE   #
#==========#

#' # Tissue

#' metab: adult - dried whole blood
#' metab: heart - dried whole blood
#' metab: sorben - cell suspension after plasma centrifugation

#' gx: adult  - whole blood
#' gx: heart - PBMC
#' gx: sorben - PBMC

# Think about the direction of association ----

# Fallzahlen müssen gleich sein
# Testlauf mit starker assoc

# pheno ~ gx + metab
# 
# pheno ~ gx
# pheno ~ metab
# 
# gx    ~ metab
# metab ~ gx
# 
# gx    -> metab -> pheno
# metab -> gx    -> pheno

# Load data ----
#' # Load data

# Gx probe and sample Annotations ----
#' ## Gx probe and sample Annotations

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

# Load Gx Data and Metabolite Data ----
#' ## Load Gx Data and Metabolite Data

what.dat <- load("res/dat_for_STEP2.RData")
what.dat
load("res/dat_for_STEP2.RData")

# Load Covars ----
#' ## Load Covars

#'  load the covariate names to be adjusted to
directory <- "/net/ifs1/san_projekte/projekte/genstat/02_projekte/1703_ge_metab_a1_b3_sorbs/170829_ConfounderAnalysis/190207_ExtraAnalysis/results"
nf <- newest_file(look_for = "FullModelWithoutQuotientsInclHighMissings", 
                  directory = directory, print_full = T)
adjust <- fread(nf)
cm <- adjust[covariate %nin% c("diabetes.status.tri", "log.bmi"), ]
cm <- cm$covariate

# ap = all phenos
ap <- c("diabetes.status.tri", "log.bmi")

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
dat.covar$cohort %>% table

# Load Meta-Analysis results ----
#' ## Load Meta-Analysis results

nf <- newest_file(look_for = "STEP4_meta_resultsAnnotated$",subfolder = "res", print_full = TRUE)
dat <- fread(nf)

# Merge phenotypical data to esets ----
#' ## Merge phenotypical data to esets

# merge to eset
for(i in c(cm, ap) ) {
  pData(eset.a1.f)[,i] <- dat.covar[match_hk(pData(eset.a1.f)$sampleID, dat.covar$id), get(i)]
  pData(eset.b3.ami.f)[,i] <- dat.covar[match_hk(pData(eset.b3.ami.f)$sampleID, dat.covar$id), get(i)]
  pData(eset.b3.heart.f)[,i] <- dat.covar[match_hk(pData(eset.b3.heart.f)$sampleID, dat.covar$id), get(i)]
  pData(eset.sorb.f)[,i] <- dat.covar[match_hk(rownames(pData(eset.sorb.f)), dat.covar$id), get(i)]
}

# is a factor!
pData(eset.a1.f)$sex <- as.character(pData(eset.a1.f)$sex)
pData(eset.b3.ami.f)$sex <- as.character(pData(eset.b3.ami.f)$sex)
pData(eset.b3.heart.f)$sex <- as.character(pData(eset.b3.heart.f)$sex)
pData(eset.sorb.f)$sex <- as.character(pData(eset.sorb.f)$sex)
pData(eset.a1.f)$diabetes.status.tri <- as.character(pData(eset.a1.f)$diabetes.status.tri)
pData(eset.b3.ami.f)$diabetes.status.tri <- as.character(pData(eset.b3.ami.f)$diabetes.status.tri)
pData(eset.b3.heart.f)$diabetes.status.tri <- as.character(pData(eset.b3.heart.f)$diabetes.status.tri)
pData(eset.sorb.f)$diabetes.status.tri <- as.character(pData(eset.sorb.f)$diabetes.status.tri)

# sorb row names
pData(eset.sorb.f)$sampleID <- rownames(pData(eset.sorb.f))

# Gx-Overlap ----
#' # Gx-Overlap

v1 <- venn4(
  featureNames(eset.a1.f),
  featureNames(eset.b3.ami.f),
  featureNames(eset.b3.heart.f),
  featureNames(eset.sorb.f)
)

# Gx-Metab pairs to test ----
#' # Gx-Metab pairs to test

# test how many one-study results are significant
dat[hfdr.sigREM==T, .N]
dat[hfdr.sigREM==T & numberStudies>1, .N]
dat[hfdr.sigREM==T & numberStudies==1, .N]

#' Filter the Eset and Metabolite data down to the metabolites and gx-probes that associated significantly

# check all significant probes per metabolite, those to test in the Sobel test
# ptt = pairs to test
ptt <- dat[significant==TRUE, .(metabolite, markerID)]

# number of runs per gx-probe (to get correct sample-size)
ptt[, .N,by=markerID][order(N,decreasing = T)] %>% hh

# restrict ptt to probes found in all cohorts to restrict mediation analysis to fully analyzable probes for consistent sample sizes

# overview
v2 <- venn2(ptt$markerID, v1$q1,mytitle = "Probes with data in all 4 studies to be analysed")

# restrict the ptt probes
# ptt <- ptt[markerID %in% v1$q1, ]

# Imputation ----
#' # Imputation

#' To guarantee identical sample sizes across marginal associations, mean-impute missings.

# bmi - sorb, ami, adult
# pData(eset.a1.f)$log.bmi[is.na(pData(eset.a1.f)$log.bmi)] <- mean(pData(eset.a1.f)$log.bmi,na.rm=T)
# pData(eset.sorb.f)$log.bmi[is.na(pData(eset.sorb.f)$log.bmi)] <- mean(pData(eset.sorb.f)$log.bmi,na.rm=T)
# pData(eset.b3.ami.f)$log.bmi[is.na(pData(eset.b3.ami.f)$log.bmi)] <- mean(pData(eset.b3.ami.f)$log.bmi,na.rm=T)
# 
# # metabolites
# metab.a1[,unique(ptt$metabolite) := lapply(.SD, function(x){
#   x[is.na(x)] <- mean(x, na.rm=TRUE)
#   return(x)
# }),.SDcols=unique(ptt$metabolite)]
# metab.b3.ami[,unique(ptt$metabolite) := lapply(.SD, function(x){
#   x[is.na(x)] <- mean(x, na.rm=TRUE)
#   return(x)
# }),.SDcols=unique(ptt$metabolite)]
# metab.b3.heart[,unique(ptt$metabolite) := lapply(.SD, function(x){
#   x[is.na(x)] <- mean(x, na.rm=TRUE)
#   return(x)
# }),.SDcols=unique(ptt$metabolite)]
# metab.sorb[,unique(ptt$metabolite) := lapply(.SD, function(x){
#   x[is.na(x)] <- mean(x, na.rm=TRUE)
#   return(x)
# }),.SDcols=unique(ptt$metabolite)]

# Gx Standardisation ----
#' # Gx Standardisation 

# center-scale data
tmp1 <- apply(exprs(eset.a1.f), 1, function(x) {((x - mean(x)) / sd(x))})
exprs(eset.a1.f) <- t(tmp1)

exprs(eset.b3.ami.f) <- t(apply(exprs(eset.b3.ami.f), 1, function(x) {((x - mean(x)) / sd(x))}))
exprs(eset.b3.heart.f) <- t(apply(exprs(eset.b3.heart.f), 1, function(x) {((x - mean(x)) / sd(x))}))
exprs(eset.sorb.f) <- t(apply(exprs(eset.sorb.f), 1, function(x) {((x - mean(x)) / sd(x))}))

# Associations ----
#' # Associations

# Pheno~Gx Associations ----
#' ## Pheno~Gx Associations

# create a nested list containing all objects for easier looping
dat.list <- list(
  "Adult" = list(
    "gx.dat" = eset.a1.f[featureNames(eset.a1.f) %in% unique(ptt$markerID),],
    "metab.dat" = metab.a1[,.SD,.SDcols=c("id",unique(ptt$metabolite))],
    "gx.annot" = probe.adult,
    "covariates" = cm
  ),
  "AMI" = list(
    "gx.dat" = eset.b3.ami.f[featureNames(eset.b3.ami.f) %in% unique(ptt$markerID),],
    "metab.dat" = metab.b3.ami[,.SD,.SDcols=c("id",unique(ptt$metabolite))],
    "gx.annot" = probe.heart,
    "covariates" = cm
  ),
  "Heart" = list(
    "gx.dat" = eset.b3.heart.f[featureNames(eset.b3.heart.f) %in% unique(ptt$markerID),],
    "metab.dat" = metab.b3.heart[,.SD,.SDcols=c("id",unique(ptt$metabolite))],
    "gx.annot" = probe.heart,
    "covariates" = cm
  ),
  "Sorb" = list(
    "gx.dat" = eset.sorb.f[featureNames(eset.sorb.f) %in% unique(ptt$markerID),],
    "metab.dat" = metab.sorb[,.SD,.SDcols=c("id",unique(ptt$metabolite))],
    "gx.annot" = probe.adult,
    "covariates" = setdiff(cm, "fasting.hours")
  )
)

# check NAs in phenotypes
lapply(dat.list,function(x){
  tmp1 <- pData(x[["gx.dat"]])
  m1 <- tmp1$log.bmi %>% is.na %>% sum
  m2 <- tmp1$diabetes.status.tri %>% is.na %>% sum
  return(c("bmi" = m1,
           "t2d" = m2))
})

# check NAs in metabolites
lapply(dat.list,function(x){
  tmp1 <- x[["metab.dat"]]
  tmp2 <- tmp1[,sapply(.SD,function(x){
    sum(is.na(x))}),.SDcols=unique(ptt$metabolite)]
  return(tmp2[tmp2!=0])
})

# Data-Overlap ----
#' # Data-Overlap

#' I need to ensure a consistent sample size over each test pheno~metab in each cohort. Gx has no missings and can be disregarded in this. I'll need to get the sample IDs of each combination in each cohort and save the IDs for use in the test associations.

# ass = all sample sizes
ass <- expand.grid(
  metabolite = unique(ptt$metabolite),
  pheno = ap, 
  stringsAsFactors = FALSE
)
setDT(ass)
ass[, index:=1:.N]


# loop through each cohort
res <- lapply(
  list("Adult",
       "AMI",
       "Heart",
       "Sorb"), function(i){
         
         # get the required data for the column
         pheno <- pData(dat.list[[i]][["gx.dat"]])
         setDT(pheno)
         metab <- dat.list[[i]][["metab.dat"]]
         covars <- dat.list[[i]][["covariates"]]
         var3 <- covars # needed for the loop below
         
         # loop through the indices
         res.cohort <- lapply(ass$index, function(j){
           
           # get a copy to return as result
           ass.cohort <- copy(ass[index==(j)])
           
           # get the specific combination of pheno-metabolite
           var1 <- ass[index==j, pheno]
           var2 <- ass[index==j, metabolite]
           
           # get ids of pheno + covariates in the cohort
           sample.size.pheno <- pheno[, na.omit(.SD),.SDcols=c("sampleID", var1, var3)]$sampleID
           
           # get the ids of the metabolite
           sample.size.metab <- metab[, na.omit(.SD),.SDcols=c("id", var2)]$id
           
           # get the overlap
           v1 <- venn2(sample.size.pheno, sample.size.metab, plotte = FALSE)
           
           # annotate in result table
           ass.cohort[, cohort := (i)]
           ass.cohort[index == j, N := length(v1$q1)]
           # ass.cohort[index == j, ids := list(v1$q1)]
         })
         
         # concatenate
         res.cohort <- rbindlist(res.cohort)
         return(res.cohort)
         
       })

# use these IDs for the individual associations
res.ids <- rbindlist(res)

# Loop through each association ----
#' # Loop through each association

# loop through each cohort
res <- lapply(
  list(
    "Adult",
    "AMI",
    "Heart",
    "Sorb"), function(i){
      
      # Pheno ~ Gx ----
      
      # put gx into the right order
      gx <- as.data.frame(exprs(dat.list[[i]][["gx.dat"]]))
      setDT(gx,keep.rownames = TRUE)
      # hh(gx)
      gx <- dcast(melt(gx, id.vars = "rn",
                       variable.name = "id",
                       value.name = "measurement"),
                  id~rn, # transpose here
                  value.var = "measurement")
      # hh(gx)
      pheno <- as.data.table(pData(dat.list[[i]][["gx.dat"]]))
      metab <- dat.list[[i]][["metab.dat"]]
      covars <- dat.list[[i]][["covariates"]]
      
      # get the correct sample overlap for each pheno ~ metab + covariate combination
      # res.ids[cohort==(i), ]
      
      # put everything into the right order
      pheno.glm <- pheno[,.SD,.SDcols=c("sampleID", ap)]
      confounders <- pheno[,.SD,.SDcols=c("sampleID", covars)]
      setnames(pheno.glm,"sampleID","id")
      setnames(confounders,"sampleID","id")
      
      # fit models for binomial phenotypes
      res.pheno.gx.glm <- lapply(ap, function(j){
        
        res <- glm_loop(
          y = pheno.glm[,.SD,.SDcols=c("id",j)],
          x = gx,
          c = confounders,
          type = ifelse(j=="diabetes.status.tri", "binomial", "gaussian"),
          id = "id"
        )
        return(res)
        
      })
      
      # RESULT: these are all pheno~gx+conf associations
      res.pheno.gx.glm <- rbindlist(res.pheno.gx.glm)
      res.pheno.gx.glm[, study := i]
      res.pheno.gx.glm[, type:="pheno~gx"]
      
      # provide some feedback
      message(paste0("Pheno~Gx associations done in ", i))
      
      # Pheno ~ Metab ----
      
      # fit models for both phenotypes
      res.pheno.metab.glm <- lapply(ap, function(j){
        
        res <- glm_loop(
          y = pheno.glm[,.SD,.SDcols=c("id",j)],
          x = metab,
          c = confounders,
          type = ifelse(j=="diabetes.status.tri", "binomial", "gaussian"),
          id = "id"
        )
        return(res)
        
      })
      
      # RESULT: these are all pheno~metab+conf associations
      res.pheno.metab.glm <- rbindlist(res.pheno.metab.glm)
      res.pheno.metab.glm[, study := i]
      res.pheno.metab.glm[, type:="pheno~metab"]
      
      # provide some feedback
      message(paste0("Pheno~Metab associations done in ", i))
      
      # Gx ~ Metab ----
      
      #' Important note: Only use the probes available in this study! Some are only avaiable in (at least) two studies!
      
      res.gx.metab.glm <- lapply(unique(ptt$markerID)[unique(ptt$markerID)%in%names(gx)], function(j){
        
        # get the metabolites for this marker
        each.x <- ptt[markerID == (j), metabolite]
        
        # loop through the metabolites of this marker
        res <- glm_loop(
          y = gx[,.SD,.SDcols=c("id",j)],
          x = metab[, .SD,.SDcols=c("id", each.x)],
          c = confounders,
          type = "gaussian",
          id = "id"
        )
        return(res)
        
      })
      
      # RESULT: these are all gx~metab+conf associations
      res.gx.metab.glm <- rbindlist(res.gx.metab.glm)
      res.gx.metab.glm[, study := i]
      res.gx.metab.glm[, type:="gx~metab"]
      
      # provide some feedback
      message(paste0("Gx~Metab associations done in ", i))
      
      # Metab ~ Gx ----
      
      res.metab.gx.glm <- lapply(unique(ptt$metabolite), function(j){
        
        # get the metabolites for this marker
        each.x <- ptt[metabolite == (j), markerID]
        x.in.study <- each.x[each.x %in% names(gx)]
        
        # loop through the metabolites of this marker
        res <- glm_loop(
          y = metab[, .SD,.SDcols=c("id", j)],
          x = gx[,.SD,.SDcols=c("id",x.in.study)],
          c = confounders,
          type = "gaussian",
          id = "id"
        )
        return(res)
        
      })
      
      # RESULT: these are all gx~metab+conf associations
      res.metab.gx.glm <- rbindlist(res.metab.gx.glm)
      res.metab.gx.glm[, study := i]
      res.metab.gx.glm[, type:="metab~gx"]
      
      # provide some feedback
      message(paste0("Metab~Gx associations done in ", i))
      
      # return all the results
      res <- rbindlist(
        list(
          res.pheno.gx.glm,
          res.pheno.metab.glm,
          res.gx.metab.glm,
          res.metab.gx.glm
        )
      )
      
      return(res)
      
    # },mc.cores = 4,mc.cleanup = TRUE)
    })

# these are the results from all studies
res <- rbindlist(res)

fwrite(res,dpx("mediationParts1","res/"),sep = "\t")
# res <- fread(newest_file("mediationParts1","res",print_full = TRUE))

# CHECKS ----

# check wheter gx-metab and metab-gx are reciprocal assocs by comparing t-statistic
dat1 <- res[type=="gx~metab", .(id = paste(x,y,study), gx.metab = test.statistic)]
dat2 <- res[type=="metab~gx", .(id = paste(y,x,study), metab.gx = test.statistic)]
v1 <- venn2(dat1$id,dat2$id)
dat3 <- cbind(dat1[match(v1$q1, id), .(gx.metab)], dat2[match(v1$q1, id), .(metab.gx)])
plot(dat3$gx.metab,dat3$metab.gx, main="comparison of t-statistics gx-metab and metab-gx assoc")
abline(0,1,col="red",lty="dashed")

# Meta Analysis of single associations ----
#' # Meta Analysis of single associations

#' This loop takes every sub-type for the mediation analysis (e.g. pheno~gx, metab~gx, etc.) and runs the meta-analysis.

# test plan
test.plan <- res[,.(pheno=unique(y)), by = type]
test.plan[, index:=1:.N]

# loop through each thing
res.meta <- mclapply(unique(test.plan$type), function(i){
  
  # get the subset of the data
  # tmptmptmp <- dat
  dat <- res[type==(i), ]
  
  # set names as used in the function
  setnames(dat, 
           old = c("x", "y", "beta.se", "p.value", "study"),
           new = c("markerID", "var", "se", "p", "cohort"),skip_absent = T
  )
  
  # dummy columns for the meta-analysis
  dat[, `:=`(
    codedAll = "A",
    noncodedAll = "B",
    eaf = 0.2,
    imputed = FALSE,
    infoscore = 0.9
  )]
  
  # get the subset of the test plan
  test.plan.sub <- test.plan[type== (i), ]
  
  # run the meta-analysis for each phenotype in the type of analysis
  res.type <- lapply(test.plan.sub$pheno, function(j){
    
    # run the meta-analysis
    res.meta <- meta_analysis(
      phenotype = j,
      assocData = dat
    )
    return(res.meta)
  })
  
  # these are the results for a single type
  res.type <- rbindlist(res.type,use.names = T,fill = T)
  
  # remove columns that are of no use
  colfilt <- grepl(pattern = "nWeightedMAF|minMAF|nWeightedInfoScore|minInfoScore|eaf.|maf.|imputed.|infoscore.",
                   x = names(res.type))
  res.type <- res.type[,.SD,.SDcols=names(res.type)[!colfilt]]
  
  # set back to correct names
  setnames(dat, 
           old = c("markerID", "var", "se", "p", "cohort"),
           new = c("x", "y", "beta.se", "p.value", "study")
  )
  setnames(res.type,
           old=c("metabolite","markerID"),
           new=c("y","x"))
  
  # loop through each cohort and the fixed and random effects result
  cohort.cols <- paste0(".", unique(dat$study))
  for(k in c(cohort.cols, "FEM", "REM")){
    
    pcol <- paste0("p", k)
    bcol <- paste0("beta", k)
    secol <- paste0("se", k)
    
    # calculate the confidence Intervals/hierarchical fdr
    set(x = res.type, j = paste0("ci.l", k), value = res.type[[bcol]] - 1.96 * res.type[[secol]])
    set(x = res.type, j = paste0("ci.u", k), value = res.type[[bcol]] + 1.96 * res.type[[secol]])
    
    # add hierarchical fdr
    hfdr <- addHierarchFDR(pvalues = res.type[!is.na(get(pcol))][[pcol]],
                           # categs = res.type[!is.na(get(pcol)), y],
                           # "pheno~gx"    "pheno~metab" "gx~metab"    "metab~gx" 
                           categs = if(i == "pheno~metab" | i == "pheno~gx"){
                             res.type[!is.na(get(pcol)),x]} else {
                               res.type[!is.na(get(pcol)),y]},
                           fdrmethod_level1 = "BH",
                           fdrmethod_level2 = "BH",
                           correctionLevel1 = "BB")
    
    # filter for NA in pcol T/F
    na.filt <- which(!is.na(res.type[[pcol]]))
    set(x = res.type, i = na.filt, j = paste0("hfdr.1", k), value = hfdr$fdr_level1)
    set(x = res.type, i = na.filt, j = paste0("hfdr.2", k), value = hfdr$fdr_level2)
    set(x = res.type, i = na.filt, j = paste0("hfdr.sig", k), value = hfdr$hierarch_fdr5proz)
  }
  
  # RESULT
  res.type <- move_col_front(res.type, "y")
  res.type[ , type := i]
  res.type <- move_col_front(res.type, "type")
  return(res.type)
  
},mc.cores = 4,mc.cleanup = TRUE)

# RESULTS
res.meta <- rbindlist(res.meta)

# save 
fwrite(res.meta,dpx("mediationMetaParts1","res/"),sep = "\t")
# res.meta <- fread(newest_file("mediationMetaParts1","res",print_full = T))

# CHECKS ----

#' Compare t-statistics of gx~metab and metab~gx based on fixed and random effects model

tmp1 <- res.meta[type=="gx~metab",.(id=paste(x,y,sep="__"), tREM = betaREM/seREM, tFEM=betaFEM/seFEM)]
tmp2 <- res.meta[type=="metab~gx",.(id=paste(y,x,sep="__"), tREM = betaREM/seREM, tFEM=betaFEM/seFEM)]
v2 <- venn2(tmp1$id,tmp2$id)

tiff(dpx("metaTStatREMFEMComparison.tiff","plt/"),res=300,unit="in",height=5,width = 10)
par(mfrow=c(1,2))
plot(
  tmp1[match(v2$q1,id),tREM],
  tmp2[match(v2$q1,id),tREM],
  pch=19,
  col=alpha("grey40",0.5),
  main="REM gx~metab, metab~gx comparison"
)
abline(0,1,lwd=2,col="red",lty="dashed")
plot(
  tmp1[match(v2$q1,id),tFEM],
  tmp2[match(v2$q1,id),tFEM],
  pch=19,
  col=alpha("grey40",0.5),
  main="FEM gx~metab, metab~gx comparison"
)
abline(0,1,lwd=2,col="red",lty="dashed")
par(mfrow=c(1,1))
dev.off()

# Select clinically relevant associations ----
#' # Select clinically relevant associations

# clinically relevant gx-metab:
ptt %>% hh

uniqueN(ptt$metabolite)
uniqueN(ptt$markerID)

# gx associating with phenotype
res.meta[type== "pheno~gx", ] %>% hh
res.meta[type== "pheno~gx", .N,by=y]
res.meta[type== "pheno~gx" & numberStudies >=2, .N,by=y]
res.meta[type== "pheno~gx" & hfdr.sigREM == T & numberStudies >=2, .N,by=y]
res.meta[type== "pheno~gx" & numberStudies >=2,sum(hfdr.sigREM == T)/.N,by=y]
res.meta[type== "pheno~gx", table(hfdr.sigREM),by=y]

# metab associating with phenotype
res.meta[type== "pheno~metab", ] %>% hh()
res.meta[type== "pheno~metab", .N,by=y]
res.meta[type== "pheno~metab" & numberStudies >=2, .N,by=y]
res.meta[type== "pheno~metab" & hfdr.sigREM == T & numberStudies >=2, .N,by=y]
res.meta[type== "pheno~metab" & numberStudies >=2,sum(hfdr.sigREM == T)/.N,by=y]
res.meta[type== "pheno~metab", table(hfdr.sigREM),by=y]

# get significant gx
significant.gx <- res.meta[type== "pheno~gx" & hfdr.sigREM == T & numberStudies >=2, x, by=y]
v.gx <- venn2(significant.gx[y=="log.bmi",x], significant.gx[y=="diabetes.status.tri",x])
length(v.gx$q1)/res.meta[type== "pheno~gx", uniqueN(x)] # PAPER

# get significant metab
significant.metab <- res.meta[type== "pheno~metab" & hfdr.sigREM == T & numberStudies >=2, x, by=y]
v.metab <- venn2(significant.metab[y=="log.bmi",x], significant.metab[y=="diabetes.status.tri",x])

#' This means I'll only get the metabolite-gx pairs where either the metabolite OR the gx probe associate significantly with the phenotype. I already have the gx-metab associations

# Pheno ~ Gx + Metab Associations ----
#' ## Pheno ~ Gx + Metab Associations

#' I need to extract all the eligible combinations where either gx OR metab associates with any of the phenotypes

# amt = all mediation tests
amt <- lapply(ap, function(ii){
  
  sme <- significant.metab[y==(ii), x]
  sgx <- significant.gx[y==(ii), x]
  
  by.me <- ptt[metabolite %in% sme, ]
  by.gx <- ptt[markerID %in% sgx, ]
  
  # add info on direction, because I only want to mediate pairs with a sig. raw effect
  by.me[,exposure:="metab"]
  by.gx[,exposure:="gx"]
  
  # get overlapping pairs
  v <- venn2(by.me[,paste0(markerID,metabolite)],
             by.gx[,paste0(markerID,metabolite)],
             plotte = F)
  res <- rbindlist(list(by.gx,by.me))
  res[, pheno := (ii)]
  
  # mark pairs to be tested in both directions
  res[paste0(markerID,metabolite) %in% v$q1, exposure:="both"]
  res <- res[!duplicated(paste0(markerID,metabolite))]
  
}) %>% rbindlist
setnames(amt, c("metabolite", "gx.probe", "exposure", "pheno"))

# PAPER
table(amt[pheno=="log.bmi",exposure])
table(amt[pheno=="diabetes.status.tri",exposure])
amt[,.N,by=pheno]
amt[,.N,by=exposure]
amt[,id:=paste(metabolite,gx.probe,sep = "__")]
v.med <- venn2(amt[pheno=="diabetes.status.tri",id],amt[pheno=="log.bmi", id])

# save for use in part 2
save_csv_carl(file = amt,file_name = "allMediationsList","res")

# loop through each cohort
res <- mclapply(
  list(
    "Adult",
    "AMI",
    "Heart",
    "Sorb"), function(i){
      
      # Pheno ~ Gx + Metab ----
      
      # put gx into the right order
      gx <- as.data.frame(exprs(dat.list[[i]][["gx.dat"]]))
      setDT(gx,keep.rownames = TRUE)
      hh(gx)
      gx <- dcast(melt(gx, id.vars = "rn",
                       variable.name = "id",
                       value.name = "measurement"),
                  id~rn, # transpose here
                  value.var = "measurement")
      hh(gx)
      pheno <- as.data.table(pData(dat.list[[i]][["gx.dat"]]), keep.rownames = T)
      setnames(pheno,"rn","sampleID",skip_absent = TRUE)
      metab <- dat.list[[i]][["metab.dat"]]
      covars <- dat.list[[i]][["covariates"]]
      
      # put everything into the right order
      pheno.glm <- pheno[,.SD,.SDcols=c("sampleID", ap)]
      confounders <- pheno[,.SD,.SDcols=c("sampleID", covars)]
      setnames(pheno.glm,"sampleID","id")
      setnames(confounders,"sampleID","id")
      
      # get the correct sample overlap for each pheno ~ metab + covariate combination
      # res.ids[cohort==(i), ]
      
      # fit models for binomial phenotypes
      res.mediation.glm <- lapply(ap, function(j){
        
        # select the correct gx
        
        res <- glm_mediation_loop(
          y = pheno.glm[,.SD,.SDcols=c("id",j)],
          x1 = gx,
          x2 = metab,
          testPlan = amt[pheno==(j) & 
                           gx.probe %in% names(gx)[-1], 
                         ], # IMPORTANT: Only test available markers per study!
          c = confounders,
          type = ifelse(j=="diabetes.status.tri", "binomial", "gaussian"),
          id = "id"
        )
        return(res)
        
      })
      
      # RESULT: these are all pheno~gx+conf associations
      res.mediation.glm <- rbindlist(res.mediation.glm)
      res.mediation.glm[, study := i]
      res.mediation.glm[, type:="pheno~gx+Metab"]
      
      # provide some feedback
      message(paste0("Pheno~Gx+Metab associations done in ", i))
      
      return(res.mediation.glm)
      
    },mc.cores = 4,mc.cleanup = TRUE)

# these are the results from all studies
res <- rbindlist(res)

fwrite(res,dpx("mediationParts2","res/"),sep = "\t")
# res <- fread(newest_file(look_for = "mediationParts2",subfolder = "res",print_full = TRUE))
#############################################=

# Meta Analysis of mediation associations ----
#' # Meta Analysis of mediation associations

#' This loop takes the adjusted betas of pheno ~ gx + metab for the mediation triangles and runs the meta-analysis. That means, that the effects of the metabolites  adjusted for the gx.probe (and vice versa) on each phenotype are to be meta-analysed. 

# loop through each thing
res.meta <- lapply(unique(res$x.type), function(i){
  
  # the data to be meta-analysed
  dat <- res[x.type==(i), ]
  
  # set names as used in the function
  setnames(dat, 
           old = c("id.name", "y", "beta.se", "p.value", "study"),
           new = c("markerID", "var", "se", "p", "cohort"),
           skip_absent = T
  )
  
  # dummy columns for the meta-analysis
  dat[, `:=`(
    codedAll = "A",
    noncodedAll = "B",
    eaf = 0.2,
    imputed = FALSE,
    infoscore = 0.9
  )]
  
  # run the meta-analysis for each phenotype in the type of analysis
  res.type <- lapply(unique(dat$var), function(j){
    
    # run the meta-analysis
    res.meta <- meta_mediation_analysis(
      phenotype = j,
      assocData = dat
    )
    return(res.meta)
  })
  
  # these are the results for a single type
  res.type <- rbindlist(res.type,use.names = T,fill = T)
  
  # remove columns that are of no use
  colfilt <- grepl(pattern = "nWeightedMAF|minMAF|nWeightedInfoScore|minInfoScore|eaf.|maf.|imputed.|infoscore.",
                   x = names(res.type))
  res.type <- res.type[,.SD,.SDcols=names(res.type)[!colfilt]]
  
  # set back to correct names
  setnames(dat, 
           old = c("markerID", "var", "se", "p", "cohort"),
           new = c("id.name", "y", "beta.se", "p.value", "study")
  )
  setnames(res.type,
           old=c("metabolite","markerID"),
           new=c("y","id.name"))
  
  # loop through each cohort and the fixed and random effects result
  cohort.cols <- paste0(".", unique(dat$study))
  for(k in c(cohort.cols, "FEM", "REM")){
    
    pcol <- paste0("p", k)
    bcol <- paste0("beta", k)
    secol <- paste0("se", k)
    
    # calculate the confidence Intervals/hierarchical fdr
    set(x = res.type, j = paste0("ci.l", k), value = res.type[[bcol]] - 1.96 * res.type[[secol]])
    set(x = res.type, j = paste0("ci.u", k), value = res.type[[bcol]] + 1.96 * res.type[[secol]])
    
    # add hierarchical fdr
    hfdr <- addHierarchFDR(pvalues = res.type[!is.na(get(pcol))][[pcol]],
                           categs = res.type[!is.na(get(pcol)), x],
                           fdrmethod_level1 = "BH",
                           fdrmethod_level2 = "BH",
                           correctionLevel1 = "BB")
    
    # filter for NA in pcol T/F
    na.filt <- which(!is.na(res.type[[pcol]]))
    set(x = res.type, i = na.filt, j = paste0("hfdr.1", k), value = hfdr$fdr_level1)
    set(x = res.type, i = na.filt, j = paste0("hfdr.2", k), value = hfdr$fdr_level2)
    set(x = res.type, i = na.filt, j = paste0("hfdr.sig", k), value = hfdr$hierarch_fdr5proz)
  }
  
  # RESULT
  res.type <- move_col_front(res.type, "y")
  res.type[ , type := i]
  res.type <- move_col_front(res.type, "type")
  return(res.type)
  
})

# RESULTS
res.meta <- rbindlist(res.meta)

# remove column name inconsistency with the Part1 mediation results, where "type" indicates the association (e.g. gx->metab)
res.meta[, mediation.target := type]
res.meta[, type := "pheno~gx+metab"]
res.meta <- move_col_front(res.meta, "mediation.target")
res.meta <- move_col_front(res.meta, "type")

# extract and save interaction seperately
res.interaction <- res.meta[mediation.target=="interaction"]

# no sig. interaction
res.interaction[,table(hfdr.sigREM)]

# save 
fwrite(res.interaction,dpx("mediationMetaInteraction","res/"),sep = "\t")
fwrite(res.meta[mediation.target!="interaction"],dpx("mediationMetaParts2","res/"),sep = "\t")

# Wrap-Up ----
#' # Wrap-Up

sessionInfo()