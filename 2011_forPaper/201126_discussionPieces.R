#' ---
#' title: "Paper Discussion content"
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
  if(grepl(x= getwd(),pattern =  "/home/carl")){
    bp <- "/home/carl/imise/projekte/genstat/"
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
  "ggplot2",
  "readxl"
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
# source(here::here("../functions/all_annotation_tables.R"))
source(here::here("../functions/option_setup.R"))
source(here::here("../functions/theme_gxMetab.R"))
setDTthreads(8)
option_setup()

#+ include=F
#==========#
#   DONE   #
#==========#

# Load data ----
#' # Load data

# results from the main meta-analysis
nf <- newest_file(look_for = "STEP4_meta_resultsAnnotated$",subfolder = "res", print_full = TRUE)
res.init.meta <- fread(nf)

# load the joint mediation results
res.joint <- fread(newest_file("mediationStatisticsAnnotatedJoint", "res",print_full = TRUE))

# annotated sobel results
res.sobel <- fread(newest_file("mediationStatisticsAnnotated$","res",print_full = T))
r2.mediation <- fread(newest_file("mediationStatisticsMatchingMarginalInfo", "res/",print_full = T))

# metabolite names
load("/net/ifs1/san_projekte/projekte/genstat/02_projekte/1703_ge_metab_a1_b3_sorbs/1805_gxMetaboliteAssociation/obj/200306_metaboliteNames.RData")

# merge long names to results

m1 <- match(res.init.meta$metabolite, agn.all$rid)
res.init.meta$metab_abbr <- agn.all[m1, abbr]
res.init.meta$metab_long <- agn.all[m1, full.name]

# paper discussion

# gene aus Burkhardt/Kirsten 2015 ----
res.init.meta[symbol_INGENUITY == "ARG1" & significant == T,]
res.init.meta[symbol_INGENUITY == "SLC22A16" & significant == T]
res.init.meta[symbol_INGENUITY %in% c("ANUBL1", "FAM21C") & significant == T, .(metabolite, betaREM)]
res.init.meta[symbol_INGENUITY %in% c("HLCS") & significant == T, .(metabolite, betaREM)]
res.init.meta[symbol_INGENUITY %in% c("PPP1R16A","LRRC14") & significant == T, .(metabolite, betaREM)]

# auch explizit genannt
res.init.meta[symbol_INGENUITY %in% c("JAM3") & significant == T, .(metabolite, betaREM)]

res.init.meta[symbol_INGENUITY %in% c("MCCC1", "ETFDH", "SLC22A4","SLC22A5", "ACADM", "ACADS", "CPS1", "CRAT") & significant == T, .(metabolite, betaREM, symbol_INGENUITY)]

res.init.meta[significant==T & grepl("SLC22A", symbol_INGENUITY), table(symbol_INGENUITY)]

# speziell genannt in disc von burkhardt - C8 & C10
res.init.meta[symbol_INGENUITY %in% c("ACADM") & significant == T, .(metabolite, betaREM, symbol_INGENUITY)]

# burkhardt kirsten dat
bkd <- read_xlsx("dat/journal.pgen.1005510.s014.XLSX")

# Other Paper ----
# genannt in https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4745136/

v1 <- venn2(
  c("SLC22A16","SLC16A10","ARG1","AGPS","CENPU","ACSL1"),
  unique(res.init.meta$symbol_INGENUITY)
)
v2 <- venn2(
  c("SLC22A16","SLC16A10","ARG1","AGPS","CENPU","ACSL1"),
  unique(res.init.meta[significant==T,]$symbol_INGENUITY)
)

# how many assocs with which metabolite?
res.init.meta[symbol_INGENUITY %in% c("SLC22A16","SLC16A10","ARG1","AGPS","CENPU","ACSL1") & significant == T, .(metabolite, betaREM, symbol_INGENUITY)]

# Summary of results ----

# first paragraph

# summary for total assocs and mediations
res.init.meta[,sum(significant)]
res.init.meta[significant==T & symbol_INGENUITY != "",.(uniqueN(metabolite),uniqueN(symbol_INGENUITY))]

relative.mediation[pheno=="log.bmi" & mediation.group != "None", uniqueN(paste0(gx.probe, metabolite))]
relative.mediation[pheno=="log.bmi" & mediation.group != "None", uniqueN(gx.probe)]
relative.mediation[pheno=="log.bmi" & mediation.group != "None", uniqueN(metabolite)]
relative.mediation[pheno=="log.bmi" & mediation.group != "None" & gene != "", uniqueN(gene)]

r2.mediation[pheno == "log.bmi" & mediation.sig==T, .N, type]
r2.mediation[pheno == "log.bmi" & mediation.sig==T, .N]
r2.mediation[pheno == "log.bmi" & mediation.sig==T, sum(duplicated(paste0(gx.probe,metabolite)))]

# Barthel Theis 2015 ----
# https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1005274#sec012

# barthel theis data
btd <- read_xlsx("dat/journal.pgen.1005274.s013.XLSX")
setDT(btd)

# barthel theiss mapping
btm <- fread("dat/170815_metabolitesMapping_BarthelTheiss.csv")
btm[corresponding_metabolite == "-", corresponding_metabolite := NA]
btm[corresponding_metabolite == "", corresponding_metabolite := NA]

# example Isoleucine
btd[Metabolite=="isoleucine", ]
res.init.meta[metabolite=="LeuIle" & symbol_INGENUITY %in% btd$mRNA & significant == T, .(metabolite, betaREM, symbol_INGENUITY)]


# overlap with associations
  v3 <- venn2(
  btd$mRNA,
  unique(res.init.meta$symbol_INGENUITY)
)

v4 <- venn2(
  btd$mRNA,
  unique(res.init.meta[significant==T,]$symbol_INGENUITY)
)

btd$mRNA %>% uniqueN
btd %>% nrow
res.init.meta[symbol_INGENUITY %in% btd$mRNA & significant == T, .(metabolite, betaREM, symbol_INGENUITY)]
res.init.meta[symbol_INGENUITY %in% btd$mRNA & significant == T, .(uniqueN(metabolite), uniqueN(symbol_INGENUITY))]
res.init.meta[symbol_INGENUITY %in% btd$mRNA & significant == T, .(uniqueN(symbol_INGENUITY)),metabolite][order(-V1)]
r2.mediation[mediation.sig==T & gene %in% v4$q1]

# bcaas
btd[Metabolite %in% c("valine", "leucine", "isoleucine"), unique(Metabolite)]
res.init.meta[significant==T & metabolite %in% c("Val", "LeuIle"), unique(metabolite)]

# isoleucine ABCG1
btd[Metabolite %in% c("valine", "leucine", "isoleucine") & mRNA == "ABCG1", ]
res.init.meta[significant==T & metabolite %in% c("Val", "LeuIle") &
                symbol_INGENUITY == "ABCG1", .(metabolite, symbol_INGENUITY, betaREM, pREM)]

# overlap in sig genes
v1 <- venn2(
  btd[Metabolite %in% c("valine", "leucine", "isoleucine"), unique(mRNA)],
  res.init.meta[symbol_INGENUITY != "" &
                  significant==T & 
                  metabolite %in% c("Val", "LeuIle"), unique(symbol_INGENUITY)]
)

# Mediation Summary ----

# Strongest PM ----

r2.mediation[, pm :=  (beta.mediation) / (beta.mediation + tauAdj)]
r2.mediation[, total.effect :=  (beta.mediation + tauAdj)]

r2.mediation[mediation.sig==TRUE & pheno == "log.bmi"][
  order(-pm),head(.SD,5),type,
  .SDcols=c("type",
            "gene",
            "metabolite",
            "beta.mediation",
            "ci.95.lower",
            "ci.95.upper",
            "pm",
            "total.effect"
            )
  ]

# sortet by mediation beta
r2.mediation[mediation.sig==TRUE & pheno == "log.bmi"][
  order(-abs(beta.mediation)),head(.SD,5),type,
  .SDcols=c("type",
            "gene",
            "metabolite",
            "beta.mediation",
            "ci.95.lower",
            "ci.95.upper",
            "pm",
            "total.effect"
            )
  ]

# C4/C6 -> ABCG1 -> BMI
res.joint[type=="pheno~metab" & pheno == "log.bmi" & metabolite %in% c("C4", "C6") , 
          .(metabolite, round(betaREM,3), round(seREM,3), hfdr.sigREM, beta.alpha.beta + tauAdj)]


res.joint[type=="pheno~metab" & pheno == "log.bmi" & metabolite %in% c("C4", "C6") , 
          .(metabolite, round(betaREM,3), round(seREM,3), hfdr.sigREM, beta.alpha.beta + tauAdj)]
agn.all[rid %in% c("C4OH", "C6")]


# Shin 2014 - An atlas of genetic influences on human blood metabolites ----

# supp 4 
shin6 <- read_xlsx("dat/41588_2014_BFng2982_MOESM50_ESM.xlsx", 
                   sheet = 7,
                   skip = 3)
setDT(shin6)

shin_loci <- unlist(strsplit(shin6$`Locus name`,split = "|", fixed = T))
shin_loci <- unique(shin_loci[shin_loci != "unknown"])

v5 <- venn2(
  shin_loci,
  unique(res.init.meta$symbol_INGENUITY)
)

v6 <- venn2(
  shin_loci,
  unique(res.init.meta[significant==T,]$symbol_INGENUITY)
)

# clean the list
shin_causal_genes <- unlist(strsplit(shin6$`Predicted causal gene`,split = "|", fixed = T))
shin_causal_genes <- unique(unlist(strsplit(shin_causal_genes,split = " ", fixed = T)))
shin_causal_genes <- gsub("[^a-zA-z0-9]","",x = shin_causal_genes)
shin_causal_genes <- shin_causal_genes[shin_causal_genes != "unknown"] %>% sort %>% unique

# check predicted causal genes
v5 <- venn2(
  shin_causal_genes,
  unique(res.init.meta$symbol_INGENUITY)
)

v6 <- venn2(
  shin_causal_genes,
  unique(res.init.meta[significant==T,]$symbol_INGENUITY)
)

uniqueN(c(shin_loci,shin_causal_genes))

v7 <- venn3(
  shin_loci,
  shin_causal_genes,
  unique(res.init.meta$symbol_INGENUITY)
)

v8 <- venn3(
  shin_loci,
  shin_causal_genes,
  unique(res.init.meta[significant==T,]$symbol_INGENUITY)
)

# assocs with shin genes
res.init.meta[significant==T & symbol_INGENUITY %in% c(shin_loci,shin_causal_genes),]

# check against bcaas
shin_bcaa <- c(
  "isoleucine/ X-11529",
  "leucine",
  "valine/ isovalerylcarnitine",
  "valine/ proline"
  )
shin6[
  `Most associated metabolite or ratio` %in% (shin_bcaa), 
  .(`Predicted causal gene`, `Most associated metabolite or ratio`)
  ]
shin_bcaa_genes <- shin6[`Most associated metabolite or ratio` %in% (shin_bcaa), `Predicted causal gene`]
res.init.meta[symbol_INGENUITY %in% shin_bcaa_genes, unique(symbol_INGENUITY)]
res.init.meta[significant==T & 
                metabolite %in% c("Val", "LeuIle") &
                symbol_INGENUITY %in% c(shin_bcaa_genes),
              .(symbol_INGENUITY, metabolite, betaREM, pREM)]

# other associations with these genes
res.init.meta[symbol_INGENUITY %in% shin_bcaa_genes & significant==T, 
              .(metabolite, betaREM,pREM, metab.super.path, metab.path),symbol_INGENUITY]

agn.all[rid %in% c("Q27","Q38", "MeHis"), full.name]

# Metabolite ratios with new genes ----

# Most of the objects for this are created in STEP5_[...]_primaryStatistic

base.plot.dat[order(-new.perc),] %>% head(10)

# # check q14, q15, q28

# quotientensinn:
# important ratio
res.init.meta[significant==TRUE & metabolite == "Q28"]
base.plot.dat[metabolite=="Q28"]

# check q28
res.init.meta[significant==TRUE & metabolite == "Q28"]


q1 <- res.init.meta[significant==TRUE & metabolite == "Gly", unique(symbol_INGENUITY)]
q2 <- res.init.meta[significant==TRUE & metabolite == "Ser", unique(symbol_INGENUITY)]
q3 <- unique(c(q1,q2))
q3 <- q3[q3 != ""]
res.init.meta[significant==TRUE & metabolite == "Q28", uniqueN(symbol_INGENUITY)]
res.init.meta[significant==TRUE & metabolite == "Q28", .(metabolite, symbol_INGENUITY, betaREM, pREM)]
res.init.meta[significant==TRUE & metabolite == "Q28" & symbol_INGENUITY %nin% (q3), 
              .(metabolite, symbol_INGENUITY, betaREM, pREM, ci.lREM,ci.uREM)]

# Association hubs ----

assoc_hubs <- res.init.meta[
  significant==T & 
    symbol_INGENUITY!="", 
  .(uniqueN(metabolite),
    N = .N,
    mclass = paste(
      unlist(
        table(
          metab.class
          )
        ), 
      collapse = ", "
      )
    ),
  .(symbol_INGENUITY, 
    types_INGENUITY)][
    order(-V1)][
      1:50,
      .(sig_metabs = V1, 
        N,
        mclass,
        symbol_INGENUITY, 
        types_INGENUITY
        )
      ]

v9 <- venn4(
  res.init.meta[significant==T & symbol_INGENUITY!="", uniqueN(metabolite),symbol_INGENUITY][order(-V1)][1:50,symbol_INGENUITY],
  shin_causal_genes,
  shin_loci,
  btd$mRNA, mylabels = c("meta", "shin1", "shin2", "theis")
)

# mcccc1
v9$q1
v9$q14

# discuss top five
assoc_hubs %>% head(10)

# assoc hub genes

# ccdc50
res.init.meta[significant==T & symbol_INGENUITY=="CCDC50",uniqueN(metabolite)]

res.init.meta[significant==T & symbol_INGENUITY=="CCDC50",
              .(metabolite, 
                symbol_INGENUITY, 
                betaREM,
                pREM,
                metab.path, 
                metab.super.path)
              ] %>% data.frame

# abcg1
res.init.meta[
  significant==T & symbol_INGENUITY=="ABCG1",
  .(metabolite, metab.class, metab.path, metab.super.path, betaREM, pREM)
  ][order(-abs(betaREM))] %>% as.data.frame


res.init.meta[significant==T & symbol_INGENUITY=="ABGC1",unique(metabolite)]

# pcca 
res.init.meta[significant==T & symbol_INGENUITY=="PCCA",unique(metabolite)]
res.init.meta[
  significant==T & symbol_INGENUITY=="PCCA",
  .(metabolite, metab_abbr, metab_long, metab.class, metab.path, metab.super.path, betaREM, pREM)
  ][order(-abs(betaREM))] %>% as.data.frame

# acadvl
res.init.meta[significant==T & symbol_INGENUITY=="ACADVL",uniqueN(metabolite)]
res.init.meta[
  significant==T & symbol_INGENUITY=="ACADVL",
  .(metabolite, metab_abbr, metab_long, metab.class, metab.path, metab.super.path, betaREM, pREM)
  ][order(-abs(betaREM))] %>% as.data.frame

tmp1 <- res.init.meta[
  significant==T & symbol_INGENUITY=="ACADVL",
  .(unique(metabolite),
    unique(metab.path), 
    unique(metab.super.path), 
    unique(metab.class),
    unique(metab.class.q)),by=metabolite]
tmp1$V4 %>% table %>% sort
tmp1$V3 %>% table %>% sort
tmp1$V2 %>% table %>% sort

# alas2
res.init.meta[significant==T & symbol_INGENUITY=="ALAS2",unique(metabolite)]
res.init.meta[
  significant==T & symbol_INGENUITY=="ALAS2",
  .(metabolite, metab_abbr, metab_long, metab.class, metab.path, metab.super.path, betaREM, pREM)
  ][order(-abs(betaREM))] %>% as.data.frame


# Mediation hubs ----

# mediation hub type gene
mhg <- r2.mediation[pheno=="log.bmi" & mediation.sig == TRUE, uniqueN(metabolite), .(gene)][order(-V1)]

# mediation hub gene
mhtg <- r2.mediation[pheno=="log.bmi" & mediation.sig == TRUE, uniqueN(metabolite), .(type, gene)][order(-V1)]

# look at results by mediation type 

# mediated genes
mhtg[type=="pheno~gx(+metab)"] %>% head(5)

# mediating genes
mhtg[type=="pheno~metab(+gx)"] %>% head(5)

# interesting genes
res.init.meta[symbol_INGENUITY=="NRCAM" & significant==T, unique(metabolite)]
res.init.meta[symbol_INGENUITY=="MCCC1" & significant==T, unique(metabolite)]

# Top PM Mediations ----

r2.mediation[, pm := beta.mediation / (beta.mediation + tauAdj)]

r2.mediation[pheno == "log.bmi" & mediation.sig == T][
  order(-pm), head(.SD,5),type, 
  .SDcols = c("type", "pheno", "gx.probe",
              "gene", "metabolite", "pm", 
              "beta.mediation", "raw.sig.group",
              "metab.abbr", "metab.long")
  ]

# Bi-directional mediations ----

r2.mediation[pheno=="log.bmi" & mediation.sig == TRUE, ][duplicated(paste0(gx.probe, metabolite)), ]

# BMI genes ----

# WWOX 
r2.mediation[pheno=="log.bmi" & mediation.sig == TRUE & gene == "WWOX", .(type, metabolite, pm, total.effect)] %>% as.data.frame
r2.mediation[pheno=="log.bmi" & mediation.sig == TRUE & gene == "WWOX", .(uniqueN(metabolite))]

# IRF4 
r2.mediation[pheno=="log.bmi" & mediation.sig == TRUE & gene == "IRF4", .(type, metabolite, pm, total.effect)] %>% as.data.frame
r2.mediation[pheno=="log.bmi" & mediation.sig == TRUE & gene == "IRF4", .(uniqueN(metabolite))]


# BMP6 
r2.mediation[pheno=="log.bmi" & mediation.sig == TRUE & gene == "BMP6", .(type, metabolite, pm, total.effect)] %>% as.data.frame
r2.mediation[pheno=="log.bmi" & mediation.sig == TRUE & gene == "BMP6", .(uniqueN(metabolite))]
agn.all[rid %in% c("Q19", "Q30")]

res.init.meta[, .(unique(metab.class), unique(metab.class.q)),metabolite]$V2 %>% table
res.init.meta[symbol_INGENUITY!="", uniqueN(symbol_INGENUITY)]
res.init.meta[symbol_INGENUITY!="", uniqueN(markerID)]
