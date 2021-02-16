#' ---
#' title: "Mediation Results"
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
#'      number_sections: FALSE
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
  "ggplot2"
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
source(here("../functions/all_annotation_tables.R"))
source(here("../functions/option_setup.R"))
an <- all_annotation_tables(mount.point = "/net/ifs2/san_projekte/projekte/genstat/")
theme_set(theme_light(base_size = 15, base_family = "Helvetica"))
setDTthreads(1)
option_setup()

#+ include=F
#==========#
#   DONE   #
#==========#

# Load data ----

#' # Load data

# joint association results
res.joint <- fread(newest_file("mediationStatisticsAnnotatedJoint", "res",print_full = TRUE))

# results from the main meta-analysis
nf <- newest_file(look_for = "STEP4_meta_resultsAnnotated$",subfolder = "res", print_full = TRUE)
res.init.meta <- fread(nf)

# mediation results
r2.mediation <- fread(newest_file("mediationStatisticsMatchingMarginalInfo", "res/",print_full = T))

# re-formatted mediation results
relative.mediation <- fread(newest_file("mediatonGroupList","res/",print_full = T))

# all associating pairs
amt <- fread(newest_file("allMediationsList","res",print_full = TRUE))

# proper metabolite names
load("/net/ifs1/san_projekte/projekte/genstat/02_projekte/1703_ge_metab_a1_b3_sorbs/1805_gxMetaboliteAssociation/obj/200306_metaboliteNames.RData")

# annotate genes in res.joint
m <- match(res.joint$gx.probe, res.init.meta$markerID)
res.joint$symbol_INGENUITY <- res.init.meta[m,symbol_INGENUITY]

# match proper names
r2.mediation$rid <- r2.mediation$metabolite
m <- match(r2.mediation$metabolite, agn.all$rid)
r2.mediation$metabolite <- agn.all[m,abbr]

res.joint$rid <- res.joint$metabolite
m <- match(res.joint$rid, agn.all$rid)
res.joint$metabolite <- agn.all[m, abbr]

relative.mediation$rid <- relative.mediation$metabolite
m <- match(relative.mediation$rid, agn.all$rid)
relative.mediation$metabolite <- agn.all[m, abbr]

# metab.class reference
mcref <- res.init.meta[, .(
  class = unique(metab.class), 
  class.q = unique(metab.class.q),
  m.sup = unique(metab.super.path),
  m.path = unique(metab.path)
  ),by=metabolite]

# merge to data
m <- match(r2.mediation$rid, mcref$metabolite)
r2.mediation$mclass <- mcref[m,class.q]
r2.mediation$mclassq <- mcref[m,class]
r2.mediation$mpath <- mcref[m,m.path]
r2.mediation$mspath <- mcref[m,m.sup]

#' # Create Info
# Create Info ----


# Paper: How much mediation? ----

# Paper: Who mediates? ----


empty.genes <- r2.mediation[gene=="" & mediation.sig==T,.N]
d <- r2.mediation[gene!="" & pheno=="log.bmi"]

#' # Text
# Text ----

#' # 4.3	Mediation of effects on BMI
#'
#' From the `r res.init.meta[significant==T,uniqueN(markerID)]` probes (mapped to `r res.init.meta[significant==T & symbol_INGENUITY!="",uniqueN(symbol_INGENUITY)]` genes) that associated significantly  (hierarchical FDR=5%) with `r res.init.meta[significant==T,uniqueN(metabolite)]` metabolites (`r res.init.meta[significant==T & symbol_INGENUITY!="",uniqueN(metabolite),by=metab.class.q][metab.class.q == "aa", V1]` amino acids, `r res.init.meta[significant==T & symbol_INGENUITY!="",uniqueN(metabolite),by=metab.class.q][metab.class.q == "ac", V1]` acylcarnitines and `r res.init.meta[significant==T & symbol_INGENUITY!="",uniqueN(metabolite),by=metab.class.q][metab.class.q == "quotient", V1]` quotients) in the meta-analysis of probe-metabolite associations from at least two single studies, `r res.joint[type=="pheno~gx" & pheno == "log.bmi" & hfdr.sigREM==TRUE, uniqueN(gx.probe)]` (`r round(res.joint[type=="pheno~gx" & pheno == "log.bmi" & hfdr.sigREM==TRUE, uniqueN(gx.probe)]/res.joint[type=="pheno~gx" & pheno == "log.bmi", uniqueN(gx.probe)]*100,1)`%) probes (`r res.joint[type=="pheno~gx" & pheno == "log.bmi" & hfdr.sigREM==TRUE, uniqueN(symbol_INGENUITY)]` genes; `r round(res.joint[type=="pheno~gx" & pheno == "log.bmi" & hfdr.sigREM==TRUE, uniqueN(symbol_INGENUITY)]/res.joint[type=="pheno~gx" & pheno == "log.bmi", uniqueN(symbol_INGENUITY)]*100,1)`%) and `r res.joint[type=="pheno~metab" & pheno == "log.bmi" & hfdr.sigREM==TRUE, uniqueN(metabolite)]` (`r round(res.joint[type=="pheno~metab" & pheno == "log.bmi" & hfdr.sigREM==TRUE, uniqueN(metabolite)]/res.joint[type=="pheno~metab" & pheno == "log.bmi", uniqueN(metabolite)]*100,1)`%) Metabolites associated directly with BMI at hierarchical FDR=5%. Mediation of effects on BMI was tested for each associating probe-metabolite pair where either metabolite or gene expression associated directly with BMI, with the associated feature as exposure, resulting in `r r2.mediation[pheno=="log.bmi",uniqueN(paste0(gx.probe,metabolite))]` gene expression probe-metabolite combinations of `r r2.mediation[pheno=="log.bmi", uniqueN(metabolite)]` metabolites and `r r2.mediation[pheno=="log.bmi", uniqueN(gx.probe)]` gene expression probes (`r r2.mediation[pheno=="log.bmi", uniqueN(gene)]` genes). Tests were conducted both with gene expression and metabolite as exposure when both associated with BMI. Thus, a total of `r r2.mediation[pheno=="log.bmi",.N]` mediations were tested.
#'
#' Hierarchical correction of p-values on global and family level (the exposure is considered the family) results in a total of`r r2.mediation[pheno=="log.bmi", sum(mediation.sig==T)]` significant mediations on BMI at hierarchical FDR=5%. Of those, `r r2.mediation[pheno=="log.bmi" & type=="pheno~gx(+metab)", sum(mediation.sig==T)]` (`r (r2.mediation[pheno=="log.bmi" & type=="pheno~gx(+metab)", sum(mediation.sig==T)] / r2.mediation[pheno=="log.bmi",.N]*100) %>% round(1)`% of the tested `r r2.mediation[pheno=="log.bmi",.N]`) are with  gene expression as exposures and `r r2.mediation[pheno=="log.bmi" & type=="pheno~metab(+gx)", sum(mediation.sig==T)]` (`r (r2.mediation[pheno=="log.bmi" & type=="pheno~metab(+gx)", sum(mediation.sig==T)] / r2.mediation[pheno=="log.bmi",.N]*100) %>% round(1)`% of the `r r2.mediation[pheno=="log.bmi",.N]` tested) metabolites. `r r2.mediation[pheno=="log.bmi" & mediation.sig == TRUE, sum(duplicated(paste0(gx.probe, metabolite)))]` (`r round(r2.mediation[pheno=="log.bmi" & mediation.sig == TRUE, sum(duplicated(paste0(gx.probe, metabolite)))] / r2.mediation[pheno=="log.bmi", uniqueN(paste0(gx.probe, metabolite))]*100, 1)`% of all `r r2.mediation[pheno=="log.bmi", uniqueN(paste0(gx.probe, metabolite))]` tested pairs) mediations of `r r2.mediation[pheno=="log.bmi" & mediation.sig == TRUE][duplicated(paste0(gx.probe, metabolite)), uniqueN(metabolite)]` metabolites and `r r2.mediation[pheno=="log.bmi" & mediation.sig == TRUE][duplicated(paste0(gx.probe, metabolite)), uniqueN(gx.probe)]` gene expression probes (`r r2.mediation[pheno=="log.bmi" & mediation.sig == TRUE & gene != ""][duplicated(paste0(gx.probe, metabolite)), uniqueN(gene)]` genes) on BMI were significant for both probe- and metabolite-mediated effects. An overview of available features and tested mediations is available in Figure X.
#'
# Plot: Barplot here ----
#'
#' ## 4.3.1	Linkage of mediated and mediating genes and metabolites
#'
#' Due easier interpretability of the results, only  mediations involving mapped genes were considered in the following. `r empty.genes` significant mediations involving non-mapped probes were removed. Overall, among the significant mediations on BMI, the five most frequently mediating/mediated metabolites were `r r2.mediation[mediation.sig==TRUE,.N,by=metabolite][order(-N)]$metabolite %>% head(5) %>% paste(collapse = ", ")`. `r d[mediation.sig==TRUE, uniqueN(gene), by=.(mclassq)][mclassq=="aa",V1]` genes act as either exposure or mediator for effects on BMI with `r d[mediation.sig==TRUE&mclassq=="aa",uniqueN(metabolite)]` amino acids (or ratios of amino acids) and `r d[mediation.sig==TRUE, uniqueN(gene), by=.(mclassq)][mclassq=="ac",V1]` genes with `r d[mediation.sig==TRUE&mclassq=="ac",uniqueN(metabolite)]` acylcarnitines (or ratios of acylcarnitines). `r d[mediation.sig==TRUE, uniqueN(gene), by=.(mclassq)][mclassq=="mix",V1]` genes are exposure or mediator for effects of `r d[mediation.sig==TRUE&mclassq=="mix",uniqueN(metabolite)]` mix quotient (ratios of amino acids and acylcarnitines), namely the ratio of `r agn.all[rid=="Q19", full.name]`, representing activity of the BCAA (#TAG-Abbreviation#) Isoleucine pathway. Conversely, the 5 most frequent mediated/medating genes were `r r2.mediation[mediation.sig==TRUE,.N,by=gene][order(-N)]$gene %>% head(5) %>% paste(collapse = ", ")`. Summarizing mediation results based on the mediator (either genes or metabolites) emphasizes the metabolites role in mediating effects of gene expression on BMI. `r r2.mediation[mediation.sig==TRUE,uniqueN(metabolite),by=.(type, mclassq)][type=="pheno~gx(+metab)" & mclassq=="aa",V1]` Amino acids were mediators for `r r2.mediation[mediation.sig==TRUE,uniqueN(gene),by=.(type, mclassq)][type=="pheno~gx(+metab)" & mclassq=="aa",V1]` genes, while effects of `r r2.mediation[mediation.sig==TRUE,uniqueN(metabolite),by=.(type, mclassq)][type=="pheno~metab(+gx)" & mclassq=="aa",V1]` amino acids were exposures for mediation by `r r2.mediation[mediation.sig==TRUE,uniqueN(gene),by=.(type, mclassq)][type=="pheno~metab(+gx)" & mclassq=="aa",V1]` genes. `r r2.mediation[mediation.sig==TRUE,uniqueN(metabolite),by=.(type, mclassq)][type=="pheno~gx(+metab)" & mclassq=="ac",V1]` acylcarnitines were mediators of effects of `r r2.mediation[mediation.sig==TRUE,uniqueN(gene),by=.(type, mclassq)][type=="pheno~gx(+metab)" & mclassq=="ac",V1]` genes, while effects of `r r2.mediation[mediation.sig==TRUE,uniqueN(metabolite),by=.(type, mclassq)][type=="pheno~metab(+gx)" & mclassq=="ac",V1]` acylcarnitines wer exposures in mediations by `r r2.mediation[mediation.sig==TRUE,uniqueN(gene),by=.(type, mclassq)][type=="pheno~metab(+gx)" & mclassq=="ac",V1]` genes. 
#' 
#' Acylcarnitines were mediators for up to `r d[mediation.sig==T & type=="pheno~gx(+metab)" & mclassq=="ac",uniqueN(gene),by=.(metabolite,mclassq)][order(-V1)][1,V1]` genes (`r d[mediation.sig==T & type=="pheno~gx(+metab)" & mclassq=="ac",uniqueN(gene),by=.(metabolite,mclassq)][order(-V1)][1,metabolite]`),  and exposures in mediations involving more unique genes than any amino acid (`r d[mediation.sig==T & type=="pheno~metab(+gx)" & mclassq=="ac",uniqueN(gene),by=.(metabolite,mclassq)][order(-V1)][1,paste(V1, metabolite, sep=" mediations: ")]`). The top five acylcarnitine mediators were `r d[mediation.sig==T & type=="pheno~gx(+metab)" & mclassq=="ac",uniqueN(gene),by=.(metabolite,mclassq)][order(-V1)][1:5,paste(metabolite, collapse=", ")]`, mediating effects of `r d[mediation.sig==T & type=="pheno~gx(+metab)" & mclassq=="ac",uniqueN(gene),by=.(metabolite,mclassq)][order(-V1)][1:5,paste(V1, collapse=", ")]` genes. The top five most frequent acylcarnitine exposures were `r d[mediation.sig==T & type=="pheno~metab(+gx)" & mclassq=="ac",uniqueN(gene),by=.(metabolite,mclassq)][order(-V1)][, head(metabolite,5)] %>% paste(collapse = ", ")`, being mediated by `r d[mediation.sig==T & type=="pheno~metab(+gx)" & mclassq=="ac",uniqueN(gene),by=.(metabolite,mclassq)][order(-V1)][, head(V1,5)] %>% paste(collapse = ", ")` unique genes, each. 

# top aa mediators
tmp1 <- d[mediation.sig==T & type=="pheno~gx(+metab)" & mclassq=="aa",uniqueN(gene),by=.(metabolite,mclassq)][order(-V1)]
tmp2 <- d[mediation.sig==T & type=="pheno~metab(+gx)" & mclassq=="aa",uniqueN(gene),by=.(metabolite,mclassq)][order(-V1)]
tmp3 <- d[mediation.sig==T,uniqueN(gene),by=.(type,metabolite,mclassq)][order(-V1)]

#' Top five amino acid mediators were `r tmp1$metabolite %>% head(5) %>% paste(collapse=", ")`, mediating effects of `r tmp1$V1 %>% head(5) %>% paste(collapse=", ")` unique genes. Reversely, `r tmp2$metabolite %>% head(5) %>% paste(collapse=", ")` were the top five exposures whose effects on BMI were mediated by the most unique genes (`r tmp2$V1 %>% head(5) %>% paste(collapse=", ")`). The top four acylcarnitines are mediators of more genes on BMI (from 435 to 298 genes) than the top amino acid, namely `r tmp3[mclassq=="aa", head(metabolite, 1)]`.

tmp1 <- d[mediation.sig==T,uniqueN(metabolite),by=.(gene,type)][order(-V1)]
tmp1[,type:=ifelse(type=="pheno~gx(+metab)", "gx","metab")]
tmp2 <- d[mediation.sig==T,uniqueN(metabolite),by=.(gene)][order(-V1)]

#' The genes `r tmp2$gene[1]` (being the exposure for mediation of `r tmp1[gene==tmp2$gene[1]&type=="gx",V1]` and mediator of effects of `r tmp1[gene==tmp2$gene[1]&type=="metab",V1]` metabolites) and `r tmp2$gene[2]` (exposure for `r tmp1[gene==tmp2$gene[2]&type=="gx",V1]` unique metabolite effects and mediator of `r tmp1[gene==tmp2$gene[2]&type=="metab",V1]` metabolite effects) being the most connected overall as well as mediating and mediated genes (Figure X). Pathways of metabolites and class of genes involved in mediation are displayed in Figure X. Mostly, mediating/mediated metabolites are generally involved in the metabolism of (branched-chain) amino acids and the energy/fatty acid metabolism. Most mediated/mediating genes were not attributed to any class. The next most frequent gene classes involved in mediation were enzymes, transcription regulators and transporters (Figure X).
#'
#' Plot: Treeplot here ----
#'

tmp1 <- res.init.meta[significant==TRUE & symbol_INGENUITY!="",uniqueN(metabolite),by=symbol_INGENUITY][order(-V1)]
tmp2 <- d[mediation.sig==T,uniqueN(metabolite),by=gene][order(-V1)]
d[,beta.reduct:= beta * alpha / (tauAdj + beta * alpha)]
d[mediation.sig==T,mean(beta.reduct),by=type][type=="pheno~metab(+gx)",V1]

d[mediation.sig==T,max(beta.reduct),by=type]
g1 <- d[mediation.sig==T&type=="pheno~gx(+metab)",beta.reduct]
g2 <- d[mediation.sig==T&type=="pheno~metab(+gx)",beta.reduct]
wilcox.test(x=g1,y=g2,paired = F)

#'
#' Genes that associated with the highest number of unique metabolites (top five: `r tmp1$symbol_INGENUITY[1:5] %>% paste(collapse=", ")`) in the association analysis were also involved in the mediations with the most metabolites (`r tmp2$gene[1:5] %>% paste(collapse=", ")`). Effects of gene expression on BMI were on average mediated to a higher extend by metabolites ($\mu_{1-\frac{\beta_{direct}}{\beta_{raw}}}=`r d[mediation.sig==T,mean(beta.reduct),by=type][type=="pheno~gx(+metab)",V1] %>% round(2)`$, $(1-\frac{\beta_{direct}}{\beta_{raw}})_{max}=`r d[mediation.sig==T,max(beta.reduct),by=type][type=="pheno~gx(+metab)",V1] %>% round(2)`$), as well as to a higher maximal proportion compared to gene-expression mediated effects ($\mu_{1-\frac{\beta_{direct}}{\beta_{raw}}}=`r d[mediation.sig==T,mean(beta.reduct),by=type][type=="pheno~metab(+gx)",V1] %>% round(2)`$, $(1-\frac{\beta_{direct}}{\beta_{raw}})_{max}=`r d[mediation.sig==T,max(beta.reduct),by=type][type=="pheno~metab(+gx)",V1] %>% round(2)`$; **Figure X**).
#'
# Plot: Histogram here ----
#'

# **HIER WEITER**

gx.top <- d[mediation.sig==T,][order(-abs(beta.mediation))][which(type == "pheno~gx(+metab)")[1],]
gx.top2 <- res.joint[type=="pheno~gx(+metab)" & 
            metabolite==gx.top$metabolite & 
            symbol_INGENUITY==gx.top$gene & 
            pheno=="log.bmi" & 
            gx.probe==gx.top$gx.probe,.(se,p.value)] %>% round(3)

# aa/ac distribution
d[mediation.sig==T,.N,mclassq]

# strongest effects not aa
d[mediation.sig==T,][order(-abs(beta.mediation))][,which(mclassq != "aa")[1]]

top.med <- res.joint[!is.na(beta.alpha.beta) & pheno=="log.bmi"][order(-beta.alpha.beta)] %>% head(5)

#'
#' Overall , mediation of metabolite effects by gene expression exhibit the highest overall mediation effect estimates. Sorting significant mediation results by their absolute effect estimate ($|\beta_{mediation}|$) shows that the `r d[mediation.sig==T,][order(-abs(beta.mediation))][,which(type == "pheno~gx(+metab)")[1]]` strongest mediations of metabolite effects exhibit larger effect estimates than the strongest mediation with a gene expression exposure (mediation of the effect of `r gx.top$gene` via `r gx.top$metabolite`; $\beta_{mediation}=`r gx.top$beta.mediation %>% round(3)`$, $SE_{mediation}=`r gx.top2$se`$, $p_{mediation}=`r gx.top2$p.value`$). While acylcarnitines are involved in mediations with the most genes, the highest mediation effects were estimated for mediations involving amino acids as mediators as well as mediated exposures (Table X). The strongest `r d[mediation.sig==T,][order(-abs(beta.mediation))][,which(mclassq != "aa")[1]]` mediations involve amino acids. Overall, the strongest mediation of effects were observed for the mediation of the effect of the quotient `r top.med[1, metabolite]` on BMI via the gene  `r top.med[1, symbol_INGENUITY]` ($\beta_{mediation}=`r  `$, $SE{mediation}=`r top.med[1, se] %>% round(4)`$, $p_{mediation}=`r top.med[1, p.value] %>% round(4)`$). Similarly strong mediations were observed for the metabolite `r top.med[2, metabolite]` and the gene  `r top.med[2, symbol_INGENUITY]`, `r top.med[3, metabolite]` and the genes `r top.med[3, symbol_INGENUITY]` and `r top.med[4, symbol_INGENUITY]` and `r top.med[5, metabolite]` and the gene `r top.med[5, symbol_INGENUITY]` (**Table X**).
#'
#' ## 4.3.2 Categorization of mediation results
#' --> Methods
#' To single out potentially interesting mediations, we categorized significantly (hierarchical FDR=5%) mediating pairs of genes and metabolites by the ratio of direct versus raw effect of the exposure on the outcome (i.e. the standardized regression coefficients $\beta$, see *Methods*). Mediations were classified as “gene expression mediated metabolite effects” when the mediation of metabolite effect on the outcome was significantly mediated by the gene expression probe, the ratio of the direct effect of the metabolite on the outcome (adjusted for the probe) and the raw effect of the metabolite of the outcome was below 0.8 and the ratio of the direct effect of the gene expression probe on the outcome (adjusted for the metabolite) and the raw effect of the gene expression on the outcome was between 0.8 and 1.25. Similarly, mediations were classified as “metabolite mediated gene expression effects” when the mediation of the gene expression effect on the outcome was significantly mediated by the metabolite, the ratio of the direct effect of the gene expression and the raw effect of the gene expression was below 0.8 and the direct effect of the ratio of the direct metabolite effect and the raw effect of the metabolite was between 0.8 and 1.25. “Bi-directional mediations” exhibit significant mediation effects for both metabolite- and gene expression-mediated effect and a ratio of direct and raw metabolite and gene expression effect on the outcome below 0.8. Mediations were classified as “other” when it was significant but the ratio of direct to raw effect did not meet the above specifications. Non-significant mediations were classified as “None”.
#'
#' Plot: Scatterplot ----
#'

relative.mediation[,.N,mediation.group][mediation.group=="None",N]

#' Thus, `r ifelse(length(relative.mediation[pheno=="log.bmi",.N,mediation.group][mediation.group=="Bi-directional mediation",N])==0,0,relative.mediation[pheno=="log.bmi",.N,mediation.group][mediation.group=="Bi-directional mediation",N])` mediations on BMI were classified as bi-directional, `r relative.mediation[pheno=="log.bmi",.N,mediation.group][mediation.group=="Gene expression mediated metabolite effects",N]` as gene expression-mediated, `r relative.mediation[pheno=="log.bmi",.N,mediation.group][mediation.group=="Metabolite mediated gene expression effects",N]` as metabolite-mediated and `r relative.mediation[pheno=="log.bmi",.N,mediation.group][mediation.group=="Other",N]` as “Other” (Figure X). **Table X** shows the genes and metabolites that are mediated/mediating in each of these groups. Using this grouping, the most frequent mediator genes were `r relative.mediation[pheno=="log.bmi" & mediation.group=="Gene expression mediated metabolite effects", .N,by=gene][order(-N),gene] %>% head(5) %>% paste(collapse=", ")`. The most freqent exposure genes were `r relative.mediation[pheno=="log.bmi" & mediation.group=="Metabolite mediated gene expression effects", .N,by=gene][order(-N),gene] %>% head(5) %>% paste(collapse=", ")`. Conversely, the most mediating metabolites were `r relative.mediation[pheno=="log.bmi" & mediation.group=="Metabolite mediated gene expression effects", .N,by=metabolite][order(-N),metabolite] %>% head(5) %>% paste(collapse=", ")` and the most mediated metabolites were `r relative.mediation[pheno=="log.bmi" & mediation.group=="Gene expression mediated metabolite effects", .N,by=metabolite][order(-N),metabolite] %>% head(5) %>% paste(collapse=", ")`.
#'
#' ## 4.3.3	Mediation of obesity-related genes

# vorher bmi.genes & d2f aus 3.5 laden!
bmi.genes <- c("FTO", "MC4R", "TMEM18", "KCTD15", "GNPDA2", "SH2B1", "MTCH2",
               "NEGR1", "BDNF", "OLFM4", "ADCY3", "PFKP", "NPC1", "MAF", "PTER",
               "NRXN3", "TNN13K", "SEC16B", "POMC", "PRKD1", "CADM2", "QPCTL",
               "FANCL", "EFEMP1", "BMP6", "MIR-129-2", "HSD17B12", "PRDM11",
               "WWOX", "KCNJ2", "GSTCD", "PTCH1", "LINC01122", "NLRC3-ADCY9",
               "GPRC5B-GP2", "AGBL4-ELAVL4", "ATP2A1-SBK1", "TCF7L2", "GIPR",
               "IRS1", "FOXO3", "ASB4", "RPTOR", "CREB1", "FAM57B", "APOBR",
               "PTBP2", "ELAVL4", "CELF1", "RALYL", "MAP2K5", "MAPK3", "FAIM2",
               "PARK2", "IRX3", "IRX5", "ARID5B",
               # thermogenesis genes
               "IRF4", "PGC1α", "PRDM16", "TBX15", "UCP1",
               # non-syndromic monogenetic variants
               "LEPR", "LEP", "PCSK1", "NTRK2", "SIM1","KSR2", "TUB","MRAP2","AMY1",
               # Inflammatory genes hinzufügen?
               "PBRM1",
               # Lotta 2016
               "PPM1K"
)

bmi.genes <- res.init.meta[symbol_INGENUITY %in% bmi.genes, unique(symbol_INGENUITY)]
bmi.genes %>% length
sig.bmi.genes<-res.init.meta[symbol_INGENUITY %in% bmi.genes & significant==T, unique(symbol_INGENUITY)]
sig.bmi.genes %>% length


bmi.genes %>% length
sig.bmi.genes %>% length
d2f <- fread("obj/d2f.csv")
d2f[,unique(metabolite),by=m.class1][,table(m.class1)]
tmp <- d2f[,.(.N, m=paste(metabolite,collapse = ", ")),by=gene]


#'
#' We explore the range of mediated effects by investigating a variety of obesity-associated genes. Categorized mediation results on obesity were filtered for mediations of/with `r bmi.genes %>% length` obesity-associated genes reported in European populations (Table X). Of those, `r sig.bmi.genes %>% length` associated with at least one metabolite in the meta-analysis. We found a total of `r d2f %>% nrow` mediations of gene expression effects of `r d2f$gene %>% uniqueN` (`r d2f$gene %>% unique %>% paste(collapse=", ")`) obesity-associated genes by `r d2f$metabolite %>% uniqueN()` metabolites (`r d2f[,unique(metabolite),by=.(m.class1,m.class2)][m.class1=="aa",.N]` amino acids including `r d2f[,unique(metabolite),by=.(m.class1,m.class2)][m.class1=="aa" & m.class2=="quotient",.N]` ratios of amino acids and `r d2f[,unique(metabolite),by=.(m.class1,m.class2)][m.class1=="ac",.N]` acylcarnitines, including `r d2f[,unique(metabolite),by=.(m.class1,m.class2)][m.class1=="ac" & m.class2=="quotient",.N]` ratios of acylcarnitines and `r d2f[,unique(metabolite),by=.(m.class1,m.class2)][m.class1=="mix",.N]` mix quotient; **Table X**). 
#' 
#' The genes `r d2f[,.N,by=gene][N==1,gene] %>% paste(collapse=", ")` were involved in mediations with one metabolite each (`r paste(tmp[N==1]$m,collapse=", ")`, respectively), wile the genes `r d2f[,.N,by=gene][N>1,gene] %>% paste(collapse=", ")` were involved with `r paste(tmp[N>1]$N,collapse=", ")` metabolites, respectively. `r tmp[N>1][1,gene]` was involved in mediations with the metabolites `r tmp[N>1][1,m]`, `r tmp[N>1][2,gene]` with `r tmp[N>1][2,m]` and `r tmp[N>1][3,gene]` with `r tmp[N>1][3,m]` (**Figure X**).
#'
#' Plot: BMI genes network ----
#'
#' 