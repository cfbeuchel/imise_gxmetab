#' ---
#' title: "Network of the Mediation"
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
  bp <- "/net/ifs2/san_projekte/projekte/genstat/"
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
  "igraph",
  "magrittr",
  "ggplot2",
  "colorspace"
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
source(here::here("../functions/all_annotation_tables.R"))
source(here::here("../functions/option_setup.R"))
source(here::here("../functions/theme_gxMetab.R"))
an <- all_annotation_tables(mount.point = "/net/ifs2/san_projekte/projekte/genstat/")
theme_set(theme_light(base_size = 15, base_family = "Helvetica"))
setDTthreads(1)
option_setup()

#+ include=F
#==========#
#   DONE   #
#==========#

#=================================================================================#
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#=================================================================================#

# igraph mediation overview
 
# tutorial here: https://kateto.net/networks-r-igraph

# results from the main meta-analysis
nf <- newest_file(look_for = "STEP4_meta_resultsAnnotated$",subfolder = "res", print_full = TRUE)
res.init.meta <- fread(nf)

# joint results of mediation analysis
res.joint <- fread(newest_file("mediationStatisticsAnnotatedJoint", "res",print_full = TRUE))

# reformatted mediation analysis results
r2.mediation <- fread(newest_file("mediationStatisticsMatchingMarginalInfo", "res/",print_full = T))
relative.mediation <- fread(newest_file("mediatonGroupList","res/",print_full = T))

# for agn.all
load("/net/ifs2/san_projekte/projekte/genstat/02_projekte/1703_ge_metab_a1_b3_sorbs/1805_gxMetaboliteAssociation/obj/200306_metaboliteNames.RData")

#=================================================================================#
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#=================================================================================#

# Mediation Table ----

# supp table with mediation summary gene-centric
table2 <- r2.mediation[pheno=="log.bmi"]
table2[gene == "", gene := NA]

# join long metabolite names
table2$metab.abbr <- agn.all[match(table2$metabolite, rid), abbr]

# remove empty entries
table2 <- table2[!is.na(gene)]

# type column in line with the text!
table2[type=="pheno~gx(+metab)", type := "Gene Expression->Metabolite->Phenotype"]
table2[type=="pheno~metab(+gx)", type := "Metabolite->Gene Expression->Phenotype"]

# sum the results by gene
table2_genes <- table2[
  ,.(
     metabs_tested = uniqueN(metabolite),
     all_metabs_tested = paste(metab.abbr, collapse=", "),
     sum_sigs = sum(mediation.sig), 
     prop_sigs = sum(mediation.sig)/.N,
     metab_sigs = paste(metab.abbr[mediation.sig==TRUE], collapse=", "),
     strongest_med_metab = .SD[mediation.sig == T][abs(beta.mediation) == max(abs(beta.mediation), na.rm = T)][,metabolite],
     strongest_med_beta = .SD[mediation.sig == T][abs(beta.mediation) == max(abs(beta.mediation), na.rm = T)][,beta.mediation],
     strongest_med_ci_l = .SD[mediation.sig == T][abs(beta.mediation) == max(abs(beta.mediation), na.rm = T)][,ci.95.lower],
     strongest_med_ci_u = .SD[mediation.sig == T][abs(beta.mediation) == max(abs(beta.mediation), na.rm = T)][,ci.95.upper],
     strongest_med_p = .SD[mediation.sig == T][abs(beta.mediation) == max(abs(beta.mediation), na.rm = T)][,p.value.prod],
     strongest_med_fdr1 = .SD[mediation.sig == T][abs(beta.mediation) == max(abs(beta.mediation), na.rm = T)][,fdr.local],
     strongest_med_fdr2 = .SD[mediation.sig == T][abs(beta.mediation) == max(abs(beta.mediation), na.rm = T)][,fdr.global]
  ),
  .(type, gene)]

# order to highest sigs
table2_genes <- table2_genes[order(-sum_sigs)]

table2_genes <- table2_genes[,
  .(
    Mediation_Type = type,
    Mediation_Outcome = "log-BMI",
    Gene = gene,
    Number_Of_Tested_Metabolites = metabs_tested,
    Tested_Metabolites = all_metabs_tested,
    Significant_Mediations = sum_sigs,
    Proportion_Significant_To_Tested_Metabolites = round(prop_sigs, 2),
    Sigificant_Metabolites = metab_sigs,
    Strongest_Mediation_Metabolites = strongest_med_metab,
    Strongest_Mediation_Effect_Estimate = strongest_med_beta,
    `Strongest_Mediation_CI95%_Lower` = strongest_med_ci_l,
    `Strongest_Mediation_CI95%_Upper` = strongest_med_ci_u,
    Strongest_Mediation_P_Value = strongest_med_p,
    Strongest_Mediation_FDR_Level_1 = strongest_med_fdr1,
    Strongest_Mediation_FDR_Level_2 = strongest_med_fdr2
  )
]

# supp table with mediation summary metabolite-centric
table2_metab <- table2[
  ,.(
    genes_tested = uniqueN(gene),
    all_genes_tested = paste(gene, collapse=", "),
    sum_sigs = sum(mediation.sig), 
    prop_sigs = sum(mediation.sig)/.N,
    gene_sigs = paste(gene[mediation.sig==TRUE], collapse=", "),
    strongest_med_gene = .SD[mediation.sig == T][abs(beta.mediation) == max(abs(beta.mediation), na.rm = T)][,gene],
    strongest_med_beta = .SD[mediation.sig == T][abs(beta.mediation) == max(abs(beta.mediation), na.rm = T)][,beta.mediation],
    strongest_med_ci_l = .SD[mediation.sig == T][abs(beta.mediation) == max(abs(beta.mediation), na.rm = T)][,ci.95.lower],
    strongest_med_ci_u = .SD[mediation.sig == T][abs(beta.mediation) == max(abs(beta.mediation), na.rm = T)][,ci.95.upper],
    strongest_med_p = .SD[mediation.sig == T][abs(beta.mediation) == max(abs(beta.mediation), na.rm = T)][,p.value.prod],
    strongest_med_fdr1 = .SD[mediation.sig == T][abs(beta.mediation) == max(abs(beta.mediation), na.rm = T)][,fdr.local],
    strongest_med_fdr2 = .SD[mediation.sig == T][abs(beta.mediation) == max(abs(beta.mediation), na.rm = T)][,fdr.global]
  ),
  .(type, metab.abbr)]

# order
table2_metab <- table2_metab[order(-sum_sigs)]

# proper format
table2_metab <- table2_metab[
  ,
  .(
    Mediation_Type = type,
    Mediation_Outcome = "log-BMI",
    Metabolite = metab.abbr,
    Number_Of_Tested_Genes = genes_tested,
    Tested_Genes = all_genes_tested,
    Significant_Mediations = sum_sigs,
    Proportion_Significant_To_Tested_Genes = round(prop_sigs, 2),
    Sigificant_Genes = gene_sigs,
    Strongest_Mediation_Genes = strongest_med_gene,
    Strongest_Mediation_Effect_Estimate = strongest_med_beta,
    `Strongest_Mediation_CI95%_Lower` = strongest_med_ci_l,
    `Strongest_Mediation_CI95%_Upper` = strongest_med_ci_u,
    Strongest_Mediation_P_Value = strongest_med_p,
    Strongest_Mediation_FDR_Level_1 = strongest_med_fdr1,
    Strongest_Mediation_FDR_Level_2 = strongest_med_fdr2
  )
  ]

# supp table with all mediation results

# join long metabolite names
r2.mediation$metab.abbr <- agn.all[match(r2.mediation$metabolite, rid), abbr]
r2.mediation$metab.long <- agn.all[match(r2.mediation$metabolite, rid), full.name]

# format all mediation results table
supp_table_all <- r2.mediation[pheno=="log.bmi", .(
  metab = metabolite, # for matching, remove later
  Mediation_Type = type, 
  Mediation_Outcome = pheno, 
  Gene_Expression_Probe = gx.probe, 
  Gene = gene, 
  Metabolite_Abbreviation = metab.abbr, 
  Metabolite = metab.long, 
  Mediation_Effect_Estimate = beta.mediation, 
  `Mediation_Effect_CI95%_Lower` = ci.95.lower,
  `Mediation_Effect_CI95%_Lower` = ci.95.upper, 
  Proportion_Mediated = beta.mediation / (tauAdj + beta.mediation),
  Mediation_P_Value = p.value.prod,
  Mediation_FDR_Level_1 = fdr.local,
  Mediation_FDR_Level_2 = fdr.global, 
  Mediation_Significance = mediation.sig)]

# corrections as above
supp_table_all[Gene == "", Gene := NA]
supp_table_all[,Mediation_Outcome := "log-BMI"]
supp_table_all[Mediation_Type=="pheno~gx(+metab)", Mediation_Type := "Gene Expression->Metabolite->Phenotype"]
supp_table_all[Mediation_Type=="pheno~metab(+gx)", Mediation_Type := "Metabolite->Gene Expression->Phenotype"]

# remove mediation directions that were not tested!
supp_table_all <- supp_table_all[!is.na(Mediation_Effect_Estimate)]

# merge mediation group to the table

m1 <- match(
  supp_table_all[,paste0(Gene_Expression_Probe,metab)], 
  relative.mediation[pheno=="log.bmi",paste0(gx.probe,metabolite)]
)
supp_table_all$Mediation_Group <- relative.mediation[pheno=="log.bmi",][m1, mediation.group]
supp_table_all$metab <- NULL

# save everything

fwrite(
  table2_genes, 
  dpx("geneCentricMediationResults.tsv","ppr/supp_tables/"), 
  sep = "\t",
  quote = FALSE
)

fwrite(
  table2_metab,
  dpx("metaboliteCentricMediationResults.tsv","ppr/supp_tables/"), 
  sep = "\t",
  quote = FALSE
)

fwrite(
  supp_table_all,
  dpx("mediationAnalysisSumStats.tsv","ppr/supp_tables/"), 
  sep = "\t",
  quote = FALSE
)

#=================================================================================#
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#=================================================================================#

# BMI Network ----

#=================================================================================#
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#=================================================================================#

# Paper: Mediation of BMI genes? ----

#' get BMI genes

# PMID: 27605733
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

uniqueN(bmi.genes)

# create table with gene sources
supp_tab <- data.table(
  "Obesity-associated genes" = bmi.genes,
  "Source PMID" = "PLACEHOLDER",
  "Source DOI" = "PLACEHOLDER"
)

supp_tab[`Obesity-associated genes` == "PPM1K", 
         `:=`(
           `Source PMID` = "27898682",
           `Source DOI` = "https://doi.org/10.1371/journal.pmed.1002179"
         )
         ]
supp_tab[`Obesity-associated genes` != "PPM1K", 
         `:=`(
           `Source PMID` = "27605733",
           `Source DOI` = "https://doi.org/10.1007/s12291-015-0541-x"
         )]

WriteXLS_hk(x = supp_tab,
            ExcelFileName = "ppr/supp_tables/literatureObesityGenes.xls",
            SheetNames = "Obesity-associated genes")

# x in data
res.init.meta[symbol_INGENUITY %in% bmi.genes, uniqueN(symbol_INGENUITY)]
res.init.meta[symbol_INGENUITY %in% bmi.genes & significant == TRUE, .N]

# check mediations
r2.mediation[pheno=="log.bmi" & gene %in% bmi.genes, .N]
r2.mediation[pheno=="log.bmi" & gene %in% bmi.genes, uniqueN(gene)]
r2.mediation[pheno=="log.bmi" & gene %in% bmi.genes & mediation.sig == TRUE, ]
r2.mediation[pheno=="log.bmi" & gene %in% bmi.genes & mediation.sig == TRUE, uniqueN(gene)]
r2.mediation[pheno=="log.bmi" & gene %in% bmi.genes & mediation.sig == TRUE, uniqueN(metabolite)]
tmp1 <- r2.mediation[pheno=="log.bmi" & gene %in% bmi.genes & mediation.sig == TRUE, unique(metabolite)]

# mediated metabolite class
res.init.meta[metabolite %in% (tmp1), unique(metab.class),metabolite]$V1 %>% table

r2.mediation[pheno=="log.bmi" & gene %in% bmi.genes & mediation.sig == TRUE, .N]
r2.mediation[pheno=="log.bmi" & gene %in% bmi.genes & mediation.sig == TRUE, .N, by = .(type)]
r2.mediation[pheno=="log.bmi" & gene %in% bmi.genes & mediation.sig == TRUE, .N, by = .(type, metabolite)]
r2.mediation[pheno=="log.bmi" & gene %in% bmi.genes & mediation.sig == TRUE, .N, by = .(gene)]


# number of categorized mediations
relative.mediation[pheno=="log.bmi" & gene %in% bmi.genes & mediation.group %nin% c("Other", "None"), ]
metab <- relative.mediation[pheno=="log.bmi" & gene %in% bmi.genes & mediation.group %nin% c("Other", "None"), metabolite]
agn.all[rid %in% (metab),]

# number of mediations per direction
relative.mediation[pheno=="log.bmi" & gene %in% bmi.genes & mediation.group != c("None"), 
                   .(metab = sum(`mediation.sig_pheno~metab(+gx)`),
                     gx = sum(`mediation.sig_pheno~gx(+metab)`))]

# number of bi-directional mediations
relative.mediation[pheno=="log.bmi" & gene %in% bmi.genes & mediation.group != c("None") & 
                   (`mediation.sig_pheno~metab(+gx)` == T &
                     `mediation.sig_pheno~gx(+metab)` == T),]

# bi-directional strong mediations
relative.mediation[pheno=="log.bmi" & gene %in% bmi.genes & mediation.group %nin% c("None", "Other") & 
                   (`mediation.sig_pheno~metab(+gx)` == T &`mediation.sig_pheno~gx(+metab)` == T),]

# metabolites per genes
relative.mediation[pheno=="log.bmi" & gene %in% bmi.genes & mediation.group != c("None"), .N, gene]

# mehrfach mediierte metabolite
relative.mediation[pheno=="log.bmi" & gene %in% bmi.genes & mediation.group != c("None"), .N, metabolite][N==2,]
relative.mediation[pheno=="log.bmi" & gene %in% bmi.genes & mediation.group != c("None"), .(.N, paste(gene, collapse=", ")), metabolite][N==2,]
agn.all[rid %in% c("Q30","Q38","Q27"),]


# table of mediation results
tablex <- r2.mediation[pheno=="log.bmi" & gene %in% bmi.genes & mediation.sig == TRUE, ]

# formatting of table
tablex[, `:=`(
  Exposure = ifelse(type=="pheno~gx(+metab)", "Gene Expression", "Metabolite"),
  Mediator = ifelse(type=="pheno~gx(+metab)", "Metabolite", "Gene Expression")
)]

# match proper metabolite names
m1 <- match(tablex$metabolite, agn.all$rid)
tablex$`Metabolite Abbreviation` <- agn.all[m1, abbr]
tablex$Metabolite <- agn.all[m1, full.name]

tablex <- tablex[,.(
  Gene = gene,
  `Gene Expression Probe` = gx.probe,
  Metabolite,
  `Metabolite Abbreviation`,
  Exposure,
  Mediator,
  Outcome = "BMI",
  `Proportion Mediated` = round(alpha*beta/(alpha*beta+tauAdj),3),
  `Mediation effect size` = round(alpha*beta, 4),
  `Total effect` = round(alpha*beta+tauAdj, 3),
  `Explained varianc of Exposure` = ifelse(type=="pheno~gx(+metab)",round(r2.gx.raw,3),round(r2.metab.raw,3))
)]

WriteXLS_hk(tablex, ExcelFileName = dpx("bmiMediationTable.xls", "ppr/main_tables/"), SheetNames = "Table")

#' Display: raw assoc to phenotype, quotient, mediated r2, sign of assoc, significance t/F, direction of mediation

#=================================================================================#
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#=================================================================================#

#' # Colorscheme Metabolites/Genes
# Colorscheme Metabolites/Genes ----

pal.metabs <- sequential_hcl(palette = "OrYel",n=100)
pal.genes <- sequential_hcl(palette = "Teal",n=100)
pal.pairs <- sequential_hcl(palette = "Greens2",n=100)

#=================================================================================#
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#=================================================================================#

# BMI Network ----
#' ## BMI Network

# filter data including some renaming
d2 <- relative.mediation[
  gene != "" &
    pheno=="log.bmi" &
    mediation.group %nin% c("None"),
  .(pheno, 
    metabolite, 
    gene, 
    gx.probe, 
    mediation.group,
    sig.group.raw.gx = `raw.sig.group_pheno~gx(+metab)`,
    sig.group.raw.metab = `raw.sig.group_pheno~metab(+gx)`,
    sig.mediation.gx = `mediation.sig_pheno~gx(+metab)`,
    sig.mediation.metab = `mediation.sig_pheno~metab(+gx)`,
    beta.mediation.gx = `beta.mediation_pheno~gx(+metab)`,
    beta.mediation.metab = `beta.mediation_pheno~metab(+gx)`,
    quot.dir.vs.raw.gx = `prop.mediated_pheno~gx(+metab)`,
    quot.dir.vs.raw.metab = `prop.mediated_pheno~metab(+gx)`,
    sig.raw.gx = `sig.gx.raw_pheno~gx(+metab)`,
    sig.raw.metab = `sig.metab.raw_pheno~metab(+gx)`,
    r2.raw.gx = `r2.gx.raw_pheno~gx(+metab)`,
    r2.raw.metab = `r2.metab.raw_pheno~metab(+gx)`,
    r2.dir.gx = `r2.arnd_pheno~gx(+metab)`,
    r2.dir.metab = `r2.arnd_pheno~metab(+gx)`,
    beta.raw.gx = `beta.gx.raw_pheno~gx(+metab)`,
    beta.raw.metab = `beta.metab.raw_pheno~metab(+gx)`)]

#### FUNCTION ####

source("../functions/create_mediation_network_data.R")
source("../functions/mediation_network.R")

# gene.filter
my.filter <- bmi.genes

# create data
d2f <- create_mediation_network_data(
  d2 = d2,
  res.init.meta = res.init.meta,
  agn.all = agn.all, 
  gene.filter = my.filter, save.table = F
)

# fwrite(d2f, "obj/d2f.csv")

# plot network ----

tiff(dpx("r2BMIMediationNetwork.tiff","plt/"),res=300,width = 15,height = 15,unit="in")

# plot
mediation_network(d2f = d2f)

dev.off()

#' # TODO: Histogram of distribution of effect sizes in nodes/edges
# TODO: Histogram of distribution of effect sizes in nodes/edges ----

#=================================================================================#
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#=================================================================================#

# Full Network ----


# get all genes - filtering is done below based on max beta per gene
my.filter <- unique(d2$gene)

# additional optional filtering
d2[, max.mediation.beta := ifelse(abs(beta.mediation.gx) > abs(beta.mediation.metab), beta.mediation.gx, beta.mediation.metab)]

# create data
d2f2 <- create_mediation_network_data(
  d2 = d2,
  res.init.meta = res.init.meta,
  agn.all = agn.all, 
  gene.filter = my.filter, save.table = F
)

# filter to only display the 250 strongest each
d2f2 <- d2f2[(sig.mediation.gx == T | sig.mediation.metab == T) & mediation.group %nin% "Other", ][
  order(-abs(max.mediation.beta)), head(.SD,250),by=mediation.group]
d2f2[allDuplicatedEntries(gene)]

# save plot
tiff(dpx("r2BMIMediationObesityGenesIncludingOtherNetwork.tiff","plt/"),res=300,width = 20,height = 20,unit="in")

mediation_network(d2f = d2f2)

dev.off()

#=================================================================================#
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#=================================================================================#

# T2D Network ----

# General T2D mediation info ----

# mediations
r2.mediation[pheno=="diabetes.status.tri" & !is.na(mediation.sig) & mediation.sig==T & gene != "", .N]

# unique genes/metabolites
r2.mediation[pheno=="diabetes.status.tri" & !is.na(mediation.sig) & mediation.sig==T & gene != "", uniqueN(gene)]
r2.mediation[pheno=="diabetes.status.tri" & !is.na(mediation.sig) & mediation.sig==T & gene != "", uniqueN(metabolite)]

# type
r2.mediation[pheno=="diabetes.status.tri" & !is.na(mediation.sig) & mediation.sig==T & gene != "", table(type)]

# master-mediater metabolites

# as exposure
r2.mediation[pheno=="diabetes.status.tri" & !is.na(mediation.sig) & mediation.sig==T & gene != "", 
             .(uniqueN(gene)),by=.(metabolite,type)][order(-V1),][
               type=="pheno~metab(+gx)"] %>% head(5)

# as mediator
r2.mediation[pheno=="diabetes.status.tri" & !is.na(mediation.sig) & mediation.sig==T & gene != "", 
             .(uniqueN(gene)),by=.(metabolite,type)][order(-V1),][
               type=="pheno~gx(+metab)"] %>% head(5)

# master mediator genes
# as exposure
r2.mediation[pheno=="diabetes.status.tri" & !is.na(mediation.sig) & mediation.sig==T & gene != "", 
             .(uniqueN(metabolite)),by=.(gene,type)][order(-V1),][
               type=="pheno~gx(+metab)"] %>% head(5)

# as mediator
r2.mediation[pheno=="diabetes.status.tri" & !is.na(mediation.sig) & mediation.sig==T & gene != "", 
             .(uniqueN(metabolite)),by=.(gene,type)][order(-V1),][
               type=="pheno~metab(+gx)"] %>% head(7)

# get proportion  mediated
r2.mediation[!is.na(mediation.sig), pm := alpha*beta / (alpha*beta + tauAdj)]

# strongest mediations by beta
r2.mediation[pheno=="diabetes.status.tri" & !is.na(mediation.sig) & mediation.sig==T & gene != "", ][
  order(-abs(beta.mediation)), ][
    , .(
      type,
      gene,
      gx.probe,
      metabolite,
      beta=round(beta.mediation,3),
      cil = round(ci.95.lower,3),
      cil = round(ci.95.upper,3),
      pm = round(pm,2),
      pval = p.value.prod
    )] %>% as.data.frame() %>% head(15)

# strongest mediations by PM
r2.mediation[pheno=="diabetes.status.tri" & !is.na(mediation.sig) & mediation.sig==T & gene != "", ][
  order(-abs(pm)), ][
    , .(
      type,
      gene,
      gx.probe,
      metabolite,
      beta=round(beta.mediation,3),
      cil = round(ci.95.lower,3),
      cil = round(ci.95.upper,3),
      pm = round(pm,2),
      pval = p.value.prod
    )] %>% as.data.frame() %>% head(15)


ggplot(r2.mediation[!is.na(mediation.sig) & mediation.sig==T & gene != "",],
       aes(pm)
       ) +
  geom_histogram(col="white") +
  facet_grid(pheno~type) +
  theme_gxMetab()

#=================================================================================#
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#=================================================================================#

# Mediation of T2D genes ----

# 2017: PMID: 28741265  DOI: 10.1007/s11892-017-0908-x
t2d.genes1 <- c("TCF7L2", "MTNR1B", "PPARG", "FTO", "ANK1", "CENTD2", "CAMK1D", "ZFAND3", "JAZF1", "KCNQ1", "KCNK17", "WFS1")

# 2015: PMID: 26119161  DOI: 10.1111/1753-0407.12323
t2d.genes2 <- c("NOTCH2", "PROX1", "THADA", "IRS1", "GCKR", "RBMS1", "BCL11A", "IGF2BP2", "PPARG", "ADAMTS9", "ADCY5", "WFS1", "ZBED3", "ANKRD55","CDKAL1", "JAZF1", "KLF14", "GCK", "DGKB/TMEM195", "SLC30A8", "TP53INP1", "CDKN2A/B", "TLE4", "TLE1", "TCF7L2", "HHEX", "IDE","CDC123", "CAMK1D", "ZMIZ1", "KCNJ11", "ARAP1", "MTNR1B", "DUSP8", "LGR5", "HMGA2", "HNF1A", "KLHDC5", "CCND2", "ZFAND6", "PRC1", "FTO", "BCAR1", "HNF1B", "MC4R", "GATAD2A", "PBX4", "DUSP9")

# get all genes
t2d.genes <- unique(c(t2d.genes1,t2d.genes2))
t2d.genes %>% length

# genes in meta-analysis
res.init.meta[symbol_INGENUITY %in% t2d.genes, uniqueN(symbol_INGENUITY)]

# sig genes in meta-analysis
res.init.meta[symbol_INGENUITY %in% t2d.genes & significant == T, uniqueN(symbol_INGENUITY)]

# genes in mediation
r2.mediation[gene %in% (t2d.genes) & pheno == "diabetes.status.tri", uniqueN(gene)]

# sig mediating/-ed genes
r2.mediation[gene %in% (t2d.genes) & pheno == "diabetes.status.tri" & mediation.sig == TRUE,]
r2.mediation[gene %in% (t2d.genes) & pheno == "diabetes.status.tri" & mediation.sig == TRUE, uniqueN(gene)]

# number of mediations tested
r2.mediation[pheno=="diabetes.status.tri" & gene %in% t2d.genes, .N]

# number of sig mediations
r2.mediation[pheno=="diabetes.status.tri" & gene %in% t2d.genes & mediation.sig == TRUE, ]

# number of sig mediation genes
r2.mediation[pheno=="diabetes.status.tri" & gene %in% t2d.genes & mediation.sig == TRUE, unique(gene)]
r2.mediation[pheno=="diabetes.status.tri" & gene %in% t2d.genes & mediation.sig == TRUE, uniqueN(gene)]

# number of sig mediation metabolites
r2.mediation[pheno=="diabetes.status.tri" & gene %in% t2d.genes & mediation.sig == TRUE, uniqueN(metabolite)]

# mediated metabolite class
tmp1 <- r2.mediation[pheno=="diabetes.status.tri" & gene %in% t2d.genes & mediation.sig == TRUE, unique(metabolite)]
res.init.meta[metabolite %in% (tmp1), unique(metab.class),metabolite]$V1 %>% table

r2.mediation[pheno=="diabetes.status.tri" & gene %in% t2d.genes & mediation.sig == TRUE, .N]
r2.mediation[pheno=="diabetes.status.tri" & gene %in% t2d.genes & mediation.sig == TRUE, .N, by = .(type)]

# mediations & metabolites per gene
r2.mediation[pheno=="diabetes.status.tri" & gene %in% t2d.genes & mediation.sig == TRUE, .(.N,uniqueN(metabolite)), by = .(gene)][
  order(-N)
]

# number of categorized mediations
relative.mediation[pheno=="diabetes.status.tri" & gene %in% t2d.genes & mediation.group %nin% c("Other", "None"), ]
metab <- relative.mediation[pheno=="diabetes.status.tri" & gene %in% t2d.genes & mediation.group %nin% c("Other", "None"), metabolite]
agn.all[rid %in% (metab),]

# number of bi-directional mediations
relative.mediation[pheno=="diabetes.status.tri" & gene %in% t2d.genes & mediation.group != c("None") & 
                     (`mediation.sig_pheno~metab(+gx)` == T &
                        `mediation.sig_pheno~gx(+metab)` == T),]

# who is bi-directional
relative.mediation[pheno=="diabetes.status.tri" & gene %in% t2d.genes & mediation.group != c("None") & 
                     (`mediation.sig_pheno~metab(+gx)` == T &
                        `mediation.sig_pheno~gx(+metab)` == T), table(gene)]

# bi-directional strong mediations
relative.mediation[pheno=="diabetes.status.tri" & gene %in% t2d.genes & mediation.group %nin% c("None", "Other") & 
                     (`mediation.sig_pheno~metab(+gx)` == T &`mediation.sig_pheno~gx(+metab)` == T),]

# mehrfach mediierte metabolite
relative.mediation[pheno=="diabetes.status.tri" & gene %in% t2d.genes & mediation.group != c("None") & metabolite =="Tyr", ]
relative.mediation[pheno=="diabetes.status.tri" & gene %in% t2d.genes & mediation.group != c("None"), .N, metabolite][N==2,]
relative.mediation[pheno=="diabetes.status.tri" & gene %in% t2d.genes & mediation.group != c("None"), .(.N, paste(gene, collapse=", ")), metabolite][N==2,][order(-N)]
agn.all[rid %in% c("Q27"),]

# table of mediation results
tablex <- r2.mediation[pheno=="diabetes.status.tri" & gene %in% t2d.genes & mediation.sig == TRUE, ]

# formatting of table
tablex[, `:=`(
  Exposure = ifelse(type=="pheno~gx(+metab)", "Gene Expression", "Metabolite"),
  Mediator = ifelse(type=="pheno~gx(+metab)", "Metabolite", "Gene Expression")
)]

# match proper names
m1 <- match(tablex$metabolite, agn.all$rid)
tablex$`Metabolite Abbreviation` <- agn.all[m1, abbr]
tablex$Metabolite <- agn.all[m1, full.name]

# create table
tablex <- tablex[,.(
  Gene = gene,
  `Gene Expression Probe` = gx.probe,
  Metabolite,
  `Metabolite Abbreviation`,
  Exposure,
  Mediator,
  Outcome = "Type 2 Diabetes",
  `Proportion Mediated` = round(alpha*beta/(alpha*beta+tauAdj),3),
  `Mediation effect size` = round(alpha*beta, 4),
  `Total effect` = round(alpha*beta+tauAdj, 3),
  `Explained varianc of Exposure` = ifelse(type=="pheno~gx(+metab)",round(r2.gx.raw,3),round(r2.metab.raw,3))
)]
tablex[order(-`Proportion Mediated`)]

# save table
WriteXLS_hk(tablex, ExcelFileName = dpx("t2dMediationTable.xls", "ppr/main_tables/"), SheetNames = "Table")

#=================================================================================#
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#=================================================================================#


# filter data including some renaming
d2 <- relative.mediation[
  gene != "" &
    pheno=="diabetes.status.tri" &
    mediation.group %nin% c("None"),
  .(pheno, 
    metabolite, 
    gene, 
    gx.probe, 
    mediation.group,
    sig.group.raw.gx = `raw.sig.group_pheno~gx(+metab)`,
    sig.group.raw.metab = `raw.sig.group_pheno~metab(+gx)`,
    sig.mediation.gx = `mediation.sig_pheno~gx(+metab)`,
    sig.mediation.metab = `mediation.sig_pheno~metab(+gx)`,
    beta.mediation.gx = `beta.mediation_pheno~gx(+metab)`,
    beta.mediation.metab = `beta.mediation_pheno~metab(+gx)`,
    quot.dir.vs.raw.gx = `prop.mediated_pheno~gx(+metab)`,
    quot.dir.vs.raw.metab = `prop.mediated_pheno~metab(+gx)`,
    sig.raw.gx = `sig.gx.raw_pheno~gx(+metab)`,
    sig.raw.metab = `sig.metab.raw_pheno~metab(+gx)`,
    r2.raw.gx = `r2.gx.raw_pheno~gx(+metab)`,
    r2.raw.metab = `r2.metab.raw_pheno~metab(+gx)`,
    r2.dir.gx = `r2.arnd_pheno~gx(+metab)`,
    r2.dir.metab = `r2.arnd_pheno~metab(+gx)`,
    beta.raw.gx = `beta.gx.raw_pheno~gx(+metab)`,
    beta.raw.metab = `beta.metab.raw_pheno~metab(+gx)`)]

# gene.filter
my.filter <- t2d.genes

# create data
d2f <- create_mediation_network_data(
  d2 = d2,
  res.init.meta = res.init.meta,
  agn.all = agn.all, 
  gene.filter = my.filter, save.table = T
)

# plot network ----

tiff(dpx("r2T2DMediationNetwork.tiff","plt/"),res=300,width = 20,height = 20,unit="in")

# plot
mediation_network(d2f = d2f)

dev.off()

