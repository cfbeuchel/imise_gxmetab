#' ---
#' title: "Mendelian Randomisation"
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
  "magrittr"
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
setDTthreads(8)
option_setup()

#+ include=F
#==========#
#   DONE   #
#==========#

# Mediation Results ----
#' # Mediation Results 

# res_bmi <- readxl::read_xls("ppr/main_tables/201001_bmiMediationTable.xls")
# setDT(res_bmi)
res_bmi <- fread(newest_file("mediationBMIGenesTable","res",print_full = T))
# res_t2d <- fread(newest_file("mediationBMIGenesTable","res",print_full = T))

# Data for Bi-Directional MR ----
#' # Data for Bi-Directional MR

#' This means that I'll first establish the causal link gx<->metab

# load metabolite IDs
load("/net/ifs2/san_projekte/projekte/genstat/02_projekte/1703_ge_metab_a1_b3_sorbs/1805_gxMetaboliteAssociation/obj/200306_metaboliteNames.RData")

# mQTL ----
#' ## mQTL

# mapping: grch37
top.mqtl <- fread("/net/ifs2/san_projekte/projekte/genstat/02_projekte/1512_meta_meta/03_Meta/03_metaAnnotation/gwasresults1.7_alle_I2_90/synopsis/topliste_tabdelim/topliste_2019-09-19_alle_I2_90.txt")

# GTEx eQTL ----
#' ### GTEx eQTL

# mapping: grch38
gtex.all.eqtl <- fread("/net/ifs2/san_projekte/projekte/genstat/01_daten/2007_GTEx_v8/GTEx_Analysis_v8_QTLs_GTEx_Analysis_v8_eQTL_all_associations_Whole_Blood.allpairs.txt.gz")


gtex.all.eqtl %>% head
setkey(gtex.all.eqtl, variant_id)

# get position
variant.split <- strsplit(x = gtex.all.eqtl$variant_id, split = "_")
gtex.all.eqtl$chr <- unlist(lapply(variant.split, `[[`, 1))
gtex.all.eqtl[, chr := gsub(pattern = "chr", replacement = "", chr)]
gtex.all.eqtl$pos <- as.integer(unlist(lapply(variant.split, `[[`, 2)))
gtex.all.eqtl$effect.allele <- unlist(lapply(variant.split, `[[`, 3))
gtex.all.eqtl$other.allele <- unlist(lapply(variant.split, `[[`, 4))
rm(variant.split)
gtex.all.eqtl %>% head %>% str

# Match ENSG ID ----

# object genes
load("/net/ifs2/san_projekte/projekte/genstat/07_programme/rtools/140813_neighbouring_genes/data/ensembl_genes_GRCh37.rel91_1712hk.RData")
setDT(genes)
genes38 <- fread("/net/ifs2/san_projekte/projekte/genstat/02_projekte/1703_ge_metab_a1_b3_sorbs/1805_gxMetaboliteAssociation/dat/201022_GRCh38.p13_ensembleGenes.gz")

# annotate genes

# removethe version indicator
gtex.all.eqtl[, gene_id_no_version := gsub(x=gene_id, "\\.[0-9]+$", "")]

# check gene overlap
tmp1 <- venn2(gtex.all.eqtl$gene_id_no_version,
             genes38$`Gene stable ID`, 
             plotte = T)
tmp2 <- venn2(gtex.all.eqtl$gene_id,
             genes38$`Gene stable ID version`)

m1 <- match(gtex.all.eqtl$gene_id_no_version, genes$ensembl_gene_id)
matched.ids <- genes[m1,associated_gene_name]
gtex.all.eqtl[, gene.id:=(matched.ids)]

# large objects - remove
rm(m1, matched.ids)

# how to match IDs

# Step1: top eQTL (grch38) suchen, in mQTL (grch37) nachschlagen
# -> da mQTL viel mehr ist, convert eqtl von 38 zu 37

# reduce in size for lighter lift over
allsnp <- unique(gtex.all.eqtl[, .(variant_id, chr, pos)])
allsnp[, snps := variant_id]
allsnp$variant_id <- NULL
allsnp[chr == "X", chr := "23"]
allsnp[, chr := as.numeric(chr)]
allsnp$chr %>% unique %>% as.numeric %>% sort

# all ensbl Ids of GTEX (in grch38)
allgenes <- unique(gtex.all.eqtl[, .(gene_id)])

# LIFTOVER GTEx SNP to mQTL mapping ----

# use file provided by gtex
gtex8annot <- fread("../../../01_daten/2007_GTEx_v8/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.lookup_table.txt.gz")

# match b37-id to gtex
setkey(gtex8annot, variant_id)
setkey(allsnp, snps)
allsnp$variant_id_b37 <- gtex8annot[allsnp, variant_id_b37]

# check
allsnp[, variant_id_b38 := snps]
allsnp[, snps := NULL]

# match ENSG gene id to HGNC Identifier

# remove id version
allgenes[, gene := gsub("\\.\\d+" , "",gene_id)]

# highest overlap in gtex and mapped ENSG IDs?
v1 <- venn2(allgenes$gene, genes38$`Gene stable ID`)
v2 <- venn2(allgenes$gene_id, genes38$`Gene stable ID version`)

# match HGNC symbol to allgenes
m1 <- match(allgenes$gene, genes38$`Gene stable ID`)
allgenes[, `:=`(
  hgnc = genes38[m1, `Gene name`],
  hgnc_description = genes38[m1, `Gene description`],
  gene_start = genes38[m1, `Gene start (bp)`],
  gene_end = genes38[m1, `Gene end (bp)`]
)]


# Data for MR ----
#' # Data for 2-Step MR

# BMI Hits ----
#' ## BMI Hits

# res.bmi <- fread("/net/ifs2/san_projekte/projekte/genstat/02_projekte/1703_ge_metab_a1_b3_sorbs/1805_gxMetaboliteAssociation/dat/bmi_t2d_sum_stats/200508_smr_bmi.txt")

#=#################################################################=#
#===================================================================#
#=#################################################################=#

# Bi-Directional MR ----

#' # Select Instruments from top mediations
# Select Instruments ----

# unique is the same!
lookup <- unique(res_bmi[pheno == "log.bmi" & mediation.group %nin% c("None"), .(gx.probe,gene,metabolite)])
lookup[, metab.search := metabolite]
lookup[grepl("\\:", metabolite), metab.search := gsub(pattern = "\\:.*", "", metabolite)]
lookup[,index:=1:.N]
# lookup <- unique(res.med[pheno == "log.bmi" & mediation.group %nin% c("None","Other"), .(gx.probe,gene,metabolite)])

#=#################################################################=#
#===================================================================#
#=#################################################################=#

# INSTRUMENTS G->Gx ----
#' ## Instruments G->Gx

#' Find instruments for these probes:

sig.probes <- lookup[, unique(gx.probe)]
sig.genes <- lookup[, unique(gene)]
length(sig.probes)
length(sig.genes)

# G->Gx SNP in GTEx ----
#' ### G->Gx SNP in GTEx

# TODO: check r2 between SNP via plink ----
# TODO: Sind das nur tag-SNP? ----

# overlap gtex genes and my sig probes
v6 <- venn2(unique(gtex.all.eqtl$gene.id),sig.genes)

# get genome wide sig eqtl
fdr.gtex <- addHierarchFDR(pvalues = gtex.all.eqtl$pval_nominal,
                           categs = gtex.all.eqtl$gene_id
                           )
gtex.all.eqtl[, `:=`(
  fdr.local = fdr.gtex$fdr_level1,
  fdr.global = fdr.gtex$fdr_level2,
  fdr.hierarch.5perc = fdr.gtex$hierarch_fdr5proz
)]

# how many sig genes
gtex_sig <- gtex.all.eqtl[,sum(fdr.hierarch.5perc==T),gene_id][order(-V1),][V1 != 0, ]
gtex.all.eqtl[,sum(fdr.hierarch.5perc==T),gene_id][, table(V1>0)]
# gtex.all.eqtl[,sum(fdr.hierarch.5perc==T),gene_id][, boxplot(V1)]

gtex.all.eqtl[ , sig.genomewide := ifelse(pval_nominal <= 5e-8, T, F)]

# how many bonferroni/fdr sigs?
# this shows that all bonferroni sigs are also 
gtex.all.eqtl[,.N]
gtex_sig_fdr_bf <- gtex.all.eqtl[
  , .(fdr = sum(fdr.hierarch.5perc), 
  bonferroni = sum(sig.genomewide))
  ]
gtex_sig_fdr_bf
gtex.all.eqtl[, xtabs_hk(~fdr.hierarch.5perc+sig.genomewide)]

# are all genes represented? 
v7 <- venn2(
  gtex.all.eqtl[, gene.id],
  sig.genes
  )

# are sig snps for all genes found?
v7 <- venn3(
  gtex.all.eqtl[fdr.hierarch.5perc == T, gene.id],
  gtex.all.eqtl[sig.genomewide == T, gene.id], 
  sig.genes)

# which is the non sig gene
v7 %>% str

# no sig hit in gtex!
gtex.all.eqtl[gene.id == "FANCL", min(pval_nominal)]


# instruments!
# sig.eqtl <- gtex.eqtl[gene_name %in% sig.genes, ] # empty!
sig.eqtl <- gtex.all.eqtl[sig.genomewide == TRUE & gene.id %in% sig.genes, ]

# do I have sig snp for all my genes?
sig.eqtl$gene.id %>% unique %>% length == length(sig.genes) # do I have all my genes

# overview
sig.eqtl[,.N,gene.id]

# match the liftover b37 variant IDs 
m1 <- match(sig.eqtl$variant_id, allsnp$variant_id_b38)
sig.eqtl$variant_id_b37 <- allsnp[m1, variant_id_b37]


b37id <- strsplit(sig.eqtl$variant_id_b37, split = "_",fixed = T)
sig.eqtl$chr_b37 <- unlist(lapply(b37id, `[[`, 1))
sig.eqtl$pos_b37 <- unlist(lapply(b37id, `[[`, 2))
sig.eqtl$effect_allele_b37 <- unlist(lapply(b37id, `[[`, 3))
sig.eqtl$other_allele_b37 <- unlist(lapply(b37id, `[[`, 4))

#=#################################################################=#
#===================================================================#
#=#################################################################=#

# Corresponding G->metab SNP in MetaMeta ----
#' ### Corresponding G->metab SNP in MetaMeta

# lookup only for eqtl I have isntruments for
# lookup eqtl IN mqtl data
lookup.eqtl <- lookup[gene %in% unique(sig.eqtl$gene.id), ]

# file and location of step 10 objects
step.10.loc <- "/net/ifs2/san_projekte/projekte/genstat/02_projekte/1512_meta_meta/03_Meta/03_metaAnnotation/gwasresults1.7_alle_I2_90/"
step.10.name <- "step10_2_primaryANDsecondary_"


# metab.lookup <- lookup.eqtl[, .(metab=unique(metabolite))]
# 
# # clean the quotient names in the file
# metab.lookup[, metab.search := metab]
# metab.lookup[grepl("\\:", metab), metab.search := gsub(pattern = "\\:.*", "", metab)]

# replace name in search
lookup.eqtl[metab.search == "AC-total",metab.search := "acges"]

lookup.eqtl[, file := grep(
  pattern = paste0(step.10.name, metab.search, ".RData"), 
  x = list.files(step.10.loc),
  value = T), by=metab.search]

# check for missing files
stopifnot(lookup.eqtl[is.na(file), .N] == 0)
lookup.eqtl[is.na(file)]

# indexing
lookup.eqtl[,index:=1:.N]

# loop through each metameta object and find snp for which i have an eqtl
all.metab.eqtl.metab.snp <- lapply(lookup.eqtl$index, function(i){
  
  my.metab <- lookup.eqtl[index==(i),metabolite]
  my.file <- lookup.eqtl[index==(i),file]
  
  res.snp <- tryCatch({
    # load file
    load(paste0(step.10.loc,my.file))
    # stopifnot(any(ls()=="erg1"))
    # erg1 %>% head
    
    # which genes to look for snp for
    my.genes <- lookup.eqtl[metabolite==(my.metab), unique(gene)]
    my.snps <- sig.eqtl[gene.id %in% my.genes, 
                        .(gene.id,
                          # rs.id = rs_id_dbSNP147_GRCh37p13,
                          id = variant_id_b37,
                          
                          # nach id_meta (ref/alt wie gtex)
                          snp = paste(paste0("chr",chr_b37),effect_allele_b37,other_allele_b37,sep=":"),
                          snp.reverse = paste(paste0("chr",chr_b37),pos,other_allele_b37,effect_allele_b37,sep=":"),
                          
                          # nach markername (alt/ref - verkehrt zu gtex)
                          # snp = paste(rs_id_dbSNP147_GRCh37p13,pos,ref,alt,sep=":"),
                          # snp.reverse = paste(rs_id_dbSNP147_GRCh37p13,pos,alt,ref,sep=":"),
                          effect_allele_b37,
                          other_allele_b37,
                          index = 1:.N)]
    
    # look through ech gene
    res.snp <- lapply(my.snps$index, function(ii){
      
      ii.gene <- my.snps[index==(ii),gene.id]
      ii.rs <- my.snps[index==(ii),id]
      ii.snp <- my.snps[index==(ii),snp]
      ii.snp.rev <- my.snps[index==(ii),snp.reverse]
      
      # search for snp
      res.snp <- erg1[id_meta %in% c(ii.snp, ii.snp.rev), ]
      res.snp[,ref.allele.like.gtex := ifelse(id_meta == ii.snp, TRUE, FALSE)]
      
      # annotate gene
      res.snp[, instrument.gene := (ii.gene)]
      
    }) %>% rbindlist
    
    # coment
    res.snp[,comment:=NA]
    
    # remove huge object
    rm(erg1)
    
    return(res.snp)
    
  },error = function(cond){
    res.snp <- data.table(
      pheno=my.metab,
      comment=as.character(cond))
    return(res.snp)
  })
  
  return(res.snp)
})

# save!
all.metab.eqtl.metab.snp <- all.metab.eqtl.metab.snp %>% rbindlist(fill = T,use.names = T)
fwrite(all.metab.eqtl.metab.snp, dpx(suffix = "allSNPeQTLInstrumentsForMetab.csv", "res/"),sep = ",")
# all.metab.eqtl.metab.snp <- newest_file("llSNPeQTLInstrumentsForMetab","/res/",print_full = T) %>% fread()

# rename
eqtl.gx <- sig.eqtl
eqtl.metab <- all.metab.eqtl.metab.snp

# overlap
eqtl.metab$ref.allele.like.gtex %>% table # all other direction in mqtl data
venn3(
  eqtl.gx[, unique(paste(paste0("chr",chr_b37),pos_b37,effect_allele_b37,other_allele_b37,sep=":"))],
  eqtl.gx[, unique(paste(paste0("chr",chr_b37),pos_b37,other_allele_b37,effect_allele_b37,sep=":"))],
  eqtl.metab[, unique(id_meta)],
  mylabels = c("gtex", "gtex reverse", "mqtl")
)

# TODO: Match results (unified ID & beta) ----

#=#################################################################=#
#===================================================================#
#=#################################################################=#

# INSTRUMENTS G->Metabolite ----
#' ## Instruments G->Metabolite

# G->Metabolite SNP in MetaMeta ----
#' ### G->Metabolite SNP in MetaMeta

#' Find instruments for these metabolites:

sig.metabs <- res_bmi[mediation.group %nin% c("None"), unique(metabolite)]
length(sig.metabs)
sig.metabs <- agn.all[abbr %in% sig.metabs, rid]


# get colnames from mqtl data
v1 <- venn2(
  unlist(lapply(strsplit(names(top.mqtl)[grep("beta_score_", names(top.mqtl))],split = "_",fixed = T),`[[`,3)),
  sig.metabs,
  mylabels = c("mqtl","mediation")
)

# filter unnecessary metabolites
top.mqtlf <-  top.mqtl[,.SD,.SDcols=c(
  "markername", "duplicated_position", "invalidAssoc", "chr", 
  "cyto", "chunk", "pos", "tagger", "r2_tagger", "tagsnp", "id_meta", 
  "numberStudies", "cochQ", "cochQpval", "I2", "minMAF", "minInfoScore", 
  "Direction", "nearestgenes", "gene_biotype", "nearestgene", "Eigen", 
  "EigenPC", "CADD_scaled", "DANN", "GWAVA_Region", "GWAVA_TSS", 
  "GWAVA_Unmatched", "regulome_score", "regulome_score_numeric", 
  "regulome_details", "corinfo_gwas2", "cisgene", "transgene", 
  "coremine_genes", "hgnc4pathway", "entrez4pathway", "KEGG", "reactome", 
  "DOSE", "GO", "tissues", "markertoppheno", "toppheno", "toplogp_score", 
  "toplogp_expect", "r2_snp_toppheno", "strand", "effect_allele", 
  "other_allele", "topmaf", "topeaf", "topn", "topinfo", "topimputation_type",
  grep(
    pattern = paste(paste0(sig.metabs,"$"),collapse = "|"),
    x = names(top.mqtl),
    value = T)
)]

# match all the measure variables to 
name.list <- list(
  unlist(lapply(strsplit(grep(x = names(top.mqtlf),pattern = "^beta_score_",value = T),"_"),`[[`,3)),
  unlist(lapply(strsplit(grep(x = names(top.mqtlf),pattern = "^se_score_",value = T),"_"),`[[`,3)),
  unlist(lapply(strsplit(grep(x = names(top.mqtlf),pattern = "^I2_",value = T),"_"),`[[`,2)),
  unlist(lapply(strsplit(grep(x = names(top.mqtlf),pattern = "^logp_score_",value = T),"_"),`[[`,3))
)

# check wheter the names are equal
for(i in name.list){
  x <- all(
    sapply(
      name.list,identical,i
    )
  )
  stopifnot(x)
}

# melt

# CAVEAT: Ist das richtig? ----

top.mqtlfm <- melt(top.mqtlf,  
                   id.vars = c("markername", "duplicated_position", "invalidAssoc", "chr", 
                               "cyto", "chunk", "pos", "tagger", "r2_tagger", "tagsnp", "id_meta", 
                               "numberStudies", "cochQ", "cochQpval", "I2", "minMAF", "minInfoScore", 
                               "Direction", "nearestgenes", "gene_biotype", "nearestgene", "Eigen", 
                               "EigenPC", "CADD_scaled", "DANN", "GWAVA_Region", "GWAVA_TSS", 
                               "GWAVA_Unmatched", "regulome_score", "regulome_score_numeric", 
                               "regulome_details", "corinfo_gwas2", "cisgene", "transgene", 
                               "coremine_genes", "hgnc4pathway", "entrez4pathway", "KEGG", "reactome", 
                               "DOSE", "GO", "tissues", "markertoppheno", "toppheno", "toplogp_score", 
                               "toplogp_expect", "r2_snp_toppheno", "strand", "effect_allele", 
                               "other_allele", "topmaf", "topeaf", "topn", "topinfo", "topimputation_type"),
                   measure.vars = list(
                     grep(x = names(top.mqtlf),pattern = "^beta_score_",value = T),
                     grep(x = names(top.mqtlf),pattern = "^se_score_",value = T),
                     grep(x = names(top.mqtlf),pattern = "^I2_",value = T),
                     grep(x = names(top.mqtlf),pattern = "^logp_score_",value = T)
                   ),
                   value.name = c("beta","se","i2","logp"),
                   variable.factor = F,
                   value.factor = F
)

# match metabolite names to results
matchit <- data.table(
  metab=name.list[[1]],
  num=1:length(name.list[[1]])
)
m1 <- match(top.mqtlfm$variable,matchit$num)
top.mqtlfm$metabolite <- matchit[m1,metab]

# testing = looks good!
top.mqtlfm[markername=="rs10684465:131886278:C:CTT" & metabolite=="Pro",.(markername,chunk,pos,metabolite,beta,se,logp)]
top.mqtl[markername=="rs10684465:131886278:C:CTT",.(markername,chunk,pos,beta_score_Pro,se_score_Pro,logp_score_Pro)]

# check number of studies on tagsnp
top.mqtlfm[tagsnp==TRUE, table(numberStudies)]
top.mqtlfm[tagsnp==TRUE, table(toppheno)]

# get the top snp per phenotype
# pos chro
tag.mqtl <- top.mqtlfm[metabolite %in% sig.metabs & logp > -log10(5e-8),
                       .(tagsnp, markername, metabolite, toppheno,
                         markertoppheno, beta,se,i2,logp, 
                         effect_allele, other_allele,
                         topmaf, topinfo, pos, chr
                         )
                       ]
tag.mqtl[,hist(logp)]
tag.mqtl[,range(logp)]
tag.mqtl[,table(toppheno)]
tag.mqtl[,table(markername)][tag.mqtl[,table(markername)]>1]
tag.mqtl[,.N,by=toppheno]

# merge rid to res_bmi
m1 <- match(res_bmi$metabolite, agn.all$abbr)
res_bmi$rid <- agn.all[m1, rid]

# Overlap of tagsnp 
v2 <- venn2(
  tag.mqtl[,unique(metabolite)],
  res_bmi[mediation.group %nin% c("None"),unique(rid)]
)

# check for strong snp on remaining phenos -> DISREGARD
top.mqtl.nontag <- top.mqtlfm[metabolite %in% v2$q3 & 
                                logp>-log10(5e-8), 
                              .(markername,r2_tagger,chr,cyto,chunk,pos,toppheno,metabolite,beta,se,i2,logp,tagger,cisgene)][
                                , .SD[logp == max(logp)],by=.(metabolite,chr,cyto)]

# get list of SNP 
venn2(top.mqtl.nontag$markername,tag.mqtl$markername)
top.mqtl.nontag$markername %>% length
top.mqtl.nontag$markername %>% uniqueN
tag.mqtl$markername %>% length
tag.mqtl$markername %>% uniqueN

# get all mqtl
length(tag.mqtl$markername)
uniqueN(tag.mqtl$markername)

# check epistatic snp
table(tag.mqtl$markername) %>% sort(decreasing = T)
tag.mqtl[markername=="rs12210538:110760008:A:G", metabolite]

# extract rs ID
tag.mqtl$rs.id <- unlist(lapply(strsplit(tag.mqtl$markername,split = ":",fixed = T),`[[`,1))

mqtl <- unique(tag.mqtl[,.(chr,markername, pos)])
length(mqtl)
mqtl$rs <- unlist(lapply(strsplit(mqtl$markername,split = ":",fixed = T),`[[`,1))

#=#################################################################=#
#===================================================================#
#=#################################################################=#

# Corresponding G->Gx SNP in GTEx ----
#' ### Corresponding G->Gx SNP in GTEx

# via lookup
lookup.mqtl <- lookup[metabolite %in% unique(tag.mqtl$metabolite), ]

v1<-venn3(
  tag.mqtl[,paste0("chr", chr, "_",pos,"_",effect_allele,"_",other_allele, "_b38")],
  tag.mqtl[,paste0("chr", chr, "_",pos,"_",other_allele,"_",effect_allele, "_b38")],
  gtex.all.eqtl$variant_id,
  mylabels = c("mQTL", "mQTL rev", "Gtex v8")
  )

# check gene overlap with lookup genes
lookup$gene %>% uniqueN
lookup.mqtl$gene %>% uniqueN
gtex.all.eqtl$gene.id %>% uniqueN

v7 <- venn2(
  unique(gtex.all.eqtl$gene.id),
  unique(lookup.mqtl$gene)
)

# remove entries with genes I dont have
lookup.mqtl <- lookup.mqtl[gene %in% v7$q1]
lookup.mqtl[,index:=1:.N]

# reduce gtex data set
gtex.all.eqtl <- gtex.all.eqtl[gene.id %in% v7$q1, ]

# look for snp for each gene-metab pair
all.mqtl.gx.snp <- lapply(lookup.mqtl$index, function(i){
  
  my.gene <- lookup.mqtl[index==(i),gene]
  my.metabolite <- lookup.mqtl[index==(i),metabolite]
  my.snp.id <- tag.mqtl[metabolite == (my.metabolite), rs.id]
  my.snp <- tag.mqtl[metabolite == (my.metabolite), paste(
      paste0("chr",chr), pos, effect_allele,other_allele,"b38", sep="_")]
  my.snp.rev <- tag.mqtl[metabolite == (my.metabolite), paste(
      paste0("chr",chr), pos, other_allele,effect_allele,"b38", sep="_")]
  
  # check in gtex results
  res <- gtex.all.eqtl[variant_id %in% c(my.snp,my.snp.rev) & gene.id %in% (my.gene)]
  res[, ref.allele.like.meta.meta := ifelse(variant_id %in% my.snp, TRUE, FALSE)]
  res[,instrument.metab := (my.metabolite)]
  res[ref.allele.like.meta.meta==T,rs.id := my.snp.id[match(variant_id, my.snp)]]
  res[ref.allele.like.meta.meta==F,rs.id := my.snp.id[match(variant_id, my.snp.rev)]]
  return(res)
}) %>% rbindlist(fill = T,use.names = T)

fwrite(all.mqtl.gx.snp,dpx("allSNPmQTLInstrumentsForGx.csv", "res/"))

# overlap of genes in results
uniqueN(all.mqtl.gx.snp$gene_id)
uniqueN(genes$ensembl_gene_id_version)
vtmp <- venn2(all.mqtl.gx.snp$gene_id,genes$ensembl_gene_id_version)

# TODO: Check Ref-Allele ----

mqtl.metab <- tag.mqtl
mqtl.gx <- all.mqtl.gx.snp

# check overlap
table(mqtl.gx$ref.allele.like.meta.meta)
venn2(mqtl.metab[,unique(paste(chr,pos,other_allele,effect_allele,"b37",sep="_"))],
      mqtl.gx[,unique(variant_id)])

#=#################################################################=#
#===================================================================#
#=#################################################################=#

# Check total overlap ----
#' ## Check total overlap ----

# unify ID
eqtl.gx[, mr.id := paste(chr, pos, ref, alt, sep="_")]
eqtl.metab[, mr.id := ifelse(ref.allele.like.gtex==TRUE, 
                             paste(chr, pos, effect_allele, other_allele, sep="_"),
                             paste(chr, pos, other_allele, effect_allele, sep="_")
                             )
           ]
mqtl.metab[, mr.id := paste(chr, pos, effect_allele, other_allele, sep="_")]

var.split <- mqtl.gx$variant_id %>% strsplit(split="_")
mqtl.gx$chr <- lapply(var.split,`[[`,1)
mqtl.gx$pos <- lapply(var.split,`[[`,2)
mqtl.gx$ref <- lapply(var.split,`[[`,3)
mqtl.gx$alt <- lapply(var.split,`[[`,4)
mqtl.gx[, mr.id := ifelse(ref.allele.like.meta.meta==TRUE, 
                          paste(chr, pos, ref, alt, sep="_"),
                          paste(chr, pos, alt, ref, sep="_")
)
]

# check total overlap
v.total <- venn4(
  eqtl.gx$mr.id,
  eqtl.metab$mr.id,
  mqtl.metab$mr.id,
  mqtl.gx$mr.id
)

# check number of instruments
eqtl.gx[]
eqtl.metab$mr.id
mqtl.metab$mr.id
mqtl.gx$mr.id


#=#################################################################=#
#===================================================================#
#=#################################################################=#

#' # MR
# MR ----

# get everyone a variant.id | variant.id.rev based on match in data
eqtl.gx$variant_id_mr <- eqtl.gx$variant_id
eqtl.metab[, variant_id_mr := ifelse(ref.allele.like.gtex == TRUE, 
                                     paste(chr, pos, effect_allele, other_allele, "b37",sep="_"),
                                     paste(chr, pos, other_allele, effect_allele, "b37",sep="_")
)
]

# mqtl.gx$ref.allele.like.meta.meta %>% table
mqtl.metab[, variant_id_mr := paste(chr, pos, effect_allele, other_allele, "b37",sep="_")]
variant.split <- strsplit(mqtl.gx$variant_id, "_")
mqtl.gx[, variant_id_mr := ifelse(
  ref.allele.like.meta.meta == TRUE, 
  paste(
    lapply(variant.split,`[[`,1),
    lapply(variant.split,`[[`,2),
    lapply(variant.split,`[[`,3),
    lapply(variant.split,`[[`,4),
    lapply(variant.split,`[[`,5),
    sep="_"
  ),paste(
    lapply(variant.split,`[[`,1),
    lapply(variant.split,`[[`,2),
    lapply(variant.split,`[[`,4),
    lapply(variant.split,`[[`,3),
    lapply(variant.split,`[[`,5),
    sep="_"
  )
)]

# check on correct thingy
venn2(mqtl.metab$rs.id, mqtl.gx$rs.id)
venn2(mqtl.metab$variant_id_mr, mqtl.gx$variant_id_mr)

venn2(eqtl.gx$rs_id_dbSNP147_GRCh37p13, unlist(lapply(strsplit(eqtl.metab$markername,":"),`[[`,1)))
venn2(eqtl.gx$variant_id_mr, eqtl.metab$variant_id_mr)

# check if maf>1% & info>0.8 in (backwards) snp
eqtl.gx <- eqtl.gx[maf>=0.01] # TODO maf filer ok? ----
eqtl.metab <- eqtl.metab[maf>=0.01 & info >= 0.8]
# eqtl.metab$reason4exclusion_topl %>% table
# eqtl.metab$reason4exclusion_xtrct %>% table
mqtl.gx <- mqtl.gx[maf>=0.01]

# reverse beta for snp with reverese effect all
# TODO: decide on REM/FEM ----
eqtl.metab[, mr.beta := ifelse(ref.allele.like.gtex == TRUE, beta_expect, -1 * beta_expect)]
mqtl.gx[, mr.beta := ifelse(ref.allele.like.meta.meta == TRUE, slope, -1 * slope)]

# use this file to do all the checks
lookup[,index := 1:.N]

mr.snps <- mclapply(lookup$index, function(i){
  
  # get metab/gx pair
  my.probe <- lookup[index==(i),gx.probe]
  my.gene <- lookup[index==(i),gene]
  my.metab <- lookup[index==(i),metabolite]
  
  # get SNP from eqtl
  my.eqtl.gx.snp <- eqtl.gx[gene_name == (my.gene), variant_id_mr]
  my.eqtl.metab.snp <- eqtl.metab[pheno == (my.metab) & instrument.gene == (my.gene), variant_id_mr]
  eqtl.instruments <- intersect(my.eqtl.gx.snp, my.eqtl.metab.snp)
  
  # get SNP from mqtl
  my.mqtl.metab.snp <- mqtl.metab[metabolite == (my.metab), variant_id_mr]
  my.mqtl.gx.snp <- mqtl.gx[gene.id == (my.gene) & instrument.metab == (my.metab), variant_id_mr]
  mqtl.instruments <- intersect(my.mqtl.gx.snp, my.mqtl.metab.snp)
  
  mr.instruments <- rbindlist(
    list(
      data.table(
        mr.direction = "gx.to.metab", # gx->metab
        probe = my.probe,
        gene = my.gene,
        metabolite = my.metabolite,
        instrument = eqtl.instruments
      ),
      res.mqtl <- data.table(
        mr.direction = "metab.to.gx", # gx<-metab
        probe = my.probe,
        gene = my.gene,
        metabolite = my.metabolite,
        instrument = mqtl.instruments
      )
    ))
      
  return(mr.instruments)
},mc.cores = 20,mc.cleanup = T) %>% rbindlist  

  # Check ld<0.1 between Instruments in in both directions
  
  # MR gx->metab
  
  # MR gx<-metab
  
venn4(
  mqtl.metab$variant_id_mr,
  mqtl.gx$variant_id_mr,
  eqtl.gx$variant_id_mr,
  eqtl.metab$variant_id_mr
)

mr.snps.c <- dcast(mr.snps[!is.na(instrument)], 
      probe + gene + metabolite ~ mr.direction, 
      value.var = "instrument")
mr.snps.c[gx.to.metab > 0 & metab.to.gx > 0]
mr.snps[gene == "SLC22A16",]


# MendelianRandomization::mr_ivw()

# Instrumente für Metabolit --> GX aus der Meta-Meta GWAS und Effekte dieser Instrumente in Holger eQTL Datensatz nachschlagen
# Instrumente für GX --> Metabolit aus dem eQTL Datensatz, und Effekte dieser Instrumente auf Metabolite in der step10 nachschlagen



##########################

source("/net/ifs2//san_projekte//projekte/genstat/07_programme/rtools/1409_priority_pruning/priorityPruningGenedoses_190508.R")
readLines("/net/ifs2//san_projekte//projekte/genstat/07_programme/rtools/1409_priority_pruning/priorityPruningGenedoses_190508.R")



# toplist_maxR2prune=0.1

# Oder Plink:
# Als Ref 1000Genomes P3
"R:\genstat\01_daten\1503_1000genomes_phase3_vs5_okt14_IMPUTE\plinkformat\results"
"/net/ifs2/san_projekte/projekte/genstat/01_daten/1503_1000genomes_phase3_vs5_okt14_IMPUTE/plinkformat/results/s05_1_503EUR_maf0_5_1kg_phase3v5_from_imputeRef_autoANDchrXnonPAR"

# https://www.cog-genomics.org/plink/2.0/ld#indep

datafn<-"--bfile ../data/1KG_PCA"
snplistfn<-"--extract results_PCA/mySnps_filtered.txt"
samplelistfn<-"--keep results_PCA/mySamples.txt"
ld<-"--indep-pairphase 50 5 0.2"
outfn<-"--out results_PCA/pruning_filter"
call1<-paste(plink_call,datafn,snplistfn,samplelistfn, ld,outfn, sep=" ")
system(call1)
plink_call<-"/net/ifs2/san_projekte/projekte/genstat/07_programme/plink1.9/vs180103_stable_beta_51/unix64/plink"
