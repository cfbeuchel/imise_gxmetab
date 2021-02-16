#' ---
#' title: "GTEx v8 LiftOver"
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
# root.dir <- paste0(bp, "02_projekte/1703_ge_metab_a1_b3_sorbs/1805_gxMetaboliteAssociation/")
# setwd(root.dir)

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
source(paste0(bp, "02_projekte/1703_ge_metab_a1_b3_sorbs/functions/option_setup.R"))
setDTthreads(8)
option_setup()


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
