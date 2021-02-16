#' ---
#' title: "Enrichment of metabolite pathways for master regulator genes"
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
  "ggplot2",
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
source(here("../functions/hyper_test.R"))
source("/net/ifs2/san_projekte/projekte/genstat/07_programme/rtools/1409_pathway_enrichment/hypergeometricEnrichment_KEGGneu_GO_Reactome_DOSE_viaEntrez_180608.R")
an <- all_annotation_tables(mount.point = "/net/ifs2/san_projekte/projekte/genstat/")
option_setup()

#+ include=F
#==========#
#   DONE   #
#==========#

#' I evaluated annotated data in the previous step.

nf <- newest_file(look_for = "STEP4_meta_resultsAnnotated",subfolder = "res", print_full = TRUE)
dat <- fread(nf)

# get the annotated IDs
nf <- newest_file(look_for = "metabOverview", subfolder = "obj", print_full = TRUE)
st1 <- fread(nf)

# load metabolite data base
# current.msetlib -> these are the pathways and each metabolite in it 
load("../170829_ConfounderAnalysis/190123_pathwayenrichExample/pathway.rda")
setDT(current.msetlib)

current.msetlib %>% str

# pdb = pathway data base
pdb <- current.msetlib

# basierende auf MetPa https://academic.oup.com/bioinformatics/article/26/18/2342/208464
# compund data base -> This has the matched HMDB etc IDs!
cdb <- readRDS("../170829_ConfounderAnalysis/190123_pathwayenrichExample/compound_db.rds")
setDT(cdb)
cdb %>% str

# Load Supp table for most recent metabolite hmdb/pubchem IDs
ids <- read_excel2(fn = "../170829_ConfounderAnalysis/paper/submission/revised_molecular_metabolism/Supplementary_tables.xlsx", sheet = 2, skip = 1)
ids$Abbreviation <- gsub(pattern = "*",replacement = "",x = ids$Abbreviation, fixed = TRUE)
ids$Abbreviation[ids$Abbreviation == "AC-total ()"] <- "AC-total"

#' # Find the overlap of the ids in cbd and pdb

# full overlap?
all(ids$Abbreviation %in% st1$abbr) # -> TRUE

m1 <- match(st1$abbr, ids$Abbreviation)

# hmdb2 is the one from the supp table
st1 <- cbind(st1, ids[m1, .(PubChem)])

# use only identified compounds
all.ids <- st1[, .(rid, abbr, full.name, hmdb, pubchem = PubChem, metab.class, metab.class.q, metab.super.path)]
all.ids$hmdb[all.ids$hmdb == ""] <- NA
all.ids$hmdb[all.ids$hmdb == "-"] <- NA
all.ids$pubchem[all.ids$pubchem=="-"] <- NA

# easy fix
all.ids$pubchem[is.na(all.ids$pubchem)] <- "---"
all.ids$hmdb[is.na(all.ids$hmdb)] <- "---"

#' Everything except compound 37 looks OK. Is it GABA or Abscisic acid. The latter is a plant hormone, so it must be GABA! I changed that before running the search!

# change wrong pubchem ID
all.ids$pubchem[all.ids$pubchem == 5702609] <- 119

# map all metabolites based on the data-base
found.hmdb <- lapply(1:nrow(all.ids), function(i){
  
  if(grepl(pattern = "/", x = all.ids$hmdb[i], fixed = TRUE)){
    
    # there are sometimes multiple hmdb ids per metabolite (isoform etc), split these and search for
    pat <- gsub(pattern = "/", replacement = "|", x = all.ids$hmdb[i], fixed = TRUE)
    x1 <- grep(pattern = pat, x = cdb$hmdb_id, value = TRUE)
    
    pat2 <- gsub(pattern = "/", replacement = "$|^", x = all.ids$pubchem[i], fixed = TRUE)
    pat2 <- paste0("^", pat2, "$")
    x2 <- grep(pattern = pat2, x = cdb$pubchem_id, value = TRUE)
    
    # get metabolite
    m1 <- cdb[hmdb_id %in% x1, name]
    m2 <- cdb[pubchem_id %in% x2, name]
    
    # get additional ids
    # id1 <- cdb[hmdb_id %in% x1, .(name, chebi_id,kegg_id,metlin_id)]
    # id2 <- cdb[pubchem_id %in% x2, .(name, chebi_id,kegg_id,metlin_id)]
    
    # id <- rbindlist(list(id1,id2))
    x <- c(m1,m2)
    # x <- c(x1,x2)
    
  } else{  
    
    # in case only one compund matches
    x1 <- grep(pattern = all.ids$hmdb[i], x = cdb$hmdb_id, value = TRUE, fixed = TRUE)
    
    x2 <- cdb$pubchem_id[cdb$pubchem_id %in% all.ids$pubchem[i]]
    # x2 <- grep(pattern = paste0("^", all.ids$pubchem[i], "$"), x = cdb$pubchem_id, value = TRUE)
    
    # get metabolite
    m1 <- cdb[hmdb_id %in% x1, name]
    m2 <- cdb[pubchem_id %in% x2, name]
    
    # get additional ids
    # id1 <- cdb[hmdb_id %in% x1, .(chebi_id,kegg_id,metlin_id)]
    # id2 <- cdb[pubchem_id %in% x2, .(chebi_id,kegg_id,metlin_id)]
    
    # id <- c(id1,id2)
    x <- c(m1,m2)
    # x <- c(x1, x2)
  }
  
  return(x)
})

# get the unique compunds
unique.found.hmdb <- lapply(found.hmdb, unique)
unique.found.hmdb <- lapply(unique.found.hmdb, function(x){
  paste0(x, collapse=", ")
})

# Correct mapping according to
all.ids[, mapped.metabs := unlist(unique.found.hmdb)]
all.ids[mapped.metabs == "", mapped.metabs := NA]



# expand quotients
source("../../1703_ge_metab_a1_b3_sorbs/functions/quotient_calc_table.R")
qt <- quotient_calc_table()
qt <- melt(qt,
           id.vars = "Quotient", 
           measure.vars = c("Numerator1","Numerator2","Denominator1","Denominator2"), 
           value.name = "metabolite",variable.factor = FALSE)
qt <- qt[metabolite != "0"]

# amm = all mapped metabolites
amm <- unlist(strsplit(x = unlist(unique.found.hmdb), split = ", "))

#' # Enrichment of metabolite pathways in master-regulators

#' ## Function for metabolite pathway enrichment

# get the metabolites per gene
mrg.num <- dat[symbol_INGENUITY != "" & significant == TRUE, .(sig.metabs = uniqueN(metabolite)),by=symbol_INGENUITY][order(sig.metabs, decreasing = TRUE), ]

#' First check each gene that has >=10 significantly associated metabolites

# gtt = genes to test
# mrg.num[1:50, range(sig.metabs)]
# gtt <- mrg.num[, symbol_INGENUITY][1:50]
gtt <- mrg.num[sig.metabs>=5, symbol_INGENUITY]

# what are the quantiles
mrg.num$sig.metabs %>% quantile(.,probs=c(0,.25,.4,.75,.95,1))

# get the list with drawn metabolites per gene
mrg <- dat[symbol_INGENUITY %in% gtt & significant == TRUE, .(sig.metabs = unique(metabolite)),by=symbol_INGENUITY]

# map the quotient parts to the metabolite names
# only take quotients in the table
v.q <- venn2(qt$Quotient, mrg$sig.metabs)
qt <- qt[Quotient %in% v.q$q1, ]

# match mrg info to quotients
m.q <- match(qt[,Quotient], mrg$sig.metabs )
qt <- cbind(qt, mrg[m.q, ])
qt.remove <- unique(qt$Quotient)
qt <- qt[,.(symbol_INGENUITY, sig.metabs = metabolite)]

mrg <- rbindlist(list(mrg[sig.metabs %nin% qt.remove, ], qt))

# map the correct metabolites to this list
m1 <- match(mrg$sig.metabs, all.ids$rid)
mrg$compound <- all.ids[m1, mapped.metabs]

# remove all non-mapped metabolites
mrg <- mrg[!is.na(compound), ]

# apm = all pathway metabolites
apm <- unlist(strsplit(pdb$member, split = "; "))

# should I include q3?
v1 <- venn2(apm, amm)

# this should be the background
# metabolite.bg <- c(v1$q1, v1$q2)
metabolite.bg <- c(v1$q1)

# gpe = gene-wise pathway enrichment
gpe <- lapply(unique(mrg$symbol_INGENUITY), function(i){
  
  # i <- unique(mrg$symbol_INGENUITY)[1]
  gezogene <- mrg[symbol_INGENUITY == (i), compound]
  
  # ppe = per pathway enrichment
  ppe <- lapply(1:nrow(pdb), function(ii){
    
    # ii <- 1
    pathway.metabolites <- pdb[ii, unlist(strsplit(member, "; "))]
    
    # use holgers function
    ht <- suppressMessages(
      hyper_test(
        pathway_gene = pathway.metabolites,
        gezogene = gezogene,
        alle_gene = metabolite.bg
      )
    )
    
    # create results 
    one.pathway <- data.table(
        id = pdb[ii, id],
        name = pdb[ii, name],
        in.pathway = pdb[ii, member],
        in.drawn = paste(gezogene, collapse = "; "),
        in.drawn = ht$in_gezogen,
        in.bg = ht$in_bk,
        enrichment = ht$enrichment,
        or = ht$or,
        or.lower = ht$or_lower,
        or.upper = ht$or_upper,
        pval.hyper = ht$pval,
        pva.fisher = ht$pval_fisher
      )
      
    
    return(one.pathway)
    
  })
  
  one.gene.all.pathways <- rbindlist(ppe)
  one.gene.all.pathways[, gene := (i)]
  return(one.gene.all.pathways)
})

# consolidate all enrichments
all.res <- rbindlist(gpe)

all_res_hfdr <- addHierarchFDR(pvalues = all.res$pval.hyper,categs = all.res$gene, quiet = F)
all_res <- cbind(all.res, all_res_hfdr[,.(fdr_level1,fdr_level2, hierarch_fdr5proz)])
all_res[,table(hierarch_fdr5proz)]

all.res[, p.hyper.fdr := p.adjust(pval.hyper, method = "BH")]
all.res[, .SD[pval.hyper==min(pval.hyper)], by = gene][, head(.SD,1),by = gene][order(p.hyper.fdr,decreasing = F)]
all.res[order(p.hyper.fdr)]

# prepare table with results

# annotate top 50 association hubs + those that are not

# hubs
gtt <- mrg.num[, symbol_INGENUITY][1:50]

# sub-hubs
gtt <- mrg.num[sig.metabs>=5, symbol_INGENUITY]




save_csv_carl(file = all.res, 
              file_name = "masterRegulatorEnrichment",
              subfolder = "obj", 
              quote = TRUE)

# check significancies
all.res[p.hyper.fdr<=0.05, unique(gene)]
all.res[p.hyper.fdr<=0.05, uniqueN(gene)]
all.res[p.hyper.fdr<=0.05, table(name)]
all.res[p.hyper.fdr<=0.05 & (name=="Phenylalanine and Tyrosine Metabolism")]

sig.genes <- all.res[p.hyper.fdr<=0.05, gene]

mrg.num[,max(sig.metabs)]
mrg.num[symbol_INGENUITY %in% (sig.genes), range(sig.metabs)]


#'  Based on Holgers Hypergeometric Test function and the pathway info used also in the MetaboAnalyst.

#'  # END

devtools::session_info()