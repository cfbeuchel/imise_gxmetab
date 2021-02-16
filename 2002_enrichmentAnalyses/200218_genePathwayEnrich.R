#' ---
#' title: "Genetic pathway enrichment"
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
  "magrittr",
  "ggplot2",
  "clusterProfiler",
  "ReactomePA",
  "DOSE",
  "org.Hs.eg.db"
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
source(here("../functions/all_annotation_tables.R"))
source(here("../functions/option_setup.R"))
an <- all_annotation_tables(mount.point = "/net/ifs2/san_projekte/projekte/genstat/")
setDTthreads(1)
option_setup()

#+ include=F
#==========#
#   DONE   #
#==========#

# Load data ----
#' # Load data

# sobel test results
res.sobel <- fread(newest_file("mediationStatisticsAnnotatedJoint","res",print_full = T))

# mediation results annoted for the groups based the quotient of direct and raw effect
res.mediation <- fread(newest_file("mediatonGroupList","res",print_full = T))

# meta-analysis results
res.meta <- fread(newest_file("STEP4_meta_resultsAnnotated","res",print_full = T))

# gx mapping
ilmn.map <- fread(paste0(bp, "/07_programme/rtools/1807_gx_tools/mappingHT12/results/102_8_remappingHT12_INGENUITY_positionalinfo_ilmnCentric.txt"))
entrez.map <- fread(paste0(bp, "/07_programme/rtools/1807_gx_tools/mappingHT12/results/102_8_remappingHT12_INGENUITY_positionalinfo_entrezCentric.txt"))

# probe annotations
probe.adult  <- fread(paste0(bp, "02_projekte/1705_ge_lifea1v2/180719_a1_gePreProForMetab/tosend_NoBMI_NoDiabetes/s401_1_probeannot_HT12v4.txt"))
probe.heart  <- fread(paste0(bp, "02_projekte/1102_ge_lifeb3/07_prepro_ami_allb3_for_metab/tosend_NoBMI_NoDiabetes/probeannot_HT12v4.txt"))
probe.sorb  <- fread(paste0(bp, "02_projekte/1706_ge_sorben_v2/181009_sorb_gePreProForMetab/tosend_NoBMI_NoDiabetes/s401_1_probeannot_HT12v4.txt"))
probe <- rbindlist(list(probe.adult,probe.heart,probe.sorb),use.names = T,fill = T)
probe <- probe[perfectprobe_PBMC == TRUE | 
                 perfectprobe_ami == TRUE | 
                 perfectprobe_heart == TRUE | 
                 perfectprobe_all == TRUE | 
                 perfectprobe_allsubgroups == TRUE, .(Entrez_Gene_ID, ilmn, ILMN_Gene, Symbol)]

# any duplicated entrez?
probe[,uniqueN(Entrez_Gene_ID),by=ilmn][V1>1, ]

# just use unique IDs
probe <- probe[, .SD[1], by = ilmn]

grep(pattern = ",",fixed = TRUE,x =  probe$Entrez_Gene_ID,value = T)
grep(pattern = ",",fixed = TRUE,x =  probe$Entrez_Gene_ID,value = T) %>% length

# match entrez to meta results
entrez.qc <- unique(entrez.map$ilmn_mapping_hg19v2) %>% sort

# check gene overlap
v1 <- venn2(
  entrez.map[ilmn_mapping_hg19v2 %in% entrez.qc[1:8], ilmn],
  res.meta$markerID
)
v2 <- venn2(
  res.meta$symbol_INGENUITY,entrez.map[ilmn_mapping_hg19v2 %in% entrez.qc[1:4], 
                                       symbol_INGENUITY]
)

# Matching ----
#' ## Matching

res.meta$entrez.id <- NULL
m1 <- match(res.meta$symbol_INGENUITY, entrez.map[ilmn_mapping_hg19v2 %in% entrez.qc[1:8],symbol_INGENUITY])
res.meta <- cbind(res.meta, entrez.map[ilmn_mapping_hg19v2 %in% entrez.qc[1:8],][m1, .(entrez.id = ilmn_entrezID_INGENUITY)])
# m1 <- match(res.meta$markerID, probe$ilmn)
# res.meta <- cbind(res.meta, probe[(m1), .(entrez.id = Entrez_Gene_ID)])

# match to sobel results
m2 <- match(res.sobel$symbol_INGENUITY, entrez.map[ilmn_mapping_hg19v2 %in% entrez.qc[1:8],symbol_INGENUITY])
res.sobel <- cbind(res.sobel, entrez.map[ilmn_mapping_hg19v2 %in% entrez.qc[1:8],][m2, .(entrez.id = ilmn_entrezID_INGENUITY)])

# match to sobel group results
m3 <- match(res.mediation$gene, entrez.map[ilmn_mapping_hg19v2 %in% entrez.qc[1:8],symbol_INGENUITY])
res.mediation <- cbind(res.mediation, entrez.map[ilmn_mapping_hg19v2 %in% entrez.qc[1:8],][m3, .(entrez.id = ilmn_entrezID_INGENUITY)])

#=================================================================================#
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#=================================================================================#

# Pathway enrichment ----
#' # Pathway enrichment

# Enrichment of Gx-Metabolite Associations ----
#' ## Enrichment of Gx-Metabolite Associations

# get fg
assoc.fg <- res.meta[!is.na(entrez.id) & significant==TRUE, 
                              unique(entrez.id), 
                              by=metabolite]

# seperate rows with multiple IDs
assoc.fg <-sapply(unique(assoc.fg$metabolite), function(x){
  res <- assoc.fg[metabolite==(x), unique(V1)]
  res <- unlist(strsplit(x = res, split = ", ", fixed = T))
  res <- list(res)
  names(res) <- x
  return(res)
},USE.NAMES=FALSE)

# get bg - all IDs
assoc.bg <- res.meta[!is.na(entrez.id), unique(entrez.id)]
assoc.bg <- as.character(unlist(strsplit(x = assoc.bg, split = ", ", fixed = T)))

# enrichment analysis
assoc.enrich.go <- compareCluster(
  geneCluster = assoc.fg, 
  fun = "enrichGO", 
  OrgDb = org.Hs.eg.db,
  universe = assoc.bg, 
  minGSSize = 2,
  pvalueCutoff = 0.2, 
  qvalueCutoff = 0.2
)

assoc.enrich.go <- simplify(assoc.enrich.go)
res.go <- as.data.table(assoc.enrich.go@compareClusterResult)
res.go[,mediation.group := (i)]
res.go[,data.base := "GO"]

# viz
dotplot(assoc.enrich.go)

# kegg
enrich.kegg <- compareCluster(
  geneCluster = assoc.fg, 
  fun = "enrichKEGG", 
  pvalueCutoff = 0.3, 
  universe = assoc.bg, 
  minGSSize = 2,
  use_internal_data=F,
  qvalueCutoff = 0.3
)

# save results
res.kegg <- as.data.table(enrich.kegg@compareClusterResult)
res.kegg[,mediation.group := (i)]
res.kegg[,data.base := "KEGG"]
res.kegg[!is.na(qvalue) & qvalue<=0.05][order(qvalue), head(.SD,10)]

# viz
dotplot(enrich.kegg)







# Enrichment of Mediation Quotient Groups ----
#' ## Enrichment of Mediation Quotient Groups

#' Only possible for log-BMI

# loop through mediation groups
mediation.groups <- unique(
  subset(res.mediation, 
         pheno == "log.bmi" & 
           mediation.group != "Other")$mediation.group)
res.enrich.groups <- mclapply(mediation.groups,
  function(i){
    
    # define background
    bg.mediation <- res.mediation[!is.na(entrez.id) & 
                                    pheno == "log.bmi", 
                                  unique(entrez.id)]
    bg.mediation <- as.character(unlist(strsplit(x = bg.mediation, split = ", ", fixed = T)))
    
    # define foreground
    fg.mediation <- res.mediation[!is.na(entrez.id) & 
                                    pheno == "log.bmi" &
                                    mediation.group == (i), 
                                  entrez.id, 
                                  by=metabolite]
    fg.mediation <-sapply(unique(fg.mediation$metabolite), function(x){
      res <- fg.mediation[metabolite==(x), unique(entrez.id)]
      res <- unlist(strsplit(x = res, split = ", ", fixed = T))
      res <- list(res)
      names(res) <- x
      return(res)
    },USE.NAMES=FALSE)
    
    #' ### GO Data Base
    # GO Data Base ----
    
    enrich.go <- compareCluster(
      geneCluster = fg.mediation, 
      fun = "enrichGO", 
      OrgDb = org.Hs.eg.db,
      universe = bg.mediation, 
      minGSSize = 2,
      pvalueCutoff = 0.2, 
      qvalueCutoff = 0.2
    )
    
    # save results
    enrich.go <- simplify(enrich.go)
    res.go <- as.data.table(enrich.go@compareClusterResult)
    res.go[,mediation.group := (i)]
    res.go[,data.base := "GO"]
    
    # viz
    dotplot(enrich.go)
    
    #' ### KEGG Data Base
    # KEGG Data Base ----
    
    enrich.kegg <- compareCluster(
      geneCluster = fg.mediation, 
      fun = "enrichKEGG", 
      pvalueCutoff = 0.3, 
      universe = bg.mediation, 
      minGSSize = 2,
      use_internal_data=F,
      qvalueCutoff = 0.3
    )
    
    # save results
    res.kegg <- as.data.table(enrich.kegg@compareClusterResult)
    res.kegg[,mediation.group := (i)]
    res.kegg[,data.base := "KEGG"]
    
    # viz
    dotplot(enrich.kegg)
    
    #' ### REACTOME Data Base
    # REACTOME Data Base ----
    
    enrich.path <- compareCluster(
      geneCluster = fg.mediation, 
      fun = "enrichPathway", 
      universe = bg.mediation)
    
    # save results
    res.reactome <- as.data.table(enrich.path@compareClusterResult)
    res.reactome[,mediation.group := (i)]
    res.reactome[,data.base := "REACTOME"]
    
    # plot
    dotplot(enrich.path)
    
    #' ### DO Data Base
    # DO Data Base ----
    
    enrich.do <- compareCluster(
      geneCluster = fg.mediation, 
      fun = "enrichDO", 
      universe = bg.mediation
    )
    
    # save results
    res.do <- as.data.table(enrich.path@compareClusterResult)
    res.do[,mediation.group := (i)]
    res.do[,data.base := "DO"]
    
    # plot
    dotplot(enrich.do,includeAll=TRUE)
    
    # Join the results and return
    res <- rbindlist(list(
      res.go,
      res.kegg,
      res.reactome,
      res.do
    ))
    
    # bundle with the single results
    res <- list(joint.results = res,
                cluster.profiler.objects = list(
                  enrich.do,
                  enrich.go,
                  enrich.kegg,
                  enrich.path
                ))
    names(res)[2] <- i
    
    # return results
    return(res)
    
  },mc.cores = 3)

# save results
res.enrich.tables <- rbindlist(lapply(res.enrich.groups, `[[`, 1))
res.enrich.objects <- lapply(res.enrich.groups, `[[`, 2)

#' ## Enrichment of joint group gx->metab & metab->gx
# Enrichment of joint group gx->metab & metab->gx ----

# define background
bg.mediation <- res.mediation[!is.na(entrez.id) & 
                                pheno == "log.bmi", 
                              unique(entrez.id)]
bg.mediation <- as.character(unlist(strsplit(x = bg.mediation, split = ", ", fixed = T)))

# define foreground
fg.mediation <- res.mediation[!is.na(entrez.id) & 
                                pheno == "log.bmi" &
                                mediation.group %in% c(
                                  "Metabolite mediated gene expression effects",
                                  "Gene expression mediated metabolite effects"), 
                              entrez.id, 
                              by=metabolite]
i<-paste0(c(
  "Metabolite mediated gene expression effects",
  "Gene expression mediated metabolite effects"), collapse = " and ")
fg.mediation <-sapply(unique(fg.mediation$metabolite), function(x){
  res <- fg.mediation[metabolite==(x), unique(entrez.id)]
  res <- unlist(strsplit(x = res, split = ", ", fixed = T))
  res <- list(res)
  names(res) <- x
  return(res)
},USE.NAMES=FALSE)

#' ### GO Data Base
# GO Data Base ----

enrich.go <- compareCluster(
  geneCluster = fg.mediation, 
  fun = "enrichGO", 
  OrgDb = org.Hs.eg.db,
  universe = bg.mediation, 
  minGSSize = 2,
  pvalueCutoff = 0.2, 
  qvalueCutoff = 0.2
)

# save results
enrich.go <- simplify(enrich.go)
res.go <- as.data.table(enrich.go@compareClusterResult)
res.go[,mediation.group := (i)]
res.go[,data.base := "GO"]

# viz
# dotplot(enrich.go)

#' ### KEGG Data Base
# KEGG Data Base ----

enrich.kegg <- compareCluster(
  geneCluster = fg.mediation, 
  fun = "enrichKEGG", 
  pvalueCutoff = 0.3, 
  universe = bg.mediation, 
  minGSSize = 2,
  use_internal_data=F,
  qvalueCutoff = 0.3
)

# save results
res.kegg <- as.data.table(enrich.kegg@compareClusterResult)
res.kegg[,mediation.group := (i)]
res.kegg[,data.base := "KEGG"]

# viz
# dotplot(enrich.kegg)

#' ### REACTOME Data Base
# REACTOME Data Base ----

enrich.path <- compareCluster(
  geneCluster = fg.mediation, 
  fun = "enrichPathway", 
  universe = bg.mediation)

# save results
res.reactome <- as.data.table(enrich.path@compareClusterResult)
res.reactome[,mediation.group := (i)]
res.reactome[,data.base := "REACTOME"]

# plot
# dotplot(enrich.path)

#' ### DO Data Base
# DO Data Base ----

enrich.do <- compareCluster(
  geneCluster = fg.mediation, 
  fun = "enrichDO", 
  universe = bg.mediation
)

# save results
res.do <- as.data.table(enrich.path@compareClusterResult)
res.do[,mediation.group := (i)]
res.do[,data.base := "DO"]

# plot
# dotplot(enrich.do,includeAll=TRUE)

# Join the results and return
res <- rbindlist(list(
  res.go,
  res.kegg,
  res.reactome,
  res.do
))

# bind together with the single results
res.enrich.tables <- rbindlist(list(res.enrich.tables,
               res))
res.enrich.objects[[4]] <- list(
  enrich.do,
  enrich.go,
  enrich.kegg,
  enrich.path
)

# save seperately
fwrite(res.enrich.tables, dpx("mediationGroupEnrichResults", "res/"), sep="\t")
# res.enrich.tables <- fread(newest_file("mediationGroupEnrichResults","res",print_full = T))
save(res.enrich.objects, file = dpx("mediationGroupEnrichObjects.RData", "res/"))

# Explore Results ----
#' ### Explore Results

# add column enrichment: found/expected
# https://yulab-smu.github.io/clusterProfiler-book/chapter2.html
# https://www.biostars.org/p/220465/
# GeneRatio = k/n (Meine Gene im Pathway / Meine Gene in allen Genen)
# BgRatio = M/N (Alle gene im Pathway / Alle Gene)
res.enrich.tables$enrichment <- NULL
fg.ratio <- strsplit(res.enrich.tables$GeneRatio, split = "/", fixed = T)
fg.ratio <- unlist(lapply(fg.ratio,function(i){
  as.numeric(i[1])/as.numeric(i[2])
}))
bg.ratio <- strsplit(res.enrich.tables$BgRatio, split = "/", fixed = T)
bg.ratio <- unlist(lapply(bg.ratio,function(i){
  as.numeric(i[1])/as.numeric(i[2])
}))
res.enrich.tables$enrichment <- fg.ratio/bg.ratio

# Add enriched Genes ----  
  
#' Add the ID of the enriched gene in order of mediation effect strength

# add a sum of the quotient for joint consideration
res.mediation[, quot.sum := `quot.direct.vs.marg_pheno~gx(+metab)` + `quot.direct.vs.marg_pheno~metab(+gx)`]


gene.list <- sapply(1:nrow(res.enrich.tables), function(index){
  
  i <- res.enrich.tables[(index), ]
  metab <- i$Cluster
  genes <- unlist(strsplit(i$geneID,"/"))
  group <- i$mediation.group
  
  # set the correct joint groups for filtering
  if (group == "Metabolite mediated gene expression effects and Gene expression mediated metabolite effects") {
    group <- c("Gene expression mediated metabolite effects", 
               "Metabolite mediated gene expression effects")
    
    # sortiere nach stärke der mediation, also nach kleinem quotienten = viel Effekt wird mediiert
    setorder(res.mediation, quot.sum)
    
  } else if(group == "Gene expression mediated metabolite effects"){
    
    setorder(res.mediation, `quot.direct.vs.marg_pheno~metab(+gx)`)
    
  } else if(group == "Metabolite mediated gene expression effects"){
    
    setorder(res.mediation, `quot.direct.vs.marg_pheno~gx(+metab)`)
    
  } else if(group == "Bi-directional mediation"){
    
    setorder(res.mediation, quot.sum)
    
  }
  
  # match the entrez ID to the Gene of interest
  gene.list <- res.mediation[pheno=="log.bmi"][
    mediation.group %in% group & metabolite %in% metab, 
    .(gene, entrez.id, gx.probe)]
  
  # split all comma-seperated genes into rows
  # unlist(strsplit(gene.list$entrez.id, split = ", ",fixed = T))
  gene.list$entrez.id.greppable <- sapply(strsplit(gene.list$entrez.id, split = ", ",fixed = T), function(ii){
    ii <- ifelse(length(ii)==1, paste0("^", ii, "$"),paste0(paste0("^",ii,"$"), collapse="|"))
   return(ii)
  })
  
  # match the gene name to the entrez ID of the found genes in the pathway
  for(iii in 1:nrow(gene.list)){
    tmp <- gene.list[iii,]
    matches <- grep(tmp$entrez.id.greppable, genes)
    names(genes)[matches] <- tmp$gene
  }
  
  # names are still in order of the entry in the enrichment result file.
  # order them by mediation effect size
  
  # reference order
  genes <- names(genes)
  genes <- factor(genes, ordered = T, levels = unique(gene.list$gene)[unique(gene.list$gene) %in% genes])
  genes <- sort(genes)
  return(paste(genes, collapse=","))
  
})

# add new column
res.enrich.tables$enriched.genes <- gene.list

# save for inspection
WriteXLS_hk(res.enrich.tables, ExcelFileName = dpx("genePathwayEnrichment.xls","res/"),SheetNames = "1")
# res.enrich.tables <- readxl::read_xls(newest_file("genePathwayEnrichment.xls","res",print_full = T))
# setDT(res.enrich.tables)

# check the results
res.enrich.tables[order(p.adjust,decreasing = FALSE), head(.SD,5), by= .(mediation.group,data.base)] %>% as.data.frame
setorder(res.enrich.tables, p.adjust, -Count)
res.enrich.tables[data.base=="GO" & p.adjust<=0.05, ]
res.enrich.tables[data.base=="GO" & p.adjust<=0.05, .N, by= Description][order(N,decreasing = T)] %>% head(10)


#' ## Obesity related genes 
# Obesity related genes ----

res.meta[symbol_INGENUITY %in% c("GLUT1")] # none

# relevant obesity gene: https://link.springer.com/article/10.1007/s00439-020-02155-1
obesity.assocs1 <- res.meta[significant==TRUE & 
                              symbol_INGENUITY %in% c("PRKCQ","IFNG","SOX6","MARK3","DNM3","CDKAL1","MPP7"), 
                            .(symbol_INGENUITY,metabolite)]

# Genetics of obesity - DOI: 10.1007/s12291-015-0541-x
# https://www.nature.com/articles/nrg2594
obesity.assocs2 <- res.meta[significant==TRUE & symbol_INGENUITY %in% c(
  "LINC01122", "NLRC3-ADCY9", "GPRC5B-GP2", "BDNF", "MC4R", "AGBL4-ELAVL4", 
  "ATP2A1-SBK1", "TCF7L2", "GIPR", "IRS1", "FOXO3", "ASB4",
  "RPTOR", "NPC1", "CREB1", "FAM57B", "APOBR", "HSD17B12", "PTBP2",
  "ELAVL4", "CELF1", "RALYL", "MAP2K5", "MAPK3", "FAIM2", "PARK2",
  "OLFM4", "FTO", "TMEM18", "KCTD15", "GNPDA2", "SH2B1", "MTCH2",
  "NEGR1","SH2B1","ADCY3", "CTNNBL1", "PCSK1", "POMC", "BDNF","ADIPOQ",
  "LEP","LEPR","INSIG2","PPARG", "GPRC5B"
), .(symbol_INGENUITY,metabolite)]

# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3014648/
obesity.assocs3 <- res.meta[significant==TRUE & symbol_INGENUITY %in% c(
  "ETV5","TNNI3K","SEC16B","TFAP2B","NRXN3","RBJ",
  "MAP2K5","LBXCOR1","QPCTL","GIPR",
  "TNNI3K","FLJ35779","HMGCR","LRRN6C",
  "FANCL","CADM2","PRKD1","LRP1B",
  "PTBP2","MTIF3","ZNF608","RPL27A",
  "TUB","NUDT3","HMGA1",
  "CRAT", #
  "ADCY3", # Gx
  "POMC", # Gx
  "IQCK", # Gx
  "GPRC5B", # Gx
  "TMEM160", # Gx
  "SLC39A8",  # Gx
  "GTF3A", # Gx
  "ZC3H4"),  # Gx
  .(symbol_INGENUITY,metabolite)]

# join all assocs related to obesity genes
obesity.assocs <- rbindlist(list(obesity.assocs1,obesity.assocs2, obesity.assocs3))
obesity.assocs[, joint := paste0(metabolite,symbol_INGENUITY)]

# check in mediation results
obesity.mediations <- res.sobel[is.mediation == "yes" & 
            fdr.sig == TRUE & 
            paste0(metabolite, symbol_INGENUITY) %in% obesity.assocs$joint, 
            .(type, pheno, metabolite, gene = symbol_INGENUITY,
              id = paste0(pheno, metabolite, symbol_INGENUITY))]

# overlap of mediations per phenotype
venn.mediation <- venn2(
  obesity.mediations[pheno=="diabetes.status.tri", paste0(metabolite, "__", gene)],
  obesity.mediations[pheno=="log.bmi", paste0(metabolite, "__", gene)]
)

# strong mediations classified by r2
res.mediation[paste0(pheno, metabolite, gene) %in% obesity.mediations$id & 
                !(mediation.group %in% c("Other", "None")), 
              .(pheno, metabolite, gene, mediation.group)][,table(gene)] %>% sort
res.mediation[paste0(pheno, metabolite, gene) %in% obesity.mediations$id & 
                !(mediation.group %in% c("Other", "None")), 
              .(pheno, metabolite, gene, mediation.group)][,table(metabolite)] %>% sort

#' ## CPT1 activity
#  CPT1 activity ----

res.mediation[metabolite %in% c("Q2","Q3"), table(mediation.group)]
res.sobel[metabolite %in% c("Q2","Q3")  & (fdr.sig==TRUE | hfdr.sigREM == TRUE), 
          .(symbol_INGENUITY, 
            metabolite, 
            pheno,
            type,
            beta.tau.tauAdj, 
            betaREM
          )][,table(sign(beta.tau.tauAdj))]
res.sobel[metabolite %in% c("Q2","Q3")  & (fdr.sig==TRUE | hfdr.sigREM == TRUE), 
          .(symbol_INGENUITY, 
            metabolite, 
            pheno,
            type,
            beta.tau.tauAdj, 
            betaREM
          )][,table(symbol_INGENUITY)] %>% sort(decreasing = T) %>% head(10)

res.mediation[metabolite %in% c("Q2","Q3"), .(beta)]
res.mediation[metabolite %in% c("Q2","Q3"), table(gene)] %>% sort(decreasing = T) %>% head(10)

# enrichment for Q2/Q3 mediation
res.enrich.tables[p.adjust<=0.05 & Cluster %in% c("Q2","Q3"), table(Description)] %>% sort

#' ## Metabolic felxibilty (C2)
#  Metabolic felxibilty (C2) ----

res.sobel[metabolite %in% c("C2")  & (fdr.sig==TRUE | hfdr.sigREM == TRUE), 
          .(symbol_INGENUITY, 
            metabolite, 
            pheno,
            type,
            beta.tau.tauAdj, 
            betaREM
          )][,table(symbol_INGENUITY)] %>% sort(decreasing = T) %>% head(10)

res.mediation[metabolite %in% c("C2"), table(mediation.group)] %>% sort(decreasing = T)
res.mediation[metabolite %in% c("C2"), table(gene)] %>% sort(decreasing = T) %>% head(10)
res.enrich.tables[p.adjust<=0.05 & Cluster %in% c("C2"), table(Description)] %>% sort(decreasing = T) %>% head(10)

# in combination with CrAT

res.meta[symbol_INGENUITY=="CRAT" & hfdr.sigREM==TRUE]
res.sobel[ symbol_INGENUITY == "CRAT", ]

res.mediation[metabolite=="C2"&gene=="CRAT"] # nüscht
res.mediation[gene=="CRAT"] # nüscht

# C4OH
res.mediation[metabolite=="C4OH",table(mediation.group)]
res.mediation[metabolite=="C4OH",table(gene)] %>% sort(decreasing = T) %>% head(10)
res.meta[symbol_INGENUITY=="PDK4" & hfdr.sigREM]
res.mediation[metabolite=="C4OH"&gene=="PDK4"]

#' ## AA derived ACs
#  AA derived ACs ----

res.sobel[type=="pheno~metab" & hfdr.sigREM==TRUE, .N,by=pheno]

venn2(
  res.sobel[type=="pheno~metab" & hfdr.sigREM==TRUE & pheno=="diabetes.status.tri", unique(metabolite)],
  res.sobel[type=="pheno~metab" & hfdr.sigREM==TRUE & pheno=="log.bmi", unique(metabolite)],
  mylabels = c("T2DM","BMI")
)

#' ## Long Chain ACs
#  Long Chain ACs ----
s
res.mediation[metabolite %in% c("C16","C18","C181","C182"), table(gene)] %>% sort(decreasing = T) %>% head(10)

# 

#' ## Some Plots
# Some Plots ----

# Paper: Overlaps between groups ----
# verlap between significant GO pathways in the 4 groups
v1 <- venn3(
  res.enrich.tables[mediation.group == unique(mediation.group)[1] & data.base=="GO" & p.adjust<=0.05, Description],
  res.enrich.tables[mediation.group == unique(mediation.group)[2] & data.base=="GO" & p.adjust<=0.05, Description],
  res.enrich.tables[mediation.group == unique(mediation.group)[4] & data.base=="GO" & p.adjust<=0.05, Description],
  mylabels = c("Bi-Directional","m->gx->bmi","gx->m->bmi") #"both",
)

# check number of enriched metabolites
my.groups <- res.enrich.tables$mediation.group %>% unique
my.groups[1]
res.enrich.tables[mediation.group == my.groups[1] & p.adjust <= 0.05, table(Description)] %>% sort
res.enrich.tables[mediation.group == my.groups[1] & p.adjust <= 0.05, table(Cluster)] %>% sort

my.groups[2]
res.enrich.tables[mediation.group == my.groups[2] & p.adjust <= 0.05, table(Description)] %>% sort
res.enrich.tables[mediation.group == my.groups[2] & p.adjust <= 0.05, table(Cluster)] %>% sort

my.groups[3]
res.enrich.tables[mediation.group == my.groups[3] & p.adjust <= 0.05, table(Description)] %>% sort
res.enrich.tables[mediation.group == my.groups[3] & p.adjust <= 0.05, table(Cluster)] %>% sort

my.groups[4]
res.enrich.tables[mediation.group == my.groups[4] & p.adjust <= 0.05, table(Description)] %>% sort
res.enrich.tables[mediation.group == my.groups[4] & p.adjust <= 0.05, table(Cluster)] %>% sort

# same thing with genes
venn3(
  res.enrich.tables[mediation.group == unique(mediation.group)[1] & 
                      data.base=="GO" & 
                      p.adjust<=0.05, 
                    unlist(strsplit(enriched.genes,split = ","))],
  res.enrich.tables[mediation.group == unique(mediation.group)[2] & 
                      data.base=="GO" & 
                      p.adjust<=0.05, 
                    unlist(strsplit(enriched.genes,split = ","))],
  res.enrich.tables[mediation.group == unique(mediation.group)[4] & 
                      data.base=="GO" & 
                      p.adjust<=0.05, 
                    unlist(strsplit(enriched.genes,split = ","))],
  mylabels = c("Bi-Directional","m->gx->bmi","gx->m->bmi") #"both",
)

# same thing with metabolites
res.enrich.tables[, `:=`(class1 = res.meta[match(Cluster, metabolite),metab.class],
                         class2 = res.meta[match(Cluster, metabolite),metab.class.q])]
venn3(
  res.enrich.tables[mediation.group == unique(mediation.group)[1] & 
                      data.base=="GO" & 
                      p.adjust<=0.05, 
                    Cluster],
  res.enrich.tables[mediation.group == unique(mediation.group)[2] & 
                      data.base=="GO" & 
                      p.adjust<=0.05, 
                    Cluster],
  res.enrich.tables[mediation.group == unique(mediation.group)[4] & 
                      data.base=="GO" & 
                      p.adjust<=0.05, 
                    Cluster],
  mylabels = c("Bi-Directional","m->gx->bmi","gx->m->bmi")
  )


#' ## Interesting genes
# Interesting genes ----

# PTPRF -> involved in insulin resistance
res.enrich.tables[data.base=="GO" & p.adjust<=0.05, ][grep("PTPRF",enriched.genes), table(Description)]
res.enrich.tables[data.base=="GO" & p.adjust<=0.05, ][grep("PTPRF",enriched.genes), table(class2)]

# with metabolite classes
res.enrich.tables[p.adjust<=0.05 & data.base == "GO", .N,by=.(mediation.group, class1)][order(mediation.group)] %>% as.data.frame
res.enrich.tables[p.adjust<=0.05 & data.base == "GO", uniqueN(Cluster),by=.(mediation.group, class1)][order(mediation.group)] %>% as.data.frame
res.enrich.tables[p.adjust<=0.05 & data.base == "GO", uniqueN(Cluster),by=.(mediation.group, class2)][order(mediation.group)] %>% as.data.frame

# Top Metabolite: Q16, Q12, Q25
res.enrich.tables[p.adjust<=0.05 & data.base == "GO", .N,by=Cluster][order(N,decreasing = T)]

# Top enriched genes
res.enrich.tables[p.adjust<=0.05 & data.base == "GO", unlist(strsplit(enriched.genes,split = ","))] %>% table %>% sort


# bi-directional & GO Terms
dotplot(res.enrich.objects[[3]][[2]], 
        showCategory=2,
        font.size=10) + 
  theme(
    # axis.text.x = element_text(angle=45,hjust=1),
    # axis.text.y = element_text(angle=45,hjust=1,size = 8)
  ) +
  theme_gxMetab()

#=================================================================================#
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#=================================================================================#

# Metabolite enrichment ----
#' ### Metabolite enrichment

source("../../../07_programme/rtools/2003_metabolitePathwayEnrichment/functions/enrich_metabolites.R")
source("../../../07_programme/rtools/2003_metabolitePathwayEnrichment/functions/hyper_test.R")
source("../../../07_programme/rtools/2003_metabolitePathwayEnrichment/functions/prepare_enrich_metabolites.R")
dat <- prepare_enrich_metabolites()

bg.metab <- res.mediation[pheno=="log.bmi", unique(metabolite)]
fg.metab.1 <- res.mediation[pheno=="log.bmi" & mediation.group == "Metabolite mediated gene expression effects", unique(metabolite),by=gene]
res.mediation[pheno=="log.bmi" & mediation.group == "Metabolite mediated gene expression effects", uniqueN(metabolite),by=gene][order(V1)]

fg.metab.2 <- res.mediation[pheno=="log.bmi" & mediation.group == "Gene expression mediated metabolite effects", unique(metabolite)]
fg.metab.3 <- res.mediation[pheno=="log.bmi" & mediation.group == "Bi-directional mediation", unique(metabolite),by=gene]
res.mediation[pheno=="log.bmi" & mediation.group == "Bi-directional mediation", uniqueN(metabolite),by=gene][order(V1)]

enrich.metab.1 <- enrich_metabolites(
  fg = fg.metab.1[gene =="SLC4A7", V1],
  bg = bg.metab,
  convertQuotients = TRUE,
  verbose=TRUE,
  infoFile = dat)
enrich.metab.1[order(pval.hyper)]

enrich.metab.3 <- enrich_metabolites(
  fg = fg.metab.3[gene =="PTPRF", V1],
  bg = bg.metab,
  convertQuotients = TRUE,
  verbose=TRUE,
  infoFile = dat)
enrich.metab.3[order(pval.hyper)]


# -----
# -----
# -----

#=================================================================================#
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#=================================================================================#

# Enrichment of Sobel results for T2D & BMI ----
#' ## Enrichment of Sobel results for T2D & BMI

res.sobel[is.mediation=="yes",.N,by=.(type,pheno)]

# select both backgrounds

bg.sobel <- res.sobel[!is.na(entrez.id), unique(as.character(entrez.id)),by=y]
bg.sobel <- sapply(unique(bg.sobel$y),function(x){
  res <- bg.sobel[y==(x), unique(V1)]
  res <- unlist(strsplit(x = res, split = ", ", fixed = T))
  res <- list(res)
  names(res) <- x
  return(res)
},USE.NAMES = FALSE)

# foreground
fg.sobel <- res.sobel[!is.na(entrez.id) & fdr.sig==TRUE, unique(as.character(entrez.id)), by=.(y,metabolite)]

# select foregrounds
fg.res <-sapply(unique(fg.sobel$y), function(x){
  
  fg.y <- sapply(fg.sobel[y==x, unique(metabolite)], function(xx){
    
    res <- fg.sobel[y== (y) & metabolite==(xx), V1]
    res <- unlist(strsplit(x = res, split = ", ", fixed = T))
    res <- list(res)
    names(res) <- xx
    return(res)
    
  },USE.NAMES=FALSE)
  
  res <- list(fg.y)
  names(res) <- x
  return(res)
  
},USE.NAMES=FALSE)
fg.sobel <- fg.res

lapply(unique(res.sobel$y), function(pheno){
  
  fg <- fg.sobel[[pheno]]
  bg <- bg.sobel[[pheno]]
  
  # KEGG ----
  
  enrich.kegg <- compareCluster(
    geneCluster = fg, 
    fun = "enrichKEGG", 
    pvalueCutoff = 0.3, 
    universe = as.character(bg), 
    minGSSize = 2,
    use_internal_data=F,
    qvalueCutoff = 0.3
  )
  dotplot(enrich.kegg)
  
  # check KEGG results
  res.kegg <- enrich.kegg@compareClusterResult
  setDT(res.kegg)
  str(res.kegg)
  res.kegg[order(p.adjust,decreasing = FALSE), ] %>% head(10)
  
  # KEGG GSEA ----
  # geneList <- res.sobel[fdr.sig==TRUE&!is.na(entrez.id) &y=="diabetes.status.tri"&metabolite=="C2", 2^beta.tau.tauAdj]
  # names(geneList) <- res.sobel[fdr.sig==TRUE&!is.na(entrez.id) &y=="diabetes.status.tri"&metabolite=="C2", entrez.id]
  # geneList <- geneList  %>% sort(decreasing = T)
  # kk2 <- gseKEGG(geneList     = geneList,
  #                organism     = 'hsa',
  #                nPerm        = 1000,
  #                minGSSize    = 120,
  #                pvalueCutoff = 0.05,
  #                verbose      = FALSE)
  
  
  # GO ----
  
  enrich.go <- compareCluster(
    geneCluster = fg, 
    fun = "enrichGO", 
    OrgDb = org.Hs.eg.db,
    universe = as.character(bg), 
    minGSSize = 2,
    pvalueCutoff = 0.3, 
    qvalueCutoff = 0.3
  )
  dotplot(enrich.go,showCategory=2)
  res.go <- enrich.go@compareClusterResult
  setDT(res.go)
  res.go[order(p.adjust,decreasing = FALSE), ] %>% head(10)
  
  
  # check GO results
  res.go <- enrich.go@compareClusterResult
  setDT(res.go)
  str(res.go)
  res.go[order(p.adjust,decreasing = FALSE), ] %>% head(5)
  
  # some GO group comparison
  group.go <- compareCluster(
    geneCluster = fg, 
    fun = "groupGO", 
    OrgDb = org.Hs.eg.db
  )
  
  res.group <- group.go@compareClusterResult
  setDT(res.group)
  
  # geneList <- res.sobel[!is.na(entrez.id) & y == (pheno) & metabolite == "OHProl" & x.type == "gx.probe", 2^beta.tau.tauAdj]
  # names(geneList) <- res.sobel[!is.na(entrez.id) & y == (pheno) & metabolite == "OHProl" & x.type == "gx.probe", entrez.id] 
  # geneList <- sort(geneList,decreasing = TRUE)
  # res.go.gsea <- gseGO(geneList = geneList,
  #               OrgDb        = org.Hs.eg.db,
  #               ont          = "CC",
  #               nPerm        = 1000,
  #               minGSSize    = 100,
  #               maxGSSize    = 500,
  #               pvalueCutoff = 0.05,
  #               verbose      = FALSE)
  
  
})




















#=================================================================================#
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#=================================================================================#

# Enrichment of Meta-analysis results ----
#' ## Enrichment of Meta-analysis results

bg.meta <- res.meta[!is.na(entrez.id), unique(entrez.id)]
bg.meta <- unlist(strsplit(x = bg.meta, split = ", ", fixed = T))

fg.meta <- res.meta[!is.na(entrez.id) & significant==TRUE, entrez.id, by=metabolite]
fg.meta <-sapply(unique(fg.meta$metabolite), function(x){
  res <- fg.meta[metabolite==(x), unique(entrez.id)]
  res <- unlist(strsplit(x = res, split = ", ", fixed = T))
  res <- list(res)
  names(res) <- x
  return(res)
},USE.NAMES=FALSE)

# KEGG ----
#' ### KEGG ----

#' Feingliedrig, funktion, verständlichkeit

# enrichment
# c2.kegg <- enrichKEGG(gene = fg.meta[["C2"]],
#                       universe = as.character(bg.meta),
#                       pvalueCutoff = .3,
#                       minGSSize = 2,
#                       qvalueCutoff = .3,
#                       use_internal_data = FALSE)
# c2.res <- c2.kegg@result
# setDT(c2.res)
# c2.res[order(p.adjust)] %>% head(10)

enrich.kegg <- compareCluster(
  geneCluster = fg.meta, 
  fun = "enrichKEGG", 
  pvalueCutoff = 0.3, 
  universe = as.character(bg.meta), 
  minGSSize = 2,
  use_internal_data=F,
  qvalueCutoff = 0.3
)

# dotplot
hh(enrich.kegg@compareClusterResult,8)
hh(summary(enrich.kegg),7)
dotplot(enrich.kegg, showCategory=1)
res.kegg <- enrich.kegg@compareClusterResult
setDT(res.kegg)
res.kegg[order(p.adjust)]

# check for
metab <- "Sarc"
edo <- enrichDGN(gene         = as.character(fg.meta[[metab]]), 
                 universe     = as.character(bg.meta),
                 qvalueCutoff = 0.3)

# another plot
emapplot(edo)

# beta (==logFC) to foldChange = 2^logFC
AspList <- res.meta[metabolite==metab & entrez.id %in% fg.meta[[metab]], 2^betaREM]
names(AspList) <- res.meta[metabolite==metab & entrez.id %in% fg.meta[[metab]], entrez.id]

# gene concept network
# convert gene ID to Symbol
edox <- setReadable(edo, 'org.Hs.eg.db', 'ENTREZID')
p1 <- cnetplot(edox, foldChange=AspList)
p2 <- cnetplot(edox, categorySize="pvalue", foldChange=geneList)
p3 <- cnetplot(edox, foldChange=geneList, circular = TRUE, colorEdge = TRUE)

# view kegg pathway
# browseKEGG(kk, 'hsa04110')


# GO ----
#' ### GO

#' Umfangreich, allgemeiner

library("org.Hs.eg.db")
enrich.go <- compareCluster(
  geneCluster = fg.meta, 
  fun = "enrichGO", 
  OrgDb = org.Hs.eg.db,
  universe = as.character(bg.meta), 
  minGSSize = 2,
  pvalueCutoff = 0.2, 
  qvalueCutoff = 0.2
)
dotplot(enrich.go, showCategory=2)

# GO Gene Set Enrichment Analysis
# meta.gsea <- GSEA(geneList = fg.meta,nPerm = 1000)
# ego3 <- gseGO(geneList     = geneList,
#               OrgDb        = org.Hs.eg.db,
#               ont          = "CC",
#               nPerm        = 1000,
#               minGSSize    = 100,
#               maxGSSize    = 500,
#               pvalueCutoff = 0.05,
#               verbose      = FALSE)

# DO ----
#' ### DO

#' Krankheitsbezogen

enrich.do <- compareCluster(
  geneCluster = fg.meta, 
  fun = "enrichDO", 
  universe = as.character(bg.meta)
)
dotplot(enrich.do,includeAll=TRUE)

# Reactome ----
#' ### Reactome

library("ReactomePA")
enrich.path <- compareCluster(
  geneCluster = fg.meta, 
  fun = "enrichPathway", 
  universe = as.character(bg.meta)
)
dotplot(enrich.do,includeAll=FALSE)

# library("pathview")
# hsa04110 <- pathview(gene.data  = geneList,
#                      pathway.id = "hsa04110",
#                      species    = "hsa",
#                      limit      = list(gene=max(abs(geneList)), cpd=1))



save(list = c("enrich.kegg","enrich.go","enrich.do", "enrich.path"), file = "res/GeneEnrichment.RData")
