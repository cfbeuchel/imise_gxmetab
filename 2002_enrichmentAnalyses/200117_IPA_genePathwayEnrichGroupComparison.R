#' ---
#' title: "IPA Group Comparison Summaries"
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
  "ComplexHeatmap"
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
option_setup()

#+ include=F
#==========#
#   DONE   #
#==========#

# Load data ----

# afs = all files
afs <- list.files(path = "obj/IPA/200120_ipa_comparisonResults/")

# no csv
afs <- afs[!grepl(pattern = "csv$", afs)]

# read reference data
nf <- newest_file(look_for = "metabOverview",subfolder = "obj",print_full = TRUE)
st1 <- fread(nf)

# Read single files ----

# spi = super path identifier
# spi <- st1[,.(rid, metab.super.path)]
spi <- data.table(
  long = st1[,unique(metab.super.path)]
  )
spi[,short:=c(
  "CarnitineTransport",
  "FattyAcidMetabolism",
  "EnergyMetabolism",
  "BCAAMetabolism",
  "AminoAcidMetabolismOther",
  "Other",
  "ureaCycle",
  "AmmoniaRecycling",
  "CollagenSynthesis",
  "FolicAcidCycle"
  )]

# IPA table types
# tbi = table identifier
tbi <- lapply(strsplit(afs, split = "_",fixed = TRUE), function(x){
  x[2]
})
tbi <- unique(unlist(tbi))
# order the files for each group into a list

# Annotate files ----

# get all the analyses of each type
all.files <- sapply(tbi, function(x){
  
  files <- grep(pattern = paste0("_", x, "_"), afs, value = T, fixed = T)
  
  # extract and clean file name
  tbl.type <- lapply(strsplit(files, split = "_", fixed = TRUE), `[[`, 3)
  tbl.type <- gsub(tbl.type, pattern = ".txt", replacement = "", fixed = TRUE)
  tbl.type <- gsub(tbl.type, pattern = "[1-2]", replacement = "")
  # tbl.type <- gsub(tbl.type, pattern = "\\d", replacement = "")
  
  # assign correct name of superpathway
  m1 <- match(tbl.type, spi$short)
  names(files) <- spi[m1, long]
  
  # try to read the files
  files2 <- lapply(files, function(y){
    tbl <- try(read.delim(paste0("obj/IPA/200120_ipa_comparisonResults/", y), sep = "\t", skip = 2))
    return(tbl)
  })
  
  # return this including the name of the table type
  res <- list(files2)
  
  return(res)

}, USE.NAMES = TRUE)

# check
names(all.files)
lapply(all.files, names)

# create tables for each pathway-type
all.files$NetworksTable %>% str

# Canonical Pathways ----

# get new file for the first element
all.files[[names(all.files)[1]]]$`Amino acid metabolism, other` %>% str
can.path <- all.files[[names(all.files)[1]]]

# create molten long table of canonical enrichment
can.path.long <- lapply(seq_along(can.path), function(x){
  
  # format and melt
  dat <- can.path[[x]]
  setDT(dat)
  res <- melt.data.table(data = dat,
                  id.vars = "Canonical.Pathways", 
                  variable.name = "metabolite", 
                  value.name = "enrichment", 
                  variable.factor = FALSE)
  res$enrichment <- as.numeric(res$enrichment)
  res[, super.path := names(can.path)[[x]]]
  
})

# acp = all canonical pathway info
acp <- rbindlist(can.path.long)

# remove non-enriched Pathways/metabolites
na.paths <- acp[, all(is.na(enrichment)), by = Canonical.Pathways][V1==TRUE,Canonical.Pathways]
na.metab <- acp[, all(is.na(enrichment)), by = metabolite][V1==TRUE,metabolite]

acp <- acp[Canonical.Pathways %nin% na.paths & metabolite %nin% na.metab, ]

# dcast to display as heatmap
acp2 <- dcast(acp, formula = super.path + metabolite ~ Canonical.Pathways, value.var = "enrichment")
acp3 <- copy(acp2)
acp.super.paths <- acp3$super.path
acp3$super.path <- NULL
setDF(acp3, rownames = acp3$metabolite)
acp3$metabolite <- NULL
acp3 <- as.matrix(acp3)

# looks uninteresting
Heatmap(acp3, 
        cluster_rows = FALSE, 
        cluster_columns = FALSE, 
        row_names_gp = gpar(fontsize = 6),
        column_names_gp = gpar(fontsize = 6)
        )

# Regulator Effects ----

# extract info
which.table <- names(all.files)[6]
reg.eff.tbl <- all.files[[which.table]]

# quick look into the data
names(reg.eff.tbl[[1]])
str(reg.eff.tbl)

# annotate single sets with super.pathway
reg.eff.tbl2 <- lapply(seq_along(reg.eff.tbl), function(x){
  
  # format and melt
  dat <- reg.eff.tbl[[x]]
  setDT(dat)
  dat[, Consistency.Score := as.numeric(sub(pattern = ",",
                                            replacement = ".",
                                            x = Consistency.Score,
                                            fixed = TRUE))]
  dat[, super.path := names(reg.eff.tbl)[[x]]]
  return(dat)
  
})

# combine
reg.eff.tbl3 <- rbindlist(reg.eff.tbl2)

# make the column with the overall connection info machine readable
known.relationships.perc <- lapply(strsplit(reg.eff.tbl3$Known.Regulator.Disease.Function.Relationship,split = " ", fixed = T),`[[`, 1)
known.relationships.perc <- as.numeric(sub(unlist(known.relationships.perc), pattern = "%",replacement = "",fixed = TRUE))/100
known.relationships.num <- lapply(strsplit(reg.eff.tbl3$Known.Regulator.Disease.Function.Relationship,split = " ", fixed = T),`[[`, 2)
known.relationships.num <- unlist(known.relationships.num)
known.relationships.num <- as.integer(regmatches(known.relationships.num,regexpr("[0-9]+",known.relationships.num)))

# enter variables
reg.eff.tbl3[,`:=`(
  known.regulator.rel.perc = known.relationships.perc,
  known.regulator.rel.num = known.relationships.num
)]

# preview
reg.eff.tbl3[order(Consistency.Score,decreasing = TRUE),]
reg.eff.tbl3[order(known.regulator.rel.num,decreasing = TRUE),] %>% head()

# Upstream Regulators ----

which.table <- names(all.files)[10]
ups.reg.tbl <- all.files[[which.table]]

# quick look into the data
names(ups.reg.tbl[[1]])
str(ups.reg.tbl[[1]])

# annotate single sets with super.pathway
ups.reg.tbl2 <- lapply(seq_along(ups.reg.tbl), function(x){
  
  # format
  dat <- ups.reg.tbl[[x]]
  setDT(dat)
  dat[, super.path := names(ups.reg.tbl)[[x]]]
  return(dat)
})

# flags col is not present in every instance
lapply(ups.reg.tbl2,names)

# bind together, fill the flags col with NAs if necessary
ups.reg.tbl3 <- rbindlist(ups.reg.tbl2,fill = TRUE)

# several comma-separated cols -> turn numeric
num.cols <- c("Expr.False.Discovery.Rate..q.value.",
              "Expr.Log.Ratio", 
              "Expr.p.value",
              "Activation.z.score",
              "p.value.of.overlap")

# turn numeric
ups.reg.tbl3[, (num.cols) := lapply(.SD, function(x){
  as.numeric(sub(pattern = ",",
                 replacement = ".",
                 x = x,
                 fixed = TRUE))
  }),.SDcols=(num.cols)]

setorder(ups.reg.tbl3, Expr.False.Discovery.Rate..q.value.,-Activation.z.score, -Expr.Log.Ratio, na.last = TRUE)

# number of regulators
ups.reg.tbl3$Upstream.Regulator %>% uniqueN






# Tox Lists ----

