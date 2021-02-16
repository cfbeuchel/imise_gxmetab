#' ---
#' title: "Visualization of metabolites and gx-probe PCs"
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
  if(grepl(x = getwd(), pattern =  "mnt")){
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
  "limma",
  "ggplot2",
  "magrittr",
  "uwot"
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
source(here("../functions/beta_plot.R"))
source(here("../functions/all_annotation_tables.R"))
an <- all_annotation_tables(mount.point = "/net/ifs2/san_projekte/projekte/genstat/")
option_setup()

#+ include=F
#==========#
#   DONE   #
#==========#

#' # Load metabolite data

#' The data was cleaned in STEP1
load("res/dat_for_STEP2.RData", verbose = TRUE)

# all metabolites
am <- an$metab.annot$metabolite
an$metab.annot[,pheno.class.q:=pheno.class]
an$metab.annot[grepl("^Q",metabolite),pheno.class.q:="quotient"]

# load supplementary table 1
st1 <- fread("obj/191209_metabOverview.csv")

#' # Metabolites

# amd = all metab data
am <- c("id", am)
amd <- rbindlist(
  list(
    metab.nc.a1[,..am],
    metab.nc.b3.heart[,..am],
    metab.nc.b3.ami[,..am],
    metab.nc.sorb[,..am]
  ) 
) %>% na.omit

# add cohort identifier
covar.nc.a1[,cohort:="a1"]
covar.nc.b3.heart[,cohort:="heart"]
covar.nc.b3.ami[,cohort:="ami"]
covar.nc.sorb[,cohort:="sorb"]

# acd = all covar data
acd <- rbindlist(
  list(
    covar.nc.a1,
    covar.nc.b3.heart,
    covar.nc.b3.ami,
    covar.nc.sorb
  ),fill = TRUE 
)

# filter down to metabolite data
acd <- acd[id %in% amd$id, ]

# reorder to metabolite IDs
m1 <- match(amd$id, acd$id)
acd <- acd[m1, ]

# plot PCs and or UMaps
# pcm = PCs of the metabolites
pcm <- amd[,.SD, .SDcols = names(amd)[-1]] %>% prcomp()
summary(pcm)


tiff(dpx("metabPC_300dpi.tiff", "plt/"), 
     width = 7, 
     height = 7, 
     res=300, 
     unit="in")
plot(pcm$x[,1]~pcm$x[,2], pch=20, col=alpha(
  factor(
    acd$cohort, 
    labels = c(
      "a1" = "blue",
      "heart" = "red", 
      "ami" = "green", 
      "sorb" = "yellow"
    )
  ),
  .7))
dev.off()

#' ## Plot via UMAP
metab.umap <- umap(
  # X = pcm$x[,1:50],
  X = amd[,.SD, .SDcols = names(amd)[-1]],
  min_dist = .001, #  default: min_dist = 0.1
  spread = 1.5 # default: spread = 1
  # a = 0.1, # default: a = 0.0001
  # b = 0.9 # default: b = 100
) # tweak the parameters to get proper separation

tiff(dpx("metabUMAP_300dpi.tiff", "plt/"), 
     width = 7, 
     height = 7, 
     res=300, 
     unit="in")
plot(metab.umap[,1]~metab.umap[,2],
     col=alpha(
       factor(
         acd$cohort, 
         labels = c(
           "a1" = "blue",
           "heart" = "red", 
           "ami" = "green", 
           "sorb" = "yellow"
         )
       ),
       .7),
     pch=20)
dev.off()

#' # Gene Expression Probes

# get the overlapping probe IDs

# melt for correct PC calc
# sorben
gx.sorb <- as.data.frame(exprs(eset.nc.sorb.f))
setDT(gx.sorb,keep.rownames = "probe")
gx.sorb <- dcast.data.table(melt(gx.sorb,
                                 id.vars="probe",
                                 variable.name="id", 
                                 variable.factor = FALSE),
                            formula=id~probe)

# adult
gx.a1 <- as.data.frame(exprs(eset.nc.a1.f))
setDT(gx.a1,keep.rownames = "probe")
gx.a1 <- dcast.data.table(melt(gx.a1,
                               id.vars="probe",
                               variable.name="id",
                               variable.factor = FALSE),
                          formula=id~probe)

# ami
gx.b3.ami <- as.data.frame(exprs(eset.nc.b3.ami.f))
setDT(gx.b3.ami,keep.rownames = "probe")
gx.b3.ami <- dcast.data.table(melt(gx.b3.ami,
                                   id.vars="probe",
                                   variable.name="id", 
                                   variable.factor = FALSE),
                              formula=id~probe)

# heart
gx.b3.heart <- as.data.frame(exprs(eset.nc.b3.heart.f))
setDT(gx.b3.heart,keep.rownames = "probe")
gx.b3.heart <- dcast.data.table(melt(gx.b3.heart,
                                     id.vars="probe",
                                     variable.name="id", 
                                     variable.factor = FALSE),
                                formula=id~probe)

v1 <- venn4(
  featureNames(eset.nc.a1.f),
  featureNames(eset.nc.b3.ami.f),
  featureNames(eset.nc.b3.heart.f),
  featureNames(eset.nc.sorb.f)
)

# check 2
v3 <- venn4(
  colnames(gx.a1)[-1],
  colnames(gx.b3.ami)[-1],
  colnames(gx.b3.heart)[-1],
  colnames(gx.sorb)[-1]
)

# check
v2 <- venn4(
  as.character(gx.sorb$id),
  as.character(gx.a1$id),
  as.character(gx.b3.ami$id),
  as.character(gx.b3.heart$id)
)

# get the overlap of IDs 
# agx = all gx ids
agx <- v1$q1

# all gx sample id
gsi <- c(as.character(gx.sorb$id),
         as.character(gx.a1$id),
         as.character(gx.b3.ami$id),
         as.character(gx.b3.heart$id))

# jgx = joint gx data
jgx <- rbindlist(
  list(
    gx.sorb[,..agx],
    gx.a1[,..agx],
    gx.b3.ami[,..agx],
    gx.b3.heart[,..agx]
  ))

# nag = NAs in gx data
nag <- lapply(jgx,function(x){
  which(is.na(x))
})
str(nag)

# pci = PC ids 
pci <- which(gsi %nin% nag[[1]])
jgx[nag[[1]], .SD,.SDcols = agx[1:10]]

# Whats up with these samples????
eset.nc.sorb.f[, gsi[nag[[1]]]]

# They will probably be removed?

# calculate PCs of joint gx data
pcg <- rbindlist(
  list(
    gx.sorb[,..agx],
    gx.a1[,..agx],
    gx.b3.ami[,..agx],
    gx.b3.heart[,..agx]
  )
) %>% na.omit %>% prcomp(center=TRUE)

# save for later use
fwrite(pcg$x, dpx("gxUnadjusted_all_PCs.csv", "obj/"))
# pcg <- fread("obj/191211_gxUnadjusted_all_PCs.csv")
sadasdasd
# matrix for handling
pcg <- pcg %>% as.matrix

#' ## Plot PCs

plot(pcg[,1]~pcg[,2], 
     pch=20, 
     col= alpha(
       factor(
         acd[pci, cohort], 
         labels = c("a1" = "blue",
                    "heart" = "red", 
                    "ami" = "green", 
                    "sorb" = "yellow"
         )
       ),
       .7)
)

# calculate UMAP 
gx.umap <- umap(
  X = pcg$x[,1:50],
  min_dist = .001, #  default: min_dist = 0.1
  spread = 1.5) # default: spread = 1) # tweak the parameters to get proper separation
plot(gx.umap[,1]~gx.umap[,2],
     col=alpha(
       factor(
         acd[pci, cohort], 
         labels = c("a1" = "blue",
                    "heart" = "red", 
                    "ami" = "green", 
                    "sorb" = "yellow"
         )
       ),
       .7),
     pch=20)