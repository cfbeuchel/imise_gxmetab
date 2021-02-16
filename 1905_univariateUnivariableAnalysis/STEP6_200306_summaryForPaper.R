#' ---
#' title: "Summary statistics from primary analysis"
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
  nc <- parallel::detectCores()-6
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
  "ComplexHeatmap",
  "gridExtra",
  "corrplot",
  "magrittr",
  "colorspace",
  "igraph"
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
source("../functions/all_annotation_tables.R")
source("../functions/option_setup.R")
source("../functions/theme_gxMetab.R")
an <- all_annotation_tables(mount.point = bp)
option_setup()

#+ include=F
#==========#
#   DONE   #
#==========#

# Load data ----

#' I evaluated annotated data in the previous step.

theme_set(theme_light(base_size = 15, base_family = "Helvetica"))
setDTthreads(8)

nf <- newest_file(look_for = "STEP4_meta_resultsAnnotated",subfolder = "res", print_full = TRUE)
dat <- fread(nf)

#' These are the correct names of all metabolites and quotients
load("obj/200107_metaboliteIDs.RData")

# load supplementary table 1
nf <- newest_file(look_for = "metabOverview",subfolder = "obj",print_full = TRUE)
st1 <- fread(nf)

# metabolite-metablite beta correlation
mat <- read.csv(file = "obj/191203_fullBetaCorrMatrix.csv", row.names = 1, header = TRUE)
mat <- as.matrix(mat) # was ist mit den spalten namen los???
colnames(mat) <- rownames(mat)

# get the long beta-corr table including p-values
nf <- newest_file("meta_betaCorrelationR2", subfolder = "res",print_full = TRUE)
all.combinations <- fread(nf)

# enrichment data
nf <- newest_file(look_for = "masterRegulatorEnrichment",subfolder = "obj", print_full = TRUE)

# gwe = gene-wise enrichment
gwe <- fread(nf)

# rename metabolites agn etc...
load(paste0(bp,"02_projekte/1703_ge_metab_a1_b3_sorbs/1805_gxMetaboliteAssociation/obj/200306_metaboliteNames.RData"))

#++++++++++++++++++++++++++++++++++++++++++++#
#============================================#

# Supp Table Meta-Analysis ----

#' I want a supplementary Table of the Meta-Analysis Results

# create a copy I can change and format & remove cols from
supp_table <- copy(dat)

# remove empty entries for correct display
supp_table[symbol_INGENUITY == "", symbol_INGENUITY := NA]
supp_table[location_INGENUITY == "", location_INGENUITY := NA]
supp_table[types_INGENUITY == "", types_INGENUITY := NA]
supp_table[biomarkers_INGENUITY == "", biomarkers_INGENUITY := NA]
supp_table[drugs_INGENUITY == "", drugs_INGENUITY := NA]

# add correct metabolite names
supp_table[, Metabolite_Abbreviation := agn.all[match(supp_table$metabolite, agn.all$rid), abbr]]
supp_table[, Metabolite := agn.all[match(supp_table$metabolite, agn.all$rid), full.name]]

# give the class a more descriptive name
supp_table[metab.class == "aa", metab.class := "Amino Acid"]
supp_table[metab.class == "ac", metab.class := "Acylcarnitine"]
supp_table[metab.class == "mix", metab.class := "Mix Quotient"]
supp_table[metab.class.q == "aa", metab.class.q := "Amino Acid"]
supp_table[metab.class.q == "ac", metab.class.q := "Acylcarnitine"]
supp_table[metab.class.q == "quotient", metab.class.q := "Quotient"]
supp_table$metab.class %>% mytable
supp_table$metab.class.q %>% mytable

supp_table <- supp_table[,.(
  
  Probe_ID = markerID,
  HGNC_ID = symbol_INGENUITY,
  
  Gene_Type = types_INGENUITY,
  Gene_Expression_Location = location_INGENUITY,
  Gene_Biomarker_Application = biomarkers_INGENUITY,
  Gene_Drug_Application = drugs_INGENUITY,
  
  Metabolite,
  Metabolite_Abbreviation,
  
  Metabolite_Class_1 = metab.class,
  Metabolite_Class_2 = metab.class.q,
  Metabolite_Path_1 = metab.super.path,
  Metabolite_Path_2 = metab.path,
  
  Effect_Size = betaREM,
  `Lower_CI_95%` = ci.lREM,
  `Upper_CI_95%` = ci.uREM,
  Standard_Error = seREM,
  P_Value = pREM,
  Q_Value = qREM,
  hierarchical_FDR_Level_1 = hfdr.1REM,
  hierarchical_FDR_Level_2 = hfdr.2REM,
  `Significance_hierarchical_FDR=5%` = hfdr.sigREM,
  
  Studies_Analyzed = numberStudies,
  Total_N = totalN,
  Heterogeneity_I2 = I2,
  
  Analysis_In_Single_Study = analyzed.adult.heart.ami.sorb,
  Significant_Single_Study_Effect = significant.adult.heart.ami.sorb
)]

fwrite(supp_table,
       dpx(suffix = "metaAnalysisSumStats.tsv",
           folder = "ppr/supp_tables/"), sep = "\t", 
       quote = FALSE)

#============================================#
#++++++++++++++++++++++++++++++++++++++++++++#

# total number of assocs in single study
dat[!is.na(p.Adult) | !is.na(p.Heart) | !is.na(`p.Heart-AMI`) | !is.na(p.Sorb),
    uniqueN(markerID)]
dat[!is.na(p.Adult) | !is.na(p.Heart) | !is.na(`p.Heart-AMI`) | !is.na(p.Sorb),
    uniqueN(metabolite)]
dat[(!is.na(p.Adult) | !is.na(p.Heart) | !is.na(`p.Heart-AMI`) | !is.na(p.Sorb)) & 
      symbol_INGENUITY != "",
    uniqueN(symbol_INGENUITY)]


st1[,.N]
st1[,.N,metab.class]

# Table of single-study results 

# single study results table
strest <- data.table(
  "Study" = c("LIFE Adult", "LIFE Heart", "LIFE AMI","Sorb Study")
)

# total number of associations
totassoc = c(
  dat[!is.na(p.Adult), .(adult = .N)],
  dat[!is.na(p.Heart), .(heart = .N)],
  dat[!is.na(`p.Heart-AMI`), .(ami = .N)],
  dat[!is.na(p.Sorb), .(sorb = .N)]
)

# estimate eta1
library(fdrtool)
eta1 = c(
  dat[!is.na(p.Adult), .(adult = 1-fdrtool::pval.estimate.eta0(p.Adult,diagnostic.plot = F))],
  dat[!is.na(p.Heart), .(heart = 1-fdrtool::pval.estimate.eta0(p.Heart,diagnostic.plot = F))],
  dat[!is.na(`p.Heart-AMI`), .(ami = 1-fdrtool::pval.estimate.eta0(`p.Heart-AMI`,diagnostic.plot = F))],
  dat[!is.na(p.Sorb), .(sorb = 1-fdrtool::pval.estimate.eta0(p.Sorb,diagnostic.plot = F))]
)

# total sigs fdr
totalsigs = c(
  dat[!is.na(p.Adult), .(adult = sum(hfdr.sig.Adult))],
  dat[!is.na(p.Heart), .(heart = sum(hfdr.sig.Heart))],
  dat[!is.na(`p.Heart-AMI`), .(ami = sum(`hfdr.sig.Heart-AMI`))],
  dat[!is.na(p.Sorb), .(sorb = sum(hfdr.sig.Sorb))]
)

# associating metabolites
umetabolites = c(
  dat[!is.na(p.Adult)& hfdr.sig.Adult == TRUE, .(adult = uniqueN(metabolite))],
  dat[!is.na(p.Heart)& hfdr.sig.Heart == TRUE, .(heart = uniqueN(metabolite))],
  dat[!is.na(`p.Heart-AMI`)& `hfdr.sig.Heart-AMI` == TRUE, .(ami = uniqueN(metabolite))],
  dat[!is.na(p.Sorb)& hfdr.sig.Sorb == TRUE, .(sorb = uniqueN(metabolite))]
)

# associating probes
uprobes = c(
  dat[!is.na(p.Adult)& hfdr.sig.Adult == TRUE, .(adult = uniqueN(markerID))],
  dat[!is.na(p.Heart)& hfdr.sig.Heart == TRUE, .(heart = uniqueN(markerID))],
  dat[!is.na(`p.Heart-AMI`)& `hfdr.sig.Heart-AMI` == TRUE, .(ami = uniqueN(markerID))],
  dat[!is.na(p.Sorb)& hfdr.sig.Sorb == TRUE, .(sorb = uniqueN(markerID))]
)

# associating genes
ugenes = c(
  dat[!is.na(p.Adult) & symbol_INGENUITY != "" & hfdr.sig.Adult == TRUE, .(adult = uniqueN(symbol_INGENUITY))],
  dat[!is.na(p.Heart) & symbol_INGENUITY != "" & hfdr.sig.Heart == TRUE, .(heart = uniqueN(symbol_INGENUITY))],
  dat[!is.na(`p.Heart-AMI`) & symbol_INGENUITY != "" & `hfdr.sig.Heart-AMI` == TRUE, .(ami = uniqueN(symbol_INGENUITY))],
  dat[!is.na(p.Sorb) & symbol_INGENUITY != "" & hfdr.sig.Sorb == TRUE, .(sorb = uniqueN(symbol_INGENUITY))]
)

# exclusively associating probes
vtmp <- venn4(
  dat[!is.na(p.Adult) & hfdr.sig.Adult == TRUE, unique(markerID)],
  dat[!is.na(p.Heart) & hfdr.sig.Heart == TRUE, unique(markerID)],
  dat[!is.na(`p.Heart-AMI`) & `hfdr.sig.Heart-AMI` == TRUE, unique(markerID)],
  dat[!is.na(p.Sorb) & hfdr.sig.Sorb == TRUE, unique(markerID)],
  mylabels = c("adult","heart","ami","sorb")
)
eprobes <- c(length(vtmp$q10), length(vtmp$q11), length(vtmp$q12), length(vtmp$q13))

# exclusively associating genes
vtmp <- venn4(
  dat[!is.na(p.Adult) & symbol_INGENUITY != "" & hfdr.sig.Adult == TRUE, unique(symbol_INGENUITY)],
  dat[!is.na(p.Heart) & symbol_INGENUITY != "" & hfdr.sig.Heart == TRUE, unique(symbol_INGENUITY)],
  dat[!is.na(`p.Heart-AMI`) & symbol_INGENUITY != "" & `hfdr.sig.Heart-AMI` == TRUE, unique(symbol_INGENUITY)],
  dat[!is.na(p.Sorb) & symbol_INGENUITY != "" & hfdr.sig.Sorb == TRUE, unique(symbol_INGENUITY)],
  mylabels = c("adult","heart","ami","sorb")
)
egenes <- c(length(vtmp$q10), length(vtmp$q11), length(vtmp$q12), length(vtmp$q13))

# exclusively associating metabolites
vtmp <- venn4(
  dat[!is.na(p.Adult) &  hfdr.sig.Adult == TRUE, unique(metabolite)],
  dat[!is.na(p.Heart) &  hfdr.sig.Heart == TRUE, unique(metabolite)],
  dat[!is.na(`p.Heart-AMI`)  & `hfdr.sig.Heart-AMI` == TRUE, unique(metabolite)],
  dat[!is.na(p.Sorb) & hfdr.sig.Sorb == TRUE, unique(metabolite)],
  mylabels = c("adult","heart","ami","sorb")
)
emetabolites <- c(length(vtmp$q10), length(vtmp$q11), length(vtmp$q12), length(vtmp$q13))

# enter into table
strest[["Total Number of associations"]] <- unlist(totassoc)
strest[["Eta 1"]] <- round(unlist(eta1),2)
strest[["Total significant associations"]] <- unlist(totalsigs)
strest[["Associating probes"]] <- unlist(uprobes)
strest[["Associating genes"]] <- unlist(ugenes)
strest[["Associating metabolites"]] <- unlist(umetabolites)
strest[["Exclusively associating probes"]] <- (eprobes)
strest[["Exclusively associating genes"]] <- (egenes)
strest[["Exclusively associating metabolites"]] <- (emetabolites)
strest.md <- dcast(melt(strest, id.vars = "Study"), variable~Study, value.var = "value")
fwrite(strest.md,file = dpx("singleStudySummary.csv","res/"))

# template for more cols
# tmp = c(
#   dat[!is.na(p.Adult), .(adult = )],
#   dat[!is.na(p.Heart), .(heart = )],
#   dat[!is.na(`p.Heart-AMI`), .(ami = )],
#   dat[!is.na(p.Sorb), .(sorb = )]
# )

# overlapping gene-metabolite assocs
dat[,assoc.pairs := paste0(metabolite,"__",symbol_INGENUITY)]
vtmp <- venn4(
  dat[!is.na(p.Adult) & symbol_INGENUITY != "" & hfdr.sig.Adult == TRUE, unique(assoc.pairs)],
  dat[!is.na(p.Heart) & symbol_INGENUITY != "" & hfdr.sig.Heart == TRUE, unique(assoc.pairs)],
  dat[!is.na(`p.Heart-AMI`) & symbol_INGENUITY != "" & `hfdr.sig.Heart-AMI` == TRUE, unique(assoc.pairs)],
  dat[!is.na(p.Sorb) & symbol_INGENUITY != "" & hfdr.sig.Sorb == TRUE, unique(assoc.pairs)],
  mylabels = c("adult","heart","ami","sorb")
)

# total all 4 overlap
pairtoverlap <- length(vtmp$q1)

# check these
dat[assoc.pairs %in% vtmp$q1, .(genes = uniqueN(symbol_INGENUITY),
                                metabolite = uniqueN(metabolite))]
  
dat[assoc.pairs %in% vtmp$q1, xtabs( ~symbol_INGENUITY + metabolite)]
dat[assoc.pairs %in% vtmp$q1, unique(assoc.pairs)]

# minus sorb overlapping
pairminussorboverlap <- length(vtmp$q4)
dat[assoc.pairs %in% vtmp$q4, .(genes = uniqueN(symbol_INGENUITY),
                                metabolite = uniqueN(metabolite))]

# at least two studies overlappint
pairmin2overlap <- sapply(c(1:9, 14, 15),function(i){
  length(vtmp[[paste0("q",i)]])
}) %>% sum
pairmin2overlapnames <- sapply(c(1:9, 14, 15),function(i){
  vtmp[[paste0("q",i)]]
}) %>% unlist
dat[assoc.pairs %in% (pairmin2overlapnames),
    .(genes = uniqueN(symbol_INGENUITY),
      metabolite = uniqueN(metabolite))]

# significancies with sign
dat[,table(analyzed.adult.heart.ami.sorb)]
dat[,table(significant.adult.heart.ami.sorb)]

# effect signs
dat[assoc.pairs %in% vtmp$q1,table(significant.adult.heart.ami.sorb)]
dat[assoc.pairs %in% vtmp$q4 & (
  hfdr.sig.Sorb ==T | hfdr.sig.Adult == T | hfdr.sig.Heart ==T | `hfdr.sig.Heart-AMI` == T
  ),.N]
dat[assoc.pairs %in% vtmp$q4 & (
  hfdr.sig.Sorb ==T | hfdr.sig.Adult == T | hfdr.sig.Heart ==T | `hfdr.sig.Heart-AMI` == T
  ),table(significant.adult.heart.ami.sorb)]
9+6+59+12+10+55+1+4+14+2+1


#++++++++++++++++++++++++++++++++++++++++++++#
#============================================#

#============================================#
#++++++++++++++++++++++++++++++++++++++++++++#

# Beta-Correlation Plot ----

#' # Beta-Correlation Plot

# preprate annotation for color labels based on class
m1 <- match(rownames(mat), st1$abbr)

# super pathways matching the matrix rows/cols
spw <- st1[m1, metab.super.path]

# create a color mapping
cls <- c("#E16A86", "#CB7F2F", "#9F9400", "#50A315", "#00AC79", "#00AAB7", 
         "#009ADE", "#A87BE4", "#DA65C3", "#E16A86")

mapping <- data.table(
  group = unique(spw),
  col = cls
  )
m2 <- match(spw, mapping$group)

# TODO: Add Legend for name coloring ----

tiff(dpx("betaCorrelationHeatmap_300dpi.tiff", "plt/"), 
     width = 10, 
     height = 10, 
     res=300, 
     unit="in")

# simple heatmap of this
ComplexHeatmap::Heatmap(
  mat, 
  col = colorRamp2(c(-1,0,1), 
                   c("#023FA5", "#E2E2E2", "#8E063B")),
  row_names_gp = gpar(
    col = mapping[m2, col],
    cex = 0.5
  ),
  column_names_gp = gpar(
    col = mapping[m2, col],
    cex = 0.5
  ),
  cluster_rows = TRUE, 
  cluster_columns = TRUE, 
  clustering_method_columns = "ward.D2", 
  clustering_method_rows = "ward.D2",
  clustering_distance_rows = "euclidean", 
  clustering_distance_columns = "euclidean", 
  heatmap_legend_param = list(title="Pearson's r"))

dev.off()

#++++++++++++++++++++++++++++++++++++++++++++#
#============================================#

#============================================#
#++++++++++++++++++++++++++++++++++++++++++++#

# PLOT: Correlation Plot ----

# load corr. matrix
# cam
tmp <- load("obj/metabCorrelation.RData")
tmp
load("obj/metabCorrelation.RData")

# match correlation matrix names to right names
m1 <- match(rownames(cam), st1$rid)
rownames(cam) <- st1[m1, abbr]
colnames(cam) <- st1[m1, abbr]

# reorder
m2 <- match(rownames(mat), rownames(cam))
cam <- cam[m2,m2]

# calculate difference
# diff.mat <- abs(mat-cam)
diff.mat <- cam
m2 <- match(spw, mapping$group)

# tiff(dpx("betaCorrelationMetabCorrelationDiffHeatmap_300dpi.tiff", "plt/"), 
tiff(dpx("metabCorrelation_300dpi.tiff", "plt/"), 
     width = 10, 
     height = 10, 
     res=300, 
     unit="in")


# simple heatmap of this
ComplexHeatmap::Heatmap(
  cam, 
  col = colorRamp2(c(-1,0,1), 
                   c("#023FA5", "#E2E2E2", "#8E063B")),
  # col = colorRamp2(c(-1,1),
  #                  c("#E2E2E2", "#8E063B")),
  row_names_gp = gpar(
    col = mapping[m2, col],
    cex = 0.5
  ),
  column_names_gp = gpar(
    col = mapping[m2, col],
    cex = 0.5
  ),
  cluster_rows = TRUE, 
  cluster_columns = TRUE, 
  clustering_method_columns = "ward.D2", 
  clustering_method_rows = "ward.D2",
  clustering_distance_rows = "euclidean", 
  clustering_distance_columns = "euclidean", 
  heatmap_legend_param = list(title="Metabolite Correlation\nPearson's r"))

dev.off()

#++++++++++++++++++++++++++++++++++++++++++++#
#============================================#

#============================================#
#++++++++++++++++++++++++++++++++++++++++++++#

# MAYBE: Supp graphic of effect size distribution over the metaboliltes

#' # Density of effect sizes per metabolite

dat.beta.dens <- dat[significant==TRUE, betaREM, by = .(metabolite,
                                                        symbol_INGENUITY, 
                                                        metab.class, 
                                                        metab.class.q, 
                                                        metab.super.path)]

# reorder by number of associations and by maximum abs beta
ref <- dat.beta.dens[,.(sum = length(symbol_INGENUITY),
                 max = max(abs(betaREM))
                 ),by=metabolite]
setorder(ref, -sum, -max)

# metabolite order
mord <- ref$metabolite

# effect size distribution per metabolite
dat.beta.dens$metabolite.f <- factor(dat.beta.dens$metabolite,
                                     ordered = T,
                                     levels = rev(mord))
p.tmp <- ggplot(dat.beta.dens, 
       aes(
         group=metabolite.f,
         y = metabolite.f,
         x = betaREM, 
         color = ..x..>0)) +
  geom_point(pch="|", cex = 2, alpha=1) +
  scale_color_manual(values = c("#023FA5","#8E063B")) +
  # scale_color_gradient2(low = "#023FA5",
  #                       mid =  "#E2E2E2",
  #                       high =  "#8E063B") +
  guides(xlab = "Standardized effect size",
         color="none") +
  theme_gxMetab() +
  theme(
    axis.title.y = element_blank(),
    panel.grid.major.x   = element_line(
      color = "grey30",
      linetype = "dotted")
  )

tiff(dpx("metaBetaDistributionPerMetabolite_300dpi.tiff", "plt/"), 
     width = 6, 
     height = 10, 
     res=300, 
     unit="in")
p.tmp
dev.off()

# Metabolite Correlation ----

#' # Metabolite Correlation 

#' This is done in a seperate Script

# PC Visualization ----

#' Visualize overall similaritie between studies as PCs of complete phenotypic
#' data. This is done in a seperate script under additional analyses.

# Master Regulater Genes Plot ----

#' # Look for master regulators 

#' ## get the top genes

#' The top regugulated/regulating genes all map to several probes!

# total metabolites per gene include duplications!
mrg.tot <- dat[significant==TRUE, .(sig.metabs = .N),by=.(symbol_INGENUITY)]
mrg.tot <- mrg.tot[order(sig.metabs,decreasing = TRUE),]
mrg.tot <- mrg.tot[symbol_INGENUITY!="", ]

# mrg = master regulator genes
mrg.probe <- dat[ , .(probes = uniqueN(markerID)),by=.(symbol_INGENUITY)][order(probes,decreasing = TRUE)]
hh(mrg.probe,10)
mrg.probe[symbol_INGENUITY %in% mrg.tot$symbol_INGENUITY[1:5], ]

# the unique metabolite count is the important one
mrg.unique <- dat[significant==TRUE, .(sig.metabs = uniqueN(metabolite),
                                       sig.probes = uniqueN(markerID),
                                       max.beta = max(abs(betaREM))),
                  by=.(symbol_INGENUITY, 
                       location_INGENUITY, 
                       types_INGENUITY)]

# order by number of sig.metabs AND max abs. beta
setorder(mrg.unique, -sig.metabs, -max.beta)
# mrg.unique <- mrg.unique[order(sig.metabs,decreasing = TRUE),]

#' Interesting to see how many unmapped probes associate
hh(mrg.unique, 10)

#' ## Show distribution of number of associations across lists

# 33 metabolites associate with unmapped probes! Remove those for the further analysis
mrg.unique <- mrg.unique[symbol_INGENUITY!="", ]

# show the cumulative associations over number of genes ordered by strongest genes
cum.genes <- cumsum(mrg.unique$sig.metabs)/sum(mrg.unique$sig.metabs)
mrg.unique[,cumulative.metab.assoc := cum.genes]
mrg.unique[,perc.sig.genes := (1:length(mrg.unique$sig.metabs))/length(mrg.unique$sig.metabs)]

# PAPER: % of metabolites associating top 5 genes ----

mrg.unique[perc.sig.genes <= 0.01, .(.N, cumulative.metab.assoc)]
mrg.unique[1:50, .(perc.sig.genes, cumulative.metab.assoc)]

# PLOT: Power Law - gene centric ----

# summarize frequency of each sig.metabs
totn <- nrow(mrg.unique)
pwlaw <- mrg.unique[, .(freq = length(symbol_INGENUITY)/totn),by=sig.metabs]

# colorschemes for plotting
pal.metabs <- sequential_hcl(palette = "OrYel",n=100)
pal.genes <- sequential_hcl(palette = "Teal",n=100)
pal.pairs <- sequential_hcl(palette = "Greens2",n=100)

# power law plot focused on genes
pwl_p1 <- ggplot(pwlaw,
       aes(
         x = sig.metabs,
         y = freq
       )) +
  scale_y_log10() +
  scale_x_log10() +
  annotation_logticks() +
  geom_point(col=pal.genes[1]) +
  ylab("Frequency") +
  xlab("Significant metabolites per gene") +
  annotate("rect", 
           xmin = 17.5, 
           xmax = 28.5, 
           ymin = 1*10^-4, 
           ymax = 1*10^-2.5,
           alpha = .2,
           fill="grey45") +
  theme_gxMetab() +
  theme(text=element_text(color="grey30"),
        panel.grid.major.y = element_line(
          color = "grey30",
          linetype = "dotted"
          ),
        panel.grid.major.x = element_line(
          color = "grey30",
          linetype = "dotted"
          )
        )

#' How many genes associate with more than one metabolite?

# see question above
mrg.unique[sig.metabs>1, .N] / nrow(mrg.unique)

#' Find the percentage of associations that is defined by a number of genes
dat.beta.dens$gene.f <- factor(dat.beta.dens$symbol_INGENUITY,ordered = T,levels = rev(mrg.unique$symbol_INGENUITY))

# plot of most associating genes and their effect sizes
p3 <- ggplot(dat.beta.dens[symbol_INGENUITY %in% mrg.unique$symbol_INGENUITY[1:50]], 
       aes(
         y = gene.f,
         x = betaREM, 
         col = metab.class)
       ) +
  scale_color_manual(values = c("aa"="darkgoldenrod",
                                "ac"="cadetblue4",
                                "mix"="darkolivegreen4"),
                     labels=c("Amino Acid",
                              "Acylcarnitine",
                              "Mix Quotient")) +
  scale_x_continuous(limits = c(-.21,.21)) + 
  geom_point(lwd=3, 
             stroke=1.1, 
             alpha=0.8, 
             pch="|") +
  ylab(NULL) +
  xlab("Standardized effect size") +
  theme_gxMetab() +
  theme(legend.position = "top",
        text = element_text(colour = "grey30"),
        legend.box = "vertical",
        legend.title = element_blank(),
        legend.text = element_text(),
        legend.direction = "horizontal",
        legend.justification = "left",
        panel.grid.major.y = element_line(
          color = "grey30",
          linetype = "dotted"
        ),
        panel.grid.major.x = element_line(
          color = "grey30",
          linetype = "dotted"
        ))

# arrange plot
tiff(dpx("masterRegulaterPlot_300dpi.tiff", "plt/"), 
     width = 12, 
     height = 10, 
     res=300, 
     unit="in")

grid.arrange(p3 + labs(tag = "C"), 
             pwl_p1 + labs(tag = "A"), 
             ggplot() + labs(tag = "B") + theme(panel.border = element_blank()), 
             layout_matrix = matrix(c(2,1,3,1), 2, 2, byrow = TRUE))

dev.off()

# Network of master regulators ----

# pmr = plot master regulators

# proper metabolite names
load("/net/ifs2/san_projekte/projekte/genstat/02_projekte/1703_ge_metab_a1_b3_sorbs/1805_gxMetaboliteAssociation/obj/200306_metaboliteNames.RData")

# ndt = network data
plot.genes <- mrg.unique$symbol_INGENUITY[1:50]
r.squared.cutoff <- 0.005
# perc.genes.cutoff <- 0.01
# pmr <- mrg.unique[ perc.sig.genes <= (perc.genes.cutoff), symbol_INGENUITY]

ndt <- dat[significant==TRUE & symbol_INGENUITY %in% (plot.genes) & r.squared >= r.squared.cutoff, .(symbol_INGENUITY, metabolite, r.squared, sign = sign(betaREM))]
# ndt <- dat[significant==TRUE & symbol_INGENUITY %in% (pmr) & r.squared >= r.squared.cutoff, .(symbol_INGENUITY, metabolite, r.squared)]

# simplify dt
ndt <- ndt[, .(gene = symbol_INGENUITY, metabolite, rsquared = r.squared, sign)]

# match correct names
m1 <- match(ndt$metabolite, agn.all$rid)
ndt$metabolite <- agn.all[m1, abbr]


# build network

source("../functions/master_regulator_network.R")

tiff(dpx("masterRegulaterNetwork_300dpi.tiff", "plt/"), 
     width = 12, 
     height = 10, 
     res=300, 
     unit="in")

master_regulator_network(dat = ndt)

dev.off()


# Primary Tables ----

#' # Primary table

#' ## Table Superpathway-centric Overview

# Supplementary Tables ----

#' # Supplementary Table

# Supp Table 1 - Metabolite Overview ----

#' ## Supp Table 1 - Metabolite overview

#' Add the more specific pathway information
minf <- copy(st1)
m1 <- match(st1$rid, an$metab.annot$metabolite)
st1$metab.path <-  an$metab.annot$pathway[m1]

#' ### Add Pubchem from supp table

# Load Supp table for most recent metabolite hmdb/pubchem IDs
ids <- read_excel2(fn = "../170829_ConfounderAnalysis/paper/submission/revised_molecular_metabolism/Supplementary_tables.xlsx",
                   sheet = 2,
                   skip = 1)
ids$Abbreviation <- gsub(pattern = "*",replacement = "",x = ids$Abbreviation, fixed = TRUE)
ids$Abbreviation[ids$Abbreviation == "AC-total ()"] <- "AC-total"

#' # Find the overlap of the ids in cbd and pdb

# full overlap?
all(ids$Abbreviation %in% st1$abbr) # -> TRUE

m1 <- match(st1$abbr, ids$Abbreviation)

# hmdb2 is the one from the supp table
st1 <- cbind(st1, ids[m1, .(PubChem)])

# change wrong ID
st1[rid=="Aba", PubChem := "119"]

#' ### Re-Order and set labels

#' Include: Superpath, top gene, most enriched pathway (IPA), eta1, # sig genes, maybe overall unique gene?, weighted zero-inflation,
#' new percentage quotients

setnames(st1, 
         old = c("rid", "abbr", "full.name", "hmdb", "metab.class", "metab.class.q", 
                "metab.super.path", "i2.q99", "zero.inf", "fdr5hier", "eta1", 
                "markerID", "symbol_INGENUITY", "location_INGENUITY", "types_INGENUITY", 
                "betaREM", "ci.lREM", "ci.uREM", "seREM", "pREM", "r.squared", 
                "hfdr.1REM", "hfdr.2REM", "analyzed.adult.heart.ami.sorb", "metab.path"),
         new = c("R-ID", 
                "Metabolite Abbreviation", 
                "Full metabolite name",
                "HMDB", 
                "Metabolite Class", 
                "Metabolite Class 2", 
                "Metabolite Super Pathway", 
                "I2, 99th Percentile",
                "Weighted Metabolite Zero-Inflation",
                "Significant genes at hierarchical FDR=5%",
                "Eta1",
                "Top-Hit ILMN-ID",
                "Top-Hit Mapped Gene",
                "Top-Hit Mapped Gene Location",
                "Top-Hit Mapped Gene Function",
                "Top-Hit Effect Size",
                "Top-Hit CI-95% Lower Bound",
                "Top-Hit CI-95% Upper Bound",
                "Top-Hit Standard Error",
                "Top-Hit P-Value",
                "Top-Hit Explained Variance (R2)",
                "Top-Hit FDR Level 1",
                "Top-Hit FDR Level 2",
                "Single Study Effect (Order: Adult, Heart, AMI, Sorb)",
                "Metabolite Pathway"
                ))

# set new order
st1.order  <- c(
  "Metabolite Abbreviation", 
  "Full metabolite name",
  "HMDB", 
  "PubChem",
  "Metabolite Class", 
  "Metabolite Class 2", 
  "Metabolite Pathway",
  "Metabolite Super Pathway", 
  "Significant genes at hierarchical FDR=5%",
  "Eta1",
  "I2, 99th Percentile",
  "Weighted Metabolite Zero-Inflation",
  "Top-Hit Mapped Gene",
  "Top-Hit Mapped Gene Function",
  "Top-Hit Mapped Gene Location",
  "Top-Hit ILMN-ID",
  "Top-Hit Effect Size",
  "Top-Hit Standard Error",
  "Top-Hit CI-95% Lower Bound",
  "Top-Hit CI-95% Upper Bound",
  "Top-Hit Explained Variance (R2)",
  "Top-Hit P-Value",
  "Top-Hit FDR Level 1",
  "Top-Hit FDR Level 2",
  "Single Study Effect (Order: Adult, Heart, AMI, Sorb)"
)

# reorder
st1 <- st1[,..st1.order]

# set correct labels for metabolite class
st1[,`Metabolite Class`:=factor(
  `Metabolite Class`, 
  labels = c(
    "aa" = "Amino Acid",
    "ac" = "Acylcarnitine",
    "mix" = "Mix Quotient")
  )]
st1[,`Metabolite Class 2`:=factor(`Metabolite Class 2`, labels = c(
  "aa" = "Amino Acid",
  "ac" = "Acylcarnitine",
  "quotient" = "Quotient"
))]

WriteXLS_hk(st1, "ppr/supp_tables/metaboliteInfo.xls", SheetNames = "Supplementary Table 2")

# Supp Table 2 - Gene-centric overview ----

#' ## Supp Table 2 - Gene-centric overview

# the unique metabolite count is the important one
mrg <- dat[symbol_INGENUITY!="" & significant==TRUE, .(
  sig.metabs = uniqueN(metabolite),
  all.sig.metabs = paste0(unique(metabolite), collapse = ", "),
  sig.probes = uniqueN(markerID),
  all.probes = paste(unique(markerID), collapse=", ")
),
by=.(symbol_INGENUITY, 
     location_INGENUITY, 
     types_INGENUITY,
     biomarkers_INGENUITY,
     drugs_INGENUITY)]

# order from largest to smalles (top ones are my master-regulators)
mrg <- mrg[order(sig.metabs,decreasing = TRUE),]

# add the cumulative information again
cum.genes <- cumsum(mrg$sig.metabs)/sum(mrg$sig.metabs)
mrg[,cumulative.metab.assoc := cum.genes]
mrg[,perc.sig.genes := (1:length(mrg$sig.metabs))/length(mrg$sig.metabs)]

# get the UNIQUE number of AA/AC/Mix metabolites per master regulator
tmp1 <- dat[symbol_INGENUITY!="" & significant==TRUE, .(
  sig.metabs = unique(metabolite)),
by=.(symbol_INGENUITY)]

tmp1$class <- minf[match(tmp1$sig.metabs, minf$rid), metab.class]
tmp1[, num.aa.ac.mix := paste(sum(class == "aa"),
                              sum(class == "ac"),
                              sum(class == "mix"), sep = "/"),by=symbol_INGENUITY]

# now match to the mrg.unique
mrg$num.aa.ac.mix <- tmp1[match(mrg$symbol_INGENUITY, tmp1$symbol_INGENUITY), num.aa.ac.mix]

#' ### Add enriched pathways to table

setorder(gwe, p.hyper.fdr, pval.hyper)
top.enrich <- gwe[, .SD[pval.hyper==min(pval.hyper)], by = gene][, head(.SD,1),by = gene]
m1 <- match(mrg$symbol_INGENUITY, top.enrich$gene)
top.enrich <- top.enrich[m1, .(
  `Enriched Pathway` = name, 
  Enrichment = enrichment, 
  `P-Value (Hypergeom. Test)` = pval.hyper, 
  `Q-Value (FDR=5%)` = p.hyper.fdr
  )]

# add the enrichment info

mrg <- cbind(mrg, top.enrich)
# rename and clean
setnames(mrg, old = c(
  "symbol_INGENUITY", "location_INGENUITY", "types_INGENUITY", 
  "biomarkers_INGENUITY", "drugs_INGENUITY", "sig.metabs", "all.sig.metabs", 
  "sig.probes", "all.probes", "cumulative.metab.assoc", "perc.sig.genes", 
  "num.aa.ac.mix", "Enriched Pathway"),
  new = c(
    "Gene", 
    "Gene Expression Location", 
    "Gene Type", 
    "Biomarker Function of Gene",
    "Druggability of Gene",
    "Number of Significantly associating Metabolites",
    "Significantly associating Metabolites",
    "Number of mapped Probes",
    "Mapped Probes",
    "Cumulative fraction of total significant metabolite associations",
    "Cumulative fraction of total significantly associating genes",
    "Number of associating Amino Acids/Acylcarnitines/Mix-quotients",
    "Top Enriched Pathway"
    )
  )

# reorder for saving
st2.order <- c("Gene", 
              "Number of mapped Probes", 
              "Mapped Probes", 
               "Gene Expression Location", 
               "Gene Type", 
               "Biomarker Function of Gene", 
               "Druggability of Gene", 
               "Number of Significantly associating Metabolites",
              "Significantly associating Metabolites", 
              "Number of associating Amino Acids/Acylcarnitines/Mix-quotients",
              "Cumulative fraction of total significant metabolite associations", 
              "Cumulative fraction of total significantly associating genes",
              "Top Enriched Pathway", 
              "Enrichment", 
              "P-Value (Hypergeom. Test)", 
              "Q-Value (FDR=5%)"
              )
st2 <- mrg[, ..st2.order]
WriteXLS_hk(st2, "ppr/supp_tables/Supplementary_Table_2.xls", SheetNames = "Supplementary Table 2")

# In the Results Text mentioned Numbers ----

#' # In the Results Text mentioned Numbers

# total number of associations
dat[significant==TRUE, .N]

# total number of associating metabolites
dat[significant==TRUE, unique(metabolite)]
dat[significant==TRUE, unique(metabolite)]

# total number of associating genes
dat[symbol_INGENUITY != "", significant==TRUE, uniqueN(symbol_INGENUITY)]

# Pathway Enrichment ----

#' # Pathway Enrichment

st2[`Q-Value (FDR=5%)`<=0.05, uniqueN(Gene)]
st2[`Q-Value (FDR=5%)`<=0.05, mytable(`Gene Expression Location`)]
st2[`Q-Value (FDR=5%)`<=0.05, mytable(`Gene Type`)]

st2[`P-Value (Hypergeom. Test)`<=0.05, uniqueN(Gene)]
st2[`P-Value (Hypergeom. Test)`<=0.05, uniqueN(`Top Enriched Pathway`)]
st2[`P-Value (Hypergeom. Test)`<=0.05, table(`Top Enriched Pathway`)] %>% sort(decreasing = T)

#=============================================#
#+++++++++++++++++++++++++++++++++++++++++++++#

# Metabolite-pathway enrichment ----

st2[`Q-Value (FDR=5%)`<=0.05, ][order(`Q-Value (FDR=5%)`),table(`Gene Expression Location`)]
st2[`Q-Value (FDR=5%)`<=0.05, ][order(`Q-Value (FDR=5%)`),table(`Biomarker Function of Gene`)]
st2[`Q-Value (FDR=5%)`<=0.05, ][order(`Q-Value (FDR=5%)`),table(`Gene Type`)]

# PAPER: Amount of enriched pathways ----
st2[`P-Value (Hypergeom. Test)`<=0.05, .(Gene, `Enriched Pathway`)][,uniqueN(Gene)]
st2[`P-Value (Hypergeom. Test)`<=0.05, .(Gene, `Enriched Pathway`)][,unique(`Enriched Pathway`)]
st2[`P-Value (Hypergeom. Test)`<=0.05, .(Gene, `Enriched Pathway`)][,uniqueN(`Enriched Pathway`)]

# check main regulators
setorder(st2, -`Number of Significantly associating Metabolites`, `Q-Value (FDR=5%)`)
st2[1:50, sum(`Q-Value (FDR=5%)`<=0.05)]

# PAPER: number of significant metabolites in enriched genes ----
st2[ ,.(num.enrich = sum(`Q-Value (FDR=5%)`<=0.05), 
        `Number of Significantly associating Metabolites`), 
     by=Gene][
       !is.na(num.enrich), ][
         order(num.enrich, decreasing = T), ][
           num.enrich>=1, range(`Number of Significantly associating Metabolites`)
         ]

#+++++++++++++++++++++++++++++++++++++++++++++#
#=============================================#

#=============================================#
#+++++++++++++++++++++++++++++++++++++++++++++#

# Gene-wise-pathway enrichment ----

#' Using KEGG/GO function from Andreas/Holger -> in other script

#+++++++++++++++++++++++++++++++++++++++++++++#
#=============================================#

#=============================================#
#+++++++++++++++++++++++++++++++++++++++++++++#

# Role of BCL11A ----

#' # Investigate BCL11A associations

# Association to BCAA metabolism intermediaries
# https://www.sciencedirect.com/science/article/pii/S1079979615000947
dat[significant==TRUE & symbol_INGENUITY == "BCL11A", uniqueN(metabolite)]
dat[significant==TRUE & symbol_INGENUITY == "BCL11A", .N]

# how many associated metabolites are part of BCAA-metabolism? -> many
dat[significant==TRUE & symbol_INGENUITY == "BCL11A", table(metab.path)]
dat[significant==TRUE & symbol_INGENUITY == "BCL11A", table(metab.super.path)]

# https://www.sciencedirect.com/science/article/pii/S1079979615000947
# repression of BCL11A in definitive erythroid cells
# effect direction
dat[significant==TRUE & symbol_INGENUITY == "BCL11A", table(sign(betaREM))]

# explained variance of metabolite 
dat[significant==TRUE & symbol_INGENUITY == "BCL11A", r.squared]


#+++++++++++++++++++++++++++++++++++++++++++++#
#=============================================#

#=============================================#
#+++++++++++++++++++++++++++++++++++++++++++++#

# Strongest associating gene ----

#' # Describe the strongest association found

# ALAS2
dat[significant==TRUE,.SD, .SDcols= c(
  "metabolite",
  "symbol_INGENUITY",
  "metab.path",
  "metab.super.path",
  "r.squared",
  "betaREM",
  "pREM", 
  "hfdr.1REM", 
  "hfdr.2REM", 
  "hfdr.sigREM")][
    order(abs(betaREM), decreasing = TRUE)][
      1:10, ] %>% as.data.frame

# https://ghr.nlm.nih.gov/gene/ALAS2
# codes for ALA-Synthase, an important actor in heme production as the 
# first of eight steps

#' The five strongest associations were all with the gene ALAS2

#+++++++++++++++++++++++++++++++++++++++++++++#
#=============================================#

#=============================================#
#+++++++++++++++++++++++++++++++++++++++++++++#


# Co-Expression Metabolites ----

# sig correlations overall
all.combinations <- all.combinations[m1!=m2, ]
metab.coexp <- all.combinations[m1 != m2, .(
  total.pairs          = .N,
  total.sigs           = sum(q.rho <= 0.05),
  total.sigs.rel       = sum(q.rho <= 0.05)/.N,
  high.total.sigs      = sum(q.rho <= 0.05 & abs(rho)>=.9),
  high.total.sigs.rel  = sum(q.rho <= 0.05 & abs(rho)>=.9)/.N,
  within.sig.abs       = sum(q.rho <= 0.05 & super.path1==super.path2),
  within.sig.rel       = sum(q.rho <= 0.05 & super.path1==super.path2)/.N,
  between.sig.abs      = sum(q.rho <= 0.05 & super.path1!=super.path2),
  between.sig.rel      = sum(q.rho <= 0.05 & super.path1!=super.path2)/.N,
  high.within.sig.abs  = sum(q.rho <= 0.05 & super.path1==super.path2 & abs(rho)>=.9),
  high.within.sig.rel  = sum(q.rho <= 0.05 & super.path1==super.path2 & abs(rho)>=.9)/.N,
  high.between.sig.abs = sum(q.rho <= 0.05 & super.path1!=super.path2 & abs(rho)>=.9),
  high.between.sig.rel = sum(q.rho <= 0.05 & super.path1!=super.path2 & abs(rho)>=.9)/.N
  ), by=.(m1, path1, super.path1)]

# save for easier lookup
WriteXLS_hk(metab.coexp, "ppr/supp_tables/TMP_Table_3.xls", SheetNames = "TMP Table 3")

# high coexpression between pathways
setorder(metab.coexp, -high.between.sig.rel, -high.within.sig.rel)
setorder(metab.coexp, -high.total.sigs)
metab.coexp[,.(m1,path1,super.path1, high.total.sigs, high.between.sig.abs, high.within.sig.abs, high.total.sigs)]

# metabolites that act across pathways
c("Glut", "Q35", "Q15")

all.combinations[m1!=m2 & m1=="Q34" & q.rho<=0.05 & abs(rho)>=0.9, .(m1,m2,super.path2, rho)]
st1[rid%in%c("C0","C2", "C14", "C181", "C201", "acges"),.(rid,abbr,full.name)]

all.combinations[m1!=m2 & m1=="Glut" & q.rho<=0.05 & abs(rho)>=0.9, .(m1,m2,super.path2, rho)]
st1[rid%in%c("C10", "C12", "C141", "C201", "Q33"),.(rid,abbr)]

st1[rid%in%c("C12", "C141", "Q34"),.(rid,abbr)]

# coexpression summarised by super-pathway
super.path.coexp <- all.combinations[m1!=m2, .(
  total.pairs          = .N,
  total.sigs           = sum(q.rho <= 0.05),
  total.sigs.rel       = sum(q.rho <= 0.05)/.N,
  high.total.sigs      = sum(q.rho <= 0.05 & abs(rho)>=.9),
  high.total.sigs.rel  = sum(q.rho <= 0.05 & abs(rho)>=.9)/.N,
  within.sig.abs       = sum(q.rho <= 0.05 & super.path1==super.path2),
  within.sig.rel       = sum(q.rho <= 0.05 & super.path1==super.path2)/.N,
  between.sig.abs      = sum(q.rho <= 0.05 & super.path1!=super.path2),
  between.sig.rel      = sum(q.rho <= 0.05 & super.path1!=super.path2)/.N,
  high.within.sig.abs  = sum(q.rho <= 0.05 & super.path1==super.path2 & abs(rho)>=.9),
  high.within.sig.rel  = sum(q.rho <= 0.05 & super.path1==super.path2 & abs(rho)>=.9)/.N,
  high.between.sig.abs = sum(q.rho <= 0.05 & super.path1!=super.path2 & abs(rho)>=.9),
  high.between.sig.rel = sum(q.rho <= 0.05 & super.path1!=super.path2 & abs(rho)>=.9)/.N
), by=.(super.path1, m2)]

WriteXLS_hk(super.path.coexp, "ppr/supp_tables/TMP_Table_4.xls", SheetNames = "TMP Table 4")

# high coexpression within pathways

# same on super pathway 

super.path.coexp <- all.combinations[m1!=m2, .(
  total.pairs          = .N,
  total.sigs           = sum(q.rho <= 0.05),
  total.sigs.rel       = sum(q.rho <= 0.05)/.N,
  high.total.sigs      = sum(q.rho <= 0.05 & abs(rho)>=.9),
  high.total.sigs.rel  = sum(q.rho <= 0.05 & abs(rho)>=.9)/.N,
  within.sig.abs       = sum(q.rho <= 0.05 & super.path1==super.path2),
  within.sig.rel       = sum(q.rho <= 0.05 & super.path1==super.path2)/.N,
  between.sig.abs      = sum(q.rho <= 0.05 & super.path1!=super.path2),
  between.sig.rel      = sum(q.rho <= 0.05 & super.path1!=super.path2)/.N,
  high.within.sig.abs  = sum(q.rho <= 0.05 & super.path1==super.path2 & abs(rho)>=.9),
  high.within.sig.rel  = sum(q.rho <= 0.05 & super.path1==super.path2 & abs(rho)>=.9)/.N,
  high.between.sig.abs = sum(q.rho <= 0.05 & super.path1!=super.path2 & abs(rho)>=.9),
  high.between.sig.rel = sum(q.rho <= 0.05 & super.path1!=super.path2 & abs(rho)>=.9)/.N
  ), by=.(super.path1, super.path2)]

# heatmap
coexp.matrix <- dcast(super.path.coexp, super.path1~super.path2, value.var = "high.total.sigs.rel")
setDF(coexp.matrix, rownames(coexp.matrix))
cm <- as.matrix(coexp.matrix[,-1])
rownames(cm) <- coexp.matrix$super.path1

# complex heatmap
ComplexHeatmap::Heatmap(
  matrix = cm, 
  col=c("#E2E2E2", "#8E063B"),
  clustering_method_rows    = "complete", 
  clustering_method_columns = "complete", 
  heatmap_legend_param = list(title="Overlapping significant\nbeta-correlation")
)

# save test-wise pic
tiff(dpx("TEST_heatmap.tiff", "plt/"),
     width = 7, 
     height = 5, 
     res=300, 
     unit="in")
# par(mar=c(8.1, 4.1, 8.1, 2.1))
heatmap(cm, margins=c(15,10))
dev.off()

# Network ----

# coexpression summarised by super-pathway
network.dat <- all.combinations[m1!=m2, .(
  total.pairs          = .N,
  total.sigs           = sum(q.rho <= 0.05),
  total.sigs.rel       = sum(q.rho <= 0.05)/.N,
  high.total.sigs      = sum(q.rho <= 0.05 & abs(rho)>=.9),
  high.total.sigs.rel  = sum(q.rho <= 0.05 & abs(rho)>=.9)/.N,
  within.sig.abs       = sum(q.rho <= 0.05 & super.path1==super.path2),
  within.sig.rel       = sum(q.rho <= 0.05 & super.path1==super.path2)/.N,
  between.sig.abs      = sum(q.rho <= 0.05 & super.path1!=super.path2),
  between.sig.rel      = sum(q.rho <= 0.05 & super.path1!=super.path2)/.N,
  high.within.sig.abs  = sum(q.rho <= 0.05 & super.path1==super.path2 & abs(rho)>=.9),
  high.within.sig.rel  = sum(q.rho <= 0.05 & super.path1==super.path2 & abs(rho)>=.9)/.N,
  high.between.sig.abs = sum(q.rho <= 0.05 & super.path1!=super.path2 & abs(rho)>=.9),
  high.between.sig.rel = sum(q.rho <= 0.05 & super.path1!=super.path2 & abs(rho)>=.9)/.N
), by=.(super.path1, m2)]

# get the coexpression by superpathway data and exclude all connections with 0 overlapping sig metabolites
fp <- network.dat[high.total.sigs.rel!=0, ]

# remove duplicates as in A-B = B-A etc
for(i in 1:nrow(fp)){
  dbs <- fp[i, paste(sort(c(super.path1, m2)),collapse = "___")]
  data.table::set(x = fp, i = i, j = "doubles", value=dbs)
}

# create filter col
fp[super.path1 != super.path2 & !duplicated(doubles), filter.remove := TRUE]
fp[is.na(filter.remove), filter.remove:=FALSE]

# keep only the non-duplicated
fp <- fp[filter.remove != T, ]

# create unique colnames so the function doesn't get confused
# fp[,`:=`(
#   node1=paste0("M1: ", super.path1),
#   node2=paste0("M2: ", super.path2)
# )]

# in case of Super-path to path connection
# fp[,`:=`(
#   node1=paste0("SP: ", super.path1),
#   node2=paste0("P: ", path2)
# )]

source("../functions/coexpression_network_plot.R")
network_plot(assocResults = fp,
             node1 = "super.path1",
             node2 = "m2",
             rSquaredColumn = "high.total.sigs.rel"
) %>% visPhysics(enabled = TRUE)


visExport(graph = n1,
  name = dpx("superPathCoExpressionNetwork.png", "plt/"))


# Direction of Co-Regulation ----

# effect sign all significant assocs
all.combinations[m1!=m2 & q.rho <= 0.05, mytable(sign(rho))]

# nur in den starken korrelationen
top.corr <- all.combinations[abs(rho)>=.9 & m1!=m2 & q.rho <= 0.05, ]

# sign of the top metabolite
top.corr[, mytable(sign(rho))] 

# Pathway-Co-Regulation-Network 


#+++++++++++++++++++++++++++++++++++++++++++++#
#=============================================#
