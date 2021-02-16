#' ---
#' title: "Explore primary analysis"
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
  "ggplot2",
  "plotly",
  "circlize",
  "ComplexHeatmap",
  "magrittr",
  # "dendextend",
  "RColorBrewer",
  "corrplot"
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
source(here("../functions/correlate_meta_statistics.R"))
source(here("../functions/option_setup.R"))
source(here("../functions/beta_plot.R"))
source(here("../functions/all_annotation_tables.R"))
source(here("../functions/calc_quotients.R"))
source(here("../functions/quotient_calc_table.R"))
source(here("../functions/theme_gxMetab.R"))
an <- all_annotation_tables(mount.point = "/net/ifs1/san_projekte/projekte/genstat/")
option_setup()

#+ include=F
#==========#
#   DONE   #
#==========#

# Load metabolites ----

#' I evaluated annotated data in the previous step.

nf <- newest_file(look_for = "STEP4_meta_resultsAnnotated",subfolder = "res", print_full = TRUE)
dat <- fread(nf)

# Fix wrong mix quotient label!
dat[metabolite=="Q6", metab.class := "mix"]

#' # Annotate quotient numerator/denominator

# qt = quotient table
qt <- quotient_calc_table()

# read in supp table 3 from previous paper for proper metabolite names
# st = supplementary table
st <- read_excel2("../170829_ConfounderAnalysis/paper/submission/revised_molecular_metabolism/Supplementary_tables.xlsx", sheet = 2, skip = 1)

# cn = clean names
cn <- gsub(st$Abbreviation, pattern = "*", replacement = "", fixed = TRUE) 
cn <- gsub(x=cn, pattern = " ()", replacement = "", fixed = TRUE)

# agn = all good names
agn <- data.table(
  abbr = cn,
  full.name = st$Name,
  hmdb = st$HMDB,
  rid = c("C0", "C10", "C101", "C12", "C14", "C141", "C14OH", "C16", 
          "C161", "C161OH", "C16OH", "C18", "C181", "C181OH", "C182", 
          "C182OH", "C18OH", "C2", "C201", "C202", "C203", "C3", "C3DC", 
          "C4", "C4OH", "C5", "C51", "C5OHHMG", "C6", "C6DC", "C8", "C81", 
          "Glut", "MeGlut", "MMA", "acges", "Aba", "Ala", "Arg", "Asn", 
          "Asp", "Carnosin", "Cit", "Gln", "Glu", "Gly", "His", "LeuIle", 
          "Lys", "MeHis", "Met", "OHProl", "Orn", "Phe", "PiPA", "Pro", 
          "Sarc", "Ser", "Tau", "Thr", "Trp", "Tyr", "Val"))

# Do I have the full overlap?
venn2(unique(dat$metabolite),agn$rid) %>% invisible

# Add quotient definitions to agn
agn.q <- an$metab.annot[grep("^Q",metabolite), .(metabolite, quotient.def, full.name, pheno.class, super.pathway)]
agn.q[pheno.class=="",pheno.class:="mix"]

# merge the name dts
setnames(agn.q, 
         old = c("metabolite", "quotient.def"),
         new = c("rid", "abbr"))

# rbind
agn.all <- rbindlist(list(agn[,.(rid, abbr, full.name, hmdb)], agn.q[,.(rid, abbr, full.name)]), fill = TRUE)

# save for use in step6
# save(list = c("agn.all", "agn", "agn.q"),file = dpx("metaboliteNames.RData", "obj/"))
load("obj/200306_metaboliteNames.RData")

#' # How does the r squared distribute over the metabolites

# plot template
pd <- position_dodge(.5)
theme_set(theme_light(base_size = 20, base_family = "Helvetica"))

# Calculate correlation of association lists ----
#' # Compare Effect-size distributions per metabolite pair

#' first Plot the betas of metabolites

all.combinations <- 
  expand.grid(m1=unique(dat$metabolite),
              m2=unique(dat$metabolite),
              stringsAsFactors = FALSE)
setDT(all.combinations)
all.combinations[,index:=1:.N]

#+ eval=F
tiff("plt/meta_betaPlots.tiff", width=8, height=7, res=300, units = "in")
# effect size correlations? by gene or metabolite
for(i in all.combinations$index){
  i1 <- all.combinations[index==(i),m1]
  i2 <- all.combinations[index==(i),m2]
  if(length(dat[numberStudies>1 & metabolite == (i1), .N]) == 
     length(dat[numberStudies>1 & metabolite == (i2), .N])){
    # all betas
    # source("../functions/beta_plot.R")
    p.beta <- beta_plot(
      x     = dat[numberStudies>1 & metabolite == (i1), betaREM],
      y     = dat[numberStudies>1 & metabolite == (i2), betaREM],
      xQval = dat[numberStudies>1 & metabolite == (i1), qREM],
      yQval = dat[numberStudies>1 & metabolite == (i2), qREM],
      showPlot=FALSE)
    # plot
    print(
      p.beta$plot + ggtitle(paste0(i1, "~",i2)) + xlab(i1) + ylab(i2)
    )
  }
}
dev.off()

#' The beta plots will not be computed for every combination, but I need to calculate the statistics

#+ eval=F
# this function only correlates FDR 5% (hierarchical) significant associations in either x or y
all.combinations <- correlate_meta_statistics(
  assocDat=dat,
  statistic="betaREM",
  # controlFDR = 0.05,
  useCores = 10,
  verbose = TRUE)

#' MeGlut~MeGlut Fails due to only having two significant associations. However,
#' the associations are identical, therefore I'll set them manually. This also
#' applies to every other m1==m2 situation.

#+ eval=F
all.combinations[m1==m2 & !is.na(comment), `:=`(
  rho=1,
  p.rho=0,
  beta=1,
  r2=1,
  signed.r2=1,
  n=0,
  p.lm=0,
  q.rho=0,
  q.lm=0)]

#+ eval=F
# set everything that failed due to sample size etc to 0 (will be filtered later by circlize)
all.combinations[!is.na(comment) & m1 != m2, `:=`(
  rho=0,
  p.rho=1,
  beta=0,
  r2=0,
  signed.r2=0,
  n=2,
  p.lm=1,
  q.rho=1,
  q.lm=1
)]

#+ eval=F
# save for later use
save_csv_carl(file = all.combinations,
              ile_name = "meta_betaCorrelationR2",
              ubfolder = "res",
              sep = ",",
              quote = TRUE)

#' Load the results from a previous run due to long runtime and plotting.
nf <- newest_file("meta_betaCorrelationR2", subfolder = "res",print_full = TRUE)
all.combinations <- fread(nf)

# ============================================= #

# Prepare data for the circular plot ----
#' # Circular Plot of r-squared of each metabolite-combination beta-comparison

# define groups by metabolite class aa, ac, q
group <- dat[,unique(metabolite), by=.(metab.class, metab.class.q,metab.super.path)]

#' ## Set correct names and groups all objects

# match the metabolite names
group[,metabolite:=agn.all[match(group$V1, agn.all$rid), abbr]]

# match the correct names to the results from the beta correlation
all.combinations[, `:=`(
  metabolite1 = group[match(m1, V1), metabolite],
  metabolite2 = group[match(m2, V1), metabolite]
)]

dat$metab.path %>% unique
dat$metab.super.path %>% unique

# Change ACGes to sum as an extra group
group[metab.super.path=="",metab.super.path:="Other"]

# Better group label
group$metab.class.long <- factor(
  group$metab.class, 
  levels = c("aa",
             "ac",
             "mix"),
  label = c("Amino Acids", 
            "Acylcarnitines",
            "Mix"))

#' Filter step: which correlations are to be plotted! Set everything zero prior to building the correlation matrix.

# decide on FDR cutoff (or rho cutoff, see what works)
all.combinations[, rho.filtered := rho]

#' Set the correlation cutoff
rho.cutoff <- .90

#' Version 1 is filtering by q-value, e.g. for FDR at 5%
all.combinations[q.rho>0.05, rho.filtered := 0]

#' Version 2 is filtering by absolute correlation
all.combinations[abs(rho)<(rho.cutoff), rho.filtered := 0]

# remove all combinations where rho.filtered=0
# ac.backup <- all.combinations
# all.combinations.filtered <- all.combinations[m1!=m2 &
# q.rho<=0.05 &
# abs(rho.filtered)>=rho.cutoff, ]

# sanity check
all.combinations[,range(abs(rho.filtered))]
all.combinations[,max(q.rho)]
all.combinations[m1==m2,.N]

#' Info is taken from here: https://jokergoo.github.io/blog/html/large_matrix_circular.html

#==========================================================#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

# External info for circos plot ----
# Build object containing the information for the outer circle

# get the zero-inflation weighted by total N
n.sorb <- 935 # extracted from ESet prior to LIMMA Analysis
n.heart <- 1271 + 2355 # ami + heart 
n.adult <- 3145
n.total <- n.sorb+n.heart+n.adult

# get the weighted zero-inflation for all cohorts
zinf <- an$metab.annot[, .(metabolite, percent.0.a1, percent.0.b3, percent.0.sorb)]
zinf[, zinf.weighted := ((n.adult / n.total) * percent.0.a1 + (n.heart / n.total) * percent.0.b3 + (n.sorb / n.total) * percent.0.sorb )]

# check heterogeneity of significant assocs
plot.i2 <- dat[significant==TRUE & metabolite %in% group$V1, .(
  i2.q50 = quantile(I2, probs = 0.50),
  i2.q75 = quantile(I2, probs = 0.75),
  i2.q95 = quantile(I2, probs = 0.95),
  i2.q99 = quantile(I2, probs = 0.99),
  fdr5hier = uniqueN(symbol_INGENUITY)
), by=.(metabolite,metab.class, metab.super.path)]

# get eta1 of each metabolite
plot.info <- dat[metabolite %in% group$V1, .(
  eta1 = 1-fdrtool::pval.estimate.eta0(pREM, diagnostic.plot = FALSE)
), by=.(metabolite,metab.class,metab.super.path)]

# match I2 info to the rest
venn2(plot.info$metabolite, plot.i2$metabolite) # all metabolites with significant assocs
m1 <- match(plot.info$metabolite, plot.i2$metabolite)
new.cols <-plot.i2[m1, .(i2.q50, i2.q75,i2.q95,i2.q99, fdr5hier)]
plot.info[, names(new.cols):=new.cols]

# set the missing info to zero (data.table removes empty entries when summing)
plot.info[is.na(fdr5hier), fdr5hier:= 0]

# Change ACGes to sum as an extra group
# plot.info[metabolite=="acges",metab.class:="sum"]
plot.info[metab.super.path=="",metab.super.path:="Other"]
plot.info$metab.class <- factor(
  plot.info$metab.class,
  levels = c("aa",
             "ac",
             "mix"),
  label = c("Amino Acids",
            "Acylcarnitines",
            "Mix"))

# order data
setorder(plot.info, metab.super.path, -fdr5hier, metab.class)

# set set correct names
plot.info$rid <- plot.info$metabolite

# first get the group dt the missing quotient names in the V1 colum
# group[metab.class.q=="quotient", metabolite:=V1]
# group[, metabolite:=V1]

# set the proper names insteat of the RID
plot.info[, metabolite := group[match(rid, V1), metabolite]]

# add zero inflation to plot.info
m2 <- match(plot.info$rid, zinf$metabolite)
plot.info$zero.inf <- zinf[(m2), zinf.weighted]

# prepare plot.info as Supp Table
supp.table1 <- plot.info[,.(metabolite, metab.class, metab.super.path, i2.q99, zero.inf, fdr5hier, eta1)]

# get the top gene for each metabolite 
top.gene <- dat[significant==TRUE, .SD[pREM==min(pREM)],
                by=metabolite, 
                .SDcols = c("markerID", 
                            "symbol_INGENUITY",
                            "location_INGENUITY",
                            "types_INGENUITY", 
                            "betaREM", 
                            "ci.lREM",
                            "ci.uREM",
                            "seREM",
                            "pREM",
                            "r.squared",
                            "hfdr.1REM",
                            "hfdr.2REM",
                            "analyzed.adult.heart.ami.sorb"
                            )]

# match the three dt agn.all, top.gene, supp.table1
m1 <- match(agn.all$abbr, supp.table1$metabolite)
m2 <- match(agn.all$rid, top.gene$metabolite)
m3 <- match(agn.all$rid, dat$metabolite)

# join the info
supp.table.1 <- cbind(
  agn.all,
  dat[m3, .(metab.class, metab.class.q)],
  supp.table1[m1, .SD,.SDcols=names(supp.table1)[-2:-1]],
  top.gene[m2, .SD,.SDcols=names(top.gene)[-1]]
)

# Save data for Supp. Table ----

fwrite(supp.table.1, "obj/201112_metabOverview.csv")

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#==========================================================#

# Create matrix from correlations ----
# first create a matrix of observations for each pair
corr.mat <- dcast(all.combinations, metabolite1~metabolite2, value.var = "rho.filtered")
setDF(corr.mat,rownames=corr.mat$metabolite1)
corr.mat$metabolite1<-NULL
mat<-as.matrix(corr.mat)

#' I need to reorder the matrix for it to correctly display all the groups

# match to ordered reference and reorder
m1 <- match(plot.info$metabolite, rownames(mat))
mat <- mat[m1,m1]
rownames(mat)

# also order the group reference by the plot info
m2 <- match(plot.info$metabolite, group$metabolite)
group <- group[m2, ]

# remove the diagonal and lower triangle, because I only need the pairwise
# correlations once.
mat[is.na(mat)] <- 0

#' Check a regular heatmap
corrplot(mat,
         method="circle",
         #method="number",
         number.cex=.3,
         tl.cex = .5,
         type="full",
         order="hclust")

diag(mat) = 0
mat[lower.tri(mat)] <- 0

# get index of each row and number of 
n = nrow(mat)
rn = rownames(mat)

#===============================================# 
#' Create the matrix with the full matrix and save for use in the next step

# first create a matrix of observations for each pair
corr.mat.2 <- dcast(all.combinations, metabolite1~metabolite2, value.var = "rho")
setDF(corr.mat.2,rownames=corr.mat.2$metabolite1)
corr.mat.2$metabolite1<-NULL
mat.2<-as.matrix(corr.mat.2)

#' I need to reorder the matrix for it to correctly display all the groups

# match to ordered reference and reorder
m1 <- match(plot.info$metabolite, rownames(mat.2))
mat.2 <- mat.2[m1,m1]
write.csv(mat.2, file = dpx("fullBetaCorrMatrix.csv","obj/"), quote = TRUE)

#===============================================# 

#' Create the sector groups for the metabolites

# create list object with group names and metabolites in each group
gl <-sapply(unique(group$metab.super.path), function(x){
  
  # sda
  # res <- list(group[metab.class.q==(x), V1])
  res <- group[metab.super.path==(x), metabolite]
  res <- list(res)
  names(res) <- x
  return(res)
},USE.NAMES=FALSE)

# preview
str(gl)
gl <- gl[
  c("Amino acid metabolism, other", 
    "Folic acid cycle", 
    "Ammonia recycling", 
    "BCAA metabolism", 
    "Collagen synthesis", 
    "Energy metabolism", 
    "Carnitine transport", 
    "Fatty acid metabolism", 
    "Other", 
    "Urea cycle"
  )]

# get entries per group
group_size <- sapply(gl, length)
# group_size <- group[,uniqueN(metabolite), by=metab.super.path]$V1

# reformat the gl list into named vector and generate color for each group
gd <- structure(rep(names(gl), times = sapply(gl, length)), names = unlist(gl))

# create a group indicating class ownership for each metabolite
gd.class <- as.character(group$metab.class.long)
names(gd.class) <- group$metabolite

# color for metabolite groups
mclass_color <- structure(
  c("darkgoldenrod", "darkolivegreen4", "cadetblue4"), # stay in line with quotient plot
  # c("cornflowerblue", "darkmagenta", "aquamarine3"), 
  # c("#027C1E", "#97A753", "#DEC79D", "#E2E2E2"),
  names = unique(as.character(group$metab.class.long)))

# color for superpathway --- stay in line with heatmap coloring
superpath_color <- structure(
  c("#E16A86", "#CB7F2F", "#9F9400", 
    "#50A315", "#00AC79", "#00AAB7", 
    "#009ADE", "#A87BE4", "#DA65C3", 
    "#E16A86"),
  names=names(gl))

# this is prettier, but harder to read later on!
# superpath_color <- structure(
#   c("#7D0112", "#91331B", "#A55229", 
#     "#B76F39", "#C78A4D", "#D5A463", 
#     "#E1BD7C", "#EAD396", "#F0E6B3", 
#     "#F2F1E4"),
#   names=names(gl))

# set the points for each metabolite indicating heterogeneity
# heterogeneity: Non, good: 0-.25, ok: 0.25-0.75, bad: 0.75-1 as
col.het <- c("green4", "gold", "gold3", "indianred4")
col.fun.het <- colorRamp2(c(0, 0.25, 0.75, 1), 
                          c("green4", "gold", "gold3", "indianred4"), 
                          transparency = 0.0)

# number of sectors
n_group <- length(gl)

# set color for the connections (e.g. like in corrplot)
col_fun <- colorRamp2(c(-1, -0.5, 0.5, 1), 
                      c("dodgerblue4", "white", "white", "indianred4"), 
                      transparency = 0.4)


# give the points the coloring and size by eta1
col.eta <- c("indianred4","gold", "gold3", "green4")
col.fun.eta <- colorRamp2(c(0, .2, 0.4), 
                          c("indianred4","gold", "green4"), 
                          transparency = 0)

# give the points the coloring and size by eta1
col.zinf <- c("green4", "gold","gold3", "indianred4")
col.fun.zinf <- colorRamp2(c(0, 0.5, 1), 
                          c("green4", "gold", "indianred4"), 
                          transparency = 0)

#' # Build circular plot

col.border.light <- "grey50"
col.border.dark <- "grey45"
col.text.dark <- "grey26"
col.text <- "grey30"
col.bg.dark <- "grey97"
col.bg.light <- "grey98"

# check plot parameters
# cp.back <- circos.par()
circos.info()

# Circos Plot function ----
#' # Build the Plot from the function

#' I put the code for the plot into a function for easier debugging

# build plot
source("../functions/build_circos.R")

### Start Plot
tiff(dpx("circos_600dpi.tiff", "plt/"), 
     width = 10, 
     height = 10, 
     res=600, 
     unit="in")

# Plotting function to be found in the functions-folder
build_circos()

dev.off() ### End Plot


#=============================================#
#+++++++++++++++++++++++++++++++++++++++++++++#

# Prepare data for quotient plot ----
#' # Do quotients add information to the regular AAs/ACs

#' Add quotient definition to the dat in order to connetct them in the plot. For this, I wil create groups for each quotient (named after the Q) that contain all the metabolites that make out the quotient in order to tell ggplot which metabolites to link to which quotients

# get the total number of significant genes
all.metabs <- dat[,.N,by=.(metabolite, metab.class, metab.class.q, metab.super.path)]
all.metabs$N <- NULL
sig.genes <- dat[significant==TRUE,.(sig.genes = uniqueN(symbol_INGENUITY)),by=.(metabolite, metab.class,metab.super.path)]

# merge
m1 <- match(all.metabs$metabolite, sig.genes$metabolite)
all.metabs[, sig.genes := sig.genes[(m1), sig.genes]]
all.metabs[is.na(sig.genes), sig.genes:= 0]

# set zeros in the base table to NA (syntax is easier then)
change_in_dt(qt,from = 0,to = NA,change_in_dat = TRUE) %>% invisible()

# melt the data 
quotient.info <- melt(qt, 
                      id.vars = "Quotient", 
                      measure.vars = c("Quotient", "Numerator1", "Numerator2", "Denominator1", "Denominator2"), 
                      variable.name = "quotient.def",
                      value.name = "metabolite")

# remove all NA entries (empty second numerator/denominator etc)
quotient.info <- quotient.info[!is.na(metabolite), ]

# join the quotient info to the all.metabs info
m2 <- match(quotient.info$metabolite, all.metabs$metabolite)

# add no sig genes to quotient info
quotient.info[, c("sig.genes", 
                  "metab.class",
                  "metab.class.q",
                  "super.path") := all.metabs[(m2), .(sig.genes, 
                                                      metab.class, 
                                                      metab.class.q, 
                                                      metab.super.path)]]

# create object for remaining metabolites not in any quotient
# remaining metabolites
rem.met <- quotient.info$metabolite %>% unique
quotient.add <- all.metabs[metabolite %nin% rem.met, ]
quotient.add[, Quotient := metabolite]
quotient.add[, quotient.def := "none"]
setnames(quotient.add, "metab.super.path", "super.path")

# join all the information
plot.dat <- rbindlist(list(quotient.info, quotient.add), use.names = TRUE)
setnames(plot.dat, old = "Quotient", "group")

# add quotient definition
m3 <- match(plot.dat$metabolite, agn.q$rid)
plot.dat[, full.name := agn.q[m3, abbr]]
m4 <- match(plot.dat[is.na(full.name), metabolite], agn$rid)
plot.dat[is.na(full.name), full.name := agn[m4, abbr]]

#' It would be interesting whether the quotients associate with different genes than the metabolites they are computed from.

# get the overlap of significant genes
for(i in plot.dat[quotient.def=="Quotient", unique(group)]){
  
  # get the different parts of the quotient for comparison of unique genes
  # quotient = qtt 
  qtt <- i
  # prt = parts
  prt <- plot.dat[group == (qtt) & quotient.def != "Quotient", metabolite]
  
  # get the unique genes in each of the quotient parts
  # msg = metabolite sig genes
  msg <- dat[significant==TRUE & metabolite %in% (prt), symbol_INGENUITY]
  
  # remove empty entries
  msg <- msg[!(msg %in% "")]
  
  length(msg)
  
  # mug = metabolite unique genes
  mug <- unique(msg)
  
  # qsg = quotient sig gines
  qsg <- dat[significant==TRUE & metabolite %in% (qtt), symbol_INGENUITY]
  qsg <- qsg[!(qsg %in% "")]
  
  # overlap of quotient/metabolite genes
  vn2 <- venn2(mug, qsg, plotte = FALSE)
  
  # percentage of metabolite genes that do not associate with the metabolite parts
  # nqg = new quotient genes
  nqg <- length(vn2$q3) / uniqueN(qsg)

  # enter into plot.dat
  plot.dat[quotient.def == "Quotient" & group == qtt, new.perc := nqg]
  plot.dat[quotient.def == "Quotient" & group == qtt, sigs.sum := uniqueN(mug)]

  # plot.dat[quotient.def != "Quotient" & group == qtt, joint.genes := length(vn2$q1)]
  # plot.dat[quotient.def != "Quotient" & group == qtt, total.genes := length(qsg)]
  # plot.dat[quotient.def != "Quotient" & group == qtt, total.genes := length(qsg)]
  
}

# set all non quotient entries to 0
plot.dat[is.na(new.perc), new.perc:=0]
plot.dat

# check the sigs.sum
# plot.dat[, unique(sigs.sum[!is.na(sigs.sum)]),group] %>% as.data.frame

# create a new group for correct plotting: color for each metabolite, proper links of each Quotient
# line-link group with duplicated Quotient for each quotient part for separate lines to be joined to the quotient
plot.dat[, link.group := ifelse(group==metabolite, NA, paste(quotient.def,group,sep="_"))]
plot.dat <- rbind(plot.dat, plot.dat[quotient.def =="Quotient"])
setorder(plot.dat, group, link.group, na.last = T)

# long quotients
lq <- plot.dat[quotient.def %in% c("Denominator2", "Numerator2"),unique(group)]
plot.dat[group %nin% (lq), link.group2 := ifelse(is.na(link.group)==T, link.group[is.na(link.group)==F], link.group), by = group]

# bind enough rows of each remaining quotient to plot.dat
plot.dat[group %in% lq,.N,by=group]

# twice for these quotients
plot.dat <- rbind(plot.dat, 
                  plot.dat[group == lq[1] & quotient.def=="Quotient", ][1],
                  plot.dat[group == lq[1] & quotient.def=="Quotient", ][1])
plot.dat <- rbind(plot.dat, 
                  plot.dat[group == lq[2] & quotient.def=="Quotient", ][1])
plot.dat <- rbind(plot.dat, 
                  plot.dat[group == lq[3] & quotient.def=="Quotient", ][1])
plot.dat[group %in% (lq), link.group2 := ifelse(is.na(link.group)==T, link.group[is.na(link.group)==F], link.group), by = group]

# non-empty group for everything else that prevents from linking NAs as a group
plot.dat[is.na(link.group2), link.group2 := paste0("dummy_", 1:.N)]

# should be at most 2
plot.dat[,.N,by=link.group2] %>% as.data.frame

# set empty supercass to other
plot.dat[super.path=="", super.path:="Other"]

# color the links when quotient has more associations than the seperate metabolites 

# get the sum of the individual metabs for each quotient

# 201215 moved this into the loop to get unique sum of genes
# plot.dat[quotient.def %nin% c("Quotient", "none"), sigs.sum := sum(sig.genes,na.rm = TRUE),by=group]

plot.dat[,unique(sigs.sum),by=group] %>% as.data.frame

# set sigs.sum the same value for everything in the group
plot.dat[,sigs.sum:=unique(sigs.sum[!is.na(sigs.sum)]),by=group]
plot.dat[quotient.def=="Quotient", group.col := ifelse(sig.genes > sigs.sum, "more","less")]
plot.dat[,group.col:=unique(group.col[!is.na(group.col)]),by=group]

# new column for more/less in at least one metabolite

quot.sigs <- plot.dat[quotient.def=="Quotient", .(metabolite, sig.genes)][!duplicated(metabolite), ]
plot.dat[, quotient.sigs := quot.sigs$sig.genes[match(plot.dat$group, quot.sigs$metabolite)]]

# T/F wheter the the quotient has more sig. genes than the single metabolite
plot.dat[quotient.def != "Quotient", quot.sigs := quotient.sigs > sig.genes]

#' ## Plot significant number of associations

plot.dat$metab.class.qf <- factor(plot.dat$metab.class.q, 
                                  ordered = TRUE, 
                                  levels = c("aa", "quotient", "ac"),
                                  labels= c("Amino Acid", "Quotient", "Acylcarnitine"))
plot.dat$metab.class <- factor(plot.dat$metab.class, 
                               levels = c("aa", "mix", "ac"),
                               labels = c("Amino Acid", "Mix", "Acylcarnitine"))

#' ## How large is the percentage of quotient genes that is not significant in the numerator/denomiator

# 201123 New genes ----

#' I need to select quotients with highest new.perc

# show those with >50% new genes
plot.dat[new.perc>=0.5, .(full.name,sig.genes, new.perc, metab.class)][order(-new.perc),] %>% unique %>%  as.data.frame


#' Create a base-R plot

base.plot.dat <- plot.dat[!duplicated(metabolite), ]

# numeric x axis
base.plot.dat[, group.num := ifelse(metab.class.q == "aa", 0.5, ifelse(metab.class.q == "ac", 2.5, 1.5))]
base.plot.dat[metab.class.q=="quotient", group.num := ifelse(metab.class == "Amino Acid", group.num - 0.3, 
                                                             ifelse(metab.class == "Acylcarnitine", group.num + 0.3,
                                                                    group.num))]
# create new col with jittered values
base.plot.dat[ , group.num.jitter := group.num]
base.plot.dat[metab.class.q != "quotient", group.num.jitter := jitter(group.num, .2)]
base.plot.dat[metab.class.q == "quotient", group.num.jitter := jitter(group.num, 1)]

# get log values, set log(1)=0 to 1
base.plot.dat[, sig.genes.log := log10(sig.genes)]

# set log(0)=-Inf to 0
base.plot.dat[is.infinite(sig.genes.log), sig.genes.log := -0.5]

# log values also needed from this file for drawing of the lines
plot.dat[, sig.genes.log := log10(sig.genes)]

# set the log(0) -> -Inf to Zero 
plot.dat[is.infinite(sig.genes.log), sig.genes.log := -0.5]

# Quotient plot ----
#' ### Try base-plotting

# cool way to show palette:
# scales::show_col(pal)


# Plotting base
source("../functions/quotient_plot.R")

tiff(
  paste0("plt/", format(Sys.Date(), "%g%m%d"), "_quotient_300dpi.tiff"), 
  width = 10, 
  height = 10, 
  res=300, 
  unit="in"
)

# draw the plot
quotient_plot()

dev.off()

# Quotient Paper Info ----

#' ## Whats special about the quotient that has hundreds of significant genes more? (the red ones in the plot)

# Which ones are those?
base.plot.dat[group %in% c("Q26", "Q22"), ]

# Q26
base.plot.dat[metabolite %in% c("Q26", "Glu", "Orn")]

# Q22
base.plot.dat[metabolite %in% c("Q22", "Ser", "Met")]

base.plot.dat[group.col =="more", ]
plot.dat[group.col =="more", ]

plot.dat[quot.sigs==TRUE, .(group, metabolite, sigs.sum, quotient.sigs, metab.class)] %>% as.data.frame
plot.dat[quot.sigs==TRUE, .(group, metabolite, metab.class)]$metab.class %>% table 

#' ## Which quotients add more info than at least one of their parts

# metabolites that have less sig genes than their quotients
plot.dat[quot.sigs==TRUE & group.col!="more",metabolite]
plot.dat[quot.sigs==TRUE & group.col!="more", .(group, metab.class)]
plot.dat[quot.sigs==TRUE & group.col!="more",uniqueN(group), by=metab.class]

#' ## Additional Info added by quotients on AC 

plot.dat[group=="Q34", ] 

# groups in which more genes associate than in at least one part
more.groups <- plot.dat[quot.sigs==TRUE & group.col!="more", unique(group)] 
base.plot.dat[quotient.def=="Quotient" & group %in% more.groups, .N, by=(metab.class)] 

#' ## Info on new percentage of Quotients

base.plot.dat[metab.class.q == "quotient", mean(new.perc),by=metab.class]
base.plot.dat[metab.class.q == "quotient", median(new.perc),by=metab.class]

base.plot.dat[order(-new.perc) & new.perc==1,] %>% head(10)
agn.all[rid%in%c("Q14", "Q15", "Q28", "Q6"), abbr]

#+++++++++++++++++++++++++++++++++++++++++++++#
#=============================================#

#' # Session info

devtools::session_info()

