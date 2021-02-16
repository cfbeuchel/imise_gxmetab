#' ---
#' title: "Explore meta-analysis results"
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
  "ggplot2",
  "ggridges",
  "plotly",
  "visNetwork",
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
source(here("../functions/option_setup.R"))
source(here("../functions/beta_plot.R"))
source(here("../functions/all_annotation_tables.R"))
source(here("../functions/calc_quotients.R"))
source(here("../functions/quotient_calc_table.R"))
an <- all_annotation_tables(mount.point = "/net/ifs2/san_projekte/projekte/genstat/")
option_setup()

#+ include=F
#==========#
#   DONE   #
#==========#

#' I evaluated meta-analysis approaches in the previous step. I decided on the double-adjusted approach for the main analysis.  

#' # Load Data from STEP3

load("obj/190729_STEP3.RData")

#' How many unique probe Ids are in the result-table
res.meta$markerID %>% uniqueN

#' # Check on REM/FEM differences

#' This looks like the heterogeneity for REM significant genes is often lower than .75

#' ## REM Significancies

# plot I2 over metabolites for different significance filters
i2.probe.sum <- res.meta[hfdr.sigREM==TRUE & numberStudies >= 2, .(mean.i2 = mean(I2),
                             median.i2 = median(I2),
                             sd.i2 = sd(I2)), by = .(markerID)]

head(i2.probe.sum)

i2.metab.sum <- res.meta[hfdr.sigREM==TRUE & numberStudies >= 2, .(mean.i2 = mean(I2),
                             median.i2 = median(I2),
                             sd.i2 = sd(I2)), by = .(metabolite)]

head(i2.metab.sum)

#' ## I2 plotting

#' Venn Diagram of Significancies (probes) based on REM, FEM & FEM with I2 cutoff

m1 <- res.meta[hfdr.sigFEM==TRUE, uniqueN(markerID), by=metabolite]
m20 <- res.meta[hfdr.sigFEM==TRUE & I2 < .75, uniqueN(markerID), by=metabolite]
m21 <- res.meta[hfdr.sigFEM==TRUE & I2 < .5, uniqueN(markerID), by=metabolite]
m22 <- res.meta[hfdr.sigFEM==TRUE & I2 < .25, uniqueN(markerID), by=metabolite]
m3 <- res.meta[hfdr.sigREM==TRUE, uniqueN(markerID), by=metabolite]

# start a new table to match the above summaries to
sigs <- data.table(
  metabolite= unique(res.meta$metabolite)
  )

# match
sigs$FEM.noI2 <- m1[match(sigs$metabolite, m1$metabolite), V1]
sigs$FEM.75.I2 <- m20[match(sigs$metabolite, m20$metabolite), V1]
sigs$FEM.50.I2 <- m21[match(sigs$metabolite, m21$metabolite), V1]
sigs$FEM.25.I2 <- m22[match(sigs$metabolite, m22$metabolite), V1]
sigs$REM <- m3[match(sigs$metabolite, m3$metabolite), V1]

# set zeros
change_in_dt(dat = sigs, from = NA, to = 0, change_in_dat = TRUE) %>% invisible()

# show in final doc
kable(sigs)

# save
save_csv_carl(file = sigs, file_name = "meta_remFemI2Comparison",subfolder = "obj")

# annotate significant by metabolite fem/rem
sigs.m <- melt(sigs, id.vars = "metabolite", variable.name = "approach", value.name = "number.sigs")
sigs.m[,ord:=max(number.sigs),by=approach]
sigs.m[, class := an$metab.annot[match(sigs.m$metabolite, an$metab.annot$metabolite), pheno.class]]
sigs.m[grepl(x = metabolite,pattern =  "^Q"), class:="quotient"]


# plot
theme_set(theme_light(base_size = 15, base_family = "Helvetica"))
pd <- position_dodge(0.5)
p1 <- ggplot(sigs.m, aes(y=number.sigs,
                         x=reorder(approach, ord)
                         ))+
  geom_violin(alpha=.5,col="lightblue4",draw_quantiles = c(0.25, 0.5, 0.75),fill="grey90", lwd=1) +
  geom_point(position = pd, aes(col=metabolite, pch=class),cex=1.5) +
  geom_line(position = pd, aes(group=metabolite, col=metabolite), alpha=.5,lty="dotted") +
  guides(color=FALSE,line=FALSE) + 
  ylab("Significant associations\nper metabolite") +
  theme(
    axis.text.x = element_text(hjust=1,angle = 45,size=10), 
    axis.title.x = element_blank(),
    panel.border = element_rect(fill=NA),
    legend.position = "bottom"
  )

# show normal and interactive
p1

#' check the unique Gx-Metabolite associations

# unique ID for metab-gx assoc
res.meta[, metab.marker := paste(metabolite, markerID, sep="__")]
uniqueN(res.meta$metab.marker)
length(res.meta$metab.marker)

# venn of significancies based on model type and i2 filter
v1 <- venn3(
  res.meta[numberStudies>=2 & hfdr.sigFEM==TRUE, metab.marker],
  res.meta[numberStudies>=2 & hfdr.sigFEM==TRUE & I2 < .25, metab.marker],
  res.meta[numberStudies>=2 & hfdr.sigREM==TRUE, metab.marker]
)

#' Check some of the associatoins in the FEM-I2<.25 & REM scenario 


# complete distribution over each metabolite
plot.meta <- res.meta

#' How many genes are REM significant per metabolite?

hist(plot.meta[,  sum(hfdr.sigREM), by = metabolite]$V1, breaks =20,main="Significant genes per metabolite - REM") 
hist(plot.meta[,  sum(hfdr.sigFEM), by = metabolite]$V1, breaks =20,main="Significant genes per metabolite - FEM") 
hist(plot.meta[I2<.75,  sum(hfdr.sigFEM), by = metabolite]$V1, breaks =20,main="Significant genes per metabolite - FEM + .75 I2 cutoff") 
hist(plot.meta[I2<.5,  sum(hfdr.sigFEM), by = metabolite]$V1, breaks =20,main="Significant genes per metabolite - FEM + .5 I2 cutoff") 
hist(plot.meta[I2<.25,  sum(hfdr.sigFEM), by = metabolite]$V1, breaks =20,main="Significant genes per metabolite - FEM + .25 I2 cutoff") 

# mean significancies (based on hierarchical FDR)
plot.meta[, mean.sig := sum(hfdr.sigFEM & I2<=.25), by = metabolite]

#' Plot the I2 distribution over each metabolite based on the 50 most FEM significant
#' genes ordered by the total number of significant FEM genes

# identify high heterogeneity assocs
plot.meta[ , i2.cutoff := I2 <= 0.25]

#' This table shows the number of associations that have high heterogeneity compared to their significance. i2.cutoff TRUE means the i2 is >=.75.

# show
plot.meta[,xtabs(~ hfdr.sigFEM + i2.cutoff)]

#' check distribution based on their REM significance

#+ ridge.plot, fig.height = 15, fig.width = 5
ggplot(plot.meta[hfdr.sigREM==TRUE & numberStudies>=2, ],
             # order(hfdr.1FEM, decreasing = F), head(.SD, 50 ), by = metabolite], 
             aes(y = reorder(metabolite, mean.sig),
                 x = I2,
                 fill = ..x..)) +
  stat_density_ridges(
    geom = "density_ridges_gradient",
    calc_ecdf = T,
    panel_scaling = T
  ) +
  scale_fill_gradientn(
    colours = rev(c("#8E063B", "#AB4147", "#C56551", 
                    "#DA8459", "#E99F61", "#F2B669",
                    "#F6C971", "#F4D97B", "#EDE388", 
                    "#E2E6BD"))
  ) + 
  geom_vline(xintercept = .75, lty="dashed",lwd = .2) +
  guides(fill = FALSE) + 
  ylab(label = "") +
  xlab(label = "I2 of REM significant assocs") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

#' Heterogeneity of FEM-significant associations

#+ ridge.plot.fem, fig.height = 15, fig.width = 5
ggplot(plot.meta[hfdr.sigFEM==TRUE & numberStudies>=2, ],
             # order(hfdr.1FEM, decreasing = F), head(.SD, 50 ), by = metabolite], 
             aes(y = reorder(metabolite, mean.sig),
                 x = I2,
                 fill = ..x..)) +
  stat_density_ridges(
    geom = "density_ridges_gradient",
    calc_ecdf = T,
    panel_scaling = T
  ) +
  scale_fill_gradientn(
    colours = rev(c("#8E063B", "#AB4147", "#C56551", 
                    "#DA8459", "#E99F61", "#F2B669",
                    "#F6C971", "#F4D97B", "#EDE388", 
                    "#E2E6BD"))
  ) + 
  geom_vline(xintercept = .75, lty="dashed",lwd = .5) +
  geom_vline(xintercept = .5, lty="dashed",lwd = .5) +
  geom_vline(xintercept = .25, lty="dashed",lwd = .5) +
  guides(fill = FALSE) + 
  ylab(label = "") +
  xlab(label = "I2 of FEM significant assocs") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

# don't need this object any longer
rm(plot.meta)

# add qval for FEM
res.meta$qFEM <- addKorbinianFDR(showplot = FALSE, pvals = res.meta$pFEM)

# Only REM significant
res.meta[
  numberStudies>=2,  .(markerID, metabolite, I2, betaREM, seREM, pREM, qREM, hfdr.sigREM)
  ][
    order(qREM), 
    ][,plot(-log(pREM)~betaREM, 
            col=ifelse(hfdr.sigREM==TRUE, "red","black"),
            pch=20,
            main="REM significancies"
            )]

#' ## Tables of top genes per metabolite

annot <- fread("../../../07_programme/rtools/1807_gx_tools/mappingHT12/results/102_8_remappingHT12_INGENUITY_positionalinfo_ilmnCentric.txt")

# filter annot
matched.annot <- match(res.meta$markerID, annot[ilmn_mapping_hg19_num %in% 1:8, ilmn])
dat <- cbind(res.meta, annot[(matched.annot), .(symbol_INGENUITY,
                                                           location_INGENUITY,  
                                                           types_INGENUITY, 
                                                           biomarkers_INGENUITY,
                                                           drugs_INGENUITY)])

#' ## Annotate number of significant associations

#' Create Columns to indicate 1.) presence in study and effect direction, 2.)
#' per cohort significance and effect direction and 3.) significance overall in
#' meta study using FEM significance, I2<=0.25 & at least two studies present

#' clean and annotate result table
names(dat) <- gsub(names(dat), pattern = " w/ covariates", replacement = "", fixed = TRUE)

# annotate effect directions
dat[, ':='(
  significant.adult.heart.ami.sorb = paste0(
    # test adult
    ifelse(is.na(beta.Adult) | hfdr.sig.Adult == FALSE,
           "n",
           ifelse(sign(beta.Adult)==1, "+", "-")),
    # test heart
    ifelse(is.na(beta.Heart) | hfdr.sig.Heart == FALSE,
           "n",
           ifelse(sign(beta.Heart)==1, "+", "-")),
    # test ami
    ifelse(is.na(`beta.Heart-AMI`) | `hfdr.sig.Heart-AMI` == FALSE,
           "n",
           ifelse(sign(`beta.Heart-AMI`)==1, "+", "-")),
    # test sorb
    ifelse(is.na(beta.Sorb) | hfdr.sig.Sorb == FALSE,
           "n",
           ifelse(sign(beta.Sorb)==1, "+", "-"))
  ),
  analyzed.adult.heart.ami.sorb = paste0(
    # test adult
    ifelse(is.na(beta.Adult),
           "n",
           ifelse(sign(beta.Adult)==1, "+", "-")),
    # test heart
    ifelse(is.na(beta.Heart),
           "n",
           ifelse(sign(beta.Heart)==1, "+", "-")),
    ifelse(is.na(`beta.Heart-AMI`),
           "n",
           ifelse(sign(`beta.Heart-AMI`)==1, "+", "-")),
    # test sorb
    ifelse(is.na(beta.Sorb),
           "n",
           ifelse(sign(beta.Sorb)==1, "+", "-"))
  )
)]

#' ## Annotate metabolite class

#' clean and annotate result table

#' Annotate group as per metabolite annotation
venn2(an$metab.annot$metabolite,unique(dat$metabolite))

m1 <- match(dat$metabolite, an$metab.annot$metabolite)
dat[, `:=`(
  metab.class = an$metab.annot[m1, pheno.class],
  metab.path = an$metab.annot[m1, pathway],
  metab.super.path = an$metab.annot[m1, super.pathway]
)]

# annotate aa/ac mix quotients
dat[metab.class=="", metab.class := "mix"]

# different annotation for quotients
dat[, metab.class.q := metab.class]
dat[grepl(pattern = "^Q", metabolite), metab.class.q := "quotient"]

#' ## Clean INGENUITY annotation

# set dummy instead of NA for better plotting
dat[is.na(location_INGENUITY), location_INGENUITY := "no_data"]
dat[is.na(types_INGENUITY), types_INGENUITY := "no_data"]
dat[location_INGENUITY=="", location_INGENUITY := "no_data"]
dat[types_INGENUITY=="", types_INGENUITY := "no_data"]

#' ## Simple significance indicator

#' It's easier to use one column instead of always selecting two

dat[, significant := ifelse(hfdr.sigREM==TRUE & numberStudies>=2, TRUE, FALSE)]

#' ## Calculate R²

#' PLOS reference used by Arnd: https://doi.org/10.1371/journal.pone.0120758
#' calculate using beta and se(beta)
#' Calculate se as beta/t-statistic

#' Variant one
#' beta^2 / ( beta^2 + N * se(beta)^2 )

dat[, r.squared :=  betaREM^2 / ( betaREM^2 + totalN * seREM^2 )]

#' Variant two -> seems to be a problem? I don't know
#' seen here: http://www.bios.unc.edu/research/genomic_software/Matrix_eQTL/faq.html
#' Discussed here: https://www.quora.com/Whats-the-difference-between-t-test-and-correlation
#' r^2 = tstat^2 / dfFull + tstat^2
#' I do not have the full df nor the t-statistic of the meta GWAS, therefore I'll use the other formula!

#' ## Save

# save for later
# nf <- newest_file(look_for = "meta_resultsAnnotated", subfolder = "res", print_full = TRUE)
# dat <- fread(nf)
fwrite(dat, file = "res/STEP4_meta_resultsAnnotated.csv", quote = TRUE)
# save_csv_carl(file = dat, file_name = "meta_resultsAnnotated", subfolder = "res", )

#' # Get New INGENUITY Gx Mapping

#' The new mapping for the Gx probes is available, use it for probe-gene annotation

# load the mapping file
gx.map <- fread(paste0(bp, "/07_programme/rtools/1807_gx_tools/mappingHT12/results/102_8_remappingHT12_INGENUITY_positionalinfo_ilmnCentric.txt"))

#' There are several QC columns I need to consider before matching these to my result. There are duplications and badly mapped probes.
#' I need to explore these a little longer.

# check names
names(gx.map)

# quality col
table(gx.map$ilmn_mapping_hg19)

# corresponds to 
table(gx.map$ilmn_mapping_hg19_num)

# special case for y-chromosome
table(gx.map$ychr_2nd_uniquemap)

#' # Prepare Data for IPA

ipa <- dcast.data.table(
  data = dat[numberStudies>=2, .(markerID, metabolite, betaREM, pREM, qREM)], 
  formula = markerID~metabolite,
  value.var = c("betaREM", "pREM", "qREM"), 
  sep = "__"
)

# paste together names for IPA
m1 <- grepl(pattern = paste("__", an$metab.annot$metabolite[1:20], "$", collapse = "|", sep = ""), x = names(ipa))
m2 <- grepl(pattern = paste("__", an$metab.annot$metabolite[21:40], "$", collapse = "|", sep = ""), x = names(ipa))
m3 <- grepl(pattern = paste("__", an$metab.annot$metabolite[41:60], "$", collapse = "|", sep = ""), x = names(ipa))
m4 <- grepl(pattern = paste("__", an$metab.annot$metabolite[61:80], "$", collapse = "|", sep = ""), x = names(ipa))
m5 <- grepl(pattern = paste("__", an$metab.annot$metabolite[81:97], "$", collapse = "|", sep = ""), x = names(ipa))

# get the 97 metabolites in five sets with max 20 metabolites
set1 <- ipa[, .SD, .SDcols= c("markerID", names(ipa)[m1])]
set2 <- ipa[, .SD, .SDcols= c("markerID", names(ipa)[m2])]
set3 <- ipa[, .SD, .SDcols= c("markerID", names(ipa)[m3])]
set4 <- ipa[, .SD, .SDcols= c("markerID", names(ipa)[m4])]
set5 <- ipa[, .SD, .SDcols= c("markerID", names(ipa)[m5])]

# save them
write.delim(set1, "obj/IPA/set1.txt", sep="\t")
write.delim(set2, "obj/IPA/set2.txt", sep="\t")
write.delim(set3, "obj/IPA/set3.txt", sep="\t")
write.delim(set4, "obj/IPA/set4.txt", sep="\t")
write.delim(set5, "obj/IPA/set5.txt", sep="\t")

#' # Forest Plots

# what are the columns called
study <- c("FEM", "REM", ".Sorb", ".Heart", ".Heart-AMI", ".Adult") 
study.pattern <- paste(study, collapse = "|")

# get all the columns with test statistics I want to melt
use.cols <- grepl(names(dat), pattern = study.pattern)

# get the columns I do not want to melt
hfdr.sig.cols <- !grepl(names(dat), pattern = "hfdr.sig")

# these are all the columns I want to melt
melt.vars <- names(dat)[use.cols & hfdr.sig.cols]

# melt the data for better reordering
res.meta.melt <- melt(dat,
                      id.vars = c("markerID", "metabolite", "symbol_INGENUITY"),
                      measure.vars = c(melt.vars, "I2", "totalN"),
                      variable.name = "statistic",
                      value.name = "stat.val")

# get a function to extract cohort ID from the data
str_extract <- function(string, extract){
  str.extract <- substring(
    string,
    regexpr(extract, string),
    regexpr(extract, string) + attr(regexpr(extract, string), "match.length") - 1
  )
  return(str.extract)
}

# example:
# str_extract(string = c("skdfSorbskfjsd", "ksdfjdsfAdultklsdjfsd"),extract = "Adult|Sorb")

# get a study identifier
setDTthreads(4)
res.meta.melt[, statistic.clean := gsub(pattern = ".Adult|.Sorb|.Heart|.Heart-AMI|FEM|REM", replacement = "", x = statistic)]
# res.meta.melt[, statistic.clean := gsub(pattern = " w/ covariates", replacement = "", x = statistic.clean, fixed=TRUE)]
res.meta.melt[,.N,by=statistic.clean]

# clean the study indicator
res.meta.melt[, study := str_extract(string = statistic, extract = "Adult|Sorb|Heart|Heart-AMI|FEM|REM")]
res.meta.melt[,.N,by=study]
res.meta.melt[, study.clean := gsub(pattern = " w/ covariates", replacement = "", x = study, fixed = TRUE)]
res.meta.melt[,.N,by=study.clean]

# sort out I2
res.meta.melt[statistic.clean=="I2", study.clean:="REM"]
res.meta.melt[statistic.clean=="totalN", study.clean:="REM"]

# re-cast the data 
res.meta.cast <- dcast(
  data = res.meta.melt,
  formula = markerID + metabolite + study.clean + symbol_INGENUITY ~ statistic.clean, 
  value.var = "stat.val",
  sep = "."
)

# create col indicating study or summary
res.meta.cast[, summary := ifelse(study.clean %in% c("REM", "FEM"), TRUE, FALSE)]

# save cleaned results used for plotting
# nf <- newest_file(look_for = "a1_b3_sorb_metaAnalysisResultsClean", subfolder = "res", print_full = TRUE)
# res.meta.cast <- fread(nf)
save_csv_carl(file = res.meta.cast, file_name = "a1_b3_sorb_metaAnalysisResultsClean", subfolder = "res")

# draw example forest plot
library(forestplot)
forest_plot <- function(
  plotProbe,
  plotMetabolite,
  results
) {
  example.plot <- results[markerID == (plotProbe) & metabolite == (plotMetabolite), ][order(summary), ]
  i2 <- example.plot[study.clean=="REM", I2]
  total.n <- example.plot[study.clean=="REM", totalN]
  example.plot[, forestplot(labeltext = paste0(study.clean, ", p=", signif(p,3)),
                            mean = beta,
                            lower = ci.l,
                            upper = ci.u, 
                            is.summary = summary,
                            title = paste0(unique(example.plot$symbol_INGENUITY), 
                                           " (", unique(example.plot$markerID), ")",
                                           " ~ ",
                                           unique(example.plot$metabolite),
                                           ", I2 = ", 
                                           signif(i2,2),
                                           ", N = ", total.n))]
}

#' Get some example associations
plot.these <- dat[hfdr.sigFEM == T & I2<=0.25 & numberStudies >=2, ][
  order(abs(betaFEM), decreasing = T), .(markerID, metabolite, betaFEM)][
    1:10, ]

# plot some examples
plot.these[, index := 1:.N]
pdf("plt/190912_exampleForestPlots.pdf", width = 7, height = 7)
for(i in plot.these$index){
  
  my.probe <- plot.these[index == (i), markerID]
  my.metab <- plot.these[index == (i), metabolite]
  
  forest_plot(plotProbe = my.probe, plotMetabolite = my.metab, results = res.meta.cast)
  
}
dev.off()

#' # Finalize script
devtools::session_info()
