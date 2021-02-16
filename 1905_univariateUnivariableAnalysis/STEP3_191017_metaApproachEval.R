#' ---
#' title: "Comparison of three meta-analysis approaches for Gene Expression ~ Metabolite association"
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
#'          smooth_scroll: FALSE
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
  computer <- "forostar"
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
  "fdrtool",
  "here",
  "ggplot2",
  "scales",
  "plotly",
  "magrittr",
  "ggridges",
  "readxl"
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
source(here("../functions/compare_limma_matrixeqtl.R"))
source(here("../functions/all_annotation_tables.R"))
annot <- all_annotation_tables(mount.point = "/net/ifs2/san_projekte/projekte/genstat/")
option_setup()

#+ include=F
#==========#
#   DONE   #
#==========#

#' # Load Data from previous Step

#' ## Load metabolite data

# these are the metabolites adjusted for the gx covariates
metab.dir <- paste0(bp, "/02_projekte/1703_ge_metab_a1_b3_sorbs/170824_MetabAnnotTab/1906_adjustMetabsForGx_NoBMI_NoDiabetes/results/")
metab.loc <- newest_file(look_for = "a1_b3_sorb_metabolites.csv",
                         directory = metab.dir,
                         print_full = TRUE)
dat.metab <- fread(metab.loc)

#' ## Load Meta Analysis Assoc Results (based on LIMMA Summary Stats)

#' Meta-Analysis with additional added covariates

nf <- newest_file(
  look_for = "a1_b3_sorb_metaAnalysisWithCovarsLimmaResults",
  subfolder = "res",
  print_full = T)
res.c.meta <- fread(nf)

#' Check dimensions of the result files

dim(res.c.meta)

#' Add FDR for plotting of approach comparison

#' Use Strimmer fdr for this
res.c.meta$qREM <- addKorbinianFDR(showplot = FALSE, pvals = res.c.meta$pREM)

#===================================================================#

#' ## Summaries of different Meta-Analysis approaches

#' compile FDR info about each of the approaches in one table

# create table by filtering for 5% and 20% FDR
sigs.meta.all <- data.table(
  metabolite = unique(res.c.meta$metabolite),
  c.fdr20 = res.c.meta[, sum(qREM<=.2), by=metabolite]$V1,
  c.fdr05 = res.c.meta[, sum(qREM<=.05), by=metabolite]$V1
)

# display results
kable(sigs.meta.all)

# draw plot
theme_set(theme_minimal(base_size = 13, base_family = "Helvetica"))

#' Try to visualize the table
sma <- melt(sigs.meta.all, id.vars = "metabolite")
sma[,`:=`(
  approach=ifelse(grepl(x = variable, pattern = "nc.", fixed = TRUE), "unadjustPrepro_adjustLimma",
                  ifelse(grepl(x = variable, pattern = "^c."), "adjustPrepro_adjustLimma",
                         "adjustPrepro_unadjustLimma")),
  fdr=ifelse(grepl(x = variable, pattern = "20", fixed = TRUE), "fdr20","fdr05")
)]

# get something to oder the approaches by
sma[,ord:=max(value),by=approach]

# plotting
pd <- position_dodge(0.5)
p1 <- ggplot(sma, aes(y=value,
                      x=reorder(approach, ord)) #, 
                      # label=paste(metabolite,  value, sep="\n"), 
                      # pch = grepl("^Q", metabolite))
             )+
  facet_grid(~fdr) +
  geom_violin(alpha=.5,col="lightblue4",draw_quantiles = c(0.25, 0.5, 0.75),fill="grey90") +
  geom_point(position = pd, aes(col=metabolite),cex=1.2) +
  geom_line(position = pd, aes(group=metabolite, col=metabolite), alpha=.5,lty="dotted") +
  guides(color=FALSE,line=FALSE) + 
  ylab("Significant associations\nper metabolite") +
  theme(
    axis.text.x = element_text(hjust=1,angle = 45,size=10), 
    axis.title.x = element_blank(),
    panel.border = element_rect(fill=NA)
  )

# p1
ggplotly(p1)

#' ## Check I2 of the different analyses

# uc = use columns
uc <- c("markerID", "metabolite", "I2", "qREM", "numberStudies", "totalN", 
        grep(x = names(res.c.meta),
             pattern = "REM",
             fixed = TRUE, 
             value = TRUE))

# concatinate
res.c.meta[,type:="adjustPrePro_adjustLimma"]
res.meta.all <- rbindlist(
  list(
    res.c.meta
  ), use.names = TRUE, fill = TRUE)

# summarise
i2.meta <- res.meta.all[hfdr.sigREM==T & numberStudies>1,
                        .(mean.i2 = mean(I2, na.rm=TRUE),
                          median.i2 = median(I2, na.rm=TRUE),
                          perc025 = quantile(I2,.25, na.rm=TRUE),
                          perc050 = quantile(I2,.50, na.rm=TRUE),
                          perc075 = quantile(I2,.75, na.rm=TRUE),
                          min = min(I2, na.rm=TRUE),
                          max = max(I2, na.rm=TRUE),
                          sd.i2 = sd(I2, na.rm=TRUE)), by = .(metabolite, type)]

save_csv_carl(file = i2.meta, file_name = "meta_i2_summary", subfolder = "obj")

#' Plot description: Red square + lines are median + IQR, grey triangle is mean, + is minimum valie and x is maximum value

# plot
ggplot(i2.meta) +
  coord_flip() +
  facet_grid(~type) + 
  geom_point(aes(x=reorder(metabolite, mean.i2),y=mean.i2), pch = 2, col="grey30") +
  geom_point(aes(x=reorder(metabolite, mean.i2),y=min), pch = 3, col = "grey60") +
  geom_point(aes(x=reorder(metabolite, mean.i2),y=max), pch = 4, col = "grey60") +
  geom_pointrange(aes(ymin=perc025,
                      y=perc050,
                      ymax=perc075, 
                      x=reorder(metabolite, mean.i2)),
                  lwd=.3, col="indianred4", pch=15) +
  ylab("I2") +
  theme(
    axis.title.y = element_blank()
  )

#' ## Check Betas/Standard-Errors

#' summarise beta/se over each metabolite

# annotate number of significant assocs per metabolite
res.meta.all[, num.sigs := sum(hfdr.sigREM==T & numberStudies>1), by = .(metabolite, type)]

bse <- res.meta.all[hfdr.sigREM==T & numberStudies>1, 
                    .(beta.perc25 = quantile(betaREM,.25, na.rm=TRUE),
                      beta.perc50 = quantile(betaREM,.50, na.rm=TRUE),
                      beta.perc75 = quantile(betaREM,.75, na.rm=TRUE),
                      beta.min = min(betaREM, na.rm=TRUE),
                      beta.max = max(betaREM, na.rm=TRUE),
                      beta.mean = mean(betaREM, na.rm=TRUE),
                      se.perc25 = quantile(seREM,.25, na.rm=TRUE),
                      se.perc50 = quantile(seREM,.50, na.rm=TRUE),
                      se.perc75 = quantile(seREM,.75, na.rm=TRUE),
                      se.min = min(seREM, na.rm=TRUE),
                      se.max = max(seREM, na.rm=TRUE),
                      se.mean = mean(seREM, na.rm=TRUE),
                      num.sigs = mean(num.sigs)
                      ), by=.(metabolite,type)]

# check betas
ggplot(bse, aes(ymin = beta.perc25, 
                y = beta.perc50, 
                ymax = beta.perc75,
                x = reorder(metabolite, beta.perc50))) +
  coord_flip() + 
  scale_y_continuous() + 
  facet_wrap(~type) +
  geom_hline(aes(yintercept=0), col="grey4") +
  geom_pointrange(fatten=1, size=1, col="indianred4") +
  geom_point(aes(y=beta.mean), shape=2, col="grey15") +
  geom_point(aes(y=beta.min), shape=3, col="grey60") +
  geom_point(aes(y=beta.max), shape=4, col="grey60") +
  ylab("Beta") +
  geom_text(aes(label=num.sigs, y = .8)) + 
  theme(
    axis.title.y = element_blank()
    )

#' ## Look at standard errors
 
ggplot(bse, aes(ymin = se.perc25, 
                y = se.perc50, 
                ymax = se.perc75,
                x = reorder(metabolite, se.perc50))) +
  coord_flip() + 
  facet_grid(~type) +
  scale_y_continuous() + 
  geom_hline(aes(yintercept=0), col="grey4") +
  geom_pointrange(fatten=1, size=1, col="indianred4") +
  geom_point(aes(y=se.mean), shape=2, col="grey15") +
  geom_point(aes(y=se.min), shape=3, col="grey60") +
  geom_point(aes(y=se.max), shape=4, col="grey60") +
  ylab("Standard Error") +
  theme(
    axis.title.y = element_blank()
  )

#' look at a few strong associations and how they relate over the different approaches

# get all significant assocs for each type
all.sigs <- res.meta.all[hfdr.sigREM==TRUE & numberStudies>1, ][order(hfdr.1REM, decreasing = FALSE), .(metabolite, type, markerID, seREM, hfdr.1REM, betaREM)]

# asc =all sigs cast
asc <- dcast.data.table(data = all.sigs, 
                        formula = metabolite+markerID~type, value.var = c("hfdr.1REM", "betaREM"), sep = "__")

# this is weird: adding T + T = 2, T + !T = 1 but !T + T = FALSE,
# however (!T) + T = 1
# T + T
# (!T) + T
# T + !T
# is.na(1) + is.na(1)
# !is.na(1) + is.na(1)
# (!is.na(1)) + is.na(1)

#' Look at eta1 of each analysis approach
res.eta1 <- res.meta.all[, .(eta1 = 1 - fdrtool::pval.estimate.eta0(pREM, diagnostic.plot = FALSE)), by = .(metabolite, type)]
mean.eta1 <- res.eta1[, .(eta1 = mean(eta1)), by = metabolite]
mean.eta1[, type:="mean"]

#' Plot eta 1 by mean over each type and individually

ggplot()+
  geom_point(data=mean.eta1, aes(y = reorder(metabolite, eta1), x=eta1), col="grey5") + 
  geom_point(data=res.eta1,
             aes(y= metabolite, col=type, x=eta1), 
             alpha=.75) +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    axis.title.y = element_blank()
  )

ggplot(data=res.eta1,
       aes(fill=type, x=eta1))+
  geom_density(alpha=.5) +
  geom_rug(alpha=0.5, aes(col=type)) +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    axis.title.y = element_blank()
  )

#' ## Plot metabolite-wise eta1 

# use a better theme for the plot
theme_set(theme_light(base_size = 13, base_family = "Helvetica"))

#' summarise results over metabolites
res.eta1[,overall.mean:=mean(eta1), by=type]

#' This plot shows the distribution of eta1 for each metabolite in the three approaches
pd <- position_dodge(.5)

# plot metabolite eta1
eta.plot <- ggplot(res.eta1, aes(reorder(type, overall.mean), eta1, color = type)) +
  coord_flip() +
  scale_y_continuous(limits = c(0, .6), expand = c(0.005, 0.005), breaks = pretty_breaks()) +
  
  # metabolites and connecting lines
  stat_ydensity(col="grey60",fill="grey90", alpha=.1, draw_quantiles = c(.25,.5,.75)) +
  scale_color_manual(values = c("#001889", "#91008D", "#D24E71")) +
  geom_line(position = pd, aes(group=metabolite), alpha=.5,lty="dotted", col="grey50") +
  geom_point(position = pd, aes(group=metabolite), size = 2, alpha = 0.7) +
  
  # over all mean and type-mean
  geom_point(aes(type, overall.mean), size = 5) +
  geom_hline(aes(yintercept = mean(eta1)), color = "grey70", size = 0.6) +
  geom_segment(aes(x = type,
                   xend = type,
                   y = mean(eta1),
                   yend = overall.mean),
               size = 0.8) +
  labs(x = NULL, y = "Eta1") +
  theme(legend.position = "none",
        axis.title = element_text(size = 12),
        axis.text.x = element_text(size = 10),
        panel.grid = element_blank())

# interactive version 
ggplotly(eta.plot)

save_csv_carl(file = res.eta1, file_name = "meta_eta1_summary", subfolder = "obj")

# show the metabolites with the highest eta1
kable(res.eta1[order(eta1, decreasing = T), .(type, metabolite)][1:15,])

res.meta <- res.meta.all[type=="adjustPrePro_adjustLimma",
                         .(markerID, metabolite, numberStudies, num.sigs, totalN, pFisher, Z, pZ, 
                           betaFEM, seFEM, pFEM, I2, betaREM, seREM, pREM, ci.lFEM, 
                           ci.uFEM, hfdr.1FEM, hfdr.2FEM, hfdr.sigFEM, ci.lREM, ci.uREM, 
                           hfdr.1REM, hfdr.2REM, hfdr.sigREM, qREM, 
                           `beta.Adult w/ covariates`, `se.Adult w/ covariates`, `p.Adult w/ covariates`, 
                           `n.Adult w/ covariates`, `beta.Heart w/ covariates`, `se.Heart w/ covariates`, 
                           `p.Heart w/ covariates`, `n.Heart w/ covariates`, `beta.Heart-AMI w/ covariates`, 
                           `se.Heart-AMI w/ covariates`, `p.Heart-AMI w/ covariates`, `n.Heart-AMI w/ covariates`, 
                           `beta.Sorb w/ covariates`, `se.Sorb w/ covariates`, `p.Sorb w/ covariates`, 
                           `n.Sorb w/ covariates`, `ci.l.Adult w/ covariates`, `ci.u.Adult w/ covariates`, 
                           `hfdr.1.Adult w/ covariates`, `hfdr.2.Adult w/ covariates`, `hfdr.sig.Adult w/ covariates`, 
                           `ci.l.Heart w/ covariates`, `ci.u.Heart w/ covariates`, `hfdr.1.Heart w/ covariates`, 
                           `hfdr.2.Heart w/ covariates`, `hfdr.sig.Heart w/ covariates`, 
                           `ci.l.Heart-AMI w/ covariates`, `ci.u.Heart-AMI w/ covariates`, 
                           `hfdr.1.Heart-AMI w/ covariates`, `hfdr.2.Heart-AMI w/ covariates`, 
                           `hfdr.sig.Heart-AMI w/ covariates`, `ci.l.Sorb w/ covariates`, 
                           `ci.u.Sorb w/ covariates`, `hfdr.1.Sorb w/ covariates`, `hfdr.2.Sorb w/ covariates`, 
                           `hfdr.sig.Sorb w/ covariates`)]

# save for STEP4
save( 
  i2.meta,
  res.eta1,
  bse,
  res.meta,
  file = "obj/190729_STEP3.RData")

devtools::session_info()
