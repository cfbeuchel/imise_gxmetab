#' ---
#' title: "Vizualize the Mediation"
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
  "igraph",
  "magrittr",
  "ggplot2",
  # "treemap",
  "Vennerable",
  "ggExtra",
  "ggrepel",
  "cowplot",
  "grid",
  "colorspace",
  "parallel",
  "plyr"
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
# source(here::here("../functions/all_annotation_tables.R"))
source(here::here("../functions/option_setup.R"))
source(here::here("../functions/theme_gxMetab.R"))
an <- all_annotation_tables(mount.point = bp)
theme_set(theme_light(base_size = 15, base_family = "Helvetica"))
setDTthreads(8)
option_setup()

#+ include=F
#==========#
#   DONE   #
#==========#

# Load data ----
#' # Load data

# results from the main meta-analysis
nf <- newest_file(look_for = "STEP4_meta_resultsAnnotated$",subfolder = "res", print_full = TRUE)
res.init.meta <- fread(nf)

# load the joint mediation results
res.joint <- fread(newest_file("mediationStatisticsAnnotatedJoint", "res",print_full = TRUE))

# annotated sobel results
res.sobel <- fread(newest_file("mediationStatisticsAnnotated$","res",print_full = T))

# metabolite names
load("/net/ifs2/san_projekte/projekte/genstat/02_projekte/1703_ge_metab_a1_b3_sorbs/1805_gxMetaboliteAssociation/obj/200306_metaboliteNames.RData")

#=================================================================================#
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#=================================================================================#

#' # Colorscheme Metabolites/Genes
# Colorscheme Metabolites/Genes ----

pal.metabs <- sequential_hcl(palette = "OrYel",n=100)
pal.genes <- sequential_hcl(palette = "Teal",n=100)
pal.pairs <- sequential_hcl(palette = "Greens2",n=100)

#=================================================================================#
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#=================================================================================#

# Checking Plots ----
#' # Checking Plots 

# Gx->M vs. M->Gx ----

# match both
m1 <- match(
  res.joint[type=="gx~metab",paste0(metabolite,gx.probe)],
  res.joint[type=="metab~gx",paste0(metabolite,gx.probe)])

# plot betas
plot(
  res.joint[type=="gx~metab",betaREM],
  res.joint[type=="metab~gx",][m1,betaREM],
  col=alpha("grey35",0.5),
  pch=20,
  xlab="gx~metab",
  ylab="metab~gx"
)
abline(0,1,col="indianred4")

# plot se
plot(
  res.joint[type=="gx~metab",seREM],
  res.joint[type=="metab~gx",][m1,seREM],
  col=alpha("grey35",0.5),
  pch=20,
  xlab="gx~metab",
  ylab="metab~gx"
)
abline(0,1,col="indianred4")

# \tau vs \beta

m2 <- match(
  res.joint[type=="pheno~gx(+metab)" & pheno == "log.bmi",paste0(metabolite, gx.probe)],
  res.joint[type=="pheno~metab(+gx)" & pheno == "log.bmi",paste0(metabolite, gx.probe)]
)

plot(
  res.joint[type=="pheno~gx(+metab)" & pheno == "log.bmi", betaREM],
  res.joint[type=="pheno~metab(+gx)" & pheno == "log.bmi",][m2, betaREM],
  col=alpha("grey35",0.5),
  pch=20,
  xlab="gx effect",
  ylab="metab effect"
)

# R2 mediation yes/no ----
#' ## R2 mediation yes/no

#' Plot the mediation effect alpha*beta

p.alpha.beta <- ggplot(res.joint[is.mediation=="yes"],aes(x=beta.alpha.beta, fill=fdr.sig))+
  geom_histogram(position="stack",col="grey35",bins=60)+ # ,fill="grey80"
  geom_vline(xintercept = 0,lwd=0.6,col="indianred4",lty="dashed") +
  facet_wrap(type~pheno,scales="free")+
  scale_x_continuous(limits = c(-0.05,0.05)) +
  theme_gxMetab() + 
  guides(fill=guide_legend(title = "Significance at hierarchical FDR=5%")) + 
  labs(title = expression("Distribution of mediation effects calculated as " * beta[alpha] * beta[beta]),
       y=NULL,
       x=expression(beta[alpha] - beta[beta])) + 
  theme(legend.position = "bottom", 
        legend.justification = c(0,0)) 

print(p.alpha.beta)

# tiff(dpx("alphaTimesBetaHistogram.tiff","plt/"),res=300,unit="in",height=8,width=8)
# print(p.alpha.beta)
# dev.off()

#' Plot different mediation betas

# tiff(dpx("tauTauAdjVsAlphaBetaEstimate.tiff","plt/"),res=300,unit="in",height=8,width=8)
ggplot(res.joint[is.mediation=="yes" & fdr.sig == TRUE & !is.na(beta.tau.tauAdj),],
       aes(x=beta.tau.tauAdj,
           y=beta.alpha.beta))+
  geom_point(col="grey35",pch=21)+ # ,fill="grey80"
  geom_abline(slope = 1,intercept = 0,lwd=1.2,col="indianred4",lty="dashed") +
  facet_wrap(type~pheno,scales="free")+
  theme_gxMetab()
# dev.off()

# Joint vs. raw effect sizes ----
#' ## Joint vs. raw effect sizes

# distribution of all combinations

r2.mediation <- res.joint[type%in%c("pheno~gx(+metab)","pheno~metab(+gx)"), 
                          .(r2.arnd = r.squared.REM, # this is not the r2 of the mediation, but of the direct effect! misleading!
                            # beta.mediation = beta.tau.tauAdj, # beta of the sobel test 
                            beta.mediation = beta.alpha.beta, # beta of the sobel test 
                            ci.95.lower,
                            ci.95.upper,
                            alpha,
                            seAlpha,
                            beta,
                            seBeta,
                            tau,
                            tauAdj,
                            p.value.prod,
                            fdr.local = fdr1,
                            fdr.global = fdr2,
                            mediation.sig = fdr.sig
                          ),
                          by=.(type,pheno,gx.probe,gene = symbol_INGENUITY,metabolite)]

# NEW 200608 remove non-mediation tests from r2.mediation: It still contains ALL assocs, but not all assocs were taken fo mediation

#' Add the single assoc r2 to each adjusted pair
m.gx <- match(r2.mediation[type=="pheno~gx(+metab)",
                           paste0(pheno,gx.probe)],
              res.joint[type=="pheno~gx",paste0(pheno,gx.probe)])
m.metab <- match(r2.mediation[type=="pheno~metab(+gx)",
                              paste0(pheno,metabolite)],
                 res.joint[type=="pheno~metab",paste0(pheno,metabolite)])

# add the marginal explained variance of the main effect (not the mediator!)
r2.mediation[type=="pheno~gx(+metab)", `:=`(
  r2.gx.raw = res.joint[type=="pheno~gx",][m.gx, r.squared.REM], # not marginal, but "raw"
  beta.gx.raw = res.joint[type=="pheno~gx",][m.gx, betaREM],
  sig.gx.raw = res.joint[type=="pheno~gx",][m.gx, hfdr.sigREM]
)]
r2.mediation[type=="pheno~metab(+gx)", `:=`(
  r2.metab.raw = res.joint[type=="pheno~metab",][m.metab, r.squared.REM],
  beta.metab.raw = res.joint[type=="pheno~metab",][m.metab, betaREM],
  sig.metab.raw = res.joint[type=="pheno~metab",][m.metab, hfdr.sigREM]
)]

# add mediator raw significance 

m.mediator.gx <- match(r2.mediation[,paste0(pheno,metabolite)],
                       res.joint[type=="pheno~metab",paste0(pheno,metabolite)])
m.mediator.metab <- match(r2.mediation[,paste0(pheno,gx.probe)],
                          res.joint[type=="pheno~gx",paste0(pheno,gx.probe)])

r2.mediation[, sig.metab.raw := res.joint[type=="pheno~metab",][m.mediator.gx, hfdr.sigREM]]
r2.mediation[, r2.metab.raw := res.joint[type=="pheno~metab",][m.mediator.gx, r.squared.REM]]
r2.mediation[, beta.metab.raw := res.joint[type=="pheno~metab",][m.mediator.gx, betaREM]]
r2.mediation[, sig.gx.raw    := res.joint[type=="pheno~gx",][m.mediator.metab, hfdr.sigREM]]
r2.mediation[, r2.gx.raw    := res.joint[type=="pheno~gx",][m.mediator.metab, r.squared.REM]]
r2.mediation[, beta.gx.raw    := res.joint[type=="pheno~gx",][m.mediator.metab, betaREM]]

# classify mediation by their marginal significance
# r2.mediation[, both.marginal.sig  := ifelse(marginal.gx.sig & marginal.metab.sig, TRUE, FALSE)]
# r2.mediation$both.marginal.sig %>% table
r2.mediation[, raw.sig.group  := ifelse(sig.gx.raw & sig.metab.raw, "both", 
                                        ifelse((type=="pheno~gx(+metab)" & sig.gx.raw==TRUE) | 
                                                 (type=="pheno~metab(+gx)" & sig.metab.raw==TRUE), "main.effect",
                                               ifelse((type=="pheno~gx(+metab)" & sig.metab.raw==TRUE) | 
                                                        (type=="pheno~metab(+gx)" & sig.gx.raw==TRUE), "mediating.effect", 
                                                      "CHECK!!!")))]

# checks
r2.mediation$raw.sig.group %>% table
r2.mediation$raw.sig.group %>% is.na %>% any
any(r2.mediation$raw.sig.group == "CHECK!!!", na.rm = T)

# check - in all cases, no row should be FALSE for both marginal associations! (prerequisite of analysis)
r2.mediation[, all(sig.metab.raw | sig.gx.raw, na.rm = T)] # -> all TRUE

# Save Mediation Info with matching marginal effects 
fwrite(r2.mediation, dpx("mediationStatisticsMatchingMarginalInfo",folder = "res/"),quote = TRUE, sep="\t")
# r2.mediation <- fread(newest_file("mediationStatisticsMatchingMarginalInfo", "res/",print_full = T))

#=================================================================================#
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#=================================================================================#

# Plot: Venn Sobel unique features ----
#' ## Plot: Venn Sobel unique features 

# get the unique genes/metabolites per phenotype/mediation direction

# metab
w1 <- compute.Venn(Venn(Sets=list(
  "Metabolites as Exposure" = res.sobel[fdr.sig==TRUE & y=="log.bmi" & x.type=="metabolite",unique(metabolite)],
  "Metabolites as Mediator" = res.sobel[fdr.sig==TRUE & y=="log.bmi" & x.type=="gx.probe",unique(metabolite)]
)), type = "circles")

gp1 <- VennThemes(w1)
names(gp1)
gp1[["Face"]][["11"]]$fill <-  pal.metabs[50]
gp1[["Face"]][["01"]]$fill <-  pal.metabs[1]
gp1[["Face"]][["10"]]$fill <-  pal.metabs[100]
gp1[["Set"]][["Set1"]][["col"]] <- pal.pairs[50]
gp1[["Set"]][["Set2"]][["col"]] <- pal.pairs[1]
gp1[["SetText"]][["Set1"]][["col"]] <- pal.pairs[50]
gp1[["SetText"]][["Set2"]][["col"]] <- pal.pairs[1]

tiff(dpx("sobelVennBMIMetabolites.tiff","plt/"),res=300,unit="in",height=7,width=7.9)
plot(w1, gp = gp1)
dev.off()

# gx
w2 <- compute.Venn(Venn(Sets=list(
  "Genes as Exposure" = res.sobel[fdr.sig==TRUE & y=="log.bmi" & x.type=="gx.probe" & symbol_INGENUITY != "",symbol_INGENUITY],
  "Genes as Mediator" = res.sobel[fdr.sig==TRUE & y=="log.bmi" & x.type=="metabolite" & symbol_INGENUITY != "",symbol_INGENUITY]
)), type = "circles")

gp2 <- VennThemes(w2)
names(gp2)
gp2[["Face"]][["11"]]$fill <-  pal.genes[50]
gp2[["Face"]][["01"]]$fill <-  pal.genes[1]
gp2[["Face"]][["10"]]$fill <-  pal.genes[100]
gp2[["Set"]][["Set1"]][["col"]] <- pal.pairs[50]
gp2[["Set"]][["Set2"]][["col"]] <- pal.pairs[1]
gp2[["SetText"]][["Set1"]][["col"]] <- pal.pairs[50]
gp2[["SetText"]][["Set2"]][["col"]] <- pal.pairs[1]

tiff(dpx("sobelVennBMIGenes.tiff","plt/"),res=300,unit="in",height=7,width=7.9)
plot(w2, gp = gp2)
dev.off()

# how many metabolites are in which 

tiff(dpx("sobelVennMetabolites.tiff","plt/"),res=300,unit="in",height=7,width=7.9)
v1 <- venn4(
  res.sobel[fdr.sig==TRUE & y=="log.bmi" & x.type=="metabolite",x.gene],
  res.sobel[fdr.sig==TRUE & y=="log.bmi" & x.type=="gx.probe",m.gene],
  res.sobel[fdr.sig==TRUE & y=="diabetes.status.tri" & x.type=="metabolite",x.gene],
  res.sobel[fdr.sig==TRUE & y=="diabetes.status.tri" & x.type=="gx.probe",m.gene],
  mylabels = c("BMI mediated\nmetabolites", "BMI mediating\nmetabolites",
               "T2D mediated\nmetabolites", "T2D mediating\nmetabolites"),
  mytitle = "Overview of unique significant metabolite features\nin 4 mediation scenarios"
)
dev.off()

# gene expression

tiff(dpx("sobelVennGenes.tiff","plt/"),res=300,unit="in",height=7,width=7.9)
v2 <- venn4(
  res.sobel[fdr.sig==TRUE & y=="log.bmi" & x.type=="metabolite",m.gene],
  res.sobel[fdr.sig==TRUE & y=="log.bmi" & x.type=="gx.probe",x.gene],
  res.sobel[fdr.sig==TRUE & y=="diabetes.status.tri" & x.type=="metabolite",m.gene],
  res.sobel[fdr.sig==TRUE & y=="diabetes.status.tri" & x.type=="gx.probe",x.gene],
  mylabels = c("BMI mediating\ngenes", "BMI mediated\ngenes",
               "T2D mediating\ngenes", "T2D mediated\ngenes"),
  mytitle = "Overview of unique significant genes\nin 4 mediation scenarios"
)
dev.off()

#=================================================================================#
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#=================================================================================#

#' ## Plot: Barplot of mediated metabolites
# Plot: Barplot of mediated metabolites ----

tmp1 <- res.init.meta[, .(analysis = "meta.all", metabolites = unique(metabolite))]
tmp2 <- res.init.meta[significant==TRUE, .(analysis = "meta.sig", metabolites = unique(metabolite))]

tmp3 <- res.joint[pheno=="log.bmi" & type %in% c("pheno~metab") & hfdr.sigREM == TRUE, 
                  .(analysis = "bmi.assoc", 
                    metabolites = unique(metabolite)
                  )
                  ]

tmp4 <- res.joint[pheno=="log.bmi" & type %in% c("pheno~gx(+metab)","pheno~metab(+gx)") & fdr.sig == TRUE, 
                  .(analysis = "mediation.sig", 
                    metabolites = unique(metabolite)
                  )
                  ]
tmp5 <- res.joint[pheno=="log.bmi" & type %in% c("pheno~metab(+gx)") & fdr.sig == TRUE, 
                  .(analysis = "mediated.by.gx", 
                    metabolites = unique(metabolite)
                  )
                  ]
tmp6 <- res.joint[pheno=="log.bmi" & type %in% c("pheno~gx(+metab)") & fdr.sig == TRUE, 
                  .(analysis = "attenuating.gx", 
                    metabolites = unique(metabolite)
                  )
                  ]
bmi.metabolites <- rbindlist(
  list(
    tmp1,tmp2,tmp3,tmp4,tmp5,tmp6
  ))
bmi.metabolites[,phenotype:="BMI"]
bmi.metabolites[, analysis := factor(
  analysis,
  levels=c("meta.all", "meta.sig", "bmi.assoc", "mediation.sig", "mediated.by.gx","attenuating.gx"),
  labels=c("meta.all" = "Available ", 
           "meta.sig" = "Associating with\nGene Expression", 
           "bmi.assoc" = "Associating with\nBMI",
           "mediation.sig" = "Exposure or\nMediator",
           "mediated.by.gx" = "Exposure",
           "attenuating.gx" = "Mediator"
  )
)]

# first part for metabolites
bp1 <- ggplot(bmi.metabolites,
              aes(x=analysis)) + #, fill=analysis
  geom_bar(position = position_dodge(),
           fill=pal.metabs[1]) + 
  geom_text(stat='count', aes(label=..count..), vjust=-.5,cex=5) +
  guides(fill=FALSE) + 
  scale_y_continuous(expand = expansion(mult = c(0, .1))) +
  # scale_fill_discrete_sequential("OrYel",rev=F)+
  ylab("Number of unique metabolites") +
  theme_gxMetab() +
  ggtitle("Metabolites") +
  theme(
    text=element_text(size=13),
    axis.text.x = element_text(angle=45,hjust=1,size=13),
    axis.title.x = element_blank(),
    plot.title = element_text(hjust = .5), 
    panel.grid.major.x = element_blank()
  )

#=================================================================================#
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#=================================================================================#

#' ## Plot: Barplot of mediated genes
# Plot: Barplot of mediated genes ----

tmp1 <- res.init.meta[symbol_INGENUITY!="", .(analysis = "meta.all", genes = unique(symbol_INGENUITY))]
tmp2 <- res.init.meta[significant==TRUE & symbol_INGENUITY!="", .(analysis = "meta.sig", genes = unique(symbol_INGENUITY))]
tmp3 <- res.joint[pheno=="log.bmi" & type %in% c("pheno~gx") & hfdr.sigREM == TRUE & symbol_INGENUITY!="", 
                  .(analysis = "bmi.assoc", 
                    genes = unique(symbol_INGENUITY)
                  )]
tmp4 <- res.joint[pheno=="log.bmi" & type %in% c("pheno~gx(+metab)","pheno~metab(+gx)") & fdr.sig == TRUE & symbol_INGENUITY!="", 
                  .(analysis = "mediation.sig", 
                    genes = unique(symbol_INGENUITY)
                  )]
tmp5 <- res.joint[pheno=="log.bmi" & type %in% c("pheno~gx(+metab)") & fdr.sig == TRUE & symbol_INGENUITY!="", 
                  .(analysis = "mediated.by.metab", 
                    genes = unique(symbol_INGENUITY)
                  )]
tmp6 <- res.joint[pheno=="log.bmi" & type %in% c("pheno~metab(+gx)") & fdr.sig == TRUE & symbol_INGENUITY!="", 
                  .(analysis = "attenuating.metab", 
                    genes = unique(symbol_INGENUITY)
                  )]
bmi.genes <- rbindlist(
  list(
    tmp1,tmp2,tmp3,tmp4,tmp5,tmp6
  ))
bmi.genes[,phenotype:="BMI"]
bmi.genes[, analysis := factor(
  analysis,
  levels=c("meta.all", "meta.sig","bmi.assoc", "mediation.sig", "mediated.by.metab", "attenuating.metab"),
  labels=c("meta.all" = "Available", 
           "meta.sig" = "Associating with\nMetabolites", 
           "bmi.assoc" = "Associating with\nBMI", 
           "mediation.sig" = "Exposure or\nMediator",
           "mediated.by.metab" = "Exposure",
           "attenuating.metab" = "Mediator"
  )
)]

bp2 <- ggplot(bmi.genes,aes(analysis)) + # , fill=analysis
  geom_bar(position = position_dodge(),
           fill=pal.genes[1]) + 
  geom_text(stat='count', aes(label=..count..), vjust=-.5, cex=5) +
  guides(fill=FALSE) + 
  scale_y_continuous(expand = expansion(mult = c(0, .1))) +
  # scale_fill_discrete_sequential("Teal",rev=F)+
  ylab("Number of unique genes") +
  theme_gxMetab() +
  ggtitle("Gene Expression") +
  theme(
    text=element_text(size=13),
    axis.text.x = element_text(angle=45,hjust=1, size=13),
    axis.title.x = element_blank(),
    plot.title = element_text(hjust = .5),
    panel.grid.major.x = element_blank()
  )

#=================================================================================#
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#=================================================================================#

#' ## Plot: Barplot of mediations
# Plot: Barplot of mediations ----

tmp1 <- res.init.meta[symbol_INGENUITY!="", .(analysis = "meta.all", assocs = unique(paste0(symbol_INGENUITY, metabolite)))]
tmp2 <- res.init.meta[symbol_INGENUITY!="" & significant==TRUE, .(analysis = "meta.sig", assocs = unique(paste0(symbol_INGENUITY, metabolite)))]

# paare, die mindestens 1 mal mit bmi assoziiern 
# das sind die, für die ich mediation getestet habe
tmp3 <- res.joint[pheno=="log.bmi" & type %in% c("pheno~gx(+metab)","pheno~metab(+gx)") & symbol_INGENUITY!="",
                  .(analysis = "bmi.assoc",
                    assocs = unique(paste0(symbol_INGENUITY, metabolite))
                  )]
tmp4 <- res.joint[pheno=="log.bmi" & type %in% c("pheno~gx(+metab)","pheno~metab(+gx)") & fdr.sig==TRUE & symbol_INGENUITY!="", 
                  .(analysis = "mediations.sig", 
                    assocs = unique(paste0(symbol_INGENUITY, metabolite))
                  )]
tmp5 <- res.joint[pheno=="log.bmi" & type %in% c("pheno~gx(+metab)") & fdr.sig == TRUE & symbol_INGENUITY!="", ][
  ,.(analysis = "mediated.gx", 
     assocs = unique(paste0(symbol_INGENUITY, metabolite))
  )]
tmp6 <- res.joint[pheno=="log.bmi" & type %in% c("pheno~metab(+gx)") & fdr.sig == TRUE & symbol_INGENUITY!="", ][
  ,.(analysis = "mediated.metab", 
     assocs = unique(paste0(symbol_INGENUITY, metabolite))
  )]
# this seems tricky. get all mediations that are duplicated and sig
tmp7 <- res.joint[pheno=="log.bmi" & 
                    type %in% c("pheno~gx(+metab)","pheno~metab(+gx)") & fdr.sig == TRUE & symbol_INGENUITY!="", ][
                      ,.(id = unique(paste0(symbol_INGENUITY,metabolite))),type][
  duplicated(id),
  .(analysis = "mediated.both", 
    assocs = id
  )]

bmi.pairs <- rbindlist(
  list(
    tmp1,tmp2,tmp3,tmp4,tmp5,tmp6,tmp7
  ))
bmi.pairs[,phenotype:="BMI"]
bmi.pairs[, analysis := factor(
  analysis,
  levels=c("meta.all","meta.sig", "bmi.assoc","mediations.sig", "mediated.gx", "mediated.metab","mediated.both"), 
  labels=c(
    "meta.all" = "Available",
    "meta.sig" = "Associating\nPairs", 
    "bmi.assoc" = "Associating with\nBMI",
    "mediations.sig" = "Mediating Pairs",
    "mediated.gx" = "Gene Expression is\nExposure",
    "mediated.metab" = "Metabolite is\nExposure",
    "mediated.both" = "Bi-directional\nMediation"
  )
)]
bmi.pairs[,dummy:=ifelse(analysis=="Available", "Available Pairs","")]
bmi.pairs$dummy <- as.factor(bmi.pairs$dummy)

# second plot for reduction in pairs by analysis
bp3.2 <- ggplot(bmi.pairs[dummy==""],
                aes(analysis)) + #, fill=analysis
  geom_bar(position = position_dodge(),
           fill=pal.pairs[1]) + 
  geom_text(stat='count', aes(label=..count..), vjust=-.5, cex=5) +
  scale_y_continuous(expand = expansion(mult = c(0, .1))) +
  guides(fill=FALSE) + 
  # facet_wrap(~dummy, scales="free",width=)+
  # scale_fill_discrete_sequential("Greens2",rev=F)+
  # ylab("Number of metabolite-probe pairs") +
  ggtitle("Metabolite-Gene Pairs") +
  theme_gxMetab() +
  theme(
    text=element_text(size=13),
    axis.text.x = element_text(angle=45,hjust=1),
    axis.title.x = element_blank(),
    plot.title = element_text(hjust = .5),
    axis.title.y = element_blank(),
    panel.grid.major.x = element_blank()
  )

#=================================================================================#
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#=================================================================================#

# Build Plot
plot_a <- plot_grid(
  ggplot(bmi.pairs[dummy==""]) + theme(panel.border = element_blank()) + labs(tag="A"),
  bp2 + labs(tag="C"),
  ncol = 1,
  nrow = 2,
  axis="b",
  align="v"
)
plot_b <- plot_grid(
  ggplot(bmi.pairs[dummy==""]) + theme(panel.border = element_blank()) + labs(tag="B"),
  bp1 + labs(tag="D"),
  ncol = 1,
  nrow = 2,
  axis="b",
  align="v"
)
plot_c <- plot_grid(
  bp3.2 + labs(tag="E"), 
  ncol = 1,
  nrow = 1,
  axis="b",
  align="v"
)

plot_all <- plot_grid(plot_a,plot_b,plot_c,axis="b", align="vh",ncol=3, nrow=1)


ggsave2(plot = plot_all,
        filename = dpx("mediationMetaboliteProbeOverviewBarplot.tiff","plt/"),
        dpi = 300,
        units = "in",
        width = 13,
        height = 13)

#=================================================================================#
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#=================================================================================#

# plot by gene/metabolite group

# Plot: Treemap of mediator classes ----

# no_data group is not informative -> change to "other"
res.init.meta[, uniqueN(symbol_INGENUITY),by=types_INGENUITY][order(-V1)] %>% as.data.frame

# add super path and ingenuity type
metabpathref <- res.init.meta[,.(path = unique(metab.super.path)),by=metabolite]
genepathref <- res.init.meta[,.(path = unique(types_INGENUITY)),by=symbol_INGENUITY]
mm <- metabpathref[match(r2.mediation$metabolite, metabolite),path]
gm <- genepathref[match(r2.mediation$gene, symbol_INGENUITY),path]
r2.mediation[, metab.path := mm]
r2.mediation[, gene.path := gm]

# data for metabolite centric plot
bp1 <- r2.mediation[pheno == "log.bmi" & mediation.sig == TRUE, .(metabolite,metab.path,type)]
bp1 <- bp1[,.(metabolite = unique(metabolite)),.(metab.path,type)]
bp1[, type := factor(type, labels = c("pheno~gx(+metab)" = "Metabolite-mediated Gx Effects",
                                      "pheno~metab(+gx)" = "Gx-mediated Metabolite Effects"))]

# 201112 Add 4th bar with unique mediations for each type (from meeting with markus)
bp1_all <- res.init.meta[, .(metabolite = unique(metabolite),
                            type = "Available Metabolites"),metab.super.path]
bp1_all$metab.path <- bp1_all$metab.super.path
bp1_all$metab.super.path <- NULL

# data for gene centric plot
bp2 <- r2.mediation[gene != "" & pheno == "log.bmi" & mediation.sig == TRUE, .(gene, gene.path, type)]
bp2 <- bp2[,.(gene = unique(gene)),.(gene.path,type)]
bp2[, type := factor(type, labels = c("pheno~gx(+metab)" = "Metabolite-mediated Gx Effects",
                                      "pheno~metab(+gx)" = "Gx-mediated Metabolite Effects"))]

# 201112 Add 4th bar with unique mediations for each type (from meeting with markus)
bp2_all <- res.init.meta[symbol_INGENUITY!="", .(gene = unique(symbol_INGENUITY),
                            type = "Available Genes"),types_INGENUITY]
bp2_all$gene.path <- bp2_all$types_INGENUITY
bp2_all$types_INGENUITY <- NULL

# third plot with pairs (and metabolite pathway)
bp3_all <- r2.mediation[!is.na(mediation.sig) & gene != "" & pheno == "log.bmi", 
                              .(pair = unique(paste0(metabolite, gene))),
                              .(gene.path,metab.path,type)]
bp3_all$type <- "All tested pairs"

bp3 <- r2.mediation[gene != "" & pheno == "log.bmi" & mediation.sig == TRUE, 
                    .(pair = unique(paste0(metabolite,gene))),
                    .(metab.path, gene.path, type)]
bp3[, type := factor(type, labels = c("pheno~gx(+metab)" = "Metabolite-mediated Gx Effects",
                                      "pheno~metab(+gx)" = "Gx-mediated Metabolite Effects"))]
bp3_intersect <- bp3[, .(
  pair=intersect(
    pair[type=="Metabolite-mediated Gx Effects"],
    pair[type=="Gx-mediated Metabolite Effects"]
  ),
  type = "Intersect"
), by = .(metab.path,gene.path)]
bp3_intersect <- bp3_intersect[!is.na(pair)]

# add common metabolites/genes to plot

bp1 <- rbindlist(
  use.names = T,
  list(
    bp1_all,
    bp1,
    # nice i thought that wouldnt work
    bp1[, .(
      metabolite=intersect(
        metabolite[type=="Metabolite-mediated Gx Effects"],
        metabolite[type=="Gx-mediated Metabolite Effects"]
      ),
      type = "Intersect"
    ), by = metab.path]
  )
)

bp2 <- rbindlist(
  use.names = TRUE,
  list(
    bp2_all,
    bp2,
    bp2[, .(
      gene=intersect(
        gene[type=="Metabolite-mediated Gx Effects"],
        gene[type=="Gx-mediated Metabolite Effects"]
      ),
      type = "Intersect"
    ), by = gene.path]
  ))

# the pair object
bp3 <- rbindlist(
  use.names = T,
  fill = T,
  list(
    bp3,
    bp3_all,
    bp3_intersect
  )
)

# reorder pathway by most features
bp1.ord <- bp1[, .(ord = uniqueN(metabolite)), by=metab.path][order(-ord)]
bp1[,metab.path.fac := factor(metab.path, levels = bp1.ord$metab.path)]
bp2.ord <- bp2[, .(ord = uniqueN(gene)), by=gene.path][order(-ord)]
bp2[,gene.path.fac := factor(gene.path, levels = bp2.ord$gene.path)]

bp3.ord.met <- bp3[, .(ord = uniqueN(pair)), by=.(metab.path)][order(-ord)]
bp3[,metab.path.fac := factor(metab.path, levels = bp3.ord.met$metab.path)]
bp3.ord.gene <- bp3[, .(ord = uniqueN(pair)), by=.(gene.path)][order(-ord)]
bp3[,gene.path.fac := factor(gene.path, levels = bp3.ord.gene$gene.path)]

#=================================================================================#
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#=================================================================================#

# Plot: Barplot of mediator classes ----

# superpathway colors (from circus plot)
gl <- c("Amino acid metabolism, other", 
        "Folic acid cycle", 
        "Ammonia recycling", 
        "BCAA metabolism", 
        "Collagen synthesis", 
        "Energy metabolism", 
        "Carnitine transport", 
        "Fatty acid metabolism", 
        "Other", 
        "Urea cycle")
superpath_color <- structure(
  c("#E16A86", "#CB7F2F", "#9F9400", 
    "#50A315", "#00AC79", "#00AAB7", 
    "#009ADE", "#A87BE4", "#DA65C3", 
    "#E16A86"),
  names=gl)

# create an ordered factor to reorder the y axis in the plot below (weird behaviour due to flip)
bp1[, type.f := factor(
  type, 
  levels = c(
    "Intersect",
    "Gx-mediated Metabolite Effects", 
    "Metabolite-mediated Gx Effects",
    "Available Metabolites"
             ),
  labels = c(
    "Intersect" = "Both",
    "Gx-mediated Metabolite Effects" = "Metabolite->Gene expression->BMI", 
    "Metabolite-mediated Gx Effects" = "Gene expression->Metabolite->BMI",
    "Available Metabolites" = "Available Metabolites"
  ),
  ordered = T)]

p1 <- ggplot(
  bp1,
  aes(
    x = type.f,
    fill = metab.path.fac
  )) +
  geom_bar(position="fill",lwd=.5,col="grey30") + 
  scale_fill_manual(values = superpath_color) +
  guides(fill = guide_legend(title = "Metabolite pathway")) +
  scale_y_continuous(
    breaks = c(0,0.25,0.5,0.75,1),
    labels = c(0,25,50,75,100)
  ) +
  labs(y = "Unique metabolites [%]") +
  theme_gxMetab() +
  theme(
    plot.margin = unit(c(0,1,0.3,1),"cm"),
    axis.title.y = element_blank(),
    # axis.text.x = element_text(angle = 45, 
    #                            hjust = 1),
    # legend.title = element_blank(),
    legend.justification = c(0,1), 
    legend.position = "right",
    legend.direction = "vertical",
    panel.grid.major.x = element_line(color = "grey30",
                                      linetype = "dotted")
    
  ) +
  coord_flip()
  
# custom palette for all the gene classes
gene_pal <- colorspace::qualitative_hcl(n = uniqueN(bp2$gene.path), h = c(60, 240), c = 50, l = 70, register = "Custom-Palette")
gene_pal <- sample(gene_pal, size = length(gene_pal))

# new ordering of axis
bp2[, type.f := factor(
  type, 
  levels = c(
    "Intersect",
    "Gx-mediated Metabolite Effects", 
    "Metabolite-mediated Gx Effects",
    "Available Genes"
  ),
  labels = c(
    "Intersect" = "Both",
    "Gx-mediated Metabolite Effects" = "Metabolite->Gene expression->BMI", 
    "Metabolite-mediated Gx Effects" = "Gene expression->Metabolite->BMI",
    "Available Genes"
  ),
  ordered = T)]

p2 <- ggplot(
  bp2,
  aes(
    x = type.f,
    fill = gene.path.fac
  )) +
  geom_bar(position="fill",lwd=.5,col="grey30") + 
  scale_fill_manual(values = gene_pal) +
  guides(fill = guide_legend(title = "Gene pathway")) +
  scale_y_continuous(
    breaks = c(0,0.25,0.5,0.75,1),
    labels = c(0,25,50,75,100)
  ) +
  labs(y = "Unique genes [%]") +
  theme_gxMetab() +
  theme(
    plot.margin = unit(c(0,1,0.3,1),"cm"),
    axis.title.y = element_blank(),
    # axis.text.x = element_text(angle = 45, 
    #                            hjust = 1),
    # legend.title = element_blank(),
    legend.justification = c(0,1), 
    legend.position = "right",
    legend.direction = "vertical",
    panel.grid.major.x = element_line(color = "grey30",
                                      linetype = "dotted")
    
  ) +
  coord_flip()

# Pair plot
bp3[, type.f := factor(
  type, 
  levels = c(
    "Intersect",
    "Gx-mediated Metabolite Effects", 
    "Metabolite-mediated Gx Effects",
    "All tested pairs"
  ),
  labels = c(
    "Intersect" = "Both",
    "Gx-mediated Metabolite Effects" = "Metabolite->Gene expression->BMI", 
    "Metabolite-mediated Gx Effects" = "Gene expression->Metabolite->BMI",
    "All tested pairs"
  ),
  ordered = T)]

# try plotting...
p3.1 <- ggplot(
  bp3,
  aes(
    x = type.f,
    fill = gene.path.fac
  )) +
  geom_bar(position="fill",lwd=.5,col="grey30") +
  scale_fill_manual(values = gene_pal) +
  scale_y_continuous(
    breaks = c(0,0.25,0.5,0.75,1),
    labels = c(0,25,50,75,100)
  ) +
  labs(y = "Unique pairs [%]") +
  theme_gxMetab() +
  theme(
    plot.margin = unit(c(0,1,0.3,1),"cm"),
    axis.title.y = element_blank(),
    # axis.text.x = element_text(angle = 45, 
    #                            hjust = 1),
    legend.title = element_blank(),
    legend.justification = c(0,1), 
    legend.position = "right",
    legend.direction = "vertical",
    panel.grid.major.x = element_line(color = "grey30",
                                      linetype = "dotted")
    
  ) +
  coord_flip()

p3.2 <- ggplot(
  bp3,
  aes(
    x = type.f,
    fill = metab.path.fac
  )) +
  geom_bar(position="fill",lwd=.5,col="grey30") +
  scale_fill_manual(values = superpath_color) +
  scale_y_continuous(
    breaks = c(0,0.25,0.5,0.75,1),
    labels = c(0,25,50,75,100)
  ) +
  labs(y = "Unique pairs [%]") +
  theme_gxMetab() +
  theme(
    plot.margin = unit(c(0,1,0.3,1),"cm"),
    axis.title.y = element_blank(),
    # axis.text.x = element_text(angle = 45, 
    #                            hjust = 1),
    legend.title = element_blank(),
    legend.justification = c(0,1), 
    legend.position = "right",
    legend.direction = "vertical",
    panel.grid.major.x = element_line(color = "grey30",
                                      linetype = "dotted")
    
  ) +
  coord_flip()

# arrange the plot
pg1 <- plot_grid(p1 + theme(legend.position="none"),
          p2 + theme(legend.position="none"),
          p3.2 + theme(legend.position="none"),
          p3.1 + theme(legend.position="none"),
          labels = "AUTO", align = "v", nrow = 4, ncol =1)

legend2 <- get_legend(p2)
legend1 <- get_legend(p1)
legends <- plot_grid(
  legend1,
  legend2,
  align = "v",
  nrow = 2,
  ncol=1,
  rel_heights = c(10,15)
)

pg <- plot_grid(
  pg1,
  legends,
  rel_widths = c(1,.3)
)

tiff(dpx("uniqueFeaturesPerPath.tiff", "plt/"), 
     width = 13,
     height = 7,
     res=300,
     units = "in")
 
pg

dev.off()

#=================================================================================#
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#=================================================================================#

# Difference in pathway proportions ----

# metabolites
path_names_metab <- as.character(unique(bp1$metab.path.fac))
chi_metab <- dcast(
  bp1[, .N, .(type.f, metab.path.fac)], 
  type.f ~ metab.path.fac, value.var = "N")
change_in_dt(chi_metab, from = NA, to = 0, change_in_dat = T)

# genes
chi_genes <- dcast(
  bp2[, .N, .(type.f, gene.path.fac)], 
  gene.path.fac ~ type.f)
change_in_dt(chi_genes, from = NA, to = 0, change_in_dat = T)

# pairs - metabolites
chi_pair_metab <- dcast(bp3[, .N, .(type.f, metab.path.fac)], 
                        metab.path.fac ~ type.f)
change_in_dt(chi_pair_metab, from = NA, to = 0, change_in_dat = T)

# pairs - genes
chi_pair_genes <- dcast(bp3[, .N, .(type.f, gene.path.fac)], 
                        gene.path.fac ~ type.f)
change_in_dt(chi_pair_genes, from = NA, to = 0, change_in_dat = T)

# calculations for metabolites
n_observed <- chi_metab[1:3, as.matrix(.SD), .SDcols = path_names_metab]
row_sums <- chi_metab[, rowSums(.SD), .SDcols = path_names_metab]
names(row_sums) <- chi_metab$type.f
n_total <- row_sums["Available Metabolites"]
row_sums <- row_sums[1:3]
col_sums <- unlist(chi_metab[, lapply(.SD, sum),.SDcols=path_names_metab])

# calc the expected occurence
n_expected <- matrix(nrow = length(row_sums),
              ncol = length(col_sums), 
              dimnames = list(names(row_sums),
                              names(col_sums)),
              data = NA)
res2 <- n_expected

# erwartete häufigkeiten
for(i in 1:length(row_sums)){
  for(j in 1:length(col_sums)){
    n_expected[i,j] <- ( row_sums[i]*col_sums[j] ) / n 
  }
}

# get chisq statistic
chisq <- sum((n_observed - n_expected) / n_expected)


# df
df <- (length(row_sums)-1) * (length(col_sums)-1)
pchisq(abs(chisq), df = df, lower.tail = F)


#=================================================================================#
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#=================================================================================#

# Plot: Treemap of most mediat* genes/metabolites ----

# mediation results
d <- res.joint[type%in%c("pheno~gx(+metab)","pheno~metab(+gx)") & pheno=="log.bmi" & symbol_INGENUITY!=""]

# add super path and ingenuity type
metabclassref <- res.init.meta[,.(class = unique(metab.class)),by=metabolite]
metabpathref <- res.init.meta[,.(path = unique(metab.super.path)),by=metabolite]
genepathref <- res.init.meta[,.(path = unique(types_INGENUITY)),by=symbol_INGENUITY]
mm <- metabpathref[match(d$metabolite, metabolite),path]
mc <- metabclassref[match(d$metabolite, metabolite),class]
gm <- genepathref[match(d$symbol_INGENUITY, symbol_INGENUITY),path]
d[, metab.class := mc]
d[, metab.path := mm]
d[, gene.path := gm]

# data for gene centric plot
tm.gbm <- d[fdr.sig == TRUE, uniqueN(symbol_INGENUITY), by = .(type, metabolite)][order(-V1)][,head(.SD,25),by=type]
tm.gbm[, type := factor(type, labels = c("pheno~gx(+metab)" = "Metabolite-mediated Gx Effects",
                                         "pheno~metab(+gx)" = "Gx-mediated Metabolite Effects"))]

tm.mbg <- d[fdr.sig == TRUE, uniqueN(metabolite), by = .(type, symbol_INGENUITY)][order(-V1)][,head(.SD,25),by=type]
tm.mbg[, type := factor(type, labels = c("pheno~gx(+metab)" = "Metabolite-mediated Gx Effects",
                                         "pheno~metab(+gx)" = "Gx-mediated Metabolite Effects"))]

# Plot
tiff(dpx("treemapTopMetabGeneMediationBMI.tiff","plt/"),res=300,unit="in",height=9,width=13)

grid.newpage()
pushViewport(viewport(layout=grid.layout(1,2)))
vp <- viewport(layout.pos.col=1, layout.pos.row=1)

treemap(dtf = tm.gbm,
        index=c("type", "metabolite"),
        vSize="V1",
        vColor="V1",
        palette = rev(pal.metabs),
        overlap.labels = 1,
        border.col=c("grey70","grey35"),
        border.lwds = c(7,2),
        type="manual",
        title = "Top 25 metabolites involved in mediation on BMI",
        title.legend = "Number of unique genes involved with metabolite",
        vp=vp)

vp <- viewport(layout.pos.col=2, layout.pos.row=1)
treemap(dtf = tm.mbg,
        index=c("type", "symbol_INGENUITY"),
        vSize="V1",
        vColor="V1",
        palette = rev(pal.genes),
        overlap.labels = 1,
        border.col=c("grey70","grey35"),
        border.lwds = c(7,2),
        type="manual",
        title = "Top 25 genes involved in mediation on BMI",
        title.legend = "Number of unique metabolites involved with gene",
        vp=vp)

dev.off()

#=================================================================================#
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#=================================================================================#

# Plot: Histogram of mediation betas ----
#' ## Plot: Histogram of mediation betas

p.hist <- ggplot(r2.mediation[pheno=="log.bmi" & mediation.sig == TRUE,],
                 aes(x=beta.mediation, fill=type)) +
  geom_histogram(col="grey35",
                 alpha=0.9,
                 bins=50
  )+
  geom_vline(xintercept = 0,lwd=0.6,col="indianred4",lty="dashed") +
  scale_fill_manual(values = c("pheno~gx(+metab)" = pal.genes[1],
                               "pheno~metab(+gx)" = pal.metabs[1]),
                    name=NULL,
                    labels = c(
                      "pheno~gx(+metab)" = "Mediation of gene expression effects",
                      "pheno~metab(+gx)" = "Mediation of metabolite effects"
                    )) +
  theme_gxMetab() + 
  labs(title = expression("Distribution of mediation effect estimates "*(beta["mediation"])),
       subtitle = "All significant mediation pairs (hierarchical FDR=5%)", 
       y=NULL,
       x=expression(beta["mediation"])) +
  theme(
    legend.position = "bottom", 
    legend.direction = "vertical",legend.justification = c(0,0)
  )

print(p.hist)

tiff(dpx("distributionMediationBetaBMI.tiff","plt/"),res=300,unit="in",height=8,width=8)
print(p.hist)
dev.off()

#=================================================================================#
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#=================================================================================#

# Plot: Histogram of relative beta-reduction ----
#' ## Plot: Histogram of relative beta-reduction

# Measure of proportion mediated from: The intermediate endpoint effect in logistic and probit regression paper
# 1 - (tauAdj/tau)
# gamma_m * alpha_z / beta_z
# gamma_m * alpha_z / (gamma_z + gamma_m * alpha_z) = beta * alpha / (tauAdj + beta * alpha)
# -> vermeide nutzung von tau (siehe paper!!!)

r2.mediation[mediation.sig == TRUE, prop.mediated.old := 1-(tauAdj/tau)]
r2.mediation[mediation.sig == TRUE, prop.mediated.new := beta * alpha / (tauAdj + beta * alpha)]
r2.mediation[mediation.sig == TRUE, prop.mediated.tmp := beta * alpha / (tau)]
with(r2.mediation, lm(prop.mediated.new ~ prop.mediated.old-1)) %>% summary

par(mfrow=c(1,3))
with(r2.mediation, plot(prop.mediated.new ~ prop.mediated.old, 
                        xlab="1-(tau'/tau)", 
                        ylab="alpha*beta/(alpha*beta+tau')",
                        col=alpha("grey35",0.5),
                        pch=20))
abline(0,1,col="indianred4")
with(r2.mediation, plot(prop.mediated.new ~ prop.mediated.tmp, 
                        xlab="alpha*beta/tau",
                        ylab="alpha*beta/(alpha*beta+tau')",
                        col=alpha("grey35",0.5),
                        pch=20))
abline(0,1,col="indianred4")
with(r2.mediation, plot(prop.mediated.old ~ prop.mediated.tmp, 
                        xlab="alpha*beta/tau", 
                        ylab="1-(tau'/tau)",
                        col=alpha("grey35",0.5),
                        pch=20))
abline(0,1,col="indianred4")
par(mfrow=c(1,1))

# consolidate
r2.mediation$prop.mediated <- r2.mediation$prop.mediated.new
r2.mediation$prop.mediated.new <- NULL
r2.mediation$prop.mediated.old <- NULL
r2.mediation$prop.mediated.tmp <- NULL

# plot
p.hist <- ggplot(r2.mediation[pheno=="log.bmi" & mediation.sig == TRUE,],
                 aes(x=prop.mediated, fill=type)) +
  geom_histogram(col="white",lwd=.3,
                 alpha=0.9,
                 bins=50
  )+
  geom_vline(xintercept = r2.mediation[pheno=="log.bmi" & 
                                         mediation.sig == TRUE & 
                                         type=="pheno~gx(+metab)",max(prop.mediated)],
             col=pal.genes[1],lwd=0.6,lty="dashed") +
  geom_vline(xintercept = r2.mediation[pheno=="log.bmi" & 
                                         mediation.sig == TRUE & 
                                         type=="pheno~gx(+metab)",median(prop.mediated)],
             col=pal.genes[1],lwd=0.6,lty="solid") +
  geom_vline(xintercept = r2.mediation[pheno=="log.bmi" & 
                                         mediation.sig == TRUE & 
                                         type=="pheno~metab(+gx)",max(prop.mediated)],
             col=pal.metabs[1],lwd=0.6,lty="dashed") +
  geom_vline(xintercept = r2.mediation[pheno=="log.bmi" & 
                                         mediation.sig == TRUE & 
                                         type=="pheno~metab(+gx)",median(prop.mediated)],
             col=pal.metabs[1],lwd=0.6,lty="solid") +
  
  scale_color_manual(values = c("pheno~gx(+metab)" = pal.genes[1],
                                "pheno~metab(+gx)" = pal.metabs[1]),
                     guide=FALSE)+
  scale_fill_manual(values = c("pheno~gx(+metab)" = pal.genes[1],
                               "pheno~metab(+gx)" = pal.metabs[1]),
                    name=NULL,
                    labels = c(
                      "pheno~gx(+metab)" = "Gene expression->Metabolite->BMI",
                      "pheno~metab(+gx)" = "Metabolite->Gene expression->BMI"
                    )) +
  scale_x_continuous(limits = c(0,NA))+
  theme_gxMetab() + 
  labs(
    # title = "Distribution of raw effect reduction by mediation",
    # subtitle = "All significant mediation pairs (hierarchical FDR=5%)", 
    y=NULL,
    x=expression("Proportion mediated effect")
    ) +
  theme(
    legend.position = "bottom", 
    legend.direction = "vertical",
    legend.justification = c(0,0)
  )

print(p.hist)

tiff(dpx("MediatedEffectProportionHistogramBMI.tiff","plt/"),res=300,unit="in",height=8,width=8)
print(p.hist)
dev.off()


#=================================================================================#
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#=================================================================================#


# Plot: Boxplot of reduction of direct effect ----
#' ## Plot: Boxplot of reduction of direct effect

# quotient for coloring

# new col for raw effect
r2.mediation[, beta.tot := alpha*beta+tauAdj]

beta.mediation.m <- melt(
  r2.mediation, 
  id.vars = c("type","pheno","gx.probe", "metabolite", "mediation.sig", "raw.sig.group","prop.mediated"),
  measure.vars = c("beta.tot", "tauAdj"), 
  variable.name = "beta.type",
  value.name = "beta" 
)

# pd <- position_dodge(width = .2, height = 0, seed = 1)
beta.mediation.m[, beta.type.num := as.numeric(as.factor(beta.type))]
beta.mediation.m[, beta.type.num := jitter(beta.type.num, 1)]

p.box <- ggplot(beta.mediation.m[mediation.sig==TRUE & 
                                   pheno=="log.bmi", ],
                aes(x=beta.type.num,
                    y=beta,
                    group=paste0(pheno,gx.probe,metabolite)))+
  geom_point(pch=21,alpha=.8,col="white",
             aes(fill=type,
                 size=abs(beta))) +
  scale_size_continuous(range = c(0.7,3),guide=FALSE) +
  scale_fill_manual(values = c("pheno~gx(+metab)" = pal.genes[1],
                               "pheno~metab(+gx)" = pal.metabs[1]),
                    guide=FALSE) +
  geom_line(col="grey60",alpha=.1) +
  facet_grid(~type, labeller = labeller(
    "type" = c("pheno~gx(+metab)" = "Gene Expression",
               "pheno~metab(+gx)" = "Metabolites"))) +
  scale_x_continuous(breaks = c(1,2),
                     labels = c("Raw Effect", "Direct Effect")) +
  geom_hline(yintercept = 0,lwd=1,col="grey35") +
  labs(y = "Standardized Effect estimate",
       x = NULL,
       title = "Comparison of direct (adjusted for a mediator)\nand raw effect of gene expression\nand metabolite levels on log-BMI"
  ) + 
  theme_gxMetab() + 
  theme(legend.position = "bottom", 
        legend.justification = c(0,0)) 

# Plot
# p.box

# Save for posteriority
# tiff(dpx("PropMediatedBetaBMIBoxplot.tiff","plt/"),res=300,unit="in",height=9,width=7)
# print(p.box)
# dev.off()

#=================================================================================#
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#=================================================================================#

# Mediation Groups - BMI ----
#' # Mediation Groups - BMI

# cast relevant information
relative.mediation <- dcast(r2.mediation, 
                            pheno + gx.probe + gene + metabolite ~ type,
                            value.var = c("r2.arnd",
                                          "beta.mediation",
                                          "mediation.sig", 
                                          "r2.gx.raw",
                                          "r2.metab.raw",
                                          "beta.gx.raw",
                                          "beta.metab.raw",
                                          "sig.gx.raw",
                                          "sig.metab.raw",
                                          "raw.sig.group",
                                          "prop.mediated"
                            ))

# set the NA significancies to FALSE
relative.mediation[is.na(`mediation.sig_pheno~gx(+metab)`), `mediation.sig_pheno~gx(+metab)` := FALSE]
relative.mediation[is.na(`mediation.sig_pheno~metab(+gx)`), `mediation.sig_pheno~metab(+gx)` := FALSE]
relative.mediation[is.na(`sig.metab.raw_pheno~gx(+metab)`), `sig.metab.raw_pheno~gx(+metab)` := FALSE]
relative.mediation[is.na(`sig.metab.raw_pheno~metab(+gx)`), `sig.metab.raw_pheno~metab(+gx)` := FALSE]
relative.mediation[is.na(`sig.gx.raw_pheno~gx(+metab)`), `sig.gx.raw_pheno~gx(+metab)` := FALSE]
relative.mediation[is.na(`sig.gx.raw_pheno~metab(+gx)`), `sig.gx.raw_pheno~metab(+gx)` := FALSE]

# small check
relative.mediation[is.na(`mediation.sig_pheno~gx(+metab)`),]
relative.mediation[is.na(`mediation.sig_pheno~metab(+gx)`),]

# Define the 3 mediatoin groups

# Update: Groups based on proportion mediated

# based on proportion mediated!
relative.mediation[((`prop.mediated_pheno~gx(+metab)` < 0.2) |
                      is.na(`prop.mediated_pheno~gx(+metab)`) == TRUE) & (
                        (`mediation.sig_pheno~metab(+gx)` == TRUE) &
                          `prop.mediated_pheno~metab(+gx)` >= 0.2), 
                   mediation.group := "Gene expression mediated metabolite effects"]

relative.mediation[((`prop.mediated_pheno~metab(+gx)` < 0.2) | 
                      is.na(`prop.mediated_pheno~metab(+gx)`) == TRUE) & (
                        (`mediation.sig_pheno~gx(+metab)` == TRUE) &
                          `prop.mediated_pheno~gx(+metab)` >= 0.2), 
                   mediation.group := "Metabolite mediated gene expression effects"]

relative.mediation[(`prop.mediated_pheno~metab(+gx)` >= 0.2 & `prop.mediated_pheno~gx(+metab)` >= .2) & 
                     (`mediation.sig_pheno~gx(+metab)` == TRUE & `mediation.sig_pheno~metab(+gx)` == TRUE),
                   mediation.group := "Bi-directional mediation"]

relative.mediation[(is.na(`mediation.sig_pheno~gx(+metab)`) | `mediation.sig_pheno~gx(+metab)` == FALSE) & 
                     (is.na(`mediation.sig_pheno~metab(+gx)`) | `mediation.sig_pheno~metab(+gx)` == FALSE),
                   mediation.group := "None"]

relative.mediation[is.na(mediation.group) & pheno=="log.bmi", mediation.group := "Other"]
relative.mediation[is.na(mediation.group) & pheno=="diabetes.status.tri", mediation.group := "Other"]
relative.mediation[pheno=="log.bmi", table(mediation.group)]
# relative.mediation[pheno=="diabetes.status.tri", mediation.group := "None"]

# check
relative.mediation[, unique(mediation.group),by=pheno]
relative.mediation$mediation.group %>% unique
relative.mediation$mediation.group %>% is.na %>% any

# # match genes to info
# m1 <- match(relative.mediation$gx.probe, res.init.meta$markerID)
# relative.mediation[, gene := res.init.meta[m1, symbol_INGENUITY]]

# Check results for BMI, not T2D
relative.mediation[pheno=="log.bmi", .(.N)]
relative.mediation[pheno=="log.bmi", .(.N), by=mediation.group]
relative.mediation[pheno=="log.bmi", .(genes = uniqueN(gene),
                                       metabolites = uniqueN(metabolite)), 
                   by=mediation.group]


# save results
fwrite(relative.mediation, dpx("mediatonGroupList","res/"),sep = "\t")
# relative.mediation <- fread(newest_file("mediatonGroupList","res/",print_full = T))

#' ## Characterize Groups
# Characterize Groups ----

# Overlap of genes/metabolites in group
my.groups <- relative.mediation$mediation.group %>% unique

# remove empty gene mediations
relative.mediation[,sum(gene == "")]
relative.mediation <- relative.mediation[gene != "",]

# basic numbers
relative.mediation[pheno=="log.bmi", .N]
relative.mediation[pheno=="log.bmi" & (`mediation.sig_pheno~gx(+metab)` == TRUE |  `mediation.sig_pheno~metab(+gx)` == TRUE),.N]
relative.mediation[pheno=="log.bmi", table(mediation.group)]

# total genes/metabolites in all groups
relative.mediation[pheno=="log.bmi" & mediation.group %in% my.groups[c(3,4)], .(genes=uniqueN(gene),
                                                                                metabolites=uniqueN(metabolite))]

# number of genes/metabolites per group
relative.mediation[pheno=="log.bmi", .(genes=uniqueN(gene),
                                       metabolites=uniqueN(metabolite)),by=(mediation.group)]

# genes
venn3(
  relative.mediation[pheno=="log.bmi" & mediation.group == my.groups[4], unique(gene)],
  relative.mediation[pheno=="log.bmi" & mediation.group == my.groups[3], unique(gene)],
  relative.mediation[pheno=="log.bmi" & mediation.group == my.groups[2], unique(gene)],
  mylabels = c("m->gx->bmi","gx->m->bmi","Other"),
  mytitle = "Unique Genes"
)

# metabolites
venn3(
  relative.mediation[pheno=="log.bmi" & mediation.group == my.groups[4], unique(metabolite)],
  relative.mediation[pheno=="log.bmi" & mediation.group == my.groups[3], unique(metabolite)],
  relative.mediation[pheno=="log.bmi" & mediation.group == my.groups[2], unique(metabolite)],
  mylabels = c("m->gx->bmi","gx->m->bmi","Other"),
  mytitle = "Unique Metabolites"
)

# TODO: Supp table with bmi categorization ----

#=================================================================================#
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#=================================================================================#

# Plot: Relative Mediation   ----
#' # Plot: Relative Mediation 

# checking
r2.mediation[,.(
  sig = sum(mediation.sig==T,na.rm = T),
  all = .N,
  range = paste(signif(range(prop.mediated,na.rm = T),3),collapse = ",")),by=type]
relative.mediation[!is.na(`mediation.sig_pheno~gx(+metab)`)][
  `mediation.sig_pheno~gx(+metab)`==T,
  .(sig= .N,
    range=paste(signif(range(`prop.mediated_pheno~gx(+metab)`,na.rm = T),3),collapse = ","))]
relative.mediation[!is.na(`mediation.sig_pheno~metab(+gx)`)][
  `mediation.sig_pheno~metab(+gx)`==T,
  .(sig = .N,
    range=paste(signif(range(`prop.mediated_pheno~metab(+gx)`,na.rm = T),3),collapse = ","))]

# range per group
relative.mediation[mediation.group %in% my.groups[c(2,3,4)], .(
  range=paste(signif(range(`prop.mediated_pheno~gx(+metab)`,na.rm = T),3),collapse = ","),
  range=paste(signif(range(`prop.mediated_pheno~metab(+gx)`,na.rm = T),3),collapse = ",")
),by = mediation.group] 

#=================================================================================#
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#=================================================================================#

# melt for plotting
rel.med.m <- melt(relative.mediation[mediation.group != "None"], 
                  id.vars = c("pheno","gx.probe","gene", "metabolite", "mediation.group",
                              "r2.gx.raw_pheno~gx(+metab)",
                              "r2.metab.raw_pheno~metab(+gx)"
                  ),
                  measure.vars = c("prop.mediated_pheno~gx(+metab)","prop.mediated_pheno~metab(+gx)"), 
                  variable.name = "type",
                  value.name = "prop.mediated" 
)

# color the 'other' groups by metab/gx
rel.med.m[, mediation.group.col := mediation.group]
rel.med.m[mediation.group.col == "Other", mediation.group.col := ifelse(
  type == "prop.mediated_pheno~gx(+metab)", 
  "Other (Exposure: Gene Expression)",
  "Other (Exposure: Metabolite)")]

# turn label into numeric for plotting
rel.med.m[,r2 := ifelse(type=="prop.mediated_pheno~gx(+metab)", `r2.gx.raw_pheno~gx(+metab)`, `r2.metab.raw_pheno~metab(+gx)`)]

# set mediation.group as char->factor->numeric 
rel.med.m[, med.group.num := as.numeric(as.factor(mediation.group))]

# create nice verbose groups with offset for plotting
rel.med.m[med.group.num == 1 & type == "prop.mediated_pheno~gx(+metab)", med.group.num := 0.8]
rel.med.m[med.group.num == 1 & type == "prop.mediated_pheno~metab(+gx)", med.group.num := 1.2]
rel.med.m[med.group.num == 2 & type == "prop.mediated_pheno~gx(+metab)", med.group.num := 1.8]
rel.med.m[med.group.num == 2 & type == "prop.mediated_pheno~metab(+gx)", med.group.num := 2.2]
rel.med.m[med.group.num == 3 & type == "prop.mediated_pheno~gx(+metab)", med.group.num := 2.8]
rel.med.m[med.group.num == 3 & type == "prop.mediated_pheno~metab(+gx)", med.group.num := 3.2]


# jitter points for visibility 
rel.med.m[, med.group.num := jitter(med.group.num, 2)]
rel.med.m[is.na(prop.mediated)] %>% head
setkey(rel.med.m, prop.mediated)

# merge long name
# rel.med.m$metab.long <- agn.all[match(rel.med.m$metabolite, rid), abbr]
rel.med.m[,id := paste0(gene,gx.probe,metabolite,type)]

# label the top mediations
ids <- rel.med.m[mediation.group != "Other" & pheno == "log.bmi",tail(.SD,5), by = mediation.group][
  , .(id=paste0(gene,gx.probe,metabolite,type))]
rel.med.m[pheno=="log.bmi" & id %in% (ids$id), label := paste0("Metabolite: ", metabolite,"\nGene: ", gene)]

# plot
p.box <- ggplot(rel.med.m[mediation.group%nin% c("None","Other") & pheno=="log.bmi", ],
                aes(x=med.group.num,
                    y=prop.mediated,
                    col=type,
                    group=paste0(pheno,gx.probe,metabolite)
                ))+
  geom_line(
    data = subset(rel.med.m[mediation.group!="None" & pheno=="log.bmi", ],
                  mediation.group %nin% c("Other")),
    col="grey60",alpha=.1) +
  geom_point(alpha=.6,aes(size=r2)) +
  geom_text_repel(aes(label = label),
                  col           = "grey29",
                  cex           = 3.4,
                  force         = 10,
                  box.padding   = 0.5, 
                  nudge_y = .3,
                  point.padding = 0.5,
                  ylim=c(.40,1)) +
  scale_size_continuous(range = c(0.7,3),
                        guide=guide_legend(
                          title = "Explained variance of\nexposure in log-BMI"
                        )) +
  scale_color_manual(values = c(
    "prop.mediated_pheno~gx(+metab)"=pal.genes[1],
    "prop.mediated_pheno~metab(+gx)"=pal.metabs[1]
  ),
  labels = c(
    "prop.mediated_pheno~gx(+metab)"="Gene expression->Metabolite->BMI",
    "prop.mediated_pheno~metab(+gx)"="Metabolite->Gene expression->BMI"
    ),
  guide=guide_legend(title = "Mediation Analysis Direction")) +
  scale_y_continuous(limits = c(0,1)) + 
  scale_x_continuous(breaks = c(1,2),
                     labels = c(
                       "Category:\nGene expression mediated\nmetabolite effects",
                       "Category:\nMetabolite mediated gene\nexpression effects"
                       # "Category:\nOther"
                       )) +
  labs(y = "Proportion Mediated",
       x = NULL
  ) +
  theme_gxMetab() + 
  theme(legend.position = "bottom", 
        legend.justification = c(0,0),
        legend.direction = "vertical")

# tiff(dpx("bothDirectionPropmediatedEffectBetaBMIBoxplot.tiff","plt/"),res=300,unit="in",height=10,width=9)
# print(p.box)
# dev.off()

# Joint Plot of histogram of PM and boxplot of effect reductions in both directions
tiff(dpx("propMediatedHistAndBoxplot.tiff","plt/"),res=300,unit="in",height=10,width=14)
plot_grid(
  p.hist + guides(fill = guide_legend(title = "Mediation Analysis Direction")),
  p.box + guides(color = "none"),
          labels = "AUTO",
          align = "hv",
          axis = "b",
          ncol = 2, 
          nrow=1, 
          rel_widths = c(1,1))
dev.off()

#=================================================================================#
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#=================================================================================#

# summarize 
med.sum <- relative.mediation[!is.na(mediation.group),.(
  unique.probes = uniqueN(gx.probe),
  unique.genes = uniqueN(gene),
  unique.metabolites = uniqueN(metabolite)
),by=mediation.group]

# most frequently mediating genes
# ABCG1 -> Cholesterintransport

# Paper: list! ----

#' show strongest, most often mediating genes etc.

relative.mediation[pheno=="log.bmi" & !is.na(gene), .N]
relative.mediation[pheno=="log.bmi" & !is.na(gene), uniqueN(gene)]
relative.mediation[pheno=="log.bmi" & !is.na(gene), uniqueN(metabolite)]
relative.mediation[pheno=="log.bmi" & !is.na(gene) & mediation.group == "None", uniqueN(gene)]
relative.mediation[pheno=="log.bmi" & !is.na(gene) & mediation.group == "None", uniqueN(metabolite)]
relative.mediation[pheno=="log.bmi" & !is.na(gene) & mediation.group == "None", .N]
relative.mediation[pheno=="log.bmi" & !is.na(gene) & mediation.group != "None", .N]

relative.mediation[pheno=="log.bmi" & gene !="" & mediation.group != "None" & mediation.group != "Other", table(mediation.group)]


relative.mediation[!(mediation.group %in% c("None", "Other")), .N, by=gene][order(N,decreasing = T)] %>% head(10)
relative.mediation[!(mediation.group %in% c("None", "Other")) & gene == "ATM", unique(metabolite)]

# get most mediating genes by group
tmp1 <- dcast(relative.mediation[(`mediation.sig_pheno~gx(+metab)` == TRUE |
                                    `mediation.sig_pheno~metab(+gx)` == TRUE) &
                                   gene != "", 
                                 uniqueN(metabolite),by=.(pheno,gene,mediation.group)][
                                   pheno=="log.bmi" & mediation.group != "Other"][
                                     order(V1,decreasing=T)],
              pheno + gene ~ mediation.group,value.var =  "V1")

# replace na with 0 to be able to sum it all up

tmp1[is.na(`Gene expression mediated metabolite effects`), `Gene expression mediated metabolite effects` := 0]
tmp1[is.na(`Metabolite mediated gene expression effects`), `Metabolite mediated gene expression effects` := 0]
tmp1[,sig.sum := (`Gene expression mediated metabolite effects` + 
                    `Metabolite mediated gene expression effects`)]
tmp1[order(sig.sum,decreasing=T),] %>% head(10)

# general number of sig mediations
res.joint[is.mediation=="yes", .N,by=.(pheno, type)]
res.joint[is.mediation=="yes", sum(fdr.sig),by=.(pheno, type)]
res.joint[is.mediation=="yes", sum(fdr.sig)/.N,by=.(pheno, type)]

# occurences of mediating genes on both phenotypes
res.joint[is.mediation=="yes" & fdr.sig==T & pheno=="log.bmi", table(symbol_INGENUITY)] %>% sort(decreasing = T) %>% head(10)
res.joint[is.mediation=="yes" & fdr.sig==T & pheno=="diabetes.status.tri", table(symbol_INGENUITY)] %>% sort(decreasing = T) %>% head(10)

res.joint[is.mediation=="yes" & fdr.sig==TRUE, .(uniqueN(metabolite), uniqueN(symbol_INGENUITY))]
res.init.meta[significant==TRUE & symbol_INGENUITY != "",
              .(uniqueN(metabolite),
                unique(types_INGENUITY)),
              by=.(symbol_INGENUITY)][
                order(V1,decreasing = T)] %>% head(10)

# check genes in mediation result groups
relative.mediation$mediation.group %>% unique
venn3(
  relative.mediation[mediation.group == unique(mediation.group)[3], unique(gene)],
  relative.mediation[mediation.group == unique(mediation.group)[4], unique(gene)],
  relative.mediation[mediation.group == unique(mediation.group)[5], unique(gene)]
)

#=============================================================================#
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#=============================================================================#

devtools::session_info()