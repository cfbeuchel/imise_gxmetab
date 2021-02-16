source("functions/plot_multivar.R")
source("functions/generic_multiple_testing_correction.R")
library(ggplot2)
library(scales)
library(data.table)

res_sig <- fread("../../170829_ConfounderAnalysis/190207_ExtraAnalysis/results/190212_a1_b3_sorb_allPartialRsquaredWithoutQuotientsInclHighMissings.csv")

res_remaining  <- fread("../../170829_ConfounderAnalysis/190207_ExtraAnalysis/results/190212_a1_b3_sorb_partialRSquaredRemainingNotSelectedCovars.csv")

res_sig_m <- melt(res_sig, measure.vars = c("age", "diabetes.status.tri", "fasting.hours", "hematocrit", "log.bmi", "monocytes.percent", "neutrophils.percent", "sex"))

setnames(res_sig_m, "variable","covariate2")
setnames(res_sig_m, "value","term.r.squared")
res_sig_m[,r.squared.cutoff := 0.025]

res_sig_m$cohort <- factor(
    res_sig_m$cohort,
    levels = c("a1","b3","sorb"),
    labels = c("LIFE Adult", "LIFE Heart", "Sorb Study"))

res_sig_m$covariate <- factor(
    res_sig_m$covariate, 
    levels = unique(res_sig_m$covariate),
    labels = c(
    "Age",
    "Diabetes status (yes/no)",
    "Fasting hours (hours)",
    "Hematocrit (%)",
    "log-BMI (kg/m^2)",
    "Monocytes (%)",
    "Neutrophils (%)",
    "Sex (female/male)"
    ))

p1  <- plot_multivar(
    res_sig_m
) + ylab("Partial explained variance of factor")

tiff(
    "../../1805_gxMetaboliteAssociation/ppr/supp_fig/covarSelection.tiff",
    width = 7,
    height = 7, 
    res = 300,
    unit = "in",
)
p1 
dev.off()