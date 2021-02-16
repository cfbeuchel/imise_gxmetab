#' ---
#' title: "Statistics necessary for mediation analysis"
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
  "here",
  "magrittr",
  "ggplot2",
  "ggExtra",
  "parallel",
  "RMediation"
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
source(here("../functions/sobel_test.R"))
source(here("../functions/theme_gxMetab.R"))
an <- all_annotation_tables(mount.point = "/net/ifs2/san_projekte/projekte/genstat/")
theme_set(theme_light(base_size = 15, base_family = "Helvetica"))
setDTthreads(1)
option_setup()

#+ include=F
#==========#
#   DONE   #
#==========#

# Load data ----
#' # Load data

#' Load the meta-analysed effect estimates for each mediation-relevant association.

# these are the meta-analysed effect estimates for each mediation variant to be tested
nf <- newest_file("mediationMetaParts1","res",print_full = T)
res1 <- fread(nf)
nf <- newest_file("mediationMetaParts2","res",print_full = T)
res2 <- fread(nf)

# concatinate
res.meta <- rbindlist(list(res1,res2),use.names = TRUE,fill = TRUE )
rm(res1,res2) # clean up

# results from the main meta-analysis
nf <- newest_file(look_for = "STEP4_meta_resultsAnnotated$",subfolder = "res", print_full = TRUE)
res.init.meta <- fread(nf)

# cleaning
res.init.meta[metab.super.path =="",metab.super.path:="Other"]
res.init.meta[symbol_INGENUITY=="",symbol_INGENUITY:=NA]
res.init.meta[types_INGENUITY == "no_data" ,types_INGENUITY := "other"]
res.init.meta[location_INGENUITY == "no_data" ,location_INGENUITY := "Other"]

# save corrected meta analysis results
fwrite(res.init.meta, dpx("STEP4_meta_resultsAnnotated", "res/"), sep = "\t")

# these are the meta-analysed effect estimates for each mediation variant to be tested
nf <- newest_file("mediationParts1","res",print_full = T)
res1 <- fread(nf)
nf <- newest_file("mediationParts2","res",print_full = T)
res2 <- fread(nf)
res.single <- rbindlist(list(res1,res2),use.names = TRUE,fill = TRUE)

# remove all interaction associations
res.single <- res.single[x.type %nin% "interaction"]

rm(res1,res2) # clean up

# load the list with all relevant association tests
amt <- fread(newest_file("allMediationsList","res",print_full = TRUE))

# Checks ----
#' # Checks 

#' I need to check whether what I do was correct. Check the sample sizes of various mediation triangles, check the re-analyzed gx-metab associations, whether the results are the same.

# Identify the relevant statistics for a single mediation test

# identify a random assication to test
res.init.meta[hfdr.sigREM==TRUE & metabolite == "Lys" & markerID == "ILMN_1729980", 
         .(beta.Adult,`beta.Heart-AMI`,beta.Sorb,beta.Heart, betaFEM,betaREM,pREM, hfdr.sigREM, totalN)]

#' This is the starting point. Lys~ILMN_1729980 define an interesting pair for further investigation. Next, we check, wheter any of the two associate with any of my phenotypes.

# check for association with metabolite
res.meta[type=="pheno~gx" & x == "ILMN_1729980" & hfdr.sigREM==TRUE, .(type, y, x, m, pREM, totalN,hfdr.sigREM)]
res.meta[type=="pheno~metab" & x == "Lys" & hfdr.sigREM==TRUE, .(type, y, x, m, pREM, totalN,hfdr.sigREM)]

#' Both, metabolite and probe associate with the phenotype "diabetes.status.tri". Therefore, this association pair is considered clinically relevant. Subsequently, mediation effects are to be estimated for both directions (gx->metab->pheno and metab->gx->pheno). I have association statistics for both cases, each with x and m switched

#' First check whether the effect of the gene expression is mediated through the metabolite.

# Visualize R2 ----
#' # Visualize R2

res.single[type=="gx~metab", .SD[r.squared.1==max(r.squared.1)]]

p1 <- ggplot(res.single,aes(x=r.squared.1,
                            y=r.squared.2,
                            col=r.squared.type,
                            pch=study)
             )+
  facet_grid(r.squared.type~type) +
  geom_point() +
  ylim(c(0,1)) +
  xlim(c(0,1)) +
  geom_abline(slope = 1,intercept = 0) +
  theme_gxMetab()

tiff(dpx("r2Comparison.tiff","plt/"),res=300,unit="in",height=10,width=12)
p1
dev.off()

# Calculate Mediation ----
#' # Calculate Mediation

#' Calculate mediation effects for all associations.

# get index
amt[, index := 1:.N]

# how many tests are there?
nrow(amt)

# were checking mediation in 2 directions, so every combination is tested twice (gx->metab->pheno & metab->gx->pheno)
nrow(amt)*2

res.sobel <- mclapply(amt$index, function(i){
  
  # Exposure: Gx ----
  
  # get the relevant info on each association statistic
  exposure.type <- amt[index==(i),exposure]
  gx.probe <- amt[index==(i),gx.probe]
  metabolite <- amt[index==(i),metabolite]
  pheno <- amt[index==(i),pheno]
  
  if(exposure.type %in% c("both", "gx")) {
    
    # get the relevant statistics
    tau = res.meta[type == "pheno~gx" & y == (pheno) & x == (gx.probe), betaREM]
    
    tauAdj = res.meta[type == "pheno~gx+metab" & 
                        y == (pheno) & 
                        x == (gx.probe) & 
                        m == (metabolite), betaREM]
    alpha = res.meta[type=="metab~gx" & y == (metabolite) & x == (gx.probe), betaREM]
    seAlpha = res.meta[type=="metab~gx" & y == (metabolite) & x == (gx.probe), seREM]
    beta = res.meta[type == "pheno~gx+metab" & 
                      y == (pheno) & 
                      m == (gx.probe) & 
                      x == (metabolite), betaREM]
    seBeta = res.meta[type == "pheno~gx+metab" & 
                        y == (pheno) & 
                        m == (gx.probe) & 
                        x == (metabolite), seREM]
    
    # get unique sample sizes
    n1 = res.meta[type == "pheno~gx" & y == (pheno) & x == (gx.probe), totalN]
    n2 = res.meta[type == "pheno~gx+metab" & 
                    y == (pheno) & 
                    x == (gx.probe) & 
                    m == (metabolite), totalN]
    n3 = res.meta[type=="metab~gx" & y == (metabolite) & x == (gx.probe), totalN]
    n4 = res.meta[type == "pheno~gx+metab" & 
                    y == (pheno) & 
                    m == (gx.probe) & 
                    x == (metabolite), totalN]
    
    # calcuate the statistics
    res1 <- sobel_test(
      tau = tau,
      tauAdj = tauAdj,
      alpha = alpha,
      seAlpha = seAlpha,
      beta = beta,
      seBeta = seBeta)
    
    res1[,`:=`(
      tau = tau,
      tauAdj = tauAdj,
      alpha = alpha,
      seAlpha = seAlpha,
      beta = beta, 
      seBeta = seBeta,
      y = pheno,
      x = gx.probe,
      x.type = "gx.probe",
      m = metabolite,
      m.type = "metabolite",
      exposure = exposure.type,
      n = list(c("tau" = n1, 
                 "tauAdj" = n2,
                 "alpha" = n3, 
                 "beta" = n4))
    )]
    
  }
  
  # Exposure: Metab ----
  
  if(exposure.type %in% c("both", "metab")){

    # calculate the mediation effect in the other direction
    # get the relevant statistics
    tau = res.meta[type == "pheno~metab" & y == (pheno) & x == (metabolite), betaREM]
    tauAdj = res.meta[type == "pheno~gx+metab" & 
                        y == (pheno) & 
                        x == (metabolite) & 
                        m == (gx.probe), betaREM]
    alpha = res.meta[type=="gx~metab" & y == (gx.probe) & x == (metabolite), betaREM]
    seAlpha = res.meta[type=="gx~metab" & y == (gx.probe) & x == (metabolite), seREM]
    beta = res.meta[type == "pheno~gx+metab" & 
                      y == (pheno) & 
                      x == (gx.probe) & 
                      m == (metabolite), betaREM]
    seBeta = res.meta[type == "pheno~gx+metab" & 
                        y == (pheno) & 
                        x == (gx.probe) & 
                        m == (metabolite), seREM]
    
    # get the sample sizes
    n1 = res.meta[type == "pheno~metab" & y == (pheno) & x == (metabolite), totalN]
    n2 = res.meta[type == "pheno~gx+metab" & 
                    y == (pheno) & 
                    x == (metabolite) & 
                    m == (gx.probe), totalN]
    n3 = res.meta[type=="gx~metab" & y == (gx.probe) & x == (metabolite), totalN]
    n4 = res.meta[type == "pheno~gx+metab" & 
                    y == (pheno) & 
                    x == (gx.probe) & 
                    m == (metabolite), totalN]
    
    # calcuate the statistics
    res2 <- sobel_test(
      tau = tau,
      tauAdj = tauAdj,
      alpha = alpha,
      seAlpha = seAlpha,
      beta = beta,
      seBeta = seBeta)
    
    res2[,`:=`(
      tau = tau,
      tauAdj = tauAdj,
      alpha = alpha,
      seAlpha = seAlpha,
      beta = beta, 
      seBeta = seBeta,
      y = pheno,
      x = metabolite,
      x.type = "metabolite",
      m = gx.probe,
      m.type = "gx.probe",
      exposure = exposure.type,
      n = list(c("tau" = n1, 
                 "tauAdj" = n2, 
                 "alpha" = n3, 
                 "beta" = n4))
    )]
    
  }
  
  # return results based on what was computed
  if(exposure.type == "gx"){
    
    return(res1)
    
  }else if(exposure.type == "metab"){
    
    return(res2)
    
  }else if(exposure.type == "both"){
    
    res <- rbindlist(list(res1,res2))
    return(res)
    
  } else {
    stop("ERROR! Unclear 'exposure.type'")
  }
  
},mc.cores = 20,mc.cleanup = TRUE) %>% rbindlist()

# Testing
range(res.sobel$p.value)
range(res.sobel$p.value.prod)
par(mfrow=c(2,1))
hist(res.sobel$p.value, main=NULL)
hist(res.sobel$p.value.prod, main=NULL)
par(mfrow=c(1,1))
res.sobel[,plot(p.value~p.value.prod,col=alpha("grey40",5),pch=20)]
abline(0,1,col="indianred4",lwd=1.2,lty="dashed")

# save results
fwrite(res.sobel, dpx("mediationStatistics",folder = "res/"),quote = TRUE, sep="\t")
# res.sobel <- fread(newest_file("mediationStatistics$","res",print_full = T))

# Annotate Mediation ----
#' # Annotate Mediation

# what is being tested
res.sobel[,.N,by=x.type]
res.sobel[,.N,by=.(y,x.type)]

# get summary info on metabolite
res.sobel[, metabolite := ifelse(x.type=="metabolite",x,m)]
res.sobel[, gx.probe := ifelse(x.type=="gx.probe",x,m)]

# mediation through metab->gx
sobel.fdr <- res.sobel[,addHierarchFDR(
  pvalues = p.value,
  categs = x,
  fdrmethod_level1 = "BH",
  fdrmethod_level2 = "BH",
  correctionLevel1 = "BB")]

# quick summary
sobel.fdr[,table(hierarch_fdr5proz)]
sobel.fdr[hierarch_fdr5proz==TRUE,.(sig.medi = .N),by=.(category)][
  order(sig.medi,decreasing = TRUE)] %>% head(10)

# merge to results
res.sobel[, `:=`(
  fdr1 = sobel.fdr$fdr_level1,
  fdr2 = sobel.fdr$fdr_level2,
  fdr.sig = sobel.fdr$hierarch_fdr5proz
)]

# match annotation info
m1 <- match(res.sobel$gx.probe, res.init.meta$markerID)
res.sobel[, c("symbol_INGENUITY",
              "location_INGENUITY",
              "types_INGENUITY",
              "biomarkers_INGENUITY",
              "drugs_INGENUITY") := res.init.meta[m1,.SD,.SDcols=c(
                "symbol_INGENUITY",
                "location_INGENUITY",
                "types_INGENUITY",
                "biomarkers_INGENUITY",
                "drugs_INGENUITY")]
          ]

# get the type to indicate the mediation direction 
res.sobel[, type := ifelse(x.type=="gx.probe", "pheno~gx(+metab)", "pheno~metab(+gx)")]

# remove empty gene description
res.sobel[, sum(symbol_INGENUITY == "",na.rm = T)]
res.sobel[symbol_INGENUITY == "", symbol_INGENUITY := NA]

# get a col with gene names instead of probe IDs
res.sobel[,x.gene := x]
res.sobel[x.type=="gx.probe",x.gene := symbol_INGENUITY]
res.sobel[,m.gene := m]
res.sobel[x.type=="metabolite",m.gene := symbol_INGENUITY]

# Check - Plot the different mediation effect estimates tau-tau' and alpha*beta

ggplot(res.sobel, aes(beta.tau.tauAdj, beta.alpha.beta, col=fdr.sig)) +
  geom_point() +
  geom_abline(slope = 1,intercept = 0, col="indianred4",lwd=1.2) +
  facet_wrap(y~x.type, scales="free") +
  theme_gxMetab()

# check total effect estimate tau vs alpha*beta + tauAdj
res.sobel[fdr.sig==T, plot(y = tau, x = alpha * beta + tauAdj, col=alpha("grey35",.8),pch=16)]
abline(0,1, col="indianred4", lwd=1.2)

# save results
fwrite(res.sobel, dpx("mediationStatisticsAnnotated",folder = "res/"),quote = TRUE, sep="\t")
# res.sobel <- fread(newest_file("mediationStatisticsAnnotated","res",print_full = T))

#=============================================================================#
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#=============================================================================#

# Check Sobel Results ----
#' # Check Sobel Results

# Create table to compare results with single study results ----
#' ## Create table to compare results with single study results

res.single[type=="pheno~gx", x.type := "gx.probe"]
res.single[type=="pheno~metab", x.type := "metabolite"]
res.single[type=="gx~metab", x.type:= "metabolite"]
res.single[type=="metab~gx", x.type:= "gx.probe"]
res.single[type=="pheno~gx+Metab" & x.type == "metabolite", type := "pheno~metab(+gx)"]
res.single[type=="pheno~gx+Metab" & x.type == "gx.probe", type := "pheno~gx(+metab)"]
res.single[x.type=="metabolite", `:=`(metabolite = x)]
res.single[x.type=="gx.probe", `:=`(gx.probe = x)]
res.single[type=="gx~metab", `:=`(gx.probe = y, metabolite = x)]
res.single[type=="metab~gx", `:=`(gx.probe = x, metabolite = y)]
res.single[type=="pheno~metab(+gx)", gx.probe := m]
res.single[type=="pheno~metab(+gx)", pheno := y]
res.single[type=="pheno~gx(+metab)", metabolite := m]
res.single[type=="pheno~gx(+metab)", pheno := y]
res.single[type=="pheno~gx" | type == "pheno~metab", pheno := y]

# check
all.types <- res.single$type %>% unique
all.types
res.single[type==all.types[1]] %>% hh
res.single[type==all.types[2]] %>% hh
res.single[type==all.types[3]] %>% hh
res.single[type==all.types[4]] %>% hh
res.single[type==all.types[5]] %>% hh
res.single[type==all.types[6]] %>% hh

# cast
#[type == "pheno~gx(+metab)" | 
# type == "pheno~metab(+gx)"]
res.single.cast <- dcast(res.single, 
                    pheno + metabolite + gx.probe + type ~ study, 
                    value.var = c("beta",
                                  "beta.se", 
                                  "test.statistic", 
                                  "p.value", 
                                  "r.squared.1",
                                  "r.squared.2",
                                  "r.squared.3"
                                  ),
                    sep = ".")

#=============================================================================#
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#=============================================================================#

# Create table to compare results meta analyzed study results ----
#' ## Create table to compare results with meta analyzed study results

res.meta[type=="pheno~gx", x.type := "gx.probe"]
res.meta[type=="pheno~metab", x.type := "metabolite"]
res.meta[type=="gx~metab", x.type:= "metabolite"]
res.meta[type=="metab~gx", x.type:= "gx.probe"]
res.meta[type=="pheno~gx+metab" & x.type == "metabolite", type := "pheno~metab(+gx)"]
res.meta[type=="pheno~gx+metab" & x.type == "gx.probe", type := "pheno~gx(+metab)"]
res.meta[x.type=="metabolite", `:=`(metabolite = x)]
res.meta[x.type=="gx.probe", `:=`(gx.probe = x)]
res.meta[type=="gx~metab", `:=`(gx.probe = y, metabolite = x)]
res.meta[type=="metab~gx", `:=`(gx.probe = x, metabolite = y)]
res.meta[type=="pheno~metab(+gx)", gx.probe := m]
res.meta[type=="pheno~metab(+gx)", pheno := y]
res.meta[type=="pheno~gx(+metab)", metabolite := m]
res.meta[type=="pheno~gx(+metab)", pheno := y]
res.meta[type=="pheno~gx" | type == "pheno~metab", pheno := y]

# get the explained variance for each metab-gx mediation pair

# Calculate R2 - Arnd Formel ----

res.meta[, r.squared.REM :=  betaREM^2 / ( betaREM^2 + totalN * seREM^2 )]

res.meta[type=="pheno~gx(+metab)" | type=="pheno~metab(+gx)" , 
    total.r.squared.REM := sum(r.squared.REM), by = .(id.name, y)]
res.meta.all <- res.meta[,.(pheno, 
              metabolite, 
              gx.probe,
              type, 
              betaREM,
              seREM,
              I2,
              n.REM = totalN,
              numberStudies, 
              hfdr.1REM,
              hfdr.2REM,
              hfdr.sigREM, 
              r.squared.REM,
              total.r.squared.REM
              )]

# merge res.single.cast and res.sobel
v2 <- venn3(
  res.single.cast[, paste(type,pheno,gx.probe,metabolite)],
  
  # edit 200902: changed y to pheno, dont know why I had to..
  res.sobel[, paste(type,y,gx.probe,metabolite)],
  # res.sobel[, paste(type,pheno,gx.probe,metabolite)],
  res.meta.all[, paste(type,pheno,gx.probe,metabolite)])

#=============================================================================#
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#=============================================================================#

# Create Joint Result table ----
#' # Create Joint Result table 

# join the results table for easier comparison
mid <- res.single.cast[, paste(type, pheno,gx.probe,metabolite)]
res.joint <- cbind(res.single.cast,
      res.sobel[match(mid,paste(type,y,gx.probe,metabolite))],
      # 200902 - changes: res.sobel[match(mid,paste(type,pheno,gx.probe,metabolite))],
      res.meta.all[match(mid, paste(type,pheno,gx.probe,metabolite))]
      )

# check and remove double gx.probe and metabolite column
names(res.joint)[names(res.joint)%in%c(
  "gx.probe",
  "metabolite",
  "type")] <- c("metabolite", 
                "gx.probe",
                "type",
                "metabolite2",
                "gx.probe2",
                "type2",
                "metabolite3",
                "gx.probe3",
                "type3")

# remove unnecessary cols
res.joint$type2 <- NULL
res.joint$metabolite2 <- NULL
res.joint$gx.probe2 <- NULL
res.joint$type3 <- NULL
res.joint$metabolite3 <- NULL
res.joint$gx.probe3 <- NULL
res.joint$y <- NULL
res.joint$x <- NULL
res.joint$m <- NULL
res.joint[, which(duplicated(names(res.joint))) := NULL]
# set(res.joint, j = names(res.joint)[duplicated(names(res.joint))], value = NULL)
res.joint$n <- NULL


# check for duplicated colnames
stopifnot(!any(duplicated(names(res.joint))))

# to avoid filtering
res.joint[is.na(pheno),pheno:="none"]
res.joint[,is.mediation:=ifelse(type=="pheno~gx(+metab)" | type=="pheno~metab(+gx)", "yes", "no")]
res.joint$is.mediation %>% mytable

# save results
fwrite(res.joint, 
       dpx("mediationStatisticsAnnotatedJoint",folder = "res/"), 
       quote = TRUE, 
       sep="\t")
# res.joint <- fread(newest_file("mediationStatisticsAnnotatedJoint","res",print_full = T))

#=============================================================================#
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#=============================================================================#

# Summaries ----
#' # Summaries 

#' Following are basic counts of significant mediations per type etc.

# number of probes and metabolites in the sobel test
amt$metabolite %>% uniqueN
amt$gx.probe %>% uniqueN

# significancies for each x.type and response
res.sobel[x.type=="gx.probe",.(gx = uniqueN(x),
                               metab = uniqueN(m)
)]

# number of signifcant mediations
res.sobel[fdr.sig==TRUE,.N]
res.sobel[,.N,by=x.type]
res.sobel[fdr.sig==TRUE,.N,by=.(x.type)]
res.sobel[,sum(fdr.sig==TRUE)/.N,by=.(x.type)]
res.sobel[fdr.sig==TRUE,.N,by=.(x.type,y)][order(x.type)]
res.sobel[,.(SIG = sum(fdr.sig==TRUE),
             N=.N, 
             REL = sum(fdr.sig==TRUE)/.N),
          by=.(x.type,y)][order(x.type)]

# number of probes and metabolites
res.sobel[fdr.sig==TRUE,
          .(metab=uniqueN(metabolite),
            gene=uniqueN(symbol_INGENUITY),
            probes=uniqueN(gx.probe))] # minus empty genes
res.sobel[,
          .(metab=uniqueN(metabolite),
            gene=uniqueN(symbol_INGENUITY),
            probes=uniqueN(gx.probe))] # minus empty genes
res.sobel[fdr.sig==TRUE,.(metab=uniqueN(metabolite),
                          gene=uniqueN(gx.probe)),by=x.type]
res.sobel[fdr.sig==TRUE,.(metab=uniqueN(metabolite),
                          gene=uniqueN(gx.probe)),by=y]
res.sobel[fdr.sig==TRUE,.(metab=uniqueN(metabolite),
                          gene=uniqueN(gx.probe)),by=.(x.type,y)]

# Compare
res.sobel[, id:=paste(y,gx.probe,metabolite,sep="__")]
res.sobel[x.type=="gx.probe",.(.N,uniqueN(id))]
res.sobel[x.type=="gx.probe" & fdr.sig==TRUE,.(.N,uniqueN(id))]
res.sobel[x.type=="metabolite",.(.N,uniqueN(id))]
res.sobel[x.type=="metabolite" & fdr.sig==TRUE,.(.N,uniqueN(id))]

# significant mediations in both directions
v.type1 <- venn2(
  res.sobel[x.type=="gx.probe" & fdr.sig==TRUE,id],
  res.sobel[x.type=="metabolite" & fdr.sig==TRUE,id]
)
res.sobel[fdr.sig==TRUE, uniqueN(id)]
res.sobel[fdr.sig==TRUE, length(v.type1$q1)/uniqueN(id)]

# significant bi-directional mediations for bmi and diabets 
v.type2 <- venn2(
  res.sobel[y == "log.bmi" & x.type=="gx.probe" & fdr.sig==TRUE,id],
  res.sobel[y == "log.bmi" & x.type=="metabolite" & fdr.sig==TRUE,id]
)
v.type3 <- venn2(
  res.sobel[y == "diabetes.status.tri" & x.type=="gx.probe" & fdr.sig==TRUE,id],
  res.sobel[y == "diabetes.status.tri" & x.type=="metabolite" & fdr.sig==TRUE,id]
)
res.sobel[fdr.sig==TRUE, uniqueN(id),by=y]
res.sobel[fdr.sig==TRUE, length(v.type2$q1)/uniqueN(id)]
res.sobel[fdr.sig==TRUE, length(v.type3$q1)/uniqueN(id)]

#=============================================================================#
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#=============================================================================#

sessionInfo()
