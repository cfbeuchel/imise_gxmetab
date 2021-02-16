require(data.table)


calcR2 = function(estimate, std.error, FallzahlODERFreiheitsgrad) {
  # message("Calculating explained variance accounting for covariates according to PMID:25898129")
  estimate^2 / (estimate^2 + FallzahlODERFreiheitsgrad *(std.error^2))
  
}

getR2fromMod = function(mod_x1, var = 'x1', n = 1000, modname = "") {
  res = c()
  res$var = var
  res$modname = modname
  res$estimate = mod_x1$coefficients[var,"Estimate"]
  res$se = mod_x1$coefficients[var,"Std. Error"]
  res$pval = mod_x1$coefficients[var,"Pr(>|t|)"]
  
  res$r2 = calcR2(estimate = res$estimate, std.error = res$se, FallzahlODERFreiheitsgrad = n)
  as.data.table(res)
}



compareR2 = function(mod_x1, mod_x2, mod_x1x2) {
  r2_modx1_x1=  getR2fromMod(mod_x1, var= "x1", modname = "modx1_x1")
  r2_modx2_x2=  getR2fromMod(mod_x2, var= "x2", modname = "modx2_x2")
  r2_modx1x2_x1=  getR2fromMod(mod_x1x2, var= "x1", modname = "modx1x2_x1")
  r2_modx1x2_x2=  getR2fromMod(mod_x1x2, var= "x2", modname = "modx1x2_x2")
  allres = rbindlist(list(r2_modx1_x1, r2_modx2_x2, r2_modx1x2_x1, r2_modx1x2_x2))  
  allres2 = dcast.data.table(allres, .~modname, value.var = c("estimate", "se", "r2", "pval"))
  allres2[,diff_x1 := r2_modx1_x1 -r2_modx1x2_x1]
  allres2[,diff_x2 := r2_modx2_x2 -r2_modx1x2_x2]
  allres2
}


set.seed(1902)

simu = 1:1000
allsim = lapply(simu, function (mysimu) {
  # mysimu=1
  # message(mysimu)
  x1 = rnorm(1000)
  x2 = rnorm(1000) + 0.1*x1
  y = rnorm(1000) +0.1*x1+0.2*x2
  
  mod_x1 = summary(lm(y~x1))
  mod_x2 = summary(lm(y~x2))
  mod_x1x2 = summary(lm(y~x1+x2))
  
  erg = compareR2(mod_x1, mod_x2, mod_x1x2)
  erg$simu = mysimu
  erg
})




allsim_noassoc = lapply(simu, function (mysimu) {
  # mysimu=1
  # message(mysimu)
  x1 = rnorm(1000)
  x2 = rnorm(1000) + 0.1*x1
  y = rnorm(1000) +0.2*x2
  
  mod_x1 = summary(lm(y~x1))
  mod_x2 = summary(lm(y~x2))
  mod_x1x2 = summary(lm(y~x1+x2))
  
  erg = compareR2(mod_x1, mod_x2, mod_x1x2)
  erg$cor_x1x2 = cor(x1, x2)
  erg$cor_x1x2_pval = cor.test(x1, x2)$p.value
  erg$simu = mysimu
  erg
})
allsim_noassoc2 = rbindlist(allsim_noassoc)


# tiff(dpx("simulationR2.tiff","plt/"),res=300,unit="in",height=10,width=10)

par(mfrow = c(2,2))

allsim2 = rbindlist(allsim)
allsim2[, plot(r2_modx1x2_x1 ~r2_modx1_x1, main = "x1 and x2 associated with y")]
abline(0,1, col = "red")
allsim2[, plot(r2_modx1x2_x2 ~r2_modx2_x2)]
abline(0,1, col = "red")

allsim_noassoc2[, plot(r2_modx1x2_x1 ~r2_modx1_x1, main = "x2, but not x1 associated with y")]
abline(0,1, col = "red")
allsim_noassoc2[, plot(r2_modx1x2_x2 ~r2_modx2_x2)]
abline(0,1, col = "red")

# dev.off()

tiff(dpx("simulationR2DiffHist.tiff","plt/"),res=300,unit="in",height=10,width=10)

par(mfrow = c(2,2))

hist(allsim2$diff_x1, main = "x1 and x2 associated with y")
abline(v = 0,lty="dashed",col="indianred4",lwd=2)
hist(allsim2$diff_x2)
abline(v = 0,lty="dashed",col="indianred4",lwd=2)
hist(allsim_noassoc2$diff_x1, main = "x2, but not x1 associated with y")
abline(v = 0,lty="dashed",col="indianred4",lwd=2)
hist(allsim_noassoc2$diff_x2)
abline(v = 0,lty="dashed",col="indianred4",lwd=2)

dev.off()



### Diff 

plan <- data.table(
  x1.in.x2 = c(0.5,0.5,.5,.5, 0),
  x1.in.y  = c(0.5,0,.5,0, 0.5),
  x2.in.y  = c(0,0.5,.5,0, 0.5)
)

res <- list()

for(i in 1:nrow(plan)){
  
  allsim = lapply(simu, function (mysimu) {
    
    # mysimu=1
    # message(mysimu)
    
    # if(i==1){
    #   # S1 scheinkorrelation - 
    #   x1 = rnorm(1000)
    #   x2 = rnorm(1000) + 0.5*x1 # gx macht metab
    #   y = rnorm(1000) +0.5*x1+0.0*x2 # nur gx macht bmi
    # }
    # 
    # if(i==2){
    #   # S2
    #   x1 = rnorm(1000)
    #   x2 = rnorm(1000) + 0.5*x1 # gx macht metab
    #   y = rnorm(1000) +0.0*x1+0.5*x2 # metab macht bmi
    # }
    # 
    # if(i==3){
    #   # S3 - bei metab bleibt mehr übrig 
    #   x1 = rnorm(1000)
    #   x2 = rnorm(1000) + 0.5*x1 # gx macht metab
    #   y = rnorm(1000) +0.5*x1+0.5*x2 # gx & metab machen bmi
    # }
    # 
    # if(i==4){
    #   # S4 - no association between y and x1|x2
    #   x1 = rnorm(1000)
    #   x2 = rnorm(1000) + 0.5*x1 #  metab macht gx
    #   y = rnorm(1000) + 0.0*x1+0.0*x2 # metab macht bmi
    # }
    # 
    # if(i==5){
    #   # S5 
    #   x1 = rnorm(1000)
    #   x2 = rnorm(1000) + 0*x1 
    #   y = rnorm(1000) +0.5*x1+0.5*x2 # gx & metab machen bmi
    # }
    # 

    x1.in.x2 <- plan[(i), x1.in.x2]
    x1.in.y  <- plan[(i), x1.in.y]
    x2.in.y <- plan[(i), x2.in.y]
    
    x1 = rnorm(1000) # gx
    x2 = rnorm(1000) + x1.in.x2*x1 # metab
    y = rnorm(1000) + x1.in.y*x1+x2.in.y*x2 # bmi
    
    mod_x1 = summary(lm(y~x1))
    mod_x2 = summary(lm(y~x2))
    mod_x1x2 = summary(lm(y~x1+x2))
    
    erg = compareR2(mod_x1, mod_x2, mod_x1x2)
    erg$simu = mysimu
    erg
    
  }) %>% rbindlist
  
  allsim[, index := (i)]
  res[[i]] <- allsim
  
}


res <- rbindlist(res)
res[,`:=`(rel.diffx1 = diff_x1/r2_modx1_x1,
          rel.diffx2 = diff_x2/r2_modx2_x2,
          gx.proz.dir.of.marg = r2_modx1x2_x1 / r2_modx1_x1,
          metab.proz.dir.of.marg = r2_modx1x2_x2 / r2_modx2_x2)]

# 
p <- ggplot(res,
            aes(gx.proz.dir.of.marg,
                metab.proz.dir.of.marg,
                fill=as.factor(index))) +
  geom_point(alpha=0.5, pch=21,col="grey35") + 
  geom_hline(yintercept = 1,lwd=0.6,col="indianred4",lty="dashed") +
  geom_vline(xintercept = 1,lwd=0.6,col="indianred4",lty="dashed") +
  scale_x_continuous(limits = c(.0,2)) +
  scale_y_continuous(limits = c(.0,2)) +
  labs(
    x=expression("direct-r"[Gx]^2 %/% "raw-r"[Gx]^2),
    y=expression("direct-r"[Metabolite]^2 %/% "raw-r"[Metabolite]^2)
    ) +
  theme_gxMetab() +
  theme(legend.position = "none")

# 
p


tiff(dpx("simulationMediation.tiff","plt/"),res=300,unit="in",height=8,width=8)
print(p)
dev.off()

allsim[,.(r2_modx1_x1,gx.proz.dir.of.marg, r2_modx2_x2,metab.proz.dir.of.marg, r2_modx1x2_x1,  r2_modx1x2_x2)]

toolboxH::finalizeSkript()
