
#' Hi! THis script is basically an archive. I might want to do some dendrogram related stuff in the future, that's what this script is for!

########################################################
############## NOT USED RIGHT NOW ######################

#' ## Draw dendrogram

#' Maybe draw a dendrogram of both correlation of the phenotypes and beta correlations of the meta-analysis and compare the clustering

#' Info here: https://cran.r-project.org/web/packages/dendextend/vignettes/introduction.html

dend <- dist(mat,method="euclidean") %>% 
  hclust(method = "average") %>% 
  as.dendrogram()
order.dendrogram(dend)

dend_height = attr(dend, "height")

labels(dend) %>% sort

labels(reorder(dend, x))

circos.trackPlotRegion(
  ylim = c(0, dend_height), 
  bg.border = NA, 
  track.index=4,
  track.height = 0.2,
  panel.fun = function(x, y) {
    
    # draw circ
    circos.dendrogram(
      dend = dend,
      facing = "outside",
      max_height = dend_height
    )
  })

# tracks are the rings ordered from the most outer circle (1) to the most inner circle 1+n
circos.trackPlotRegion(
  ylim = c(0, 1),
  track.index = 3,
  bg.border = NA, 
  track.height = 0.1,
  panel.fun = function(x, y) {
    
    # this only gives the *sector* name (in my case, AA, AC)
    nm = get.cell.meta.data("sector.index")
    
    # this returns all the regions of one sector
    r = gl[[nm]]
    
    # this returns the number of regions/indices per sector
    n = length(r)
    
    # draw the content of each sector (i.e. the metabolites)
    circos.text(x = 1:n - 0.5, 
                y = rep(0.5, n), 
                labels = r, 
                facing = "clockwise", 
                niceFacing = TRUE, 
                cex = 0.5)
    
  })


#' ### Phylogram of effect-size clustering

# use dendextend with circlize
# https://cran.r-project.org/web/packages/dendextend/vignettes/introduction.html

#' calculate distance for each metabolite combination then do hierarchical clustering and create the dendrogram object

# turn clustering object into a dendrogram
dend <- dist(mat,method="euclidean") %>% 
  hclust(method = "average") %>% 
  as.dendrogram()
dend %>% highlight_branches_col %>% plot()

# get colored bar based on class
bars <- rep(NA,length(labels(dend)))
names(bars) <- labels(dend)
bars[names(bars) %in% group[metab.class.q=="aa", V1]] <- group_color[names(group_color)=="aa"]
bars[names(bars) %in% group[metab.class.q=="ac", V1]] <- group_color[names(group_color)=="ac"]
bars[names(bars) %in% group[metab.class.q=="quotient", V1]] <- group_color[names(group_color)=="quotient"]
bars[names(bars) %in% group[metab.class.q=="sum", V1]] <- group_color[names(group_color)=="sum"]

# plot with bars indicating class
par.backup <- par(no.readonly = TRUE)
par(mar = c(10,2,1,1))
plot(dend)
colored_bars(colors = bars, dend = dend, rowLabels = "Class",text_shift = -.5)
par(par.backup)

# plot as cicular plot
par(mar = rep(0,4))
circlize_dendrogram(dend)