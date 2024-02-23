

library(phytools); library(geiger); library(ape)

mtrees_M <- readRDS("~/GitHub/Evo-Legumes-Herbivory/LegumesArmature/output/RDS/Mammals_4traits_simmap_100.rds")

### Modified function to plot symbols instead of lines =============
markChanges1<-function(mtree,colors=cols1,cex=1,lwd=2){
  states<-sort(unique(getStates(mtree)))
  if(is.null(colors)) colors<-setNames(palette()[1:length(states)],
                                       states)
  obj<-get("last_plot.phylo",envir=.PlotPhyloEnv)
  nc<-sapply(mtree$maps,length)-1
  ii<-which(nc>0)
  nc<-nc[ii]
  h<-vector()
  for(i in 1:length(ii)){
    for(j in 1:nc[i]){
      ss<-names(tree$maps[[ii[i]]])[j+1]
      mm<-tree$edge[ii[i],1]
      dd<-tree$edge[ii[i],2]
      x<-rep(obj$xx[mm]+cumsum(tree$maps[[ii[i]]])[j],2)
      y<-c(obj$yy[dd]-0.000001*mean(strheight(LETTERS)*cex),
           obj$yy[dd]+0.000001*mean(strheight(LETTERS)*cex))
      #lines(x,y,lwd=5, col=cols1[[ss]], fill=cols2[[ss]],lend=2)
      points(x,y, pch = syms[ss], col=cols1[[ss]], bg=cols1[ss], cex=0.8)
      #h<-c(h,x[1])
    }
  }
  #invisible(h)
}
#### ===========================================

# Color Palette ========================================================
states <- sort(unique(getStates(mtrees_M[[1]], "tips")))
cols1<-setNames(c("#476653", "#87b09a", "darkgrey", "lightgrey"), states)
syms <- setNames(c(24, 23, NA, NA), states)   # Put NA where the symbol should be removed
### Check this link for the symbols: 
### http://www.sthda.com/english/wiki/r-plot-pch-symbols-the-different-point-shapes-available-in-r

## Large Plot ===========================================================
pdf(file="output/Figs/Phylo_states_mammals_large.pdf", heigh = 80, width = 8)
plotSimmap(mtrees_M[[1]],ftype="i",fsize=0.1,ylim=c(-1,Ntip(mtrees_M[[1]])), col=cols1)
obj <- markChanges1(mtrees_M[[1]])
add.simmap.legend(x=0,y=-1,colors=cols1,prompt=FALSE,vertical=FALSE)
points(obj,pch=8)
dev.off()
