library(dplyr)
library(rstatix)
library(phytools)
### Functions =================================================================

getChanges<-function(tree){
  states<-sort(unique(getStates(tree)))
  nc<-sapply(tree$maps,length)-1
  b<-which(nc>0)
  nc<-nc[b]
  xx<-vector()
  H<-nodeHeights(tree)
  for(i in 1:length(b)){
    for(j in 1:nc[i]){
      ss<-names(tree$maps[[b[i]]])[j+1]
      x<-rep(H[b[i],1]+cumsum(tree$maps[[b[i]]])[j],2)
      xx<-c(xx,setNames(x[1],
                        paste(names(tree$maps[[b[i]]])[j:(j+1)],
                              collapse="->")))
    }
  }
  xx
}



#### Data =====================================================================
changes_table_Mammals <- changes_table_wo_edge
changes_table_Mammals_full <- changes_table_w_edge
#saveRDS(changes_table_Mammals_full, "~/GitHub/Evo-Legumes-Herbivory/LegumesArmature/output/RDS/changes_table_mammals_with_number_of_species.rds")
#saveRDS(changes_table_Mammals, "~/GitHub/Evo-Legumes-Herbivory/LegumesArmature/output/RDS/changes_table_mammals.rds")


changes_table_Legumes <- changes_table_wo_edge
#saveRDS(changes_table_Legumes, "~/GitHub/Evo-Legumes-Herbivory/LegumesArmature/output/RDS/changes_table_legumes.rds")


#### Legumes ==================================================================
mtrees_L <- readRDS("~/GitHub/Evo-Legumes-Herbivory/LegumesArmature/output/RDS/Legumes_mtrees_multiSimmap.rds")
changes_table_Legumes <- readRDS("~/GitHub/Evo-Legumes-Herbivory/LegumesArmature/output/RDS/changes_table_legumes.rds")
tree <- readRDS("output/RDS/Legumes_MCC.rds")


#### Subset to test script::
set.seed(123)
species_samp <- sample(tree$tip.label, 30)

tree <- keep.tip(tree, species_samp)
mtrees_L <-keep.tip(mtrees_L, species_samp)

######


plotTree(tree,type="fan",ftype="i",lwd=1, fsize = 0.2)
nodelabels(cex=0.5)


as.data.frame(ltt_L$times)


# Transform the original changes list of matrices from sapply() for all simulations to dataframe 
# add trans_no to each transition based on grouping by transition 
# (i.e, there are two trans_no = 1, one for each transition type)
changes<-sapply(mtrees_L,getChanges)
i_changes_list <- list() 
for (i in 1:length(changes)) {
  
  i_changes <- changes[[i]] # for each simulation from multiSimmap
  
  
  # extract transition types and times and assign ID column for simmap simulation
  i_changes_df <- data.frame(transition = names(i_changes), time = i_changes, N_sim = i)
  i_changes_df$transition <- as.factor(i_changes_df$transition)
  
  # Add trans_no to each transition from a group (transition type)
  i_changes_df <- i_changes_df %>%  arrange(time) %>%
    group_by(transition) %>%
    mutate(trans_no = row_number())
  
  i_changes_list[[i]] <- i_changes_df
  
}

# Final Table without edge.lengths:
changes_table <- plyr::rbind.fill(i_changes_list, fill=T)
changes_table_wo_edge <- changes_table



h<-max(nodeHeights(tree))
b<-20
segs<-cbind(seq(0,h-h/b,h/b),
            seq(1/b*h,h,h/b))
nchanges<-rep(0,b)
for(i in 1:length(changes)){
  for(j in 1:length(changes[[i]])){
    ind<-which((changes[[i]][j]>segs[,1])+
                 (changes[[i]][j]<=segs[,2])==2)
    nchanges[ind]<-nchanges[ind]+1/length(changes)
  }
}


# Estimate number of lineages through time ===================================
ltt_L <- phytools::ltt(tree)
# transform LTT data so that we can use it:
LTT <- cbind(ltt_L$ltt[2:(length(ltt_L$ltt)-1)],            # Number of species at each node minus the root
             ltt_L$times[2:(length(ltt_L$ltt)-1)],          # "start" time for each node minus the root
             ltt_L$times[3:length(ltt_L$ltt)])              # "end" time for each node minus the root

tail(LTT)
head(LTT)

## Edge Lengths (by which nchanges will be divided)  =========================
ii<-1
edge.length<-rep(0,b)
for(i in 1:nrow(segs)){
  done.seg<-FALSE
  while(LTT[ii,2]<=segs[i,2]&&done.seg==FALSE){
    edge.length[i]<-edge.length[i]+
      LTT[ii,1]*(min(segs[i,2],LTT[ii,3])-
                   max(segs[i,1],LTT[ii,2]))
    if(LTT[ii,3]>=segs[i,2]) done.seg<-TRUE
    if(LTT[ii,3]<=segs[i,2]) ii<-if(ii<nrow(LTT)) ii+1 else ii
  }
}

# merge together to ctt() object =============================================
ctt_obj<-list(segments=segs,nchanges=nchanges,edge.length=edge.length,tree=tree)
class(ctt_obj)<-"ctt"
ctt_obj

ctt_df <- data.frame(ctt_obj$segments,
      ctt_obj$nchanges,
      ctt_obj$edge.length)
names(ctt_df) <- c("seg_start_time", "seg_end_time", "n_changes", "edge.length")


## Adjusted Plot CTT function for rate

plot.ctt_optim<-function(x,...){
  h<-max(nodeHeights(x$tree))
  args<-list(...)
  type<-"rate"
  if(!is.null(args$show.tree)){
    show.tree<-args$show.tree
    args$show.tree<-NULL
  } else show.tree<-FALSE
  if(!is.null(args$add)){ 
    add<-args$add
    args$add<-NULL
  } else add<-FALSE
  if(is.null(args$ylim)) 
    args$ylim<-c(0,max(x$nchanges/x$edge.length))
  if(is.null(args$xlim))
    args$xlim<-c(max(x$segments),min(x$segments))
  if(is.null(args$lwd)) args$lwd<-2
  if(is.null(args$xlab)) args$xlab<-"time since the present"
  if(is.null(args$ylab)) 
    args$ylab<-"mean number of changes / total edge length"
  args$type<-"l"
  args$x<-h-as.vector(t(x$segments))
  args$y<-rbind(x$nchanges/x$edge.length,x$nchanges/x$edge.length)
  if(!add) do.call(plot,args)
  else do.call(lines,args)
  if(show.tree) plotTree(x$tree,add=TRUE,ftype="off",lwd=1,
                         color=make.transparent("blue",0.1),mar=par()$mar,
                         direction="leftwards",xlim=args$xlim)
}


ctt_plot <- plot.ctt_optim(ctt_obj)















ageL <- max(nodeHeights(tree))# root age: 41.14834

changes_table_Legumes <- changes_table_Legumes %>% 
  mutate(time_new = ageL - time) %>%  # reshape time column so the 0 is Present and 41 is the root
  mutate(time_new = round(time_new,5)) %>% 
  select(-time)

changes_table_Legumes %>% 
  group_by(N_sim, transition) %>% 
  summarise(N = n(), avg_time = mean(time_new))

changes_table_Legumes %>% 
  group_by(transition, trans_no) %>% 
  summarise(avg_time=mean(time_new)) %>% 
  View()

mean_plotdf <- changes_table_Legumes %>% 
  group_by(transition, trans_no) %>% 
  summarise(N=n(), avg_time = mean(time_new)) 

mean_cols <- setNames( c("darkgrey", "#FFBB00"), c(unique(mean_plotdf$transition)))


### Good Figure: 

plot_df <- changes_table_Legumes %>% 
  mutate(time_rounded = round(time_new, 0)) %>% 
  group_by(transition, time_rounded) %>% summarise(N = n()) %>% 
  distinct() 

### Export via Pane -> Export -> as PDF -> Cairo device -> width 6 x height 3.5: Name : Legumes_CumTrans_Bar.pdf
Fig1a_L <- ggplot()+
  geom_col(aes(x = time_rounded, y = N), fill = "lightgrey", col = "darkgrey",  
           data = plot_df %>% filter(transition == unique(plot_df$transition)[1]))+
  geom_col(aes(x = time_rounded, y = N), fill = "#ffd500", col = "#FFBB00", alpha = 0.2,
           data = plot_df %>% filter(transition == unique(plot_df$transition)[2]))+
 ylab("Cumulative Number of Transitions") +
  xlab("Time since Present (Ma) [Resolution 1 Ma]") +
  theme_classic()+
  scale_x_reverse()+
  scale_fill_manual(values = mean_cols)+
  labs(title ="Number of transitions per 1 million years (Legumes)")


plot_df <- changes_table_Legumes %>% 
  mutate(time_rounded = round(time_new, 0)) %>% 
  group_by(transition, N_sim, time_rounded) %>%
  summarise(mean_N = mean(length(trans_no))) %>%
  ungroup() %>%
  group_by(transition, time_rounded) %>%
  summarise(N = mean(mean_N))

Fig1b_L <- ggplot()+
  geom_col(aes(x = time_rounded, y = N), fill = "lightgrey", col = "darkgrey",  
           data = plot_df %>% filter(transition == unique(plot_df$transition)[1]))+
  geom_col(aes(x = time_rounded, y = N), fill = "#ffd500", col = "#FFBB00", alpha = 0.2,
           data = plot_df %>% filter(transition == unique(plot_df$transition)[2]))+
  ylab("Mean Cumulative Number of Transitions") +
  xlab("Time since Present (Ma) [Resolution 1 Ma]") +
  theme_classic()+
  scale_x_reverse()+
  scale_fill_manual(values = mean_cols)+
  labs(title ="Mean Number of transitions per 1 million years (Legumes)")


gridExtra::grid.arrange(Fig1a_M, Fig1a_L, nrow = 2)

#### Mammals ================================================================== 
mtrees_M <- readRDS("~/GitHub/Evo-Legumes-Herbivory/LegumesArmature/output/RDS/Mammals_mtrees_multiSimmap.rds")
changes_table_Mammals <- readRDS("~/GitHub/Evo-Legumes-Herbivory/LegumesArmature/output/RDS/changes_table_mammals.rds")




ageM <- max(nodeHeights(mtrees_M[[1]])) # root age: 217.84 Ma years old

changes_table_Mammals <- changes_table_Mammals %>% 
  mutate(time_new = ageM - time) %>% 
  mutate(time_new = round(time_new,5)) %>% 
  select(-time)

changes_table_Mammals %>% 
  group_by(N_sim, transition) %>% 
  summarise(N = n(), avg_time = mean(time_new))

changes_table_Mammals %>% 
  group_by(transition, trans_no) %>% 
  summarise(avg_time=mean(time_new)) %>% 
  View()

mean_plotdf <- changes_table_Mammals %>% 
  group_by(transition, trans_no) %>% 
  summarise(N=n(), avg_time = mean(time_new)) 

# Plots ====
#  transition times (cummulative)
mean_cols <- setNames( c("darkgrey", "#FFBB00"), c(unique(mean_plotdf$transition)))
# --------------------------------------- #

plot_df <- changes_table_Mammals %>% 
  mutate(time_rounded = round(time_new, 0)) %>% 
  group_by(transition, time_rounded) %>% summarise(N = n()) %>% 
  distinct() 

Fig1a_M <- ggplot()+
  geom_col(aes(x = time_rounded, y = N), fill = "lightgrey", col = "darkgrey",
           data = plot_df %>% filter(transition == unique(plot_df$transition)[1]))+
  geom_col(aes(x = time_rounded, y = N), fill = "#ffd500", col = "#FFBB00", alpha = 0.2,
           data = plot_df %>% filter(transition == unique(plot_df$transition)[2]))+
  ylab("Cumulative Number of Transitions") +
  xlab("Time since Present (Ma) [Resolution 1 Ma]") + xlim(45, 0) +
  theme_classic()+
  scale_fill_manual(values = mean_cols)+
  labs(title ="Number of transitions per 1 million years (Mammals)")

## Mean cumulative rates (averaged across simulations)
plot_df <- changes_table_Mammals %>% 
  mutate(time_rounded = round(time_new, 0)) %>% 
  group_by(transition, N_sim, time_rounded) %>%
  summarise(mean_N = mean(length(trans_no))) %>%
  ungroup() %>%
  group_by(transition, time_rounded) %>%
  summarise(N = mean(mean_N))

Fig1b_M <- ggplot()+
  geom_col(aes(x = time_rounded, y = N), fill = "lightgrey", col = "darkgrey",
           data = plot_df %>% filter(transition == unique(plot_df$transition)[1]))+
  geom_col(aes(x = time_rounded, y = N), fill = "#ffd500", col = "#FFBB00", alpha = 0.2,
           data = plot_df %>% filter(transition == unique(plot_df$transition)[2]))+
  ylab("Mean Cumulative Number of Transitions") +
  xlab("Time since Present (Ma) [Resolution 1 Ma]") + xlim(45, 0) +
  theme_classic()+
  scale_fill_manual(values = mean_cols)+
  labs(title ="Mean Number of transitions per 1 million years (Mammals)")


gridExtra::grid.arrange(Fig1b_L, Fig1b_M, nrow = 2)


## Together ====================== ###
### Export via Pane -> Export -> as PDF -> Cairo device -> width 6 x height 3.5: Name : Mammals_CumTrans_Bar.pdf
ggplot()+
  geom_col(aes(x = time_rounded, y = N), fill = "lightgrey", col = "darkgrey",  
           data = plot_df %>% filter(transition == unique(mean_plotdf$transition)[1]))+
  geom_col(aes(x = time_rounded, y = N), fill = "#ffd500", col = "#FFBB00", alpha = 0.2,
           data = plot_df %>% filter(transition == unique(mean_plotdf$transition)[2]))+
  ylim(0, 50)+ ylab("Cumulative Number of Transitions") +
  xlim(230, -2) + xlab("Time since Present (Ma) [Resolution 1 Ma]") +
  theme_classic()+
  scale_fill_manual(values = mean_cols)+
  labs(title ="Number of transitions per 1 million years (Megaherbivores)")




# average and mean transition times (cumulative)
p1 <- ggplot()+
  geom_point(aes(x = time_new, y = trans_no), shape = 21, col = "lightgrey", fill = "lightgrey", alpha = 0.3, 
             data = changes_table_Mammals %>% filter(transition == unique(mean_plotdf$transition)[1]))+
  geom_point(aes(x = time_new, y = trans_no), shape = 24, col = "#FFBB00", fill= "#FFBB00", alpha = 0.3, 
             data = changes_table_Mammals%>% filter(transition == unique(mean_plotdf$transition)[2]))+
  scale_x_reverse()+
  ylim(0, 50)+
  xlim(230, -5) #+theme_classic()


p1 + 
  geom_point(aes(x = avg_time, y = trans_no), shape = 21, col = "black",  fill = "darkgrey", 
             data = mean_plotdf %>% filter(transition == unique(mean_plotdf$transition)[1])) +
  geom_point(aes(x = avg_time, y = trans_no), shape = 24, col = "black", fill= "#ffd500",
             data = mean_plotdf %>% filter(transition == unique(mean_plotdf$transition)[2]))


#########
changes_table_Mammals_full <- changes_table_Mammals_full %>% 
  mutate(time_new = age - time) %>% 
  mutate(time_new = round(time_new,5)) %>% 
  select(-time)


changes_table_Mammals_full %>% 
  mutate(time_rounded = round(time_new, 0)) %>% 
  group_by(transition, time_rounded) %>% summarise(N = n()) %>% 
  distinct() %>%
  
  ggplot()+
  geom_col(aes(x = time_rounded, y = N, fill = transition, alpha = 0.5))+
  ylim(0, 50)+
  xlim(230, -5) #+theme_classic()


changes_table_Mammals_full
