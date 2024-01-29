library(dplyr)
library(rstatix)
#### Data =====================================================================
mtrees_M <- readRDS("~/GitHub/Evo-Legumes-Herbivory/LegumesArmature/output/RDS/Mammals_mtrees_multiSimmap.rds")
changes_table_Mammals <- changes_table_wo_edge
changes_table_Mammals_full <- changes_table_w_edge
saveRDS(changes_table_Mammals_full, "~/GitHub/Evo-Legumes-Herbivory/LegumesArmature/output/RDS/changes_table_mammals_with_number_of_species.rds")
saveRDS(changes_table_Mammals, "~/GitHub/Evo-Legumes-Herbivory/LegumesArmature/output/RDS/changes_table_mammals.rds")


mtrees_L <- readRDS("~/GitHub/Evo-Legumes-Herbivory/LegumesArmature/output/RDS/Legumes_mtrees_multiSimmap.rds")
#changes_table_Legumes <- changes_table_wo_edge
changes_table_Legumes <- readRDS("~/GitHub/Evo-Legumes-Herbivory/LegumesArmature/output/RDS/changes_table_legumes.rds")
#saveRDS(changes_table_Legumes, "~/GitHub/Evo-Legumes-Herbivory/LegumesArmature/output/RDS/changes_table_legumes.rds")


#### Legumes ==================================================================
ageL <- max(nodeHeights(mtrees_L[[1]]))# root age: 41.14834

changes_table_Legumes <- changes_table_Legumes %>% 
  mutate(time_new = ageL - time) %>% 
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




### Good Figure: 

plot_df <- changes_table_Legumes %>% 
  mutate(time_rounded = round(time_new, 0)) %>% 
#  select(transition, time_rounded, trans_no) %>%
  group_by(transition, time_rounded) %>% summarise(N = n()) %>% 
  distinct() 



### Export via Pane -> Export -> as PDF -> Cairo device -> width 6 x height 3.5: Name : Legumes_CumTrans_Bar.pdf
ggplot()+
  geom_col(aes(x = time_rounded, y = N), fill = "lightgrey", col = "darkgrey",  
           data = plot_df %>% filter(transition == unique(plot_df$transition)[1]))+
  geom_col(aes(x = time_rounded, y = N), fill = "#ffd500", col = "#FFBB00", alpha = 0.2,
           data = plot_df %>% filter(transition == unique(plot_df$transition)[2]))+
  ylim(0, 50)+ ylab("Cumulative Number of Transitions") +
  xlab("Time since Present (Ma) [Resolution 1 Ma]") +
  theme_classic()+
  scale_x_reverse()+
  scale_fill_manual(values = mean_cols)+
  labs(title ="Number of transitions per 1 million years (Legumes)")

#### Mammals ================================================================== 
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
# average transition times (cummulative)
mean_cols <- setNames( c("darkgrey", "#FFBB00"), c(unique(mean_plotdf$transition)))


plot_df <- changes_table_Mammals %>% 
  mutate(time_rounded = round(time_new, 0)) %>% 
  group_by(transition, time_rounded) %>% summarise(N = n()) %>% 
  distinct() 
  
ggplot()+
  # geom_col(aes(x = time_rounded, y = N), fill = "lightgrey", col = "darkgrey",  
  #          data = plot_df %>% filter(transition == unique(mean_plotdf$transition)[1]))+
  geom_col(aes(x = time_rounded, y = N), fill = "#ffd500", col = "#FFBB00", alpha = 0.2,
           data = plot_df %>% filter(transition == unique(mean_plotdf$transition)[2]))+
  ylim(0, 50)+ ylab("Cumulative Number of Transitions") +
  xlim(230, -2) + xlab("Time since Present (Ma) [Resolution 1 Ma]") +
  theme_classic()+
scale_fill_manual(values = mean_cols)+
  labs(title ="Number of transitions per 1 million years")


ggplot()+
  geom_col(aes(x = time_rounded, y = N), fill = "lightgrey", col = "darkgrey",  
           data = plot_df %>% filter(transition == unique(mean_plotdf$transition)[1]))+
  # geom_col(aes(x = time_rounded, y = N), fill = "#ffd500", col = "#FFBB00", alpha = 0.2,
  #          data = plot_df %>% filter(transition == unique(mean_plotdf$transition)[2]))+
  ylim(0, 50)+ ylab("Cumulative Number of Transitions") +
  xlim(230, -2) + xlab("Time since Present (Ma) [Resolution 1 Ma]") +
  theme_classic()+
  scale_fill_manual(values = mean_cols)+
  labs(title ="Number of transitions per 1 million years")


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
