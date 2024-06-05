# Final script: Evolution Armature Legumes-Megaherbivores #
# 1. Accumulation of transition rates through time
# 2. Relative transition rates through time (corrected for species numbers)

rm(list=ls())
source("Script/functions.R")


# Libraries:
my_pcks <- c("phytools", "dplyr", "geiger", "tictoc", "tidyr", "ggplot2", "prettyGraphs")

# Packages: install (if needed) and load:
install_and_load(my_pcks)


## 1. Legumes ========================================================================================
legumes_simmap_100 <- readRDS("output/RDS/Legumes_simmap_100.rds") 
legumes_simmap_100_v2 <- legumes_simmap_100  #%>% lapply(ladderize)
class(legumes_simmap_100_v2) <- class(legumes_simmap_100)
rm(legumes_simmap_100)

mtrees <- legumes_simmap_100_v2

## 1.1 Data processing ==== 
# apply function across list of multiple simulations:
armature_transitions <- data.frame()
for(i in 1:length(mtrees)){
  temp <- cbind(i, transition_times(mtrees[[i]]))
  armature_transitions <- rbind(armature_transitions, temp)
}
rm(temp, i)

# Number of transitions across all simulations
table(armature_transitions$transition)

# build new data frame with cumulative number of transitions
armature_trans_cumul <- data.frame()
for(n in 1:100){
  trans <- armature_transitions %>% 
    dplyr::filter(i == n) %>% 
    dplyr::group_by(transition) %>%
    arrange(desc(time), by.group = T) %>%
    dplyr::mutate(trans_no = row_number(i))
  
  armature_trans_cumul <- rbind(armature_trans_cumul, trans)
}
rm(n, trans)


# copied from somewhere else in the same script....:
# first need to know average number of transitions across simulations, rounded
trans_avg_length <- armature_trans_cumul %>%
  dplyr::group_by(transition, i) %>%
  dplyr::mutate(no_trans = max(trans_no)) %>%
  dplyr::ungroup() %>%
  dplyr::filter(trans_no == no_trans) %>%
  dplyr::group_by(transition) %>%
  dplyr::mutate(avg_length = round(mean(no_trans), digits = 0)) %>%
  dplyr::select(transition, avg_length) %>%
  dplyr::ungroup() %>%
  dplyr::distinct()

# rearrange data and average times across rows by average time across simulations:
avg_trans_times <- armature_trans_cumul %>%
  dplyr::group_by(transition, trans_no) %>%
  dplyr::mutate(avg_time = mean(time)) %>%
  dplyr::mutate(SE_time = sqrt(var(time) / length(time))) %>%
  dplyr::select(transition, trans_no, avg_time, SE_time) %>%
  dplyr::distinct()

# reduce avg_trans_times to trans_avg_length
avg_trans_times_a2na <- avg_trans_times %>%
  dplyr::filter(transition == "armature->no_armature") %>%
  dplyr::filter(trans_no <= trans_avg_length[1,2])

avg_trans_times_na2a <- avg_trans_times %>%
  dplyr::filter(transition == "no_armature->armature") %>%
  dplyr::filter(trans_no <= trans_avg_length[2,2])
avg_trans_times <- rbind(avg_trans_times_a2na, avg_trans_times_na2a)
rm(avg_trans_times_a2na, avg_trans_times_na2a)



## 1.2 Plotting: Cumulative rates (yellow ~ grey) =============================================

# Colors
myColours = c("lightgrey", "#FFBB00")
myColoursAlpha <- add.alpha(myColours, alpha=0.2) # transparency
my_cols <- setNames(myColoursAlpha, c("armature->no_armature", "no_armature->armature"))

# Prepare data for simulations
simulations <- armature_transitions %>%
  group_by(i, transition) %>%
  mutate(trans_no = row_number()) %>%
  ungroup()

# Prepare data for average lines
avg_data <- avg_trans_times %>%
  mutate(transition = factor(transition, levels = c("armature->no_armature", "no_armature->armature")))

# Define the plot
Fig1a_Legumes_CumRates <- ggplot() +
  # Plot points and lines for 1000 simulations
  geom_point(data = simulations %>% filter(transition == "armature->no_armature"),
             aes(x = time, y = trans_no), color = my_cols["armature->no_armature"], shape = 15, alpha = 0.2) +
  geom_line(data = simulations %>% filter(transition == "armature->no_armature"),
            aes(x = time, y = trans_no), color = my_cols["armature->no_armature"], alpha = 0.2) +
  geom_point(data = simulations %>% filter(transition == "no_armature->armature"),
             aes(x = time, y = trans_no), color = my_cols["no_armature->armature"], shape = 17, alpha = 0.2) +
  geom_line(data = simulations %>% filter(transition == "no_armature->armature"),
            aes(x = time, y = trans_no), color = my_cols["no_armature->armature"], alpha = 0.2) +
  
  # Plot average points and lines
  geom_point(data = avg_data %>% filter(transition == "armature->no_armature"),
             aes(x = avg_time, y = trans_no), color = "darkgrey", shape = 15) +
  geom_line(data = avg_data %>% filter(transition == "armature->no_armature"),
            aes(x = avg_time, y = trans_no), color = "darkgrey", size = 1) +
  geom_point(data = avg_data %>% filter(transition == "no_armature->armature"),
             aes(x = avg_time, y = trans_no), color = "yellow", shape = 17) +
  geom_line(data = avg_data %>% filter(transition == "no_armature->armature"),
            aes(x = avg_time, y = trans_no), color = "yellow", size = 1) +
  
  # Customize plot
  xlim(220,0) +
  labs(x = "Time of transitions (mya)", 
       y = "Cumulative number of transitions",
       title = "Cumulative Transition Rates (Legumes)") +
  theme_classic(base_size = 18) +
  theme(legend.position = "top") +
  
  # Add legend
  scale_shape_manual(values = c(15, 17), name = "Transition",
                     labels = c("armature->no_armature", "no_armature->armature")) +
  scale_color_manual(values = myColours, name = "Transition",
                     labels = c("armature->no_armature", "no_armature->armature"))

# Save to powerpoint vector file =====
esquisse:::ggplot_to_ppt("Fig1a_Legumes_CumRates")







rm(list = ls())








## 2. Mammals =============================================
source("Script/functions.R")
mtrees <- readRDS("output/RDS/Mammals_simmap_100.rds")
# mtrees <- readRDS("output/RDS/Mammals_4traits_simmap_100.rds")

## 2.1 Data processing =================================

# Ref[0] -> Ref[3]


# apply function across list of multiple simulations:

megaherbivore_transitions <- data.frame()
for(i in 1:length(mtrees)){
  temp <- cbind(i, transition_times(mtrees[[i]]))
  megaherbivore_transitions <- rbind(megaherbivore_transitions, temp)
}
rm(temp, i)

# Number of transitions across all simulations
table(megaherbivore_transitions$transition)

# build new data frame with cumulative number of transitions
megaherbivore_trans_cumul <- data.frame()
for(n in 1:100){
  trans <- megaherbivore_transitions %>% 
    dplyr::filter(i == n) %>% 
    dplyr::group_by(transition) %>%arrange(desc(time), by.group = T) %>%
    dplyr::mutate(trans_no = row_number(i))
  
  megaherbivore_trans_cumul <- rbind(megaherbivore_trans_cumul, trans)
}
rm(n, trans)


##### copied from somewhere else in the same script....:

# first need to know average number of transitions across simulations, rounded
trans_avg_length <- megaherbivore_trans_cumul %>%
  dplyr::group_by(transition, i) %>%
  dplyr::mutate(no_trans = max(trans_no)) %>%
  dplyr::ungroup() %>%
  dplyr::filter(trans_no == no_trans) %>%
  dplyr::group_by(transition) %>%
  dplyr::mutate(avg_length = round(mean(no_trans), digits = 0)) %>%
  dplyr::select(transition, avg_length) %>%
  dplyr::ungroup() %>%
  dplyr::distinct()

avg_trans_times <- megaherbivore_trans_cumul %>%
  dplyr::group_by(transition, trans_no) %>%
  dplyr::mutate(avg_time = mean(time)) %>%
  dplyr::mutate(SE_time = sqrt(var(time) / length(time))) %>%
  dplyr::select(transition, trans_no, avg_time, SE_time) %>%
  dplyr::distinct()

avg_trans_list <- list()
for (i in seq_along(unique(trans_avg_length$transition))){ # i in 10
  type <- unique(trans_avg_length$transition)[i]
  avg_trans_trait <- avg_trans_times %>%
    dplyr::filter(transition == type) %>%
    dplyr::filter(trans_no <= trans_avg_length[i,2]) 
  
  avg_trans_list[[i]] <- avg_trans_trait 
}

avg_trans_df <- do.call(rbind, avg_trans_list)


## 2.2 Plotting: Cumulative Rates Plot =========================

# Colors
myColours = c("lightgrey", "#FFBB00")
myColoursAlpha <- add.alpha(myColours, alpha=0.2) # transparency
my_cols <- setNames(myColoursAlpha, c("megaherbivore->other", "other->megaherbivore"))

# Prepare data for simulations
simulations <- megaherbivore_transitions %>%
  group_by(i, transition) %>%
  mutate(trans_no = row_number()) %>%
  ungroup()

# Prepare data for average lines
avg_data <- avg_trans_times %>%
  filter(transition %in% c("megaherbivore->other", "other->megaherbivore"))

# Define the plot
Fig1b_Mammals_CumRates <- ggplot() +
  # Plot points and lines for 100 simulations
  geom_point(data = simulations %>% filter(transition == "megaherbivore->other"),
             aes(x = time, y = trans_no), color = my_cols["megaherbivore->other"], shape = 15, alpha = 0.2) +
  geom_line(data = simulations %>% filter(transition == "megaherbivore->other"),
            aes(x = time, y = trans_no), color = my_cols["megaherbivore->other"], alpha = 0.2) +
  geom_point(data = simulations %>% filter(transition == "other->megaherbivore"),
             aes(x = time, y = trans_no), color = my_cols["other->megaherbivore"], shape = 17, alpha = 0.2) +
  geom_line(data = simulations %>% filter(transition == "other->megaherbivore"),
            aes(x = time, y = trans_no), color = my_cols["other->megaherbivore"], alpha = 0.2) +
  
  # Plot average points and lines
  geom_point(data = avg_data %>% filter(transition == "megaherbivore->other"),
             aes(x = avg_time, y = trans_no), color = "darkgrey", shape = 15) +
  geom_line(data = avg_data %>% filter(transition == "megaherbivore->other"),
            aes(x = avg_time, y = trans_no), color = "darkgrey", size = 1) +
  geom_point(data = avg_data %>% filter(transition == "other->megaherbivore"),
             aes(x = avg_time, y = trans_no), color = "yellow", shape = 17) +
  geom_line(data = avg_data %>% filter(transition == "other->megaherbivore"),
            aes(x = avg_time, y = trans_no), color = "yellow", size = 1) +
  
  # Customize plot
  scale_x_reverse(limits = c(220, 0)) +
  ylim(0,110) +
  labs(x = "Time of transitions (mya)", y = "Cumulative number of transitions",
       title = "Cumulative Transition Rates (Megaherbivores)") +
  theme_classic(base_size = 18) +
  theme(legend.position = "right") +
  
  # Add legend
  scale_shape_manual(values = c(15, 17), name = "Transition",
                     labels = c("megaherbivore->other", "other->megaherbivore")) +
  scale_color_manual(values = myColours, name = "Transition",
                     labels = c("megaherbivore->other", "other->megaherbivore"))


Fig1b_Mammals_CumRates
# Save to powerpoint ===============
esquisse:::ggplot_to_ppt("Fig1b_Mammals_CumRates")
