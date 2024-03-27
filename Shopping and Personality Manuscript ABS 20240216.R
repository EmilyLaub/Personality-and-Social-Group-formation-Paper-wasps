# Analysis for manuscript "Personality and body mass impact social group formation and function in Paper wasps"
# Code written by Emily C. Laub (eclaub@ucla.edu; eclaub@umich.edu)

############# load libraries for analysis 
library(lme4)
library(car)
library(ggplot2)
library(dplyr)
library(igraph)
library(MASS)
library(geepack)
library(performance)

############# ggplot theme
theme_mine2 <- function(base_size = 20, base_family = "Helvetica") {
  # Starts with theme_grey and then modify some parts
  theme_bw(base_size = base_size, base_family = base_family) %+replace%
    theme(
      strip.background = element_blank(),
      strip.text.x = element_text(size = 20),
      strip.text.y = element_text(size = 20),
      axis.text.x = element_text(size=22),
      axis.text.y = element_text(size=22,hjust=1),
      axis.ticks =  element_line(colour = "black"), 
      axis.title.x= element_text(size=28),
      axis.title.y= element_text(size=28,angle=90),
      panel.background = element_blank(), 
      panel.border =element_blank(), 
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank(), 
      panel.margin = unit(1.0, "lines"), 
      plot.background = element_blank(), 
      plot.margin = unit(c(1,  1, 1, 1), "lines"),
      axis.line.x = element_line(color="black", size = 1),
      axis.line.y = element_line(color="black", size = 1)
    )
}




############ load data 
## Wasp individual attributes (body weight, personality traits), partner assessment, nest size, nest site assessment 
shop <- read.csv("Wasp_Attributes_PersonalityShopping_20240217.csv")

head(shop)
hist(shop$Days_observed)

mean(shop$Days_observed)
sd(shop$Days_observed)

nrow(shop)
hist(shop$max_stable_group)
##################################################
# Personality measures:
# Explore = exploratory personality
# body_contact = affiliation personality
# Log_agg = aggression personality
# neutral = investigation personality, neutral investigation of dummy (antennation)

##############################################################
####### Shopping and individual attributes

hist(shop$Nest_sites_sampled)
### Use histograms to examine data distribution, choose poisson for nest sites

head(shop)

# Nest sites sampled
boxes1 <- glm( Nest_sites_sampled ~ body_contact + Explore + Log_agg +neutral + Weight, family = "poisson", data = shop)

summary(boxes1)
Anova(boxes1)

model_performance(boxes1)
check_model(boxes1)

#### create graphs

ggplot(shop, aes(body_contact, Nest_sites_sampled)) +
  geom_point() +
  xlab("Affiliation") + ylab("Roosting sites sampled") +
  theme_mine2() +
  theme(axis.title.x = element_text(margin=margin(t= 6)), #add margin to x-axis title
        axis.title.y = element_text(margin=margin(r= 6))) +       #
  scale_y_continuous(name="Roosting sites sampled", limits=c(0, 25)) +
  stat_smooth(method="lm")


ggplot(shop, aes(Weight, Nest_sites_sampled)) +
  geom_point() +
  xlab("Body mass") + ylab("Roosting sites sampled") +
  theme_mine2() +
  theme(axis.title.x = element_text(margin=margin(t= 6)), #add margin to x-axis title
        axis.title.y = element_text(margin=margin(r= 6))) +
  scale_y_continuous(name="Roosting sites sampled", limits=c(0, 25)) +
  stat_smooth(method="lm", se=TRUE)


#################### filter only wasps that chose nests
att_nests = 
  shop %>%
  filter(nest_choice == "y") %>%
  as.data.frame()

####### Days to choose a nest
hist(att_nests$Individual.choice.day..partners...3.time.seeing.together.consequetive..members.of.same.ultimate.group..)
##### negative binomial distribution


sitechoiceday1 <- glm.nb(Individual.site.choice.day..3.day.of.same.consequetive.location. ~ 
                           body_contact + Explore + Log_agg +neutral + Weight, link = "log", data = att_nests)

summary(sitechoiceday1)
Anova(sitechoiceday1)


ggplot(att_nests, aes(Explore , Individual.site.choice.day..3.day.of.same.consequetive.location.)) +
  geom_point() +
  xlab("Exploration") + ylab("Days to choose nest site") +
  theme_mine2() +
  scale_y_continuous(name="Days to choose nest site", limits=c(0, 50)) +
  stat_smooth(method="lm", se=TRUE)

##########Partner Choice Day
hist(att_nests$Individual.choice.day..partners...3.time.seeing.together.consequetive..members.of.same.ultimate.group..)
# negative binomial distribution

partners4 <- glm.nb(Individual.choice.day..partners...3.time.seeing.together.consequetive..members.of.same.ultimate.group.. ~ 
                      body_contact + Explore + Log_agg +neutral + Weight, link = "log", data = att_nests)
summary(partners4)
Anova(partners4)

check_model(partners4)


ggplot(att_nests, aes(Explore , Individual.choice.day..partners...3.time.seeing.together.consequetive..members.of.same.ultimate.group..)) +
  geom_point() +
  xlab("Exploration") + ylab("Days to choose partners") +
  theme_mine2() +
  scale_y_continuous(name="Days to choose partners", limits=c(0, 50)) +
  stat_smooth(method="lm", se=TRUE)



############################### Social network analysis 
# Load interaction matrix 
int2021 = read.csv("interactionmatrixall2021_cleaned_20230209.csv")

### Get data in symmetrical matrix for analysis
row.names(int2021) <- int2021[,1]
head(int2021)
int2021<- int2021[, -1]
shopping_mat2021<- data.matrix(int2021, rownames.force = NA)
isSymmetric(shopping_mat2021)
# Check that matrix is symmetric before making graph

wasps2021 = graph_from_adjacency_matrix(adjmatrix = shopping_mat2021, mode = "undirected",
                                        weighted = T, diag = F) 

#### Get attributes in same order as node names
wasp_attributes2021 = data.frame(matrix(nrow = 52))
wasp_attributes2021$names <- data.frame(names = V(wasps2021)$name)
wasp_attributes2021 <- wasp_attributes2021[,-1]



str(wasp_attributes2021)

wasp_attributes2021$names = as.character(wasp_attributes2021$names)


wasp_attributes2021$Weight = shop$Weight[match(wasp_attributes2021$names, 
                                               shop$Color)]

wasp_attributes2021$Explore = shop$Explore[match(wasp_attributes2021$names, 
                                                 shop$Color)]
wasp_attributes2021$Log_agg = shop$Log_agg[match(wasp_attributes2021$names, 
                                                 shop$Color)]
wasp_attributes2021$body_contact = shop$body_contact[match(wasp_attributes2021$names, 
                                                           shop$Color)]
wasp_attributes2021$neutral = shop$neutral[match(wasp_attributes2021$names, 
                                                 shop$Color)]

####### Assign attributes to nodes

V(wasps2021)$Weight <- wasp_attributes2021$Weight
V(wasps2021)$Explore <- wasp_attributes2021$Explore
V(wasps2021)$Log_agg <- wasp_attributes2021$Log_agg
V(wasps2021)$neutral <- wasp_attributes2021$neutral
V(wasps2021)$body_contact <- wasp_attributes2021$body_contact

#####################################################################
# Calculate Assortativity - likelihood of wasps forming ties with wasps that share attribute
#####################################################################

#######################################
# Explore Assortativity
#######################################
V(wasps2021)$Explore <- wasp_attributes2021$Explore
obs_ass_Explore <- assortativity(wasps2021, types = as.numeric((V(wasps2021)$Explore)), types2 = NULL, directed = FALSE) ### set to variable so it won't be recalculated

itr = 1000 ### set iteration count. This will set how many null models are created.

all_refass_Explore <- rep(NA, itr) #### data frame that will get filled.

##### now run for loop that will shuffle node attributes and calculate assortativity therefrom
for (i in 1:itr) {
  shuffnet <- wasps2021
  V_Explore <- (V(shuffnet)$Explore) ### get list of vertex values
  shuff <- sample(V_Explore) #### get random sample of vertex values
  
  V(shuffnet)$Explore = shuff #### assign  random sample to vertecies
  ref_ass <- assortativity(shuffnet, types = as.numeric((V(shuffnet)$Explore)), types2 = NULL, directed = FALSE) #### calculate assortativity
  
  all_refass_Explore[i] <- ref_ass ### add to all_refass list
}

### again, VERY IMPORTANT to make sure that the observed assortativity didn't get recalculated


rank_obs_ass_Explore <- sum(all_refass_Explore <= obs_ass_Explore)

# calculate the quantile based on the rank of the observed
quantile_value_ass_Explore <- rank_obs_ass_Explore / length(all_refass_Explore)

# calculate 2-tailed p-value
p_ass_Explore <- min(c(quantile_value_ass_Explore,(1-quantile_value_ass_Explore)))*2

p_ass_Explore

#######################################
# Log_agg Assortativity
#######################################
V(wasps2021)$Log_agg <- wasp_attributes2021$Log_agg
obs_ass_Log_agg <- assortativity(wasps2021, types = as.numeric((V(wasps2021)$Log_agg)), types2 = NULL, directed = FALSE) ### set to variable so it won't be recalculated

itr = 1000 ### set iteration count. This will set how many null models are created.

all_refass_Log_agg <- rep(NA, itr) #### data frame that will get filled.

##### now run for loop that will shuffle node attributes and calculate assortativity therefrom
for (i in 1:itr) {
  shuffnet <- wasps2021
  V_Log_agg <- (V(shuffnet)$Log_agg) ### get list of vertex values
  shuff <- sample(V_Log_agg) #### get random sample of vertex values
  
  V(shuffnet)$Log_agg = shuff #### assign  random sample to vertecies
  ref_ass <- assortativity(shuffnet, types = as.numeric((V(shuffnet)$Log_agg)), types2 = NULL, directed = FALSE) #### calculate assortativity
  
  all_refass_Log_agg[i] <- ref_ass ### add to all_refass list
}

### again, VERY IMPORTANT to make sure that the observed assortativity didn't get recalculated




rank_obs_ass_Log_agg <- sum(all_refass_Log_agg <= obs_ass_Log_agg)

# calculate the quantile based on the rank of the observed
quantile_value_ass_Log_agg <- rank_obs_ass_Log_agg / length(all_refass_Log_agg)

# calculate 2-tailed p-value
p_ass_Log_agg <- min(c(quantile_value_ass_Log_agg,(1-quantile_value_ass_Log_agg)))*2


#############################################################################
# body contact assortatvity
#############################################################################

V(wasps2021)$body_contact <- wasp_attributes2021$body_contact
obs_ass_body_contact <- assortativity(wasps2021, types = as.numeric((V(wasps2021)$body_contact)), types2 = NULL, directed = FALSE) ### set to variable so it won't be recalculated

itr = 1000 ### set iteration count. This will set how many null models are created.

all_refass_body_contact <- rep(NA, itr) #### data frame that will get filled.

##### now run for loop that will shuffle node attributes and calculate assortativity therefrom
for (i in 1:itr) {
  shuffnet <- wasps2021
  V_body_contact <- (V(shuffnet)$body_contact) ### get list of vertex values
  shuff <- sample(V_body_contact) #### get random sample of vertex values
  
  V(shuffnet)$body_contact = shuff #### assign  random sample to vertecies
  ref_ass <- assortativity(shuffnet, types = as.numeric((V(shuffnet)$body_contact)), types2 = NULL, directed = FALSE) #### calculate assortativity
  
  all_refass_body_contact[i] <- ref_ass ### add to all_refass list
}

### again, VERY IMPORTANT to make sure that the observed assortativity didn't get recalculated

rank_obs_ass_body_contact <- sum(all_refass_body_contact <= obs_ass_body_contact)

# calculate the quantile based on the rank of the observed
quantile_value_ass_body_contact <- rank_obs_ass_body_contact / length(all_refass_body_contact)

# calculate 2-tailed p-value
p_ass_body_contact <- min(c(quantile_value_ass_body_contact,(1-quantile_value_ass_body_contact)))*2

#########################################################################
# Assortativity investigation
########################################################################

V(wasps2021)$neutral <- wasp_attributes2021$neutral
obs_ass_neutral <- assortativity(wasps2021, types = as.numeric((V(wasps2021)$neutral)), types2 = NULL, directed = FALSE) ### set to variable so it won't be recalculated

itr = 1000 ### set iteration count. This will set how many null models are created.

all_refass_neutral <- rep(NA, itr) #### data frame that will get filled.

##### now run for loop that will shuffle node attributes and calculate assortativity therefrom
for (i in 1:itr) {
  shuffnet <- wasps2021
  V_neutral <- (V(shuffnet)$neutral) ### get list of vertex values
  shuff <- sample(V_neutral) #### get random sample of vertex values
  
  V(shuffnet)$neutral = shuff #### assign  random sample to vertecies
  ref_ass <- assortativity(shuffnet, types = as.numeric((V(shuffnet)$neutral)), types2 = NULL, directed = FALSE) #### calculate assortativity
  
  all_refass_neutral[i] <- ref_ass ### add to all_refass list
}

### again, VERY IMPORTANT to make sure that the observed assortativity didn't get recalculated

### Now graph reference and then add in observed w/ abline
hist(all_refass_neutral,
     main = "Assortativity of body contact, 2021",
     xlab = "Reference Assortativity")
abline(v=obs_ass_neutral,lwd=2, col='purple')


rank_obs_ass_neutral <- sum(all_refass_neutral <= obs_ass_neutral)

# calculate the quantile based on the rank of the observed
quantile_value_ass_neutral <- rank_obs_ass_neutral / length(all_refass_neutral)

# calculate 2-tailed p-value
p_ass_neutral <- min(c(quantile_value_ass_neutral,(1-quantile_value_ass_neutral)))*2


#########################################################################
# Assortativity Weight
########################################################################

V(wasps2021)$Weight <- wasp_attributes2021$Weight
obs_ass_Weight <- assortativity(wasps2021, types = as.numeric((V(wasps2021)$Weight)), types2 = NULL, directed = FALSE) ### set to variable so it won't be recalculated

itr = 1000 ### set iteration count. This will set how many null models are created.

all_refass_Weight <- rep(NA, itr) #### data frame that will get filled.

##### now run for loop that will shuffle node attributes and calculate assortativity therefrom
for (i in 1:itr) {
  shuffnet <- wasps2021
  V_Weight <- (V(shuffnet)$Weight) ### get list of vertex values
  shuff <- sample(V_Weight) #### get random sample of vertex values
  
  V(shuffnet)$Weight = shuff #### assign  random sample to vertecies
  ref_ass <- assortativity(shuffnet, types = as.numeric((V(shuffnet)$Weight)), types2 = NULL, directed = FALSE) #### calculate assortativity
  
  all_refass_Weight[i] <- ref_ass ### add to all_refass list
}

### again, VERY IMPORTANT to make sure that the observed assortativity didn't get recalculated

### Now graph reference and then add in observed w/ abline
hist(all_refass_Weight,
     main = "Assortativity of body contact, 2021",
     xlab = "Reference Assortativity")
abline(v=obs_ass_Weight,lwd=2, col='purple')


rank_obs_ass_Weight <- sum(all_refass_Weight <= obs_ass_Weight)

# calculate the quantile based on the rank of the observed
quantile_value_ass_Weight <- rank_obs_ass_Weight / length(all_refass_Weight)

# calculate 2-tailed p-value
p_ass_Weight <- min(c(quantile_value_ass_Weight,(1-quantile_value_ass_Weight)))*2






###############################################################################################
# Analysis of how different attributes influence Degree: number of unique individuals sampled
###############################################################################################

V(wasps2021)$degree <- degree(wasps2021, v = V(wasps2021), mode = "all",  loops = FALSE, normalized = FALSE)

wasp_attributes2021$degree <- V(wasps2021)$degree

degree1 <- glm( degree ~ Explore + body_contact + Log_agg + neutral + Weight, 
                family = "poisson", data = wasp_attributes2021)

degsum <- summary(degree1)

obs_degree_p_Explore = degsum$coefficients[2]
obs_degree_p_bodycontact = degsum$coefficients[3]
obs_degree_p_logagg = degsum$coefficients[4]
obs_degree_p_neutral = degsum$coefficients[5]
obs_degree_p_weight = degsum$coefficients[6]

itr <- 1000

all_cor_p_degree_Explore <- rep(NA, itr)
all_cor_p_degree_Weight <- rep(NA, itr)
all_cor_p_degree_bodycontact <- rep(NA, itr)
all_cor_p_degree_logagg <- rep(NA, itr)
all_cor_p_degree_neutral <- rep(NA, itr)

for (i in 1:itr) {
  wasps_shuff <- graph_from_adjacency_matrix(adjmatrix = shopping_mat2021, mode = "undirected",
                                             weighted = T, diag = F) 
  names1 <- V(wasps_shuff)$name ### get list of vertex values
  shuff <- sample(names1) #### get random sample of vertex values
  V(wasps_shuff)$name <- shuff #### assign  random sample to vertecies
  V(wasps_shuff)$degree <- degree(wasps_shuff, v = V(wasps_shuff), mode = "all", loops = FALSE, normalized = FALSE )
  wasp_attributes_shuff = data.frame(matrix(nrow = 52))
  wasp_attributes_shuff$names <- as.factor(V(wasps_shuff)$name)
  wasp_attributes_shuff$degree <- V(wasps_shuff)$degree
  wasp_attributes_shuff$Weight = shop$Weight[match(wasp_attributes_shuff$names, 
                                                   shop$Color)]
  wasp_attributes_shuff$Explore = shop$Explore[match(wasp_attributes_shuff$names, 
                                                     shop$Color)]
  wasp_attributes_shuff$Log_agg = shop$Log_agg[match(wasp_attributes_shuff$names, 
                                                     shop$Color)]
  wasp_attributes_shuff$body_contact = shop$body_contact[match(wasp_attributes_shuff$names, 
                                                               shop$Color)]
  wasp_attributes_shuff$neutral = shop$neutral[match(wasp_attributes_shuff$names, 
                                                     shop$Color)]
  degree_shuff <- glm( degree ~  Explore + body_contact + Log_agg + neutral + Weight, 
                       family = "poisson", data = wasp_attributes_shuff)
  degsum_shuff <- summary(degree_shuff)
  all_cor_p_degree_Explore[i] <- degsum_shuff$coefficients[2]
  all_cor_p_degree_bodycontact[i] <- degsum_shuff$coefficients[3]
  all_cor_p_degree_logagg[i] <- degsum_shuff$coefficients[4]
  all_cor_p_degree_neutral[i] <- degsum_shuff$coefficients[5]
  all_cor_p_degree_Weight[i] <- degsum_shuff$coefficients[6]
  
}



####### Explore

rank_obs_degsum_degree_Explore <- sum(all_cor_p_degree_Explore <= obs_degree_p_Explore)

# calculate the quantile based on the rank of the observed
quantile_value_degree_Explore <- rank_obs_degsum_degree_Explore / length(all_cor_p_degree_Explore)

# calculate 2-tailed p-value
p_Explore <- min(c(quantile_value_degree_Explore,(1-quantile_value_degree_Explore)))*2


###### Weight
rank_obs_degsum_degree_Weight <- sum(all_cor_p_degree_Weight <= obs_degree_p_Weight)

# calculate the quantile based on the rank of the observed
quantile_value_degree_Weight <- rank_obs_degsum_degree_Weight / length(all_cor_p_degree_Weight)

# calculate 2-tailed p-value
p_Weight <- min(c(quantile_value_degree_Weight,(1-quantile_value_degree_Weight)))*2


###### neutral
rank_obs_degsum_degree_neutral <- sum(all_cor_p_degree_neutral <= obs_degree_p_neutral)

# calculate the quantile based on the rank of the observed
quantile_value_degree_neutral <- rank_obs_degsum_degree_neutral / length(all_cor_p_degree_neutral)

# calculate 2-tailed p-value
p_neutral <- min(c(quantile_value_degree_neutral,(1-quantile_value_degree_neutral)))*2

###### bodycontact
rank_obs_degsum_degree_bodycontact <- sum(all_cor_p_degree_bodycontact <= obs_degree_p_bodycontact)

# calculate the quantile based on the rank of the observed
quantile_value_degree_bodycontact <- rank_obs_degsum_degree_bodycontact / length(all_cor_p_degree_bodycontact)

# calculate 2-tailed p-value
p_bodycontact <- min(c(quantile_value_degree_bodycontact,(1-quantile_value_degree_bodycontact)))*2

###### aggression
rank_obs_degsum_degree_logagg <- sum(all_cor_p_degree_logagg <= obs_degree_p_logagg)

# calculate the quantile based on the rank of the observed
quantile_value_degree_logagg <- rank_obs_degsum_degree_logagg / length(all_cor_p_degree_logagg)

# calculate 2-tailed p-value
p_logagg <- min(c(quantile_value_degree_logagg,(1-quantile_value_degree_logagg)))*2




library(ggplot2)
ggplot(wasp_attributes2021, aes(Explore , degree )) +
  geom_point() +
  xlab("Exploration") + ylab("Degree") +
  theme_mine2() +
  scale_y_continuous(name="Degree", limits=c(0, 30)) +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18))+
  stat_smooth(method="lm", se=FALSE)

########################################################################################################
###### Betweenness - global centrality of wasps
########################################################################################################

V(wasps2021)$betweenness <- betweenness(wasps2021, v = V(wasps2021), directed = FALSE, weights = NULL,
                                        nobigint = TRUE, normalized = FALSE)
wasp_attributes2021$betweenness <- V(wasps2021)$betweenness

between1 <- glm( betweenness ~  Explore + body_contact + Log_agg + neutral + Weight, 
                 family = "gaussian", data = wasp_attributes2021)

betsum <- summary(between1)
check_model(between1)
obs_between_p_Explore = betsum$coefficients[2]
obs_between_p_bodycontact = betsum$coefficients[3]
obs_between_p_logagg = betsum$coefficients[4]
obs_between_p_neutral = betsum$coefficients[5]
obs_between_p_weight = betsum$coefficients[6]

itr <- 1000

all_cor_p_between_Explore <- rep(NA, itr)
all_cor_p_between_Weight <- rep(NA, itr)
all_cor_p_between_bodycontact <- rep(NA, itr)
all_cor_p_between_logagg <- rep(NA, itr)
all_cor_p_between_neutral <- rep(NA, itr)

for (i in 1:itr) {
  wasps_shuff <- graph_from_adjacency_matrix(adjmatrix = shopping_mat2021, mode = "undirected",
                                             weighted = T, diag = F) 
  names1 <- V(wasps_shuff)$name ### get list of vertex values
  shuff <- sample(names1) #### get random sample of vertex values
  V(wasps_shuff)$name <- shuff #### assign  random sample to vertecies
  V(wasps_shuff)$between <- betweenness(wasps_shuff, v = V(wasps_shuff), directed = FALSE, weights = NULL,
                                        nobigint = TRUE, normalized = FALSE)
  wasp_attributes_shuff = data.frame(matrix(nrow = 52))
  wasp_attributes_shuff$names <- as.factor(V(wasps_shuff)$name)
  wasp_attributes_shuff$between <- V(wasps_shuff)$between
  wasp_attributes_shuff$Weight = shop$Weight[match(wasp_attributes_shuff$names, 
                                                   shop$Color)]
  wasp_attributes_shuff$Explore = shop$Explore[match(wasp_attributes_shuff$names, 
                                                     shop$Color)]
  wasp_attributes_shuff$Log_agg = shop$Log_agg[match(wasp_attributes_shuff$names, 
                                                     shop$Color)]
  wasp_attributes_shuff$body_contact = shop$body_contact[match(wasp_attributes_shuff$names, 
                                                               shop$Color)]
  wasp_attributes_shuff$neutral = shop$neutral[match(wasp_attributes_shuff$names, 
                                                     shop$Color)]
  between_shuff <- glm( between ~ Explore + body_contact + Log_agg + neutral + Weight, 
                        family = "gaussian", data = wasp_attributes_shuff)
  betsum_shuff <- summary(between_shuff)
  all_cor_p_between_Explore[i] <- betsum_shuff$coefficients[2]
  all_cor_p_between_bodycontact[i] <- betsum_shuff$coefficients[3]
  all_cor_p_between_logagg[i] <- betsum_shuff$coefficients[4]
  all_cor_p_between_neutral[i] <- betsum_shuff$coefficients[5]
  all_cor_p_between_Weight[i] <- betsum_shuff$coefficients[6]
  
}



####### Explore

rank_obs_betsum_between_Explore <- sum(all_cor_p_between_Explore <= obs_between_p_Explore)

# calculate the quantile based on the rank of the observed
quantile_value_between_Explore <- rank_obs_betsum_between_Explore / length(all_cor_p_between_Explore)

# calculate 2-tailed p-value
p_Explore <- min(c(quantile_value_between_Explore,(1-quantile_value_between_Explore)))*2


###### Weight
rank_obs_betsum_between_Weight <- sum(all_cor_p_between_Weight <= obs_between_p_weight)

# calculate the quantile based on the rank of the observed
quantile_value_between_Weight <- rank_obs_betsum_between_Weight / length(all_cor_p_between_Weight)

# calculate 2-tailed p-value
p_Weight <- min(c(quantile_value_between_Weight,(1-quantile_value_between_Weight)))*2


###### neutral
rank_obs_betsum_between_neutral <- sum(all_cor_p_between_neutral <= obs_between_p_neutral)

# calculate the quantile based on the rank of the observed
quantile_value_between_neutral <- rank_obs_betsum_between_neutral / length(all_cor_p_between_neutral)

# calculate 2-tailed p-value
p_neutral <- min(c(quantile_value_between_neutral,(1-quantile_value_between_neutral)))*2

###### bodycontact
rank_obs_betsum_between_bodycontact <- sum(all_cor_p_between_bodycontact <= obs_between_p_bodycontact)

# calculate the quantile based on the rank of the observed
quantile_value_between_bodycontact <- rank_obs_betsum_between_bodycontact / length(all_cor_p_between_bodycontact)

# calculate 2-tailed p-value
p_bodycontact <- min(c(quantile_value_between_bodycontact,(1-quantile_value_between_bodycontact)))*2

###### aggression
rank_obs_betsum_between_logagg <- sum(all_cor_p_between_logagg <= obs_between_p_logagg)

# calculate the quantile based on the rank of the observed
quantile_value_between_logagg <- rank_obs_betsum_between_logagg / length(all_cor_p_between_logagg)

# calculate 2-tailed p-value
p_logagg <- min(c(quantile_value_between_logagg,(1-quantile_value_between_logagg)))*2


###########################################################################################################
########## Number of stable days on a nest
###########################################################################################################
hist(att_nests$Stable_days_individual)
# best matches poisson distribution

stable_days2 <- glm(Stable_days_individual ~ 
                      body_contact + Log_agg + neutral + Explore + Weight, family = "poisson", data = att_nests)

summary(stable_days2)
Anova(stable_days2)

check_model(stable_days2)

ggplot(att_nests, aes(Explore , Stable_days_individual)) +
  geom_point() +
  xlab("Exploration") + ylab("Stable days on nest") +
  theme_mine2() +
  theme(axis.title.x = element_text(margin=margin(t= 6)), #add margin to x-axis title
        axis.title.y = element_text(margin=margin(r= 6))) +
  scale_y_continuous(name="Stable days on nest", limits=c(0, 50)) +
  stat_smooth(method="lm", se=TRUE)



################################################################################
###### On-nest behavior and personality
################################################################################

ons <- read.csv("On nest 2021 Personality 20240217.csv")
head(ons)

## Personality measures:

# Explore = exploratory personality
# body_contact = affiliation personality
# Log_agg = aggression personality
# neutral = investigation personality, neutral investigation of dummy (antennation)

## Nest ID column = X.2

########## On-Nest Aggression
str(ons)

### add personality to columns
aggression_on_nest <- glmer.nb(sum_agg_initiate_nodart_1 ~ wasp.weight +  Explore + neutral + body_contact + Log_agg + temp+ + number.wasps.on.nest + Vid_Day + log_mins_on_nest+ Big_larvae+ (1|X.2/wasp_ID), 
                               data = ons, verbose = FALSE)
summary(aggression_on_nest)
Anova(aggression_on_nest)
check_model(aggression_on_nest)


ggplot(ons, aes(Explore, sum_agg_initiate_nodart_1)) +
  geom_point() +
  xlab("Exploration") + ylab("Aggression") +
  theme_mine2() +
  theme(axis.title.x = element_text(margin=margin(t= 6)), #add margin to x-axis title
        axis.title.y = element_text(margin=margin(r= 6))) +
  scale_y_continuous(name="On-nest aggression", limits=c(0, 70)) +
  stat_smooth(method="lm", se=FALSE)

########### Darting


darts <- glmer.nb(dart.initiated ~ wasp.weight + Explore + neutral + body_contact + Log_agg+ temp + number.wasps.on.nest + Vid_Day + log_mins_on_nest+ Big_larvae+ (1|X.2/wasp_ID), 
                  data = ons, verbose = FALSE)
summary(darts)
Anova(darts)

str(ons)
ggplot(ons, aes(Log_agg, dart.initiated)) +
  geom_point() +
  xlab("Log-Aggression") + ylab("Darting") +
  theme_mine2() +
  theme(axis.title.x = element_text(margin=margin(t= 6)), #add margin to x-axis title
        axis.title.y = element_text(margin=margin(r= 6))) +
  scale_y_continuous(name="Darting", limits=c(0, 105)) +
  stat_smooth(method="lm", se=FALSE)

ggplot(ons, aes(wasp.weight, dart.initiated)) +
  geom_point() +
  xlab("Body mass") + ylab("Darting") +
  theme_mine2() +
  theme(axis.title.x = element_text(margin=margin(t= 6)), #add margin to x-axis title
        axis.title.y = element_text(margin=margin(r= 6))) +
  scale_y_continuous(name="Darting", limits=c(0, 105)) +
  stat_smooth(method="lm", se=FALSE)


########### Trophallaxis


trophs<- glmer.nb(all_troph_interactions_positive ~
                    wasp.weight + Explore + neutral + body_contact + Log_agg+ temp + number.wasps.on.nest + Vid_Day + log_mins_on_nest+ Big_larvae+ (1|X.2/wasp_ID), 
                  data = ons, verbose = FALSE) 
summary(trophs)
Anova(trophs)


##########################################################################################
############# Dominance on nests 
##########################################################################################
library(gee)
library(geepack)
##### get only wasps who are dominant or subordinate on nests
att_nests_domsub =
  att_nests %>%
  filter(!(dominance == "a")) %>%
  mutate(dominance_num = case_when(
    dominance == "d" ~ "1",
    dominance == "s" ~ "0"
  )) %>%
  mutate(Dominance_Rank = case_when(
    dominance == 'd' ~ "Dominant", #### Added for ease of graphing
    dominance == 's' ~ "Subordinate"
  )) %>%
  as.data.frame()

att_nests_domsub$dominance_num <- as.numeric(att_nests_domsub$dominance_num)
att_nests_domsub$NestID <- as.factor(att_nests_domsub$NestID)

mod <- geeglm(dominance_num ~ Explore +neutral + Weight+ body_contact  , 
              id = NestID, data = att_nests_domsub, family=binomial)

summary(mod)

check_model(mod)
head(att_nests_domsub)


ggplot(att_nests_domsub, aes(x=as.factor(Dominance_Rank), y=Explore)) + 
  geom_boxplot(alpha= 0.8) + 
  theme_mine2() +
  theme(axis.title.x = element_text(margin=margin(t= 6)), #add margin to x-axis title
        axis.title.y = element_text(margin=margin(r= 6))) +
  xlab("Dominance Rank") + ylab("Exploration")


###########################################################################################
############ Nest size and personality 
###########################################################################################

cells_all <- glmer(Cell_Count ~  Explore + Log_agg + neutral + body_contact + Weight + max_stable_group  + (1| NestID), data = att_nests, family = poisson)
summary(cells_all)
Anova(cells_all)

check_model(cells_all)


###### just wasps who are dominant and alone - should be responsible for most reproduction, 

att_nest_dom =
  att_nests %>%
  filter(Cell_Count > 0) %>% ###### remove wasp that chose site, but had not built a nest
  filter(dominance == "a" | dominance == "d") %>%
  filter(ID != "MN15") %>% ###### need to get rid of MN15 - lost nest due to usurpation.
  as.data.frame()


cells_dom <- glm(Cell_Count ~  Explore + Log_agg + neutral + body_contact + Weight + max_stable_group, data = att_nest_dom, family = poisson)
summary(cells_dom)
Anova(cells_dom)


check_model(cells_dom)



ggplot(att_nest_dom, aes(body_contact, Cell_Count)) +
  geom_point() +
  xlab("Affiliation") + ylab("Cell count") +
  theme_mine2() +
  theme(axis.title.x = element_text(margin=margin(t= 6)), #add margin to x-axis title
        axis.title.y = element_text(margin=margin(r= 6))) +
  scale_y_continuous(name="Cell count", limits=c()) +
  stat_smooth(method="lm", se=TRUE)

ggplot(att_nest_dom, aes(Log_agg, Cell_Count)) +
  geom_point() +
  xlab("Log-Aggression") + ylab("Cell count") +
  theme_mine2() +
  theme(axis.title.x = element_text(margin=margin(t= 6)), #add margin to x-axis title
        axis.title.y = element_text(margin=margin(r= 6))) +
  scale_y_continuous(name="Cell count", limits=c()) +
  stat_smooth(method="lm", se=TRUE)





