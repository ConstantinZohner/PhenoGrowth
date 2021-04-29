##############################################################################################################################################
############################################################# Rscript for: ###################################################################
##############################################################################################################################################
################################################# Zohner et al. 2021 Journal of Ecology ######################################################
#### How changes in spring and autumn phenology translate into growth - experimental evidence of asymmetric effects ##########################
##############################################################################################################################################
##############################################################################################################################################
################################ This script creates Figures 1-5, S2, and S5 #################################################################
##############################################################################################################################################

  

#required packages
require(ggplot2)
require(dplyr)
require(patchwork)
require(lubridate)
require(data.table)
require(git2r)
require(stringr)
require(broom)
require(viridis)
require(car)
require(tidyverse)



##############################################################################################################################################
##############################################################################################################################################

  

############### 
## Load data ##
###############  



# define paths
setwd(".../Zohner_et_al_JEcol_Analysis_data_code")
data.dir   = "Data"
output.dir = "Output"

#load data
clim.data    = data.frame(fread(paste(data.dir,"S1_climate_data.csv",sep="/")))     # Climate 
pheno.data   = data.frame(fread(paste(data.dir,"S2_phenology_data.csv",sep="/")))   # Phenology
volume.data  = data.frame(fread(paste(data.dir,"S3_stem_volume_data.csv",sep="/"))) # Stem volume
biomass.data = data.frame(fread(paste(data.dir,"S4_biomass_data.csv",sep="/")))     # Biomass



##############################################################################################################################################
##############################################################################################################################################



########################
## Define plot themes ##
########################



plotTheme = theme(legend.position   = "right",
                  legend.background = element_rect(fill=NA, size=0.5, linetype="solid"),
                  legend.text       = element_text(color="black"),
                  panel.grid.major  = element_blank(),
                  panel.grid.minor  = element_blank(),
                  panel.background  = element_rect(colour = NA, size=1, fill="grey97"),
                  axis.line         = element_line(color = "black"),
                  strip.background  = element_rect(fill=NA),
                  strip.text        = element_text(colour = 'black'))

plotTheme2 = theme(legend.position   = "none",
                   legend.background = element_rect(fill=NA, size=0.5, linetype="solid"),
                   legend.text       = element_text(color="black"),
                   panel.grid.major  = element_blank(),
                   panel.grid.minor  = element_blank(),
                   panel.background  = element_rect(colour = NA, size=1, fill="grey97"),
                   axis.text         = element_blank(),
                   axis.ticks        = element_blank(),
                   strip.background  = element_rect(fill=NA),
                   strip.text        = element_text(colour = 'black'))

plotTheme3 = theme(legend.position   = "none",
                   legend.background = element_rect(fill=NA, size=0.5, linetype="solid"),
                   legend.text       = element_text(color="black"),
                   panel.grid.major  = element_blank(),
                   panel.grid.minor  = element_blank(),
                   panel.background  = element_rect(colour = NA, size=1, fill="grey97"),
                   axis.line         = element_line(color = "black"),
                   strip.background  = element_rect(fill=NA),
                   strip.text        = element_text(colour = 'black'))

plotTheme4 = theme(legend.position   = "right",
                   legend.background = element_rect(fill=NA, size=0.5, linetype="solid"),
                   legend.text       = element_text(color="black"),
                   panel.grid.major  = element_blank(),
                   panel.grid.minor  = element_blank(),
                   panel.background  = element_rect(colour = NA, size=1, fill="grey97"),
                   axis.text.x       = element_text(angle = 45, hjust = 1),
                   axis.line.y       = element_line(color = "black"),
                   strip.background  = element_rect(fill=NA),
                   strip.text        = element_text(colour = 'black'))



##############################################################################################################################################
##############################################################################################################################################



######################  
## Data preparation ##	
######################
 


# Climate data
##############

#Convert dates to DOY
clim.data$date = dmy(clim.data$date)

#get median deviations from control + interquartile ranges
clim.data1 = data.frame(clim.data %>%
                          filter(!is.na(inside_tmp)) %>% #delete summer period
                          mutate(diffOutIn  = inside_tmp-outside_tmp)) %>% #add deviation column
  group_by(treatment) %>%
  summarize(medianDiff = median(diffOutIn),
            IQrange    = IQR(diffOutIn, na.rm = FALSE, type = 7))

#get new dataframe for deviation plot
clim.data2 = data.frame(clim.data %>%
                          filter(!is.na(inside_tmp)) %>% #delete NAs
                          mutate(diffOutIn  = inside_tmp-outside_tmp)) #get deviation from control

#reshape initial dataframe
#clim.data3 = melt(clim.data, id.vars = c("date","treatment"), measure.vars = c("outside_tmp","inside_tmp"))
clim.data3 = pivot_longer(clim.data, c(outside_tmp, inside_tmp), names_to = "variable", values_to = "value")

#create treatment column
clim.data3$variable2 = ifelse(clim.data3$variable=="outside_tmp", "control", 
                              ifelse(clim.data$treatment=="autumn", "autumn", "spring"))


##############################################################################################################################################


# Stem volume
#############

#order treatments
Treatment = c("spring", "autumn", "autumn-spring", "control")
volume.data$Treatment = factor(volume.data$Treatment, levels=Treatment)

#calculate volumes
volume.data$Volume = 1/3 * pi * volume.data$Height *
  ((volume.data$StemDiameter/2)^2 + 
     (volume.data$StemDiameter/2)*(volume.data$StemDiameterTop/2) + 
     volume.data$StemDiameterTop^2)

#sum stem and twig volume measurements within individuals for each year
VariablesVector = c("Volume")
volume.data = data.frame(volume.data %>% 
                           group_by(Species,ID,Treatment,Year) %>% 
                           summarize_at(VariablesVector, sum, na.rm = TRUE))

#sample sizes
table(volume.data$Treatment, volume.data$Species)/4

  
##############################################################################################################################################


# Phenology data
################
  
#get RMFs
biomass.data$RSR = biomass.data$Root/biomass.data$Shoot

#delete species epithet in all tables
DeleteEpithet <- function(df) {
  df$Species = gsub(' [A-z ]*', '' , df$Species)
  return(df)
}
dfnames    = list(pheno.data,biomass.data,volume.data)
dfs        = lapply(dfnames, DeleteEpithet)
names(dfs) = c("pheno.data","biomass.data","volume.data")

#Ordering of factors
Treatments      = c("control", "autumn", "autumn-spring", "spring") #order treatments
dfs$pheno.data$Treatment = factor(pheno.data$Treatment, levels=Treatments, ordered=T)
Species         = c("Quercus", "Fagus", "Lonicera") #order species
dfs$pheno.data$Species   = factor(dfs$pheno.data$Species, levels=Species, ordered=T)

#reshape phenology table to long format
pheno.data = pivot_longer(dfs$pheno.data, -c(Species, ID, Treatment), names_to = "variable", values_to = "value")

#Convert dates to DOY
pheno.data$value = dmy(pheno.data$value)
pheno.data$DOY   = yday(pheno.data$value) # Jan 1 = day 1 

#create year and phenophase columns
pheno.data$Year  = str_sub(pheno.data$variable, str_length(pheno.data$variable)-3, str_length(pheno.data$variable))
pheno.data$Pheno = str_sub(pheno.data$variable, 1, str_length(pheno.data$variable)-5)
pheno.data       = select(pheno.data, -c(variable, value))
pheno.data       = na.omit(pheno.data)

#from long (many rows) to short (multiple columns) format
pheno.data = pivot_wider(pheno.data, names_from = Pheno, values_from = DOY)

#get growing season length information (leaf-out to senescence)
pheno.data$GSL   = pheno.data$Senescence - pheno.data$Leaf_out

#constrain late leafers
pheno.data$Leaf_out = ifelse(pheno.data$Leaf_out>140,140,pheno.data$Leaf_out)


##############################################################################################################################################


# Biomass data
##############
  
#Melt data to a long format
biomass.data = pivot_longer(dfs$biomass.data, -c(Species, ID, Treatment), names_to = "organ", values_to = "biomass")

#Ordering of factors
Treatments     = c("control", "autumn", "autumn-spring", "spring") #order treatments
biomass.data$Treatment = factor(biomass.data$Treatment, levels=Treatments, ordered=T)
Species        = c("Quercus", "Fagus", "Lonicera") #order species
biomass.data$Species   = factor(biomass.data$Species, levels=Species, ordered=T)
variable       = c("Total", "Shoot", "Root", "RSR") #order organ
biomass.data$organ     = factor(biomass.data$organ, levels=variable, ordered=T)



##############################################################################################################################################
##############################################################################################################################################



##############  
## Analysis ##	
##############



################
# 1. Phenology #
################


####################################
# 1.1. Treatment effect on phenology
####################################


########################################
# 1.1.1. Summarize annual phenology data
########################################


# get average phenology anomalies per individual
################################################

#reshape phenology columns to long format
pheno.data = pivot_longer(pheno.data, -c(Species, ID, Treatment, Year), names_to = "phenology", values_to = "day")

#get means per species and year
phenoMean = pheno.data %>% 
  group_by(Species,phenology,Year) %>% 
  summarize(Mean = mean(day, na.rm=TRUE),
            SD = sd(day, na.rm=TRUE))

#merge data
pheno.data <-merge(pheno.data, phenoMean, all.x=T, by=c("Species","phenology","Year"))

#get change in phenology relative to mean
pheno.data = pheno.data %>% 
  mutate(phenologyAnomaly = ifelse(phenology=="Leaf_out", Mean-day, day-Mean)) %>% #add column
  dplyr::select(!c(Mean)) #delete columns

#from long (many rows) to short (multiple columns) format
pheno.data3 = pivot_wider(pheno.data %>% select(Species, ID, Treatment, Year, phenology, phenologyAnomaly), 
                          names_from = phenology, values_from = phenologyAnomaly)

#summarize years
VariablesVector = c("Leaf_out","Senescence","GSL")
pheno.data3 = data.frame(pheno.data3 %>%
                           filter(!Year %in% c("2017")) %>% #delete rows based on condition
                           group_by(Species,ID,Treatment) %>% 
                           summarize_at(VariablesVector, mean, na.rm = TRUE))


# Add biomass information
#########################

#merge phenology and biomass tables
pheno.data3 = merge(pheno.data3, dfs$biomass.data, all.x=T, by=c("Species","ID","Treatment"))

#reshape phenology columns to long format
pheno.data3 = pivot_longer(pheno.data3, -c(Species, ID, Treatment, Total, Root, Shoot, RSR), 
                           names_to = "phenology", values_to = "phenologyAnomaly")
#Ordering of factors
phenology     = c("Leaf_out", "Senescence", "GSL") #order treatments
pheno.data3$phenology = factor(pheno.data3$phenology, levels=phenology, ordered=T)


#get average phenology dates per individual
###########################################

#from long (many rows) to short (multiple columns) format
pheno.data1 = pivot_wider(pheno.data %>% select(Species, ID, Treatment, Year, phenology, day), 
                          names_from = phenology, values_from = day)

#summarize years
VariablesVector = c("Leaf_out","Senescence","GSL")
pheno.data1 = data.frame(pheno.data1 %>%
                           filter(!Year %in% c("2017")) %>% #delete rows based on condition
                           group_by(Species,ID,Treatment) %>% 
                           summarize_at(VariablesVector, mean, na.rm = TRUE))

#merge phenology and biomass tables
pheno.data1 = merge(pheno.data1, dfs$biomass, all.x=T, by=c("Species","ID","Treatment"))

#reshape phenology columns to long format
pheno.data1 = pivot_longer(pheno.data1, -c(Species, ID, Treatment, Total, Root, Shoot, RSR), 
                           names_to = "phenology", values_to = "day")


##############################################################################################################################################


##################################
# 1.1.2. Multivariate linear model
##################################


#Create spring and autumn treatments
####################################

pheno.data4        = pheno.data3
pheno.data4$spring = ifelse(pheno.data4$Treatment %in% c("control","autumn"), "no","warming")
pheno.data4$autumn = ifelse(pheno.data4$Treatment %in% c("control","spring"), "no","warming")

#Multivariate linear model
##########################

resultsLM.pheno = pheno.data4 %>% 
  group_by(Species,phenology) %>% 
  do({model = lm(phenologyAnomaly ~ spring*autumn, data=.)  # create your model
  data.frame(tidy(model))}) %>%           # get coefficient info
  filter(!term %in% c("(Intercept)")) %>% # delete rows
  dplyr::select(Species, phenology, term, p.value, estimate, std.error) %>% #delete columns
  mutate_if(is.numeric, round, 2) %>%   
  mutate(sig = ifelse(p.value<0.05, "**", ifelse(p.value<=0.1, "*", ""))) #add significance
#Ordering of factors
term                 = c("springwarming", "autumnwarming", "springwarming:autumnwarming") #order treatments
resultsLM.pheno$term = factor(resultsLM.pheno$term, levels=term, ordered=T)
Species                 = c("Quercus", "Fagus", "Lonicera") #order species
resultsLM.pheno$Species = factor(resultsLM.pheno$Species, levels=Species, ordered=T)
#order statistics table to match dataframe
resultsLM.pheno = resultsLM.pheno[with(resultsLM.pheno, order(resultsLM.pheno$Species, resultsLM.pheno$phenology)), ]
data.frame(resultsLM.pheno[,-c(7)])


##############################################################################################################################################


#############################################################
# 1.1.3. Check for deviation from control (one-sample T-test)
#############################################################


#get means of control group
pheno.data2 = pheno.data1 %>% 
  filter(Treatment %in% c("control")) %>%
  group_by(Species,phenology) %>% 
  summarize(Mean = mean(day, na.rm=TRUE))
#merge data
pheno.data2 <-merge(pheno.data1[!pheno.data1$Treatment=="control",], 
                    pheno.data2, all.x=T, by=c("Species","phenology"))
#get change in phenology relative to control
pheno.data2$phenologyChange = pheno.data2$day - pheno.data2$Mean

#Ordering of factors
Phenology     = c("Leaf_out", "Senescence", "GSL") #order treatments
pheno.data2$phenology = factor(pheno.data2$phenology, levels=Phenology, ordered=T)
Treatments     = c("spring", "autumn", "autumn-spring") #order treatments
pheno.data2$Treatment = factor(pheno.data2$Treatment, levels=Treatments, ordered=T)
Species         = c("Quercus", "Fagus", "Lonicera") #order species
pheno.data2$Species   = factor(pheno.data2$Species, levels=Species, ordered=T)

#T-test
resultsTT = pheno.data2 %>% 
  group_by(Species,phenology,Treatment) %>% 
  do({model = t.test(.$phenologyChange) # create your model
  data.frame(tidy(model),
             tidy(shapiro.test(x = .$phenologyChange))[1,2])}) %>%           
  dplyr::select(Species, phenology, Treatment, p.value, estimate, p.value.1) %>% #delete columns
  rename(shapiroTest=p.value.1) %>%      
  mutate(significance = ifelse(p.value<0.05, "**", ifelse(p.value<=0.1, "*", ""))) %>% #add significance
  mutate_if(is.numeric, round, 2) 
#order statistics table to match dataframe
resultsTT = resultsTT[with(resultsTT, order(resultsTT$Species, resultsTT$phenology)), ]


##############################################################################################################################################


###################################################################
# 1.2. Phenology effect on total biomass (Univariate linear models)
###################################################################


#reshape biomass columns to long format
pheno.data3 = pivot_longer(pheno.data3, -c(Species, ID, Treatment, phenology, phenologyAnomaly), 
                           names_to = "organ", values_to = "biomass")

#Transform biomass to percentage anomaly
########################################

#get means per species
pheno.data5 = pheno.data3 %>% 
  group_by(Species,organ) %>% 
  summarize(Mean = mean(biomass, na.rm=TRUE))
#merge data
pheno.data3 <-merge(pheno.data3, pheno.data5, all.x=T, by=c("Species","organ"))
#get precentage change in biomass
pheno.data3$relativeBiomass = (pheno.data3$biomass / pheno.data3$Mean -1) * 100

#order species
Species             = c("Quercus", "Fagus", "Lonicera") 
pheno.data3$Species = factor(pheno.data3$Species, levels=Species, ordered=T)

#Extract linear model info
resultsLM.pheno2 = pheno.data3 %>% 
  group_by(Species,phenology,organ) %>% 
  do({model = lm(relativeBiomass ~ phenologyAnomaly, data=.)  # create model
  data.frame(tidy(model),                # get coefficient info
             glance(model),
             tidy(shapiro.test(x = residuals(object = model)))[1,2]
  )}) %>%                     # get model info
  filter(term != "(Intercept)") %>%      # delete intercept info
  dplyr::select(Species, phenology, organ, estimate, std.error, p.value, r.squared, p.value.2) %>% # delete columns
  rename(shapiro=p.value.2) %>%          #rename columns
  mutate(sig = ifelse(p.value<0.05, "**", ifelse(p.value<=0.1, "*", ""))) %>% # add column
  mutate_if(is.numeric, round, 2)

#Order phenophases
Organs            = c("Total","Shoot","Root","RSR") #order treatments
resultsLM.pheno2$organ   = factor(resultsLM.pheno2$organ, levels=Organs, ordered=T)
#order statistics table to match dataframe
resultsLM.pheno2 = resultsLM.pheno2[with(resultsLM.pheno2, 
                                         order(resultsLM.pheno2$Species, 
                                               resultsLM.pheno2$organ, 
                                               resultsLM.pheno2$phenology)), ]
data.frame(resultsLM.pheno2)

#Visually inspect model assumptions 
data.assumptions = pheno.data3
data.assumptions$category = paste(data.assumptions$Species,
                                  data.assumptions$organ,
                                  data.assumptions$phenology, sep="_") #create identifier column
category.list = as.factor(unique(data.assumptions$category))       #create category vector
par(mfrow=c(2,2))                                       #set plot layout
for (category in category.list){                        #loop over categories
  tab.subset=data.assumptions[data.assumptions$category==category, ] #create table subset
  plot(lm(relativeBiomass ~ phenologyAnomaly, data=tab.subset), main=category)
}



##############################################################################################################################################
##############################################################################################################################################



##############
# 2. Biomass #
##############



####################################################
# 2.1. Total and per treatment biomass distributions
####################################################
  

#Species-level plot
ggplot(data = biomass.data, mapping = aes(biomass, fill=organ)) +
  geom_histogram(bins=10, position="identity", color="black") +
  xlab("Biomass (g) [panels 1-3]                     
                                                                                           Biomass ratio [panel 4]") +
  ylab("Count") +
  scale_fill_viridis(option="E",discrete=TRUE) +
  facet_grid(Species~organ, scales="free") + 
  plotTheme4

#Treatment-level plot
ggplot(data = biomass.data[biomass.data$organ %in% c("Shoot","Root"),], mapping = aes(biomass, fill=organ)) +
  geom_histogram(bins=10, color="black") +
  xlab("Biomass (g)") +
  ylab("Count") +
  scale_fill_viridis(option="C",discrete=TRUE) +
  facet_grid(organ+Species~Treatment, scales="free") + 
  plotTheme4


##############################################################################################################################################


#########################
# 2.2. Treatment analysis
#########################


#######
# ANOVA
#######

#Extract model info and check ANOVA assumptions
resultsLM = biomass.data %>% 
  group_by(Species,organ) %>% 
  do({model = aov(biomass ~ Treatment, data=.)        # create your model
  data.frame(tidy(model)[1,],                         # get coefficient info
             glance(model)[1,],                       # get model info
             leveneTest=tidy(leveneTest(model))[1,4], # check Homogeneity of variances
             tidy(shapiro.test(x = residuals(object = model)))[1,2] # check for Normality
  )}) %>%           
  #select(Species, organ, p.value, r.squared, adj.r.squared, p.value.2, p.value.3) %>% #delete columns
  rename(leveneTest=p.value.2, shapiroTest=p.value.3) %>%                             #rename columns
  mutate(sig = ifelse(p.value<0.05, "**", ifelse(p.value<0.1, "*", "")),
         Treatment="control")     #add significance and dummy column for plotting
data.frame(resultsLM)

# Plots to check ANOVA assumptions
biomass.data$category = paste(biomass.data$Species,biomass.data$organ, sep="_") #create identifier column
category.list = as.factor(unique(biomass.data$category))        #create category vector
par(mfrow=c(2,2))                                       #set plot layout
for (category in category.list){                        #loop over categories
  tab.subset=biomass.data[biomass.data$category==category, ]            #create table subset
  plot(aov(biomass ~ Treatment, data=tab.subset), 1, main=category) # 1. Homogeneity of variances
  plot(aov(biomass ~ Treatment, data=tab.subset), 2, main=category) # 2. Normality
}


#####
#Plot
#####

Boxp <- ggplot(data = biomass.data, mapping = aes(x = Treatment, y = biomass, fill=Treatment)) +
  geom_boxplot(position=position_dodge(0.8),coef=1e30) +
  geom_text(data = resultsLM,
            mapping = aes(x = Inf, y = Inf, hjust = 2, vjust = 2, 
                          label = paste("R2 = ", round(r.squared,2), sig, sep="")), 
            size=3.5, color="black")+
  xlab("Treatment") +
  ylab("Biomass (g)") +
  scale_fill_viridis(option="E",discrete=TRUE) +
  facet_grid(organ~Species, scales="free") + 
  plotTheme4
Boxp


##############################################################################################################################################


################################
# 2.3. Relative biomass analysis
################################


###########################
# Multivariate linear model
###########################


#Create spring and autumn treatments
####################################

biomass.data$spring = ifelse(biomass.data$Treatment %in% c("control","autumn"), "no","warming")
biomass.data$autumn = ifelse(biomass.data$Treatment %in% c("control","spring"), "no","warming")


#Transform biomass to percentage anomaly
########################################

#get means per species
biomass.data1 = biomass.data %>% 
  group_by(Species,organ) %>% 
  summarize(Mean = mean(biomass, na.rm=TRUE))
#merge data
biomass.data1 <-merge(biomass.data, biomass.data1, all.x=T, by=c("Species","organ"))
#get precentage change in biomass
biomass.data1$relativeBiomass = (biomass.data1$biomass / biomass.data1$Mean -1) * 100


#Multivariate linear model
##########################

resultsLM.biomass = biomass.data1 %>% 
  group_by(Species,organ) %>% 
  do({model = lm(relativeBiomass ~ spring*autumn, data=.)  # create your model
  data.frame(tidy(model))}) %>%           # get coefficient info
  filter(!term %in% c("(Intercept)")) %>% # delete rows
  dplyr::select(Species, organ, term, p.value, estimate, std.error) %>% #delete columns
  mutate_if(is.numeric, round, 2) %>%   
  mutate(sig = ifelse(p.value<0.05, "**", ifelse(p.value<=0.1, "*", ""))) #add significance
data.frame(resultsLM.biomass)


########
# T-test
########

#get means of control group
###########################

biomass.data1 = biomass.data %>% 
  filter(Treatment %in% c("control")) %>%
  group_by(Species,organ) %>% 
  summarize(Mean = mean(biomass, na.rm=TRUE))
#merge data
biomass.data1 <-merge(biomass.data[!biomass.data$Treatment=="control",], biomass.data1, all.x=T, by=c("Species","organ"))
#get precentage change in biomass
biomass.data1$relativeBiomass = ifelse(biomass.data1$organ!="RSR", (biomass.data1$biomass / biomass.data1$Mean -1) * 100,
                                       (biomass.data1$biomass - biomass.data1$Mean)*100)
#Ordering of factors
Treatments     = c("spring", "autumn", "autumn-spring") #order treatments
biomass.data1$Treatment = factor(biomass.data1$Treatment, levels=Treatments, ordered=T)
Species     = c("Quercus", "Fagus", "Lonicera") #order treatments
biomass.data1$Species = factor(biomass.data1$Species, levels=Species, ordered=T)

#T-test
#######

resultsTTbiomass = biomass.data1 %>% 
  group_by(Species,organ,Treatment) %>% 
  do({model = t.test(.$relativeBiomass)        # create your model
  data.frame(tidy(model),
             tidy(shapiro.test(x = .$relativeBiomass))[1,2])}) %>%           
  select(Species, organ, Treatment, p.value, estimate, p.value.1) %>% #delete columns
  rename(shapiroTest=p.value.1) %>%   
  mutate_if(is.numeric, round, 2) %>%   
  mutate(sig = ifelse(p.value<0.05, "**", ifelse(p.value<=0.1, "*", ""))) #add significance
#order statistics table to match dataframe
resultsTTbiomass = resultsTTbiomass[with(resultsTTbiomass, 
                                         order(resultsTTbiomass$Species, resultsTTbiomass$organ)), ]
data.frame(resultsTTbiomass)



##############################################################################################################################################
##############################################################################################################################################



#############
## Figures ##	
#############



###########
# Figure1 #
###########
  
 
#2D density plot
Dens2D = ggplot(clim.data, aes(x=outside_tmp, y=inside_tmp-outside_tmp) ) +
  stat_density_2d(aes(fill = ..level..), geom = "polygon") +
  scale_fill_distiller(palette=4, direction=1) +
  geom_hline(aes(yintercept=4.7))+
  coord_cartesian(ylim = c(0, 10),xlim = c(-4, 19))+
  xlab("Ambient temperature (°C)") +
  ylab("Deviation from control (°C)") +
  plotTheme

#Boxplot
boxp = ggplot(clim.data2, aes(x=treatment, y=diffOutIn, fill=treatment)) + 
  geom_boxplot(outlier.shape=NA,alpha=.6) +
  coord_flip(ylim = c(0.5, 10)) +
  xlab("") +
  ylab("") +
  scale_fill_manual(values=c("orange","green4"))+
  plotTheme2

#Histogram
histo = ggplot(clim.data2, aes(x=diffOutIn, fill=treatment)) +
  geom_histogram(binwidth=.5, alpha=.6, position="identity") +
  geom_vline(data=clim.data1, aes(xintercept=medianDiff, colour=c("green4","orange3"),alpha=.8),
             linetype="dashed", size=1)+
  xlab("Deviation from control (°C)") +
  ylab("Count (number of days)") +
  scale_fill_manual(values=c("orange","green4"))+
  scale_colour_manual(values=c("orange","green3"))+
  coord_cartesian(xlim = c(0.5, 10), ylim = c(3, 60))+
  plotTheme3

#Time plot

#extract date information
phenoMean$dateMin = as.Date(phenoMean$Mean-phenoMean$SD, origin = paste0(phenoMean$Year,"-01-01"))
phenoMean$dateMax = as.Date(phenoMean$Mean+phenoMean$SD, origin = paste0(phenoMean$Year,"-01-01"))
phenoMean$Y = ifelse(phenoMean$Species=="Lonicera",30,ifelse(phenoMean$Species=="Quercus",28,26))

TempPlot = ggplot() +
  
  geom_rect(data=phenoMean[phenoMean$phenology=='Leaf_out',], 
            mapping=aes(xmin=dateMin, xmax=dateMax, ymin=-Inf, ymax=Y),
            fill=rep(c('green4','green1',"green3"),each=3), alpha=0.4) +
  
  geom_rect(data=phenoMean[phenoMean$phenology=='Senescence',], 
            mapping=aes(xmin=dateMin, xmax=dateMax, ymin=-Inf, ymax=Y),
            fill=rep(c('orange4','orange1',"orange3"),each=3), alpha=0.4) +
  
  geom_line(data=clim.data3, aes(x = date, y= value, group = variable, colour = variable2)) +
  
  geom_segment(data=phenoMean[phenoMean$phenology=='Leaf_out',], 
               aes(x = dateMin, y = Y, xend = dateMax, yend = Y),
               size=3,
               col=rep(c('green4','green1',"green3"),each=3) )+
  
  geom_segment(data=phenoMean[phenoMean$phenology=='Senescence',], 
               aes(x = dateMin, y = Y, xend = dateMax, yend = Y),
               size=3,
               col=rep(c('orange4','orange1',"orange3"),each=3) )+
  
  scale_color_manual(values=c("orange3","black","green4"))+
  coord_cartesian(ylim=c(-7,29))+
  xlab("") +
  ylab("Temperature (°C)") +
  plotTheme

#define plot layout
layout <- "
AA
AA
AA
BD
CD
CD"

#Merge plots
Fig1_ClimPlot = TempPlot + boxp + histo + Dens2D + plot_layout(design = layout)
Fig1_ClimPlot

#Save PDF
#########

pdf(paste(output.dir,"Fig1_Temperature deviation.pdf",sep="/"), width=9, height=6, useDingbats=FALSE)
Fig1_ClimPlot
dev.off()


##############################################################################################################################################

  
###########
# Figure2 #
###########


#create panel labels
dat_text       = data.frame(
  label        = c("a","b","c","d","e","f","g","h","i"),
  phenology   = rep(c("Leaf_out","Senescence","GSL"),3),
  Species      = c(rep("Quercus",3), rep("Fagus",3), rep("Lonicera",3)),
  term         = "springwarming")

#Plot
Fig2_LMplot = ggplot(data = resultsLM.pheno, mapping = aes(x = term, y = estimate, 
                                                           fill=phenology, alpha=term)) +
  geom_bar(position=position_dodge(), stat="identity", color="black") +
  geom_errorbar(aes(ymin=estimate-std.error, ymax=estimate+std.error),
                width=.2,                    # Width of the error bars
                position=position_dodge(.9))+
  geom_hline(yintercept=0) +
  geom_text(data = dat_text, mapping = aes(x = Inf, y = Inf, 
                                           hjust = 1, vjust = 1,
                                           label = label,
                                           fontface = "bold"))+
  stat_summary(geom = 'text', label = resultsLM.pheno$sig, 
               fun.y = mean, vjust = 2) +
  scale_alpha_discrete(range = c(1, 0.7)) +
  xlab("Phenophase") +
  ylab("Effect size (days)") +
  scale_fill_manual(values = c("green3", "orange", "blue")) +
  facet_grid(Species~phenology) + 
  plotTheme4
Fig2_LMplot

#Save PDF
#########

pdf(paste(output.dir,"Fig2_PhenologyLinearTreatmentModel.pdf",sep="/"), width=7, height=7, useDingbats=FALSE)
Fig2_LMplot
dev.off()


##############################################################################################################################################


###########
# Figure3 #
###########


dfs$volume.data$Species = factor(dfs$volume.data$Species, levels=c("Quercus", "Fagus", "Lonicera"))
pd  = position_dodge(0.2) 
Fig3_LinePlot = ggplot(dfs$volume.data, aes(x=Year, y=Volume, colour=Treatment, group=Treatment)) + 
  stat_summary(fun = mean, geom="line", position=pd) +
  stat_summary(fun.data = "mean_se", size = 0.5, position=pd) +
  labs(x = "Year", y = "Stem volume (cm3)") +
  scale_color_manual(values=c("green3", "orange", "red3", "black")) +
  coord_cartesian(ylim = c(8, 152))+
  facet_grid(~Species) +
  plotTheme 
Fig3_LinePlot

#########
#Save PDF
#########

pdf(paste(output.dir,"Fig3_LinePlot.pdf",sep="/"), width=9, height=4, useDingbats=FALSE)
Fig3_LinePlot
dev.off()


##############################################################################################################################################


###########
# Figure4 #
###########


#Ordering of factors
term           = c("springwarming", "autumnwarming", "springwarming:autumnwarming")
resultsLM.biomass$term = factor(resultsLM.biomass$term, levels=term, ordered=T)
#create panel labels
dat_text       = data.frame(
  label        = c("a","b","c","d","e","f","g","h","i","j","k","l"),
  organ        = rep(c("Total","Shoot","Root","RSR"),3),
  Species      = c(rep("Quercus",4), rep("Fagus",4), rep("Lonicera",4)),
  term         = "springwarming")
#Plot
LMplot = ggplot(data = resultsLM.biomass, mapping = aes(x = term, y = estimate, 
                                                        fill=organ, alpha=term)) +
  geom_bar(position=position_dodge(), stat="identity", color="black") +
  geom_errorbar(aes(ymin=estimate-std.error, ymax=estimate+std.error),
                width=.2,                    # Width of the error bars
                position=position_dodge(.9))+
  geom_hline(yintercept=0) +
  geom_text(data = dat_text, mapping = aes(x = -Inf, y = -Inf, 
                                           hjust = -0.3, vjust = -1,
                                           label = label,
                                           fontface = "bold"))+
  stat_summary(geom = 'text', label = resultsLM.biomass$sig, 
               fun.y = mean, vjust = 2) +
  scale_fill_manual(values=c("orange","yellow2","red3","grey70"))+
  scale_alpha_discrete(range = c(1, 0.7)) +
  xlab("Phenophase") +
  ylab("Effect size (%)") +
  facet_grid(Species~organ) + 
  plotTheme4
LMplot

#Save PDF
#########

pdf(paste(output.dir,"Fig4_BiomassLinearTreatmentModel.pdf",sep="/"), width=8, height=7, useDingbats=FALSE)
LMplot
dev.off()


##############################################################################################################################################


###########
# Figure5 #
###########


#Effect sizes
#############

#create panel labels
dat_text      = data.frame(
  label       = c("a","b","c","d","e","f","g","h","i","j","k","l"),
  organ       = rep(c("Total","Shoot","Root","RSR"),3),
  Species     = c(rep("Quercus",4), rep("Fagus",4), rep("Lonicera",4)),
  phenology   = "Leaf_out")

#Plot
LMplot = ggplot(data = resultsLM.pheno2, mapping = aes(x = phenology, y = estimate, 
                                                       fill=organ, alpha=phenology)) +
  geom_bar(position=position_dodge(), stat="identity", color="black") +
  geom_errorbar(aes(ymin=estimate-std.error, ymax=estimate+std.error),
                width=.2,                    # Width of the error bars
                position=position_dodge(.9))+
  geom_hline(yintercept=0) +
  geom_text(data = dat_text, mapping = aes(x = -Inf, y = -Inf, 
                                           hjust = -0.1, vjust = -1,
                                           label = label,
                                           fontface = "bold"))+
  stat_summary(geom = 'text', label = resultsLM.pheno2$sig, 
               fun.y = mean, vjust = 2) +
  scale_fill_manual(values=c("orange","yellow2","red3","grey70"))+
  scale_alpha_discrete(range = c(1, 0.7)) +
  xlab("Phenophase") +
  ylab("Effect size (%/day)") +
  facet_grid(Species~organ) + 
  plotTheme4
LMplot

#Biomass
########

#order
Organs      = c("Total","Shoot","Root","RSR") #order treatments
pheno.data3$organ = factor(pheno.data3$organ, levels=Organs, ordered=T)
#create panel labels
dat_text       = data.frame(
  label        = c("a","b","c","d","e","f","g","h","i"),
  phenology    = rep(c("Leaf_out","Senescence","GSL"),3),
  Species      = c(rep("Quercus",3), rep("Fagus",3), rep("Lonicera",3)),
  organ        = "Root")
BiomassPlot = ggplot(pheno.data3[pheno.data3$organ %in% c("Total","Shoot","Root"),], 
                     aes(x=phenologyAnomaly, y=relativeBiomass, colour=organ)) + 
  geom_hline(yintercept=0) +
  geom_vline(xintercept=0) +
  geom_smooth(method=lm, fullrange = T, alpha = 0.2) +
  labs(x = "Phenology anomaly (days)", y = "Biomass anomaly (%)") +
  scale_color_manual(values=c("orange","yellow2","red3"))+
  geom_text(data = dat_text, mapping = aes(x = -Inf, y = Inf, 
                                           hjust = -0.2, vjust = 1.3,
                                           label = label,
                                           fontface = "bold"), color="black")+
  geom_text(data    = resultsLM.pheno2[resultsLM.pheno2$organ=="Total",],
            mapping = aes(x = Inf, y = Inf, hjust = 1.5, vjust = 2.5, 
                          label = paste(round(estimate,1), " %/day", sig, sep="")), 
            size=3.5, color="orange")+
  geom_text(data    = resultsLM.pheno2[resultsLM.pheno2$organ=="Shoot",],
            mapping = aes(x = Inf, y = Inf, hjust = 1.5, vjust = 4.5, 
                          label = paste(round(estimate,1), " %/day", sig, sep="")), 
            size=3.5, color="yellow2")+
  geom_text(data    = resultsLM.pheno2[resultsLM.pheno2$organ=="Root",],
            mapping = aes(x = Inf, y = Inf, hjust = 1.5, vjust = 6.5, 
                          label = paste(round(estimate,1), " %/day", sig, sep="")), 
            size=3.5, color="red3")+
  facet_grid(Species~phenology, scales = "free") +
  plotTheme
BiomassPlot

#RSR
RSRplot = ggplot(pheno.data3[pheno.data3$organ %in% c("RSR"),], aes(x=phenologyAnomaly, y=biomass)) + 
  geom_point(size=0.2) +  
  geom_smooth(method=lm, fullrange = F) +
  labs(x = "Day", y = "Root:shoot ratio") +
  scale_color_viridis(option="viridis",discrete=TRUE) +
  geom_text(data    = resultsLM.pheno2[resultsLM.pheno2$organ=="RSR",],
            mapping = aes(x = Inf, y = Inf, hjust = 1.5, vjust = 2.5, 
                          label = paste("Slope = ", round(estimate,3), " day-1", sig, sep="")), 
            size=3.5, color="black")+
  facet_grid(Species~phenology, scales = "free") +
  plotTheme
RSRplot

#Save PDF
#########

pdf(paste(output.dir,"Fig5_PhenologyBiomass.pdf",sep="/"), width=8, height=7, useDingbats=FALSE)
BiomassPlot
RSRplot
dev.off()


##############################################################################################################################################


#############
# Figure S2 #
#############


#create panel labels
dat_text <- data.frame(
  label       = c("a","b","c","d","e","f","g","h","i"),
  phenology   = rep(c("Leaf_out","Senescence","GSL"),3),
  Species     = c(rep("Quercus",3), rep("Fagus",3), rep("Lonicera",3)),
  Treatment   = "spring")

#Plot
FigS2_PhenologyChangePlot = ggplot(data = pheno.data2, mapping = aes(x = Treatment, y = phenologyChange, 
                                                                     fill=phenology, alpha=Treatment)) +
  stat_summary(fun.y = mean, 
               geom = "bar", color="black") + 
  stat_summary(fun.data = mean_se,  
               geom = "errorbar", width=0.2) + 
  geom_hline(yintercept=0) +
  geom_text(data = dat_text, mapping = aes(x = -Inf, y = Inf, 
                                           hjust = -0.2, vjust = 1.5,
                                           label = label,
                                           fontface = "bold"))+
  stat_summary(geom = 'text', label = resultsTT$significance, 
               fun.y = mean, vjust = 2) +
  scale_alpha_discrete(range = c(1, 0.7)) +
  xlab("Treatment") +
  ylab("Phenological change (days)") +
  scale_fill_manual(values = c("green3", "orange", "blue")) +
  facet_grid(Species~phenology) + 
  plotTheme4
FigS2_PhenologyChangePlot

#Save PDF
#########

pdf(paste(output.dir,"FigS2_PhenologyChange.pdf",sep="/"), width=7, height=7, useDingbats=FALSE)
FigS2_PhenologyChangePlot
dev.off()


##############################################################################################################################################


#############
# Figure S5 #
#############


#create panel labels
dat_text <- data.frame(
  label = c("a","b","c","d","e","f","g","h","i","j","k","l"),
  organ   = rep(c("Total","Shoot","Root","RSR"),3),
  Species     = c(rep("Quercus",4), rep("Fagus",4), rep("Lonicera",4)),
  Treatment   = "spring")

#Plot
BiomassChangePlot = ggplot(data = biomass.data1, mapping = aes(x = Treatment, y = relativeBiomass, 
                                                               fill=organ, alpha=Treatment)) +
  stat_summary(fun.y = mean, 
               geom = "bar", color="black") + 
  stat_summary(fun.data = mean_se,  
               geom = "errorbar", width=0.2) + 
  geom_hline(yintercept=0) +
  geom_text(data = dat_text, mapping = aes(x = -Inf, y = -Inf, 
                                           hjust = -0.3, vjust = -1,
                                           label = label,
                                           fontface = "bold"))+
  stat_summary(geom = 'text', label = resultsTTbiomass$sig, 
               fun.y = max, vjust = 2.9) +
  scale_alpha_discrete(range = c(1, 0.7)) +
  coord_cartesian(ylim = c(-60, 60))+
  xlab("Treatment") +
  ylab("Biomass change (%)") +
  scale_fill_manual(values=c("orange","yellow2","red3","grey70"))+
  facet_grid(Species~organ) + 
  plotTheme4
BiomassChangePlot

#Save PDF
#########

pdf(paste(output.dir,"FigS5_BiomassChange.pdf",sep="/"), width=8, height=7, useDingbats=FALSE)
BiomassChangePlot
dev.off()



##############################################################################################################################################
##############################################################################################################################################



#####################
## Reproducibility ##	
#####################


## datetime
Sys.time()
#[1] "2021-04-29 07:42:11 CEST"

## sessioninfo
sessionInfo()
#R version 3.6.2 (2019-12-12)
#Platform: x86_64-apple-darwin15.6.0 (64-bit)
#Running under: macOS  10.16

#Matrix products: default
#LAPACK: /Library/Frameworks/R.framework/Versions/3.6/Resources/lib/libRlapack.dylib

#locale:
#  [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

#attached base packages:
#  [1] stats     graphics  grDevices utils     datasets  methods   base     

#other attached packages:
#  [1] forcats_0.4.0     purrr_0.3.3       readr_1.3.1       tidyr_1.0.2       tibble_2.1.3      tidyverse_1.3.0  
#[7] car_3.0-6         carData_3.0-3     viridis_0.5.1     viridisLite_0.3.0 broom_0.5.4       stringr_1.4.0    
#[13] git2r_0.27.1      data.table_1.12.8 lubridate_1.7.4   patchwork_1.0.0   dplyr_0.8.4       ggplot2_3.2.1    

#loaded via a namespace (and not attached):
# [1] Rcpp_1.0.3         lattice_0.20-38    assertthat_0.2.1   rprojroot_1.3-2    digest_0.6.23      R6_2.4.1          
#[7] cellranger_1.1.0   plyr_1.8.5         backports_1.1.5    reprex_0.3.0       httr_1.4.1         pillar_1.4.3      
#[13] rlang_0.4.4        lazyeval_0.2.2     curl_4.3           readxl_1.3.1       rstudioapi_0.11    labeling_0.3      
#[19] desc_1.2.0         foreign_0.8-72     munsell_0.5.0      compiler_3.6.2     modelr_0.1.6       pkgconfig_2.0.3   
#[25] tidyselect_1.0.0   gridExtra_2.3      rio_0.5.16         fansi_0.4.1        crayon_1.3.4       dbplyr_1.4.2      
#[31] withr_2.1.2        MASS_7.3-53        grid_3.6.2         nlme_3.1-142       jsonlite_1.6.1     gtable_0.3.0      
#[37] lifecycle_0.1.0    DBI_1.1.0          magrittr_1.5       scales_1.1.0       zip_2.0.4          cli_2.0.1         
#[43] stringi_1.4.5      farver_2.0.3       reshape2_1.4.3     fs_1.3.1           testthat_2.3.1     xml2_1.2.2        
#[49] generics_0.0.2     vctrs_0.2.2        openxlsx_4.1.4     RColorBrewer_1.1-2 tools_3.6.2        glue_1.3.1        
#[55] hms_0.5.3          abind_1.4-5        pkgload_1.0.2      colorspace_1.4-1   rvest_0.3.5        haven_2.2.0    



##############################################################################################################################################
#############################################################THE END##########################################################################
##############################################################################################################################################


  