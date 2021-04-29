##############################################################################################################################################
############################################################# Rscript for: ###################################################################
##############################################################################################################################################
################################################# Zohner et al. 2021 Journal of Ecology ######################################################
#### How changes in spring and autumn phenology translate into growth - experimental evidence of asymmetric effects ##########################
##############################################################################################################################################
##############################################################################################################################################
################################ This script creates Figure S3 ###############################################################################
##############################################################################################################################################



#required packages
require(ggplot2)
require(dplyr)
require(data.table)
require(lubridate)
require(tidyverse)
require(gmodels)
require(wesanderson)



##############################################################################################################################################
##############################################################################################################################################



################
## Plot theme ##
################



plotTheme = theme(legend.position   = "None",
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



#################
## Import data ##
#################



# define paths
setwd(".../Zohner_et_al_JEcol_Analysis_data_code")
data.dir     = "Data"
output.dir   = "Output"

#load data
clim.data    = data.frame(fread(paste(data.dir,"S1_climate_data.csv",sep="/"), header=T, sep=",")) # Climate 
pheno.data   = data.frame(fread(paste(data.dir,"S2_phenology_data.csv",sep="/"))) %>%                   # Phenology
  mutate(Species = gsub(' [A-z ]*', '' , Species)) %>% #delete species epithet 
  dplyr::select(-c(Senescence_2017, Senescence_2018, Senescence_2019))



##############################################################################################################################################
##############################################################################################################################################



####################
## Data wrangling ##
####################



#Climate data
#############

clim.data = clim.data %>%
  mutate(
    #add outside temperature to inside temperature column for summer period
    inside_tmp = ifelse(is.na(inside_tmp), outside_tmp, inside_tmp),
    #as date
    date       = dmy(date),
    #date to date of year (Jan 1 = day 1)
    DOY        = yday(date),
    #degree days >5°C
    GDD_inside  = ifelse(inside_tmp<5, 0, inside_tmp-5),
    GDD_outside = ifelse(outside_tmp<5, 0, outside_tmp-5),
    )


#Phenology data
###############

pheno.data = 
  #reshape phenology table to long format
  pivot_longer(pheno.data, -c(Species, ID, Treatment), names_to = "pheno_year", values_to = "date") %>%
  mutate(
    #year identifier (keep only numbers in string)
    year        = readr::parse_number(as.character(pheno_year)), 
    #Add Individual identifier
    Individual  = paste(Species, ID, year, sep='_'),
    #as date
    date        = dmy(date),
    #date to date of year (Jan 1 = day 1)
    DOY         = yday(date) ) %>%
  #delete columns
  dplyr::select(-c(pheno_year, date)) %>%
  #delete NAs
  na.omit() 



##############################################################################################################################################
##############################################################################################################################################



###########################
## Compare chilling days ##
###########################



#data wrangling
chill.data = 
  #long format
  pivot_longer(clim.data, c(inside_tmp, outside_tmp), names_to = "type", values_to = "temperature") %>%
  filter(DOY > 273 | DOY < 90) %>%
  #get chill days (<11°C)
  group_by(year, treatment, type) %>%
  mutate(ChillDay = ifelse(temperature<11,1,0)) %>%
  summarize(ChillSum = sum(ChillDay, na.rm = TRUE))

#short format
chill.data  = pivot_wider(chill.data, names_from = type, values_from = ChillSum) %>%
  mutate(ChillReduction = inside_tmp / outside_tmp) 
  
Control      = sum(chill.data $outside_tmp) / 3
Spring       = (sum(chill.data[chill.data$treatment=='spring',]$inside_tmp)+sum(chill.data[chill.data$treatment=='autumn',]$outside_tmp)) / 3
Autumn       = (sum(chill.data[chill.data$treatment=='spring',]$outside_tmp)+sum(chill.data[chill.data$treatment=='autumn',]$inside_tmp)) / 3
AutumnSpring = (sum(chill.data[chill.data$treatment=='spring',]$inside_tmp)+sum(chill.data[chill.data$treatment=='autumn',]$inside_tmp)) / 3

Spring/Control       # 0.90%
Autumn/Control       # 0.87%
AutumnSpring/Control # 0.77%



##############################################################################################################################################
##############################################################################################################################################



############################################
## Extract degree-days (>5°C) to leaf-out ##
############################################



#Unique individuals for loop
Individuals = unique(pheno.data$Individual)

#final dataframe
pheno.final = data.frame()

#loop
for(Individual in Individuals) {
 
  #Phenology subset
  pheno.sub = pheno.data[pheno.data$Individual == Individual,]
  
  #climate subset
  clim.sub = clim.data %>%
    filter(year == pheno.sub$year,
           DOY < pheno.sub$DOY)
  
  #get degree days
  if(pheno.sub$Treatment %in% c('control','autumn')) {
    pheno.sub = pheno.sub %>%
      mutate(GDD = sum(clim.sub$GDD_outside) )
  }
  if(pheno.sub$Treatment %in% c('spring','autumn-spring')) {
    pheno.sub = pheno.sub %>%
      mutate(GDD = sum(clim.sub$GDD_inside) )
  }
  
  pheno.final = rbind(pheno.final, pheno.sub)
}



##############################################################################################################################################
##############################################################################################################################################



##########
## Plot ##
##########



#Ordering of factors
Treatments     = c("control","autumn", "spring", "autumn-spring") #order treatments
pheno.final$Treatment = factor(pheno.final$Treatment, levels=Treatments, ordered=T)
Species         = c("Lonicera", "Quercus", "Fagus") #order species
pheno.final$Species   = factor(pheno.final$Species, levels=Species, ordered=T)

#Plot
FigS3 = ggplot(data = pheno.final, mapping = aes(x = Treatment, y = GDD, fill=Treatment)) +
  stat_summary(fun.y = mean, 
               geom = "bar", color="black") + 
  stat_summary(fun.data = mean_se,  
               geom = "errorbar", width=0.2) + 
  geom_hline(yintercept=0) +
  coord_cartesian(ylim=c(23.5,500))+
  xlab("Treatment") +
  ylab("Degree-days to leaf-out") +
  scale_fill_manual(values=c("black", "orange", "green3", "red3")) +
  facet_grid(~Species) + 
  plotTheme


#Save PDF
#########

pdf(paste(output.dir,"FigS3_DegreeDay.pdf",sep="/"), width=5, height=4, useDingbats=FALSE)
FigS3
dev.off()



##############################################################################################################################################
#############################################################THE END##########################################################################
##############################################################################################################################################


