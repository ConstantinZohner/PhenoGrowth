##############################################################################################################################################
############################################################# Rscript for: ###################################################################
##############################################################################################################################################
################################################# Zohner et al. 2021 Journal of Ecology ######################################################
#### How changes in spring and autumn phenology translate into growth - experimental evidence of asymmetric effects ##########################
##############################################################################################################################################
##############################################################################################################################################
################################ This script creates Figure S4 ###############################################################################
##############################################################################################################################################



#required packages
require(ggplot2)
require(dplyr)
require(data.table)
require(lubridate)
require(tidyverse)



##############################################################################################################################################
##############################################################################################################################################



################
## Plot theme ##
################



plotTheme = theme(legend.position   = "right",
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
  mutate(Species = gsub(' [A-z ]*', '' , Species))#delete species epithet



##############################################################################################################################################
##############################################################################################################################################



####################
## Data wrangling ##
####################



#Climate data
#############

clim.data = clim.data %>%
  mutate(
    inside_tmp = ifelse(is.na(inside_tmp), outside_tmp, inside_tmp),
    #as date
    date       = dmy(date),
    #date to date of year (Jan 1 = day 1)
    DOY        = yday(date) ) 

#--------------------------------------------------------------------------------------------------------------------------

#Phenology data
###############

pheno.data = 
  #reshape phenology table to long format
  pivot_longer(pheno.data, -c(Species, ID, Treatment), names_to = "pheno_year", values_to = "date") %>%
  mutate(
    #year identifier (keep only numbers in string)
    year        = readr::parse_number(as.character(pheno_year)), 
    #phenophase identifier 
    Pheno       = str_sub(pheno_year, 1, str_length(pheno_year)-5),
    #Add Individual identifier
    Individual  = paste(Species, ID, year, sep='_'),
    #as date
    date        = dmy(date),
    #date to date of year (Jan 1 = day 1)
    DOY         = yday(date) ) %>%
  #delete columns
  dplyr::select(-c(pheno_year, date))

#delete NAs
pheno.data       = na.omit(pheno.data)

#from long (many rows) to short (multiple columns) format
pheno.data = pivot_wider(pheno.data, names_from = Pheno, values_from = DOY)

#--------------------------------------------------------------------------------------------------------------------------

#Extract mean temperature 28 days after leaf-out / before senescence
####################################################################

#Unique individuals for loop
Individuals = unique(pheno.data$Individual)

#final dataframe
pheno.final = data.frame()

#loop
for(Individual in Individuals) {

  #Phenology subset
  pheno.sub = pheno.data[pheno.data$Individual == Individual,]
  
  #climate subset (28 days after leaf-out or before autumn)
  
  #Spring
  clim.sub.out = clim.data %>%
    filter(year == pheno.sub$year,
           DOY >= pheno.sub$Leaf_out,
           DOY < pheno.sub$Leaf_out+28)
  #Autumn
  clim.sub.off = clim.data %>%
    filter(year == pheno.sub$year,
           DOY > pheno.sub$Senescence-28,
           DOY <= pheno.sub$Senescence)
  
  #get mean temperatures
  if(pheno.sub$Treatment %in% c('control')) {
   pheno.sub = pheno.sub %>%
     mutate(spring.climate = mean(clim.sub.out$outside_tmp),
            autumn.climate = mean(clim.sub.off$outside_tmp) )
  }
  if(pheno.sub$Treatment %in% c('spring')) {
    pheno.sub = pheno.sub %>%
      mutate(spring.climate = mean(clim.sub.out$inside_tmp),
             autumn.climate = mean(clim.sub.off$outside_tmp) )
  }
  if(pheno.sub$Treatment %in% c('autumn-spring')) {
    pheno.sub = pheno.sub %>%
    mutate(spring.climate = mean(clim.sub.out$inside_tmp),
           autumn.climate = mean(clim.sub.off$inside_tmp) )
  }
  if(pheno.sub$Treatment %in% c('autumn')) {
    pheno.sub = pheno.sub %>%
    mutate(spring.climate = mean(clim.sub.out$outside_tmp),
           autumn.climate = mean(clim.sub.off$inside_tmp) )
  }
  
  pheno.final = rbind(pheno.final, pheno.sub)
}

pheno.final = pivot_longer(pheno.final, c(spring.climate, autumn.climate), names_to = "season", values_to = "temperature")

#--------------------------------------------------------------------------------------------------------------------------

#Change in temperature relative to control
##########################################

#get means of control group
pheno.mean = pheno.final %>% 
  filter(Treatment %in% c("control")) %>%
  group_by(Species,season) %>% 
  summarize(tempMean = mean(temperature, na.rm=TRUE))

#merge data
pheno.final <-merge(pheno.final[!pheno.final$Treatment=="control",], 
                    pheno.mean, all.x=T, by=c("Species","season"))

#get change relative to control
pheno.final$TempDiff = pheno.final$temperature - pheno.final$tempMean

#--------------------------------------------------------------------------------------------------------------------------

#get means and standard error
#############################

temp.mean = pheno.final %>% 
  group_by(Treatment,Species,season) %>% 
  summarize(tempMean = mean(TempDiff, na.rm=TRUE),
            tempSE   = sd(TempDiff))

as.data.frame(temp.mean)
#Treatment      Species         season          tempMean   tempSE
#spring         Quercus         spring.climate  1.94299272 3.886219
#spring         Quercus         autumn.climate  1.27974702 1.392956
#spring         Fagus           spring.climate  1.93572019 4.374276
#spring         Fagus           autumn.climate  0.52362923 1.275064
#spring         Lonicera        spring.climate  0.68258889 2.156495
#spring         Lonicera        autumn.climate  0.66861640 1.979595
#autumn         Quercus         spring.climate -0.23684537 3.338270
#autumn         Quercus         autumn.climate  1.85652976 2.194334
#autumn         Fagus           spring.climate  1.01614056 3.857421
#autumn         Fagus           autumn.climate -0.49104669 1.128538
#autumn         Lonicera        spring.climate  0.04790079 3.182524
#autumn         Lonicera        autumn.climate  2.34309706 1.620551
#autumn-spring  Quercus         spring.climate  2.40589933 4.218064
#autumn-spring  Quercus         autumn.climate  3.58681429 3.208678
#autumn-spring  Fagus           spring.climate  2.66260643 4.053740
#autumn-spring  Fagus           autumn.climate  0.56099299 1.565854
#autumn-spring  Lonicera        spring.climate  2.68345899 2.340641
#autumn-spring  Lonicera        autumn.climate  2.56213228 1.458263



##############################################################################################################################################
##############################################################################################################################################



##########
## Plot ##
##########



#Ordering of factors
season     = c("spring.climate", "autumn.climate") #order treatments
pheno.final$season = factor(pheno.final$season, levels=season, ordered=T)
Treatments     = c("spring", "autumn", "autumn-spring") #order treatments
pheno.final$Treatment = factor(pheno.final$Treatment, levels=Treatments, ordered=T)
Species         = c("Quercus", "Fagus", "Lonicera") #order species
pheno.final$Species   = factor(pheno.final$Species, levels=Species, ordered=T)

#create panel labels
dat_text <- data.frame(
  label       = c("a","b","c","d","e","f"),
  season      = rep(c("spring.climate","autumn.climate"),3),
  Species     = c(rep("Quercus",2), rep("Fagus",2), rep("Lonicera",2)),
  Treatment   = "spring")

#Plot
FigS4 = ggplot(data = pheno.final, mapping = aes(x = Treatment, y = TempDiff, fill=season)) +
  stat_summary(fun.y = mean, 
               geom = "bar", color="black") + 
  stat_summary(fun.data = mean_se,  
               geom = "errorbar", width=0.2) + 
  geom_hline(yintercept=0) +
  geom_text(data = dat_text, mapping = aes(x = -Inf, y = Inf, 
                                           hjust = -0.2, vjust = 1.5,
                                           label = label,
                                           fontface = "bold"))+
  coord_cartesian(ylim=c(-5,5))+
  xlab("Treatment") +
  ylab("Temperature difference (°C)") +
  scale_fill_manual(values = c("green3", "orange")) +
  facet_grid(Species~season) + 
  plotTheme


#Save PDF
#########

pdf(paste(output.dir,"FigS4_TemperatureChange.pdf",sep="/"), width=5, height=7, useDingbats=FALSE)
FigS4
dev.off()



##############################################################################################################################################
#############################################################THE END##########################################################################
##############################################################################################################################################


