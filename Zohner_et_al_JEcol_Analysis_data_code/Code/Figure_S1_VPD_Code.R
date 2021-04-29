##############################################################################################################################################
############################################################# Rscript for: ###################################################################
##############################################################################################################################################
################################################# Zohner et al. 2021 Journal of Ecology ######################################################
#### How changes in spring and autumn phenology translate into growth - experimental evidence of asymmetric effects ##########################
##############################################################################################################################################
##############################################################################################################################################
################################ This script creates Figure S1 ###############################################################################
##############################################################################################################################################



#required packages
require(dplyr)
require(data.table)
require(lubridate)
require(tidyverse)
require(plantecophys)
require(chillR)



##############################################################################################################################################
##############################################################################################################################################



###############
## Functions ##
###############



# Returns a dataframe for gradient backgrounds that can be plotted with geom_rect
#################################################################################

GenerateGradientData <- function(start,
                                 stop,
                                 start_colour,
                                 stop_colour,
                                 y_resolution = 100) {
  
  # define the colour palette
  colour_function <- colorRampPalette(
    c(start_colour, stop_colour),
    alpha = TRUE)
  
  # set up the rect coordinates
  y_range <- seq(start,
                 stop,
                 length.out = y_resolution + 1)
  grad_ymin <- y_range[-length(y_range)]
  grad_ymax <- y_range[c(1:y_resolution + 1)]
  
  # define colours
  grad_colours <- colour_function(y_resolution)
  
  # return data.frame
  data.frame(
    ymin = grad_ymin,
    ymax = grad_ymax,
    xmin = -Inf,
    xmax = Inf,
    grad_colours = grad_colours
  )
}


# Temperate inhibition function from LPJ-GUESS 
##############################################

temp_opt.fun <- function(temp) {
  x1        <- 1
  x2        <- 18
  x3        <- 25
  x4        <- 45
  k1        <- 2.*log((1/0.99)-1.)/(x1-x2)
  k2        <- (x1+x2)/2
  low       <- 1/(1+exp(k1*(k2-temp)))
  k3        <- log(0.99/0.01)/(x4-x3)
  high      <- 1-0.01*exp(k3*(temp-x3))
  tstress   <- low*high 
  if(tstress>=0) {
    tstress <- tstress
  } else {
    tstress <- 0
  }
  return(tstress)
}


# Vapour Pressure Deficit (VPD) function
########################################

# VPD = vapour pressure deficit [kPa]
# T_min & T_max = minimum and maximum daily temperature [C]
# VPD_min --> at low values, latent heat losses are unlikely to exceed available water 
# little effect on stomata
# VPD_max --> at high values, particularly if sustained, photosynthesis and growth are likely to be significantly limited
# complete stomatal closure

VPD.fun <- function(VPD, VPD_min, VPD_max) {
  if(VPD>=VPD_max) {
    y             <- 0
  }
  if(VPD<VPD_max & VPD>VPD_min) {
    y             <- 1-((VPD-VPD_min)/(VPD_max-VPD_min))
  }
  if(VPD<=VPD_min) {
    y             <- 1
  }
  return(y)
}


# Photoperiod function
######################

# photo = photoperiod 
# photo_min = minimum value during the growing season --> limited canopy development
# photo_max = maxmum value during the growing season --> allows canopies to develop unconstrained
photoperiod.fun <- function(photo, photo_min, photo_max) {
  if(photo<=photo_min) {
    photo_resp    <- 0
  }
  if(photo<photo_max & photo>photo_min) {
    photo_resp    <- (photo-photo_min)/(photo_max-photo_min)
  }
  if(photo>=photo_max) {
    photo_resp   <- 1
  }
  return(photo_resp)
}



##############################################################################################################################################
##############################################################################################################################################



################
## Plot theme ##
################



plotTheme = theme(legend.position   = c(0.7, 0.8),
                  legend.background = element_rect(fill=NA, size=0.5, linetype="solid"),
                  legend.text       = element_text(color="black"),
                  panel.grid.major  = element_blank(),
                  panel.grid.minor  = element_blank(),
                  panel.background  = element_rect(colour = NA, size=1, fill="grey97"),
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
data.dir   = "Data"
output.dir = "Output"

#load data
clim.data    = data.frame(fread(paste(data.dir,"S5_climate_hourly.csv", sep="/"), header=T, sep=",")) # Climate 
pheno.data   = data.frame(fread(paste(data.dir,"S2_phenology_data.csv", sep="/"))) %>%                # Phenology
  mutate(Species = gsub(' [A-z ]*', '' , Species)) #delete species epithet



##############################################################################################################################################
##############################################################################################################################################



####################
## Data wrangling ##
####################



clim.data = clim.data %>%
  #get date
  mutate(date = substr(date, start = 1, stop = 8),
         #as date
         date = mdy(date),
         #year
         year = year(date),
         #date to date of year (Jan 1 = day 1)
         DOY  = yday(date) ) %>%
  #get hour of day
  group_by(date) %>%
  mutate(hour = row_number()) %>%
  ungroup()%>%
  #remove 2020 climate
  filter(date < as.Date("2020-01-01"))
  


##############################################################################################################################################
##############################################################################################################################################



#####################################
## Impute temperature and humidity ##
#####################################



#Linear model based on air temperature
######################################

temperature_model=lm(Temp_inside~Temp_outside, data=clim.data)


#Predict inside air temperature for missing cases
#################################################

clim.data.pred = clim.data %>%
  #keep only missing cases
  filter(is.na(Temp_inside)) %>%
  #predict missing cases
  mutate(Temp_inside = predict(temperature_model, ., interval = c("none"),
                              level = 0.95, type = c("response")))


#Add predicitons to dataframe
#############################

clim.data = bind_rows(
  clim.data %>% filter(!is.na(Temp_inside)),
  clim.data.pred)


#---------------------------------------------------------------------------------------------------------


#Linear model based on air temperature and hour of day
######################################################

humidity_model=lm(RH_outside~Temp_outside+as.factor(hour), data=clim.data)


#Predict outside air humidity for missing cases
###############################################

clim.data.pred = clim.data %>%
  #keep only missing cases
  filter(is.na(RH_outside)) %>%
  #predict missing cases
  mutate(RH_outside = predict(humidity_model, ., interval = c("none"),
                              level = 0.95, type = c("response"))) %>%
  #set predicted RH > 100% to 100%
  mutate(RH_outside = ifelse(RH_outside>100,100, RH_outside))


#Add predicitons to dataframe
#############################

clim.data = bind_rows(
  clim.data %>% filter(!is.na(RH_outside)),
  clim.data.pred)
  

#---------------------------------------------------------------------------------------------------------


#Linear model based on air temperature and hour of day
######################################################

humidity_model_inside = lm(RH_inside~RH_outside*Temp_inside*Temp_outside+as.factor(hour), data=clim.data)


#Predict inside air humidity for missing cases
##############################################

clim.data.pred = clim.data %>%
  #keep only missing cases
  filter(is.na(RH_inside)) %>%
  #predict missing cases
  mutate(RH_inside = predict(humidity_model_inside, ., interval = c("none"),
                             level = 0.95, type = c("response"))) %>%
  #set predicted RH > 100% to 100%
  mutate(RH_inside = ifelse(RH_inside>100,100, RH_inside))


#Add predicitons to dataframe
#############################

clim.data = bind_rows(
  clim.data %>% filter(!is.na(RH_inside)),
  clim.data.pred)


#---------------------------------------------------------------------------------------------------------


#to long format
clim.data = bind_cols(
  pivot_longer(clim.data, c(Temp_inside, Temp_outside), names_to = "Treatment", values_to = "Temp") %>%
    dplyr::select(year, date, DOY, hour, Treatment, Temp),
  pivot_longer(clim.data, c(RH_inside, RH_outside), names_to = "Treatment", values_to = "RH")%>%
    dplyr::select(RH) )



##############################################################################################################################################
##############################################################################################################################################



##############
## Analysis ##
##############


#get VPD
clim.data =  clim.data %>%
  mutate(VPD  = RHtoVPD(RH,  Temp,  Pa = 101))

#get mean daytime parameters
clim.data = data.frame(clim.data %>%
                            group_by(year, date, DOY, Treatment) %>%
                            #delete night hours
                            filter(between(hour,
                                           daylength(latitude=48,JDay=DOY[1])$Sunrise,
                                           daylength(latitude=48,JDay=DOY[1])$Sunset))%>%
                            #summarise daytime hours
                            summarise(Temp = mean(Temp),
                                      RH   = mean(RH),
                                      VPD  = mean(VPD)) %>%
                            ungroup() %>%
                            #get day length
                            mutate(DL = daylength(latitude=48,JDay=DOY)$Daylength))


#---------------------------------------------------------------------------------------------------------


#Growing season index parameters
################################

#create vectors
iVPD_fagus_vector = vector()
iTemp_vector      = vector()
iDL_vector        = vector()

for(i in 1:nrow(clim.data)) {
  # apply vapor pressure deficit funtion
  # Reference: White MA, Thornton PE, Running SW et al. (2000) Parameterization and sensitivity analysis of 
  # the BIOME–BGC terrestrial ecosystem model: net primary production controls. Earth Interactions, 4, 1–85.
  iVPD_fagus        = VPD.fun(clim.data[i,]$VPD, 0.6, 3.0)
  iVPD_fagus_vector = c(iVPD_fagus_vector, iVPD_fagus)
  
  # iOpt_temp: response to optimal temperature (Gompertz function)
  iTemp             = temp_opt.fun(clim.data[i,]$Temp)
  iTemp_vector      = c(iTemp_vector, iTemp)
  
  # iPhoto: photoperiod response
  iDL               = photoperiod.fun(clim.data[i,]$DL, min(clim.data$DL), max(clim.data$DL))
  iDL_vector        = c(iDL_vector, iDL)
}

#add to dataframe
clim.data =  clim.data %>%
  mutate(iVPD   = iVPD_fagus_vector,
         iTemp  = iTemp_vector,
         iDL    = iDL_vector,
         iGSI   = iVPD * iTemp * iDL,
         #assign seasons
         season = ifelse(DOY > 243, 'Autumn', ifelse(DOY < 135, 'Spring', 'Summer')),
         season = factor(season, levels=c("Spring", 'Summer',"Autumn"), ordered=T)) %>%
  filter(season!='Summer')



##############################################################################################################################################
##############################################################################################################################################



##########
## Plot ##
##########



# define colours
high_colour <- c('#F21A00')
low_colour <- c('#3B9AB2')

# generate data for a one-hour sunrise gradient
transition_pd <- GenerateGradientData(start = 0.6,
                                   stop = 3,
                                   start_colour = low_colour,
                                   stop_colour = high_colour,
                                   y_resolution = 1000)

# Plot
FigS1 = ggplot(data = clim.data, mapping = aes(x = season, y = VPD, fill=Treatment)) +
  geom_rect(mapping=aes(xmin=-Inf, xmax=Inf, ymin=0, ymax=0.6), color=NA, fill='#3B9AB2',alpha=1) +
  geom_rect(data = transition_pd,
            mapping = aes(xmax = xmax,
                          xmin = xmin,
                          ymax = ymax,
                          ymin = ymin),
            fill = transition_pd$grad_colours,
            inherit.aes = FALSE) +
  geom_hline(yintercept=0)+
  geom_hline(yintercept=0.6, linetype='dashed')+
  coord_cartesian(ylim=c(0.134,2.866))+
  scale_fill_manual(values=c('grey','black'))+
  geom_boxplot(outlier.shape=NA, color='white') + 
  xlab("Season") +
  ylab("Vapor pressure deficit (kPa)") +
  plotTheme


#Save PDF
#########

pdf(paste(output.dir,"FigS1_VPD.pdf",sep="/"), width=3, height=3.5, useDingbats=FALSE)
FigS1
dev.off()



##############################################################################################################################################
#############################################################THE END##########################################################################
##############################################################################################################################################


