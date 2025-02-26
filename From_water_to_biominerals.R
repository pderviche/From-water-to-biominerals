#R Code constructed in February 25, 2025 by Patrick Derviche and Marcelo Soeth

#################################################################################### 
#From water to biominerals: what drives the chemical incorporation in the otolith matrix of juvenile snappers?
#################################################################################### 

#Patrick Derviche, Mario V. Condini, Michael A. Dance, Eduardo S. Costa, Fabian Sá, 
#Felippe A. Daros, Maurício Hostim–Silva, Marcelo Soeth

#################################################
#Summary
#################################################

# A. Finding the best climate window ............................................line 90
#1. Magnesium....................................................................line 130
#1.1 Obtain the p-value
#1.2 inspect the results of the models for linear effect
#1.3 Cross-validation 
#1.4 Save dataset
#1.5 Figures
#1.6 Mg - GLMMs 
#1.7 Mg - Dredge
#1.8 Mg - Model validation


#2. Manganese....................................................................line 640
#2.1 Obtain the p-value
#2.2 inspect the results of the models for linear effect
#2.3 Cross-validation 
#2.4 Save dataset
#2.5 Figures
#2.6 Mn - GLMMs 
#2.7 Mn - Dredge
#2.8  Mn - Model validation

#3. Copper.......................................................................line 1140
#3.1 Obtain the p-value
#3.2 inspect the results of the models for linear effect
#3.3 Cross-validation 
#3.4 Save dataset
#3.5 Figures
#3.6 Cu - GLMMs 
#3.7 Cu - Dredge
#3.8  Cu - Model validation

#4. Zinc.........................................................................line 1630
#4.1 Obtain the p-value
#4.2 inspect the results of the models for linear effect
#4.3 Cross-validation 
#4.4 Save dataset
#4.5 Figures
#4.6 Zn - GLMMs 
#4.7 Zn - Dredge
#4.8  Zn - Model validation

#5. Barium.......................................................................line 2140
#5.1 Obtain the p-value
#5.2 inspect the results of the models for linear effect
#5.3 Cross-validation 
#5.4 Save dataset
#5.5 Figures
#5.6 Ba - GLMMs 
#5.7 Ba - Dredge
#5.8  Ba - Model validation


#B. Figures......................................................................line 2670
#Figure 2 - Laser ablation  sampling transect ...................................line 2680
#Figure 3 - Water samplings......................................................line 2720
#Figure 4 - Water chemistry across the saline gradient...........................line 2820
#Figure 5 - Otolith elemental signatures.........................................line 2980
#Figure 6 - Matrix best temporal window..........................................line 3320
#Figure 7 - Elements relationship................................................line 3340
#Figure 8 - Partition coeficient.................................................line 3360
#Figure S2 - Correlation between variables.......................................line 3380

#C. General results..............................................................line 3430
#Table 1.........................................................................line 3450
#Dog snappers....................................................................line 3480
#Otolith elemental signatures and their relationships with water chemistry.......line 3600
#Kruskal–Wallis’ test and post hoc  Dunn’s test .................................line 3700
#Partition coeficient............................................................line 3770






#################################################
# A. Finding the best climate window
#################################################

###
#Packages
###

library(tidyverse)
library(climwin)
library(lme4)
library(multcomp)
library(lme4)
library(ggplot2)
library(MuMIn)
library(performance)
library(dplyr)
library(FSA)
library (emmeans)
library(ggeffects)
library(car)
library(DHARMa)
library(partR2)
library(scales)
library(ggpubr)

# List of packages to load
packages <- c("tidyverse", "climwin", "lme4", "multcomp", "ggplot2", "MuMIn", 
              "performance", "dplyr", "FSA", "emmeans", "ggeffects", "car", 
              "DHARMa", "partR2", "scales","ggpubr")

# Loop through the list and load each package
for(pkg in packages) {
  library(pkg, character.only = TRUE)}

#Set working directory
setwd("C:/Users/patri/OneDrive/Documentos/UFES/Tese/From water to biominerals/R/github")



###############
#1. Magnesium
###############

#Clean R environment 
rm(list = ls())

####
#Read datasets
####

#Climate variables

water <- read.csv2("water_chemistry.csv") %>%
  filter(tide == "low" &
           site == "Sao Mateus") %>% 
  mutate(across(c(Mn55w, Ca43w, Mn55w, Cu65w, Zn66w, Ba138w,
                  MgCa, MnCa, CuCa, ZnCa, BaCa, Temp, Sal), as.numeric),
         Month = factor(Month, levels = c('JUL', 'AUG', 'SEP', 'OCT', 'NOV', 'DEC', 'JAN', 'FEB', 'MAR', 'APR', 'MAY', 'JUN')),
         Month_num = case_when(Month == "JUL" ~ 7, Month == "AUG" ~ 8, Month == "SEP" ~ 9,
                               Month == "OCT" ~ 10, Month == "NOV" ~ 11, Month == "DEC" ~ 12,
                               Month == "JAN" ~ 1, Month == "FEB" ~ 2, Month == "MAR" ~ 3,
                               Month == "APR" ~ 4, Month == "MAY" ~ 5, Month == "JUN" ~ 6),
         Month_num = as.numeric(Month_num),
         Date = as.Date(paste(Year, Month_num, "15", sep = "/")))
str(water)
head(water)
table(unique(water$Month))

#Fish variables

otolith <- read.csv2("otolith_edge.csv") %>%
  mutate(
    across(c(Mg24, Mn55, Cu65, Zn66, Ba138, TL, Weight, Fulton), as.numeric),
    Month = factor(Month, levels = c('JUL', 'AUG', 'SEP', 'OCT', 'NOV', 'DEC', 'JAN', 'FEB', 'MAR', 'APR', 'MAY', 'JUN')),
    Month_num = case_when(
      Month == "JUL" ~ 7, Month == "AUG" ~ 8, Month == "SEP" ~ 9,
      Month == "OCT" ~ 10, Month == "NOV" ~ 11, Month == "DEC" ~ 12,
      Month == "JAN" ~ 1, Month == "FEB" ~ 2, Month == "MAR" ~ 3,
      Month == "APR" ~ 4, Month == "MAY" ~ 5, Month == "JUN" ~ 6),
    Month_num = as.numeric(Month_num),
    ID = as.factor(ID),
    Date = as.Date(paste(Year, Month_num, "15", sep = "-"))) %>%
  filter(Date > as.Date("2023/01/14")) %>%
  mutate(N = 1:195)


str(otolith)
head(otolith)
length(table(otolith$ID))
table(unique(otolith$Month))

head(otolith)
otolith <- otolith[, c("ID", "Time","Month","Year", "Date","TL", "Weight", "Fulton","Mg24")]

head(water)
water <- water[, c("Sample", "Month","Year", "Date","Temp", "Sal", "MgCa")]

# Mg
#Delete outliers
otolith <- otolith %>%  filter(!(ID == 'ID335' & Time == '70.2'))
otolith <- otolith %>%  filter(!(ID == 'ID204' & Time == '88.3'))
otolith <- otolith %>%  filter(!(ID == 'ID509' & Time == '83.8'))
otolith <- otolith %>%  filter(!(ID == 'ID496' & Time == '101.9'))
otolith <- otolith %>%  filter(!(ID == 'ID383' & Time == '92.8'))
otolith <- otolith %>%  filter(!(ID == 'ID274' & Time == '74.7'))
otolith <- otolith %>%  filter(!(ID == 'ID511' & Time == '83.8'))
otolith <- otolith %>%  filter(!(ID == 'ID487' & Time == '58.9'))
otolith <- otolith %>%  filter(!(ID == 'ID052' &  Time == '52.1'))
otolith <- otolith %>%  filter(!(ID == 'ID052' & Time == '54.3'))

#####
# Step 1. Determine a baseline model structure without weather effects as a null hypothesis
#####

model_baseline <- lmer(Mg24 ~ 1 + TL + Fulton + (1|ID) + (1 | Month), data = otolith)
model_baseline

#####
# Step 2. Create a candidate model set by identifying all competing hypotheses that require testing
#####

variables <- list(MgCa = water$MgCa, #Water chemistry
                  Temp = water$Temp, #Water temperature
                  Sal = water$Sal) #Salinity

#####
# Step 3. Run model set and select best candidate weather signals
#####

results <- slidingwin(xvar = variables,
                      cdate = water$Date, 
                      bdate = otolith$Date, 
                      baseline = model_baseline,
                      range = c(6, 0),
                      refday=c(31, 06),
                      cohort = otolith$Month,
                      type = c("relative"),
                      stat = c("mean"),
                      func = c("lin"), 
                      cmissing = FALSE, 
                      cinterval = "month")

####
#View results sorted by AICc, best supported model is on top
####


#MgCa water
head(results[[1]]$Dataset) 
summary(results[[1]]$BestModel)

#####
#Perform randomization for the models that assume a linear effect 
####

randomized <- randwin(repeats = 1000, ###  long, long time to run  
                      xvar = list(MgCa = water$MgCa),
                      cdate = water$Date, 
                      bdate = otolith$Date, 
                      baseline = model_baseline,
                      range = c(6, 0),
                      refday=c(31, 06),
                      cohort = otolith$Month,
                      type = c("relative"),
                      stat = c("mean"),
                      func = c("lin"), 
                      cmissing = FALSE, 
                      cinterval = "month",
                      window = "sliding")


###
#1.1 Obtain the p-value
###
pvalue <- climwin::pvalue

length(results[[1]]$Dataset)
head(results[[1]]$Dataset)
head(randomized[[1]])

#Mg:Ca water
pvalue(datasetrand = randomized[[1]],dataset = results[[1]]$Dataset,
       metric="AIC", sample.size=17)
#p value 0.499
pvalue(datasetrand = randomized[[1]], dataset = results[[1]]$Dataset,       
       metric="C", sample.size=17)
#C = 0.69

#A higher C value (closer to 1) indicates better discrimination, while a value closer to 0.5 suggests that the model's predictions are no better than random chance




###
#1.2 inspect the results of the models for linear effect
###


#Mg
plotall(datasetrand = randomized[[1]],
        dataset = results[[1]]$Dataset, 
        bestmodel = results[[1]]$BestModel,
        bestmodeldata = results[[1]]$BestModelData,
        title=results$combos[1,])

summary(randomized[[1]]$deltaAICc) 
head(results[[1]]$Dataset) 
#27.17065

Mg_frequency <- ggplot(randomized[[1]], aes(x = deltaAICc)) +
  geom_histogram(binwidth = 1, fill = "#009FC0", color = "#008080",alpha = 0.25) + 
  geom_vline(xintercept = 27.17065, linetype="dashed") +
  labs(x = "ΔAICc of the best model", y = "Frequency", title = (expression(Mg:Ca[otolith]))) +
  theme(plot.title = element_blank()) +
  theme_minimal()
Mg_frequency

# 600 x 350




####
#1.3 Cross-validation 
###

cross_val <- slidingwin(k=10,
                        xvar = list(MgCa = water$MgCa),
                        cdate = water$Date, 
                        bdate = otolith$Date, 
                        baseline = model_baseline,
                        range = c(6, 0),
                        cohort = otolith$Month,
                        type = c("relative"),
                        stat = c("mean"),
                        func = c("lin"), 
                        cmissing = FALSE, 
                        cinterval = "month")

cross_val$combos[1,]
results$combos[1,]

plotall(cross_val[[1]]$Dataset)
plotall(results[[1]]$Dataset)


###
#The final model looks like this: 
###

summary(results[[1]]$BestModel)

#summary(results[[1]]$BestModel)
#Linear mixed model fit by REML ['lmerMod']
#Formula: yvar ~ TL + Fulton + (1 | ID) + (1 | Month) + climate
#Data: modeldat

#REML criterion at convergence: 1746.1

#Scaled residuals: 
#  Min      1Q  Median      3Q     Max 
#-2.3790 -0.4485 -0.0591  0.4729  3.4608 

#Random effects:
#  Groups   Name        Variance Std.Dev.
#ID       (Intercept) 479.3    21.89   
#Month    (Intercept)   0.0     0.00   
#Residual             338.2    18.39   
#Number of obs: 188, groups:  ID, 65; Month, 6

#Fixed effects:
#  Estimate Std. Error t value
#(Intercept)  2.409e+02  7.222e+01   3.336
#TL           1.430e-01  1.105e-01   1.294
#Fulton      -2.057e+01  1.962e+01  -1.048
#climate      6.068e-07  1.337e-06   0.454

#Correlation of Fixed Effects:
#  (Intr) TL     Fulton
#TL      -0.334              
#Fulton  -0.651 -0.125       
#climate -0.892  0.164  0.350
#fit warnings:
#  Some predictor variables are on very different scales: consider rescaling
#optimizer (nloptwrap) convergence code: 0 (OK)
#boundary (singular) fit: see help('isSingular')


#####
#1.4 Save dataset
#####

#stores the values of the climate signal that was best supported by the data 
#row numbers correspond to the original row numbers of the df dataset.

res_MgCa <- results[[1]]$BestModelData
res_Temp <- results[[2]]$BestModelData
res_Sal <- results[[3]]$BestModelData

res_MgCa<-res_MgCa %>%
  rename(Mg_otolith = yvar, Mg_water = climate)
res_MgCa <- res_MgCa[, c("ID", "Mg_water")]
head(res_MgCa)

res_Temp<-res_Temp %>%
  rename(Mg_otolith = yvar, Temp_water = climate)
res_Temp <- res_Temp[, c("ID", "Temp_water")]
head(res_Temp)

res_Sal<-res_Sal %>%
  rename(Mg_otolith = yvar, Sal_water = climate)
res_Sal <- res_Sal[, c("ID", "Sal_water")]
head(res_Sal)

res_MgCa <- res_MgCa %>%
  group_by(ID) %>%
  summarise(Mg_water = mean(Mg_water, na.rm = TRUE)) 

res_Temp <- res_Temp %>%
  group_by(ID) %>%
  summarise(Temp_water = mean(Temp_water, na.rm = TRUE))

res_Sal <- res_Sal %>%
  group_by(ID) %>%
  summarise(Sal_water = mean(Sal_water, na.rm = TRUE))

signal_Mg <- otolith %>%
  full_join(res_MgCa, by = "ID") %>%
  full_join(res_Temp, by = "ID") %>%
  full_join(res_Sal, by = "ID")

head(signal_Mg)
sum(is.na(signal_Mg))
which(is.na(signal_Mg), arr.ind = TRUE)
str(signal_Mg)

#write.table(signal_Mg,"signal_Mg.csv", sep=";", dec=".",row.names = F)

####
#1.5 Figures
####

Mg_1 <-plotdelta(dataset = results[[1]]$Dataset, arrow = F)+
  geom_vline(xintercept = 1, linetype = "dotted")+
  geom_hline(yintercept = 6, linetype = "dotted")+ 
  scale_fill_gradient2(high = "white", mid = "aquamarine2", low = "cyan4", midpoint = 30,
                       name = "∆AICc",
                       breaks = seq(-12, 0, by = 2), labels = seq(-12, 0, by = 2))+
  scale_x_continuous(breaks = 0:6, 
                     labels = c('JUN','MAY', 'APR', 'MAR', 'FEB', 'JAN', 'DEC')) +
  scale_y_continuous(breaks = 0:6, 
                     labels = c('JUN','MAY', 'APR', 'MAR', 'FEB', 'JAN', 'DEC')) +
  theme_minimal()

Mg_1 + theme(plot.title = element_blank(), axis.text.y = element_text(angle = 90))

Mg_2 <-plotwin(dataset = results[[1]]$Dataset, cw = 0.95)+
  theme_minimal()
Mg_2


# 350 x 350

####
# Relationship
####

signal_Mg <- read.csv2("signal_Mg.csv")

signal_Mg <- signal_Mg %>%
  mutate(
    TL = as.numeric(TL),
    Weight = as.numeric(Weight),
    Fulton = as.numeric(Fulton),
    Mg24 = as.numeric(Mg24),
    Mg_water = as.numeric(Mg_water),
    Temp_water = as.numeric(Temp_water),
    Sal_water = as.numeric(Sal_water))

#Get the mean
signal_Mg_mean <- signal_Mg %>%
  group_by(Month) %>%
  summarise(Mg_otolith = mean(Mg24),
            Mg_water = mean(Mg_water))
str(signal_Mg_mean)

Mg_relationship <-  ggplot((signal_Mg_mean), aes(y=Mg_otolith, x=Mg_water/1000000))+
  geom_smooth(method=glm,color = "#00BFFF", fill = "#00BFFF",alpha=0.1)  +
  geom_point(size=3) +
  labs(x= expression(Mg:Ca[water]~(mol~mol^-1)),
       y= expression(Mg:Ca[otolith]~(mu~mol~mol^-1)),
       legend.title = element_blank())+ 
  stat_cor(aes(x = Mg_water, y = Mg_otolith),method = "pearson",label.y = 280)+
  theme_bw()+ylim(230,280)
Mg_relationship
#500 x 400

####
# Partition coeficient
####

library(rcompanion)

signal_Mg$Mg_coef <- signal_Mg$Mg24/signal_Mg$Mg_water
signal_Mg$Month <- factor(signal_Mg$Month, levels = c('JAN','FEB','MAR','APR','MAY','JUN'))

signal_Mg$Mg_coef_scaled <- signal_Mg$Mg_coef * 10000

#letters Mgsed on post hoc Dunn's test
dunn_Mg <- dunnTest(Mg_coef_scaled ~ Month, data = signal_Mg, method = "bonferroni")$res
dunn_Mg
cld_Mg <- cldList(P.adj ~ Comparison, data=dunn_Mg)
cld_Mg
cld_Mg <- as.data.frame(cld_Mg)                    
cld_Mg <- cld_Mg %>%
  rename(Month = Group)
cld_Mg

letter_Mg <- group_by(signal_Mg, Month) %>%
  summarise(mean=mean(Mg_coef_scaled), quant = quantile(Mg_coef_scaled, probs = 0.75)) %>%
  arrange(desc(mean))
letter_Mg <- as.data.frame(letter_Mg)                    
letter_Mg

letter_Mg <- merge(letter_Mg, cld_Mg, by = "Month")
letter_Mg

signal_Mg <- signal_Mg[-102,]

Mg_partition_coeficient <-ggplot(signal_Mg, aes(y=Mg_coef_scaled, x=Month)) +
  geom_boxplot(alpha = 0.5,width=0.3, fill = "#009FC0") +
  theme_bw() +
  xlab("") + ylab(expression(italic(D)[Mg:Ca]~(x~10^-4)))+
  theme(legend.position="none",
        axis.text = element_text(size = 15),
        axis.text.x = element_text(angle = -45),
        axis.title = element_text(size = 15))+
  geom_text(data = letter_Mg, aes(x = Month, y = quant, label = Letter), vjust=-2, hjust=-0.5, size = 4, color='black')
Mg_partition_coeficient


##############
#1.6 Mg - GLMMs 
##############

#Clean R environment 
rm(list = ls())

#Read dataset
signal_Mg <- read.csv2("signal_Mg.csv")

signal_Mg <- signal_Mg %>%
  mutate(
    TL = as.numeric(TL),
    Weight = as.numeric(Weight),
    Fulton = as.numeric(Fulton),
    Mg24 = as.numeric(Mg24),
    Mg_water = as.numeric(Mg_water),
    Temp_water = as.numeric(Temp_water),
    Sal_water = as.numeric(Sal_water))
str(signal_Mg)

head(signal_Mg)

#BoxCox Transformation
box = lm(Mg24 ~ Mg_water + Sal_water + Temp_water + TL + Fulton, data = signal_Mg)

bc<-boxcox(box, plotit = TRUE)
lambda <- bc$x[which.max(bc$y)]

if (lambda == 0) {
  signal_Mg$Mg24t <- log(signal_Mg$Mg24)
} else {
  signal_Mg$Mg24t <- (signal_Mg$Mg24^lambda - 1) / lambda
}


####
#1.7 Mg - Dredge
####

#Scale the variables
signal_Mg$Mg24ts     <-  rescale(signal_Mg$Mg24t)
signal_Mg$TLs        <-  rescale(signal_Mg$TL)
signal_Mg$Fultons   <-  rescale(signal_Mg$Fulton)
signal_Mg$Temp_waters <-  rescale(signal_Mg$Temp_water)
signal_Mg$Sal_waters  <-  rescale(signal_Mg$Sal_water)
signal_Mg$Mg_waters   <-  rescale(signal_Mg$Mg_water)

#Gaussian
full_model_Mg <- lmer(Mg24ts ~ TLs + Fultons + Mg_waters + Sal_waters + Temp_waters + (1 | Month) + (1 | ID), 
                      data = signal_Mg, REML = FALSE,
                      na.action =  na.fail)

dredge_Mg <- dredge(full_model_Mg)
head(dredge_Mg, 10)

#Model selection table 
#(Int)      Flt    Mg_wtr  Sal_wtr Tmp_wtr     TLs df  logLik   AICc delta weight
#1  0.4931                                          4 105.095 -202.0  0.00  0.191
#2  0.5297 -0.1101                                  5 105.632 -200.9  1.04  0.114
#17 0.4514                                 0.07581  5 105.605 -200.9  1.09  0.111
#5  0.4714                 0.05402                  5 105.539 -200.7  1.22  0.104
#21 0.4047                 0.07692         0.10470  6 106.449 -200.4  1.54  0.089
#18 0.4863 -0.1366                         0.09486  6 106.415 -200.4  1.60  0.086
#9  0.4689                         0.03695          5 105.344 -200.4  1.61  0.085
#10 0.5029 -0.1667                 0.06970          6 106.395 -200.3  1.64  0.084
#3  0.4770         0.03087                          5 105.240 -200.1  1.82  0.077

####
## Mg - Selected model
####
str(signal_Mg)

optimal_model_Mg24 <- lmer(Mg24ts ~ TLs + Fultons + (1 | Month) + (1 | ID), 
                           data = signal_Mg, REML = FALSE,
                           na.action =  na.fail)

summary(optimal_model_Mg24)
Anova(optimal_model_Mg24)

#Response: Mg24ts
#Chisq Df Pr(>Chisq)
#TLs     1.5848  1     0.2081
#Fultons 1.6383  1     0.2006

####
#1.8 Mg - Model validation
####
plot(optimal_model_Mg24)

#"DHARMa"
simulationOutput <- simulateResiduals(fittedModel = optimal_model_Mg24, plot = TRUE, use.u = FALSE)
testDispersion(optimal_model_Mg24, type='DHARMa')
#dispersion = 1.0168, p-value = 0.888
# 700 x 375

#r2 values
r.squaredGLMM(optimal_model_Mg24)
#Conditional R2: 0.6095565 
#Marginal R2: 0.02975404 

#Variance inflation factors (VIF)
check_collinearity(optimal_model_Mg24)
#Low Correlation




###############
#2. Manganese
###############

#Clean R environment 
rm(list = ls())

####
#Read datasets
####


#Climate variables

water <- read.csv2("water_chemistry.csv") %>%
  filter(tide == "low" &
           site == "Sao Mateus") %>% 
  mutate(across(c(Mn55w, Ca43w, Mn55w, Cu65w, Zn66w, Ba138w,
                  MnCa, MnCa, CuCa, ZnCa, BaCa, Temp, Sal), as.numeric),
         Month = factor(Month, levels = c('JUL', 'AUG', 'SEP', 'OCT', 'NOV', 'DEC', 'JAN', 'FEB', 'MAR', 'APR', 'MAY', 'JUN')),
         Month_num = case_when(Month == "JUL" ~ 7, Month == "AUG" ~ 8, Month == "SEP" ~ 9,
                               Month == "OCT" ~ 10, Month == "NOV" ~ 11, Month == "DEC" ~ 12,
                               Month == "JAN" ~ 1, Month == "FEB" ~ 2, Month == "MAR" ~ 3,
                               Month == "APR" ~ 4, Month == "MAY" ~ 5, Month == "JUN" ~ 6),
         Month_num = as.numeric(Month_num),
         Date = as.Date(paste(Year, Month_num, "15", sep = "/")))
str(water)
head(water)
table(unique(water$Month))

#Fish variables

otolith <- read.csv2("otolith_edge.csv") %>%
  mutate(
    across(c(Mn55, Mn55, Cu65, Zn66, Ba138, TL, Weight, Fulton), as.numeric),
    Month = factor(Month, levels = c('JUL', 'AUG', 'SEP', 'OCT', 'NOV', 'DEC', 'JAN', 'FEB', 'MAR', 'APR', 'MAY', 'JUN')),
    Month_num = case_when(
      Month == "JUL" ~ 7, Month == "AUG" ~ 8, Month == "SEP" ~ 9,
      Month == "OCT" ~ 10, Month == "NOV" ~ 11, Month == "DEC" ~ 12,
      Month == "JAN" ~ 1, Month == "FEB" ~ 2, Month == "MAR" ~ 3,
      Month == "APR" ~ 4, Month == "MAY" ~ 5, Month == "JUN" ~ 6),
    Month_num = as.numeric(Month_num),
    ID = as.factor(ID),
    Date = as.Date(paste(Year, Month_num, "15", sep = "-"))) %>%
  filter(Date > as.Date("2023/01/14")) %>%
  mutate(N = 1:195)


str(otolith)
head(otolith)
length(table(otolith$ID))
table(unique(otolith$Month))

head(otolith)
otolith <- otolith[, c("ID", "Time","Month","Year", "Date","TL", "Weight", "Fulton","Mn55")]

head(water)
water <- water[, c("Sample", "Month","Year", "Date","Temp", "Sal", "MnCa")]

#Delete outliers
otolith <- otolith[otolith$ID != c('ID165'), ]
otolith <- otolith %>%  filter(!(ID == 'ID382' & Time == '81.5'))
otolith <- otolith %>%  filter(!(ID == 'ID327' & Time == '70.2'))

#####
# Step 1. Determine a baseline model structure without weather effects as a null hypothesis
#####

model_baseline <- lmer(Mn55 ~ 1 + TL + Fulton + (1|ID) + (1 | Month), data = otolith)
model_baseline

#####
# Step 2. Create a candidate model set by identifying all competing hypotheses that require testing
#####

variables <- list(MnCa = water$MnCa, #Water chemistry
                  Temp = water$Temp, #Water temperature
                  Sal = water$Sal) #Salinity

#####
# Step 3. Run model set and select best candidate weather signals
#####

results <- slidingwin(xvar = variables,
                      cdate = water$Date, 
                      bdate = otolith$Date, 
                      baseline = model_baseline,
                      range = c(6, 0),
                      refday=c(31, 06),
                      cohort = otolith$Month,
                      type = c("relative"),
                      stat = c("mean"),
                      func = c("lin"), 
                      cmissing = FALSE, 
                      cinterval = "month")

####
#View results sorted by AICc, best supported model is on top
####


#BaCa water
head(results[[1]]$Dataset) 
summary(results[[1]]$BestModel)

#####
#Perform randomization for the models that assume a linear effect 
####

randomized <- randwin(repeats = 1000, ###  long, long time to run  
                      xvar = list(MnCa = water$MnCa),
                      cdate = water$Date, 
                      bdate = otolith$Date, 
                      baseline = model_baseline,
                      range = c(6, 0),
                      refday=c(31, 06),
                      cohort = otolith$Month,
                      type = c("relative"),
                      stat = c("mean"),
                      func = c("lin"), 
                      cmissing = FALSE, 
                      cinterval = "month",
                      window = "sliding")


###
#2.1 Obtain the p-value
###
pvalue <- climwin::pvalue

length(results[[1]]$Dataset)
head(results[[1]]$Dataset)
head(randomized[[1]])

#Mn:Ca water
pvalue(datasetrand = randomized[[1]],dataset = results[[1]]$Dataset,
       metric="AIC", sample.size=17)
#p value 0.063
pvalue(datasetrand = randomized[[1]], dataset = results[[1]]$Dataset,       
       metric="C", sample.size=17)
#C = 0.63

#A higher C value (closer to 1) indicates better discrimination, while a value closer to 0.5 suggests that the model's predictions are no better than random chance




###
#2.2 inspect the results of the models for linear effect
###


#Mn
plotall(datasetrand = randomized[[1]],
        dataset = results[[1]]$Dataset, 
        bestmodel = results[[1]]$BestModel,
        bestmodeldata = results[[1]]$BestModelData,
        title=results$combos[1,])

summary(randomized[[1]]$deltaAICc) 
head(results[[1]]$Dataset) 
#13.80729

Mn_frequency <- ggplot(randomized[[1]], aes(x = deltaAICc)) +
  geom_histogram(binwidth = 1, fill = "#009FC0", color = "#008080",alpha = 0.25) + 
  geom_vline(xintercept = 13.80729, linetype="dashed") +
  labs(x = "ΔAICc of the best model", y = "Frequency", title = (expression(Mn:Ca[otolith]))) +
  theme(plot.title = element_blank()) +
  theme_minimal()
Mn_frequency

# 600 x 350




####
#2.3 Cross-validation 
###

cross_val <- slidingwin(k=10,
                        xvar = list(MnCa = water$MnCa),
                        cdate = water$Date, 
                        bdate = otolith$Date, 
                        baseline = model_baseline,
                        range = c(6, 0),
                        cohort = otolith$Month,
                        type = c("relative"),
                        stat = c("mean"),
                        func = c("lin"), 
                        cmissing = FALSE, 
                        cinterval = "month")

cross_val$combos[1,]
results$combos[1,]

plotall(cross_val[[1]]$Dataset)
plotall(results[[1]]$Dataset)


###
#The final model looks like this: 
###

#summary(results[[1]]$BestModel)

#Linear mixed model fit by REML ['lmerMod']
#Formula: yvar ~ TL + Fulton + (1 | ID) + (1 | Month) + climate
#Data: modeldat

#REML criterion at convergence: 501.6

#Scaled residuals: 
#  Min      1Q  Median      3Q     Max 
#-1.7557 -0.5413 -0.1245  0.4791  3.0338 

#Random effects:
#  Groups   Name        Variance Std.Dev.
#ID       (Intercept) 0.5283   0.7268  
#Month    (Intercept) 0.0000   0.0000  
#Residual             0.4124   0.6422  
#Number of obs: 193, groups:  ID, 65; Month, 6

#Fixed effects:
#  Estimate Std. Error t value
#(Intercept)  3.3312159  1.7051938   1.954
#TL          -0.0026307  0.0036484  -0.721
#Fulton       0.7381565  0.6156167   1.199
#climate     -0.0003359  0.0001649  -2.037

#Correlation of Fixed Effects:
#  (Intr) TL     Fulton
#TL      -0.291              
#Fulton  -0.459 -0.194       
#climate -0.768  0.024 -0.068
#fit warnings:
#  Some predictor variables are on very different scales: consider rescaling
#optimizer (nloptwrap) convergence code: 0 (OK)
#boundary (singular) fit: see help('isSingular')


#####
#2.4 Save dataset
#####

#stores the values of the climate signal that was best supported by the data 
#row numbers correspond to the original row numbers of the df dataset.

res_MnCa <- results[[1]]$BestModelData
res_Temp <- results[[2]]$BestModelData
res_Sal <- results[[3]]$BestModelData

res_MnCa<-res_MnCa %>%
  rename(Mn_otolith = yvar, Mn_water = climate)
res_MnCa <- res_MnCa[, c("ID", "Mn_water")]
head(res_MnCa)

res_Temp<-res_Temp %>%
  rename(Mn_otolith = yvar, Temp_water = climate)
res_Temp <- res_Temp[, c("ID", "Temp_water")]
head(res_Temp)

res_Sal<-res_Sal %>%
  rename(Mn_otolith = yvar, Sal_water = climate)
res_Sal <- res_Sal[, c("ID", "Sal_water")]
head(res_Sal)

res_MnCa <- res_MnCa %>%
  group_by(ID) %>%
  summarise(Mn_water = mean(Mn_water, na.rm = TRUE)) 

res_Temp <- res_Temp %>%
  group_by(ID) %>%
  summarise(Temp_water = mean(Temp_water, na.rm = TRUE))

res_Sal <- res_Sal %>%
  group_by(ID) %>%
  summarise(Sal_water = mean(Sal_water, na.rm = TRUE))

signal_Mn <- otolith %>%
  full_join(res_MnCa, by = "ID") %>%
  full_join(res_Temp, by = "ID") %>%
  full_join(res_Sal, by = "ID")

head(signal_Mn)
sum(is.na(signal_Mn))
which(is.na(signal_Mn), arr.ind = TRUE)

#write.table(signal_Mn,"signal_Mn.csv", sep=";", dec=".",row.names = F)

####
#2.5 Figures
####

Mn_1 <-plotdelta(dataset = results[[1]]$Dataset, arrow = F)+
  geom_vline(xintercept = 1, linetype = "dotted")+
  geom_hline(yintercept = 6, linetype = "dotted")+ 
  scale_fill_gradient2(high = "white", mid = "aquamarine2", low = "cyan4", midpoint = 18,
                       name = "∆AICc",
                       breaks = seq(-12, 0, by = 2), labels = seq(-12, 0, by = 2))+
  scale_x_continuous(breaks = 0:6, 
                     labels = c('JUN','MAY', 'APR', 'MAR', 'FEB', 'JAN', 'DEC')) +
  scale_y_continuous(breaks = 0:6, 
                     labels = c('JUN','MAY', 'APR', 'MAR', 'FEB', 'JAN', 'DEC')) +
  theme_minimal()

Mn_1 + theme(plot.title = element_blank(), axis.text.y = element_text(angle = 90))

Mn_2 <-plotwin(dataset = results[[1]]$Dataset, cw = 0.95)
Mn_2

# 350 x 350

####
# Relationship
####

signal_Mn <- read.csv2("signal_Mn.csv")

signal_Mn <- signal_Mn %>%
  mutate(
    TL = as.numeric(TL),
    Weight = as.numeric(Weight),
    Fulton = as.numeric(Fulton),
    Mn55 = as.numeric(Mn55),
    Mn_water = as.numeric(Mn_water),
    Temp_water = as.numeric(Temp_water),
    Sal_water = as.numeric(Sal_water))

#Get the mean
signal_Mn_mean <- signal_Mn %>%
  group_by(Month) %>%
  summarise(Mn_otolith = mean(Mn55),
            Mn_water = mean(Mn_water))
str(signal_Mn_mean)

Mn_relationship <-  ggplot((signal_Mn_mean), aes(y=Mn_otolith, x=Mn_water))+
  geom_smooth(method=glm,color = "#00BFFF", fill = "#00BFFF",alpha=0.1)  +
  geom_point(size=3) +
  labs(x= expression(Mn:Ca[water]~(mol~mol^-1)),
       y= expression(Mn:Ca[otolith]~(mu~mol~mol^-1)),
       legend.title = element_blank())+ 
  stat_cor(aes(x = Mn_water, y = Mn_otolith),method = "pearson",label.y = 0.7)+
  theme_bw()
Mn_relationship
#500 x 400

####
# Partition coeficient
####

library(rcompanion)

signal_Mn$Mn_coef <- signal_Mn$Mn55/signal_Mn$Mn_water
signal_Mn$Month <- factor(signal_Mn$Month, levels = c('JAN','FEB','MAR','APR','MAY','JUN'))

signal_Mn$Mn_coef_scaled <- signal_Mn$Mn_coef * 10

#letters Mnsed on post hoc Dunn's test
dunn_Mn <- dunnTest(Mn_coef_scaled ~ Month, data = signal_Mn, method = "bonferroni")$res
dunn_Mn
cld_Mn <- cldList(P.adj ~ Comparison, data=dunn_Mn)
cld_Mn
cld_Mn <- as.data.frame(cld_Mn)                    
cld_Mn <- cld_Mn %>%
  rename(Month = Group)
cld_Mn

letter_Mn <- group_by(signal_Mn, Month) %>%
  summarise(mean=mean(Mn_coef_scaled), quant = quantile(Mn_coef_scaled, probs = 0.75)) %>%
  arrange(desc(mean))
letter_Mn <- as.data.frame(letter_Mn)                    
letter_Mn

letter_Mn <- merge(letter_Mn, cld_Mn, by = "Month")
letter_Mn


Mn_partition_coeficient<-ggplot(signal_Mn, aes(y=Mn_coef_scaled, x=Month)) +
  geom_boxplot(alpha = 0.5,width=0.3, fill = "#009FC0") +
  theme_bw() +
  xlab("") + ylab(expression(italic(D)[Mn:Ca]~(x~10^-1)))+
  theme(legend.position="none",
        axis.text = element_text(size = 15),
        axis.text.x = element_text(angle = -45),
        axis.title = element_text(size = 15))+
  geom_text(data = letter_Mn, aes(x = Month, y = quant, label = Letter), vjust=-2, hjust=-0.5, size = 4, color='black')
Mn_partition_coeficient


##############
#2.6 Mn - GLMMs 
##############

#Clean R environment 
rm(list = ls())

#Read dataset
signal_Mn <- read.csv2("signal_Mn.csv")

signal_Mn <- signal_Mn %>%
  mutate(
    TL = as.numeric(TL),
    Weight = as.numeric(Weight),
    Fulton = as.numeric(Fulton),
    Mn55 = as.numeric(Mn55),
    Mn_water = as.numeric(Mn_water),
    Temp_water = as.numeric(Temp_water),
    Sal_water = as.numeric(Sal_water))
str(signal_Mn)

head(signal_Mn)

#BoxCox Transformation
box = lm(Mn55 ~ Mn_water + Sal_water + Temp_water + TL + Fulton, data = signal_Mn)

bc<-boxcox(box, plotit = TRUE)
lambda <- bc$x[which.max(bc$y)]

if (lambda == 0) {
  signal_Mn$Mn55t <- log(signal_Mn$Mn55)
} else {
  signal_Mn$Mn55t <- (signal_Mn$Mn55^lambda - 1) / lambda
}


####
#2.7 Mn - Dredge
####

#Scale the variables
signal_Mn$Mn55ts     <-  rescale(signal_Mn$Mn55t)
signal_Mn$TLs        <-  rescale(signal_Mn$TL)
signal_Mn$Fultons   <-  rescale(signal_Mn$Fulton)
signal_Mn$Temp_waters <-  rescale(signal_Mn$Temp_water)
signal_Mn$Sal_waters  <-  rescale(signal_Mn$Sal_water)
signal_Mn$Mn_waters   <-  rescale(signal_Mn$Mn_water)

#Gaussian
full_model_Mn <- lmer(Mn55ts ~ TLs + Fultons + Mn_waters + Sal_waters + Temp_waters + (1 | Month) + (1 | ID), 
                      data = signal_Mn, REML = FALSE,
                      na.action =  na.fail)

dredge_Mn <- dredge(full_model_Mn)
head(dredge_Mn,10)

#Model selection table 
#(Int)      Flt    Mn_wtr  Sal_wtr Tmp_wtr     TLs df  logLik   AICc delta weight
#5  0.5778                  -0.11860                   5 46.319 -82.3  0.00  0.188
#3  0.5694         -0.09255                            5 45.951 -81.6  0.74  0.130
#21 0.6373                  -0.13430          -0.0951  6 46.814 -81.2  1.14  0.106
#9  0.5615                           -0.08393          5 45.669 -81.0  1.30  0.098
#1  0.5230                                             4 44.608 -81.0  1.32  0.097
#13 0.5896                  -0.09623 -0.04822          6 46.622 -80.8  1.53  0.088
#7  0.5884         -0.05034 -0.08693                   6 46.606 -80.8  1.56  0.086
#4  0.5275 0.13180 -0.09613                            6 46.447 -80.4  1.88  0.074
#10 0.5151 0.15620                   -0.09569          6 46.343 -80.2  2.08  0.066
#6  0.5675 0.02465          -0.11400                   6 46.334 -80.2  2.10  0.066

####
## Mn - Selected model
####
str(signal_Mn)

optimal_model_Mn55 <- lmer(Mn55ts ~ Mn_waters + Fultons + (1 | Month) + (1 | ID), 
                           data = signal_Mn, REML = FALSE,
                           na.action =  na.fail)

summary(optimal_model_Mn55)
Anova(optimal_model_Mn55)

#Response: Mn55ts
#Chisq Df Pr(>Chisq)  
#Mn_waters 2.9923  1    0.08366 .
#Fultons   0.9992  1    0.31751  

####
#2.8  Mn - Model validation
####
plot(optimal_model_Mn55)

#"DHARMa"
simulationOutput <- simulateResiduals(fittedModel = optimal_model_Mn55, plot = TRUE, use.u = FALSE)
testDispersion(optimal_model_Mn55, type='DHARMa')
#dispersion = 1.0252, p-value = 0.776
# 700 x 375

#r2 values
r.squaredGLMM(optimal_model_Mn55)
#Conditional R2: 0.5241215
#Marginal R2: 0.03800559 

#Variance inflation factors (VIF)
check_collinearity(optimal_model_Mn55)
#Low Correlation





###############
#3. Copper
###############

#Clean R environment 
rm(list = ls())

####
#Read datasets
####

#Climate variables

water <- read.csv2("water_chemistry.csv") %>%
  filter(tide == "low" &
           site == "Sao Mateus") %>% 
  mutate(across(c(Cu65w, Ca43w, Cu65w, Cu65w, Zn66w, Ba138w,
                  CuCa, CuCa, CuCa, ZnCa, BaCa, Temp, Sal), as.numeric),
         Month = factor(Month, levels = c('JUL', 'AUG', 'SEP', 'OCT', 'NOV', 'DEC', 'JAN', 'FEB', 'MAR', 'APR', 'MAY', 'JUN')),
         Month_num = case_when(Month == "JUL" ~ 7, Month == "AUG" ~ 8, Month == "SEP" ~ 9,
                               Month == "OCT" ~ 10, Month == "NOV" ~ 11, Month == "DEC" ~ 12,
                               Month == "JAN" ~ 1, Month == "FEB" ~ 2, Month == "MAR" ~ 3,
                               Month == "APR" ~ 4, Month == "MAY" ~ 5, Month == "JUN" ~ 6),
         Month_num = as.numeric(Month_num),
         Date = as.Date(paste(Year, Month_num, "15", sep = "/")))
str(water)
head(water)
table(unique(water$Month))

#Fish variables

otolith <- read.csv2("otolith_edge.csv") %>%
  mutate(
    across(c(Cu65, Cu65, Cu65, Zn66, Ba138, TL, Weight, Fulton), as.numeric),
    Month = factor(Month, levels = c('JUL', 'AUG', 'SEP', 'OCT', 'NOV', 'DEC', 'JAN', 'FEB', 'MAR', 'APR', 'MAY', 'JUN')),
    Month_num = case_when(
      Month == "JUL" ~ 7, Month == "AUG" ~ 8, Month == "SEP" ~ 9,
      Month == "OCT" ~ 10, Month == "NOV" ~ 11, Month == "DEC" ~ 12,
      Month == "JAN" ~ 1, Month == "FEB" ~ 2, Month == "MAR" ~ 3,
      Month == "APR" ~ 4, Month == "MAY" ~ 5, Month == "JUN" ~ 6),
    Month_num = as.numeric(Month_num),
    ID = as.factor(ID),
    Date = as.Date(paste(Year, Month_num, "15", sep = "-"))) %>%
  filter(Date > as.Date("2023/01/14")) %>%
  mutate(N = 1:195)


str(otolith)
head(otolith)
length(table(otolith$ID))
table(unique(otolith$Month))

head(otolith)
otolith <- otolith[, c("ID", "Time","Month","Year", "Date","TL", "Weight", "Fulton","Cu65")]
otolith <- na.omit(otolith)

head(water)
water <- water[, c("Sample", "Month","Year", "Date","Temp", "Sal", "CuCa")]

#Delete outliers
otolith <- otolith %>%  filter(!(ID == 'ID259' & Time == '40.7'))
otolith <- otolith %>%  filter(!(ID == 'ID502' & Time == '67.9'))

#####
# Step 1. Determine a baseline model structure without weather effects as a null hypothesis
#####

model_baseline <- lmer(Cu65 ~ 1 + TL + Fulton + (1|ID) + (1 | Month), data = otolith)
model_baseline

#####
# Step 2. Create a candidate model set by identifying all competing hypotheses that require testing
#####

variables <- list(CuCa = water$CuCa, #Water chemistry
                  Temp = water$Temp, #Water temperature
                  Sal = water$Sal) #Salinity

#####
# Step 3. Run model set and select best candidate weather signals
#####

results <- slidingwin(xvar = variables,
                      cdate = water$Date, 
                      bdate = otolith$Date, 
                      baseline = model_baseline,
                      range = c(6, 0),
                      refday=c(31, 06),
                      cohort = otolith$Month,
                      type = c("relative"),
                      stat = c("mean"),
                      func = c("lin"), 
                      cmissing = FALSE, 
                      cinterval = "month")

####
#View results sorted by AICc, best supported model is on top
####


#BaCa water
head(results[[1]]$Dataset) 
summary(results[[1]]$BestModel)

#####
#Perform randomization for the models that assume a linear effect 
####

randomized <- randwin(repeats = 1000, ###  long, long time to run  
                      xvar = list(CuCa = water$CuCa),
                      cdate = water$Date, 
                      bdate = otolith$Date, 
                      baseline = model_baseline,
                      range = c(6, 0),
                      refday=c(31, 06),
                      cohort = otolith$Month,
                      type = c("relative"),
                      stat = c("mean"),
                      func = c("lin"), 
                      cmissing = FALSE, 
                      cinterval = "month",
                      window = "sliding")


###
#3.1 Obtain the p-value
###
pvalue <- climwin::pvalue

length(results[[1]]$Dataset)
head(results[[1]]$Dataset)
head(randomized[[1]])

#Cu:Ca water
pvalue(datasetrand = randomized[[1]],dataset = results[[1]]$Dataset,
       metric="AIC", sample.size=17)
#p value 0.042
pvalue(datasetrand = randomized[[1]], dataset = results[[1]]$Dataset,       
       metric="C", sample.size=17)
#C = 0.439782

#A higher C value (closer to 1) indicates better discrimination, while a value closer to 0.5 suggests that the model's predictions are no better than random chance




###
#3.2 inspect the results of the models for linear effect
###


#Cu
plotall(datasetrand = randomized[[1]],
        dataset = results[[1]]$Dataset, 
        bestmodel = results[[1]]$BestModel,
        bestmodeldata = results[[1]]$BestModelData,
        title=results$combos[1,])

summary(randomized[[1]]$deltaAICc) 
head(results[[1]]$Dataset) 
#8.608880

Cu_frequency <- ggplot(randomized[[1]], aes(x = deltaAICc)) +
  geom_histogram(binwidth = 1, fill = "#009FC0", color = "#008080",alpha = 0.25) + 
  geom_vline(xintercept = 8.608880, linetype="dashed") +
  labs(x = "ΔAICc of the best model", y = "Frequency", title = (expression(Cu:Ca[otolith]))) +
  theme(plot.title = element_blank()) +
  theme_minimal()
Cu_frequency

# 600 x 350




####
#3.3 Cross-validation 
###

cross_val <- slidingwin(k=10,
                        xvar = list(CuCa = water$CuCa),
                        cdate = water$Date, 
                        bdate = otolith$Date, 
                        baseline = model_baseline,
                        range = c(6, 0),
                        cohort = otolith$Month,
                        type = c("relative"),
                        stat = c("mean"),
                        func = c("lin"), 
                        cmissing = FALSE, 
                        cinterval = "month")

cross_val$combos[1,]
results$combos[1,]

plotall(cross_val[[1]]$Dataset)
plotall(results[[1]]$Dataset)


###
#The final model looks like this: 
###

summary(results[[1]]$BestModel)

#Linear mixed model fit by REML ['lmerMod']
#Formula: yvar ~ TL + Fulton + (1 | ID) + (1 | Month) + climate
#Data: modeldat

#REML criterion at convergence: 34.4

#Scaled residuals: 
#  Min      1Q  Median      3Q     Max 
#-1.6505 -0.5019 -0.1930  0.2074  4.4669 

#Random effects:
#  Groups   Name        Variance Std.Dev.
#ID       (Intercept) 0.01913  0.1383  
#Month    (Intercept) 0.00000  0.0000  
#Residual             0.04604  0.2146  
#Number of obs: 192, groups:  ID, 65; Month, 6

#Fixed effects:
#  Estimate Std. Error t value
#(Intercept)  0.8523511  0.2819457   3.023
#TL          -0.0008486  0.0008516  -0.996
#Fulton      -0.0496669  0.1400862  -0.355
#climate     -0.0025668  0.0007959  -3.225

#Correlation of Fixed Effects:
#  (Intr) TL     Fulton
#TL      -0.468              
#Fulton  -0.689 -0.195       
#climate -0.466  0.209 -0.041
#optimizer (nloptwrap) convergence code: 0 (OK)
#boundary (singular) fit: see help('isSingular')


#####
#3.4 Save dataset
#####

#stores the values of the climate signal that was best supported by the data 
#row numbers correspond to the original row numbers of the df dataset.

res_CuCa <- results[[1]]$BestModelData
res_Temp <- results[[2]]$BestModelData
res_Sal <- results[[3]]$BestModelData

res_CuCa<-res_CuCa %>%
  rename(Cu_otolith = yvar, Cu_water = climate)
res_CuCa <- res_CuCa[, c("ID", "Cu_water")]
head(res_CuCa)

res_Temp<-res_Temp %>%
  rename(Cu_otolith = yvar, Temp_water = climate)
res_Temp <- res_Temp[, c("ID", "Temp_water")]
head(res_Temp)

res_Sal<-res_Sal %>%
  rename(Cu_otolith = yvar, Sal_water = climate)
res_Sal <- res_Sal[, c("ID", "Sal_water")]
head(res_Sal)

res_CuCa <- res_CuCa %>%
  group_by(ID) %>%
  summarise(Cu_water = mean(Cu_water, na.rm = TRUE)) 

res_Temp <- res_Temp %>%
  group_by(ID) %>%
  summarise(Temp_water = mean(Temp_water, na.rm = TRUE))

res_Sal <- res_Sal %>%
  group_by(ID) %>%
  summarise(Sal_water = mean(Sal_water, na.rm = TRUE))

signal_Cu <- otolith %>%
  full_join(res_CuCa, by = "ID") %>%
  full_join(res_Temp, by = "ID") %>%
  full_join(res_Sal, by = "ID")

head(signal_Cu)
sum(is.na(signal_Cu))
which(is.na(signal_Cu), arr.ind = TRUE)

#write.table(signal_Cu,"signal_Cu.csv", sep=";", dec=".",row.names = F)

####
#3.5 Figures
####

Cu_1 <-plotdelta(dataset = results[[1]]$Dataset, arrow = F)+
  geom_vline(xintercept = 2, linetype = "dotted")+
  geom_hline(yintercept = 6, linetype = "dotted")+ 
  scale_fill_gradient2(high = "white", mid = "aquamarine2", low = "cyan4", midpoint = 12,
                       name = "∆AICc",
                       breaks = seq(-12, 0, by = 2), labels = seq(-12, 0, by = 2))+
  scale_x_continuous(breaks = 0:6, 
                     labels = c('JUN','MAY', 'APR', 'MAR', 'FEB', 'JAN', 'DEC')) +
  scale_y_continuous(breaks = 0:6, 
                     labels = c('JUN','MAY', 'APR', 'MAR', 'FEB', 'JAN', 'DEC')) +
  theme_minimal()

Cu_1 + theme(plot.title = element_blank(), axis.text.y = element_text(angle = 90))

Cu_2 <-plotwin(dataset = results[[1]]$Dataset, cw = 0.95)
Cu_2

# 350 x 350

####
# Relationship
####

signal_Cu <- read.csv2("signal_Cu.csv")

signal_Cu <- signal_Cu %>%
  mutate(
    TL = as.numeric(TL),
    Weight = as.numeric(Weight),
    Fulton = as.numeric(Fulton),
    Cu65 = as.numeric(Cu65),
    Cu_water = as.numeric(Cu_water),
    Temp_water = as.numeric(Temp_water),
    Sal_water = as.numeric(Sal_water))

#Get the mean
signal_Cu_mean <- signal_Cu %>%
  group_by(Month) %>%
  summarise(Cu_otolith = mean(Cu65),
            Cu_water = mean(Cu_water))
str(signal_Cu_mean)

Cu_relationship <-  ggplot((signal_Cu_mean), aes(y=Cu_otolith, x=Cu_water))+
  geom_smooth(method=glm,color = "#00BFFF", fill = "#00BFFF",alpha=0.1)  +
  geom_point(size=3) +
  labs(x= expression(Cu:Ca[water]~(mu~mol~mol^-1)),
       y= expression(Cu:Ca[otolith]~(mu~mol~mol^-1)),
       legend.title = element_blank())+ 
  stat_cor(aes(x = Cu_water, y = Cu_otolith),method = "pearson",label.y = 0.12)+
  theme_bw()
Cu_relationship
#500 x 400

####
# Partition coeficient
####

library(rcompanion)

signal_Cu$Cu_coef <- signal_Cu$Cu65/signal_Cu$Cu_water
signal_Cu$Month <- factor(signal_Cu$Month, levels = c('JAN','FEB','MAR','APR','MAY','JUN'))

#letters Cused on post hoc Dunn's test
dunn_Cu <- dunnTest(Cu_coef ~ Month, data = signal_Cu, method = "bonferroni")$res
dunn_Cu
cld_Cu <- cldList(P.adj ~ Comparison, data=dunn_Cu)
cld_Cu
cld_Cu <- as.data.frame(cld_Cu)                    
cld_Cu <- cld_Cu %>%
  rename(Month = Group)
cld_Cu

letter_Cu <- group_by(signal_Cu, Month) %>%
  summarise(mean=mean(Cu_coef), quant = quantile(Cu_coef, probs = 0.75)) %>%
  arrange(desc(mean))
letter_Cu <- as.data.frame(letter_Cu)                    
letter_Cu

letter_Cu <- merge(letter_Cu, cld_Cu, by = "Month")
letter_Cu


Cu_partition_coeficient<-ggplot(signal_Cu, aes(y=Cu_coef, x=Month)) +
  geom_boxplot(alpha = 0.5,width=0.3, fill = "#009FC0") +
  theme_bw() +
  xlab("") + ylab(expression(italic(D)[Cu:Ca]))+
  theme(legend.position="none",
        axis.text = element_text(size = 15),
        axis.text.x = element_text(angle = -45),
        axis.title = element_text(size = 15))+
  geom_text(data = letter_Cu, aes(x = Month, y = quant, label = Letter), vjust=-2, hjust=-0.5, size = 4, color='black')
Cu_partition_coeficient


##############
#3.6 Cu - GLMMs 
##############

#Clean R environment 
rm(list = ls())

#Read dataset
signal_Cu <- read.csv2("signal_Cu.csv")

signal_Cu <- signal_Cu %>%
  mutate(
    TL = as.numeric(TL),
    Weight = as.numeric(Weight),
    Fulton = as.numeric(Fulton),
    Cu65 = as.numeric(Cu65),
    Cu_water = as.numeric(Cu_water),
    Temp_water = as.numeric(Temp_water),
    Sal_water = as.numeric(Sal_water))
str(signal_Cu)

head(signal_Cu)

#BoxCox Transformation
box = lm(Cu65 ~ Cu_water + Sal_water + Temp_water + TL + Fulton, data = signal_Cu)

bc<-boxcox(box, plotit = TRUE)
lambda <- bc$x[which.max(bc$y)]

if (lambda == 0) {
  signal_Cu$Cu65t <- log(signal_Cu$Cu65)
} else {
  signal_Cu$Cu65t <- (signal_Cu$Cu65^lambda - 1) / lambda
}


####
#3.7 Cu - Dredge
####

#Scale the variables
signal_Cu$Cu65ts     <-  rescale(signal_Cu$Cu65t)
signal_Cu$TLs        <-  rescale(signal_Cu$TL)
signal_Cu$Fultons   <-  rescale(signal_Cu$Fulton)
signal_Cu$Temp_waters <-  rescale(signal_Cu$Temp_water)
signal_Cu$Sal_waters  <-  rescale(signal_Cu$Sal_water)
signal_Cu$Cu_waters   <-  rescale(signal_Cu$Cu_water)

#Gaussian
full_model_Cu <- lmer(Cu65ts ~ TLs + Fultons + Cu_waters + Sal_waters + Temp_waters + (1 | Month) + (1 | ID), 
                      data = signal_Cu, REML = FALSE,
                      na.action =  na.fail)

dredge_Cu <- dredge(full_model_Cu)
head(dredge_Cu,10)

#Model selection table 
#    (Int)   Cu_wtr       Flt  Sal_wtr Tmp_wtr      TLs df logLik   AICc delta weigh
#2  0.5823 -0.10770                                      5 78.561 -146.8  0.00  0.211
#9  0.4682                             0.10960           5 78.390 -146.5  0.34  0.178
#5  0.4868                    0.105100                   5 77.862 -145.4  1.40  0.105
#18 0.6000 -0.11140                            -0.02829  6 78.648 -144.8  1.96  0.079

####
## Cu - Selected model
####
str(signal_Cu)

optimal_model_Cu65 <- lmer(Cu65ts ~ Cu_waters + TLs + (1 | Month) + (1 | ID), 
                           data = signal_Cu, REML = FALSE,
                           na.action =  na.fail)

summary(optimal_model_Cu65)
Anova(optimal_model_Cu65)

#Response: Cu65ts
#Chisq Df Pr(>Chisq)   
#Cu_waters 6.8039  1   0.009096 **
#TLs       0.1730  1   0.677432    

####
#3.8  Cu - Model validation
####
plot(optimal_model_Cu65)

#"DHARMa"
simulationOutput <- simulateResiduals(fittedModel = optimal_model_Cu65, plot = TRUE, use.u = FALSE)
testDispersion(optimal_model_Cu65, type='DHARMa')
#dispersion = 1.01, p-value = 0.904
# 700 x 375

#r2 values
r.squaredGLMM(optimal_model_Cu65)
#Conditional R2: 0.3631317
#Marginal R2: 0.05515462  

#Variance inflation factors (VIF)
check_collinearity(optimal_model_Cu65)
#Low Correlation






###############
#4. Zinc
###############

#Clean R environment 
rm(list = ls())

####
#Read datasets
####


#Climate variables

water <- read.csv2("water_chemistry.csv") %>%
  filter(tide == "low" &
           site == "Sao Mateus") %>% 
  mutate(across(c(Zn66w, Ca43w, Zn66w, Zn66w, Zn66w, Ba138w,
                  ZnCa, ZnCa, ZnCa, ZnCa, BaCa, Temp, Sal), as.numeric),
         Month = factor(Month, levels = c('JUL', 'AUG', 'SEP', 'OCT', 'NOV', 'DEC', 'JAN', 'FEB', 'MAR', 'APR', 'MAY', 'JUN')),
         Month_num = case_when(Month == "JUL" ~ 7, Month == "AUG" ~ 8, Month == "SEP" ~ 9,
                               Month == "OCT" ~ 10, Month == "NOV" ~ 11, Month == "DEC" ~ 12,
                               Month == "JAN" ~ 1, Month == "FEB" ~ 2, Month == "MAR" ~ 3,
                               Month == "APR" ~ 4, Month == "MAY" ~ 5, Month == "JUN" ~ 6),
         Month_num = as.numeric(Month_num),
         Date = as.Date(paste(Year, Month_num, "15", sep = "/")))
str(water)
head(water)
table(unique(water$Month))

#Fish variables

otolith <- read.csv2("otolith_edge.csv") %>%
  mutate(
    across(c(Zn66, Zn66, Zn66, Zn66, Ba138, TL, Weight, Fulton), as.numeric),
    Month = factor(Month, levels = c('JUL', 'AUG', 'SEP', 'OCT', 'NOV', 'DEC', 'JAN', 'FEB', 'MAR', 'APR', 'MAY', 'JUN')),
    Month_num = case_when(
      Month == "JUL" ~ 7, Month == "AUG" ~ 8, Month == "SEP" ~ 9,
      Month == "OCT" ~ 10, Month == "NOV" ~ 11, Month == "DEC" ~ 12,
      Month == "JAN" ~ 1, Month == "FEB" ~ 2, Month == "MAR" ~ 3,
      Month == "APR" ~ 4, Month == "MAY" ~ 5, Month == "JUN" ~ 6),
    Month_num = as.numeric(Month_num),
    ID = as.factor(ID),
    Date = as.Date(paste(Year, Month_num, "15", sep = "-"))) %>%
  filter(Date > as.Date("2023/01/14")) %>%
  mutate(N = 1:195)


str(otolith)
head(otolith)
length(table(otolith$ID))
table(unique(otolith$Month))

head(otolith)
otolith <- otolith[, c("ID", "Time","Month","Year", "Date","TL", "Weight", "Fulton","Zn66")]
otolith <- na.omit(otolith)

head(water)
water <- water[, c("Sample", "Month","Year", "Date","Temp", "Sal", "ZnCa")]

#####
# Step 1. Determine a baseline model structure without weather effects as a null hypothesis
#####

model_baseline <- lmer(Zn66 ~ 1 + TL + Fulton + (1|ID) + (1 | Month), data = otolith)
model_baseline

#####
# Step 2. Create a candidate model set by identifying all competing hypotheses that require testing
#####

variables <- list(ZnCa = water$ZnCa, #Water chemistry
                  Temp = water$Temp, #Water temperature
                  Sal = water$Sal) #Salinity

#####
# Step 3. Run model set and select best candidate weather signals
#####

results <- slidingwin(xvar = variables,
                      cdate = water$Date, 
                      bdate = otolith$Date, 
                      baseline = model_baseline,
                      range = c(6, 0),
                      refday=c(31, 06),
                      cohort = otolith$Month,
                      type = c("relative"),
                      stat = c("mean"),
                      func = c("lin"), 
                      cmissing = FALSE, 
                      cinterval = "month")

####
#View results sorted by AICc, best supported model is on top
####


#BaCa water
head(results[[1]]$Dataset) 
summary(results[[1]]$BestModel)

#####
#Perform randomization for the models that assume a linear effect 
####

randomized <- randwin(repeats = 1000, ###  long, long time to run  
                      xvar = list(ZnCa = water$ZnCa),
                      cdate = water$Date, 
                      bdate = otolith$Date, 
                      baseline = model_baseline,
                      range = c(6, 0),
                      refday=c(31, 06),
                      cohort = otolith$Month,
                      type = c("relative"),
                      stat = c("mean"),
                      func = c("lin"), 
                      cmissing = FALSE, 
                      cinterval = "month",
                      window = "sliding")


###
#4.1 Obtain the p-value
###
pvalue <- climwin::pvalue

length(results[[1]]$Dataset)
head(results[[1]]$Dataset)
head(randomized[[1]])

#Zn:Ca water
pvalue(datasetrand = randomized[[1]],dataset = results[[1]]$Dataset,
       metric="AIC", sample.size=17)
#p value 0.103
pvalue(datasetrand = randomized[[1]], dataset = results[[1]]$Dataset,       
       metric="C", sample.size=17)
#C = 0.50

#A higher C value (closer to 1) indicates better discrimination, while a value closer to 0.5 suggests that the model's predictions are no better than random chance




###
#4.2 inspect the results of the models for linear effect
###


#Zn
plotall(datasetrand = randomized[[1]],
        dataset = results[[1]]$Dataset, 
        bestmodel = results[[1]]$BestModel,
        bestmodeldata = results[[1]]$BestModelData,
        title=results$combos[1,])

summary(randomized[[1]]$deltaAICc) 
head(results[[1]]$Dataset) 
#9.14758

Zn_frequency <- ggplot(randomized[[1]], aes(x = deltaAICc)) +
  geom_histogram(binwidth = 1, fill = "#009FC0", color = "#008080",alpha = 0.25) + 
  geom_vline(xintercept = 9.14758, linetype="dashed") +
  labs(x = "ΔAICc of the best model", y = "Frequency", title = (expression(Zn:Ca[otolith]))) +
  theme(plot.title = element_blank()) +
  theme_minimal()
Zn_frequency

# 600 x 350




####
#4.3 Cross-validation 
###

cross_val <- slidingwin(k=10,
                        xvar = list(ZnCa = water$ZnCa),
                        cdate = water$Date, 
                        bdate = otolith$Date, 
                        baseline = model_baseline,
                        range = c(6, 0),
                        cohort = otolith$Month,
                        type = c("relative"),
                        stat = c("mean"),
                        func = c("lin"), 
                        cmissing = FALSE, 
                        cinterval = "month")

cross_val$combos[1,]
results$combos[1,]

plotall(cross_val[[1]]$Dataset)
plotall(results[[1]]$Dataset)


###
#The final model looks like this: 
###

summary(results[[1]]$BestModel)

#> summary(results[[1]]$BestModel)
#Linear mixed model fit by REML ['lmerMod']
#Formula: yvar ~ TL + Fulton + (1 | ID) + (1 | Month) + climate
#Data: modeldat

#REML criterion at convergence: 792.1

#Scaled residuals: 
#  Min      1Q  Median      3Q     Max 
#-1.5176 -0.3465 -0.1886  0.0152  6.3998 

#Random effects:
#  Groups   Name        Variance Std.Dev.
#ID       (Intercept) 0.5279   0.7265  
#Month    (Intercept) 0.0000   0.0000  
#Residual             2.6799   1.6370  
#Number of obs: 195, groups:  ID, 65; Month, 6

#Fixed effects:
#  Estimate Std. Error t value
#(Intercept)  7.527281   2.462427   3.057
#TL          -0.004348   0.005587  -0.778
#Fulton      -1.431127   0.896224  -1.597
#climate     -0.003089   0.001450  -2.130

#Correlation of Fixed Effects:
#  (Intr) TL     Fulton
#TL      -0.492              
#Fulton  -0.523 -0.184       
#climate -0.763  0.302  0.007
#optimizer (nloptwrap) convergence code: 0 (OK)
#boundary (singular) fit: see help('isSingular')


#####
#4.4 Save dataset
#####

#stores the values of the climate signal that was best supported by the data 
#row numbers correspond to the original row numbers of the df dataset.

res_ZnCa <- results[[1]]$BestModelData
res_Temp <- results[[2]]$BestModelData
res_Sal <- results[[3]]$BestModelData

res_ZnCa<-res_ZnCa %>%
  rename(Zn_otolith = yvar, Zn_water = climate)
res_ZnCa <- res_ZnCa[, c("ID", "Zn_water")]
head(res_ZnCa)

res_Temp<-res_Temp %>%
  rename(Zn_otolith = yvar, Temp_water = climate)
res_Temp <- res_Temp[, c("ID", "Temp_water")]
head(res_Temp)

res_Sal<-res_Sal %>%
  rename(Zn_otolith = yvar, Sal_water = climate)
res_Sal <- res_Sal[, c("ID", "Sal_water")]
head(res_Sal)

res_ZnCa <- res_ZnCa %>%
  group_by(ID) %>%
  summarise(Zn_water = mean(Zn_water, na.rm = TRUE)) 

res_Temp <- res_Temp %>%
  group_by(ID) %>%
  summarise(Temp_water = mean(Temp_water, na.rm = TRUE))

res_Sal <- res_Sal %>%
  group_by(ID) %>%
  summarise(Sal_water = mean(Sal_water, na.rm = TRUE))

signal_Zn <- otolith %>%
  full_join(res_ZnCa, by = "ID") %>%
  full_join(res_Temp, by = "ID") %>%
  full_join(res_Sal, by = "ID")

head(signal_Zn)
sum(is.na(signal_Zn))
which(is.na(signal_Zn), arr.ind = TRUE)

#write.table(signal_Zn,"signal_Zn.csv", sep=";", dec=".",row.names = F)

####
#4.5 Figures
####

Zn_1 <-plotdelta(dataset = results[[1]]$Dataset, arrow = F)+
  geom_vline(xintercept = 1, linetype = "dotted")+
  geom_hline(yintercept = 6, linetype = "dotted")+ 
  scale_fill_gradient2(high = "white", mid = "aquamarine2", low = "cyan4", midpoint = 14,
                       name = "∆AICc",
                       breaks = seq(-12, 0, by = 2), labels = seq(-12, 0, by = 2))+
  scale_x_continuous(breaks = 0:6, 
                     labels = c('JUN','MAY', 'APR', 'MAR', 'FEB', 'JAN', 'DEC')) +
  scale_y_continuous(breaks = 0:6, 
                     labels = c('JUN','MAY', 'APR', 'MAR', 'FEB', 'JAN', 'DEC')) +
  theme_minimal()

Zn_1 + theme(plot.title = element_blank(), axis.text.y = element_text(angle = 90))

Zn_2 <-plotwin(dataset = results[[1]]$Dataset, cw = 0.95)
Zn_2

# 350 x 350

####
# Relationship
####

signal_Zn <- read.csv2("signal_Zn.csv")

signal_Zn <- signal_Zn %>%
  mutate(
    TL = as.numeric(TL),
    Weight = as.numeric(Weight),
    Fulton = as.numeric(Fulton),
    Zn66 = as.numeric(Zn66),
    Zn_water = as.numeric(Zn_water),
    Temp_water = as.numeric(Temp_water),
    Sal_water = as.numeric(Sal_water))

#Get the mean
signal_Zn_mean <- signal_Zn %>%
  group_by(Month) %>%
  summarise(Zn_otolith = mean(Zn66),
            Zn_water = mean(Zn_water))
str(signal_Zn_mean)

Zn_relationship <-  ggplot((signal_Zn_mean), aes(y=Zn_otolith, x=Zn_water))+
  geom_smooth(method=glm,color = "#00BFFF", fill = "#00BFFF",alpha=0.1)  +
  geom_point(size=3) +
  labs(x= expression(Zn:Ca[water]~(mu~mol~mol^-1)),
       y= expression(Zn:Ca[otolith]~(mu~mol~mol^-1)),
       legend.title = element_blank())+ 
  stat_cor(aes(x = Zn_water, y = Zn_otolith),method = "pearson",label.y = 0.1)+
  theme_bw()
Zn_relationship
#500 x 400


####
# Partition coeficient
####

library(rcompanion)

signal_Zn$Zn_coef <- signal_Zn$Zn66/signal_Zn$Zn_water
signal_Zn$Month <- factor(signal_Zn$Month, levels = c('JAN','FEB','MAR','APR','MAY','JUN'))

#letters Znsed on post hoc Dunn's test
dunn_Zn <- dunnTest(Zn_coef ~ Month, data = signal_Zn, method = "bonferroni")$res
dunn_Zn
cld_Zn <- cldList(P.adj ~ Comparison, data=dunn_Zn)
cld_Zn
cld_Zn <- as.data.frame(cld_Zn)                    
cld_Zn <- cld_Zn %>%
  rename(Month = Group)
cld_Zn

letter_Zn <- group_by(signal_Zn, Month) %>%
  summarise(mean=mean(Zn_coef), quant = quantile(Zn_coef, probs = 0.75)) %>%
  arrange(desc(mean))
letter_Zn <- as.data.frame(letter_Zn)                    
letter_Zn

letter_Zn <- merge(letter_Zn, cld_Zn, by = "Month")
letter_Zn


Zn_partition_coeficient<-ggplot(signal_Zn, aes(y=Zn_coef, x=Month)) +
  geom_boxplot(alpha = 0.5,width=0.3, fill = "#009FC0") +
  theme_bw() +
  xlab("") + ylab(expression(italic(D)[Zn:Ca]))+
  theme(legend.position="none",
        axis.text = element_text(size = 15),
        axis.text.x = element_text(angle = -45),
        axis.title = element_text(size = 15))+
  geom_text(data = letter_Zn, aes(x = Month, y = quant, label = Letter), vjust=-2, hjust=-0.5, size = 4, color='black')
Zn_partition_coeficient

##############
#4.6 Zn - GLMMs 
##############

#Clean R environment 
rm(list = ls())

#Read dataset
signal_Zn <- read.csv2("signal_Zn.csv")
signal_Zn <- signal_Zn %>%
  mutate(
    TL = as.numeric(TL),
    Weight = as.numeric(Weight),
    Fulton = as.numeric(Fulton),
    Zn66 = as.numeric(Zn66),
    Zn_water = as.numeric(Zn_water),
    Temp_water = as.numeric(Temp_water),
    Sal_water = as.numeric(Sal_water))
str(signal_Zn)

head(signal_Zn)

#BoxCox Transformation
box = lm(Zn66 ~ Zn_water + Sal_water + Temp_water + TL + Fulton, data = signal_Zn)

bc<-boxcox(box, plotit = TRUE)
lambda <- bc$x[which.max(bc$y)]

if (lambda == 0) {
  signal_Zn$Zn66t <- log(signal_Zn$Zn66)
} else {
  signal_Zn$Zn66t <- (signal_Zn$Zn66^lambda - 1) / lambda
}


####
#4.7 Zn - Dredge
####

#Scale the variables
signal_Zn$Zn66ts     <-  rescale(signal_Zn$Zn66t)
signal_Zn$TLs        <-  rescale(signal_Zn$TL)
signal_Zn$Fultons   <-  rescale(signal_Zn$Fulton)
signal_Zn$Temp_waters <-  rescale(signal_Zn$Temp_water)
signal_Zn$Sal_waters  <-  rescale(signal_Zn$Sal_water)
signal_Zn$Zn_waters   <-  rescale(signal_Zn$Zn_water)

#Gaussian
full_model_Zn <- lmer(Zn66ts ~ TLs + Fultons + Zn_waters + Sal_waters + Temp_waters + (1 | Month) + (1 | ID), 
                      data = signal_Zn, REML = FALSE,
                      na.action =  na.fail)

dredge_Zn <- dredge(full_model_Zn)
dredge_Zn

#Model selection table 
#(Int)       Flt Sal_wtr  Tmp_wtr      TLs   Zn_wtr df  logLik   AICc delta weight
#3  0.5880          0.08055                             5 128.377 -246.4  0.00  0.147
#19 0.6093          0.07012                   -0.03056  6 128.734 -245.0  1.41  0.072
#11 0.6028          0.07551          -0.02324           6 128.481 -244.5  1.92  0.056
#7  0.5836          0.07500  0.01403                    6 128.461 -244.5  1.96  0.055

####
## Zn - Selected model
####
str(signal_Zn)

optimal_model_Zn66 <- lmer(Zn66ts ~ Sal_waters + (1 | Month) + (1 | ID), 
                           data = signal_Zn, REML = FALSE,
                           na.action =  na.fail)

summary(optimal_model_Zn66)
Anova(optimal_model_Zn66)

#Response: Zn66ts
#Chisq Df Pr(>Chisq)
#Sal_waters 4.7312  1    0.02962 *

####
#4.8  Zn - Model validation
####
plot(optimal_model_Zn66)

#"DHARMa"
simulationOutput <- simulateResiduals(fittedModel = optimal_model_Zn66, plot = TRUE, use.u = FALSE)
testDispersion(optimal_model_Zn66, type='DHARMa')
#dispersion = 1.0029, p-value = 0.976
# 700 x 375

#r2 values
r.squaredGLMM(optimal_model_Zn66)
#Conditional R2: 0.4046415
#Marginal R2: 0.05189988 

#Variance inflation factors (VIF)
check_collinearity(optimal_model_Zn66)
#Low Correlation

###
#Figure
###

optimal_model_Zn66_not_trans <- lmer(Zn66 ~ Zn_water + TL + (1 | Month) + (1 | ID), 
                                     data = signal_Zn, REML = FALSE,
                                     na.action =  na.fail)

#Significant variables
Zn <- ggpredict(optimal_model_Zn66_not_trans, terms = c("Zn_water"))

Zn_water<-plot(Zn) +
  labs(title = " ",
       x = expression(Zn:Ca[water]~(mu~mol~mol^-1)), y = expression(Zn:Ca[otolith]~(mu~mol~mol^-1)), color = "")+
  theme_bw()+
  theme(legend.position="none",
        title = element_text(size = 15),
        axis.text = element_text(size = 15),
        axis.title = element_text(size = 15))

Zn_water








###############
#5. Barium
###############

#Clean R environment 
rm(list = ls())

####
#Read datasets
####


#Climate variables

water <- read.csv2("water_chemistry.csv") %>%
  filter(tide == "low" &
           site == "Sao Mateus") %>% 
  mutate(across(c(Mg24w, Ca43w, Mn55w, Cu65w, Zn66w, Ba138w,
                  MgCa, MnCa, CuCa, ZnCa, BaCa, Temp, Sal), as.numeric),
         Month = factor(Month, levels = c('JUL', 'AUG', 'SEP', 'OCT', 'NOV', 'DEC', 'JAN', 'FEB', 'MAR', 'APR', 'MAY', 'JUN')),
         Month_num = case_when(Month == "JUL" ~ 7, Month == "AUG" ~ 8, Month == "SEP" ~ 9,
                               Month == "OCT" ~ 10, Month == "NOV" ~ 11, Month == "DEC" ~ 12,
                               Month == "JAN" ~ 1, Month == "FEB" ~ 2, Month == "MAR" ~ 3,
                               Month == "APR" ~ 4, Month == "MAY" ~ 5, Month == "JUN" ~ 6),
         Month_num = as.numeric(Month_num),
         Date = as.Date(paste(Year, Month_num, "15", sep = "/")))
str(water)
head(water)
table(unique(water$Month))

#Fish variables

otolith <- read.csv2("otolith_edge.csv") %>%
  mutate(
    across(c(Mg24, Mn55, Cu65, Zn66, Ba138, TL, Weight, Fulton), as.numeric),
    Month = factor(Month, levels = c('JUL', 'AUG', 'SEP', 'OCT', 'NOV', 'DEC', 'JAN', 'FEB', 'MAR', 'APR', 'MAY', 'JUN')),
    Month_num = case_when(
      Month == "JUL" ~ 7, Month == "AUG" ~ 8, Month == "SEP" ~ 9,
      Month == "OCT" ~ 10, Month == "NOV" ~ 11, Month == "DEC" ~ 12,
      Month == "JAN" ~ 1, Month == "FEB" ~ 2, Month == "MAR" ~ 3,
      Month == "APR" ~ 4, Month == "MAY" ~ 5, Month == "JUN" ~ 6),
    Month_num = as.numeric(Month_num),
    ID = as.factor(ID),
    Date = as.Date(paste(Year, Month_num, "15", sep = "-"))) %>%
  filter(Date > as.Date("2023/01/14")) %>%
  mutate(N = 1:195)


str(otolith)
head(otolith)
length(table(otolith$ID))
table(unique(otolith$Month))

head(otolith)
otolith <- otolith[, c("ID", "Time","Month","Year", "Date","TL", "Weight", "Fulton","Ba138")]

head(water)
water <- water[, c("Sample", "Month","Year", "Date","Temp", "Sal", "BaCa")]

#####
# Step 1. Determine a baseline model structure without weather effects as a null hypothesis
#####

model_baseline <- lmer(Ba138 ~ 1 + TL + Fulton + (1|ID) + (1 | Month), data = otolith)
model_baseline

#####
# Step 2. Create a candidate model set by identifying all competing hypotheses that require testing
#####

variables <- list(BaCa = water$BaCa, #Water chemistry
                  Temp = water$Temp, #Water temperature
                  Sal = water$Sal) #Salinity

#####
# Step 3. Run model set and select best candidate weather signals
#####

results <- slidingwin(xvar = variables,
                      cdate = water$Date, 
                      bdate = otolith$Date, 
                      baseline = model_baseline,
                      range = c(6, 0),
                      refday=c(31, 06),
                      cohort = otolith$Month,
                      type = c("relative"),
                      stat = c("mean"),
                      func = c("lin"), 
                      cmissing = FALSE, 
                      cinterval = "month")

####
#View results sorted by AICc, best supported model is on top
####


#BaCa water
head(results[[1]]$Dataset) 
summary(results[[1]]$BestModel)

#####
#Perform randomization for the models that assume a linear effect 
####

randomized <- randwin(repeats = 1000, ###  long, long time to run  
                      xvar = list(BaCa = water$BaCa),
                      cdate = water$Date, 
                      bdate = otolith$Date, 
                      baseline = model_baseline,
                      range = c(6, 0),
                      refday=c(31, 06),
                      cohort = otolith$Month,
                      type = c("relative"),
                      stat = c("mean"),
                      func = c("lin"), 
                      cmissing = FALSE, 
                      cinterval = "month",
                      window = "sliding")


###
#5.1 Obtain the p-value
###
pvalue <- climwin::pvalue

length(results[[1]]$Dataset)
head(results[[1]]$Dataset)
head(randomized[[1]])

#Ba:Ca water
pvalue(datasetrand = randomized[[1]],dataset = results[[1]]$Dataset,
       metric="AIC", sample.size=17)
#p value 0.094
pvalue(datasetrand = randomized[[1]], dataset = results[[1]]$Dataset,       
       metric="C", sample.size=17)
#C = 0.69

#A higher C value (closer to 1) indicates better discrimination, while a value closer to 0.5 suggests that the model's predictions are no better than random chance




###
#5.2 inspect the results of the models for linear effect
###


#Ba138w
plotall(datasetrand = randomized[[1]],
        dataset = results[[1]]$Dataset, 
        bestmodel = results[[1]]$BestModel,
        bestmodeldata = results[[1]]$BestModelData,
        title=results$combos[1,])

summary(randomized[[1]]$deltaAICc) 
head(results[[1]]$Dataset) 
#8.642875

Ba_frequency <- ggplot(randomized[[1]], aes(x = deltaAICc)) +
  geom_histogram(binwidth = 1, fill = "#009FC0", color = "#008080",alpha = 0.25) + 
  geom_vline(xintercept = 8.642875, linetype="dashed") +
  labs(x = "ΔAICc of the best model", y = "Frequency", title = (expression(Ba:Ca[otolith]))) +
  theme(plot.title = element_blank()) +
  theme_minimal()
Ba_frequency

# 600 x 350




####
#5.3 Cross-validation 
###

cross_val <- slidingwin(k=10,
                        xvar = list(BaCa = water$BaCa),
                        cdate = water$Date, 
                        bdate = otolith$Date, 
                        baseline = model_baseline,
                        range = c(6, 0),
                        cohort = otolith$Month,
                        type = c("relative"),
                        stat = c("mean"),
                        func = c("lin"), 
                        cmissing = FALSE, 
                        cinterval = "month")

cross_val$combos[1,]
results$combos[1,]

plotall(cross_val[[1]]$Dataset)
plotall(results[[1]]$Dataset)


###
#The final model looks like this: 
###

summary(results[[1]]$BestModel)

#> summary(results[[1]]$BestModel)
#Linear mixed model fit by REML ['lmerMod']
#Formula: yvar ~ TL + Fulton + (1 | ID) + (1 | Month) + climate
#Data: modeldat

#REML criterion at convergence: 1686.1

#Scaled residuals: 
#  Min      1Q  Median      3Q     Max 
#-2.9237 -0.2907 -0.1119  0.1688  5.0501 

#Random effects:
#  Groups   Name        Variance Std.Dev.
#ID       (Intercept) 356.40   18.879  
#Month    (Intercept)  88.22    9.392  
#Residual             168.09   12.965  
#Number of obs: 195, groups:  ID, 65; Month, 6

#Fixed effects:
#  Estimate Std. Error t value
#(Intercept) -97.451958  40.834818  -2.386
#TL            0.102014   0.099088   1.030
#Fulton       36.645511  17.056521   2.148
#climate       0.005354   0.002579   2.076

#Correlation of Fixed Effects:
#  (Intr) TL     Fulton
#TL      -0.462              
#Fulton  -0.642 -0.047       
#climate -0.598  0.094 -0.033
#fit warnings:
#  Some predictor variables are on very different scales: consider rescaling


#####
#5.4 Save dataset
#####

#stores the values of the climate signal that was best supported by the data 
#row numbers correspond to the original row numbers of the df dataset.

res_BaCa <- results[[1]]$BestModelData
res_Temp <- results[[2]]$BestModelData
res_Sal <- results[[3]]$BestModelData

res_BaCa<-res_BaCa %>%
  rename(Ba_otolith = yvar, Ba_water = climate)
res_BaCa <- res_BaCa[, c("ID", "Ba_water")]
head(res_BaCa)

res_Temp<-res_Temp %>%
  rename(Ba_otolith = yvar, Temp_water = climate)
res_Temp <- res_Temp[, c("ID", "Temp_water")]
head(res_Temp)

res_Sal<-res_Sal %>%
  rename(Ba_otolith = yvar, Sal_water = climate)
res_Sal <- res_Sal[, c("ID", "Sal_water")]
head(res_Sal)

res_BaCa <- res_BaCa %>%
  group_by(ID) %>%
  summarise(Ba_water = mean(Ba_water, na.rm = TRUE)) 

res_Temp <- res_Temp %>%
  group_by(ID) %>%
  summarise(Temp_water = mean(Temp_water, na.rm = TRUE))

res_Sal <- res_Sal %>%
  group_by(ID) %>%
  summarise(Sal_water = mean(Sal_water, na.rm = TRUE))

signal_Ba <- otolith %>%
  full_join(res_BaCa, by = "ID") %>%
  full_join(res_Temp, by = "ID") %>%
  full_join(res_Sal, by = "ID")

head(signal_Ba)
sum(is.na(signal_Ba))
which(is.na(signal_Ba), arr.ind = TRUE)

#write.table(signal_Ba,"signal_Ba.csv", sep=";", dec=".",row.names = F)

####
#5.5 Figures
####

Ba_1 <-plotdelta(dataset = results[[1]]$Dataset, arrow = F)+
  geom_vline(xintercept = 1, linetype = "dotted")+
  geom_hline(yintercept = 4, linetype = "dotted")+ 
  scale_fill_gradient2(high = "white", mid = "aquamarine2", low = "cyan4", midpoint = 11,
                       name = "∆AICc",
                       breaks = seq(-12, 0, by = 2), labels = seq(-12, 0, by = 2))+
  scale_x_continuous(breaks = 0:6, 
                     labels = c('JUN','MAY', 'APR', 'MAR', 'FEB', 'JAN', 'DEC')) +
  scale_y_continuous(breaks = 0:6, 
                     labels = c('JUN','MAY', 'APR', 'MAR', 'FEB', 'JAN', 'DEC')) +
  theme_minimal()

Ba_1 + theme(plot.title = element_blank(), axis.text.y = element_text(angle = 90))
#350 x 350

Ba_2 <-plotwin(dataset = results[[1]]$Dataset, cw = 0.95)
Ba_2

####
# Relationship
####

signal_Ba <- read.csv2("signal_Ba.csv")

signal_Ba <- signal_Ba %>%
  mutate(
    TL = as.numeric(TL),
    Weight = as.numeric(Weight),
    Fulton = as.numeric(Fulton),
    Ba138 = as.numeric(Ba138),
    Ba_water = as.numeric(Ba_water),
    Temp_water = as.numeric(Temp_water),
    Sal_water = as.numeric(Sal_water))

#Get the mean
mean_Ba <- signal_Ba %>%
  group_by(Month) %>%
  summarise(Ba_otolith = mean(Ba138),
            Ba_water = mean(Ba_water))
str(mean_Ba)

Ba_relationship <-  ggplot((mean_Ba), aes(y=Ba_otolith, x=Ba_water/1000))+
  geom_smooth(method=glm,color = "#00BFFF", fill = "#00BFFF",alpha=0.1)  +
  geom_point(size=3) +
  labs(x= expression(Ba:Ca[water]~(mmol~mol^-1)),
       y= expression(Ba:Ca[otolith]~(mu~mol~mol^-1)),
       legend.title = element_blank())+ 
  stat_cor(aes(x = Ba_water, y = Ba_otolith),method = "pearson",label.y = 60)+
  theme_bw()
Ba_relationship
#500 x 400

####
# Partition coeficient
####

library(rcompanion)

signal_Ba$Ba_coef <- signal_Ba$Ba138/signal_Ba$Ba_water
signal_Ba$Month <- factor(signal_Ba$Month, levels = c('JAN','FEB','MAR','APR','MAY','JUN'))

#letters based on post hoc Dunn's test
dunn_Ba <- dunnTest(Ba_coef ~ Month, data = signal_Ba, method = "bonferroni")$res
dunn_Ba
cld_Ba <- cldList(P.adj ~ Comparison, data=dunn_Ba)
cld_Ba
cld_Ba <- as.data.frame(cld_Ba)                    
cld_Ba <- cld_Ba %>%
  rename(Month = Group)
cld_Ba

letter_Ba <- group_by(signal_Ba, Month) %>%
  summarise(mean=mean(Ba_coef), quant = quantile(Ba_coef, probs = 0.75)) %>%
  arrange(desc(mean))
letter_Ba <- as.data.frame(letter_Ba)                    
letter_Ba

letter_Ba <- merge(letter_Ba, cld_Ba, by = "Month")
letter_Ba


Ba_partition_coeficient<-ggplot(signal_Ba, aes(y=Ba_coef, x=Month)) +
  geom_boxplot(alpha = 0.5,width=0.3, fill = "#009FC0") +
  theme_bw() +
  xlab("") + ylab(expression(italic(D)[Ba:Ca]))+
  theme(legend.position="none",
        axis.text = element_text(size = 15),
        axis.text.x = element_text(angle = -45),
        axis.title = element_text(size = 15))+
  geom_text(data = letter_Ba, aes(x = Month, y = quant, label = Letter), vjust=-2, hjust=-0.5, size = 4, color='black')
Ba_partition_coeficient


##############
#5.6 Ba - GLMMs 
##############

#Clean R environment 
rm(list = ls())

#Read dataset
signal_Ba <- read.csv2("signal_Ba.csv")
signal_Ba <- signal_Ba %>%
  mutate(
    TL = as.numeric(TL),
    Weight = as.numeric(Weight),
    Fulton = as.numeric(Fulton),
    Ba138 = as.numeric(Ba138),
    Ba_water = as.numeric(Ba_water),
    Temp_water = as.numeric(Temp_water),
    Sal_water = as.numeric(Sal_water))
str(signal_Ba)

#BoxCox Transformation
box = lm(Ba138 ~ Ba_water + Sal_water + Temp_water + TL + Fulton, data = signal_Ba)

bc<-boxcox(box, plotit = TRUE)
lambda <- bc$x[which.max(bc$y)]

if (lambda == 0) {
  signal_Ba$Ba138t <- log(signal_Ba$Ba138)
} else {
  signal_Ba$Ba138t <- (signal_Ba$Ba138^lambda - 1) / lambda
}


####
#5.7 Ba - Dredge
####

#Scale the variables
signal_Ba$Ba138ts     <-  rescale(signal_Ba$Ba138t)
signal_Ba$TLs        <-  rescale(signal_Ba$TL)
signal_Ba$Fultons   <-  rescale(signal_Ba$Fulton)
signal_Ba$Temp_waters <-  rescale(signal_Ba$Temp_water)
signal_Ba$Sal_waters  <-  rescale(signal_Ba$Sal_water)
signal_Ba$Ba_waters   <-  rescale(signal_Ba$Ba_water)

#Gaussian
full_model_Ba <- lmer(Ba138ts ~ TLs + Fultons + Ba_waters + Sal_waters + Temp_waters + (1 | Month) + (1 | ID), 
                      data = signal_Ba, REML = FALSE,
                      na.action =  na.fail)

dredge_Ba <- dredge(full_model_Ba)
head(dredge_Ba,11)

#Model selection table 
#     (Int)    Ba_wtr   Flt     Sal_wtr   Tmp_wtr   TLs     df   logLik    AICc    delta weight
#7  0.5305           0.2180 -0.2854                  6 149.937 -287.4  0.00  0.112
#16 0.7609 -0.160900 0.1964 -0.2716 -0.2215          8 151.881 -287.0  0.44  0.090
#15 0.5580           0.2288 -0.1859 -0.1197          7 150.740 -286.9  0.55  0.085
#23 0.4909           0.1931 -0.2869         0.08821  7 150.512 -286.4  1.00  0.068
#14 0.8697 -0.194000        -0.3171 -0.2311          7 150.462 -286.3  1.10  0.064
#5  0.6110                  -0.3066                  5 148.231 -286.1  1.28  0.059
#21 0.5476                  -0.3054         0.11450  6 149.193 -285.9  1.49  0.053
#11 0.5581           0.2639         -0.2596          6 149.182 -285.9  1.51  0.053
#31 0.5187           0.2043 -0.1888 -0.1179 0.08640  8 151.306 -285.8  1.59  0.050
#8  0.5757 -0.040450 0.2075 -0.3282                  7 150.042 -285.5  1.94  0.042

####
## Ba - Selected model
####
str(signal_Ba)

optimal_model_Ba138 <- lmer(Ba138ts ~ Ba_waters + Sal_waters + Temp_waters + Fultons + (1 | Month) + (1 | ID), 
                            data = signal_Ba, REML = FALSE,
                            na.action =  na.fail)

summary(optimal_model_Ba138)
Anova(optimal_model_Ba138)

#Response: Ba138ts
#            Chisq  Df Pr(>Chisq)  
#Ba_waters   2.3229  1    0.12748  
#Sal_waters  5.8156  1    0.01588 *
#Temp_waters 3.7838  1    0.05175 .
#Fultons     2.9005  1    0.08855 .

####
#5.8  Ba - Model validation
####
plot(optimal_model_Ba138)

#"DHARMa"
simulationOutput <- simulateResiduals(fittedModel = optimal_model_Ba138, plot = TRUE, use.u = FALSE)
testDispersion(optimal_model_Ba138, type='DHARMa')
#dispersion = 1.0052, p-value = 0.984
# 700 x 375

#r2 values
r.squaredGLMM(optimal_model_Ba138)
#Conditional R2: 0.3142307
#Marginal R2: 0.8567717 

#Variance inflation factors (VIF)
check_collinearity(optimal_model_Ba138)
#Low Correlation


###
#Figure
###

optimal_model_Ba138_not_trans <- lmer(Ba138 ~ Ba_water + Sal_water + Temp_water + Fulton + (1 | Month) + (1 | ID), 
                                      data = signal_Ba, REML = FALSE,
                                      na.action =  na.fail)

#Significant variables
Ba <- ggpredict(optimal_model_Ba138_not_trans, terms = c("Sal_water"))

Ba_sal<-plot(Ba) +
  labs(title = " ",
       x = expression(Salinity), y = expression(Ba:Ca[otolith]~(mu~mol~mol^-1)), color = "")+
  theme_bw()+
  theme(legend.position="none",
        title = element_text(size = 15),
        axis.text = element_text(size = 15),
        axis.title = element_text(size = 15))

Ba_sal






















##############
#B. Figures
##############

#Clean R environment 
rm(list = ls())




##############
#Figure 2b. Scatter plot of laser ablation transect, indicating the last three 
#laser spots (light green polygon) used in further analysis.
##############

#Set working directory
setwd("C:/Users/patri/OneDrive/Documentos/UFES/Tese/From water to biominerals/R/github")

#Read data water
data <- read.csv2("ID333.csv")

data <- data[,c('ID','Ba138', 'Time')]
data$Time <- as.numeric(data$Time)

data <- data %>%
  filter(ID == 'ID333')

data$distance <- data$Time*10

data <- mutate_at(data, vars(c('Ba138', 'Time','distance')), as.numeric)
str(data)

#Figure 2
ID333 <-  ggplot(data, aes(x=distance, y=Ba138)) +
  annotate('rect', xmin=890, xmax=995, ymin=0, ymax=90, alpha=0.3, fill='green')+
  annotate('rect', xmin=0, xmax=890, ymin=0, ymax=90, alpha=0.15, fill='gray')+
  geom_point(size=3, pch = 19) +
  theme_bw() +
  theme(legend.position="none",
        axis.text = element_text(size = 15),
        axis.title = element_text(size = 15))+ 
  ylab(expression(Ba:Ca[otolith]~(mu~mol~mol^-1))) +
  scale_x_continuous(expression(Distance~from~core~(mu~m)) )

ID333
#800 x 400




##############
#Figure 3
##############

#Clean R environment 
rm(list = ls())

#Set working directory
setwd("C:/Users/patri/OneDrive/Documentos/UFES/Tese/From water to biominerals/R/github")

#Read data water
water <- read.csv2("water_chemistry.csv")
water <- as.data.frame(water)
str(water)

#Description of 'water' data
#'Temp' is seawater temperature and 'Sal' is salinity
#'Mg24w', 'Ca43w', 'Mn55w', 'Cu65w', 'Zn66w', and 'Ba138w' are in umol L-1
#'MgCa', 'MnCa', 'CuCa', 'ZnCa', and 'BaCa' are the ratio El:Ca43 in umol mol-1

#Filter data
water <- water %>%
  filter(tide == "low") %>%
  filter(site == "Sao Mateus")

#Set as numeric
water <- mutate_at(water, vars(c('Mg24w', 'Ca43w', 'Mn55w', 'Cu65w', 'Zn66w',  'Ba138w',
                                 'MgCa', 'MnCa', 'CuCa', 'ZnCa', 'BaCa',
                                 'Temp', 'Sal')), as.numeric)
water$Month<-as.factor(water$Month)

water <- water %>%
  mutate(Order = case_when(
    Month == "JUL" ~ 1,
    Month == "AUG" ~ 2,
    Month == "SEP" ~ 3,
    Month == "OCT" ~ 4,
    Month == "NOV" ~ 5,
    Month == "DEC" ~ 6,
    Month == "JAN" ~ 7,
    Month == "FEB" ~ 8,
    Month == "MAR" ~ 9,
    Month == "APR" ~ 10,
    Month == "MAY" ~ 11,
    Month == "JUN" ~ 12))

#######
#Fig. 3. Seawater temperature (ºC, green line) and salinity (blue line) during the water samplings from July 2022 to June 2023 
#in the São Mateus estuary, Southwestern Atlantic. Shade and transparent areas indicate the distinct seasons: 
#late dry (from July to September), early wet (from October to December), late wet (from January to March), and early dry (from April to June).
#######

summary(water$Temp)
summary(water$Sal)

ylim.prim <- c(0, 30)
ylim.sec <- c(20, 32)

b <- diff(ylim.sec)/diff(ylim.prim)
a <- ylim.sec[1] - b*ylim.prim[1]

ggplot(data = water, aes(y = Temp, x = Order)) +
  annotate('rect', xmin=3.5, xmax=6.5, ymin=20, ymax=32, alpha=0.15, fill='#999999')+
  annotate('rect', xmin=9.5, xmax=12.5, ymin=20, ymax=32, alpha=0.15, fill='#999999')+
  theme_bw()+
  scale_x_continuous(breaks=seq(1, 12, 1),labels = c('JUL','AUG','SEP','OCT','NOV','DEC','JAN','FEB','MAR','APR','MAY','JUN'))+
  geom_line(color="#008080",size=1)+
  geom_point(aes(),color='#008080',size=5,shape=19)  +
  geom_line(aes(y = a+ Sal*b,),color='#00BFFF',size=1)+
  geom_point(aes(y = a+ Sal*b,),color='#00BFFF',size=5,shape=19)  +
  labs(title = "", x = " ", y = expression(Seawater~temperature~(ºC)))  +
  scale_y_continuous(breaks = seq(0, 32, by=4),sec.axis = sec_axis(~(.-a)/b, name=expression(Salinity)))+
  theme(axis.title.y = element_text(color = "#008080",size = 15),
        axis.title.y.right = element_text(color = "#00BFFF"))+
  theme(legend.position="none",
        title = element_text(size = 15),
        axis.text = element_text(size = 15),
        axis.text.x = element_text(angle = -45),
        axis.title = element_text(size = 15))
#700 x 500

water %>%
  mutate(SeasonGroup = ifelse(grepl("wet", Season, ignore.case = TRUE), "Wet", "Dry")) %>%
  group_by(SeasonGroup) %>%
  summarise(mean_temp = mean(Temp, na.rm = TRUE),
            sd_temp = sd(Temp, na.rm = TRUE))

summary(water$Sal)

#The salinity oscillated from zero during the early wet season and reached almost 30 
#in the dry seasons (Fig. 3). While the water temperature was higher during the wet 
#seasons (28.5 ± 2.5 ºC) compared to the dry seasons (24.3 ± 2.5 ºC; Fig. 3).








############
#Figure 4. Water chemistry across the saline gradient
############

#Set working directory
setwd("C:/Users/patri/OneDrive/Documentos/UFES/Tese/From water to biominerals/R/github")

#Clean R environment 
rm(list = ls())

data <- read.csv2("water_chemistry.csv")
data <- as.data.frame(data)
str(data)

gradient_e_wet <- data %>%
  filter(site == "Sao Mateus") %>%
  filter(tide == "gradient")  %>%
  filter(Season == "Early wet")

gradient_l_wet <- data %>%
  filter(site == "Sao Mateus") %>%
  filter(tide == "gradient")  %>%
  filter(Season == "Late wet")

gradient_e_wet <- gradient_e_wet[,c('Sample', 'Month', 'Year', 'Season', 
                                    'Ca43w',	'Mg24w',	'Mn55w',	'Cu65w',	'Zn66w', 'Ba138w',
                                    'MgCa',  'MnCa',	'CuCa',	'ZnCa',	'BaCa',
                                    'Temp', 'Sal')]

gradient_l_wet<- gradient_l_wet[,c('Sample', 'Month', 'Year', 'Season', 
                                   'Ca43w',	'Mg24w',	'Mn55w',	'Cu65w',	'Zn66w', 'Ba138w',
                                   'MgCa',  'MnCa',	'CuCa',	'ZnCa',	'BaCa',
                                   'Temp', 'Sal')]

#Set as numeric
gradient_e_wet <- mutate_at(gradient_e_wet, vars(c('Ca43w',	'Mg24w',	'Mn55w',	'Cu65w',	'Zn66w', 'Ba138w',
                                                   'MgCa',  'MnCa',	'CuCa',	'ZnCa',	'BaCa',
                                                   'Temp', 'Sal')), as.numeric)

gradient_l_wet <- mutate_at(gradient_l_wet, vars(c('Ca43w',	'Mg24w',	'Mn55w',	'Cu65w',	'Zn66w', 'Ba138w',
                                                   'MgCa',  'MnCa',	'CuCa',	'ZnCa',	'BaCa',
                                                   'Temp', 'Sal')), as.numeric)


###
# Elements
###

# Ca
Ca_e_wet<-ggplot(data = gradient_e_wet, aes(y = Ca43w, x = Sal)) +
  geom_point(color = "#008080", size=3)  +
  labs(title = "Ca",x = "Salinity", y = expression(Ca[water]~ (mu~mol~L^-1)))  +
  theme_bw()+
  geom_smooth(method="lm", formula = y ~ splines::bs(x, 3), se = FALSE,color='#008080')+
  scale_x_continuous(breaks = seq(0, 35, by = 5))
Ca_e_wet

Ca_l_wet<-ggplot(data = gradient_l_wet, aes(y = Ca43w, x = Sal)) +
  geom_point(color = "#00BFFF", size=3)  +
  labs(title = "Ca",x = "Salinity", y = expression(Ca[water]~ (mu~mol~L^-1)))  +
  theme_bw()+
  geom_smooth(method="lm", formula = y ~ splines::bs(x, 3), se = FALSE,color='#00BFFF')+
  scale_x_continuous(breaks = seq(0, 35, by = 5))
Ca_l_wet


###
# Mg:Ca
###

gradient_e_wet$MgCam <- gradient_e_wet$MgCa/1000000 #transform umol mol-1 to mol mol-1
gradient_l_wet$MgCam <- gradient_l_wet$MgCa/1000000 #transform umol mol-1 to mol mol-1

Mg <- ggplot(data = gradient_l_wet, aes(y = MgCam, x = Sal)) +
  geom_point(size=3,color='#00BFFF', alpha = 0.3)  +
  geom_smooth(method="lm", formula = y ~ splines::bs(x, 3), se = FALSE,color='#00BFFF')+
  geom_point(data = gradient_e_wet, aes(y = MgCam, x = Sal), size=3,color='#008080', alpha = 0.3)  +
  geom_smooth(data = gradient_e_wet, aes(y = MgCam, x = Sal), method="lm", formula = y ~ splines::bs(x, 3), se = FALSE,color='#008080')+
  labs(title = " ",x = " ", y = expression(Mg:Ca[water]~ ('mol'~'mol'^-1)))  +
  theme_bw()+
  scale_x_continuous(breaks = seq(0, 35, by = 5))
Mg

###
# Mn:Ca
###

gradient_e_wet$MnCam <- gradient_e_wet$MnCa/1000 #transform umol mol-1 to mmol mol-1
gradient_l_wet$MnCam <- gradient_l_wet$MnCa/1000 #transform umol mol-1 to mmol mol-1

Mn <- ggplot(data = gradient_l_wet, aes(y = MnCam, x = Sal)) +
  geom_point(size=3,color='#00BFFF', alpha = 0.3)  +
  geom_smooth(method="lm", formula = y ~ splines::bs(x, 3), se = FALSE,color='#00BFFF')+
  geom_point(data = gradient_e_wet, aes(y = MnCam, x = Sal), size=3,color='#008080', alpha = 0.3)  +
  geom_smooth(data = gradient_e_wet, aes(y = MnCam, x = Sal), method="lm", formula = y ~ splines::bs(x, 3), se = FALSE,color='#008080')+
  labs(title = " ",x = " ", y = expression(Mn:Ca[water]~ ('mmol'~'mol'^-1)))  +
  theme_bw()+
  scale_x_continuous(breaks = seq(0, 35, by = 5))
Mn

###
# Cu:Ca
###

Cu <- ggplot(data = gradient_l_wet, aes(y = CuCa, x = Sal)) +
  geom_point(size=3,color='#00BFFF', alpha = 0.3)  +
  geom_smooth(method="lm", formula = y ~ splines::bs(x, 3), se = FALSE,color='#00BFFF')+
  geom_point(data = gradient_e_wet, aes(y = CuCa, x = Sal), size=3,color='#008080', alpha = 0.3)  +
  geom_smooth(data = gradient_e_wet, aes(y = CuCa, x = Sal), method="lm", formula = y ~ splines::bs(x, 3), se = FALSE,color='#008080')+
  labs(title = " ",x = " ", y = expression(Cu:Ca[water]~ (mu~mol~mol^-1)))  +
  theme_bw()+
  scale_x_continuous(breaks = seq(0, 35, by = 5))
Cu

###
# Zn:Ca
###

gradient_e_wet$ZnCam <- gradient_e_wet$ZnCa/1000 #transform umol mol-1 to mmol mol-1
gradient_l_wet$ZnCam <- gradient_l_wet$ZnCa/1000 #transform umol mol-1 to mmol mol-1

Zn <- ggplot(data = gradient_l_wet, aes(y = ZnCam, x = Sal)) +
  geom_point(size=3,color='#00BFFF', alpha = 0.3)  +
  geom_smooth(method="lm", formula = y ~ splines::bs(x, 3), se = FALSE,color='#00BFFF')+
  geom_point(data = gradient_e_wet, aes(y = ZnCam, x = Sal), size=3,color='#008080', alpha = 0.3)  +
  geom_smooth(data = gradient_e_wet, aes(y = ZnCam, x = Sal), method="lm", formula = y ~ splines::bs(x, 3), se = FALSE,color='#008080')+
  labs(title = " ",x = " ", y = expression(Zn:Ca[water]~ (mmol~mol^-1)))  +
  theme_bw()+
  scale_x_continuous(breaks = seq(0, 35, by = 5))
Zn

###
# Ba:Ca
####

gradient_e_wet$BaCam <- gradient_e_wet$BaCa/1000 #transform umol mol-1 to mmol mol-1
gradient_l_wet$BaCam <- gradient_l_wet$BaCa/1000 #transform umol mol-1 to mmol mol-1

Ba <- ggplot(data = gradient_l_wet, aes(y = BaCam, x = Sal)) +
  geom_point(size=3,color='#00BFFF', alpha = 0.3)  +
  geom_smooth(method="lm", formula = y ~ splines::bs(x, 3), se = FALSE,color='#00BFFF')+
  geom_point(data = gradient_e_wet, aes(y = BaCam, x = Sal), size=3,color='#008080', alpha = 0.3)  +
  geom_smooth(data = gradient_e_wet, aes(y = BaCam, x = Sal), method="lm", formula = y ~ splines::bs(x, 3), se = FALSE,color='#008080')+
  labs(title = " ",x = "Salinity", y = expression(Ba:Ca[water]~ (mmol~mol^-1)))  +
  theme_bw()+
  scale_x_continuous(breaks = seq(0, 35, by = 5))
Ba


library(ggpubr)

ggarrange(Mg,Mn,Cu,Zn,Ba, 
          labels = c("A","B","C","D","E"),
          ncol = 2, nrow = 3,
          legend = "none")
#700 x 700




##############
#Figure 5 - Otolith elemental signatures
##############

#Set directory
setwd("C:/Users/patri/OneDrive/Documentos/UFES/Tese/From water to biominerals/R/github")

#Clean R environment 
rm(list = ls())

#Read data
data <- read.csv2("otolith_edge.csv")
length(unique(data$ID))
str(data)

#Set as numeric
data <- mutate_at(data, vars(c('Mg24', 'Mn55', 'Cu65','Zn66',  'Ba138', 'TL', 'Fulton')), as.numeric)

data$Month <- factor(data$Month, levels = c('JUL','AUG','SEP','OCT','NOV','DEC','JAN','FEB','MAR','APR','MAY','JUN'))
summary(data$Month)

water <- read.csv2("water_chemistry.csv") %>%
  filter(tide == "low" &
           site == "Sao Mateus") %>% 
  mutate(across(c(Mn55w, Ca43w, Mn55w, Cu65w, Zn66w, Ba138w,
                  MgCa, MnCa, CuCa, ZnCa, BaCa, Temp, Sal), as.numeric),
         Month = factor(Month, levels = c('JUL', 'AUG', 'SEP', 'OCT', 'NOV', 'DEC', 'JAN', 'FEB', 'MAR', 'APR', 'MAY', 'JUN')),
         Month_num = case_when(Month == "JUL" ~ 1, Month == "AUG" ~ 2, Month == "SEP" ~ 3,
                               Month == "OCT" ~ 4, Month == "NOV" ~ 5, Month == "DEC" ~ 6,
                               Month == "JAN" ~ 7, Month == "FEB" ~ 8, Month == "MAR" ~ 9,
                               Month == "APR" ~ 10, Month == "MAY" ~ 11, Month == "JUN" ~ 12),
         Month_num = as.numeric(Month_num),
         Date = as.Date(paste(Year, Month_num, "15", sep = "/")))
head(water)

###
# Mg:Ca
###

#Delete outliers
data_Mg <- data
data_Mg <- data_Mg %>%  filter(!(ID == 'ID335' & Time == '70.2'))
data_Mg <- data_Mg %>%  filter(!(ID == 'ID204' & Time == '88.3'))
data_Mg <- data_Mg %>%  filter(!(ID == 'ID509' & Time == '83.8'))
data_Mg <- data_Mg %>%  filter(!(ID == 'ID496' & Time == '101.9'))
data_Mg <- data_Mg %>%  filter(!(ID == 'ID383' & Time == '92.8'))
data_Mg <- data_Mg %>%  filter(!(ID == 'ID274' & Time == '74.7'))
data_Mg <- data_Mg %>%  filter(!(ID == 'ID511' & Time == '83.8'))
data_Mg <- data_Mg %>%  filter(!(ID == 'ID487' & Time == '58.9'))
data_Mg <- data_Mg %>%  filter(!(ID == 'ID052' &  Time == '52.1'))
data_Mg <- data_Mg %>%  filter(!(ID == 'ID052' & Time == '54.3'))

#transform umol mol-1 to mol mol-1
water$MgCam <- water$MgCa/1000000 

#set scale
scaleMg<- 7

#letters based on post hoc Dunn's test
library(rcompanion)
dunn_Mg <- dunnTest(Mg24 ~ Month, data = data_Mg, method = "bonferroni")$res
dunn_Mg
cld_Mg <- cldList(P.adj ~ Comparison, data=dunn_Mg)
cld_Mg
cld_Mg <- as.data.frame(cld_Mg)                    
cld_Mg <- cld_Mg %>%
  rename(Month = Group)
cld_Mg

letter_Mg <- group_by(data_Mg, Month) %>%
  summarise(mean=mean(Mg24), quant = quantile(Mg24, probs = 0.75)) %>%
  arrange(desc(mean))
letter_Mg <- as.data.frame(letter_Mg)                    
letter_Mg

letter_Mg <- merge(letter_Mg, cld_Mg, by = "Month")
letter_Mg

Mg24 <-  ggplot(data_Mg, aes(x=Month, y=Mg24)) +
  annotate('rect', xmin=3.5, xmax=6.5, ymin=170, ymax=500, alpha=0.15, fill='#999999')+
  annotate('rect', xmin=9.5, xmax=12.5, ymin=170, ymax=500, alpha=0.15, fill='#999999')+
  geom_boxplot(alpha = 1,width=0.4, color = "#008080") +
  geom_smooth(data = water, aes(y = MgCam*scaleMg, x = Month_num),
              alpha=1, size = 0.7,color='#00BFFF',
              method="lm", formula = y ~ splines::bs(x, 3), se = FALSE)+
  theme_bw() +
  theme(legend.position="none",
        axis.text = element_text(size = 15),
        axis.text.x = element_text(angle = -45),
        axis.title = element_text(size = 15),
        axis.title.y = element_text(color = "#008080",size = 15),
        axis.title.y.right = element_text(color = "#00BFFF")) +
  xlab("") + ylab(expression(Mg:Ca[otolith]~(mu~mol~mol^-1)))+
  scale_y_continuous(sec.axis = sec_axis(~./scaleMg, name=expression(Mg:Ca[water]~ ('mol'~'mol'^-1))))+
  geom_text(data = letter_Mg, aes(x = Month, y = quant, label = Letter), vjust=-6,  size = 4.5, color = "#008080")

Mg24

###
# Mn:Ca
###

#Delete outliers
data_Mn <- data
data_Mn <- data[data$ID != c('ID165'), ]
data_Mn <- data_Mn %>%  filter(!(ID == 'ID382' & Time == '81.5'))
data_Mn <- data_Mn %>%  filter(!(ID == 'ID327' & Time == '70.2'))

#transform umol mol-1 to mmol mol-1
water$MnCam <- water$MnCa/1000 

#set scale
scaleMn<- 2

#letters based on post hoc Dunn's test
dunn_Mn <- dunnTest(Mn55 ~ Month, data = data_Mn, method = "bonferroni")$res
dunn_Mn
cld_Mn <- cldList(P.adj ~ Comparison, data=dunn_Mn)
cld_Mn
cld_Mn <- as.data.frame(cld_Mn)                    
cld_Mn <- cld_Mn %>%
  rename(Month = Group)
cld_Mn

letter_Mn <- group_by(data_Mn, Month) %>%
  summarise(mean=mean(Mn55), quant = quantile(Mn55, probs = 0.75)) %>%
  arrange(desc(mean))
letter_Mn <- as.data.frame(letter_Mn)                    
letter_Mn

letter_Mn <- merge(letter_Mn, cld_Mn, by = "Month")
letter_Mn

Mn55 <-  ggplot(data_Mn, aes(x=Month, y=Mn55)) +
  annotate('rect', xmin=3.5, xmax=6.5, ymin=0, ymax=6.5, alpha=0.15, fill='#999999')+
  annotate('rect', xmin=9.5, xmax=12.5, ymin=0, ymax=6.5, alpha=0.15, fill='#999999')+
  geom_boxplot(alpha = 1 ,width=0.4, color = "#008080") +
  geom_smooth(data = water, aes(y =  MnCam/scaleMn, x = Month_num),
              alpha=1, size = 0.7,color='#00BFFF',
              method="lm", formula = y ~ splines::bs(x, 3), se = FALSE)+
  theme_bw() +
  theme(legend.position="none",
        axis.text = element_text(size = 15),
        axis.text.x = element_text(angle = -45),
        axis.title.y = element_text(color = "#008080",size = 15),
        axis.title.y.right = element_text(color = "#00BFFF")) +
  xlab("") + ylab(expression(Mn:Ca[otolith]~(mu~mol~mol^-1)))+
  scale_y_continuous(sec.axis = sec_axis(~.*scaleMn, name=expression(Mn:Ca[water]~ ('mmol'~'mol'^-1))))+
  geom_text(data = letter_Mn, aes(x = Month, y = quant, label = Letter), vjust=-6, size = 4.5,color = "#008080")

Mn55

###
# Cu:Ca
###

#Delete outliers
data_Cu <- data
data_Cu <- data_Cu %>%  filter(!(ID == 'ID083' & Time == '72.4'))
data_Cu <- data_Cu %>%  filter(!(ID == 'ID072' & Time == '81.5'))
data_Cu <- data_Cu %>%  filter(!(ID == 'ID076' & Time == '77'))
data_Cu <- data_Cu %>%  filter(!(ID == 'ID202' & Time == '70.2'))
data_Cu <- data_Cu %>%  filter(!(ID == 'ID202' & Time == '72.4'))
data_Cu <- data_Cu %>%  filter(!(ID == 'ID277' & Time == '72.4'))

#set scale
scaleCu <- 100

#letters based on post hoc Dunn's test
dunn_Cu <- dunnTest(Cu65 ~ Month, data = data_Cu, method = "bonferroni")$res
dunn_Cu
cld_Cu <- cldList(P.adj ~ Comparison, data=dunn_Cu)
cld_Cu
cld_Cu <- as.data.frame(cld_Cu)                    
cld_Cu <- cld_Cu %>%
  rename(Month = Group)
cld_Cu

letter_Cu <- group_by(data_Cu, Month) %>%
  summarise(mean=mean(Cu65), quant = quantile(Cu65, probs = 0.75)) %>%
  arrange(desc(mean))
letter_Cu <- as.data.frame(letter_Cu)                    
letter_Cu

letter_Cu <- merge(letter_Cu, cld_Cu, by = "Month")
letter_Cu

Cu65 <-  ggplot(data_Cu, aes(x=Month, y=Cu65)) +
  annotate('rect', xmin=3.5, xmax=6.5, ymin=0, ymax=2.1, alpha=0.15, fill='#999999')+
  annotate('rect', xmin=9.5, xmax=12.5, ymin=0, ymax=2.1, alpha=0.15, fill='#999999')+
  geom_boxplot(alpha = 1,width=0.4,color = "#008080") +
  geom_smooth(data = water, aes(y =  CuCa/scaleCu, x = Month_num),
              alpha=1, size = 0.7,color='#00BFFF',
              method="lm", formula = y ~ splines::bs(x, 3), se = FALSE)+
  theme_bw() +
  theme(legend.position="none",
        axis.text = element_text(size = 15),
        axis.text.x = element_text(angle = -45),
        axis.title = element_text(size = 15),
        axis.title.y = element_text(color = "#008080",size = 15),
        axis.title.y.right = element_text(color = "#00BFFF")) +
  xlab("") + ylab(expression(Cu:Ca[otolith]~(mu~mol~mol^-1)))+
  scale_y_continuous(sec.axis = sec_axis(~.*scaleCu, name=expression(Cu:Ca[water]~ (mu~mol~mol^-1))))+
  geom_text(data = letter_Cu, aes(x = Month, y = quant, label = Letter), vjust=-6, size = 4.5,color = "#008080")

Cu65

###
# Zn:Ca
###

#Delete outliers
data_Zn <- data
data_Zn <- data_Zn %>%  filter(!(ID == 'ID386' & Time == '72.4'))
data_Zn <- data_Zn %>%  filter(!(ID == 'ID079' & Time == '72.4'))
data_Zn <- data_Zn %>%  filter(!(ID == 'ID509' & Time == '83.8'))
data_Zn <- data_Zn %>%  filter(!(ID == 'ID511' & Time == '83.8'))
data_Zn <- data_Zn %>%  filter(!(ID == 'ID508' & Time == '47.5'))

#transform umol mol-1 to mmol mol-1
water$ZnCam <- water$ZnCa/1000 

#set scale
scaleZn<-5

#letters based on post hoc Dunn's test
dunn_Zn <- dunnTest(Zn66 ~ Month, data = data_Zn, method = "bonferroni")$res
dunn_Zn
cld_Zn <- cldList(P.adj ~ Comparison, data=dunn_Zn)
cld_Zn
cld_Zn <- as.data.frame(cld_Zn)                    
cld_Zn <- cld_Zn %>%
  rename(Month = Group)
cld_Zn

letter_Zn <- group_by(data_Zn, Month) %>%
  summarise(mean=mean(Zn66), quant = quantile(Zn66, probs = 0.75)) %>%
  arrange(desc(mean))
letter_Zn <- as.data.frame(letter_Zn)                    
letter_Zn

letter_Zn <- merge(letter_Zn, cld_Zn, by = "Month")
letter_Zn

Zn66 <-  ggplot(data_Zn, aes(x=Month, y=Zn66)) +
  annotate('rect', xmin=3.5, xmax=6.5, ymin=0, ymax=7.5, alpha=0.15, fill='#999999')+
  annotate('rect', xmin=9.5, xmax=12.5, ymin=0, ymax=7.5, alpha=0.15, fill='#999999')+
  geom_boxplot(alpha = 1,width=0.4,color = "#008080") +
  geom_smooth(data = water, aes(y = ZnCam*scaleZn, x = Month_num),
              alpha=1, size = 0.7,color='#00BFFF',
              method="lm", formula = y ~ splines::bs(x, 3), se = FALSE)+
  theme_bw() +
  theme(legend.position="none",
        axis.text = element_text(size = 15),
        axis.text.x = element_text(angle = -45),
        axis.title = element_text(size = 15),
        axis.title.y = element_text(color = "#008080",size = 15),
        axis.title.y.right = element_text(color = "#00BFFF")) +
  xlab("") + ylab(expression(Zn:Ca[otolith]~(mu~mol~mol^-1)))+
  scale_y_continuous(sec.axis = sec_axis(~./scaleZn, name=expression(Zn:Ca[water]~ (mmol~mol^-1))))+
  geom_text(data = letter_Zn, aes(x = Month, y = quant, label = Letter), vjust=-6, size = 4.5,color = "#008080")

Zn66

###
# Ba:Ca
###

data_Ba<-data

#transform umol mol-1 to mmol mol-1
water$BaCam <- water$BaCa/1000 

#set scale
scaleBa<-12

#letters based on post hoc Dunn's test
dunn_Ba <- dunnTest(Ba138 ~ Month, data = data_Ba, method = "bonferroni")$res
dunn_Ba
cld_Ba <- cldList(P.adj ~ Comparison, data=dunn_Ba)
cld_Ba
cld_Ba <- as.data.frame(cld_Ba)                    
cld_Ba <- cld_Ba %>%
  rename(Month = Group)
cld_Ba

letter_Ba <- group_by(data_Ba, Month) %>%
  summarise(mean=mean(Ba138), quant = quantile(Ba138, probs = 0.75)) %>%
  arrange(desc(mean))
letter_Ba <- as.data.frame(letter_Ba)                    
letter_Ba

letter_Ba <- merge(letter_Ba, cld_Ba, by = "Month")
letter_Ba

Ba138 <-  ggplot(data_Ba, aes(x=Month, y=Ba138)) +
  annotate('rect', xmin=3.5, xmax=6.5, ymin=0, ymax=165, alpha=0.15, fill='#999999')+
  annotate('rect', xmin=9.5, xmax=12.5, ymin=0, ymax=165, alpha=0.15, fill='#999999')+
  geom_boxplot(alpha = 1,width=0.4,color = "#008080") +
  geom_smooth(data = water, aes(y = BaCam*scaleBa, x = Month_num),
              alpha=1, size = 0.7,color='#00BFFF',
              method="lm", formula = y ~ splines::bs(x, 3), se = FALSE)+
  theme_bw() +
  theme(legend.position="none",
        axis.text = element_text(size = 15),
        axis.text.x = element_text(angle = -45),
        axis.title = element_text(size = 15),
        axis.title.y = element_text(color = "#008080",size = 15),
        axis.title.y.right = element_text(color = "#00BFFF")) +
  xlab("") + ylab(expression(Ba:Ca[otolith]~(mu~mol~mol^-1)))+
  scale_y_continuous(sec.axis = sec_axis(~./scaleBa, name=expression(Ba:Ca[water]~ (mmol~mol^-1))))+
  geom_text(data = letter_Ba, aes(x = Month, y = quant, label = Letter), vjust=-6, size = 4.5,color = "#008080")

Ba138


#Fig. 5. Otolith elemental signatures (El:Caotolith) of the juvenile dog snappers (green boxplots) 
#and water chemistry (El:Cawater, blue lines) monthly sampled at the São Mateus estuary – Southwestern Atlantic. 
#El:Cawater is presented with a one–month correction, and estimated curves (blue lines) were done by applying 
#linear models using a third–degree basis spline. Shade and transparent areas indicate the distinct seasons: 
#late dry (from July to September), early wet (from October to December), late wet (from January to March), 
#and early dry (from April to June). Lettering indicates post hoc Dunn’s test results in which months are significantly 
#different (p < 0.05), considering El:Caotolith as the response variable, and the sampled month as the predictor. 

ggarrange(Mg24,Mn55,Cu65,Zn66,Ba138, 
          labels = c("A","B","C","D","E"),
          ncol = 2, nrow =3,
          legend = "none")
#1000 x 1000











##############
#Figure 6 - Matrix best temporal window
##############

ggarrange(Mg_relationship, Mn_relationship, Cu_relationship,
          Zn_relationship, Ba_relationship, 
          labels = c("A","B","C","D","E"),
          ncol = 2, nrow = 3,
          legend = "none")
#800 x 800










##############
#Figure 7 - Elements relationship
##############

ggarrange(Mg_relationship, Mn_relationship, Cu_relationship,
          Zn_relationship, Ba_relationship, 
          labels = c("A","B","C","D","E"),
          ncol = 2, nrow = 3,
          legend = "none")
#800 x 800










##############
#Figure 8 - Partition coeficient
##############

ggarrange(Mg_partition_coeficient, Mn_partition_coeficient, Cu_partition_coeficient,
          Zn_partition_coeficient, Ba_partition_coeficient, 
          labels = c("A","B","C","D","E"),
          ncol = 2, nrow = 3,
          legend = "none")
#800 x 800










##############
# Figure S2 - Correlation between variables
##############

#Clean R environment 
rm(list = ls())

#Set working directory
setwd("C:/Users/patri/OneDrive/Documentos/UFES/Tese/From water to biominerals/R/github")

#Read data
otolith <- read.csv2("otolith_edge.csv")
water <- read.csv2("water_chemistry.csv")
water <- water[water$tide != c('gradient'), ]

data <- otolith %>%
  full_join(water, by = "Month")

data <- mutate_at(data, vars(c("Temp", "Sal", "TL", "Fulton")), as.numeric)

str(data)
matrix <- data[, c("Temp", "Sal", "TL", "Fulton")]
str(matrix)

library(corrplot)
par(mfrow = c(1, 1))
corrplot(cor(matrix[], method ='pearson'), type = "lower",
         method = c("number"), tl.col = "black", number.digits = 2)

library(GGally)
require(datasets)
data("swiss")

ggpairs(matrix)

#Fig. S2. Pairplots produced for identifying highly correlated explanatory variables during the exploratory data analysis. 

ggpairs(matrix, columns = 1:4, ggplot2::aes(alpha = 0.3),
        lower = list(continuous = "smooth", alpha = 0.3)) +   theme_bw()











############
#C. General results
############

#Clean R environment 
rm(list = ls())

data <- read.csv2("otolith_edge.csv")
data <- mutate_at(data, vars(c("Elapsed.Time","Time","TL","Weight","Fulton",
                               "Mg24","Mn55","Cu65","Zn66","Ba138")), as.numeric)
data$Month <- factor(data$Month, levels = c('JUL','AUG','SEP','OCT','NOV','DEC','JAN','FEB','MAR','APR','MAY','JUN'))
str(data)
options(pillar.sigfig=4)







##########
#Table 1. Sample size (n), mean ± standard deviation, range (minimum and maximum) of total length (TL, in mm), 
#and Fulton’s index (FI), for juvenile dog snapper individuals monthly collected at the São Mateus estuary –
#Southwestern Atlantic.
##########

# TL Mean SD
data %>%
  group_by(Month) %>%
  summarise(
    mean_TL = mean(TL, na.rm = TRUE),
    sd_TL = sd(TL, na.rm = TRUE),
    max_TL = max(TL, na.rm = TRUE),
    min_TL = min(TL, na.rm = TRUE))

# Fulton Mean SD
data %>%
  group_by(Month) %>%
  summarise(
    mean = mean(Fulton, na.rm = TRUE),
    sd = sd(Fulton, na.rm = TRUE),
    max = max(Fulton, na.rm = TRUE),
    min = min(Fulton, na.rm = TRUE))







############
#Dog snappers
############

####
#TL
####

#Assumptions
model <- aov(TL ~ Month, data = data)
leveneTest(TL ~ Month, data = data) #p-value = 0.004123 **, we cannot assume the homogeneity of variance
shapiro.test(residuals(model)) #p-value = 0.0003134, we cannot assume the normality

#Non-parametric test
kruskal.test(TL ~ Month, data = data) #Kruskal-Wallis chi-squared = 89.163, df = 11, p-value = 2.431e-14
dunnTest(TL ~ Month, data = data) #post hoc

summary(data$TL)
sd(data$TL)
#The TL of the juvenile dog snappers was 173.7 ± 32.5 mm (mean ± standard deviation), 
#ranging from 72.0 to 260.0 mm (Table 1), showing significant differences between the
#sampled months (Kruskal–Wallis, n = 120, χ2 = 89.16, p < 0.001). 


#letters based on post hoc Dunn's test
dunn_TL <- dunnTest(TL ~ Month, data = data, method = "bonferroni")$res
dunn_TL
cld_TL <- cldList(P.adj ~ Comparison, data=dunn_TL)
cld_TL
cld_TL <- as.data.frame(cld_TL)                    
cld_TL <- cld_TL %>%
  rename(Month = Group)
cld_TL

letter_TL <- group_by(data, Month) %>%
  summarise(mean=mean(TL), quant = quantile(TL, probs = 0.75)) %>%
  arrange(desc(mean))
letter_TL <- as.data.frame(letter_TL)                    
letter_TL

letter_TL <- merge(letter_TL, cld_TL, by = "Month")
letter_TL

TL <-  ggplot(data, aes(x=Month, y=TL)) +
  geom_boxplot(width=0.4,fill = "#008080", alpha=0.3) +
  theme_bw() +
  theme(legend.position="none",
        axis.text = element_text(size = 15),
        axis.text.x = element_text(angle = -45),
        axis.title = element_text(size = 15)) +
  geom_text(data = letter_TL, aes(x = Month, y = quant, label = Letter), vjust=-4, size = 4.5,color = "#008080")

TL


####
#Fulton
####

#Assumptions
model <- aov(Fulton ~ Month, data = data)
leveneTest(Fulton ~ Month, data = data) #p-value = 0.007608 **, we cannot assume the homogeneity of variance
shapiro.test(residuals(model)) #p-value = 1.948e-08, we cannot assume the normality

#Non-parametric test
kruskal.test(Fulton ~ Month, data = data) #Kruskal-Wallis chi-squared = 117.66, df = 11, p-value < 2.2e-16
dunnTest(Fulton ~ Month, data = data) #post hoc

summary(data$Fulton)
sd(data$Fulton)

#The Fulton’s index was 1.62 ± 0.15, ranging from 1.29 to 2.31 (Table 1), also 
#displaying significant differences between the sampled months (Kruskal–Wallis, n = 120, χ2 = 117.66, p < 0.001). 

#letters based on post hoc Dunn's test
dunn_Fulton <- dunnTest(Fulton ~ Month, data = data, method = "bonferroni")$res
dunn_Fulton
cld_Fulton <- cldList(P.adj ~ Comparison, data=dunn_Fulton)
cld_Fulton
cld_Fulton <- as.data.frame(cld_Fulton)                    
cld_Fulton <- cld_Fulton %>%
  rename(Month = Group)
cld_Fulton

letter_Fulton <- group_by(data, Month) %>%
  summarise(mean=mean(Fulton), quant = quantile(Fulton, probs = 0.75)) %>%
  arrange(desc(mean))
letter_Fulton <- as.data.frame(letter_Fulton)                    
letter_Fulton

letter_Fulton <- merge(letter_Fulton, cld_Fulton, by = "Month")
letter_Fulton

Fulton <-  ggplot(data, aes(x=Month, y=Fulton)) +
  geom_boxplot(width=0.4,fill = "#00BFFF", alpha=0.3) +
  theme_bw() +
  theme(legend.position="none",
        axis.text = element_text(size = 15),
        axis.text.x = element_text(angle = -45),
        axis.title = element_text(size = 15)) +
  geom_text(data = letter_Fulton, aes(x = Month, y = quant, label = Letter), vjust=-4, size = 4.5,color = "#008080") + 
  ylab(expression(Fulton~(italic(K))))

Fulton

ggarrange(TL, Fulton, 
          labels = c("A","B"),
          ncol = 1, nrow = 2,
          legend = "none")











############
# Otolith elemental signatures and their relationships with water chemistry
############

#Set directory
setwd("C:/Users/patri/OneDrive/Documentos/UFES/Tese/From water to biominerals/R/github")

#Clean R environment 
rm(list = ls())

data <- read.csv2("otolith_edge.csv")
length(unique(data$ID))
str(data)

#Set as numeric
data <- mutate_at(data, vars(c('Mg24', 'Mn55', 'Cu65', 'Zn66', 'Ba138',
                               'TL', 'Fulton')), as.numeric)

data$Month <- factor(data$Month, levels = c('JUL','AUG','SEP','OCT','NOV','DEC','JAN','FEB','MAR','APR','MAY','JUN'))
summary(data$Month)
length(unique(data$ID))

options(pillar.sigfig=5)



########
# Mean, SD, and range
########

# Mg
#Delete outliers
data_Mg <- data
data_Mg <- data_Mg %>%  filter(!(ID == 'ID335' & Time == '70.2'))
data_Mg <- data_Mg %>%  filter(!(ID == 'ID204' & Time == '88.3'))
data_Mg <- data_Mg %>%  filter(!(ID == 'ID509' & Time == '83.8'))
data_Mg <- data_Mg %>%  filter(!(ID == 'ID496' & Time == '101.9'))
data_Mg <- data_Mg %>%  filter(!(ID == 'ID383' & Time == '92.8'))
data_Mg <- data_Mg %>%  filter(!(ID == 'ID274' & Time == '74.7'))
data_Mg <- data_Mg %>%  filter(!(ID == 'ID511' & Time == '83.8'))
data_Mg <- data_Mg %>%  filter(!(ID == 'ID487' & Time == '58.9'))
data_Mg <- data_Mg %>%  filter(!(ID == 'ID052' &  Time == '52.1'))
data_Mg <- data_Mg %>%  filter(!(ID == 'ID052' & Time == '54.3'))

summary(data_Mg$Mg24)
sd(data_Mg$Mg24)
#The Mg:Caotolith ranged from 180.5 to 353.6 μmol mol−1 (248.5 ± 28.91 μmol mol−1)

# Mn
#Delete outliers
data_Mn <- data
data_Mn <- data[data$ID != c('ID165'), ]
data_Mn <- data_Mn %>%  filter(!(ID == 'ID382' & Time == '81.5'))
data_Mn <- data_Mn %>%  filter(!(ID == 'ID327' & Time == '70.2'))

summary(data_Mn$Mn55)
sd(data_Mn$Mn55)
#Mn:Caotolith ranged from 0.008 to 6.15 μmol mol−1 (1.07 ± 0.92 μmol mol−1)

# Cu
#Delete outliers
data_Cu <- data
data_Cu <- data_Cu %>%  filter(!(ID == 'ID083' & Time == '72.4'))
data_Cu <- data_Cu %>%  filter(!(ID == 'ID072' & Time == '81.5'))
data_Cu <- data_Cu %>%  filter(!(ID == 'ID076' & Time == '77'))
data_Cu <- data_Cu %>%  filter(!(ID == 'ID202' & Time == '70.2'))
data_Cu <- data_Cu %>%  filter(!(ID == 'ID202' & Time == '72.4'))
data_Cu <- data_Cu %>%  filter(!(ID == 'ID277' & Time == '72.4'))

summary(data_Cu$Cu65)
sd(data_Cu$Cu65)
#Cu:Caotolith ranged from 0.003 to 2.11 μmol mol−1 (0.25 ± 0.30 μmol mol−1)

# Zn
#Delete outliers
data_Zn <- data
data_Zn <- data_Zn %>%  filter(!(ID == 'ID386' & Time == '72.4'))
data_Zn <- data_Zn %>%  filter(!(ID == 'ID079' & Time == '72.4'))
data_Zn <- data_Zn %>%  filter(!(ID == 'ID509' & Time == '83.8'))
data_Zn <- data_Zn %>%  filter(!(ID == 'ID511' & Time == '83.8'))
data_Zn <- data_Zn %>%  filter(!(ID == 'ID508' & Time == '47.5'))

summary(data_Zn$Zn66)
sd(data_Zn$Zn66)
#Zn:Caotolith ranged from 0.001 to 7.05 μmol mol−1 (0.75 ± 1.11 μmol mol−1)

# Ba
summary(data$Ba138)
sd(data$Ba138)
#Ba:Caotolith ranged from 0.43 to 161.89 μmol mol−1 (19.32 ± 22.79 μmol mol−1)










########
# Kruskal–Wallis’ test and post hoc  Dunn’s test
########


# Mg
#Assumptions
model <- aov(Mg24 ~ Month, data = data_Mg)
leveneTest(Mg24 ~ Month, data = data_Mg) # 0.001818 **, we cannot assume the homogeneity of variance
shapiro.test(residuals(model)) #p-value = 0.03072, we cannot assume the normality

#Non-parametric test
kruskal.test(Mg24 ~ Month, data = data_Mg) #Kruskal-Wallis chi-squared = 52.335, df = 11, p-value = 2.371e-07
dunnTest(Mg24 ~ Month, data = data_Mg,method="bonferroni") #post hoc


# Mn
#Assumptions
model <- aov(Mn55 ~ Month, data = data_Mn)
leveneTest(Mn55 ~ Month, data = data_Mn) # 0.0006715 ***, we cannot assume the homogeneity of variance
shapiro.test(residuals(model)) #p-value = 2.91e-12, we cannot assume the normality

#Non-parametric test
kruskal.test(Mn55 ~ Month, data = data_Mn) #Kruskal-Wallis chi-squared = 35.545, df = 11, p-value = 0.0002014
dunnTest(Mn55 ~ Month, data = data_Mn,method="bonferroni") #post hoc


# Cu
#Assumptions
model <- aov(Cu65 ~ Month, data = data_Cu)
leveneTest(Cu65 ~ Month, data = data_Cu) # 0.001871 **, we cannot assume the homogeneity of variance
shapiro.test(residuals(model)) #p-value < 2.2e-16, we cannot assume the normality

#Non-parametric test
kruskal.test(Cu65 ~ Month, data = data_Cu) #Kruskal-Wallis chi-squared = 26.199, df = 11, p-value = 0.006061
dunnTest(Cu65 ~ Month, data = data_Cu,method="bonferroni") #post hoc


# Zn
#Assumptions
model <- aov(Zn66 ~ Month, data = data_Zn)
leveneTest(Zn66 ~ Month, data = data_Zn) # 0.1546, we can assume the homogeneity of variance
shapiro.test(residuals(model)) #p-value < 2.2e-16, we cannot assume the normality

#Non-parametric test
kruskal.test(Zn66 ~ Month, data = data_Zn) #Kruskal-Wallis chi-squared = 30.523, df = 11, p-value = 0.001311
dunnTest(Zn66 ~ Month, data = data_Zn,method="bonferroni") #post hoc


# Ba
#Assumptions
model <- aov(Ba138 ~ Month, data = data)
leveneTest(Ba138 ~ Month, data = data) # < 2.2e-16 ***, we cannot assume the homogeneity of variance
shapiro.test(residuals(model)) #p-value < 2.2e-16, we cannot assume the normality

#Non-parametric test
kruskal.test(Ba138 ~ Month, data = data) #Kruskal-Wallis chi-squared = 217.65, df = 11, p-value < 2.2e-16
dunnTest(Ba138 ~ Month, data = data,method="bonferroni") #post hoc












####
# Partition coeficient
####

signal_Mg <- read.csv2("signal_Mg.csv")
signal_Mn <- read.csv2("signal_Mn.csv")
signal_Cu <- read.csv2("signal_Cu.csv")
signal_Zn <- read.csv2("signal_Zn.csv")
signal_Ba <- read.csv2("signal_Ba.csv")

signal_Mg <- signal_Mg %>%
  mutate(Mg24 = as.numeric(Mg24),
    Mg_water = as.numeric(Mg_water))
signal_Mn <- signal_Mn %>%
  mutate(Mn55 = as.numeric(Mn55),
         Mn_water = as.numeric(Mn_water))
signal_Cu <- signal_Cu %>%
  mutate(Cu65 = as.numeric(Cu65),
         Cu_water = as.numeric(Cu_water))
signal_Zn <- signal_Zn %>%
  mutate(Zn66 = as.numeric(Zn66),
         Zn_water = as.numeric(Zn_water))
signal_Ba <- signal_Ba %>%
  mutate(Ba138 = as.numeric(Ba138),
         Ba_water = as.numeric(Ba_water))


signal_Mg$Mg_coef <- signal_Mg$Mg24/signal_Mg$Mg_water
signal_Mn$Mn_coef <- signal_Mn$Mn55/signal_Mn$Mn_water
signal_Cu$Cu_coef <- signal_Cu$Cu65/signal_Cu$Cu_water
signal_Zn$Zn_coef <- signal_Zn$Zn66/signal_Zn$Zn_water
signal_Ba$Ba_coef <- signal_Ba$Ba138/signal_Ba$Ba_water

###
#Mg
###

mean(signal_Mg$Mg_coef)
sd(signal_Mg$Mg_coef)
#D Mg:Ca = 6.85x10-6 ± 8.63x10-6

signal_Mg %>%
  group_by(Month) %>%
  summarise(mean = mean(Mg_coef, na.rm = TRUE),
            sd = sd(Mg_coef, na.rm = TRUE))

###
#Mn
###

signal_Mn %>%
  group_by(Month) %>%
  summarise(mean = mean(Mn_coef, na.rm = TRUE),
            sd = sd(Mn_coef, na.rm = TRUE))

###
#Cu
###

signal_Cu %>%
  group_by(Month) %>%
  summarise(mean = mean(Cu_coef, na.rm = TRUE),
            sd = sd(Cu_coef, na.rm = TRUE))

###
#Zn
###

signal_Zn %>%
  group_by(Month) %>%
  summarise(mean = mean(Zn_coef, na.rm = TRUE),
            sd = sd(Zn_coef, na.rm = TRUE))


###
#Ba
###

mean(signal_Ba$Ba_coef)
sd(signal_Ba$Ba_coef)
#D Ba:Ca = 0.003 ± 0.002

signal_Ba %>%
  group_by(Month) %>%
  summarise(mean = mean(Ba_coef, na.rm = TRUE),
            sd = sd(Ba_coef, na.rm = TRUE))


#END
