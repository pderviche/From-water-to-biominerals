#R Code constructed in July 31, 2024 by Patrick Derviche

#################################################################################### 
#From water to biominerals: what drives the chemical incorporation in the otolith matrix of juvenile snappers?
#################################################################################### 

#Patrick Derviche, Mario V. Condini, Michael A. Dance, Eduardo S. Costa, Fabian Sá, 
#Felippe A. Daros, Maurício Hostim–Silva, Marcelo Soeth

#################################################
#Summary
#################################################
  
  #1. Alignment between otolith and water chemistries, and partition coefficients.......................line 70
    #1.1 Seawater temperature and salinity across months............................................line 80
        #Figure 3
    #1.2 Figure 2b - Scatter plot of laser ablation ICP-MS (ID333)..................................line 170
    #1.3 Searching for the best temporal window.....................................................line 215
        #Element selection - coefficient of variation 
    #1.4 Save dataset with one month correction.....................................................line 1000
    #1.5 Figure S1 - One month correction...........................................................line 1110

  #2. Checking data / Statistical analysis..............................................................line 1180
    #2.1 General results............................................................................line 1240
      #Table 1 - TL and Fulton’s index  of the dog snappers
    #2.2 Dog snappers...............................................................................line 1270
        #Kruskal–Wallis’ nonparametric test, and a post hoc Dunn’s test 
    #2.3 Water chemistry across the saline gradient.................................................line 1320
      #Figure 4 - Water chemistry across the salinity gradient
    #2.4 Otolith elemental signatures and their relationships with water chemistry..................line 1490
      #2.4.1 Mean, SD, and range....................................................................line 1530
        #Kruskal–Wallis’ nonparametric test, and a post hoc Dunn’s test 
      #2.4.2 Figure 5 - Otolith elemental signatures ...............................................line 1650
      #2.4.3 Figure 6 - Correlation between otolith and water chemistries...........................line 1980
    #2.5 Partition coefficients as a proxy for absorption rate......................................line 2210
        #Figure 7 - Partition coefficients 

  #3. GLMMs.............................................................................................line 2650
    #3.1 Correlation between variables..............................................................line 2700
        #Fig. S2. Correlation matrix
    #3.2. Mg........................................................................................line 2755
    #3.3. Mn........................................................................................line 2875
    #3.4. Cu........................................................................................line 2985
    #3.5. Zn........................................................................................line 3100
    #3.6. Ba........................................................................................line 3210
    #3.7. Figure 8 - Significant relationships based on our GLMMs...............................line 3730






















#################################################
#1.Alignment between otolith and water chemistries, and partition coefficients
#################################################

#packages
library(ggplot2)
library(dplyr)



############
#1.1 Seawater temperature and salinity across months
############

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




############
#1.2 Figure 2b
############

#######
#Fig. 3b. Scatter plot of laser ablation transect, indicating the last three 
#laser spots (light green polygon) used in further analysis.
#######

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

#Figure 3
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





############
#1.3 Searching for the best temporal window
############

#Packages
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

#Clean R environment 
rm(list = ls())

#Set working directory
setwd("C:/Users/patri/OneDrive/Documentos/UFES/Tese/From water to biominerals/R/github")

####
# Read data of water chemistry
####

water <- read.csv2("water_chemistry.csv")
water <- as.data.frame(water)
str(water)

#Description of 'water' data
#'Temp' is seawater temperature and 'Sal' is salinity
#'Mg24w', 'Ca43w', 'Mn55w', 'Cu65w', 'Zn66w', and 'Ba138w' are in umol L-1
#'MgCa', 'MnCa', 'CuCa', 'ZnCa', and 'BaCa' are the ratio El:Ca43 in umol mol-1

#Filter water data
water <- water %>%
  filter(tide == "low") %>%
  filter(site == "Sao Mateus")

#Set as numeric
water <- mutate_at(water, vars(c('Mg24w', 'Ca43w', 'Mn55w', 'Cu65w', 'Zn66w', 'Ba138w',
                                 'MgCa', 'MnCa', 'CuCa', 'ZnCa', 'BaCa',
                                 'Temp', 'Sal')), as.numeric)
water$Month<-as.factor(water$Month)

#Change name column 'Month' to 'Sliding'
colnames(water)[colnames(water) == "Month"] <- "Sliding"

#Delete columns in water dataset to not repeat with the otolith dataset
water <- water[,c(-3,-4,-5,-6)]


####
# Read data of otolith chemistry
####

otolith <- read.csv2("otolith_edge.csv")

length(unique(otolith$ID))
#120 individuals
#360 laser spots

####
## No correction dataset (Otolith sampled in JUL correspond to the Water from JUL)
####

#Create the dataset with no correction 
otolith_0 <- otolith
otolith_0$Sliding <- otolith_0$Month #Create the sliding
otolith_0$Order <- otolith_0$Month #Create the Order
otolith_0$OrderS <- otolith_0$Month #Create the OrderS

otolith_0 <- otolith_0[,c(1,4,17:19,2:3,5:16)] #Reorder

otolith_0 <- otolith_0 %>%
  mutate(Order = case_when(
         Order == "JUL" ~ "1", Order == "AUG" ~ "2", Order == "SEP" ~ "3",
         Order == "OCT" ~ "4", Order == "NOV" ~ "5", Order == "DEC" ~ "6",
         Order == "JAN" ~ "7", Order == "FEB" ~ "8", Order == "MAR" ~ "9",
         Order == "APR" ~ "10", Order == "MAY" ~ "11", Order == "JUN" ~ "12",
           TRUE ~ Order))

otolith_0 <- otolith_0 %>%
  mutate(OrderS = case_when(
    OrderS == "JUL" ~ "1", OrderS == "AUG" ~ "2", OrderS == "SEP" ~ "3",
    OrderS == "OCT" ~ "4", OrderS == "NOV" ~ "5", OrderS == "DEC" ~ "6",
    OrderS == "JAN" ~ "7", OrderS == "FEB" ~ "8", OrderS == "MAR" ~ "9",
    OrderS == "APR" ~ "10", OrderS == "MAY" ~ "11", OrderS == "JUN" ~ "12",
    TRUE ~ OrderS))

otolith_0$Month <- factor(otolith_0$Month, levels = c('JUL','AUG','SEP','OCT','NOV','DEC','JAN','FEB','MAR','APR','MAY','JUN'))

length(unique(otolith_0$ID))
table(unique(otolith_0$Month))


####
## One month correction dataset (Otolith sampled in AUG correspond to the Water from JUL)
####

#Create the dataset with no correction 
otolith_1 <- otolith
otolith_1$Sliding <- otolith_1$Month #Create the sliding
otolith_1$Order <- otolith_1$Month #Create the Order
otolith_1$OrderS <- otolith_1$Month #Create the OrderS

otolith_1 <- otolith_1[,c(1,4,17:19,2:3,5:16)] #Reorder

otolith_1 <- otolith_1 %>% filter(!(Month == "JUL")) #Delete Months that you did not sampled the water previously

otolith_1 <- otolith_1 %>%
  mutate(Sliding = case_when(
    Sliding == "AUG" ~ "JUL", Sliding == "SEP" ~ "AUG",
    Sliding == "OCT" ~ "SEP", Sliding == "NOV" ~ "OCT", Sliding == "DEC" ~ "NOV",
    Sliding == "JAN" ~ "DEC", Sliding == "FEB" ~ "JAN", Sliding == "MAR" ~ "FEB",
    Sliding == "APR" ~ "MAR", Sliding == "MAY" ~ "APR", Sliding == "JUN" ~ "MAY",
    TRUE ~ Sliding))

otolith_1 <- otolith_1 %>%
  mutate(Order = case_when(
    Order == "AUG" ~ "2", Order == "SEP" ~ "3",
    Order == "OCT" ~ "4", Order == "NOV" ~ "5", Order == "DEC" ~ "6",
    Order == "JAN" ~ "7", Order == "FEB" ~ "8", Order == "MAR" ~ "9",
    Order == "APR" ~ "10", Order == "MAY" ~ "11", Order == "JUN" ~ "12",
    TRUE ~ Order))

otolith_1 <- otolith_1 %>%
  mutate(OrderS = case_when(
    OrderS == "AUG" ~ "1", OrderS == "SEP" ~ "2",
    OrderS == "OCT" ~ "3", OrderS == "NOV" ~ "4", OrderS == "DEC" ~ "5",
    OrderS == "JAN" ~ "6", OrderS == "FEB" ~ "7", OrderS == "MAR" ~ "8",
    OrderS == "APR" ~ "9", OrderS == "MAY" ~ "10", OrderS == "JUN" ~ "11",
    TRUE ~ OrderS))

otolith_1$Month <- factor(otolith_1$Month, levels = c('AUG','SEP','OCT','NOV','DEC','JAN','FEB','MAR','APR','MAY','JUN'))

length(unique(otolith_1$ID))
table(unique(otolith_1$Month))



####
## Two months correction dataset (Otolith sampled in SEP correspond to the Water from JUL)
####

#Create the dataset with no correction 
otolith_2 <- otolith
otolith_2$Sliding <- otolith_2$Month #Create the sliding
otolith_2$Order <- otolith_2$Month #Create the Order
otolith_2$OrderS <- otolith_2$Month #Create the OrderS

otolith_2 <- otolith_2[,c(1,4,17:19,2:3,5:16)] #Reorder

otolith_2 <- otolith_2 %>% filter(!(Month == "JUL") &
                                  !(Month == "AUG")) #Delete Months that you did not sampled the water previously

otolith_2 <- otolith_2 %>%
  mutate(Sliding = case_when(
    Sliding == "SEP" ~ "JUL", Sliding == "OCT" ~ "AUG",
    Sliding == "NOV" ~ "SEP", Sliding == "DEC" ~ "OCT", Sliding == "JAN" ~ "NOV",
    Sliding == "FEB" ~ "DEC", Sliding == "MAR" ~ "JAN", Sliding == "APR" ~ "FEB",
    Sliding == "MAY" ~ "MAR", Sliding == "JUN" ~ "APR",
    TRUE ~ Sliding))

otolith_2 <- otolith_2 %>%
  mutate(Order = case_when(
    Order == "SEP" ~ "3",
    Order == "OCT" ~ "4", Order == "NOV" ~ "5", Order == "DEC" ~ "6",
    Order == "JAN" ~ "7", Order == "FEB" ~ "8", Order == "MAR" ~ "9",
    Order == "APR" ~ "10", Order == "MAY" ~ "11", Order == "JUN" ~ "12",
    TRUE ~ Order))

otolith_2 <- otolith_2 %>%
  mutate(OrderS = case_when(
    OrderS == "SEP" ~ "1",
    OrderS == "OCT" ~ "2", OrderS == "NOV" ~ "3", OrderS == "DEC" ~ "4",
    OrderS == "JAN" ~ "5", OrderS == "FEB" ~ "6", OrderS == "MAR" ~ "7",
    OrderS == "APR" ~ "8", OrderS == "MAY" ~ "9", OrderS == "JUN" ~ "10",
    TRUE ~ OrderS))

otolith_2$Month <- factor(otolith_2$Month, levels = c('SEP','OCT','NOV','DEC','JAN','FEB','MAR','APR','MAY','JUN'))

length(unique(otolith_2$ID))
table(unique(otolith_2$Month))



####
## Three months correction dataset (Otolith sampled in OCT correspond to the Water from JUL)
####

#Create the dataset with no correction 
otolith_3 <- otolith
otolith_3$Sliding <- otolith_3$Month #Create the sliding
otolith_3$Order <- otolith_3$Month #Create the Order
otolith_3$OrderS <- otolith_3$Month #Create the OrderS

otolith_3 <- otolith_3[,c(1,4,17:19,2:3,5:16)] #Reorder

otolith_3 <- otolith_3 %>% filter(!(Month == "JUL") &
                                  !(Month == "AUG") &
                                  !(Month == "SEP")) #Delete Months that you did not sampled the water previously

otolith_3 <- otolith_3 %>%
  mutate(Sliding = case_when(
    Sliding == "OCT" ~ "JUL",
    Sliding == "NOV" ~ "AUG", Sliding == "DEC" ~ "SEP", Sliding == "JAN" ~ "OCT",
    Sliding == "FEB" ~ "NOV", Sliding == "MAR" ~ "DEC", Sliding == "APR" ~ "JAN",
    Sliding == "MAY" ~ "FEB", Sliding == "JUN" ~ "MAR",
    TRUE ~ Sliding))

otolith_3 <- otolith_3 %>%
  mutate(Order = case_when(
    Order == "OCT" ~ "4", Order == "NOV" ~ "5", Order == "DEC" ~ "6",
    Order == "JAN" ~ "7", Order == "FEB" ~ "8", Order == "MAR" ~ "9",
    Order == "APR" ~ "10", Order == "MAY" ~ "11", Order == "JUN" ~ "12",
    TRUE ~ Order))

otolith_3 <- otolith_3 %>%
  mutate(OrderS = case_when(
    OrderS == "OCT" ~ "1", OrderS == "NOV" ~ "2", OrderS == "DEC" ~ "3",
    OrderS == "JAN" ~ "4", OrderS == "FEB" ~ "5", OrderS == "MAR" ~ "6",
    OrderS == "APR" ~ "7", OrderS == "MAY" ~ "8", OrderS == "JUN" ~ "9",
    TRUE ~ OrderS))

otolith_3$Month <- factor(otolith_3$Month, levels = c('OCT','NOV','DEC','JAN','FEB','MAR','APR','MAY','JUN'))

length(unique(otolith_3$ID))
table(unique(otolith_3$Month))



####
## Four months correction dataset (Otolith sampled in NOV correspond to the Water from JUL)
####

#Create the dataset with no correction 
otolith_4 <- otolith
otolith_4$Sliding <- otolith_4$Month #Create the sliding
otolith_4$Order <- otolith_4$Month #Create the Order
otolith_4$OrderS <- otolith_4$Month #Create the OrderS

otolith_4 <- otolith_4[,c(1,4,17:19,2:3,5:16)] #Reorder

otolith_4 <- otolith_4 %>% filter(!(Month == "JUL") &
                                    !(Month == "AUG") &
                                    !(Month == "SEP") &
                                    !(Month == "OCT")) #Delete Months that you did not sampled the water previously

otolith_4 <- otolith_4 %>%
  mutate(Sliding = case_when(
    Sliding == "NOV" ~ "JUL", Sliding == "DEC" ~ "AUG", Sliding == "JAN" ~ "SEP",
    Sliding == "FEB" ~ "OCT", Sliding == "MAR" ~ "NOV", Sliding == "APR" ~ "DEC",
    Sliding == "MAY" ~ "JAN", Sliding == "JUN" ~ "FEB",
    TRUE ~ Sliding))

otolith_4 <- otolith_4 %>%
  mutate(Order = case_when(
    Order == "NOV" ~ "5", Order == "DEC" ~ "6",
    Order == "JAN" ~ "7", Order == "FEB" ~ "8", Order == "MAR" ~ "9",
    Order == "APR" ~ "10", Order == "MAY" ~ "11", Order == "JUN" ~ "12",
    TRUE ~ Order))

otolith_4 <- otolith_4 %>%
  mutate(OrderS = case_when(
    OrderS == "NOV" ~ "1", OrderS == "DEC" ~ "2",
    OrderS == "JAN" ~ "3", OrderS == "FEB" ~ "4", OrderS == "MAR" ~ "5",
    OrderS == "APR" ~ "6", OrderS == "MAY" ~ "7", OrderS == "JUN" ~ "8",
    TRUE ~ OrderS))

otolith_4$Month <- factor(otolith_4$Month, levels = c('NOV','DEC','JAN','FEB','MAR','APR','MAY','JUN'))

length(unique(otolith_4$ID))
table(unique(otolith_4$Month))



####
## Five months correction dataset (Otolith sampled in DEC correspond to the Water from JUL)
####

#Create the dataset with no correction 
otolith_5 <- otolith
otolith_5$Sliding <- otolith_5$Month #Create the sliding
otolith_5$Order <- otolith_5$Month #Create the Order
otolith_5$OrderS <- otolith_5$Month #Create the OrderS

otolith_5 <- otolith_5[,c(1,4,17:19,2:3,5:16)] #Reorder

otolith_5 <- otolith_5 %>% filter(!(Month == "JUL") &
                                  !(Month == "AUG") &
                                  !(Month == "SEP") &
                                  !(Month == "OCT") &
                                  !(Month == "NOV")) #Delete Months that you did not sampled the water previously

otolith_5 <- otolith_5 %>%
  mutate(Sliding = case_when(
    Sliding == "DEC" ~ "JUL", Sliding == "JAN" ~ "AUG",
    Sliding == "FEB" ~ "SEP", Sliding == "MAR" ~ "OCT", Sliding == "APR" ~ "NOV",
    Sliding == "MAY" ~ "DEC", Sliding == "JUN" ~ "JAN",
    TRUE ~ Sliding))

otolith_5 <- otolith_5 %>%
  mutate(Order = case_when(
    Order == "DEC" ~ "6",
    Order == "JAN" ~ "7", Order == "FEB" ~ "8", Order == "MAR" ~ "9",
    Order == "APR" ~ "10", Order == "MAY" ~ "11", Order == "JUN" ~ "12",
    TRUE ~ Order))

otolith_5 <- otolith_5 %>%
  mutate(OrderS = case_when(
    OrderS == "DEC" ~ "1",
    OrderS == "JAN" ~ "2", OrderS == "FEB" ~ "3", OrderS == "MAR" ~ "4",
    OrderS == "APR" ~ "5", OrderS == "MAY" ~ "6", OrderS == "JUN" ~ "7",
    TRUE ~ OrderS))

otolith_5$Month <- factor(otolith_5$Month, levels = c('DEC','JAN','FEB','MAR','APR','MAY','JUN'))

length(unique(otolith_5$ID))
table(unique(otolith_5$Month))



####
## Six months correction dataset (Otolith sampled in JAN correspond to the Water from JUL)
####

#Create the dataset with no correction 
otolith_6 <- otolith
otolith_6$Sliding <- otolith_6$Month #Create the sliding
otolith_6$Order <- otolith_6$Month #Create the Order
otolith_6$OrderS <- otolith_6$Month #Create the OrderS

otolith_6 <- otolith_6[,c(1,4,17:19,2:3,5:16)] #Reorder

otolith_6 <- otolith_6 %>% filter(!(Month == "JUL") &
                                  !(Month == "AUG") &
                                  !(Month == "SEP") &
                                  !(Month == "OCT") &
                                  !(Month == "NOV")&
                                  !(Month == "DEC")) #Delete Months that you did not sampled the water previously

otolith_6 <- otolith_6 %>%
  mutate(Sliding = case_when(
    Sliding == "JAN" ~ "JUL",
    Sliding == "FEB" ~ "AUG", Sliding == "MAR" ~ "SEP", Sliding == "APR" ~ "OCT",
    Sliding == "MAY" ~ "NOV", Sliding == "JUN" ~ "DEC",
    TRUE ~ Sliding))

otolith_6 <- otolith_6 %>%
  mutate(Order = case_when(
    Order == "JAN" ~ "7", Order == "FEB" ~ "8", Order == "MAR" ~ "9",
    Order == "APR" ~ "10", Order == "MAY" ~ "11", Order == "JUN" ~ "12",
    TRUE ~ Order))

otolith_6 <- otolith_6 %>%
  mutate(OrderS = case_when(
    OrderS == "JAN" ~ "1", OrderS == "FEB" ~ "2", OrderS == "MAR" ~ "3",
    OrderS == "APR" ~ "4", OrderS == "MAY" ~ "5", OrderS == "JUN" ~ "6",
    TRUE ~ OrderS))

otolith_6$Month <- factor(otolith_6$Month, levels = c('JAN','FEB','MAR','APR','MAY','JUN'))

length(unique(otolith_6$ID))
table(unique(otolith_6$Month))



#### Join dataframes
sliding0 = merge(otolith_0, water, on='Sliding', all=FALSE)
str(sliding0)
sliding0$Month <- factor(sliding0$Month, levels = c('JUL','AUG','SEP','OCT','NOV','DEC','JAN','FEB','MAR','APR','MAY','JUN'))
sliding0 <- sliding0[,c(2,3,1,4:33)] #Reorder
length(unique(sliding0$ID))

sliding1 = merge(otolith_1, water, on='Sliding', all=FALSE)
str(sliding1)
sliding1$Month <- factor(sliding1$Month, levels = c('AUG','SEP','OCT','NOV','DEC','JAN','FEB','MAR','APR','MAY','JUN'))
sliding0 <- sliding0[,c(2,3,1,4:33)] #Reorder
length(unique(sliding1$ID))

sliding2 = merge(otolith_2, water, on='Sliding', all=FALSE)
str(sliding2)
sliding2$Month <- factor(sliding2$Month, levels = c('SEP','OCT','NOV','DEC','JAN','FEB','MAR','APR','MAY','JUN'))
sliding0 <- sliding0[,c(2,3,1,4:33)] #Reorder
length(unique(sliding2$ID))

sliding3 = merge(otolith_3, water, on='Sliding', all=FALSE)
str(sliding3)
sliding3$Month <- factor(sliding3$Month, levels = c('OCT','NOV','DEC','JAN','FEB','MAR','APR','MAY','JUN'))
sliding0 <- sliding0[,c(2,3,1,4:33)] #Reorder
length(unique(sliding3$ID))

sliding4 = merge(otolith_4, water, on='Sliding', all=FALSE)
str(sliding4)
sliding4$Month <- factor(sliding4$Month, levels = c('NOV','DEC','JAN','FEB','MAR','APR','MAY','JUN'))
sliding0 <- sliding0[,c(2,3,1,4:33)] #Reorder
length(unique(sliding4$ID))

sliding5 = merge(otolith_5, water, on='Sliding', all=FALSE)
str(sliding5)
sliding5$Month <- factor(sliding5$Month, levels = c('DEC','JAN','FEB','MAR','APR','MAY','JUN'))
sliding0 <- sliding0[,c(2,3,1,4:33)] #Reorder
length(unique(sliding5$ID))

sliding6 = merge(otolith_6, water, on='Sliding', all=FALSE)
str(sliding6)
sliding6$Month <- factor(sliding6$Month, levels = c('JAN','FEB','MAR','APR','MAY','JUN'))
sliding0 <- sliding0[,c(2,3,1,4:33)] #Reorder
length(unique(sliding6$ID))

#Subset the data to compare the GLMMs with the same n size
sliding0 <- subset(sliding0, Month %in% c('JUL','AUG','SEP','OCT','NOV','DEC'))
sliding1 <- subset(sliding1, Month %in% c('AUG','SEP','OCT','NOV','DEC','JAN'))
sliding2 <- subset(sliding2, Month %in% c('SEP','OCT','NOV','DEC','JAN','FEB'))
sliding3 <- subset(sliding3, Month %in% c('OCT','NOV','DEC','JAN','FEB','MAR'))
sliding4 <- subset(sliding4, Month %in% c('NOV','DEC','JAN','FEB','MAR','APR'))
sliding5 <- subset(sliding5, Month %in% c('DEC','JAN','FEB','MAR','APR','MAY'))
sliding6 <- subset(sliding6, Month %in% c('JAN','FEB','MAR','APR','MAY','JUN'))

#sliding0 <- subset(sliding0, Month %in% c('JAN','FEB','MAR','APR','MAY','JUN'))
#sliding1 <- subset(sliding1, Month %in% c('JAN','FEB','MAR','APR','MAY','JUN'))
#sliding2 <- subset(sliding2, Month %in% c('JAN','FEB','MAR','APR','MAY','JUN'))
#sliding3 <- subset(sliding3, Month %in% c('JAN','FEB','MAR','APR','MAY','JUN'))
#sliding4 <- subset(sliding4, Month %in% c('JAN','FEB','MAR','APR','MAY','JUN'))
#sliding5 <- subset(sliding5, Month %in% c('JAN','FEB','MAR','APR','MAY','JUN'))
#sliding6 <- subset(sliding6, Month %in% c('JAN','FEB','MAR','APR','MAY','JUN'))

#Set as numeric

sliding0 <- mutate_at(sliding0, vars(c('Mg24', 'Mn55', 'Cu65', 'Zn66', 'Ba138', 'TL', 'Weight', 'Fulton',
                                       'Mg24w', 'Mn55w',	'Cu65w',	'Zn66w', 'Ba138w',
                                       'MgCa',	'MnCa',	'CuCa',	'ZnCa',	'BaCa',	'Temp', 'Sal')), as.numeric)
sliding0$Month <- factor(sliding0$Month, levels = c('JUL','AUG','SEP','OCT','NOV','DEC'))
sliding0$Sliding <- factor(sliding0$Sliding, levels = c('JUL','AUG','SEP','OCT','NOV','DEC'))
str(sliding0)
summary(sliding0$Month)
summary(sliding0$Sliding)

sliding1 <- mutate_at(sliding1, vars(c('Mg24', 'Mn55', 'Cu65', 'Zn66', 'Ba138', 'TL', 'Weight', 'Fulton',
                                       'Mg24w', 'Mn55w',	'Cu65w',	'Zn66w', 'Ba138w',
                                       'MgCa',	'MnCa',	'CuCa',	'ZnCa',	'BaCa',	'Temp', 'Sal')), as.numeric)
sliding1$Month <- factor(sliding1$Month, levels = c('AUG','SEP','OCT','NOV','DEC','JAN'))
sliding1$Sliding <- factor(sliding1$Sliding, levels = c('JUL','AUG','SEP','OCT','NOV','DEC'))
str(sliding1)
summary(sliding1$Month)
summary(sliding1$Sliding)

sliding2 <- mutate_at(sliding2, vars(c('Mg24', 'Mn55', 'Cu65', 'Zn66', 'Ba138', 'TL', 'Weight', 'Fulton',
                                       'Mg24w', 'Mn55w',	'Cu65w',	'Zn66w', 'Ba138w',
                                       'MgCa',	'MnCa',	'CuCa',	'ZnCa',	'BaCa',	'Temp', 'Sal')), as.numeric)
sliding2$Month <- factor(sliding2$Month, levels = c('SEP','OCT','NOV','DEC','JAN','FEB'))
sliding2$Sliding <- factor(sliding2$Sliding, levels = c('JUL','AUG','SEP','OCT','NOV','DEC'))
str(sliding2)
summary(sliding2$Month)
summary(sliding2$Sliding)

sliding3 <- mutate_at(sliding3, vars(c('Mg24', 'Mn55', 'Cu65', 'Zn66', 'Ba138', 'TL', 'Weight', 'Fulton',
                                       'Mg24w', 'Mn55w',	'Cu65w',	'Zn66w', 'Ba138w',
                                       'MgCa',	'MnCa',	'CuCa',	'ZnCa',	'BaCa',	'Temp', 'Sal')), as.numeric)
sliding3$Month <- factor(sliding3$Month, levels = c('OCT','NOV','DEC','JAN','FEB','MAR'))
sliding3$Sliding <- factor(sliding3$Sliding, levels = c('JUL','AUG','SEP','OCT','NOV','DEC'))
str(sliding3)
summary(sliding3$Month)
summary(sliding3$Sliding)

sliding4 <- mutate_at(sliding4, vars(c('Mg24', 'Mn55', 'Cu65', 'Zn66', 'Ba138', 'TL', 'Weight', 'Fulton',
                                       'Mg24w', 'Mn55w',	'Cu65w',	'Zn66w', 'Ba138w',
                                       'MgCa',	'MnCa',	'CuCa',	'ZnCa',	'BaCa',	'Temp', 'Sal')), as.numeric)
sliding4$Month <- factor(sliding4$Month, levels = c('NOV','DEC','JAN','FEB','MAR','APR'))
sliding4$Sliding <- factor(sliding4$Sliding, levels = c('JUL','AUG','SEP','OCT','NOV','DEC'))
str(sliding4)
summary(sliding4$Month)
summary(sliding4$Sliding)

sliding5 <- mutate_at(sliding5, vars(c('Mg24', 'Mn55', 'Cu65', 'Zn66', 'Ba138', 'TL', 'Weight', 'Fulton',
                                       'Mg24w', 'Mn55w',	'Cu65w',	'Zn66w', 'Ba138w',
                                       'MgCa',	'MnCa',	'CuCa',	'ZnCa',	'BaCa',	'Temp', 'Sal')), as.numeric)
sliding5$Month <- factor(sliding5$Month, levels = c('DEC','JAN','FEB','MAR','APR','MAY'))
sliding5$Sliding <- factor(sliding5$Sliding, levels = c('JUL','AUG','SEP','OCT','NOV','DEC'))
str(sliding5)
summary(sliding5$Month)
summary(sliding5$Sliding)

sliding6 <- mutate_at(sliding6, vars(c('Mg24', 'Mn55', 'Cu65', 'Zn66', 'Ba138', 'TL', 'Weight', 'Fulton',
                                       'Mg24w', 'Mn55w',	'Cu65w',	'Zn66w', 'Ba138w',
                                       'MgCa',	'MnCa',	'CuCa',	'ZnCa',	'BaCa',	'Temp', 'Sal')), as.numeric)
sliding6$Month <- factor(sliding6$Month, levels = c('JAN','FEB','MAR','APR','MAY','JUN'))
sliding6$Sliding <- factor(sliding6$Sliding, levels = c('JUL','AUG','SEP','OCT','NOV','DEC'))
str(sliding6)
summary(sliding6$Month)
summary(sliding6$Sliding)




####
## Select the element with the highest average monthly variability in otolith concentrations
####

#Element selection - coefficient of variation 

#Package
library(goeveg)

otolith <- mutate_at(otolith, vars(c('Mg24', 'Mn55', 'Cu65', 'Zn66', 'Ba138')), as.numeric)

Mg24 <- otolith %>%
  group_by(Month) %>%
  summarise(mean = mean(Mg24) )

Mn55 <- otolith %>%
  group_by(Month) %>%
  summarise(mean = mean(Mn55) )

Cu65 <- otolith %>%
  group_by(Month) %>%
  summarise(mean = mean(Cu65) )

Zn66 <- otolith %>%
  group_by(Month) %>%
  summarise(mean = mean(Zn66) )

Ba138 <- otolith %>%
  group_by(Month) %>%
  summarise(mean = mean(Ba138) )

cv(Mg24$mean, na.rm = TRUE) #cv = 0.11
cv(Mn55$mean, na.rm = TRUE) #cv = 0.38
cv(Cu65$mean, na.rm = TRUE) #cv = 0.30
cv(Zn66$mean, na.rm = TRUE) #cv = 0.52
cv(Ba138$mean, na.rm = TRUE) #cv = 0.90 #Selected element



####
#Ba138
####


#Transform data

box0 = lm(Ba138 ~ BaCa, data = sliding0)
box1 = lm(Ba138 ~ BaCa, data = sliding1)
box2 = lm(Ba138 ~ BaCa, data = sliding2)
box3 = lm(Ba138 ~ BaCa, data = sliding3)
box4 = lm(Ba138 ~ BaCa, data = sliding4)
box5 = lm(Ba138 ~ BaCa, data = sliding5)
box6 = lm(Ba138 ~ BaCa, data = sliding6)

bc0<-boxcox(box0, plotit = TRUE)
lambda0 <- bc0$x[which.max(bc0$y)]

if (lambda0 == 0) {
  sliding0$Ba138t <- log(sliding0$Ba138)
} else {
  sliding0$Ba138t <- (sliding0$Ba138^lambda0 - 1) / lambda0
}

bc1<-boxcox(box1, plotit = TRUE)
lambda1 <- bc1$x[which.max(bc1$y)]

if (lambda1 == 0) {
  sliding1$Ba138t <- log(sliding1$Ba138)
} else {
  sliding1$Ba138t <- (sliding1$Ba138^lambda1 - 1) / lambda1
}

bc2<-boxcox(box2, plotit = TRUE)
lambda2 <- bc2$x[which.max(bc2$y)]

if (lambda2 == 0) {
  sliding2$Ba138t <- log(sliding2$Ba138)
} else {
  sliding2$Ba138t <- (sliding2$Ba138^lambda2 - 1) / lambda2
}

bc3<-boxcox(box3, plotit = TRUE)
lambda3 <- bc3$x[which.max(bc3$y)]

if (lambda3 == 0) {
  sliding3$Ba138t <- log(sliding3$Ba138)
} else {
  sliding3$Ba138t <- (sliding3$Ba138^lambda3 - 1) / lambda3
}

bc4<-boxcox(box4, plotit = TRUE)
lambda4 <- bc4$x[which.max(bc4$y)]

if (lambda4 == 0) {
  sliding4$Ba138t <- log(sliding4$Ba138)
} else {
  sliding4$Ba138t <- (sliding4$Ba138^lambda4 - 1) / lambda4
}

bc5<-boxcox(box5, plotit = TRUE)
lambda5 <- bc5$x[which.max(bc5$y)]

if (lambda5 == 0) {
  sliding5$Ba138t <- log(sliding5$Ba138)
} else {
  sliding5$Ba138t <- (sliding5$Ba138^lambda5 - 1) / lambda5
}

bc6<-boxcox(box6, plotit = TRUE)
lambda6 <- bc6$x[which.max(bc6$y)]

if (lambda6 == 0) {
  sliding6$Ba138t <- log(sliding6$Ba138)
} else {
  sliding6$Ba138t <- (sliding6$Ba138^lambda6 - 1) / lambda6
}

#Scale the predictors
sliding0$Ba138ts <- scale(sliding0$Ba138t)
sliding0$Ba138s <- scale(sliding0$Ba138)
sliding0$BaCas <- scale(sliding0$BaCa)

sliding1$Ba138ts <- scale(sliding1$Ba138t)
sliding1$Ba138s <- scale(sliding1$Ba138)
sliding1$BaCas <- scale(sliding1$BaCa)

sliding2$Ba138ts <- scale(sliding2$Ba138t)
sliding2$Ba138s <- scale(sliding2$Ba138)
sliding2$BaCas <- scale(sliding2$BaCa)

sliding3$Ba138ts <- scale(sliding3$Ba138t)
sliding3$Ba138s <- scale(sliding3$Ba138)
sliding3$BaCas <- scale(sliding3$BaCa)

sliding4$Ba138ts <- scale(sliding4$Ba138t)
sliding4$Ba138s <- scale(sliding4$Ba138)
sliding4$BaCas <- scale(sliding4$BaCa)

sliding5$Ba138ts <- scale(sliding5$Ba138t)
sliding5$Ba138s <- scale(sliding5$Ba138)
sliding5$BaCas <- scale(sliding5$BaCa)

sliding6$Ba138ts <- scale(sliding6$Ba138t)
sliding6$Ba138s <- scale(sliding6$Ba138)
sliding6$BaCas <- scale(sliding6$BaCa)


########
#GLMMs
#Response variable  transformed
########

model_sw0t <- lmer(Ba138ts ~ BaCas + (1 | ID), 
                  data = sliding0, REML = FALSE,
                  na.action =  na.fail)

model_sw1t <- lmer(Ba138ts ~ BaCas + (1 | ID), 
                  data = sliding1, REML = FALSE,
                  na.action =  na.fail)

model_sw2t <- lmer(Ba138ts ~ BaCas + (1 | ID), 
                  data = sliding2, REML = FALSE,
                  na.action =  na.fail)

model_sw3t <- lmer(Ba138ts ~ BaCas + (1 | ID), 
                  data = sliding3, REML = FALSE,
                  na.action =  na.fail) 

model_sw4t <- lmer(Ba138ts ~ BaCas + (1 | ID), 
                  data = sliding4, REML = FALSE,
                  na.action =  na.fail)

model_sw5t <- lmer(Ba138ts ~ BaCas + (1 | ID), 
                  data = sliding5, REML = FALSE,
                  na.action =  na.fail)

model_sw6t <- lmer(Ba138ts ~ BaCas + (1 | ID), 
                  data = sliding6, REML = FALSE,
                  na.action =  na.fail)

#Sliding selection
sel_swt <- model.sel(model_sw0t, model_sw1t, model_sw2t, model_sw3t, model_sw4t, model_sw5t, model_sw6t)
sel_swt

#Model selection table 
#(Intrc)   BaCas df   logLik  AICc  delta weight
#model_sw1t  9.743e-16  0.6355  4 -111.954 232.2   0.00  0.774
#model_sw3t  1.217e-15  0.6044  4 -113.184 234.6   2.46  0.226
#model_sw2t  3.324e-15  0.6243  4 -121.607 251.5  19.31  0.000
#model_sw4t -1.647e-15  0.3189  4 -132.391 273.0  40.87  0.000
#model_sw5t  1.094e-15 -0.1201  4 -157.435 323.1  90.95  0.000
#model_sw6t  1.021e-15 -0.3026  4 -179.162 366.5 134.38  0.000
#model_sw0t -3.891e-16 -0.3051  4 -196.257 400.8 168.61  0.000
#Models ranked by AICc(x) 
#Random terms (all models): 
#  1 | ID

#r2 values
r.squaredGLMM(model_sw1t)

#R2m       R2c
#0.4053384 0.914108

#Profiled Confidence Intervals
profile_cit <- ggpredict(model_sw1t, terms = c("BaCas"))

plot(profile_cit) +
  labs(title = "One month corrected",
       x = expression(Ba:Ca~water), y = expression(Ba:Ca~otolith), color = "")+
  theme(legend.position="none",
        title = element_text(size = 15),
        axis.text = element_text(size = 15),
        axis.title = element_text(size = 15))

########
#GLMMs
#Response variable not transformed
########

model_sw0 <- lmer(Ba138s ~ BaCas + (1 | ID), 
                  data = sliding0, REML = FALSE,
                  na.action =  na.fail)

model_sw1 <- lmer(Ba138s ~ BaCas + (1 | ID), 
                  data = sliding1, REML = FALSE,
                  na.action =  na.fail)

model_sw2 <- lmer(Ba138s ~ BaCas + (1 | ID), 
                  data = sliding2, REML = FALSE,
                  na.action =  na.fail)

model_sw3 <- lmer(Ba138s ~ BaCas + (1 | ID), 
                  data = sliding3, REML = FALSE,
                  na.action =  na.fail) 

model_sw4 <- lmer(Ba138s ~ BaCas + (1 | ID), 
                  data = sliding4, REML = FALSE,
                  na.action =  na.fail)

model_sw5 <- lmer(Ba138s ~ BaCas + (1 | ID), 
                  data = sliding5, REML = FALSE,
                  na.action =  na.fail)

model_sw6 <- lmer(Ba138s ~ BaCas + (1 | ID), 
                  data = sliding6, REML = FALSE,
                  na.action =  na.fail)

#Sliding selection
sel_sw <- model.sel(model_sw0, model_sw1, model_sw2, model_sw3, model_sw4, model_sw5, model_sw6)
sel_sw

#Model selection table 
#           (Intrc)    BaCas    df   logLik  AICc  delta weight
#model_sw1  9.961e-15  0.7354  4  -67.172 142.6   0.00      1
#model_sw2 -9.538e-16  0.5796  4 -142.203 292.7 150.06      0
#model_sw3  2.265e-15  0.6028  4 -162.237 332.7 190.13      0
#model_sw4 -5.984e-16  0.2610  4 -179.482 367.2 224.62      0
#model_sw0 -1.402e-16 -0.2712  4 -186.560 381.4 238.78      0
#model_sw5 -3.701e-16 -0.1434  4 -191.249 390.7 248.14      0
#model_sw6 -2.158e-15 -0.2471  4 -211.196 430.6 288.01      0
#Models ranked by AICc(x) 
#Random terms (all models): 
#  1 | ID

#r2 values
r.squaredGLMM(model_sw1)

#R2m       R2c
#0.5423579 0.9574821

#Profiled Confidence Intervals
profile_ci <- ggpredict(model_sw1, terms = c("BaCas"))

plot(profile_ci) +
  labs(title = "One month corrected",
       x = expression(Ba:Ca~water), y = expression(Ba:Ca~otolith), color = "")+
  theme(legend.position="none",
        title = element_text(size = 15),
        axis.text = element_text(size = 15),
        axis.title = element_text(size = 15))








############
# 1.4 Save dataset with one month correction
############


#Clean R environment 
rm(list = ls())

#Set working directory
setwd("C:/Users/patri/OneDrive/Documentos/UFES/Tese/From water to biominerals/R/github")

####
# Read data of water chemistry
####

water <- read.csv2("water_chemistry.csv")
water <- as.data.frame(water)
str(water)

#Description of 'water' data
#'Temp' is seawater temperature and 'Sal' is salinity
#'Mg24w', 'Ca43w', 'Mn55w', 'Cu65w', 'Zn66w', and 'Ba138w' are in umol L-1
#'MgCa', 'MnCa', 'CuCa', 'ZnCa', and 'BaCa' are the ratio El:Ca43 in umol mol-1

#Filter water data
water <- water %>%
  filter(tide == "low") %>%
  filter(site == "Sao Mateus")

#Set as numeric
water <- mutate_at(water, vars(c('Mg24w', 'Ca43w', 'Mn55w', 'Cu65w', 'Zn66w', 'Ba138w',
                                 'MgCa', 'MnCa', 'CuCa', 'ZnCa',  'BaCa',
                                 'Temp', 'Sal')), as.numeric)
water$Month<-as.factor(water$Month)

#Change name column 'Month' to 'Sliding'
colnames(water)[colnames(water) == "Month"] <- "Sliding"

#Delete columns in water dataset to not repeat with the otolith dataset
water<-water[,c(-3,-4,-5,-6)]


####
# Read data of otolith chemistry
####

otolith <- read.csv2("otolith_edge.csv")

length(unique(otolith$ID))
#120 individuals
#360 laser spots

####
## One month correction dataset (Otolith sampled in AUG correspond to the Water from JUL)
####

#Create the dataset with no correction 
otolith_1 <- otolith
otolith_1$Sliding <- otolith_1$Month #Create the sliding
otolith_1$Order <- otolith_1$Month #Create the Order
otolith_1$OrderS <- otolith_1$Month #Create the OrderS

otolith_1 <- otolith_1[,c(1,4,17:19,2:3,5:16)] #Reorder

otolith_1 <- otolith_1 %>% filter(!(Month == "JUL")) #Delete Months that you did not sampled the water previously

otolith_1 <- otolith_1 %>%
  mutate(Sliding = case_when(
    Sliding == "AUG" ~ "JUL", Sliding == "SEP" ~ "AUG",
    Sliding == "OCT" ~ "SEP", Sliding == "NOV" ~ "OCT", Sliding == "DEC" ~ "NOV",
    Sliding == "JAN" ~ "DEC", Sliding == "FEB" ~ "JAN", Sliding == "MAR" ~ "FEB",
    Sliding == "APR" ~ "MAR", Sliding == "MAY" ~ "APR", Sliding == "JUN" ~ "MAY",
    TRUE ~ Sliding))

otolith_1 <- otolith_1 %>%
  mutate(Order = case_when(
    Order == "AUG" ~ "2", Order == "SEP" ~ "3",
    Order == "OCT" ~ "4", Order == "NOV" ~ "5", Order == "DEC" ~ "6",
    Order == "JAN" ~ "7", Order == "FEB" ~ "8", Order == "MAR" ~ "9",
    Order == "APR" ~ "10", Order == "MAY" ~ "11", Order == "JUN" ~ "12",
    TRUE ~ Order))

otolith_1 <- otolith_1 %>%
  mutate(OrderS = case_when(
    OrderS == "AUG" ~ "1", OrderS == "SEP" ~ "2",
    OrderS == "OCT" ~ "3", OrderS == "NOV" ~ "4", OrderS == "DEC" ~ "5",
    OrderS == "JAN" ~ "6", OrderS == "FEB" ~ "7", OrderS == "MAR" ~ "8",
    OrderS == "APR" ~ "9", OrderS == "MAY" ~ "10", OrderS == "JUN" ~ "11",
    TRUE ~ OrderS))

otolith_1$Month <- factor(otolith_1$Month, levels = c('AUG','SEP','OCT','NOV','DEC','JAN','FEB','MAR','APR','MAY','JUN'))

length(unique(otolith_1$ID))
table(unique(otolith_1$Month))

#Merge
data_m1_e3 = merge(otolith_1, water, on='Sliding', all=FALSE)
str(data_m1_e3)
data_m1_e3$Month <- factor(data_m1_e3$Month, levels = c('AUG','SEP','OCT','NOV','DEC','JAN','FEB','MAR','APR','MAY','JUN'))
data_m1_e3$OrderS <-as.numeric(data_m1_e3$OrderS)
data_m1_e3 <- data_m1_e3[,c(2,3,1,4:33)] #Reorder
length(unique(data_m1_e3$ID))

data_m1_e3 <- data_m1_e3 %>%
  arrange(OrderS, ID)

#write.table(data_m1_e3,"data_m1_e3.csv", sep=";", dec=".",row.names = F) #Month 1 correction (m1), Edge 3 otolith (m3)



############
#1.5 Figure S1 - One month correction
############

#Clean R environment 
rm(list = ls())

#Set working directory
setwd("C:/Users/patri/OneDrive/Documentos/UFES/Tese/From water to biominerals/R/github")

#read dataset
data <- read.csv2("window.csv")

#Description of 'water' data
#'Ba138' is the otolith signature (Ba:Ca) expressed in umol mol-1
#'BaCa' is the water chemistry (Ba:Ca) corrected in one month and expressed  umol mol-1
#'BaCa_no_correction' is the water chemistry (Ba:Ca) without correction and expressed in umol mol-1

data <- mutate_at(data, vars(c('Ba138','BaCa','BaCa_no_correction','Order','OrderS')), as.numeric)

data$Month <- factor(data$Month, levels = c('JUL','AUG','SEP','OCT','NOV','DEC','JAN','FEB','MAR','APR','MAY','JUN'))
data$Sliding <- factor(data$Sliding, levels = c('JUN','JUL','AUG','SEP','OCT','NOV','DEC','JAN','FEB','MAR','APR','MAY'))

data$BaCa_n_m <- data$BaCa_no_correction/1000 #transform umol mol-1 to mmol mol-1
data$BaCa_m <- data$BaCa/1000 #transform umol mol-1 to mmol mol-1

scaleBa<- 1

############
#Figure S1. Measured otolith elemental signatures (green line, Ba:Caotolith), water chemistry (Ba:Cawater) 
#without correction (blue dashed line), and Ba:Cawater with the best temporal correction (one month corrected, 
#blue solid line). Temporal corrections resulted in the Ba:Cawater value corresponding to a calendar date of one
#month prior to the Ba:Caotolith sampling (e.g., the Ba:Caotolith collected in August refers to the Ba:Cawater 
#sampled in July). Estimated scatterplot smoothing curves (lines) and 95% confidence intervals (shaded areas) 
#were done by applying the local polynomial regression fitting (loess).
############


correction <- ggplot(data = data, aes(y = Ba138,x = Order))+
  geom_smooth(data = data, aes(y = Ba138,x = Order), color = "#008080",fill = "#008080",alpha=0.3) +
  geom_smooth(data = data, aes(y = BaCa_m*scaleBa, x = Order),color = "#00BFFF", size = 0.5,alpha=0.2, method = 'loess')+
  geom_smooth(data = data, aes(y = BaCa_n_m*scaleBa, x = Order),color = "#00BFFF", size = 0.5,alpha=0.2,linetype = "dashed", method = 'loess')+
  theme_bw()+
  scale_y_continuous(sec.axis = sec_axis(~./scaleBa, name=(expression(Ba:Ca[water]~(mmol~mol^-1)))))+
  scale_x_continuous(breaks=seq(1, 13, 1),labels = c('JUL','AUG','SEP','OCT','NOV','DEC','JAN','FEB','MAR','APR','MAY','JUN','JUL')) +
  labs(x = " ", y = (expression(Ba:Ca[otolith]~(mu~mol~mol^-1))), color = "")+
  theme(axis.title.y = element_text(color = "#008080",size = 15),
        axis.title.y.right = element_text(color = "#00BFFF"))+
  theme(legend.position="none",
        title = element_text(size = 15),
        axis.text = element_text(size = 15),
        axis.text.x = element_text(angle = -45),
        axis.title = element_text(size = 15))+ 
  annotate("text", x = 3.5, y = 40, label = "One month corrected", size = 5.5, color = "black")

correction
#600 x 500













#################################################
# 2. Checking data / Statistical analysis
#################################################

#packages
library(FSA)
library(multcompView)
library(rcompanion)
library (tidyverse)
library(car)
library(ggpubr)
library(caret)
library(ggplot2)
library(tidyr)
library(patchwork)
library(rcompanion)

############
#Read data set
############

#Set directory
setwd("C:/Users/patri/OneDrive/Documentos/UFES/Tese/From water to biominerals/R/github")

#Clean R environment 
rm(list = ls())

data <- read.csv2("data_m1_e3.csv")
length(unique(data$ID))
str(data)

###
#Legend
###

#Otolith signatures, ratio El:Ca (umol mol-1): Mg24,	Mn55,	Cu65,	Zn66,	Ba138
#Water chemistry, El concentration (umol L-1): Mg24w,	Ca43w,	Mn55w,	Cu65w,	Zn66w,	Ba138w
#Water chemistry, ratio El:Ca (umol mol): MgCa,	MnCa,	CuCa,	ZnCa,	BaCa

#Set as numeric
data <- mutate_at(data, vars(c('Mg24', 'Mn55', 'Cu65', 'Zn66', 'Ba138',
                               'TL', 'Fulton',
                               'MgCa',	'MnCa',	'CuCa',	'ZnCa', 'BaCa',
                               'Temp', 'Sal','OrderS')), as.numeric)

data$Month <- factor(data$Month, levels = c('JUL','AUG','SEP','OCT','NOV','DEC','JAN','FEB','MAR','APR','MAY','JUN'))
data$Sliding <- factor(data$Sliding, levels = c('JUL','AUG','SEP','OCT','NOV','DEC','JAN','FEB','MAR','APR','MAY','JUN'))
summary(data$Month)
summary(data$Sliding)
length(unique(data$ID))

str(data)
options(pillar.sigfig=4)







############
# 2.1 General results
############

#Table 1. Sample size (n), mean ± standard deviation, range (minimum and maximum) of total length (TL, in mm), 
#and Fulton’s index (FI), for juvenile dog snapper individuals monthly collected at the São Mateus estuary –
#Southwestern Atlantic.

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
# 2.2 Dog snappers
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








############
# 2.3 Water chemistry across the saline gradient
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


###
# Fig. 4. Water chemistry (El:Cawater) across the salinity gradient (from zero to 35) collected during 
#the early wet (green line) and late wet (blue line) seasons at the São Mateus estuary – Southwestern Atlantic. 
#The dots indicate raw data, while estimated curves (lines) were done by applying linear models using a third–degree basis spline.
###

library(ggpubr)

ggarrange(Mg,Mn,Cu,Zn,Ba, 
          labels = c("A","B","C","D","E"),
          ncol = 2, nrow = 3,
          legend = "none")
#700 x 700








############
# 2.4. Otolith elemental signatures and their relationships with water chemistry
############

#Set directory
setwd("C:/Users/patri/OneDrive/Documentos/UFES/Tese/From water to biominerals/R/github")

#Clean R environment 
rm(list = ls())

data <- read.csv2("data_m1_e3.csv")
length(unique(data$ID))
str(data)

#Set as numeric
data <- mutate_at(data, vars(c('Mg24', 'Mn55', 'Cu65', 'Zn66', 'Ba138',
                               'TL', 'Fulton',
                               'MgCa',	'MnCa',	'CuCa',	'ZnCa',	'BaCa',
                               'Temp', 'Sal','OrderS')), as.numeric)

data$Month <- factor(data$Month, levels = c('JUL','AUG','SEP','OCT','NOV','DEC','JAN','FEB','MAR','APR','MAY','JUN'))
data$Sliding <- factor(data$Sliding, levels = c('JUL','AUG','SEP','OCT','NOV','DEC','JAN','FEB','MAR','APR','MAY','JUN'))
summary(data$Month)
summary(data$Sliding)
length(unique(data$ID))

options(pillar.sigfig=5)



########
# 2.4.1 Mean, SD, and range
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





########
# 2.4.2 Figure 5 - Otolith elemental signatures
########

#Set directory
setwd("C:/Users/patri/OneDrive/Documentos/UFES/Tese/From water to biominerals/R/github")

#Clean R environment 
rm(list = ls())

#Read data
data <- read.csv2("data_m1_e3.csv")
length(unique(data$ID))
str(data)

#Set as numeric
data <- mutate_at(data, vars(c('Mg24', 'Mn55', 'Cu65','Zn66',  'Ba138', 'TL', 'Fulton',
                               'MgCa',	'MnCa',	'CuCa',	'ZnCa', 'BaCa',	
                               'Temp', 'Sal','OrderS')), as.numeric)

data$Month <- factor(data$Month, levels = c('JUL','AUG','SEP','OCT','NOV','DEC','JAN','FEB','MAR','APR','MAY','JUN'))
data$Sliding <- factor(data$Sliding, levels = c('JUL','AUG','SEP','OCT','NOV','DEC','JAN','FEB','MAR','APR','MAY'))
summary(data$Month)
summary(data$Sliding)

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
data_Mg$MgCam <- data_Mg$MgCa/1000000 

#set scale
scaleMg<- 7

#letters based on post hoc Dunn's test
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
  geom_smooth(data = data_Mg, aes(y = MgCam*scaleMg, x = Order),
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
data_Mn$MnCam <- data_Mn$MnCa/1000 

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
  geom_smooth(data = data_Mn, aes(y =  MnCam/scaleMn, x = Order),
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
  geom_smooth(data = data_Cu, aes(y =  CuCa/scaleCu, x = Order),
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
data_Zn$ZnCam <- data_Zn$ZnCa/1000 

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
  geom_smooth(data = data_Zn, aes(y = ZnCam*scaleZn, x = Order),
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
data_Ba$BaCam <- data$BaCa/1000 

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
  geom_smooth(data = data_Ba, aes(y = BaCam*scaleBa, x = Order),
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












########
# 2.4.3 Figure 7 - Correlation between otolith and water chemistries
########

#Clean R environment
rm(list = ls())

#Read data
data <- read.csv2("data_m1_e3.csv")
length(unique(data$ID))
str(data)

#Set as numeric
data <- mutate_at(data, vars(c('Mg24', 'Mn55', 'Cu65', 'Zn66', 'Ba138', 'TL', 'Fulton',
                               'Mg24w','Mn55w',	'Cu65w', 'Zn66w', 'Ba138w',
                               'MgCa','MnCa',	'CuCa',	'ZnCa',	'BaCa','Temp', 'Sal','OrderS')), as.numeric)

data$Month <- factor(data$Month, levels = c('JUL','AUG','SEP','OCT','NOV','DEC','JAN','FEB','MAR','APR','MAY','JUN'))
data <- data[data$Month != c('JUL'), ]

length(unique(data$ID))

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

#Get the mean
mean_Mg <- data_Mg %>%
  group_by(Month) %>%
  summarise(Mg24 = mean(Mg24),
            MgCa = mean(MgCa))
str(mean_Mg)

###
# Mn:Ca
###

#Delete outliers
data_Mn <- data
data_Mn <- data[data$ID != c('ID165'), ]
data_Mn <- data_Mn %>%  filter(!(ID == 'ID382' & Time == '81.5'))
data_Mn <- data_Mn %>%  filter(!(ID == 'ID327' & Time == '70.2'))

#Get the mean
mean_Mn <- data_Mn %>%
  group_by(Month) %>%
  summarise(Mn55 = mean(Mn55),
            MnCa = mean(MnCa))
str(mean_Mn)

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

#Get the mean
mean_Cu <- data_Cu %>%
  group_by(Month) %>%
  summarise(Cu65 = mean(Cu65),
            CuCa = mean(CuCa))
str(mean_Cu)

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

#Get the mean
mean_Zn <- data_Zn %>%
  group_by(Month) %>%
  summarise(Zn66 = mean(Zn66),
            ZnCa = mean(ZnCa))
str(mean_Zn)

###
# Ba:Ca
###

#Get the mean
mean_Ba <- data %>%
  group_by(Month) %>%
  summarise(Ba138 = mean(Ba138),
            BaCa = mean(BaCa))
str(mean_Ba)

######
#Figures
######

###
# Mg:Ca
###

mean_Mg$MgCam <- mean_Mg$MgCa/1000000 #transform umol mol-1 to mol mol-1

Mg <-ggplot(mean_Mg, aes(y=Mg24, x=MgCam)) +
  theme_bw() +
  geom_smooth(method=glm,color = "#00BFFF", fill = "#00BFFF",alpha=0.1)  +
  geom_point(size=3) +
  labs(x = expression(Mg:Ca[water]~ ('mol'~'mol'^-1)), y = (expression(Mg:Ca[otolith]~(mu~mol~mol^-1))), color = "")+
  theme(legend.position="none",
        title = element_text(size = 15),
        axis.text = element_text(size = 15),
        axis.title = element_text(size = 15))+ 
  stat_cor(method = "pearson")
Mg


###
# Mn:Ca
###

mean_Mn$MnCam <- mean_Mn$MnCa/1000 #transform umol mol-1 to mmol mol-1

Mn <-ggplot(mean_Mn, aes(y=Mn55, x=MnCam)) +
  theme_bw() +
  geom_smooth(method=glm,color = "#00BFFF", fill = "#00BFFF",alpha=0.1)  +
  geom_point(size=3) +
  labs(x = expression(Mn:Ca[water]~ (mmol~mol^-1)), y = (expression(Mn:Ca[otolith]~(mu~mol~mol^-1))), color = "")+
  theme(legend.position="none",
        title = element_text(size = 15),
        axis.text = element_text(size = 15),
        axis.title = element_text(size = 15))+ 
  stat_cor(method = "pearson")
Mn


###
# Cu:Ca
###

Cu <-ggplot(mean_Cu, aes(y=Cu65, x=CuCa)) +
  theme_bw() +
  geom_smooth(method=glm,color = "#00BFFF", fill = "#00BFFF",alpha=0.1)  +
  geom_point(size=3) +
  labs(x = expression(Cu:Ca[water]~ (mu~mol~mol^-1)), y = (expression(Cu:Ca[otolith]~(mu~mol~mol^-1))), color = "")+
  theme(legend.position="none",
        title = element_text(size = 15),
        axis.text = element_text(size = 15),
        axis.title = element_text(size = 15))+ 
  stat_cor(method = "pearson")
Cu


###
# Zn:Ca
###

mean_Zn$ZnCam <- mean_Zn$ZnCa/1000 #transform umol mol-1 to mmol mol-1

Zn <-ggplot(mean_Zn, aes(y=Zn66, x=ZnCam)) +
  theme_bw() +
  geom_smooth(method=glm,color = "#00BFFF", fill = "#00BFFF",alpha=0.1)  +
  geom_point(size=3) +
  labs(x = expression(Zn:Ca[water]~ (mmol~mol^-1)), y = (expression(Zn:Ca[otolith]~(mu~mol~mol^-1))), color = "")+
  theme(legend.position="none",
        title = element_text(size = 15),
        axis.text = element_text(size = 15),
        axis.title = element_text(size = 15))+ 
  stat_cor(method = "pearson",label.y = 1.6)
Zn


###
# Ba:Ca
###

mean_Ba$BaCam <- mean_Ba$BaCa/1000 #transform umol mol-1 to mmol mol-1

Ba <- ggplot(mean_Ba, aes(y=Ba138, x=BaCam)) +
  theme_bw() +
  geom_smooth(method=glm,color = "#00BFFF", fill = "#00BFFF",alpha=0.1)  +
  geom_point(size=3) +
  labs(x = expression(Ba:Ca[water]~ (mmol~mol^-1)), y = (expression(Ba:Ca[otolith]~(mu~mol~mol^-1))), color = "")+
  theme(legend.position="none",
        title = element_text(size = 15),
        axis.text = element_text(size = 15),
        axis.title = element_text(size = 15))+ 
  stat_cor(method = "pearson",label.y = 65)
Ba


####
# Fig. 6. Correlation between the mean otolith elemental signatures (El:Caotolith) of the juvenile dog
#snappers and water chemistry (El:Cawater) monthly sampled at the São Mateus estuary – Southwestern Atlantic. 
#Pearson’s correlation (R) and p values of the relationships are indicated at the top left for each element. 
#The correlation was calculated considering El:Cawater with a one–month correction.
####

ggarrange(Mg,Mn,Cu,Zn,Ba, 
          labels = c("A","B","C","D","E"),
          ncol = 2, nrow = 3,
          legend = "none")
#1000 x 1000








########
# 2.5 Partition coefficients - D as a proxy for absorption rate
########

#Set working directory
setwd("C:/Users/patri/OneDrive/Documentos/UFES/Tese/From water to biominerals/R/github")

#Clean R environment
rm(list = ls())

#Read data
data <- read.csv2("data_m1_e3.csv")
length(unique(data$ID))
str(data)

data <- data %>%
  mutate(Season = dplyr::recode(Season,
                                "Late dry" = "LD",
                                "Early wet" = "EW",
                                "Late wet" = "LW",
                                "Early dry" = "ED"))

data$Season <- factor(data$Season, levels = c('LD','EW','LW','ED'))

#Set as numeric
data <- mutate_at(data, vars(c('Mg24', 'Mn55', 'Cu65', 'Zn66', 'Ba138', 'TL', 'Fulton',
                               'Mg24w','Mn55w',	'Cu65w', 'Zn66w', 'Ba138w',
                               'MgCa','MnCa',	'CuCa',	'ZnCa',	'BaCa','Temp', 'Sal','OrderS')), as.numeric)

data$Month <- factor(data$Month, levels = c('JUL','AUG','SEP','OCT','NOV','DEC','JAN','FEB','MAR','APR','MAY','JUN'))
data$Sliding <- factor(data$Sliding, levels = c('JUL','AUG','SEP','OCT','NOV','DEC','JAN','FEB','MAR','APR','MAY'))
summary(data$Month)
summary(data$Sliding)

data <- data[data$Month != c('JUL'), ]

length(unique(data$ID))

###
# Partition coefficients of elemental ratios (D) 
###

data$Mg_coef <- data$Mg24/data$MgCa
data$Mn_coef <- data$Mn55/data$MnCa
data$Cu_coef <- data$Cu65/data$CuCa
data$Zn_coef <- data$Zn66/data$ZnCa
data$Ba_coef <- data$Ba138/data$BaCa



###
# Summary Coeficients
###

#Mg_coef
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

summary(data_Mg$Mg_coef)
sd(data_Mg$Mg_coef)
#D Mg:Ca ratios ranged from 0.000003021 to 0.00002393 (0.000007919 ± 0.00000481)

#Mn_coef
data_Mn <- data
data_Mn <- data[data_Mn$ID != c('ID165'), ]
data_Mn <- data_Mn %>%  filter(!(ID == 'ID382' & Time == '81.5'))
data_Mn <- data_Mn %>%  filter(!(ID == 'ID327' & Time == '70.2'))

summary(data_Mn$Mn_coef)
sd(data_Mn$Mn_coef)
#D Mn:Ca ratios ranged from 0.000937 to 0.003562 (0.0003583 ± 0.0005)

#Cu_coef
data_Cu <- data
data_Cu <- data_Cu %>%  filter(!(ID == 'ID083' & Time == '72.4'))
data_Cu <- data_Cu %>%  filter(!(ID == 'ID072' & Time == '81.5'))
data_Cu <- data_Cu %>%  filter(!(ID == 'ID076' & Time == '77'))
data_Cu <- data_Cu %>%  filter(!(ID == 'ID202' & Time == '70.2'))
data_Cu <- data_Cu %>%  filter(!(ID == 'ID202' & Time == '72.4'))
data_Cu <- data_Cu %>%  filter(!(ID == 'ID277' & Time == '72.4'))
data_Cu$Cu_coef[is.infinite(data_Cu$Cu_coef)] <- NA
data_Cu <- data_Cu[!is.na(data_Cu$Cu_coef), ]

summary(data_Cu$Cu_coef)
sd(data_Cu$Cu_coef)
#D Cu:Ca ratios ranged from 0.00003216 to 0.05019 (0.004082 ± 0.007)

#Zn_coef
data_Zn <- data
data_Zn <- data_Zn %>%  filter(!(ID == 'ID386' & Time == '72.4'))
data_Zn <- data_Zn %>%  filter(!(ID == 'ID079' & Time == '72.4'))
data_Zn <- data_Zn %>%  filter(!(ID == 'ID509' & Time == '83.8'))
data_Zn <- data_Zn %>%  filter(!(ID == 'ID511' & Time == '83.8'))
data_Zn <- data_Zn %>%  filter(!(ID == 'ID508' & Time == '47.5'))

summary(data_Zn$Zn_coef)
sd(data_Zn$Zn_coef)
#D Zn:Ca ratios ranged from 5.8e-07 to 0.1223 (0.002846 ± 0.008)

#Ba_coef
data_Ba <- data

summary(data_Ba$Ba_coef)
sd(data_Ba$Ba_coef)
#D Ba:Ca ratios ranged from 0.00008705 to 0.1234 (0.01032 ± 0.015)

mean_Mg <- mean(data_Mg$Mg_coef, na.rm = TRUE)
mean_Mn <- mean(data_Mn$Mn_coef, na.rm = TRUE)
mean_Cu <- mean(data_Cu$Cu_coef, na.rm = TRUE)
mean_Zn <- mean(data_Zn$Zn_coef, na.rm = TRUE)
mean_Ba <- mean(data_Ba$Ba_coef, na.rm = TRUE)

mean_list <- list(
  Mg_coef = mean_Mg,
  Mn_coef = mean_Mn,
  Cu_coef = mean_Cu,
  Zn_coef = mean_Zn,
  Ba_coef = mean_Ba
)

mean_list <- mean_list[order(unlist(mean_list), decreasing = TRUE)]

mean_list


###
# Partition coeficients (Mean ± SD) in each season
###

data_Mg %>%
  group_by(Season) %>%
  summarise(
    mean_Mg_coef = mean(Mg_coef, na.rm = TRUE),
    sd_Mg_coef = sd(Mg_coef, na.rm = TRUE))

data_Mn %>%
  group_by(Season) %>%
  summarise(
    mean_Mn_coef = mean(Mn_coef, na.rm = TRUE),
    sd_Mn_coef = sd(Mn_coef, na.rm = TRUE))

data_Cu %>%
  group_by(Season) %>%
  summarise(
    mean_Cu_coef = mean(Cu_coef, na.rm = TRUE),
    sd_Cu_coef = sd(Cu_coef, na.rm = TRUE))

data_Zn %>%
  group_by(Season) %>%
  summarise(
    mean_Zn_coef = mean(Zn_coef, na.rm = TRUE),
    sd_Zn_coef = sd(Zn_coef, na.rm = TRUE))

data_Ba %>%
  group_by(Season) %>%
  summarise(
    mean_Ba_coef = mean(Ba_coef, na.rm = TRUE),
    sd_Ba_coef = sd(Ba_coef, na.rm = TRUE))

########
# Kruskal–Wallis’ test and post hoc  Dunn’s test
########


# Mg
#Assumptions
model <- aov(Mg_coef ~ Season, data = data_Mg)
leveneTest(Mg_coef ~ Season, data = data_Mg) #p-value < 2.2e-16 ***, we cannot assume the homogeneity of variance
shapiro.test(residuals(model)) #p-value < 2.2e-16, we cannot assume the normality

#Non-parametric test
kruskal.test(Mg_coef ~ Season, data = data_Mg) #Kruskal-Wallis chi-squared = 277.8, df = 3, p-value < 2.2e-16
dunnTest(Mg_coef ~ Season, data = data_Mg,method="bonferroni") #post hoc


# Mn
#Assumptions
model <- aov(Mn_coef ~ Season, data = data_Mn)
leveneTest(Mn_coef ~ Season, data = data_Mn) #p-value = 1.973e-05 ***, we cannot assume the homogeneity of variance
shapiro.test(residuals(model)) #p-value = < 2.2e-16, we cannot assume the normality

#Non-parametric test
kruskal.test(Mn_coef ~ Season, data = data_Mn) #Kruskal-Wallis chi-squared = 26.892, df = 3, p-value = 6.204e-06
dunnTest(Mn_coef ~ Season, data = data_Mn,method="bonferroni") #post hoc


# Cu
#Assumptions
model <- aov(Cu_coef ~ Season, data = data_Cu)
leveneTest(Cu_coef ~ Season, data = data_Cu) #p-value = 1.259e-09 ***, we cannot assume the homogeneity of variance
shapiro.test(residuals(model)) #p-value < 2.2e-16, we cannot assume the normality

#Non-parametric test
kruskal.test(Cu_coef ~ Season, data = data_Cu) #Kruskal-Wallis chi-squared = 45.336, df = 3, p-value = 7.848e-10
dunnTest(Cu_coef ~ Season, data = data_Cu,method="bonferroni") #post hoc


# Zn
#Assumptions
model <- aov(Zn_coef ~ Season, data = data_Zn)
leveneTest(Zn_coef ~ Season, data = data_Zn) #p-value  = 1.277e-07 ***, we cannot assume the homogeneity of variance
shapiro.test(residuals(model)) #p-value < 2.2e-16, we cannot assume the normality

#Non-parametric test
kruskal.test(Zn_coef ~ Season, data = data_Zn) #Kruskal-Wallis chi-squared = 39.787, df = 3, p-value = 1.182e-08
dunnTest(Zn_coef ~ Season, data = data_Zn,method="bonferroni") #post hoc


# Ba
#Assumptions
model <- aov(Ba_coef ~ Season, data = data)
leveneTest(Ba_coef ~ Season, data = data) #p-value = 7.658e-08 ***, we cannot assume the homogeneity of variance
shapiro.test(residuals(model)) #p-value < 2.2e-16, we cannot assume the normality

#Non-parametric test
kruskal.test(Ba_coef ~ Season, data = data) #Kruskal-Wallis chi-squared = 63.896, df = 3, p-value = 8.638e-14
dunnTest(Ba_coef ~ Season, data = data,method="bonferroni") #post hoc


###
# Figure 7 - Partition coefficients
###

###
# Mg
###

##transform x10^-4
data_Mg$Mg_coef_4 <- data_Mg$Mg_coef*10000

#letters based on post hoc Dunn's test
dunn_Mg <- dunnTest(Mg_coef ~ Season, data = data_Mg, method = "bonferroni")$res
dunn_Mg
cld_Mg <- cldList(P.adj ~ Comparison, data=dunn_Mg)
cld_Mg
cld_Mg <- as.data.frame(cld_Mg)                    
cld_Mg <- cld_Mg %>%
  rename(Season = Group)
cld_Mg

letter_Mg <- group_by(data_Mg, Season) %>%
  summarise(mean=mean(Mg_coef_4), quant = quantile(Mg_coef_4, probs = 0.75)) %>%
  arrange(desc(mean))
letter_Mg <- as.data.frame(letter_Mg)                    
letter_Mg

letter_Mg <- merge(letter_Mg, cld_Mg, by = "Season")
letter_Mg

Mg<-ggplot(data_Mg, aes(y=Mg_coef_4, x=Season)) +
  geom_boxplot(alpha = 0.5,width=0.3, fill = "#009FC0") +
  theme_bw() +
  xlab("") + ylab(expression(italic(D)[Mg:Ca]~(x~10^-4)))+
  theme(legend.position="none",
        axis.text = element_text(size = 15),
        axis.text.x = element_text(angle = -45),
        axis.title = element_text(size = 15))+
  scale_y_continuous(limits = c(0.03, 0.27), breaks=seq(0.03, 0.27, by = 0.06))+
  geom_text(data = letter_Mg, aes(x = Season, y = quant, label = Letter), vjust=-2, hjust=-1, size = 4, color='black')
Mg

###
# Mn
###

#letters based on post hoc Dunn's test
dunn_Mn <- dunnTest(Mn_coef ~ Season, data = data_Mn, method = "bonferroni")$res
dunn_Mn
cld_Mn <- cldList(P.adj ~ Comparison, data=dunn_Mn)
cld_Mn
cld_Mn <- as.data.frame(cld_Mn)                    
cld_Mn <- cld_Mn %>%
  rename(Season = Group)
cld_Mn

letter_Mn <- group_by(data_Mn, Season) %>%
  summarise(mean=mean(Mn_coef), quant = quantile(Mn_coef, probs = 0.75)) %>%
  arrange(desc(mean))
letter_Mn <- as.data.frame(letter_Mn)                    
letter_Mn

letter_Mn <- merge(letter_Mn, cld_Mn, by = "Season")
letter_Mn

Mn<-ggplot(data_Mn, aes(y=Mn_coef, x=Season)) +
  geom_boxplot(alpha = 0.5,width=0.3, fill = "#009FC0") +
  theme_bw() +
  xlab("") + ylab(expression(italic(D)[Mn:Ca]))+
  theme(legend.position="none",
        axis.text = element_text(size = 15),
        axis.text.x = element_text(angle = -45),
        axis.title = element_text(size = 15))+
  geom_text(data = letter_Mn, aes(x = Season, y = quant, label = Letter), vjust=-2, hjust=-1, size = 4, color='black')
Mn

###
# Cu
###

#letters based on post hoc Dunn's test
dunn_Cu <- dunnTest(Cu_coef ~ Season, data = data_Cu, method = "bonferroni")$res
dunn_Cu
cld_Cu <- cldList(P.adj ~ Comparison, data=dunn_Cu)
cld_Cu
cld_Cu <- as.data.frame(cld_Cu)                    
cld_Cu <- cld_Cu %>%
  rename(Season = Group)
cld_Cu

letter_Cu <- group_by(data_Cu, Season) %>%
  summarise(mean=mean(Cu_coef), quant = quantile(Cu_coef, probs = 0.75)) %>%
  arrange(desc(mean))
letter_Cu <- as.data.frame(letter_Cu)                    
letter_Cu

letter_Cu <- merge(letter_Cu, cld_Cu, by = "Season")
letter_Cu

Cu<-ggplot(data_Cu, aes(y=Cu_coef, x=Season)) +
  geom_boxplot(alpha = 0.5,width=0.3, fill = "#009FC0") +
  theme_bw() +
  xlab("") + ylab(expression(italic(D)[Cu:Ca]))+
  theme(legend.position="none",
        axis.text = element_text(size = 15),
        axis.text.x = element_text(angle = -45),
        axis.title = element_text(size = 15))+
  geom_text(data = letter_Cu, aes(x = Season, y = quant, label = Letter), vjust=-2, hjust=-1, size = 4, color='black')
Cu

###
# Zn
###

#delete outlier
data_Zn <- data_Zn %>%  filter(!(ID == 'ID324' & Time == '74.7'))

#letters based on post hoc Dunn's test
dunn_Zn <- dunnTest(Zn_coef ~ Season, data = data_Zn, method = "bonferroni")$res
dunn_Zn
cld_Zn <- cldList(P.adj ~ Comparison, data=dunn_Zn)
cld_Zn
cld_Zn <- as.data.frame(cld_Zn)                    
cld_Zn <- cld_Zn %>%
  rename(Season = Group)
cld_Zn

letter_Zn <- group_by(data_Zn, Season) %>%
  summarise(mean=mean(Zn_coef), quant = quantile(Zn_coef, probs = 0.75)) %>%
  arrange(desc(mean))
letter_Zn <- as.data.frame(letter_Zn)                    
letter_Zn

letter_Zn <- merge(letter_Zn, cld_Zn, by = "Season")
letter_Zn

Zn<-ggplot(data_Zn, aes(y=Zn_coef, x=Season)) +
  geom_boxplot(alpha = 0.5,width=0.3, fill = "#009FC0") +
  theme_bw() +
  xlab("") + ylab(expression(italic(D)[Zn:Ca]))+
  theme(legend.position="none",
        axis.text = element_text(size = 15),
        axis.text.x = element_text(angle = -45),
        axis.title = element_text(size = 15))+
  geom_text(data = letter_Zn, aes(x = Season, y = quant, label = Letter), vjust=-2, hjust=-1, size = 4, color='black')
Zn

###
# Ba
###

data_Ba<-data

#letters based on post hoc Dunn's test
dunn_Ba <- dunnTest(Ba_coef ~ Season, data = data_Ba, method = "bonferroni")$res
dunn_Ba
cld_Ba <- cldList(P.adj ~ Comparison, data=dunn_Ba)
cld_Ba
cld_Ba <- as.data.frame(cld_Ba)                    
cld_Ba <- cld_Ba %>%
  rename(Season = Group)
cld_Ba

letter_Ba <- group_by(data_Ba, Season) %>%
  summarise(mean=mean(Ba_coef), quant = quantile(Ba_coef, probs = 0.75)) %>%
  arrange(desc(mean))
letter_Ba <- as.data.frame(letter_Ba)                    
letter_Ba

letter_Ba <- merge(letter_Ba, cld_Ba, by = "Season")
letter_Ba

Ba<-ggplot(data_Ba, aes(y=Ba_coef, x=Season)) +
  geom_boxplot(alpha = 0.5,width=0.3, fill = "#009FC0") +
  theme_bw() +
  xlab("") + ylab(expression(italic(D)[Ba:Ca]))+
  theme(legend.position="none",
        axis.text = element_text(size = 15),
        axis.text.x = element_text(angle = -45),
        axis.title = element_text(size = 15))+
  geom_text(data = letter_Ba, aes(x = Season, y = quant, label = Letter), vjust=-2, hjust=-1, size = 4, color='black')
Ba

######
#Fig. 7. Partition coefficients of elemental ratios (DEl:Ca) at the São Mateus estuary – Southwestern Atlantic. 
#The DEl:Ca was calculated considering El:Cawater with a one–month correction. Acronyms indicate the distinct seasons: 
#LD (late dry 2022), EW (early wet 2022), LW (late wet 2023), and ED (early dry 2023). Lettering indicates post hoc Dunn’s 
#test results in which months are significantly different (p < 0.05), considering DEl:Ca as the response variable, and the 
#sampled month as the predictor. The DEl:Ca was calculated considering El:Cawater with a one–month correction.
######

ggarrange(Mg,Mn,Cu,Zn,Ba, 
          labels = c("A","B","C","D","E"),
          ncol = 2, nrow = 3,
          legend = "none")
#700 x 700
















#################################################
#3. GLMMs
#################################################

#Packages
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

####
# Read dataset
####

#Clean R environment 
rm(list = ls())

#Set working directory
setwd("C:/Users/patri/OneDrive/Documentos/UFES/Tese/From water to biominerals/R/github")

#Read data water
data <- read.csv2("data_m1_e3.csv")
data <- as.data.frame(data)
str(data)

length(unique(data$ID))
table(unique(data$Month))

data <- mutate_at(data, vars(c('Mg24', 'Mn55', 'Cu65','Zn66', 'Ba138', 
                               'TL', 'Fulton',
                               'MgCa','MnCa',	'CuCa',	'ZnCa',	'BaCa',
                               'Temp', 'Sal','OrderS')), as.numeric)

data <- data[data$Month != c('JUL'), ]

data$Month <- factor(data$Month, levels = c('AUG','SEP','OCT','NOV','DEC','JAN','FEB','MAR','APR','MAY','JUN'))
data$Sliding <- factor(data$Sliding, levels = c('JUL','AUG','SEP','OCT','NOV','DEC','JAN','FEB','MAR','APR','MAY'))

length(unique(data$ID))
str(data)

################################
# 3.1. Correlation between variables
################################

str(data)
matriz <- data[, c("Temp", "Sal", "TL", "Fulton")]
str(matriz)

library(corrplot)
par(mfrow = c(1, 1))
corrplot(cor(matriz[], method ='pearson'), type = "lower",
         method = c("number"), tl.col = "black", number.digits = 2)

library(GGally)
require(datasets)
data("swiss")

ggpairs(matriz)

#Fig. S2. Pairplots produced for identifying highly correlated explanatory variables during the exploratory data analysis. 

ggpairs(matriz, columns = 1:4, ggplot2::aes(alpha = 0.3),
        lower = list(continuous = "smooth", alpha = 0.3)) +   theme_bw()


################################
# GLMMs
################################

####
# Read dataset
####

#Clean R environment 
rm(list = ls())

#Read data water
data <- read.csv2("data_m1_e3.csv")

data <- data[data$Month != c('JUL'), ]

data <- mutate_at(data, vars(c('Mg24', 'Mn55', 'Cu65','Zn66', 'Ba138', 
                               'TL', 'Fulton',
                               'MgCa','MnCa',	'CuCa',	'ZnCa',	'BaCa',
                               'Temp', 'Sal','OrderS')), as.numeric)

data$Month <- factor(data$Month, levels = c('AUG','SEP','OCT','NOV','DEC','JAN','FEB','MAR','APR','MAY','JUN'))
data$Sliding <- factor(data$Sliding, levels = c('JUL','AUG','SEP','OCT','NOV','DEC','JAN','FEB','MAR','APR','MAY'))

length(unique(data$ID))





################################
# 3.2. Mg
################################

#Outliers
#Mg24
data_Mg <- data
data_Mg <- data_Mg %>%  filter(!(ID == 'ID335' & Time == '70.2'))
data_Mg <- data_Mg %>%  filter(!(ID == 'ID204' & Time == '88.3'))
data_Mg <- data_Mg %>%  filter(!(ID == 'ID052' & Time == '54.3'))
data_Mg <- data_Mg %>%  filter(!(ID == 'ID509' &  Time == '83.8'))
data_Mg <- data_Mg %>%  filter(!(ID == 'ID496' &  Time == '101.9'))
data_Mg <- data_Mg %>%  filter(!(ID == 'ID383' &  Time == '92.8'))
data_Mg <- data_Mg %>%  filter(!(ID == 'ID052' &  Time == '52.1'))
data_Mg <- data_Mg %>%  filter(!(ID == 'ID274' &  Time == '74.7'))
data_Mg <- data_Mg %>%  filter(!(ID == 'ID511' &  Time == '83.8'))
data_Mg <- data_Mg %>%  filter(!(ID == 'ID487' & Time == '58.9'))

#BoxCox Transformation
box = lm(Mg24 ~ Sal + Temp + TL + Fulton + MgCa, data = data_Mg)
bc<-boxcox(box, plotit = TRUE)
lambda <- bc$x[which.max(bc$y)]

if (lambda == 0) {
  data_Mg$Mg24t <- log(data_Mg$Mg24)
} else {
  data_Mg$Mg24t <- (data_Mg$Mg24^lambda - 1) / lambda
}

summary(data_Mg$Mg24t)

#Scale the variables
data_Mg$Mg24ts <- rescale(data_Mg$Mg24t)
data_Mg$Mg24s <- rescale(data_Mg$Mg24)
data_Mg$Sals <- rescale(data_Mg$Sal)
data_Mg$Temps <- rescale(data_Mg$Temp)
data_Mg$TLs <- rescale(data_Mg$TL)
data_Mg$Fultons <- rescale(data_Mg$Fulton)
data_Mg$MgCas <- rescale(data_Mg$MgCa)


####
## Mg - Dredge
####

#Gaussian
full_model_Mg <- lmer(Mg24ts ~ Sals + Temps + TLs + Fultons + MgCas + (1 | Month) + (1 | ID), 
                       data = data_Mg, REML = FALSE,
                       na.action =  na.fail)

dredge_Mg <- dredge(full_model_Mg)
dredge_Mg
subset(dredge_Mg, delta==0)
dredge_Mg[1:5,]

#Model selection table 
#(Intrc)     Fltns    MgCas   Sals  Temps      TLs  df  logLik   AICc delta weight
#10  0.4992 -0.2612                 0.1152           6 176.991 -341.7  0.00  0.327
#14  0.4685 -0.2706          0.0707 0.1162           7 177.857 -341.4  0.36  0.274
#2   0.5405 -0.2161                                  5 175.266 -340.3  1.37  0.165
#26  0.5080 -0.2626                 0.1191 -0.01614  7 177.012 -339.7  2.05  0.118
1#2  0.4948 -0.2623 0.008683        0.1155           7 177.002 -339.6  2.07  0.116

####
## Mg - Selected model
####

optimal_model_Mg24 <- lmer(Mg24ts ~ Fultons + Temps + Sals + (1 | Month) + (1 | ID), 
                   data = data_Mg, REML = FALSE,
                   na.action =  na.fail)


summary(optimal_model_Mg24)
Anova(optimal_model_Mg24)

#Analysis of Deviance Table (Type II Wald chisquare tests)
#Response: Mg24ts
#         Chisq Df Pr(>Chisq)  
#Fultons 7.6739  1   0.005602 **
#Temps   5.3371  1   0.020877 * 
#Sals    1.9000  1   0.168079 


####
#  Mg - Model validation
####
plot(optimal_model_Mg24)

#"DHARMa"
simulationOutput <- simulateResiduals(fittedModel = optimal_model_Mg24, plot = TRUE, use.u = FALSE)
testDispersion(optimal_model_Mg24, type='DHARMa')
#dispersion = 0.99789, p-value = 0.984

#r2 values
r.squaredGLMM(optimal_model_Mg24)
r2_nakagawa(optimal_model_Mg24)
#Conditional R2: 0.666
#Marginal R2: 0.093

#check_distribution
check_distribution <- check_distribution(optimal_model_Mg24)
plot(check_distribution)
check_distribution(optimal_model_Mg24)

#Variance inflation factors (VIF)
check_collinearity(optimal_model_Mg24)
#Low Correlation

#Temporal autocorrelation
res = simulateResiduals(optimal_model_Mg24)
res = recalculateResiduals(res, group = data_Mg$OrderS)
testTemporalAutocorrelation(res, time = unique(data_Mg$OrderS))
#Durbin-Watson test = 1.7373, p-value = 0.6504







################################
#3.3. Mn
################################

#Outliers
#Mn55
data_Mn <- data
data_Mn <- data_Mn %>%  filter(!(ID == 'ID165'))
data_Mn <- data_Mn %>%  filter(!(ID == 'ID382' & Time == '81.5'))
data_Mn <- data_Mn %>%  filter(!(ID == 'ID327' & Time == '70.2'))

#BoxCox Transformation
box = lm(Mn55 ~ Sal + Temp + TL + Fulton + MgCa, data = data_Mn)
bc<-boxcox(box, plotit = TRUE)
lambda <- bc$x[which.max(bc$y)]

if (lambda == 0) {
  data_Mn$Mn55t <- log(data_Mn$Mn55)
} else {
  data_Mn$Mn55t <- (data_Mn$Mn55^lambda - 1) / lambda
}

summary(data_Mn$Mn55t)

#Scale the variables
data_Mn$Mn55ts <- rescale(data_Mn$Mn55t)
data_Mn$Mn55s <- rescale(data_Mn$Mn55)
data_Mn$Sals <- rescale(data_Mn$Sal)
data_Mn$Temps <- rescale(data_Mn$Temp)
data_Mn$TLs <- rescale(data_Mn$TL)
data_Mn$Fultons <- rescale(data_Mn$Fulton)
data_Mn$MnCas <- rescale(data_Mn$MnCa)


####
## Mn - Dredge
####

#Gaussian
full_model_Mn <- lmer(Mn55ts ~ Sals + Temps + TLs + Fultons + MnCas + (1 |Month) + (1 | ID), 
                       data = data_Mn, REML = FALSE,
                       na.action =  na.fail)


dredge_Mn <- dredge(full_model_Mn)
dredge_Mn
subset(dredge_Mn, delta==0)
dredge_Mn[1:5,]

#Model selection table 
#(Intrc)   Fltns   MnCas     Sals  Temps      TLs     df logLik   AICc delta weight
#9   0.4252                          0.1350           5 74.254 -138.3  0.00  0.336
#25  0.4724                          0.1559 -0.09194  6 74.843 -137.4  0.90  0.215
#13  0.4439                 -0.04013 0.1357           6 74.601 -136.9  1.38  0.168
#10  0.4026 0.08186                  0.1252           6 74.545 -136.8  1.49  0.159
#11  0.4217         0.01137          0.1365           6 74.274 -136.3  2.04  0.121

####
## Mn - Selected model
####

optimal_model_Mn55 <- lmer(Mn55ts ~ Temps + Sals + (1 | Month) + (1 | ID), 
                          data = data_Mn, REML = FALSE,
                          na.action =  na.fail)

summary(optimal_model_Mn55)
Anova(optimal_model_Mn55)

#Analysis of Deviance Table (Type II Wald chisquare tests)
#Response: Mn55ts
#       Chisq Df Pr(>Chisq)   
#Temps 8.7886  1   0.003031 **
#Sals  0.6967  1   0.403895

####
#  Mn - Model validation
####
plot(optimal_model_Mn55)

#"DHARMa"
simulationOutput <- simulateResiduals(fittedModel = optimal_model_Mn55, plot = TRUE, use.u = FALSE)
testDispersion(optimal_model_Mn55, type='DHARMa')
#dispersion = 1.0073, p-value = 0.912

#r2 values
r.squaredGLMM(optimal_model_Mn55)
r2_nakagawa(optimal_model_Mn55)
# R2m       R2c
#0.0499088  0.44056281

#check_distribution
check_distribution <- check_distribution(optimal_model_Mn55)
plot(check_distribution)
check_distribution(optimal_model_Mn55)

#Variance inflation factors (VIF)
check_collinearity(optimal_model_Mn55)
#Low Correlation

#Temporal autocorrelation
res = simulateResiduals(optimal_model_Mn55)
res = recalculateResiduals(res, group = data_Mn$OrderS)
testTemporalAutocorrelation(res, time = unique(data_Mn$OrderS))
#Durbin-Watson test = 1.252, p-value = 0.1817






################################
#3.4. Cu
################################

#Outliers
#Cu65
data_Cu <- data
data_Cu <- data_Cu %>%  filter(!(ID == 'ID083' & Time == '72.4'))
data_Cu <- data_Cu %>%  filter(!(ID == 'ID072' & Time == '81.5'))
data_Cu <- data_Cu %>%  filter(!(ID == 'ID076' & Time == '77'))
data_Cu <- data_Cu %>%  filter(!(ID == 'ID202' & Time == '70.2'))
data_Cu <- data_Cu %>%  filter(!(ID == 'ID202' & Time == '72.4'))
data_Cu <- data_Cu %>%  filter(!(ID == 'ID277' & Time == '72.4'))

#BoxCox Transformation
box = lm(Cu65 ~ Sal + Temp + TL + Fulton + CuCa, data = data_Cu)
bc<-boxcox(box, plotit = TRUE)
lambda <- bc$x[which.max(bc$y)]

if (lambda == 0) {
  data_Cu$Cu65t <- log(data_Cu$Cu65)
} else {
  data_Cu$Cu65t <- (data_Cu$Cu65^lambda - 1) / lambda
}

summary(data_Cu$Cu65t)
data_Cu$Cu65ts <- scale(data_Cu$Cu65t)

#Scale the predictors
data_Cu$Cu65ts <- rescale(data_Cu$Cu65t)
data_Cu$Cu65s <- rescale(data_Cu$Cu65)
data_Cu$Sals <- rescale(data_Cu$Sal)
data_Cu$Temps <- rescale(data_Cu$Temp)
data_Cu$TLs <- rescale(data_Cu$TL)
data_Cu$Fultons <- rescale(data_Cu$Fulton)
data_Cu$CuCas <- rescale(data_Cu$CuCa)


####
## Cu - Dredge
####

#Gaussian
full_model_Cu <- lmer(Cu65ts ~ Sals + Temps + TLs + Fultons + CuCas + (1 |Month) + (1 | ID), 
                       data = data_Cu, REML = FALSE,
                       na.action =  na.fail)


dredge_Cu <- dredge(full_model_Cu)
dredge_Cu
subset(dredge_Cu, delta==0)
dredge_Cu[1:5,]

#Model selection table 
#(Intrc)    CuCas    Sals   Temps   df  logLik   AICc delta weight
#1   0.4938                          4 114.271 -220.4  0.00  0.356
#9   0.4757                 0.03631  5 114.685 -219.2  1.23  0.192
#2   0.4818 0.03398                  5 114.533 -218.9  1.54  0.165
#10  0.4445 0.06048         0.05588  6 115.439 -218.6  1.80  0.145
#5   0.4840         0.01999          5 114.381 -218.6  1.84  0.142

####
## Cu - Selected model
####

optimal_model_Cu65 <- lmer(Cu65ts ~ CuCas + Temps + (1 |Month) + (1 | ID), 
                   data = data_Cu, REML = FALSE,
                   na.action =  na.fail)

summary(optimal_model_Cu65)
Anova(optimal_model_Cu65)

#Analysis of Deviance Table (Type II Wald chisquare tests)
#Response: Cu65ts
#       Chisq Df Pr(>Chisq)
#CuCas 1.6303  1     0.2017
#Temps 1.9280  1     0.1650


####
#  Cu - Model validation
####
plot(optimal_model_Cu65)

#"DHARMa"
simulationOutput <- simulateResiduals(fittedModel = optimal_model_Cu65, plot = TRUE, use.u = FALSE)
testDispersion(optimal_model_Cu65, type='DHARMa')
#dispersion = 0.99853, p-value = 0.944

#r2 values
r.squaredGLMM(optimal_model_Cu65)
r2_nakagawa(optimal_model_Cu65)
#Conditional R2: 0.211
#Marginal R2: 0.013

#check_distribution
check_distribution <- check_distribution(optimal_model_Cu65)
plot(check_distribution)
check_distribution(optimal_model_Cu65)

#Variance inflation factors (VIF)
check_collinearity(optimal_model_Cu65)
#Low Correlation

#Temporal autocorrelation
res = simulateResiduals(optimal_model_Cu65)
res = recalculateResiduals(res, group = data_Cu$OrderS)
testTemporalAutocorrelation(res, time = unique(data_Cu$OrderS))
#Durbin-Watson test = 1.5221, p-value = 0.4046






################################
#3.5. Zn
################################

#Outliers
#Zn66
data_Zn <- data
data_Zn <- data_Zn %>%  filter(!(ID == 'ID386' & Time == '72.4'))
data_Zn <- data_Zn %>%  filter(!(ID == 'ID079' & Time == '72.4'))
data_Zn <- data_Zn %>%  filter(!(ID == 'ID509' & Time == '83.8'))
data_Zn <- data_Zn %>%  filter(!(ID == 'ID511' & Time == '83.8'))
data_Zn <- data_Zn %>%  filter(!(ID == 'ID508' & Time == '47.5'))

#BoxCox Transformation
box = lm(Zn66 ~ Sal + Temp + TL + Fulton + MgCa, data = data_Zn)
bc<-boxcox(box, plotit = TRUE)
lambda <- bc$x[which.max(bc$y)]

if (lambda == 0) {
  data_Zn$Zn66t <- log(data_Zn$Zn66)
} else {
  data_Zn$Zn66t <- (data_Zn$Zn66^lambda - 1) / lambda
}

#Scale the predictors
data_Zn$Zn66ts <- rescale(data_Zn$Zn66t)
data_Zn$Zn66s <- rescale(data_Zn$Zn66)
data_Zn$Sals <- rescale(data_Zn$Sal)
data_Zn$Temps <- rescale(data_Zn$Temp)
data_Zn$TLs <- rescale(data_Zn$TL)
data_Zn$Fultons <- rescale(data_Zn$Fulton)
data_Zn$ZnCas <- rescale(data_Zn$ZnCa)


####
## Zn - Dredge
####

#Gaussian
full_model_Zn <- lmer(Zn66ts ~ Sals + Temps + TLs + Fultons + ZnCas + (1 | Month) + (1 | ID), 
                       data = data_Zn, REML = FALSE,
                       na.action =  na.fail)


dredge_Zn <- dredge(full_model_Zn)
dredge_Zn
subset(dredge_Zn, delta==0)
dredge_Zn[1:5,]

#Model selection table 
#(Intrc)    Fltns     Sals   Temps    ZnCas   df  logLik   AICc delta weight
#1   0.6382                                    4 182.607 -357.1  0.00  0.343
#3   0.6559          -0.03653                  5 183.109 -356.0  1.06  0.202
#2   0.6553 -0.05139                           5 182.915 -355.6  1.45  0.166
#17  0.6308                           0.02663  5 182.776 -355.4  1.73  0.145
#5   0.6284                   0.02011          5 182.768 -355.3  1.74  0.144

####
## Zn - Selected model
####

optimal_model_Zn66 <- lmer(Zn66ts ~ ZnCas + (1 | Month) + (1 | ID), 
                   data = data_Zn, REML = FALSE,
                   na.action =  na.fail)


summary(optimal_model_Zn66)
Anova(optimal_model_Zn66)

#Analysis of Deviance Table (Type II Wald chisquare tests)
#Response: Mg24ts
#         Chisq Df Pr(>Chisq)  
#ZnCas 0.3435  1     0.5578

####
#  Zn - Model validation
####
plot(optimal_model_Zn66)

#"DHARMa"
simulationOutput <- simulateResiduals(fittedModel = optimal_model_Zn66, plot = TRUE, use.u = FALSE)
testDispersion(optimal_model_Zn66, type='DHARMa')
#dispersion = 1.0007, p-value = 0.96

#r2 values
r.squaredGLMM(optimal_model_Zn66)
r2_nakagawa(optimal_model_Zn66)
#Conditional R2: 0.246
#Marginal R2: 0.002

#check_distribution
check_distribution <- check_distribution(optimal_model_Zn66)
plot(check_distribution)
check_distribution(optimal_model_Zn66)

#Variance inflation factors (VIF)
check_collinearity(optimal_model_Zn66)
#Not enough model terms in the conditional part of the model to check for multicollinearity

#Temporal autocorrelation
res = simulateResiduals(optimal_model_Zn66)
res = recalculateResiduals(res, group = data_Zn$OrderS)
testTemporalAutocorrelation(res, time = unique(data_Zn$OrderS))
#Durbin-Watson test = 0.6675, p-value = 0.008102






################################
#3.6. Ba
################################

#No outliers

#BoxCox Transformation
box = lm(Ba138 ~ Sal + Temp + TL + Fulton + BaCa, data = data)

bc<-boxcox(box, plotit = TRUE)
lambda <- bc$x[which.max(bc$y)]

if (lambda == 0) {
  data$Ba138t <- log(data$Ba138)
} else {
  data$Ba138t <- (data$Ba138^lambda - 1) / lambda
}

#Scale the variables
data$Ba138ts <- rescale(data$Ba138t)
data$Ba138s <- rescale(data$Ba138)
data$Sals <- rescale(data$Sal)
data$Temps <- rescale(data$Temp)
data$TLs <- rescale(data$TL)
data$Fultons <- rescale(data$Fulton)
data$BaCas <- rescale(data$BaCa)


####
## Ba - Dredge
####

#Gaussian
full_model_Ba <- lmer(Ba138ts ~ Sals + Temps + TLs + Fultons + BaCas + (1 | Month) + (1 | ID), 
                       data = data, REML = FALSE,
                       na.action =  na.fail)

dredge_Ba <- dredge(full_model_Ba)
dredge_Ba
subset(dredge_Ba, delta==0)
dredge_Ba[1:5,]

#Model selection table 

#Model selection table 
# (Intrc)   BaCas Fltns     Sals  Temps    TLs  df  logLik   AICc delta weight
#25  0.3122                       0.2367 0.1698  6 334.436 -656.6  0.00  0.348
#26  0.2864 0.1115                0.2309 0.1697  7 335.026 -655.7  0.91  0.221
#17  0.4212                              0.1832  5 332.620 -655.1  1.56  0.160
#29  0.3394              -0.05556 0.2364 0.1697  7 334.553 -654.8  1.85  0.138
#27  0.3025        0.033          0.2322 0.1716  7 334.519 -654.7  1.92  0.133

####
## Ba - Selected model
####

optimal_model_Ba138 <- lmer(Ba138ts ~ BaCas + Temps + TLs + (1 | Month) + (1 | ID), 
                            data = data, REML = FALSE,
                            na.action =  na.fail)

summary(optimal_model_Ba138)
Anova(optimal_model_Ba138)

#Analysis of Deviance Table (Type II Wald chisquare tests)
#Response: Mg24ts
#         Chisq Df Pr(>Chisq)  
#BaCas 1.2723  1   0.259335   
#Temps 4.7840  1   0.028725 * 
# TLs   7.5194  1   0.006104 **


####
#  Ba - Model validation
####
plot(optimal_model_Ba138)

#"DHARMa"
simulationOutput <- simulateResiduals(fittedModel = optimal_model_Ba138, plot = TRUE, use.u = FALSE)
testDispersion(optimal_model_Ba138, type='DHARMa')
#dispersion = 0.9094, p-value = 0.832

#r2 values
r.squaredGLMM(optimal_model_Ba138)
r2_nakagawa(optimal_model_Ba138)
#Conditional R2: 0.916
#Marginal R2: 0.299

#check_distribution
check_distribution <- check_distribution(optimal_model_Ba138)
plot(check_distribution)
check_distribution(optimal_model_Ba138)

#Variance inflation factors (VIF)
check_collinearity(optimal_model_Ba138)
#Low Correlation

#Temporal autocorrelation
res = simulateResiduals(optimal_model_Ba138)
res = recalculateResiduals(res, group = data$OrderS)
testTemporalAutocorrelation(res, time = unique(data$OrderS))
#Durbin-Watson test = 1.3969, p-value = 0.2884











################################
# 3.7 Figure 8 - Significant relationships based on our GLMMs
################################

####
# Read dataset
####

#Clean R environment 
rm(list = ls())

#Read data water
data <- read.csv2("data_m1_e3.csv")

data <- data[data$Month != c('JUL'), ]

data <- mutate_at(data, vars(c('Mg24', 'Mn55', 'Cu65','Zn66', 'Ba138', 
                               'TL', 'Fulton',
                               'MgCa','MnCa',	'CuCa',	'ZnCa',	'BaCa',
                               'Temp', 'Sal','OrderS')), as.numeric)

####
# Mg - Figure
####

# Delete outliers
# Mg24
data_Mg <- data
data_Mg <- data_Mg %>%  filter(!(ID == 'ID335' & Time == '70.2'))
data_Mg <- data_Mg %>%  filter(!(ID == 'ID204' & Time == '88.3'))
data_Mg <- data_Mg %>%  filter(!(ID == 'ID052' & Time == '54.3'))
data_Mg <- data_Mg %>%  filter(!(ID == 'ID509' &  Time == '83.8'))
data_Mg <- data_Mg %>%  filter(!(ID == 'ID496' &  Time == '101.9'))
data_Mg <- data_Mg %>%  filter(!(ID == 'ID383' &  Time == '92.8'))
data_Mg <- data_Mg %>%  filter(!(ID == 'ID052' &  Time == '52.1'))
data_Mg <- data_Mg %>%  filter(!(ID == 'ID274' &  Time == '74.7'))
data_Mg <- data_Mg %>%  filter(!(ID == 'ID511' &  Time == '83.8'))
data_Mg <- data_Mg %>%  filter(!(ID == 'ID487' & Time == '58.9'))

optimal_model_Mg24 <- lmer(Mg24 ~ Fulton + Temp + Sal + (1 | Month) + (1 | ID), 
                           data = data_Mg, REML = FALSE,
                           na.action =  na.fail)

Anova(optimal_model_Mg24)

#Profiled Confidence Intervals
profile_ci_Mg_Fulton <- ggpredict(optimal_model_Mg24, terms = c("Fulton"))

Mg_Fulton <-plot(profile_ci_Mg_Fulton) +
  labs(title = " ",
       x = expression(Fulton~index), y = expression(Mg:Ca[otolith]~(mu~mol~mol^-1)))+
  theme_bw()+
  theme(legend.position="none",
        title = element_text(size = 15),
        axis.text = element_text(size = 15),
        axis.title = element_text(size = 15))
Mg_Fulton

profile_ci_Mg_Temp <- ggpredict(optimal_model_Mg24, terms = c("Temp"))

Mg_Temp <-plot(profile_ci_Mg_Temp) +
  labs(title = " ",
       x = expression(Temperature~(ºC)), y = expression(Mg:Ca[otolith]~(mu~mol~mol^-1)))+
  theme_bw()+
  theme(legend.position="none",
        title = element_text(size = 15),
        axis.text = element_text(size = 15),
        axis.title = element_text(size = 15))
Mg_Temp

####
# Mn - Figure
####

# Delete outliers
# Mn55

data_Mn <- data
data_Mn <- data_Mn %>%  filter(!(ID == 'ID165'))
data_Mn <- data_Mn %>%  filter(!(ID == 'ID382' & Time == '81.5'))
data_Mn <- data_Mn %>%  filter(!(ID == 'ID327' & Time == '70.2'))

optimal_model_Mn55 <- lmer(Mn55 ~ Temp + Sal + (1 | Month) + (1 | ID), 
                           data = data_Mn, REML = FALSE,
                           na.action =  na.fail)

Anova(optimal_model_Mn55)

#Profiled Confidence Intervals

profile_ci_Mn_Temp <- ggpredict(optimal_model_Mn55, terms = c("Temp"))

Mn_Temp<-plot(profile_ci_Mn_Temp) +
  labs(title = " ",
       x = expression(Temperature~(ºC)), y = expression(Mn:Ca[otolith]~(mu~mol~mol^-1)), color = "")+
  theme_bw()+
  theme(legend.position="none",
        title = element_text(size = 15),
        axis.text = element_text(size = 15),
        axis.title = element_text(size = 15))

Mn_Temp


####
# Ba - Figure
####

# No outliers

optimal_model_Ba138 <- lmer(Ba138 ~ BaCa + Temp + TL + (1 | Month) + (1 | ID), 
                            data = data, REML = FALSE,
                            na.action =  na.fail)

Anova(optimal_model_Ba138)

#Profiled Confidence Intervals

profile_ci_Ba_Temp <- ggpredict(optimal_model_Ba138, terms = c("Temp"))

Ba_Temp <-plot(profile_ci_Ba_Temp) +
  labs(title = " ",
       x = expression(Temperature~(ºC)), y =  expression(Ba:Ca[otolith]~(mu~mol~mol^-1)), color = "")+
  theme_bw()+
  theme(legend.position="none",
        title = element_text(size = 15),
        axis.text = element_text(size = 15),
        axis.title = element_text(size = 15))
Ba_Temp

profile_ci_Ba_TL <- ggpredict(optimal_model_Ba138, terms = c("TL"))

Ba_TL <-plot(profile_ci_Ba_TL) +
  labs(title = " ",
       x = expression(Total~length~(mm)), y =  expression(Ba:Ca[otolith]~(mu~mol~mol^-1)), color = "")+
  theme_bw()+
  theme(legend.position="none",
        title = element_text(size = 15),
        axis.text = element_text(size = 15),
        axis.title = element_text(size = 15))

Ba_TL

#############
#Fig. 8. The significant relationships between otolith elemental signatures (El:Caotolith) 
#and environmental and biological factors according to our generalized linear mixed models (GLMMs).
#############

library(ggpubr)

Mg_Fulton
Mg_Temp
Mn_Temp
Ba_TL
Ba_TL

ggarrange(Mg_Temp, Mg_Fulton, Mn_Temp, Ba_TL, Ba_Temp,
          labels = c("A","B","C","D","E"),
          ncol = 2, nrow = 3,
          legend = "none")
#800 x 650

#END