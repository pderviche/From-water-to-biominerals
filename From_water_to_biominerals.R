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
              #Figure 2
        #1.2 Figure 3b - Scatter plot of laser ablation ICP-MS (ID333)..................................line 165
        #1.3 Searching for the best temporal window.....................................................line 210
              #Element selection - coefficient of variation 
        #1.4 Save dataset with one month correction.....................................................line 1000
        #1.5 Figure 4 - One month correction............................................................line 1110

  #2. Checking data / Statistical analysis..............................................................line 1200
        #2.1 General results............................................................................line 1250
              #Table 1 - TL and Fulton’s index  of the dog snappers
        #2.2 Dog snappers...............................................................................line 1280
              #Kruskal–Wallis’ nonparametric test, and a post hoc Dunn’s test 
        #2.3 Water chemistry across the saline gradient.................................................line 1340
              #Figure 5 - Water chemistry across the salinity gradient
        #2.4 Otolith elemental signatures and their relationships with water chemistry..................line 1650
          #2.4.1 Mean, SD, and range....................................................................line 1680
              #Kruskal–Wallis’ nonparametric test, and a post hoc Dunn’s test 
          #2.4.2 Figure 6 - Otolith elemental signatures ...............................................line 1830
          #2.4.3 Figure 7 - Correlation between otolith and water chemistries...........................line 2190
        #2.5 Partition coefficients as a proxy for absorption rate......................................line 2430
              #Figure 8 - Partition coefficients 

  #3. GLMMs.............................................................................................line 2800
        #3.1 Read the data e correlation between variables..............................................line 2850
        #3.2. Mg........................................................................................line 2900
        #3.3. Mn........................................................................................line 3050
        #3.4. Cu........................................................................................line 3190
        #3.5. Zn........................................................................................line 3320
        #3.6. Sr........................................................................................line 3440
        #3.7. Ba........................................................................................line 3580
        #3.8. Figure 9 - The significant relationships based on our GLMMs...............................line 3730























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
#'Mg24w', 'Ca43w', 'Mn55w', 'Cu65w', 'Zn66w', 'Sr88w', and 'Ba138w' are in umol
#'MgCa', 'MnCa', 'CuCa', 'ZnCa', 'SrCa', and 'BaCa' are the ratio El:Ca43 in umol mol-1

#Filter data
water <- water %>%
  filter(tide == "low") %>%
  filter(site == "Sao Mateus")

#Set as numeric
water <- mutate_at(water, vars(c('Mg24w', 'Ca43w', 'Mn55w', 'Cu65w', 'Zn66w', 'Sr88w', 'Ba138w',
                                 'MgCa', 'MnCa', 'CuCa', 'ZnCa', 'SrCa', 'BaCa',
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
#Fig. 2. Seawater temperature (in green) and salinity (in blue) during the water samplings from July 2022 
#to June 2023 in the São Mateus estuary, Southwestern Atlantic. Shade and transparent areas indicate the
#distinct seasons: late dry (from July to September), early wet (from October to December), late wet (from 
#January to March), and early dry (from April to June).
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





############
#1.2 Figure 3b
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
#'Mg24w', 'Ca43w', 'Mn55w', 'Cu65w', 'Zn66w', 'Sr88w', and 'Ba138w' are in umol
#'MgCa', 'MnCa', 'CuCa', 'ZnCa', 'SrCa', and 'BaCa' are the ratio El:Ca43 in umol mol-1

#Filter water data
water <- water %>%
  filter(tide == "low") %>%
  filter(site == "Sao Mateus")

#Set as numeric
water <- mutate_at(water, vars(c('Mg24w', 'Ca43w', 'Mn55w', 'Cu65w', 'Zn66w', 'Sr88w', 'Ba138w',
                                 'MgCa', 'MnCa', 'CuCa', 'ZnCa', 'SrCa', 'BaCa',
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
## No correction dataset (Otolith sampled in JUL correspond to the Water from JUL)
####

#Create the dataset with no correction 
otolith_0 <- otolith
otolith_0$Sliding <- otolith_0$Month #Create the sliding
otolith_0$Order <- otolith_0$Month #Create the Order
otolith_0$OrderS <- otolith_0$Month #Create the OrderS

otolith_0 <- otolith_0[,c(1,4,18:20,2:3,5:17)] #Reorder

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

otolith_1 <- otolith_1[,c(1,4,18:20,2:3,5:17)] #Reorder

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

otolith_2 <- otolith_2[,c(1,4,18:20,2:3,5:17)] #Reorder

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

otolith_3 <- otolith_3[,c(1,4,18:20,2:3,5:17)] #Reorder

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

otolith_4 <- otolith_4[,c(1,4,18:20,2:3,5:17)] #Reorder

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

otolith_5 <- otolith_5[,c(1,4,18:20,2:3,5:17)] #Reorder

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

otolith_6 <- otolith_6[,c(1,4,18:20,2:3,5:17)] #Reorder

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
sliding0 <- sliding0[,c(2,3,1,4:36)] #Reorder
length(unique(sliding0$ID))

sliding1 = merge(otolith_1, water, on='Sliding', all=FALSE)
str(sliding1)
sliding1$Month <- factor(sliding1$Month, levels = c('AUG','SEP','OCT','NOV','DEC','JAN','FEB','MAR','APR','MAY','JUN'))
sliding1 <- sliding1[,c(2,3,1,4:36)] #Reorder
length(unique(sliding1$ID))

sliding2 = merge(otolith_2, water, on='Sliding', all=FALSE)
str(sliding2)
sliding2$Month <- factor(sliding2$Month, levels = c('SEP','OCT','NOV','DEC','JAN','FEB','MAR','APR','MAY','JUN'))
sliding2 <- sliding2[,c(2,3,1,4:36)] #Reorder
length(unique(sliding2$ID))

sliding3 = merge(otolith_3, water, on='Sliding', all=FALSE)
str(sliding3)
sliding3$Month <- factor(sliding3$Month, levels = c('OCT','NOV','DEC','JAN','FEB','MAR','APR','MAY','JUN'))
sliding3 <- sliding3[,c(2,3,1,4:36)] #Reorder
length(unique(sliding3$ID))

sliding4 = merge(otolith_4, water, on='Sliding', all=FALSE)
str(sliding4)
sliding4$Month <- factor(sliding4$Month, levels = c('NOV','DEC','JAN','FEB','MAR','APR','MAY','JUN'))
sliding4 <- sliding4[,c(2,3,1,4:36)] #Reorder
length(unique(sliding4$ID))

sliding5 = merge(otolith_5, water, on='Sliding', all=FALSE)
str(sliding5)
sliding5$Month <- factor(sliding5$Month, levels = c('DEC','JAN','FEB','MAR','APR','MAY','JUN'))
sliding5 <- sliding5[,c(2,3,1,4:36)] #Reorder
length(unique(sliding5$ID))

sliding6 = merge(otolith_6, water, on='Sliding', all=FALSE)
str(sliding6)
sliding6$Month <- factor(sliding6$Month, levels = c('JAN','FEB','MAR','APR','MAY','JUN'))
sliding6 <- sliding6[,c(2,3,1,4:36)] #Reorder
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

sliding0 <- mutate_at(sliding0, vars(c('Mg24', 'Mn55', 'Cu65', 'Zn66', 'Sr87', 'Ba138', 'TL', 'Weight', 'Fulton',
                                       'Mg24w', 'Mn55w',	'Cu65w',	'Zn66w', 'Sr88w',	'Ba138w',
                                       'MgCa',	'MnCa',	'CuCa',	'ZnCa',	'SrCa',	'BaCa',	'Temp', 'Sal')), as.numeric)
sliding0$Month <- factor(sliding0$Month, levels = c('JUL','AUG','SEP','OCT','NOV','DEC'))
sliding0$Sliding <- factor(sliding0$Sliding, levels = c('JUL','AUG','SEP','OCT','NOV','DEC'))
str(sliding0)
summary(sliding0$Month)
summary(sliding0$Sliding)

sliding1 <- mutate_at(sliding1, vars(c('Mg24', 'Mn55', 'Cu65', 'Zn66', 'Sr87', 'Ba138', 'TL', 'Weight', 'Fulton',
                                       'Mg24w', 'Mn55w',	'Cu65w',	'Zn66w', 'Sr88w',	'Ba138w',
                                       'MgCa',	'MnCa',	'CuCa',	'ZnCa',	'SrCa',	'BaCa',	'Temp', 'Sal')), as.numeric)
sliding1$Month <- factor(sliding1$Month, levels = c('AUG','SEP','OCT','NOV','DEC','JAN'))
sliding1$Sliding <- factor(sliding1$Sliding, levels = c('JUL','AUG','SEP','OCT','NOV','DEC'))
str(sliding1)
summary(sliding1$Month)
summary(sliding1$Sliding)

sliding2 <- mutate_at(sliding2, vars(c('Mg24', 'Mn55', 'Cu65', 'Zn66', 'Sr87', 'Ba138', 'TL', 'Weight', 'Fulton',
                                       'Mg24w', 'Mn55w',	'Cu65w',	'Zn66w', 'Sr88w',	'Ba138w',
                                       'MgCa',	'MnCa',	'CuCa',	'ZnCa',	'SrCa',	'BaCa',	'Temp', 'Sal')), as.numeric)
sliding2$Month <- factor(sliding2$Month, levels = c('SEP','OCT','NOV','DEC','JAN','FEB'))
sliding2$Sliding <- factor(sliding2$Sliding, levels = c('JUL','AUG','SEP','OCT','NOV','DEC'))
str(sliding2)
summary(sliding2$Month)
summary(sliding2$Sliding)

sliding3 <- mutate_at(sliding3, vars(c('Mg24', 'Mn55', 'Cu65', 'Zn66', 'Sr87', 'Ba138', 'TL', 'Weight', 'Fulton',
                                       'Mg24w', 'Mn55w',	'Cu65w',	'Zn66w', 'Sr88w',	'Ba138w',
                                       'MgCa',	'MnCa',	'CuCa',	'ZnCa',	'SrCa',	'BaCa',	'Temp', 'Sal')), as.numeric)
sliding3$Month <- factor(sliding3$Month, levels = c('OCT','NOV','DEC','JAN','FEB','MAR'))
sliding3$Sliding <- factor(sliding3$Sliding, levels = c('JUL','AUG','SEP','OCT','NOV','DEC'))
str(sliding3)
summary(sliding3$Month)
summary(sliding3$Sliding)

sliding4 <- mutate_at(sliding4, vars(c('Mg24', 'Mn55', 'Cu65', 'Zn66', 'Sr87', 'Ba138', 'TL', 'Weight', 'Fulton',
                                       'Mg24w', 'Mn55w',	'Cu65w',	'Zn66w', 'Sr88w',	'Ba138w',
                                       'MgCa',	'MnCa',	'CuCa',	'ZnCa',	'SrCa',	'BaCa',	'Temp', 'Sal')), as.numeric)
sliding4$Month <- factor(sliding4$Month, levels = c('NOV','DEC','JAN','FEB','MAR','APR'))
sliding4$Sliding <- factor(sliding4$Sliding, levels = c('JUL','AUG','SEP','OCT','NOV','DEC'))
str(sliding4)
summary(sliding4$Month)
summary(sliding4$Sliding)

sliding5 <- mutate_at(sliding5, vars(c('Mg24', 'Mn55', 'Cu65', 'Zn66', 'Sr87', 'Ba138', 'TL', 'Weight', 'Fulton',
                                       'Mg24w', 'Mn55w',	'Cu65w',	'Zn66w', 'Sr88w',	'Ba138w',
                                       'MgCa',	'MnCa',	'CuCa',	'ZnCa',	'SrCa',	'BaCa',	'Temp', 'Sal')), as.numeric)
sliding5$Month <- factor(sliding5$Month, levels = c('DEC','JAN','FEB','MAR','APR','MAY'))
sliding5$Sliding <- factor(sliding5$Sliding, levels = c('JUL','AUG','SEP','OCT','NOV','DEC'))
str(sliding5)
summary(sliding5$Month)
summary(sliding5$Sliding)

sliding6 <- mutate_at(sliding6, vars(c('Mg24', 'Mn55', 'Cu65', 'Zn66', 'Sr87', 'Ba138', 'TL', 'Weight', 'Fulton',
                                       'Mg24w', 'Mn55w',	'Cu65w',	'Zn66w', 'Sr88w',	'Ba138w',
                                       'MgCa',	'MnCa',	'CuCa',	'ZnCa',	'SrCa',	'BaCa',	'Temp', 'Sal')), as.numeric)
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

otolith <- mutate_at(otolith, vars(c('Mg24', 'Mn55', 'Cu65', 'Zn66', 'Sr87', 'Ba138')), as.numeric)

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

Sr87 <- otolith %>%
  group_by(Month) %>%
  summarise(mean = mean(Sr87) )

Ba138 <- otolith %>%
  group_by(Month) %>%
  summarise(mean = mean(Ba138) )

cv(Mg24$mean, na.rm = TRUE) #cv = 0.11
cv(Mn55$mean, na.rm = TRUE) #cv = 0.38
cv(Cu65$mean, na.rm = TRUE) #cv = 0.30
cv(Zn66$mean, na.rm = TRUE) #cv = 0.52
cv(Sr87$mean, na.rm = TRUE) #cv = 0.04
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



#GLMMs
#Response variable transformed
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
#           (Intrc)    BaCas    df   logLik  AICc  delta weight
#model_sw1t  1.208e-15  0.71710  4  -98.953 206.2   0.00      1
#model_sw3t -5.676e-15  0.59490  4 -113.770 235.8  29.63      0
#model_sw2t  5.325e-15  0.48660  4 -130.950 270.1  63.99      0
#model_sw4t -3.987e-15  0.14460  4 -136.205 280.7  74.50      0
#model_sw5t  2.961e-15 -0.09401  4 -157.615 323.5 117.31      0
#model_sw6t  1.762e-17 -0.34330  4 -178.091 364.4 158.24      0
#model_sw0t -2.053e-17 -0.24440  4 -197.598 403.4 197.29      0
#Models ranked by AICc(x) 
#Random terms (all models): 
#  1 | ID

#r2 values
r.squaredGLMM(model_sw1t)

#R2m       R2c
#0.5156959 0.9246235

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
#model_sw1  9.650e-15  0.78780  4  -61.538 131.3   0.00      1
#model_sw2  1.236e-15  0.45980  4 -147.661 303.6 172.25      0
#model_sw3 -1.404e-15  0.62940  4 -160.319 328.9 197.56      0
#model_sw4 -1.080e-16  0.08716  4 -181.569 371.4 240.06      0
#model_sw0  1.079e-15 -0.20890  4 -187.656 383.6 252.24      0
#model_sw5  6.024e-16 -0.11280  4 -191.527 391.3 259.96      0
#model_sw6 -1.461e-15 -0.28650  4 -210.306 428.8 297.50      0
#Models ranked by AICc(x) 
#Random terms (all models): 
#  1 | ID

#r2 values
r.squaredGLMM(model_sw1)

#R2m       R2c
#0.6220642 0.9575027

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
#'Mg24w', 'Ca43w', 'Mn55w', 'Cu65w', 'Zn66w', 'Sr88w', and 'Ba138w' are in umol
#'MgCa', 'MnCa', 'CuCa', 'ZnCa', 'SrCa', and 'BaCa' are the ratio El:Ca43 in umol mol-1

#Filter water data
water <- water %>%
  filter(tide == "low") %>%
  filter(site == "Sao Mateus")

#Set as numeric
water <- mutate_at(water, vars(c('Mg24w', 'Ca43w', 'Mn55w', 'Cu65w', 'Zn66w', 'Sr88w', 'Ba138w',
                                 'MgCa', 'MnCa', 'CuCa', 'ZnCa', 'SrCa', 'BaCa',
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

otolith_1 <- otolith_1[,c(1,4,18:20,2:3,5:17)] #Reorder

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
sliding1 = merge(otolith_1, water, on='Sliding', all=FALSE)
str(sliding1)
sliding1$Month <- factor(sliding1$Month, levels = c('AUG','SEP','OCT','NOV','DEC','JAN','FEB','MAR','APR','MAY','JUN'))
sliding1 <- sliding1[,c(2,3,1,4:36)] #Reorder
length(unique(sliding1$ID))

#write.table(sliding1,"data_m1_e3.csv", sep=";", dec=".",row.names = F) #Month 1 correction (m1), Edge 3 otolith (m3)







############
#1.5 Figure 4 - one month correction
############

#Clean R environment 
rm(list = ls())

#Set working directory
setwd("C:/Users/patri/OneDrive/Documentos/UFES/Tese/From water to biominerals/R/github")

#read dataset
data <- read.csv2("correction.csv")

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

######
#Fig. 4. Measured otolith elemental signatures (green line, Ba:Caotolith) and water chemistry (Ba:Cawater) 
#without correction (blue dashed line) and with the best temporal correction (one month corrected, blue solid line).
#Temporal corrections resulted in the Ba:Cawater value corresponding to a calendar date of one month prior to 
#the Ba:Caotolith sampling (e.g., the Ba:Caotolith collected in August refers to the Ba:Cawater sampled in July). 
#Estimated scatterplot smoothing curves (lines) and 95% confidence intervals (shaded areas) were done by applying
#the local polynomial regression fitting (loess).


correction <- ggplot(data = data, aes(y = Ba138,x = Order))+
  geom_smooth(data = data, aes(y = Ba138,x = Order), color = "#008080",fill = "#008080",alpha=0.3) +
  geom_smooth(data = data, aes(y = BaCa_m*scaleBa, x = Order),color = "#00BFFF", size = 0.5,alpha=0.2, method = 'loess')+
  geom_smooth(data = data, aes(y = BaCa_n_m*scaleBa, x = Order),color = "#00BFFF", size = 0.5,alpha=0.2,linetype = "dashed", method = 'loess')+
  theme_bw()+
  scale_y_continuous(sec.axis = sec_axis(~./scaleBa, name=(expression(Ba:Ca[water]~(mmol~mol^-1)))))+
  scale_x_continuous(breaks=seq(1, 13, 1),labels = c('JUL','AUG','SEP','OCT','NOV','DEC','JAN','FEB','MAR','APR','MAY','JUN','JUL')) +
  labs(title = "One month corrected",
       x = "Otolith sampling", y = (expression(Ba:Ca[otolith]~(mu~mol~mol^-1))), color = "")+
  theme(axis.title.y = element_text(color = "#008080",size = 15),
        axis.title.y.right = element_text(color = "#00BFFF"))+
  theme(legend.position="none",
        title = element_text(size = 15),
        axis.text = element_text(size = 15),
        axis.text.x = element_text(angle = -45),
        axis.title = element_text(size = 15))
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

#Set as numeric
data <- mutate_at(data, vars(c('Mg24', 'Mn55', 'Cu65', 'Zn66', 'Sr87', 'Ba138',
                               'TL', 'Fulton',
                               'MgCa',	'MnCa',	'CuCa',	'ZnCa',	'SrCa','BaCa',
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

#Parametric test
model <- aov(TL ~ Month, data = data)
leveneTest(TL ~ Month, data = data) #p-value = 0.004123 **, we cannot assume the homogeneity of variance
shapiro.test(residuals(model)) #p-value = 0.0003134, we cannot assume the normality

summary(model)
tukey <- TukeyHSD(model)
tukey

#Non-parametric test
kruskal.test(TL ~ Month, data = data) #Kruskal-Wallis chi-squared = 89.163, df = 11, p-value = 2.431e-14
dunnTest(TL ~ Month, data = data) #post hoc

summary(data$TL)
sd(data$TL)
#Dog snapper’s TL was 173.7 ± 32.49 mm (mean ± standard deviation), ranging from 
#72.0 to 260.0 mm (Table 1), showing significant differences between the sampled months (χ2 = 89.163, p < .001). 


####
#Fulton
####

#Parametric test
model <- aov(Fulton ~ Month, data = data)
leveneTest(Fulton ~ Month, data = data) #p-value = 0.007608 **, we cannot assume the homogeneity of variance
shapiro.test(residuals(model)) #p-value = 1.948e-08, we cannot assume the normality

summary(model)
tukey <- TukeyHSD(model)
tukey

#Non-parametric test
kruskal.test(Fulton ~ Month, data = data) #Kruskal-Wallis chi-squared = 117.66, df = 11, p-value < 2.2e-16
dunnTest(Fulton ~ Month, data = data) #post hoc

summary(data$Fulton)
sd(data$Fulton)

#Their Fulton’s condition factor was 1.624 ± 0.152, ranging from 1.293 to 2.310 (Table 1),
#also displaying significant differences between the sampled months (χ2 = 117.66, p < .001











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
                                    'Ca43w',	'Mg24w',	'Mn55w',	'Cu65w',	'Zn66w', 'Sr88w',	'Ba138w',
                                    'MgCa',  'MnCa',	'CuCa',	'ZnCa',	'SrCa',	'BaCa',
                                    'Temp', 'Sal')]

gradient_l_wet<- gradient_l_wet[,c('Sample', 'Month', 'Year', 'Season', 
                                   'Ca43w',	'Mg24w',	'Mn55w',	'Cu65w',	'Zn66w', 'Sr88w',	'Ba138w',
                                   'MgCa',  'MnCa',	'CuCa',	'ZnCa',	'SrCa',	'BaCa',
                                   'Temp', 'Sal')]

#Set as numeric
gradient_e_wet <- mutate_at(gradient_e_wet, vars(c('Ca43w',	'Mg24w',	'Mn55w',	'Cu65w',	'Zn66w', 'Sr88w',	'Ba138w',
                                                   'MgCa',  'MnCa',	'CuCa',	'ZnCa',	'SrCa',	'BaCa',
                                                   'Temp', 'Sal')), as.numeric)

gradient_l_wet <- mutate_at(gradient_l_wet, vars(c('Ca43w',	'Mg24w',	'Mn55w',	'Cu65w',	'Zn66w', 'Sr88w',	'Ba138w',
                                                   'MgCa',  'MnCa',	'CuCa',	'ZnCa',	'SrCa',	'BaCa',
                                                   'Temp', 'Sal')), as.numeric)


###
# Figures
###


# Ca
Ca_e_wet<-ggplot(data = gradient_e_wet, aes(y = Ca43w, x = Sal)) +
  geom_point(color = "black", size=3)  +
  labs(title = "Ca",x = "Salinity", y = "Ca in water")  +
  theme_bw()+
  geom_smooth(method="lm", formula = y ~ splines::bs(x, 3), se = FALSE,color='#008080')+
  scale_x_continuous(breaks = seq(0, 35, by = 5))
Ca_e_wet

Ca_l_wet<-ggplot(data = gradient_l_wet, aes(y = Ca43w, x = Sal)) +
  geom_point(color = "black", size=3)  +
  labs(title = "Ca",x = "Salinity", y = "Ca in water")  +
  theme_bw()+
  geom_smooth(method="lm", formula = y ~ splines::bs(x, 3), se = FALSE,color='#00BFFF')+
  scale_x_continuous(breaks = seq(0, 35, by = 5))
Ca_l_wet

###
# Elements
###

# Mg
gradient_e_wet$MgCam <- gradient_e_wet$MgCa/1000000 #transform umol mol-1 to mol mol-1
gradient_l_wet$MgCam <- gradient_l_wet$MgCa/1000000 #transform umol mol-1 to mol mol-1

# Mg:Ca
Mg_e_wet <- ggplot(data = gradient_e_wet, aes(y = MgCam, x = Sal)) +
  geom_point(size=3,color='#008080', alpha = 0.5)  +
  labs(title = " ",x = "Salinity", y = expression(Mg:Ca[water]~ ('mol'~'mol'^-1)))  +
  theme_bw()+
  geom_smooth(method="lm", formula = y ~ splines::bs(x, 3), se = FALSE,color='#008080')+
  scale_x_continuous(breaks = seq(0, 35, by = 5))

Mg_l_wet <- ggplot(data = gradient_l_wet, aes(y = MgCam, x = Sal)) +
  geom_point(size=3,color='#00BFFF', alpha = 0.5)  +
  labs(title = " ",x = "Salinity", y = expression(Mg:Ca[water]~ ('mol'~'mol'^-1)))  +
  theme_bw()+
  geom_smooth(method="lm", formula = y ~ splines::bs(x, 3), se = FALSE,color='#00BFFF')+
  scale_x_continuous(breaks = seq(0, 35, by = 5))

Mg_e_wet
Mg_l_wet

Mg <- ggplot(data = gradient_l_wet, aes(y = MgCam, x = Sal)) +
  geom_point(size=3,color='#00BFFF', alpha = 0.3)  +
  geom_smooth(method="lm", formula = y ~ splines::bs(x, 3), se = FALSE,color='#00BFFF')+
  geom_point(data = gradient_e_wet, aes(y = MgCam, x = Sal), size=3,color='#008080', alpha = 0.3)  +
  geom_smooth(data = gradient_e_wet, aes(y = MgCam, x = Sal), method="lm", formula = y ~ splines::bs(x, 3), se = FALSE,color='#008080')+
  labs(title = " ",x = "Salinity", y = expression(Mg:Ca[water]~ ('mol'~'mol'^-1)))  +
  theme_bw()+
  scale_x_continuous(breaks = seq(0, 35, by = 5))
Mg

###
# Mn
###
gradient_e_wet$MnCam <- gradient_e_wet$MnCa/1000 #transform umol mol-1 to mmol mol-1
gradient_l_wet$MnCam <- gradient_l_wet$MnCa/1000 #transform umol mol-1 to mmol mol-1

#55Mn:Ca
Mn_e_wet<- ggplot(data = gradient_e_wet, aes(y = MnCam, x = Sal)) +
  geom_point(size=3,color='#008080', alpha = 0.5)  +
  labs(title = " ",x = "Salinity", y = expression(Mn:Ca[water]~ ('mmol'~'mol'^-1)))  +
  theme_bw()+
  geom_smooth(method="lm", formula = y ~ splines::bs(x, 3), se = FALSE,color='#008080')+
  scale_x_continuous(breaks = seq(0, 35, by = 5))

Mn_l_wet<- ggplot(data = gradient_l_wet, aes(y = MnCam, x = Sal)) +
  geom_point(size=3,color='#00BFFF', alpha = 0.5)  +
  labs(title = " ",x = "Salinity", y = expression(Mn:Ca[water]~ ('mmol'~'mol'^-1)))  +
  theme_bw()+
  geom_smooth(method="lm", formula = y ~ splines::bs(x, 3), se = FALSE,color='#00BFFF')+
  scale_x_continuous(breaks = seq(0, 35, by = 5))

Mn_e_wet
Mn_l_wet

Mn <- ggplot(data = gradient_l_wet, aes(y = MnCam, x = Sal)) +
  geom_point(size=3,color='#00BFFF', alpha = 0.3)  +
  geom_smooth(method="lm", formula = y ~ splines::bs(x, 3), se = FALSE,color='#00BFFF')+
  geom_point(data = gradient_e_wet, aes(y = MnCam, x = Sal), size=3,color='#008080', alpha = 0.3)  +
  geom_smooth(data = gradient_e_wet, aes(y = MnCam, x = Sal), method="lm", formula = y ~ splines::bs(x, 3), se = FALSE,color='#008080')+
  labs(title = " ",x = "Salinity", y = expression(Mn:Ca[water]~ ('mmol'~'mol'^-1)))  +
  theme_bw()+
  scale_x_continuous(breaks = seq(0, 35, by = 5))
Mn

###
# Cu
###

# Cu:Ca
Cu_e_wet<- ggplot(data = gradient_e_wet, aes(y = CuCa, x = Sal)) +
  geom_point(size=3,color='#008080', alpha = 0.5)  +
  labs(title = " ",x = "Salinity", y = expression(Cu:Ca[water]~ (mu~mol~mol^-1)))  +
  theme_bw()+
  geom_smooth(method="lm", formula = y ~ splines::bs(x, 3), se = FALSE,color='#008080')+
  scale_x_continuous(breaks = seq(0, 35, by = 5))

Cu_l_wet<- ggplot(data = gradient_l_wet, aes(y = CuCa, x = Sal)) +
  geom_point(size=3,color='#00BFFF', alpha = 0.5)  +
  labs(title = " ",x = "Salinity", y = expression(Cu:Ca[water]~ (mu~mol~mol^-1)))  +
  theme_bw()+
  geom_smooth(method="lm", formula = y ~ splines::bs(x, 3), se = FALSE,color='#00BFFF')+
  scale_x_continuous(breaks = seq(0, 35, by = 5))

Cu_e_wet
Cu_l_wet

Cu <- ggplot(data = gradient_l_wet, aes(y = CuCa, x = Sal)) +
  geom_point(size=3,color='#00BFFF', alpha = 0.3)  +
  geom_smooth(method="lm", formula = y ~ splines::bs(x, 3), se = FALSE,color='#00BFFF')+
  geom_point(data = gradient_e_wet, aes(y = CuCa, x = Sal), size=3,color='#008080', alpha = 0.3)  +
  geom_smooth(data = gradient_e_wet, aes(y = CuCa, x = Sal), method="lm", formula = y ~ splines::bs(x, 3), se = FALSE,color='#008080')+
  labs(title = " ",x = "Salinity", y = expression(Cu:Ca[water]~ (mu~mol~mol^-1)))  +
  theme_bw()+
  scale_x_continuous(breaks = seq(0, 35, by = 5))
Cu

###
# Zn
###
gradient_e_wet$ZnCam <- gradient_e_wet$ZnCa/1000 #transform umol mol-1 to mmol mol-1
gradient_l_wet$ZnCam <- gradient_l_wet$ZnCa/1000 #transform umol mol-1 to mmol mol-1

# Zn:Ca
Zn_e_wet<- ggplot(data = gradient_e_wet, aes(y = ZnCam, x = Sal)) +
  geom_point(size=3,color='#008080', alpha = 0.3)  +
  labs(title = " ",x = "Salinity", y = Zn:Ca[water]~ (mmol~mol^-1))  +
  theme_bw()+
  geom_smooth(method="lm", formula = y ~ splines::bs(x, 3), se = FALSE,color='#008080')+
  scale_x_continuous(breaks = seq(0, 35, by = 5))

Zn_l_wet<- ggplot(data = gradient_l_wet, aes(y = ZnCam, x = Sal)) +
  geom_point(size=3,color='#00BFFF', alpha = 0.3)  +
  labs(title = " ",x = "Salinity", y = Zn:Ca[water]~ (mmol~mol^-1))  +
  theme_bw()+
  geom_smooth(method="lm", formula = y ~ splines::bs(x, 3), se = FALSE,color='#00BFFF')+
  scale_x_continuous(breaks = seq(0, 35, by = 5))

Zn_e_wet
Zn_l_wet

Zn <- ggplot(data = gradient_l_wet, aes(y = ZnCam, x = Sal)) +
  geom_point(size=3,color='#00BFFF', alpha = 0.3)  +
  geom_smooth(method="lm", formula = y ~ splines::bs(x, 3), se = FALSE,color='#00BFFF')+
  geom_point(data = gradient_e_wet, aes(y = ZnCam, x = Sal), size=3,color='#008080', alpha = 0.3)  +
  geom_smooth(data = gradient_e_wet, aes(y = ZnCam, x = Sal), method="lm", formula = y ~ splines::bs(x, 3), se = FALSE,color='#008080')+
  labs(title = " ",x = "Salinity", y = expression(Zn:Ca[water]~ (mmol~mol^-1)))  +
  theme_bw()+
  scale_x_continuous(breaks = seq(0, 35, by = 5))
Zn

###
# Sr
###
gradient_e_wet$SrCam <- gradient_e_wet$SrCa/1000 #transform umol mol-1 to mmol mol-1
gradient_l_wet$SrCam <- gradient_l_wet$SrCa/1000 #transform umol mol-1 to mmol mol-1

# Sr:Ca
Sr_e_wet<- ggplot(data = gradient_e_wet, aes(y = SrCam, x = Sal)) +
  geom_point(size=3,color='#008080', alpha = 0.3)  +
  labs(title = " ",x = "Salinity", y = expression(Sr:Ca[water]~ (mmol~mol^-1)))  +
  theme_bw()+
  geom_smooth(method="lm", formula = y ~ splines::bs(x, 3), se = FALSE,color='#008080')+
  scale_x_continuous(breaks = seq(0, 35, by = 5))

Sr_l_wet<- ggplot(data = gradient_l_wet, aes(y = SrCam, x = Sal)) +
  geom_point(size=3,color='#00BFFF', alpha = 0.3)  +
  labs(title = " ",x = "Salinity", y = expression(Sr:Ca[water]~ (mmol~mol^-1)))  +
  theme_bw()+
  geom_smooth(method="lm", formula = y ~ splines::bs(x, 3), se = FALSE,color='#00BFFF')+
  scale_x_continuous(breaks = seq(0, 35, by = 5))

Sr_e_wet
Sr_l_wet

Sr <- ggplot(data = gradient_l_wet, aes(y = SrCam, x = Sal)) +
  geom_point(size=3,color='#00BFFF', alpha = 0.3)  +
  geom_smooth(method="lm", formula = y ~ splines::bs(x, 3), se = FALSE,color='#00BFFF')+
  geom_point(data = gradient_e_wet, aes(y = SrCam, x = Sal), size=3,color='#008080', alpha = 0.3)  +
  geom_smooth(data = gradient_e_wet, aes(y = SrCam, x = Sal), method="lm", formula = y ~ splines::bs(x, 3), se = FALSE,color='#008080')+
  labs(title = " ",x = "Salinity", y = expression(Sr:Ca[water]~ (mmol~mol^-1)))  +
  theme_bw()+
  scale_x_continuous(breaks = seq(0, 35, by = 5))
Sr

# Sr Water

Srwater <- ggplot(data = gradient_e_wet, aes(y = Sr88w, x = Sal)) +
  geom_point(size=3,color='#008080', alpha = 0.3)  +
  geom_smooth(data = gradient_e_wet, aes(y = Sr88w), method="lm", formula = y ~ splines::bs(x, 3), se = FALSE, color = "#008080")+
  geom_point(data = gradient_l_wet, aes(y = Sr88w), color='#00BFFF',size=3,shape=19, alpha = 0.3)  +
  geom_smooth(data = gradient_l_wet, aes(y = Sr88w), method="lm", formula = y ~ splines::bs(x, 3), se = FALSE, color = "#00BFFF")+
  labs(title = " ",x = "Salinity", y = "Sr water")  +
  theme_bw()

Srwater

###
# Ba
####
gradient_e_wet$BaCam <- gradient_e_wet$BaCa/1000 #transform umol mol-1 to mmol mol-1
gradient_l_wet$BaCam <- gradient_l_wet$BaCa/1000 #transform umol mol-1 to mmol mol-1

#138Ba:Ca
Ba_e_wet<- ggplot(data = gradient_e_wet, aes(y = BaCam, x = Sal)) +
  geom_point(size=3,color='#008080', alpha = 0.3)  +
  labs(title = " ",x = "Salinity", y = expression(Ba:Ca[water]~ (mmol~mol^-1)))  +
  theme_bw()+
  geom_smooth(method="lm", formula = y ~ splines::bs(x, 3), se = FALSE,color='#008080')+
  scale_x_continuous(breaks = seq(0, 35, by = 5))

Ba_l_wet<- ggplot(data = gradient_l_wet, aes(y = BaCam, x = Sal)) +
  geom_point(size=3,color='#00BFFF', alpha = 0.3)  +
  labs(title = " ",x = "Salinity", y = expression(Ba:Ca[water]~ (mmol~mol^-1)))  +
  theme_bw()+
  geom_smooth(method="lm", formula = y ~ splines::bs(x, 3), se = FALSE,color='#00BFFF')+
  scale_x_continuous(breaks = seq(0, 35, by = 5))

Ba_e_wet
Ba_l_wet

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
# Fig. 5. Water chemistry (El:Cawater) across the salinity gradient (from zero to 35) collected during the
#early wet (green line) and late wet (blue line) seasons at the São Mateus estuary – Southwestern Atlantic. 
#The dots indicate raw data, while estimated curves (lines) were done by applying linear models using a 
#third-degree basis spline.
###

ggarrange(Mg,Mn,Cu,Zn,Sr,Ba, 
          labels = c("A","B","C","D","E","F"),
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
data <- mutate_at(data, vars(c('Mg24', 'Mn55', 'Cu65', 'Zn66', 'Sr87', 'Ba138',
                               'TL', 'Fulton',
                               'MgCa',	'MnCa',	'CuCa',	'ZnCa',	'SrCa','BaCa',
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

# Sr
summary(data$Sr87)
sd(data$Sr87)
#Sr:Caotolith ranged from 1837.13 to 4126.84 μmol mol−1 (3026.43 ± 368.33 μmol mol−1)

# Ba
summary(data$Ba138)
sd(data$Ba138)
#Ba:Caotolith ranged from 0.69 to 161.89 μmol mol−1 (19.54 ± 22.79 μmol mol−1)



########
# Kruskal–Wallis’ test and post hoc  Dunn’s test
########


# Mg
#Parametric test
model <- aov(Mg24 ~ Month, data = data_Mg)
leveneTest(Mg24 ~ Month, data = data_Mg) # 0.001818 **, we cannot assume the homogeneity of variance
shapiro.test(residuals(model)) #p-value = 0.03072, we cannot assume the normality

#Non-parametric test
kruskal.test(Mg24 ~ Month, data = data_Mg) #Kruskal-Wallis chi-squared = 52.335, df = 11, p-value = 2.371e-07
dunnTest(Mg24 ~ Month, data = data_Mg,method="bonferroni") #post hoc


# Mn
#Parametric test
model <- aov(Mn55 ~ Month, data = data_Mn)
leveneTest(Mn55 ~ Month, data = data_Mn) # 0.0006715 ***, we cannot assume the homogeneity of variance
shapiro.test(residuals(model)) #p-value = 2.91e-12, we cannot assume the normality

#Non-parametric test
kruskal.test(Mn55 ~ Month, data = data_Mn) #Kruskal-Wallis chi-squared = 35.545, df = 11, p-value = 0.0002014
dunnTest(Mn55 ~ Month, data = data_Mn,method="bonferroni") #post hoc


# Cu
#Parametric test
model <- aov(Cu65 ~ Month, data = data_Cu)
leveneTest(Cu65 ~ Month, data = data_Cu) # 0.001871 **, we cannot assume the homogeneity of variance
shapiro.test(residuals(model)) #p-value < 2.2e-16, we cannot assume the normality

#Non-parametric test
kruskal.test(Cu65 ~ Month, data = data_Cu) #Kruskal-Wallis chi-squared = 26.199, df = 11, p-value = 0.006061
dunnTest(Cu65 ~ Month, data = data_Cu,method="bonferroni") #post hoc


# Zn
#Parametric test
model <- aov(Zn66 ~ Month, data = data_Zn)
leveneTest(Zn66 ~ Month, data = data_Zn) # 0.1546, we can assume the homogeneity of variance
shapiro.test(residuals(model)) #p-value < 2.2e-16, we cannot assume the normality

#Non-parametric test
kruskal.test(Zn66 ~ Month, data = data_Zn) #Kruskal-Wallis chi-squared = 30.523, df = 11, p-value = 0.001311
dunnTest(Zn66 ~ Month, data = data_Zn,method="bonferroni") #post hoc


# Sr
#Parametric test
model <- aov(Sr87 ~ Month, data = data)
leveneTest(Sr87 ~ Month, data = data) #0.001297 **, we cannot assume the homogeneity of variance
shapiro.test(residuals(model)) #p-value = 0.199, we can assume the normality

#Non-parametric test
kruskal.test(Sr87 ~ Month, data = data) #Kruskal-Wallis chi-squared = 36.222, df = 11, p-value = 0.0001554
dunnTest(Sr87 ~ Month, data = data,method="bonferroni") #post hoc


# Ba
#Parametric test
model <- aov(Ba138 ~ Month, data = data)
leveneTest(Ba138 ~ Month, data = data) # < 2.2e-16 ***, we cannot assume the homogeneity of variance
shapiro.test(residuals(model)) #p-value < 2.2e-16, we cannot assume the normality

#Non-parametric test
kruskal.test(Ba138 ~ Month, data = data) #Kruskal-Wallis chi-squared = 217.65, df = 11, p-value < 2.2e-16
dunnTest(Ba138 ~ Month, data = data,method="bonferroni") #post hoc









########
# 2.4.2 Figure 6 - Otolith elemental signatures
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
data <- mutate_at(data, vars(c('Mg24', 'Mn55', 'Cu65','Zn66',  'Sr87', 'Ba138', 'TL', 'Fulton',
                               'MgCa',	'MnCa',	'CuCa',	'ZnCa', 'SrCa',	'BaCa',	
                               'Temp', 'Sal','OrderS')), as.numeric)

data$Month <- factor(data$Month, levels = c('JUL','AUG','SEP','OCT','NOV','DEC','JAN','FEB','MAR','APR','MAY','JUN'))
data$Sliding <- factor(data$Sliding, levels = c('JUL','AUG','SEP','OCT','NOV','DEC','JAN','FEB','MAR','APR','MAY'))
summary(data$Month)
summary(data$Sliding)


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
  annotate('rect', xmin=3.5, xmax=6.5, ymin=80, ymax=450, alpha=0.15, fill='#999999')+
  annotate('rect', xmin=9.5, xmax=12.5, ymin=80, ymax=450, alpha=0.15, fill='#999999')+
  geom_boxplot(alpha = 1,width=0.4, color = "#008080") +
  geom_line(data = data_Mg, aes(y = MgCam*scaleMg, x = Order),color = "#00BFFF", size = 0.3,alpha=1)+
  geom_point(aes(y = MgCam*scaleMg, x = Order),color = "#00BFFF")+
  theme_bw() +
  theme(legend.position="none",
        axis.text = element_text(size = 15),
        axis.text.x = element_text(angle = -45),
        axis.title = element_text(size = 15),
        axis.title.y = element_text(color = "#008080",size = 15),
        axis.title.y.right = element_text(color = "#00BFFF")) +
  xlab("") + ylab(expression(Mg:Ca[otolith]~(mu~mol~mol^-1)))+
  scale_y_continuous(sec.axis = sec_axis(~./scaleMg, name=expression(Mg:Ca[water]~ ('mol'~'mol'^-1))))+
  geom_text(data = letter_Mg, aes(x = Month, y = quant, label = Letter), vjust=-8,  size = 4.5, color = "#008080")

Mg24

# Mn
#Delete outliers
data_Mn <- data
data_Mn <- data[data$ID != c('ID165'), ]
data_Mn <- data_Mn %>%  filter(!(ID == 'ID382' & Time == '81.5'))
data_Mn <- data_Mn %>%  filter(!(ID == 'ID327' & Time == '70.2'))

#transform umol mol-1 to mmol mol-1
data_Mn$MnCam <- data_Mn$MnCa/1000 

#set scale
scaleMn<- 12

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
  geom_line(data = data_Mn, aes(y = MnCam/scaleMn, x = Order),color = "#00BFFF", size = 0.3,alpha=1)+
  geom_point(aes(y = MnCam/scaleMn, x = Order),color = "#00BFFF")+
  theme_bw() +
  theme(legend.position="none",
        axis.text = element_text(size = 15),
        axis.text.x = element_text(angle = -45),
        axis.title.y = element_text(color = "#008080",size = 15),
        axis.title.y.right = element_text(color = "#00BFFF")) +
  xlab("") + ylab(expression(Mn:Ca[otolith]~(mu~mol~mol^-1)))+
  scale_y_continuous(sec.axis = sec_axis(~.*scaleMn, name=expression(Mn:Ca[water]~ ('mmol'~'mol'^-1))))+
  geom_text(data = letter_Mn, aes(x = Month, y = quant, label = Letter), vjust=-8, size = 4.5,color = "#008080")

Mn55

# Cu
#Delete outliers
data_Cu <- data
data_Cu <- data_Cu %>%  filter(!(ID == 'ID083' & Time == '72.4'))
data_Cu <- data_Cu %>%  filter(!(ID == 'ID072' & Time == '81.5'))
data_Cu <- data_Cu %>%  filter(!(ID == 'ID076' & Time == '77'))
data_Cu <- data_Cu %>%  filter(!(ID == 'ID202' & Time == '70.2'))
data_Cu <- data_Cu %>%  filter(!(ID == 'ID202' & Time == '72.4'))
data_Cu <- data_Cu %>%  filter(!(ID == 'ID277' & Time == '72.4'))

#set scale
scaleCu <- 300

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
  geom_line(data = data_Cu, aes(y = CuCa/scaleCu, x = Order),color = "#00BFFF", size = 0.3,alpha=1)+
  geom_point(aes(y = CuCa/scaleCu, x = Order),color = "#00BFFF")+
  theme_bw() +
  theme(legend.position="none",
        axis.text = element_text(size = 15),
        axis.text.x = element_text(angle = -45),
        axis.title = element_text(size = 15),
        axis.title.y = element_text(color = "#008080",size = 15),
        axis.title.y.right = element_text(color = "#00BFFF")) +
  xlab("") + ylab(expression(Cu:Ca[otolith]~(mu~mol~mol^-1)))+
  scale_y_continuous(sec.axis = sec_axis(~.*scaleCu, name=expression(Cu:Ca[water]~ (mu~mol~mol^-1))))+
  geom_text(data = letter_Cu, aes(x = Month, y = quant, label = Letter), vjust=-8, size = 4.5,color = "#008080")

Cu65

# Zn
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
scaleZn<-1

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
  annotate('rect', xmin=3.5, xmax=6.5, ymin=0, ymax=8, alpha=0.15, fill='#999999')+
  annotate('rect', xmin=9.5, xmax=12.5, ymin=0, ymax=8, alpha=0.15, fill='#999999')+
  geom_boxplot(alpha = 1,width=0.4,color = "#008080") +
  geom_line(data = data_Zn, aes(y = ZnCam*scaleZn, x = Order),color = "#00BFFF", size = 0.3,alpha=1)+
  geom_point(aes(y = ZnCam*scaleZn, x = Order),color = "#00BFFF")+
  theme_bw() +
  theme(legend.position="none",
        axis.text = element_text(size = 15),
        axis.text.x = element_text(angle = -45),
        axis.title = element_text(size = 15),
        axis.title.y = element_text(color = "#008080",size = 15),
        axis.title.y.right = element_text(color = "#00BFFF")) +
  xlab("") + ylab(expression(Zn:Ca[otolith]~(mu~mol~mol^-1)))+
  scale_y_continuous(sec.axis = sec_axis(~./scaleZn, name=expression(Zn:Ca[water]~ (mmol~mol^-1))))+
  geom_text(data = letter_Zn, aes(x = Month, y = quant, label = Letter), vjust=-8, size = 4.5,color = "#008080")

Zn66


# Sr
data_Sr <- data

#transform umol mol-1 to mmol mol-1
data_Sr$Sr87m <- data$Sr87/1000
data_Sr$SrCam <- data$SrCa/1000

#set scale
scaleSr<-20

#letters based on post hoc Dunn's test
dunn_Sr <- dunnTest(Sr87 ~ Month, data = data_Sr, method = "bonferroni")$res
dunn_Sr
cld_Sr <- cldList(P.adj ~ Comparison, data=dunn_Sr)
cld_Sr
cld_Sr <- as.data.frame(cld_Sr)                    
cld_Sr <- cld_Sr %>%
  rename(Month = Group)
cld_Sr

letter_Sr <- group_by(data_Sr, Month) %>%
  summarise(mean=mean(Sr87m), quant = quantile(Sr87m, probs = 0.75)) %>%
  arrange(desc(mean))
letter_Sr <- as.data.frame(letter_Sr)                    
letter_Sr

letter_Sr <- merge(letter_Sr, cld_Sr, by = "Month")
letter_Sr

Sr87 <-  ggplot(data_Sr, aes(x=Month, y=Sr87m)) +
  annotate('rect', xmin=3.5, xmax=6.5, ymin=1.7, ymax=4.7, alpha=0.15, fill='#999999')+
  annotate('rect', xmin=9.5, xmax=12.5, ymin=1.7, ymax=4.7, alpha=0.15, fill='#999999')+
  geom_boxplot(alpha = 1,width=0.4,color = "#008080") +
  geom_line(data = data_Sr, aes(y = SrCam/scaleSr, x = Order),color = "#00BFFF", size = 0.3,alpha=1)+
  geom_point(aes(y = SrCam/scaleSr, x = Order),color = "#00BFFF")+
  theme_bw() +
  theme(legend.position="none",
        axis.text = element_text(size = 15),
        axis.text.x = element_text(angle = -45),
        axis.title = element_text(size = 15),
        axis.title.y = element_text(color = "#008080",size = 15),
        axis.title.y.right = element_text(color = "#00BFFF")) +
  xlab("") + ylab(expression(Sr:Ca[otolith]~(mmol~mol^-1)))+
  scale_y_continuous(limits = c(1.7, 4.7), sec.axis = sec_axis(~.*scaleSr, name=expression(Sr:Ca[water]~ (mmol~mol^-1))))+
  geom_text(data = letter_Sr, aes(x = Month, y = quant, label = Letter), vjust=-8, size = 4.5,color = "#008080")

Sr87


# Ba
data_Ba<-data

#transform umol mol-1 to mmol mol-1
data_Ba$BaCam <- data$BaCa/1000 

#set scale
scaleBa<-2.5

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
  geom_line(data = data_Ba, aes(y = BaCam*scaleBa, x = Order),color = "#00BFFF", size = 0.3,alpha=1)+
  geom_point(aes(y = BaCam*scaleBa, x = Order),color = "#00BFFF")+
  theme_bw() +
  theme(legend.position="none",
        axis.text = element_text(size = 15),
        axis.text.x = element_text(angle = -45),
        axis.title = element_text(size = 15),
        axis.title.y = element_text(color = "#008080",size = 15),
        axis.title.y.right = element_text(color = "#00BFFF")) +
  xlab("") + ylab(expression(Ba:Ca[otolith]~(mu~mol~mol^-1)))+
  scale_y_continuous(sec.axis = sec_axis(~./scaleBa, name=expression(Ba:Ca[water]~ (mmol~mol^-1))))+
  geom_text(data = letter_Ba, aes(x = Month, y = quant, label = Letter), vjust=-8, size = 4.5,color = "#008080")

Ba138


#Fig. 6. Otolith elemental signatures (El:Caotolith) of the juvenile dog snappers (green boxplots) and water 
#chemistry (El:Cawater, blue lines) monthly sampled at the São Mateus estuary – Southwestern Atlantic. 
#El:Cawater is presented with a one-month correction. Shade and transparent areas indicate the distinct seasons:
#late dry (from July to September), early wet (from October to December), late wet (from January to March),
#and early dry (from April to June). Lettering indicates post hoc Dunn’s test results in 
#which months are significantly different (p < 0.05), considering El:Caotolith as the response variable, 
#and the sampled month as the predictor.

ggarrange(Mg24,Mn55,Cu65,Zn66,Sr87,Ba138, 
          labels = c("A","B","C","D","E","F"),
          ncol = 2, nrow =3,
          legend = "none")
#1200 x 1200










########
# 2.4.3 Figure 7 - Correlation between otolith and water chemistries
########

#Clean R environment 
rm(list = ls())

water <- read.csv2("water_chemistry.csv")
water <- as.data.frame(water)
str(water)

water <- water %>%
  filter(site == "Sao Mateus") %>%
  filter(tide == "low")

water <- water[water$Month != 'JUN', ]

water <- water[,c('Month','MgCa','MnCa','CuCa',	'ZnCa','SrCa','BaCa')]

#Set as numeric
water <- mutate_at(water, vars(c('MgCa','MnCa','CuCa',	'ZnCa','SrCa','BaCa')), as.numeric)
water <- water %>%
  rename(Sliding = Month)

#Read data
data <- read.csv2("data_m1_e3.csv")
length(unique(data$ID))
str(data)

data <- data[data$Month != c('JUL'), ]

#Set as numeric
data <- data[,c('Month','Sliding','ID','Time','Mg24','Mn55','Cu65',	'Zn66','Sr87','Ba138')]
data <- mutate_at(data, vars(c('Mg24','Mn55','Cu65',	'Zn66','Sr87','Ba138')), as.numeric)
data$Month <- factor(data$Month, levels = c('JUL','AUG','SEP','OCT','NOV','DEC','JAN','FEB','MAR','APR','MAY','JUN'))
data$Sliding <- factor(data$Sliding, levels = c('JUL','AUG','SEP','OCT','NOV','DEC','JAN','FEB','MAR','APR','MAY'))
summary(data$Month)
summary(data$Sliding)
str(data)


# Mg
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
## Get the mean
mean_Mg <- data_Mg %>%
  group_by(Sliding) %>%
  summarise(Mg24 = mean(Mg24))
str(mean_Mg)
#### Join dataframes
mean_Mg = merge(mean_Mg, water, on='Sliding', all=FALSE)
str(mean_Mg)

# Mn
data_Mn <- data
data_Mn <- data[data$ID != c('ID165'), ]
data_Mn <- data_Mn %>%  filter(!(ID == 'ID382' & Time == '81.5'))
data_Mn <- data_Mn %>%  filter(!(ID == 'ID327' & Time == '70.2'))
## Get the mean
mean_Mn <- data_Mn %>%
  group_by(Sliding) %>%
  summarise(Mn55 = mean(Mn55))
#### Join dataframes
mean_Mn = merge(mean_Mn, water, on='Sliding', all=FALSE)
str(mean_Mg)

# Cu
data_Cu <- data
data_Cu <- data_Cu %>%  filter(!(ID == 'ID083' & Time == '72.4'))
data_Cu <- data_Cu %>%  filter(!(ID == 'ID072' & Time == '81.5'))
data_Cu <- data_Cu %>%  filter(!(ID == 'ID076' & Time == '77'))
data_Cu <- data_Cu %>%  filter(!(ID == 'ID202' & Time == '70.2'))
data_Cu <- data_Cu %>%  filter(!(ID == 'ID202' & Time == '72.4'))
data_Cu <- data_Cu %>%  filter(!(ID == 'ID277' & Time == '72.4'))
## Get the mean
mean_Cu <- data_Cu %>%
  group_by(Sliding) %>%
  summarise(Cu65 = mean(Cu65))
#### Join dataframes
mean_Cu = merge(mean_Cu, water, on='Sliding', all=FALSE)
str(mean_Cu)

# Zn
data_Zn <- data
data_Zn <- data_Zn %>%  filter(!(ID == 'ID386' & Time == '72.4'))
data_Zn <- data_Zn %>%  filter(!(ID == 'ID079' & Time == '72.4'))
data_Zn <- data_Zn %>%  filter(!(ID == 'ID509' & Time == '83.8'))
data_Zn <- data_Zn %>%  filter(!(ID == 'ID511' & Time == '83.8'))
data_Zn <- data_Zn %>%  filter(!(ID == 'ID508' & Time == '47.5'))
## Get the mean
mean_Zn <- data_Zn %>%
  group_by(Sliding) %>%
  summarise(Zn66 = mean(Zn66))
#### Join dataframes
mean_Zn = merge(mean_Zn, water, on='Sliding', all=FALSE)
str(mean_Zn)

# Sr
## Get the mean
mean_Sr <- data %>%
  group_by(Sliding) %>%
  summarise(Sr87 = mean(Sr87))
str(mean_Sr)
#### Join dataframes
mean_Sr = merge(mean_Sr, water, on='Sliding', all=FALSE)
str(mean_Sr)

#Ba
## Get the mean
mean_Ba <- data %>%
  group_by(Sliding) %>%
  summarise(Ba138 = mean(Ba138))
#### Join dataframes
mean_Ba = merge(mean_Ba, water, on='Sliding', all=FALSE)
str(mean_Ba)



#Figures

# Mg
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
cor(mean_Mg$Mg24, mean_Mg$MgCa, method = "pearson")

# Mn
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
cor(mean_Mn$Mn55, mean_Mn$MnCa, method = "pearson")

# Cu
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
cor(mean_Cu$Cu65, mean_Cu$CuCa, method = "pearson")

# Zn
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
cor(mean_Zn$Zn66, mean_Zn$ZnCa, method = "pearson")

# Sr
mean_Sr$SrCam <- mean_Sr$SrCa/1000 #transform umol mol-1 to mmol mol-1
mean_Sr$Sr87m <- mean_Sr$Sr87/1000 #transform umol mol-1 to mmol mol-1
Sr <- ggplot(mean_Sr, aes(y=Sr87m, x=SrCam)) +
  theme_bw() +
  geom_smooth(method=glm,color = "#00BFFF", fill = "#00BFFF",alpha=0.1)  +
  labs(x = expression(Sr:Ca[water]~ (mmol~mol^-1)), y = (expression(Sr:Ca[otolith]~(mmol~mol^-1))), color = "")+
  geom_point(size=3) +
  theme(legend.position="none",
        title = element_text(size = 15),
        axis.text = element_text(size = 15),
        axis.title = element_text(size = 15))+ 
  stat_cor(method = "pearson",label.y = 3.25)
Sr
cor(mean_Sr$Sr87, mean_Sr$SrCa, method = "pearson")

# Ba
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
cor(mean_Ba$Ba138, mean_Ba$BaCa, method = "pearson")

####
# Fig. 7. Correlation between the mean otolith elemental signatures (El:Caotolith) of the juvenile dog 
#snappers and water chemistry (El:Cawater) monthly sampled at the São Mateus estuary – Southwestern Atlantic.
#Pearson’s correlation (R) and p values of the relationships are indicated at the top left for each element. 
#The correlation was calculated considering El:Cawater with a one-month correction.
####
ggarrange(Mg,Mn,Cu,Zn,Sr,Ba, 
          labels = c("A","B","C","D","E","F"),
          ncol = 3, nrow = 2,
          legend = "none")
#1200 x 600










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

#Set as numeric
data <- mutate_at(data, vars(c('Mg24', 'Mn55', 'Cu65', 'Zn66', 'Sr87','Ba138', 'TL', 'Fulton',
                               'Mg24w','Mn55w',	'Cu65w',	'Zn66w',	'Sr88w',	'Ba138w',
                               'MgCa','MnCa',	'CuCa',	'ZnCa',	'SrCa',	'BaCa','Temp', 'Sal','OrderS')), as.numeric)

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
data$Sr_coef <- data$Sr87/data$SrCa
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
#D Mg:Ca ratios ranged from 0.03x10-4 to 0.23x10-4 (0.07x10-4 ± 0.04x10-4)

#Mn_coef
data_Mn <- data
data_Mn <- data[data_Mn$ID != c('ID165'), ]
data_Mn <- data_Mn %>%  filter(!(ID == 'ID382' & Time == '81.5'))
data_Mn <- data_Mn %>%  filter(!(ID == 'ID327' & Time == '70.2'))

summary(data_Mn$Mn_coef)
sd(data_Mn$Mn_coef)
#D Mn:Ca ratios ranged from 0.0003x10-3 to 1.33x10-3 (0.19x10-3 ± 0.24x10-3)

#Cu_coef
data_Cu <- data
data_Cu <- data_Cu %>%  filter(!(ID == 'ID083' & Time == '72.4'))
data_Cu <- data_Cu %>%  filter(!(ID == 'ID072' & Time == '81.5'))
data_Cu <- data_Cu %>%  filter(!(ID == 'ID076' & Time == '77'))
data_Cu <- data_Cu %>%  filter(!(ID == 'ID202' & Time == '70.2'))
data_Cu <- data_Cu %>%  filter(!(ID == 'ID202' & Time == '72.4'))
data_Cu <- data_Cu %>%  filter(!(ID == 'ID277' & Time == '72.4'))

summary(data_Cu$Cu_coef)
sd(data_Cu$Cu_coef)
#D Cu:Ca ratios ranged from 0.25x10-4 to 0.023 (0.002 ± 0.003)

#Zn_coef
data_Zn <- data
data_Zn <- data_Zn %>%  filter(!(ID == 'ID386' & Time == '72.4'))
data_Zn <- data_Zn %>%  filter(!(ID == 'ID079' & Time == '72.4'))
data_Zn <- data_Zn %>%  filter(!(ID == 'ID509' & Time == '83.8'))
data_Zn <- data_Zn %>%  filter(!(ID == 'ID511' & Time == '83.8'))
data_Zn <- data_Zn %>%  filter(!(ID == 'ID508' & Time == '47.5'))

summary(data_Zn$Zn_coef)
sd(data_Zn$Zn_coef)
#D Zn:Ca ratios ranged from 2.20x10-7 to 0.04 (0.001 ± 0.003)

#Sr_coef
summary(data$Sr_coef)
sd(data$Sr_coef)
#D Sr:Ca ratios ranged from 0.022 to 0.097 (0.051 ± 0.017)

#Ba_coef
summary(data$Ba_coef)
sd(data$Ba_coef)
#D Ba:Ca ratios ranged from 8.64x10-5 to 0.046 (0.005 ± 0.006)

###
# Figures
###

###
# Mg
###

##transform x10^-4
data_Mg$Mg_coef_5 <- data_Mg$Mg_coef*10000

#letters based on post hoc Dunn's test
dunn_Mg <- dunnTest(Mg_coef ~ Month, data = data_Mg, method = "bonferroni")$res
dunn_Mg
cld_Mg <- cldList(P.adj ~ Comparison, data=dunn_Mg)
cld_Mg
cld_Mg <- as.data.frame(cld_Mg)                    
cld_Mg <- cld_Mg %>%
  rename(Month = Group)
cld_Mg

letter_Mg <- group_by(data_Mg, Month) %>%
  summarise(mean=mean(Mg_coef_5), quant = quantile(Mg_coef_5, probs = 0.75)) %>%
  arrange(desc(mean))
letter_Mg <- as.data.frame(letter_Mg)                    
letter_Mg

letter_Mg <- merge(letter_Mg, cld_Mg, by = "Month")
letter_Mg

Mg<-ggplot(data_Mg, aes(y=Mg_coef_5, x=Month)) +
  annotate('rect', xmin=2.5, xmax=5.5, ymin=0.03, ymax=0.27, alpha=0.15, fill='#999999')+
  annotate('rect', xmin=8.5, xmax=11.5, ymin=0.03, ymax=0.27,alpha=0.15, fill='#999999')+
  geom_boxplot(alpha = 1,width=0.5) +
  theme_bw() +
  xlab("") + ylab(expression(italic(D)[Mg:Ca]~(x~10^-4)))+
  theme(legend.position="none",
        axis.text = element_text(size = 15),
        axis.text.x = element_text(angle = -45),
        axis.title = element_text(size = 15))+
  scale_y_continuous(limits = c(0.03, 0.27), breaks=seq(0.03, 0.27, by = 0.06))+
  geom_text(data = letter_Mg, aes(x = Month, y = quant, label = Letter), vjust=-3, size = 4, color='darkgray')
Mg

###
# Mn
###

#letters based on post hoc Dunn's test
dunn_Mn <- dunnTest(Mn_coef ~ Month, data = data_Mn, method = "bonferroni")$res
dunn_Mn
cld_Mn <- cldList(P.adj ~ Comparison, data=dunn_Mn)
cld_Mn
cld_Mn <- as.data.frame(cld_Mn)                    
cld_Mn <- cld_Mn %>%
  rename(Month = Group)
cld_Mn

letter_Mn <- group_by(data_Mn, Month) %>%
  summarise(mean=mean(Mn_coef), quant = quantile(Mn_coef, probs = 0.75)) %>%
  arrange(desc(mean))
letter_Mn <- as.data.frame(letter_Mn)                    
letter_Mn

letter_Mn <- merge(letter_Mn, cld_Mn, by = "Month")
letter_Mn

Mn<-ggplot(data, aes(y=Mn_coef, x=Month)) +
  annotate('rect', xmin=2.5, xmax=5.5, ymin=0, ymax=0.0020, alpha=0.15, fill='#999999')+
  annotate('rect', xmin=8.5, xmax=11.5, ymin=0, ymax=0.0020,alpha=0.15, fill='#999999')+
  geom_boxplot(alpha = 1,width=0.5) +
  theme_bw() +
  xlab("") + ylab(expression(italic(D)[Mn:Ca]))+
  theme(legend.position="none",
        axis.text = element_text(size = 15),
        axis.text.x = element_text(angle = -45),
        axis.title = element_text(size = 15))+
  geom_text(data = letter_Mn, aes(x = Month, y = quant, label = Letter), vjust=-3, size = 4, color='darkgray')
Mn

###
# Cu
###

#letters based on post hoc Dunn's test
dunn_Cu <- dunnTest(Cu_coef ~ Month, data = data_Cu, method = "bonferroni")$res
dunn_Cu
cld_Cu <- cldList(P.adj ~ Comparison, data=dunn_Cu)
cld_Cu
cld_Cu <- as.data.frame(cld_Cu)                    
cld_Cu <- cld_Cu %>%
  rename(Month = Group)
cld_Cu

letter_Cu <- group_by(data_Cu, Month) %>%
  summarise(mean=mean(Cu_coef), quant = quantile(Cu_coef, probs = 0.75)) %>%
  arrange(desc(mean))
letter_Cu <- as.data.frame(letter_Cu)                    
letter_Cu

letter_Cu <- merge(letter_Cu, cld_Cu, by = "Month")
letter_Cu

Cu<-ggplot(data_Cu, aes(y=Cu_coef, x=Month)) +
  annotate('rect', xmin=2.5, xmax=5.5, ymin=0, ymax=0.023, alpha=0.15, fill='#999999')+
  annotate('rect', xmin=8.5, xmax=11.5, ymin=0, ymax=0.023,alpha=0.15, fill='#999999')+
  geom_boxplot(alpha = 1,width=0.5) +
  theme_bw() +
  xlab("") + ylab(expression(italic(D)[Cu:Ca]))+
  theme(legend.position="none",
        axis.text = element_text(size = 15),
        axis.text.x = element_text(angle = -45),
        axis.title = element_text(size = 15))+
  geom_text(data = letter_Cu, aes(x = Month, y = quant, label = Letter), vjust=-3, size = 4, color='darkgray')
Cu

###
# Zn
###

#delete outlier
data_Zn <- data_Zn %>%  filter(!(ID == 'ID324' & Time == '74.7'))

#letters based on post hoc Dunn's test
dunn_Zn <- dunnTest(Zn_coef ~ Month, data = data_Zn, method = "bonferroni")$res
dunn_Zn
cld_Zn <- cldList(P.adj ~ Comparison, data=dunn_Zn)
cld_Zn
cld_Zn <- as.data.frame(cld_Zn)                    
cld_Zn <- cld_Zn %>%
  rename(Month = Group)
cld_Zn

letter_Zn <- group_by(data_Zn, Month) %>%
  summarise(mean=mean(Zn_coef), quant = quantile(Zn_coef, probs = 0.75)) %>%
  arrange(desc(mean))
letter_Zn <- as.data.frame(letter_Zn)                    
letter_Zn

letter_Zn <- merge(letter_Zn, cld_Zn, by = "Month")
letter_Zn

Zn<-ggplot(data_Zn, aes(y=Zn_coef, x=Month)) +
  annotate('rect', xmin=2.5, xmax=5.5, ymin=0, ymax=0.0125, alpha=0.15, fill='#999999')+
  annotate('rect', xmin=8.5, xmax=11.5, ymin=0, ymax=0.0125,alpha=0.15, fill='#999999')+
  geom_boxplot(alpha = 1,width=0.5) +
  theme_bw() +
  xlab("") + ylab(expression(italic(D)[Zn:Ca]))+
  theme(legend.position="none",
        axis.text = element_text(size = 15),
        axis.text.x = element_text(angle = -45),
        axis.title = element_text(size = 15))+
  geom_text(data = letter_Zn, aes(x = Month, y = quant, label = Letter), vjust=-3, size = 4, color='darkgray')
Zn

###
# Sr
###

data_Sr<- data

#letters based on post hoc Dunn's test
dunn_Sr <- dunnTest(Sr_coef ~ Month, data = data_Sr, method = "bonferroni")$res
dunn_Sr
cld_Sr <- cldList(P.adj ~ Comparison, data=dunn_Sr)
cld_Sr
cld_Sr <- as.data.frame(cld_Sr)                    
cld_Sr <- cld_Sr %>%
  rename(Month = Group)
cld_Sr

letter_Sr <- group_by(data_Sr, Month) %>%
  summarise(mean=mean(Sr_coef), quant = quantile(Sr_coef, probs = 0.75)) %>%
  arrange(desc(mean))
letter_Sr <- as.data.frame(letter_Sr)                    
letter_Sr

letter_Sr <- merge(letter_Sr, cld_Sr, by = "Month")
letter_Sr

Sr<-ggplot(data, aes(y=Sr_coef, x=Month)) +
  annotate('rect', xmin=2.5, xmax=5.5, ymin=0.02, ymax=0.12, alpha=0.15, fill='#999999')+
  annotate('rect', xmin=8.5, xmax=11.5, ymin=0.02, ymax=0.12,alpha=0.15, fill='#999999')+
  geom_boxplot(alpha = 1,width=0.5) +
  theme_bw() +
  xlab("") + ylab(expression(italic(D)[Sr:Ca]))+
  theme(legend.position="none",
        axis.text = element_text(size = 15),
        axis.text.x = element_text(angle = -45),
        axis.title = element_text(size = 15))+
  scale_y_continuous(limits = c(0.02, 0.12))+
  geom_text(data = letter_Sr, aes(x = Month, y = quant, label = Letter), vjust=-3, size = 4, color='darkgray')
Sr

###
# Ba
###

data_Ba<-data

#letters based on post hoc Dunn's test
dunn_Ba <- dunnTest(Ba_coef ~ Month, data = data_Ba, method = "bonferroni")$res
dunn_Ba
cld_Ba <- cldList(P.adj ~ Comparison, data=dunn_Ba)
cld_Ba
cld_Ba <- as.data.frame(cld_Ba)                    
cld_Ba <- cld_Ba %>%
  rename(Month = Group)
cld_Ba

letter_Ba <- group_by(data_Ba, Month) %>%
  summarise(mean=mean(Ba_coef), quant = quantile(Ba_coef, probs = 0.75)) %>%
  arrange(desc(mean))
letter_Ba <- as.data.frame(letter_Ba)                    
letter_Ba

letter_Ba <- merge(letter_Ba, cld_Ba, by = "Month")
letter_Ba

Ba<-ggplot(data, aes(y=Ba_coef, x=Month)) +
  annotate('rect', xmin=2.5, xmax=5.5, ymin=0, ymax=0.047, alpha=0.15, fill='#999999')+
  annotate('rect', xmin=8.5, xmax=11.5, ymin=0, ymax=0.047,alpha=0.15, fill='#999999')+geom_boxplot(alpha = 1,width=0.5) +
  theme_bw() +
  xlab("") + ylab(expression(italic(D)[Ba:Ca]))+
  theme(legend.position="none",
        axis.text = element_text(size = 15),
        axis.text.x = element_text(angle = -45),
        axis.title = element_text(size = 15))+
  geom_text(data = letter_Ba, aes(x = Month, y = quant, label = Letter), vjust=-3, size = 4, color='darkgray')
Ba

######
# Fig. 8. Partition coefficients of elemental ratios (DEl:Ca) monthly sampled at the São Mateus estuary – 
#Southwestern Atlantic. Shade and transparent areas indicate the distinct seasons: late dry (from July to September), 
#early wet (from October to December), late wet (from January to March), and early dry (from April to June). 
#Lettering indicates post hoc Dunn’s test results in which months are significantly different (p < 0.05), 
#considering DEl:Ca as the response variable, and the sampled month as the predictor. 
#The DEl:Ca was calculated considering El:Cawater with a one-month correction.
######

ggarrange(Mg,Mn,Cu,Zn,Sr,Ba, 
          labels = c("A","B","C","D","E","F"),
          ncol = 2, nrow = 3,
          legend = "none")
#1000 x 1000

















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

data <- mutate_at(data, vars(c('Mg24', 'Mn55', 'Cu65','Zn66', 'Sr87', 'Ba138', 
                               'TL', 'Fulton',
                               'MgCa','MnCa',	'CuCa',	'ZnCa',	'SrCa',	'BaCa',
                               'Temp', 'Sal','OrderS')), as.numeric)

data <- data[data$Month != c('JUL'), ]

data$Month <- factor(data$Month, levels = c('AUG','SEP','OCT','NOV','DEC','JAN','FEB','MAR','APR','MAY','JUN'))
data$Sliding <- factor(data$Sliding, levels = c('JUL','AUG','SEP','OCT','NOV','DEC','JAN','FEB','MAR','APR','MAY'))

length(unique(data$ID))
str(data)

################################
# 3.1. Correlation
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

ggpairs(matriz, columns = 1:4, ggplot2::aes(alpha = 0.3),
        lower = list(continuous = "smooth", alpha = 0.3)) +   theme_bw()


################################
#GLMMs
################################

####
# Read dataset
####

#Clean R environment 
rm(list = ls())

#Read data water
data <- read.csv2("data_m1_e3.csv")

data <- data[data$Month != c('JUL'), ]

data <- mutate_at(data, vars(c('Mg24', 'Mn55', 'Cu65','Zn66', 'Sr87', 'Ba138', 
                               'TL', 'Fulton',
                               'MgCa','MnCa',	'CuCa',	'ZnCa',	'SrCa',	'BaCa',
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
best_model_Mg24 <- get.models(dredge_Mg, delta==0)[[1]]
summary(best_model_Mg24)
Anova(best_model_Mg24)

#Model selection table 
#   (Intrc)   Fltns    Sals   Temps      TLs df  logLik   AICc delta weight
#10  0.4992 -0.2612                 0.1152           6 176.991 -341.7  0.00  0.327
#14  0.4685 -0.2706          0.0707 0.1162           7 177.857 -341.4  0.36  0.274
#2   0.5405 -0.2161                                  5 175.266 -340.3  1.37  0.165
#26  0.5080 -0.2626                 0.1191 -0.01614  7 177.012 -339.7  2.05  0.118
#12  0.4948 -0.2623 0.008683        0.1155           7 177.002 -339.6  2.07  0.116

####
## Mg - Selected model
####

optimal_model_Mg24 <- lmer(Mg24ts ~ Fultons + TLs + (1 | Month) + (1 | ID), 
                   data = data_Mg, REML = FALSE,
                   na.action =  na.fail)


summary(optimal_model_Mg24)
Anova(optimal_model_Mg24)

#Analysis of Deviance Table (Type II Wald chisquare tests)
#Response: Mg24ts
#         Chisq Df Pr(>Chisq)  
#Fultons 4.7552  1    0.02921 *
#TLs     0.1242  1    0.72451  


####
#  Mg - Model validation
####
plot(optimal_model_Mg24)

#"DHARMa"
simulationOutput <- simulateResiduals(fittedModel = optimal_model_Mg24, plot = TRUE, use.u = FALSE)
testDispersion(optimal_model_Mg24, type='DHARMa')
#dispersion = 0.99421, p-value = 0.96

#r2 values
r.squaredGLMM(optimal_model_Mg24)
r2_nakagawa(optimal_model_Mg24)
#Conditional R2: 0.668
#Marginal R2: 0.034

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
#Durbin-Watson test = 1.1386, p-value = 0.1189


####
# Mg - Figure
####

model_Mg24 <- lmer(Mg24 ~ Fulton + TL + (1 |Month) + (1 | ID), 
                   data = data_Mg, REML = FALSE,
                   na.action =  na.fail)
Anova(model_Mg24)

#Profiled Confidence Intervals
profile_ci <- ggpredict(model_Mg24, terms = c("Fulton"))

Mg_fulton <-plot(profile_ci) +
  labs(title = " ",
       x = expression(Fulton~index), y = expression(Mg:Ca[otolith]~(mu~mol~mol^-1)))+
  theme_bw()+
  theme(legend.position="none",
        title = element_text(size = 15),
        axis.text = element_text(size = 15),
        axis.title = element_text(size = 15))+
  scale_x_continuous(breaks = seq(0, 3, by = 0.2))
Mg_fulton











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
best_model_Mn55 <- get.models(dredge_Mn, delta==0)[[1]]
summary(best_model_Mn55)
Anova(best_model_Mn55)

#Model selection table 
#   (Intrc)   Fltns   MnCas     Sals  Temps      TLs df logLik   AICc delta weight
#17  0.4252                          0.1350           4 74.254 -140.4  0.00  0.334
#49  0.4724                          0.1559 -0.09194  5 74.843 -139.5  0.88  0.215
#25  0.4439                 -0.04013 0.1357           5 74.601 -139.0  1.37  0.168
#18  0.4026 0.08186                  0.1252           5 74.545 -138.9  1.48  0.159
#19  0.4221         0.01653          0.1349           5 74.295 -138.4  1.98  0.124

####
## Mn - Selected model
####

optimal_model_Mn55 <- lmer(Mn55ts ~ Temps + MnCas + (1 | Month) + (1 | ID), 
                          data = data_Mn, REML = FALSE,
                          na.action =  na.fail)

summary(optimal_model_Mn55)
Anova(optimal_model_Mn55)

#Analysis of Deviance Table (Type II Wald chisquare tests)
#Response: Mn55ts
#       Chisq Df Pr(>Chisq)   
#Temps 8.6437  1   0.003282 **
#MnCas 0.0834  1   0.772747

####
#  Mn - Model validation
####
plot(optimal_model_Mn55)

#"DHARMa"
simulationOutput <- simulateResiduals(fittedModel = optimal_model_Mn55, plot = TRUE, use.u = FALSE)
testDispersion(optimal_model_Mn55, type='DHARMa')
#dispersion = 1.0028, p-value = 0.992

#r2 values
r.squaredGLMM(optimal_model_Mn55)
r2_nakagawa(optimal_model_Mn55)
# R2m       R2c
#0.04662594 0.4405301

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
#Durbin-Watson test = 1.4035, p-value = 0.294


####
# Mn - Figure
####
model_Mn55 <- lmer(Mn55 ~ Temp + MnCa + (1 | Month) + (1 | ID), 
                   data = data_Mn, REML = FALSE,
                   na.action =  na.fail)
Anova(model_Mn55)

#Profiled Confidence Intervals
profile_ci <- ggpredict(model_Mn55, terms = c("Temp"))

Mn_temp<-plot(profile_ci) +
  labs(title = " ",
       x = expression(Temperature~(ºC)), y = expression(Mn:Ca[otolith]~(mu~mol~mol^-1)), color = "")+
  theme_bw()+
  theme(legend.position="none",
        title = element_text(size = 15),
        axis.text = element_text(size = 15),
        axis.title = element_text(size = 15))

Mn_temp










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
best_model_Cu65 <- get.models(dredge_Cu, delta==0)[[1]]
summary(best_model_Cu65)
Anova(best_model_Cu65)

#Model selection table 
#   (Intrc)   CuCas    Sals   Temps df  logLik   AICc delta weight
#1   0.4938                          4 114.271 -220.4  0.00  0.329
#2   0.4803 0.04519                  5 114.928 -219.7  0.75  0.226
#9   0.4757                 0.03631  5 114.685 -219.2  1.23  0.178
#10  0.4603 0.04774         0.03853  6 115.455 -218.6  1.77  0.136
#5   0.4840         0.01999          5 114.381 -218.6  1.84  0.131


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
#CuCas 1.6443  1     0.1997
#Temps 1.0749  1     0.2998


####
#  Cu - Model validation
####
plot(optimal_model_Cu65)

#"DHARMa"
simulationOutput <- simulateResiduals(fittedModel = optimal_model_Cu65, plot = TRUE, use.u = FALSE)
testDispersion(optimal_model_Cu65, type='DHARMa')
#dispersion = 1.0093, p-value = 0.92

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
#Durbin-Watson test = 1.6108, p-value = 0.4995

####
# Cu - No figure, because did not have significant effects
####













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
best_model_Zn66 <- get.models(dredge_Zn, delta==0)[[1]]
summary(best_model_Zn66)
Anova(best_model_Zn66)

#Model selection table 
#   (Intrc)    Fltns     Sals   Temps   ZnCas df  logLik   AICc delta weight
#1   0.6382                                    4 182.607 -357.1  0.00  0.303
#17  0.6274                           0.05383  5 183.428 -356.7  0.42  0.245
#3   0.6559          -0.03653                  5 183.109 -356.0  1.06  0.178
#2   0.6553 -0.05139                           5 182.915 -355.6  1.45  0.147
#5   0.6284                   0.02011          5 182.768 -355.3  1.74  0.127
#21  0.6144                   0.02555 0.05687  6 183.740 -355.2  1.87  0.049

####
## Zn - Selected model
####

optimal_model_Zn66 <- lmer(Zn66ts ~ ZnCas + Temps + (1 | Month) + (1 | ID), 
                   data = data_Zn, REML = FALSE,
                   na.action =  na.fail)


summary(optimal_model_Zn66)
Anova(optimal_model_Zn66)

#Analysis of Deviance Table (Type II Wald chisquare tests)
#Response: Mg24ts
#         Chisq Df Pr(>Chisq)  
#ZnCas 2.1444  1     0.1431
#Temps 0.6595  1     0.4167

####
#  Zn - Model validation
####
plot(optimal_model_Zn66)

#"DHARMa"
simulationOutput <- simulateResiduals(fittedModel = optimal_model_Zn66, plot = TRUE, use.u = FALSE)
testDispersion(optimal_model_Zn66, type='DHARMa')
#dispersion = 1.001, p-value = 0.936

#r2 values
r.squaredGLMM(optimal_model_Zn66)
r2_nakagawa(optimal_model_Zn66)
#Conditional R2: 0.245
#Marginal R2: 0.014

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
#Durbin-Watson test = 1.3796, p-value = 0.2741


####
# Zn - No figure, because did not have significant effects
####






################################
#3.6. Sr
################################

#No outliers

#BoxCox Transformation
box = lm(Sr87 ~ Sal + Temp + TL + Fulton + BaCa, data = data)

bc<-boxcox(box, plotit = TRUE)
lambda <- bc$x[which.max(bc$y)]

if (lambda == 0) {
  data$Sr87t <- log(data$Sr87)
} else {
  data$Sr87t <- (data$Sr87^lambda - 1) / lambda
}

#Scale the predictors
data$Sr87ts <- rescale(data$Sr87t)
data$Sr87s <- rescale(data$Sr87)
data$Sals <- rescale(data$Sal)
data$Temps <- rescale(data$Temp)
data$TLs <- rescale(data$TL)
data$Fultons <- rescale(data$Fulton)
data$SrCas <- rescale(data$SrCa)


####
## Sr - Dredge
####

#Gaussian
full_model_Sr <- lmer(Sr87s ~ Sals + Temps + TLs + Fultons + SrCas + (1 |Month) + (1 | ID), 
                       data = data, REML = FALSE,
                       na.action =  na.fail)


dredge_Sr <- dredge(full_model_Sr)
dredge_Sr
subset(dredge_Sr, delta==0)
dredge_Sr[1:5,]
best_model_Sr87 <- get.models(dredge_Sr, delta==0)[[1]]
summary(best_model_Sr87)
Anova(best_model_Sr87)

#Model selection table 
#   (Intrc)    Fltns     Sals    SrCas   Temps     TLs df  logLik   AICc delta weight
#9   0.5831                            -0.1275          5 178.208 -346.2  0.00  0.406
#25  0.5699                            -0.1333 0.02567  6 178.285 -344.3  1.92  0.155
#11  0.5876          -0.00978          -0.1274          6 178.243 -344.2  2.00  0.149
#10  0.5805 0.009194                   -0.1286          6 178.214 -344.2  2.06  0.145
#13  0.5816                   0.004707 -0.1304          6 178.214 -344.2  2.06  0.145

####
## Sr - Selected model
####

optimal_model_Sr87 <- lmer(Sr87s ~ Temps + TLs + (1 | Month) + (1 | ID), 
                   data = data, REML = FALSE,
                   na.action =  na.fail)

summary(optimal_model_Sr87)
Anova(optimal_model_Sr87)

#Analysis of Deviance Table (Type II Wald chisquare tests)
#Response: Mg24ts
#         Chisq Df Pr(>Chisq)  
#Temps 12.2246  1  0.0004716 ***
#TLs    0.1552  1  0.6936503


####
#  Sr - Model validation
####
plot(optimal_model_Sr87)

#"DHARMa"
simulationOutput <- simulateResiduals(fittedModel = optimal_model_Sr87, plot = TRUE, use.u = FALSE)
testDispersion(optimal_model_Sr87, type='DHARMa')
#dispersion = 1.008, p-value = 0.904

#r2 values
r.squaredGLMM(optimal_model_Sr87)
r2_nakagawa(optimal_model_Sr87)
# R2m       R2c
#0.07232007 0.5029846

#check_distribution
check_distribution <- check_distribution(optimal_model_Sr87)
plot(check_distribution)
check_distribution(optimal_model_Sr87)

#Variance inflation factors (VIF)
check_collinearity(optimal_model_Sr87)
#Low Correlation

#Temporal autocorrelation
res = simulateResiduals(optimal_model_Sr87)
res = recalculateResiduals(res, group = data$OrderS)
testTemporalAutocorrelation(res, time = unique(data$OrderS))
#Durbin-Watson test = 2.633, p-value = 0.2639


####
# Sr - Figure
####

data$Sr87m <- data$Sr87/1000 #transform umol mol-1 to mmol mol-1

model_Sr87 <- lmer(Sr87m ~ Temp + TL + (1 | Month) + (1 | ID), 
                   data = data, REML = FALSE,
                   na.action =  na.fail)
Anova(model_Sr87)

#Profiled Confidence Intervals
profile_ci <- ggpredict(model_Sr87, terms = c("Temp"))


Sr_temp <- plot(profile_ci) +
  labs(title = " ",
       x = expression(Temperature~(ºC)), y = expression(Sr:Ca[otolith]~(mmol~mol^-1)), color = "")+
  theme_bw()+
  theme(legend.position="none",
        title = element_text(size = 15),
        axis.text = element_text(size = 15),
        axis.title = element_text(size = 15))

Sr_temp











################################
#3.7. Ba
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
best_model_Ba138 <- get.models(dredge_Ba, delta==0)[[1]]
summary(best_model_Ba138)
Anova(best_model_Ba138)

#Month 1 Edge 3
#Model selection table 
#(Intrc)  BaCas   Fltns     Sals  Temps    TLs df  logLik   AICc delta weight
#26  0.3036 0.1734                  0.2167 0.1744  7 337.375 -660.4  0.00  0.367
#25  0.3356                         0.2306 0.1741  6 335.825 -659.4  1.01  0.221
#28  0.2930 0.1755 0.03486          0.2117 0.1763  8 337.471 -658.5  1.91  0.141
#18  0.3966 0.1918                         0.1897  6 335.363 -658.5  1.94  0.140
#30  0.3125 0.1700         -0.01692 0.2169 0.1744  8 337.389 -658.3  2.07  0.130

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
#BaCas 3.6190  1   0.057122 . 
#Temps 5.0569  1   0.024529 * 
#TLs   8.2349  1   0.004109 **


####
#  Ba - Model validation
####
plot(optimal_model_Ba138)

#"DHARMa"
simulationOutput <- simulateResiduals(fittedModel = optimal_model_Ba138, plot = TRUE, use.u = FALSE)
testDispersion(optimal_model_Ba138, type='DHARMa')
#dispersion = 0.92601, p-value = 0.8

#r2 values
r.squaredGLMM(optimal_model_Ba138)
r2_nakagawa(optimal_model_Ba138)
#Conditional R2: 0.914
#Marginal R2: 0.358

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
#Durbin-Watson test = 1.6605, p-value = 0.5567




####
# Ba - Figure
####

model_Ba138 <- lmer(Ba138 ~ Temp + TL + (1 |Month) + (1 | ID), 
                    data = data, REML = FALSE,
                    na.action =  na.fail)

#Profiled Confidence Intervals

profile_ci <- ggpredict(model_Ba138, terms = c("Temp"))

Ba_temp <-plot(profile_ci) +
  labs(title = " ",
       x = expression(Temperature~(ºC)), y =  expression(Ba:Ca[otolith]~(mu~mol~mol^-1)), color = "")+
  theme_bw()+
  theme(legend.position="none",
        title = element_text(size = 15),
        axis.text = element_text(size = 15),
        axis.title = element_text(size = 15))
Ba_temp

profile_ci <- ggpredict(model_Ba138, terms = c("TL"))

Ba_TL <-plot(profile_ci) +
  labs(title = " ",
       x = expression(Total~length~(mm)), y =  expression(Ba:Ca[otolith]~(mu~mol~mol^-1)), color = "")+
  theme_bw()+
  theme(legend.position="none",
        title = element_text(size = 15),
        axis.text = element_text(size = 15),
        axis.title = element_text(size = 15))

Ba_TL









################################
# 3.8 Figure 9
################################

#############
#Fig. 9. The significant relationships between otolith elemental signatures (El:Caotolith) 
#and environmental and biological factors based on our generalized linear mixed models.
#############

library(ggpubr)

Mg_fulton
Mn_temp
Sr_temp
Ba_temp
Ba_TL

ggarrange(Mg_fulton, Mn_temp, Sr_temp, Ba_temp, Ba_TL,
          labels = c("A","B","C","D","E"),
          ncol = 3, nrow = 2,
          legend = "none")
#1400 x 800

#END