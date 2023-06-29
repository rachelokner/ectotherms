#STAT 4600 R code
rm(list = ls())

datafolder <-  "C:/Users/rache/Downloads/STAT 4600/"
# load libraries
library(ape)
library(car)
library(caper)
library(multcomp)
library(geiger)
library(nlme)
library(phytools)
library(castor)
citation("phytools")
citation("geiger")

# read in all the data
alldata <- data.frame(read.csv(paste(datafolder, "FinalData_wEndos.csv", sep = ""), header = T))

#read in the phylogenetic tree
alltree <- read.tree(paste(datafolder, "allsupertreeBL.tre", sep = "")) 
#make a comparative object
all.data <- comparative.data(phy = alltree, data = alldata, names.col = SpeciesNameInTree, vcv = T, na.omit = F, warn.dropped = T)

#subset to ectotherms
ecto_data <- alldata[alldata$Therm=="Ecto",]
species_ecto<-as.vector(alldata$SpeciesNameInTree[c(which(alldata$Therm == "Ecto"))])
ecto_tree<-keep.tip(alltree, species_ecto)
ecto_tree<-multi2di(ecto_tree)

#ensuring dataset matches the species name in the tree
rownames(ecto_data) <- ecto_data$SpeciesNameInTree
ectodata <- ecto_data[match(ecto_tree$tip.label,rownames(ecto_data)),]

#creating comparative object
ecto.data <- comparative.data(phy = alltree, data = ectodata, names.col = SpeciesNameInTree, vcv = T, na.omit = F, warn.dropped = F)
ecto.data$data<-droplevels(ecto.data$data)


#trim to endotherms (used for categorical data)
endodata<-subset(alldata,alldata$Therm=="Endo")
endo.data <- comparative.data(phy = alltree, data = endodata, names.col = SpeciesNameInTree, vcv = T, na.omit = F, warn.dropped = F)
endo.data$data<-droplevels(endo.data$data)

#ANALYSES#

#descriptive statistics for longevity
summary(log(ectodata$quant95MortAge))
sd(log(ectodata$quant95MortAge))

#descriptive statistics for rate of aging
summary(ectodata$GompertzSlope)
sd(ectodata$GompertzSlope)

#################################################################################
#1. BODY SIZE HYPOTHESES

#a) Body Size vs Longevity

#plotting graph of log Body size vs log Longevity
plot(log(ectodata$Body_mass), log(ectodata$quant95MortAge), pch =c(15, 16, 17, 18, 19, 20)[factor(ectodata$Type)], col = c(4,5,3,1,2,6)[factor(ectodata$Type)],
     cex=1.4, xlab="Log Body Mass", ylab="Log Longevity", main="Body Mass vs Longevity in Ectotherms")
legend(0, 4.9, legend=c("OLS", "PGLS"),
       col=c("red", "blue"), lty=1:1:1, lwd=2:2:2, cex=0.8)
legend("bottomright",
       legend = levels(factor(ectodata$Type)),
       pch =c(15, 16, 17, 18, 19, 20), cex=0.7,
      col= c(4,5,3,1,2,6)[factor(levels(factor(ectodata$Type)))])

#descriptive statistics of Body Mass variable
summary(log(ectodata$Body_mass))
sd(log(ectodata$Body_mass))

#calculating phylogenetic signal for Body Mass variable
phylosig(ecto_tree, log(ectodata$Body_mass), method="lambda", test=T, nsim=1000, se=NULL, start=NULL,
         control=list())

names(ectodata$quant95MortAge)<-row.names(ectodata)
phylosig(ecto_tree, log(ectodata$quant95MortAge), method="lambda", test=T, nsim=1000, se=NULL, start=NULL,
         control=list())

#Method 1: OLS regression

model1<-lm(log(ectodata$quant95MortAge)~log(ectodata$Body_mass))
summary(model1)
#superimposing regression line onto plot
abline(model1$coefficients, col="red", lwd=2)
AIC(model1) #<-238.2476

#Method 2: PGLS Regression

model2<-pgls(log(ectodata$quant95MortAge)~log(ectodata$Body_mass), data=ecto.data, lambda = "ML")
summary(model2)#superimposing regression line onto plot
abline(model2$model$coef, col="blue", lwd=2)
AIC(model2) #<-223.3847
#_____________________________________________________

#Method 3: PICs

b.m<-log(ectodata$Body_mass)
longev<-log(ectodata$quant95MortAge)
names(b.m) <- row.names(ectodata)
names(longev) <- row.names(ectodata)

#calculating standardized PICs for x and y
x1 <- pic(b.m, ecto_tree, var.contrasts=T, scaled = T)
y1 <- pic(longev, ecto_tree, var.contrasts=T, scaled = T)
#linear regression through the origin (-1)
model3 <- lm(y1[,1]~x1[,1]-1)
summary.lm(model3)

plot(x1[,1], y1[,1], xlab="Log Body Mass Contrasts", ylab="Log Longevity Contrasts", main="Body Mass vs Longevity PICs",
     pch=18, col = c("grey35"))
legend("topleft", legend=c("PIC Regression Line"),
       col=c("green"), lty=1, lwd=2, cex=0.9)

#superimposing regression line onto the plot
abline(model3, col="green", lwd=2)
AIC(model3) #<-682.4394

#standardized residuals boxplot
model1_residuals<-model1$residuals/sd(model1$residuals)
model2_residuals<-model2$residuals/sd(model2$residuals)
model3_residuals<-model3$residuals/sd(model3$residuals)

boxplot(model1_residuals, model2_residuals, model3_residuals, names=c("OLS", "PGLS", "PIC"), 
        main="1.a) Body Mass vs Longevity")

#-----------------------------------------------------------------------------------------------------

#1.b. Body Size vs Rate of Aging
plot(log(ectodata$Body_mass), ectodata$GompertzSlope,pch =c(15, 16, 17, 18, 19, 20)[factor(ectodata$Type)], col = c(4,5,3,1,2,6)[factor(ectodata$Type)],
     cex=1.4, xlab="Log Body Mass", ylab="Rate of Aging", main="Body Mass vs Rate of Aging in Ectotherms")
legend(8.2, 2.1, legend=c("OLS", "PGLS"),
       col=c("red", "blue"), lty=1:1:1, lwd=2:2:2, cex=0.8)
legend(7.7, 1.5, cex=0.8,
       legend = levels(factor(ectodata$Type)),
       pch =c(15, 16, 17, 18, 19, 20),
       col= c(4,5,3,1,2,6)[factor(levels(factor(ectodata$Type)))])

#Method 1:
model1b<-lm(ectodata$GompertzSlope~log(ectodata$Body_mass))
summary(model1b)

abline(model1b$coefficients, col="red", lwd=2)
AIC(model1b) #<-132.9566

#Method 2: PGLS regression
model2b<-pgls(ectodata$GompertzSlope~log(ectodata$Body_mass), data=ecto.data, lambda = "ML")
summary(model2b)
abline(model2b$model$coef, col="blue", lwd=2)
AIC(model2b) #<-127.664

#Method 3: PICs

names(ectodata$GompertzSlope) <- row.names(ectodata)
#Calculating PICs for x and y
x1b <- pic(b.m, ecto_tree, var.contrasts=T, scaled = T)
y1b <- pic(ectodata$GompertzSlope, ecto_tree, var.contrasts=T, scaled = T)

model3b <- lm(y1b[,1]~ x1b[,1] -1)    # the "-1" forces line through origin
summary.lm(model3b)

plot(x1b[,1], y1b[,1], xlab="Log Body Mass Contrasts", ylab="Rate of Aging Contrasts", main="Body Mass vs Rate of Aging PICs",
     pch=18, col = c("grey35"))
legend("topleft", legend=c("PIC Regression Line"),
       col=c("green"), lty=1, lwd=2, cex=0.9)

abline(model3b, col="green", lwd=2)
AIC(model3b) #<- 590.803

#plotting standardized residuals for each model
mod1b_residual<-model1b$residuals/sd(model1b$residuals)
mod2b_residual<-model2b$residuals/sd(model2b$residuals)
mod3b_residual<-model3b$residuals/sd(model3b$residuals)

boxplot(mod1b_residual, mod2b_residual, mod3b_residual, 
names=c("OLS", "PGLS", "PIC"), main="1.b) Body Mass vs Rate of Aging")
#########################################################################################


#########################################################################################
#ANNUAL FECUNDITY HYPOTHESES

# 2.a) Annual Fecundity vs Longevity
AF<-log(ectodata$AnnualFecundity)
AF<-na.omit(AF)

# descriptive statistics for annual fecundity
summary(AF) 
sd(AF)

which(is.na(log(ectodata$AnnualFecundity)))
longevity<-log(ectodata$quant95MortAge)[-c(73, 98, 107)]

#plotting log annual fecundity vs log longevity
plot(AF, longevity,pch =c(15, 16, 17, 18, 19, 20)[factor(ectodata$Type)], col = c(4,5,3,1,2,6)[factor(ectodata$Type)],
     cex=1.4, xlab="Log Annual Fecundity", ylab="Log Longevity", main="Annual Fecundity vs Longevity in Ectotherms")
legend(4, 4.9, legend=c("OLS", "PGLS"),
       col=c("red", "blue", "green"), lty=1:1:1, lwd=2:2:2, cex=0.8)
legend("topright", cex=0.8,
       legend = levels(factor(ectodata$Type)),
       pch =c(15, 16, 17, 18, 19, 20),
       col= c(4,5,3,1,2,6)[factor(levels(factor(ectodata$Type)))])

#Method 1: OLS regression

AF_model1<-lm(longevity~AF)
summary(AF_model1)
#superimposing regression line onto plot
abline(AF_model1$coefficients, col="red", lwd=2)
AIC(AF_model1) #<-[1] 260.2278


#Method 2:PGLS regression
AF_model2<-pgls(log(ectodata$quant95MortAge)~log(ectodata$AnnualFecundity), data=ecto.data, lambda = "ML")
summary(AF_model2)
abline(AF_model2$model$coef, col="blue", lwd=2)
AIC(AF_model2) #<- 236.1343

#Method 3: PICs

ecto_tree_AF<-drop.tip(ecto_tree, c("Eurycea_tonkawae", "Lyciasalamandra_fazilae", "Rana_tavasensis"))

#calculating standardized PICs for x and y
x2 <- pic(AF, ecto_tree_AF, var.contrasts=T, scaled=T)
y2 <- pic(longevity, ecto_tree_AF,var.contrasts=T, scaled=T)

AF_model3 <- lm(y2[,1]~x2[,1]-1)    # the "-1" forces line through origin
summary.lm(AF_model3)
AIC(AF_model3)#<- 683.9009

plot(x2[,1], y2[,1], xlab="Log Annual Fecundity Contrasts", ylab="Log Longevity Contrasts", main="PICs for Annual Fecundity and Longevity",
     pch=18, col = c("grey35"))
legend("topright", legend=c("PIC Regression Line"),
       col=c("green"), lty=1, lwd=2, cex=0.9)

abline(AF_model3, col="green", lwd=2)

#plotting standardized residuals for each model
AFmod1_res<- AF_model1$residuals/sd(AF_model1$residuals)
AFmod2_res<- AF_model2$residuals/sd(AF_model2$residuals)
AFmod3_res<- AF_model3$residuals/sd(AF_model3$residuals)

boxplot(AFmod1_res, AFmod2_res, AFmod3_res, 
        names=c("OLS", "PGLS", "PIC"), main="2.a) Annual Fecundity vs Longevity")

#------------------------------------------------------------------------------------
#2.b) Annual Fecundity vs Rate of Aging

aging<-ectodata$GompertzSlope[-c(73, 98, 107)]

plot(AF, aging, pch =c(15, 16, 17, 18, 19, 20)[factor(ectodata$Type[-c(73, 98, 107)])], col = c(4,5,3,1,2,6)[factor(ectodata$Type[-c(73, 98, 107)])],
     cex=1.4, xlab="Log Annual Fecundity", ylab="Rate of Aging", main="Annual Fecundity vs Rate of Aging in Ectotherms")
legend(4.2, 2.1, legend=c("OLS", "PGLS"),
       col=c("red", "blue"), lty=1:1:1, lwd=2:2:2, cex=0.8)
legend("topright", cex=0.6,
       legend = levels(factor(ectodata$Type)[-c(73, 98, 107)]),
       pch =c(15, 16, 17, 18, 19, 20),
       col= c(4,5,3,1,2,6)[factor(levels(factor(ectodata$Type)[-c(73, 98, 107)]))])

#Method 1: OLS regression

AF_model1b<-lm(aging~AF)
summary(AF_model1b)
#superimposing regression line onto plot
abline(AF_model1b$coefficients, col="red", lwd=2)
AIC(AF_model1b) #<-[1] 132.8191


#Method 2:PGLS regression
AF_model2b<-pgls(ectodata$GompertzSlope~log(ectodata$AnnualFecundity), data=ecto.data, lambda = "ML")
summary(AF_model2b)
abline(AF_model2b$model$coef, col="blue", lwd=2)
AIC(AF_model2b) #<- 125.5458

#Method 3: PICs

ecto_tree_AF<-drop.tip(ecto_tree, c("Eurycea_tonkawae", "Lyciasalamandra_fazilae", "Rana_tavasensis"))

x2b <- pic(AF, ecto_tree_AF, var.contrasts=T, scaled=T)
y2b <- pic(aging, ecto_tree_AF,var.contrasts=T, scaled=T)

AF_model3b <- lm(y2b[,1]~x2b[,1]-1)    # the "-1" forces line through origin
summary.lm(AF_model3b)

plot(x2b[,1], y2b[,1], xlab="Log Annual Fecundity Contrasts", ylab="Rate of Aging Contrasts", main="PICs for Annual Fecundity and Rate of Aging",
     pch=18, col = c("grey35"))
legend("topleft", legend=c("PIC Regression Line"),
       col=c("green"), lty=1, lwd=2, cex=0.9)

abline(AF_model3b, col="green", lwd=2)
AIC(AF_model3b) #<-571.9254

#standardized residuals for each model
AFmod1b_res<- AF_model1b$residuals/sd(AF_model1b$residuals)
AFmod2b_res<- AF_model2b$residuals/sd(AF_model2b$residuals)
AFmod3b_res<- AF_model3b$residuals/sd(AF_model3b$residuals)


boxplot(AFmod1b_res, AFmod2b_res , AFmod3b_res,
        names=c("OLS", "PGLS", "PIC"), main="2.b) Annual Fecundity vs Rate of Aging")

################################################################################################
#3.a) AGE OF FIRST REPRODUCTION vs LONGEVITY
repro<-log(ectodata$Age.at.First.Repro)
repro<-na.omit(repro)

#descriptive statistics for age at first reproduction
summary(repro)
sd(repro)

which(is.na(log(ectodata$Age.at.First.Repro)))
longev_repro<-log(ectodata$quant95MortAge)[-73]

plot(repro, longev_repro, pch =c(15, 16, 17, 18, 19, 20)[factor(ectodata$Type)], col = c(4,5,3,1,2,6)[factor(ectodata$Type)],
     cex=1.4, xlab="Log Age at First Repro", ylab="Log Longevity", main="Age at First Reproduction vs Longevity in Ectotherms")
legend(-0.7, 4.9, legend=c("OLS", "PGLS"),
       col=c("red", "blue"), lty=1:1:1, lwd=2:2:2, cex=0.7)
legend(0.2, 4.9, cex=0.6,
       legend = levels(factor(ectodata$Type)),
       pch =c(15, 16, 17, 18, 19, 20),
       col= c(4,5,3,1,2,6)[factor(levels(factor(ectodata$Type)))])

#Method 1: OLS Regression
repromod1<-lm(longev_repro~repro)
summary(repromod1)
abline(repromod1, col="red", lwd=2)
AIC(repromod1) #<-206.3519

#Method 2: PGLS Regression
repromod2<-pgls(log(ectodata$quant95MortAge)~log(ectodata$Age.at.First.Repro), data=ecto.data, lambda = "ML")
summary(repromod2)
abline(repromod2, col="blue", lwd=2)
AIC(repromod2) #<-203.1525

#Method 3: PICS
ecto_tree_REPRO<-drop.tip(ecto_tree, c("Rana_tavasensis"))

#calculating standardized PICs for x and y
x3<-pic(repro, ecto_tree_REPRO, var.contrasts = T, scaled = T)
y3<-pic(longev_repro, ecto_tree_REPRO, var.contrasts = T, scaled = T)

repromod3<-lm(y3[,1]~x3[,1]-1)
summary(repromod3)

plot(x3[,1], y3[,1], xlab="Log Age at First Repro Contrasts", ylab="Log Longevity Contrasts", main="PICs for Age at First Repro and Longevity",
     pch=18, col = c("grey35"))
legend("topleft", legend=c("PIC Regression Line"),
       col=c("green"), lty=1, lwd=2, cex=0.9)

abline(repromod3, col="green", lwd=2)
AIC(repromod3) #<-662.4107

#standardized residuals for each model
repromod1_res<- repromod1$residuals/sd(repromod1$residuals)
repromod2_res<- repromod2$residuals/sd(repromod2$residuals)
repromod3_res<- repromod3$residuals/sd(repromod3$residuals)


boxplot(repromod1_res, repromod2_res , repromod3_res,
        names=c("OLS", "PGLS", "PIC"), main="3.a) Age at First Repro vs Longevity")

#---------------------------------------------------------------------------
#3.b) AGE AT FIRST REPRODUCTION vs RATE OF AGING

aging_repro<-ectodata$GompertzSlope[-73]

plot(repro, aging_repro, pch =c(15, 16, 17, 18, 19, 20)[factor(ectodata$Type)], col = c(4,5,3,1,2,6)[factor(ectodata$Type)],
     cex=1.4, xlab="Log Age at First Repro", ylab="Rate of Aging", main="Age at First Reproduction vs Rate of Aging in Ectotherms")
legend("topleft", legend=c("OLS", "PGLS"),
       col=c("red", "blue"), lty=1:1:1, lwd=2:2:2, cex=0.8)
legend("topright", cex=0.7,
       legend = levels(factor(ectodata$Type)),
       pch =c(15, 16, 17, 18, 19, 20),
       col= c(4,5,3,1,2,6)[factor(levels(factor(ectodata$Type)))])

#Method 1: OLS Regression

repromod1a<-lm(aging_repro~repro)
summary(repromod1a)
abline(repromod1a, col="red", lwd=2)
AIC(repromod1a) #<-118.9694

#Method 2: PGLS Regression
repromod2a<-pgls(ectodata$GompertzSlope~log(ectodata$Age.at.First.Repro), data=ecto.data, lambda = "ML")
summary(repromod2a)
abline(repromod2a, col="blue", lwd=2)
AIC(repromod2a) #<-116.8207


#Method 3: PICS
ecto_tree_REPRO<-drop.tip(ecto_tree, c("Rana_tavasensis"))
x3b<-pic(repro, ecto_tree_REPRO, var.contrasts = T, scaled = T)
y3b<-pic(aging_repro, ecto_tree_REPRO, var.contrasts = T, scaled = T)

repromod3a<-lm(y3b[,1]~x3b[,1]-1)
summary(repromod3a)

plot(x3b[,1], y3b[,1], xlab="Log Age at First Repro Contrasts", ylab="Rate of Aging Contrasts", main="PICs for Age at First Repro and Rate of Aging",
     pch=18, col = c("grey35"))
legend("topright", legend=c("PIC Regression Line"),
       col=c("green"), lty=1, lwd=2, cex=0.9)

abline(repromod3a, col="green", lwd=2)
AIC(repromod3a) #<-570.6765

#standardized residuals for each model
repromod1a_res<- repromod1a$residuals/sd(repromod1a$residuals)
repromod2a_res<- repromod2a$residuals/sd(repromod2a$residuals)
repromod3a_res<- repromod3a$residuals/sd(repromod3a$residuals)


boxplot(repromod1a_res, repromod2a_res , repromod3a_res,
        names=c("OLS", "PGLS", "PIC"), main="3.b) Age at First Repro vs Rate of Aging")
par(mfrow=c(2,2))
###################################################################################################################



#########################################################################################
#EXTRA- Leverage and Influence of PICs lines

#1.a)Body Mass and Longevity

leverage1a<-round(hatvalues(model3$fitted.values), digits=5)
max1a<-max(leverage1a) #<-0.42549

plot(model3b)

plot(x1b, y1b)
outliers1b <- hatvalues(model3b) > 3 * mean(hatvalues(model3b))



# #CATEGORICAL VARIABLES
# #################################################################################################################
# #ENDOTHERMS vs ECTOTHERMS HYPOTHESES
# 
# #3.a) Longevity in ectotherms vs endotherms
# 
# new_order <- with(alldata, reorder(alldata$Type , log(alldata$quant95MortAge), median , na.rm=T))
# boxplot(log(alldata$quant95MortAge)~new_order, xlab="Group", ylab="Log Longevity",
#      main="Ectotherm and Endotherm Longevity", col=c("#69b3a2","#69b3a2","#69b3a2", "slateblue1","slateblue1", "#69b3a2", "#69b3a2", "#69b3a2"))
# 
# #Method 1: ANOVA
# EEmodel1<-aov(log(alldata$quant95MortAge) ~ alldata$Type)
# anova(EEmodel1)
# summary(EEmodel1)
# 
# plot(TukeyHSD(EEmodel1), las=2, cex=0.3)
# text(cex=0.3)
# AIC(EEmodel1)
# 
# #Method 2: Phylogenetic ANOVA 
# EEmodel2 <- pgls(log(alldata$quant95MortAge) ~ factor(alldata$Type), data = all.data,lambda = "ML") 
# anova(EEmodel2)
# 
# summary(EEmodel2)
# #Lambda = 0
# AIC(EEmodel2)
# 
# 
# #Vector of names
# speciesnames<-alldata$SpeciesNameInTree
# 
# longevity<-log(alldata$quant95MortAge)
# longevity_matrix<-as.matrix(longevity)
# group<-as.factor(c(alldata$Type))
# group.vector<-as.vector(group)
# 
# names(group)<-alltree$tip.label
# names(longevity_matrix)<-alltree$tip.label
# names(group)=rownames(longevity_matrix)
# 
# aov.model1<-phylANOVA(alltree, group, longevity_matrix, posthoc = T)
# (aov.model1)
# 
# #______________________________________________________________________________
# #3.b) Rate of Aging in ectotherms vs endotherms
# new_order1 <- with(alldata, reorder(alldata$Type , alldata$GompertzSlope), median , na.rm=T)
# boxplot(alldata$GompertzSlope~new_order1, xlab="Group", ylab="Rate of Aging",
#         main="Ectotherm and Endotherm Aging", col=c("#69b3a2","#69b3a2","slateblue1","#69b3a2", "#69b3a2", "slateblue1", "#69b3a2", "#69b3a2"))
# 
# #Method 1: ANOVA
# EEmodel2a<-aov(alldata$GompertzSlope ~ alldata$Type)
# anova(EEmodel2a)
# summary.aov(EEmodel2a)
# 
# plot(TukeyHSD(EEmodel2a), las=2, cex=0.9)
# text(cex=0.3)
# AIC(EEmodel1)
# 
# #Method 2: Phylogenetic ANOVA
# EEmodel2b<-pgls(alldata$GompertzSlope ~ factor(alldata$Type), data = all.data,lambda = "ML")
# anova(EEmodel2b)
# 
# 
# all_aging<-as.vector(alldata$GompertzSlope)
# names(group)<-alltree$tip.label
# names(all_aging)<-alltree$tip.label
# 
# aov.model2<-phylANOVA(alltree, group, all_aging, posthoc = T)
# 
# #########################################################################################
# #package for nice ANOVA tables
# library(knitr)
# 
# #4.a) protection phenotype vs LONGEVITY
# new_order <- with(ectodata, reorder(ectodata$Protection , log(ectodata$quant95MortAge)), median , na.rm=T)
# boxplot(log(ectodata$quant95MortAge)~new_order, xlab="Protection Type", ylab="Log Longevity", las=1,
#         main="Protection Phenotypes vs Longevity in Ectotherms", col=c("gold","#69b3a2","#69b3a2", "#69b3a2", "#69b3a2", "slateblue1", "slateblue1", "slateblue1"))
# legend("bottomright", cex=1,
#        legend = levels(factor(ectodata$ProtectionCat)), pch=15,
#        col=c("#69b3a2", "gold", "slateblue1")[factor(levels(factor(ectodata$ProtectionCat)))])
# 
# #Method 1: ANOVA
# PRmodel1<-aov(log(ectodata$quant95MortAge) ~ ectodata$Protection)
# summary(PRmodel1)
# anova(PRmodel1)
# 
# #Method 2: Phylogenetic ANOVA
# phenotype<-ectodata$Protection
# names(phenotype)<-ecto_tree$tip.label
# names(longev)<-ecto_tree$tip.label
# 
# phenotypemod<-phylANOVA(ecto_tree, phenotype, longev)
# (phenotypemod)
# 
# #_____________________________________________________________________________________________
# #protection phenotype (ECTOTHERMS) vs AGING
# plot(factor(ectodata$ProtectionCat), ectodata$GompertzSlope, xlab="Protection Type", ylab="Rate of Aging", main="Protection Phenotypes vs Rate of Aging")
# 
# #Method 1: Regular ANOVA
# PRmodel1b<-aov(ectodata$GompertzSlope ~ ectodata$ProtectionCat)
# summary(PRmodel1b)
# anova(PRmodel1b)
# a1b<-aov(ectodata$GompertzSlope ~ ectodata$ProtectionCat)
# Tukey1b <- TukeyHSD(a1b, conf.level = 0.95)
# plot(Tukey1b)
# AIC(PRmodel1b) # <-122.6843
# 
# PRmodel2 <- pgls(ectodata$GompertzSlope ~ ectodata$ProtectionCat, data = ecto.data, lambda = "ML")
# summary(PRmodel2)
# anova(PRmodel2)
# 
# AIC(PRmodel2) #<-120.6585
# long<-as.vector(ectodata$GompertzSlope)
# protection<-as.factor(ectodata$ProtectionCat)
# 
# aov.phylo(long~protection, phy = ecto_tree, data=ecto.data)
# 
# #Method 2: Phylogenetic ANOVA
# names(long)<-ecto_tree$tip.label
# names(protection)<-ecto_tree$tip.label
# 
# phylANOVA(ecto_tree, protection, long, posthoc=TRUE, p.adj="holm")
# 
# 


#CODE taken from 
#https://sarahleejane.github.io/learning/r/2014/09/21/plotting-data-points-on-maps-with-r.html

# library(maps)
# library(ggplot2)
# world_map <- map_data("world")
# 
# 
# p <- ggplot() + coord_fixed() +
#   xlab("") + ylab("")
# 
# #Add map to base plot
# base_world_messy <- p + geom_polygon(data=world_map, aes(x=long, y=lat, group=group), 
#                                      colour="light green", fill="light green")
# 
# base_world_messy
# 
# cleanup <- 
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
#         panel.background = element_rect(fill = 'white', colour = 'white'), 
#         axis.line = element_line(colour = "white"), legend.position="none",
#         axis.ticks=element_blank(), axis.text.x=element_blank(),
#         axis.text.y=element_blank())
# 
# base_world <- base_world_messy + cleanup
# 
# base_world
# 
# map_data_LONGEV<- 
#   base_world +
#   geom_point(data=ectodata, 
#              aes(x=Longitude, y=Latitude, size=log(quant95MortAge)), colour="Deep Pink", 
#              fill="Pink",pch=21, alpha=I(0.7))
# map_data_LONGEV
# 
# map_data_AGING<- 
#   base_world +
#   geom_point(data=ectodata, 
#              aes(x=Longitude, y=Latitude, size=GompertzSlope), colour="Deep Pink", 
#              fill="Pink",pch=21, alpha=I(0.7)) 
# map_data_AGING
# 
# 
# 
# 
# citation("caper")
# citation("ape")
# citation("p")