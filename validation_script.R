
#Ricardo, este é o script que usei para validação do SST, pode servir como modelo para EFA e CFA.
#Coloquei no final comandos para IRT também
###

setwd("C:/Users/user/Google Drive/Projeto_Duke/SST Validation")
data<-read.csv("sstdata.csv",header=T)
data1<-data.frame(data$q11,data$q12,data$q13,data$q14,data$q15,data$q16,data$q17,data$q18,data$q19,data$q20,data$q21,data$q22)
data2<-read.csv("sstdata1.csv",header=T)
#data<-na.omit(data1)

#Instal packages needes for the analysis
lapply(c("ggplot2", "psych", "RCurl", "irr", "nortest", "moments"), library, character.only=T)

#uploading data ------------------------------------------------------------------------
#Functions to pull the dara from the internet file 
options(RCurlOptions = list(capath = system.file("CurlSSL", "cacert.pem", package = "RCurl"), ssl.verifypeer = FALSE))
#see http://goo.gl/mQwxO on how to get this link
#link below won't work until data is entered in the right format
data <- getURL("https://docs.google.com/spreadsheet/pub?hl=en&hl=en&key=0AoTReYGK49h_dEFHWXVfODR0NlZWdXFoVTZjT09oc0E&single=true&gid=3&output=csv")
dataicc<-read.csv(textConnection(data))
names(preference.tto)

options(RCurlOptions = list(capath = system.file("CurlSSL", "cacert.pem", package = "RCurl"), ssl.verifypeer = FALSE))
#see http://goo.gl/mQwxO on how to get this link
#link below won't work until data is entered in the right format
data4 <- getURL("https://docs.google.com/spreadsheet/pub?hl=en&hl=en&key=0AoTReYGK49h_dEFHWXVfODR0NlZWdXFoVTZjT09oc0E&single=true&gid=4&output=csv")
datacor<-read.csv(textConnection(data4))
names(preference.tto)


##########################################################################################################################
#Exploratory Data Anlysis
dim(data)
summary(data)


##########################################################################################################################
#EFA
#Determine the number of factors
library(nFactors)

par(mfrow=c(2,2))
ev <- eigen(cor(data)) # get eigenvalues
ev
ap <- parallel(subject=nrow(data),var=ncol(data),rep=100,cent=.05)
nS <- nScree(ev$values)
plotnScree(nS)

#KMO
kmo = function( data1 ){
  
  library(MASS)
  X <- cor(as.matrix(data1))
  iX <- ginv(X)
  S2 <- diag(diag((iX^-1)))
  AIS <- S2%*%iX%*%S2                      # anti-image covariance matrix
  IS <- X+AIS-2*S2                         # image covariance matrix
  Dai <- sqrt(diag(diag(AIS)))
  IR <- ginv(Dai)%*%IS%*%ginv(Dai)         # image correlation matrix
  AIR <- ginv(Dai)%*%AIS%*%ginv(Dai)       # anti-image correlation matrix
  a <- apply((AIR - diag(diag(AIR)))^2, 2, sum)
  AA <- sum(a)
  b <- apply((X - diag(nrow(X)))^2, 2, sum)
  BB <- sum(b)
  MSA <- b/(b+a)                        # indiv. measures of sampling adequacy
  
  AIR <- AIR-diag(nrow(AIR))+diag(MSA)  # Examine the anti-image of the
  # correlation matrix. That is the
  # negative of the partial correlations,
  # partialling out all other variables.
  
  kmo <- BB/(AA+BB)                     # overall KMO statistic
  
  # Reporting the conclusion
  if (kmo >= 0.00 && kmo < 0.50){
    test <- 'The KMO test yields a degree of common variance
unacceptable for FA.'
  } else if (kmo >= 0.50 && kmo < 0.60){
    test <- 'The KMO test yields a degree of common variance miserable.'
  } else if (kmo >= 0.60 && kmo < 0.70){
    test <- 'The KMO test yields a degree of common variance mediocre.'
  } else if (kmo >= 0.70 && kmo < 0.80){
    test <- 'The KMO test yields a degree of common variance middling.'
  } else if (kmo >= 0.80 && kmo < 0.90){
    test <- 'The KMO test yields a degree of common variance meritorious.'
  } else {
    test <- 'The KMO test yields a degree of common variance marvelous.'
  }
  
  ans <- list(  overall = kmo,
                report = test,
                individual = MSA,
                AIS = AIS,
                AIR = AIR )
  return(ans)
  
}    # end of kmo()
kmo(data)

#FACTOR EXTRACTION
library(psych)
library(GPArotation)

fa(data1,3,fm="pa",rotate="promax")
fa(data1,1,fm="pa",rotate="promax")
fa(data1,2,fm="pa",rotate="promax")

#######################################################################################################
#CFA Models
# Sem model
library (sem)

#1 Factor Model ############################################

mod.1 <- specifyModel()# Type these values that specify the model's relations (just use de Ctrl+R over each relation).

#Latent Variables
SST->data.q11,var3,NA
SST->data.q12,var4,NA
SST->data.q13,var5,NA
SST->data.q14,var6,NA
SST->data.q15,var7,NA
SST->data.q16,var8,NA
SST->data.q17,var9,NA
SST->data.q18,var10,NA
SST->data.q19,var11,NA
SST->data.q20,var12,NA
SST->data.q21,var13,NA
SST->data.q22,var14,NA

#Erros and COv
SST<->SST,NA,1
data.q11<->data.q11,err3,NA
data.q12<->data.q12,err4,NA
data.q13<->data.q13,err5,NA
data.q14<->data.q14,err6,NA
data.q15<->data.q15,err7,NA
data.q16<->data.q16,err8,NA
data.q17<->data.q17,err9,NA
data.q18<->data.q18,err10,NA
data.q19<->data.q19,err11,NA
data.q20<->data.q20,err12,NA
data.q21<->data.q21,err13,NA
data.q22<->data.q22,err14,NA

# Insert de covariance matrix
cov <- cov(data2, y = NULL, use = "everything", method = c("pearson", "kendall", "spearman"))
cov
cor(data1)
sem.1 <- sem(mod.1, cov, N=100)
fscore(sem.wh.2)

summary(sem.1)

modIndices(sem.1)

#2 Factor Model ##################################################

mod.wh.2 <- specifyModel()# Type these values that specify the model's relations (just use de Ctrl+R over each relation).

#Latent Variables
SST3->data.q11,var3,NA
SST3->data.q12,var4,NA
SST3->data.q13,var5,NA
SST3->data.q14,var6,NA
SST1->data.q15,var7,NA
SST1->data.q16,var8,NA
SST1->data.q17,var9,NA
SST2->data.q18,var10,NA
SST1->data.q19,var11,NA
SST2->data.q20,var12,NA
SST1->data.q21,var13,NA
SST2->data.q22,var14,NA

#Erros and COv
SST1<->SST1,NA,1
SST2<->SST2,NA,1
SST3<->SST3,NA,1
data.q11<->data.q11,err3,NA
data.q12<->data.q12,err4,NA
data.q13<->data.q13,err5,NA
data.q14<->data.q14,err6,NA
data.q15<->data.q15,err7,NA
data.q16<->data.q16,err8,NA
data.q17<->data.q17,err9,NA
data.q18<->data.q18,err10,NA
data.q19<->data.q19,err11,NA
data.q20<->data.q20,err12,NA
data.q21<->data.q21,err13,NA
data.q22<->data.q22,err14,NA
data.q14<->data.q13,cov1,NA
data.q20->data.q14,cov3,NA
SST1<->SST2,lat1,NA
SST1<->SST3,lat2,NA
SST2<->SST3,lat3,NA
data.q22<->data.q21,cov4,NA
SST<->SST,lat4,NA

# Insert de covariance matrix
cov <- cov(data1, y = NULL, use = "everything", method = c("pearson", "kendall", "spearman"))
cov
cor(data1)
sem.wh.2 <- sem(mod.wh.2, cov, N=100)
fscore(sem.wh.2)

summary(sem.wh.2)

mod.indices(sem.wh.2)

# 3 Factor Model ################################################

mod.3 <- specifyModel()# Type these values that specify the model's relations (just use de Ctrl+R over each relation).

#Latent Variables
SST3->data.q11,var3,NA
SST3->data.q12,var4,NA
SST3->data.q13,var5,NA
SST3->data.q14,var6,NA
SST1->data.q15,var7,NA
SST1->data.q16,var8,NA
SST1->data.q17,var9,NA
SST2->data.q18,var10,NA
SST1->data.q19,var11,NA
SST2->data.q20,var12,NA
SST1->data.q21,var13,NA
SST2->data.q22,var14,NA

#Erros and COv
SST1<->SST1,NA,1
SST2<->SST2,NA,1
SST3<->SST3,NA,1
data.q11<->data.q11,err3,NA
data.q12<->data.q12,err4,NA
data.q13<->data.q13,err5,NA
data.q14<->data.q14,err6,NA
data.q15<->data.q15,err7,NA
data.q16<->data.q16,err8,NA
data.q17<->data.q17,err9,NA
data.q18<->data.q18,err10,NA
data.q19<->data.q19,err11,NA
data.q20<->data.q20,err12,NA
data.q21<->data.q21,err13,NA
data.q22<->data.q22,err14,NA
data.q14<->data.q13,cov1,NA
data.q20->data.q14,cov3,NA
SST1<->SST2,lat1,NA
SST1<->SST3,lat2,NA
SST2<->SST3,lat3,NA
data.q22<->data.q21,cov4,NA

# Insert de covariance matrix
cov <- cov(data1, y = NULL, use = "everything", method = c("pearson", "kendall", "spearman"))
cov
cor(data1)
sem.3 <- sem(mod.3, cov, N=100)
effects(sem.3)
standardizedCoefficients(sem.3)
standardizedResiduals(sem.3)
summary(sem.3)

mod.indices(sem.wh.2)

#2o Order Model#######################################################
# 3 Factor Model

mod.4 <- specifyModel()# Type these values that specify the model's relations (just use de Ctrl+R over each relation).

#Latent Variables
SST3->data.q11,var3,NA
SST3->data.q12,var4,NA
SST3->data.q13,var5,NA
SST3->data.q14,var6,NA
SST1->data.q15,var7,NA
SST1->data.q16,var8,NA
SST1->data.q17,var9,NA
SST2->data.q18,var10,NA
SST1->data.q19,var11,NA
SST2->data.q20,var12,NA
SST1->data.q21,var13,NA
SST2->data.q22,var14,NA
SST->SST1,NA,1
SST->SST2,lat1,NA
SST->SST3,lat2,NA

#Erros and COv
SST1<->SST1,NA,1
SST2<->SST2,NA,1
SST3<->SST3,NA,1
SST<->SST,erro1,NA
data.q11<->data.q11,err3,NA
data.q12<->data.q12,err4,NA
data.q13<->data.q13,err5,NA
data.q14<->data.q14,err6,NA
data.q15<->data.q15,err7,NA
data.q16<->data.q16,err8,NA
data.q17<->data.q17,err9,NA
data.q18<->data.q18,err10,NA
data.q19<->data.q19,err11,NA
data.q20<->data.q20,err12,NA
data.q21<->data.q21,err13,NA
data.q22<->data.q22,err14,NA
data.q14<->data.q13,cov1,NA
data.q20->data.q14,cov3,NA
#SST1<->SST2,lat1,NA
#SST1<->SST3,lat2,NA
#SST2<->SST3,lat3,NA
#data.q22<->data.q21,cov4,NA

# Insert de covariance matrix
cov <- cov(data1, y = NULL, use = "everything", method = c("pearson", "kendall", "spearman"))
cov
cor(data1)
sem.4 <- sem(mod.4, cov, N=100)
effects(sem.3)
standardizedCoefficients(sem.4)
standardizedResiduals(sem.3)
summary(sem.4)

modIndices(sem.4)

#################################################################################################

#Alpha de Cronbach by psych package
alpha(data1, keys=NULL,title=NULL,na.rm = TRUE)
alpha(data2, keys=NULL,title=NULL,na.rm = TRUE)


#Alpha de Cronbach by ltm package
install.packages("ltm")
library(ltm)
cronbach.alpha(data2$PA, standardized = TRUE, CI = TRUE, 
               probs = c(0.025, 0.975), B = 1000, na.rm = FALSE)

#################################################################################################

# ICC by psy package

install.packages("psy")
library(psy)

x<-data(expsy)
icc(expsy[,c(12,14,16)])

# ICC by irr package
attach(dataicc)
summary(dataicc)
PAICC<-data.frame(PA,PA2)
ADLICC<-data.frame(ADL,ADL2)
CRICC<-data.frame(CR,CR2)
VGICC<-data.frame(PA,PA2)

install.packages("irr")
library(irr)

cc(ratings, model = c("oneway", "twoway"), 
   type = c("consistency", "agreement"), 
   unit = c("single", "average"), r0 = 0, conf.level = 0.95)

data(dataicc)
icc(PAICC, model="twoway", type="agreement")
icc(ADLICC, model="twoway", type="agreement")
icc(CRICC, model="twoway", type="agreement")
icc(VGICC, model="twoway", type="agreement")

###########################################################
#Correlation with SF-36

summary(datacor)
detach(dataicc)
detach(datacor)
attach(datacor)       
sst<-data.frame(PA,ADL,CR,VG)
sf36<-data.frame(CF,Limit,Dor,SG,VIT,AS,Emo,SM)

cor(sst,sf36)
cor<-as.table(cor(sst,sf36))
cor.test(PA,CF,method=c("spearman"))
cor.test(PA,Limit,method=c("spearman"))
cor.test(PA,Dor,method=c("spearman"))
cor.test(PA,SG,method=c("spearman"))
cor.test(PA,VIT,method=c("spearman"))
cor.test(PA,AS,method=c("spearman"))
cor.test(PA,Emo,method=c("spearman"))
cor.test(PA,SM,method=c("spearman"))

#### IRT

install.packages("ltm")
library(ltm)

LSAT
descript(LSAT) #descriptives


fit1 <- rasch(LSAT, constraint = cbind(length(LSAT) + 1, 1)) #fitting rasch model
summary(fit1) #Show the summary for the rasch model fitting

coef(fit1, prob = TRUE, order = TRUE) #Transform de difficulties into probabilites

GoF.rasch(fit1, B = 199) #Check the fit of the model

margins(fit1) #Alternative fit check, using two or threee way X2 residuals
margins(fit1, type = "three-way", nprint = 2) #fitting for 3way ckech 

fit2 <- rasch(LSAT) #Unconstrained method/ 1 slope for the whole data
summary(fit2)

anova(fit1,fit2)#comparing discrimination rates

margins(fit2, type = "three-way", nprint = 2)

fit3 <- ltm(LSAT ~ z1) #a Two=parameter logistic model - LSAT=dichotomus variavle/z1=latent variable
summary (fit3)
anova(fit2, fit3)

fit4 <- tpm(LSAT, type = "rasch", max.guessing = 1) #Three Parameter logistic model with a guessing parameter
summary(fit4)
anova(fit2, fit4)

par(mfrow = c(2, 2))
plot(fit2, legend = TRUE, cx = "bottomright", lwd = 3,
     cex.main = 1.5, cex.lab = 1.3, cex = 1.1)
plot(fit2, type = "IIC", annot = FALSE, lwd = 3, cex.main = 1.5,
     cex.lab = 1.3)
plot(fit2, type = "IIC", items = 0, lwd = 3, cex.main = 1.5,
     cex.lab = 1.3)
plot(0:1, 0:1, type = "n", ann = FALSE, axes = FALSE)
info1 <- information(fit2, c(-4, 0))
info2 <- information(fit2, c(0, 4))
text(0.5, 0.5, labels = paste("Total Information:", round(info1$InfoTotal, 3),
                              "\n\nInformation in (-4, 0):", round(info1$InfoRange, 3),
                              paste("(", round(100 * info1$PropRange, 2), "%)", sep = ""),
                              "\n\nInformation in (0, 4):", round(info2$InfoRange, 3),
                              paste("(", round(100 * info2$PropRange, 2), "%)", sep = "")), cex = 1.5)

factor.scores(fit2) #ability to estimate 
factor.scores(fit2, resp.patterns = rbind(c(0,1,1,0,0), c(0,1,0,1,0))) #ability for nonspecific responde patterns

#IRT for ordinal variables
rcor.test(Environment, method = "kendall") #non parametric correlation, just ans alternativo to de descriptive function

fit1 <- grm(Environment, constrained = TRUE)
fit1

margins(fit1)
margins(fit1, type = "three")

fit2 <- grm(Environment)
fit2

anova(fit1, fit2)

information(fit2, c(-4, 4)) #check numerically the probabilities

information(fit2, c(-4, 4), items = c(1, 6))

par(mfrow = c(2, 2))
plot(fit2, lwd = 2, cex = 1.2, legend = TRUE, cx = "left",
     + xlab = "Latent Trait", cex.main = 1.5, cex.lab = 1.3, cex.axis = 1.1)
plot(fit2, type = "IIC", lwd = 2, cex = 1.2, legend = TRUE, cx = "topleft",
     + xlab = "Latent Trait", cex.main = 1.5, cex.lab = 1.3, cex.axis = 1.1)
plot(fit2, type = "IIC", items = 0, lwd = 2, xlab = "Latent Trait",
     + cex.main = 1.5, cex.lab = 1.3, cex.axis = 1.1)
info1 <- information(fit2, c(-4, 0))
info2 <- information(fit2, c(0, 4))
text(-1.9, 8, labels = paste("Information in (-4, 0):",
                             + paste(round(100 * info1$PropRange, 1), "%", sep = ""),
                             + "\n\nInformation in (0, 4):",
                             + paste(round(100 * info2$PropRange, 1), "%", sep = "")), cex = 1.2)

fit1 <- ltm(WIRS ~ z1 + z2)
fit2 <- ltm(WIRS ~ z1 * z2) #interaction between the latent variables improved the fit
fit2
anova(fit1, fit2)

par(mfrow = c(2, 2))
plot(fit2, category = 1, lwd = 2, cex = 1.2, legend = TRUE, cx = -4.5,
     + cy = 0.85, xlab = "Latent Trait", cex.main = 1.5, cex.lab = 1.3,
     + cex.axis = 1.1)
for (ctg in 2:3) {
  + plot(fit2, category = ctg, lwd = 2, cex = 1.2, annot = FALSE,
         + xlab = "Latent Trait", cex.main = 1.5, cex.lab = 1.3,
         + cex.axis = 1.1)
  + }
