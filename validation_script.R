
<<<<<<< HEAD
=======
#commit test

#Ricardo, este é o script que usei para validação do SST, pode servir como modelo para EFA e CFA.
#Coloquei no final comandos para IRT também
###

setwd("C:/Users/user/Google Drive/Projeto_Duke/SST Validation")
data<-read.csv("sstdata.csv",header=T)
data1<-data.frame(data$q11,data$q12,data$q13,data$q14,data$q15,data$q16,data$q17,data$q18,data$q19,data$q20,data$q21,data$q22)
data2<-read.csv("sstdata1.csv",header=T)
#data<-na.omit(data1)

>>>>>>> 3ea635490a1fcd36bcc2d9110a512b3b9d2134fc
#Instal packages needes for the analysis
lapply(c("ggplot2", "psych", "RCurl", "irr", "nortest", "moments"), library, character.only=T)

#uploading data ------------------------------------------------------------------------
#Functions to pull the dara from the internet file 
options(RCurlOptions = list(capath = system.file("CurlSSL", "cacert.pem", package = "RCurl"), ssl.verifypeer = FALSE))
#see http://goo.gl/mQwxO on how to get this link
data1 <- getURL("https://docs.google.com/spreadsheet/pub?key=0AoTReYGK49h_dEFHWXVfODR0NlZWdXFoVTZjT09oc0E&single=true&gid=1&output=csv
")
data<-read.csv(textConnection(data1))

#data for ICC analysis
data2 <- getURL("https://docs.google.com/spreadsheet/pub?hl=en&hl=en&key=0AoTReYGK49h_dEFHWXVfODR0NlZWdXFoVTZjT09oc0E&single=true&gid=3&output=csv")
dataicc<-read.csv(textConnection(data2))

#data for correlation analysis
data4 <- getURL("https://docs.google.com/spreadsheet/pub?hl=en&hl=en&key=0AoTReYGK49h_dEFHWXVfODR0NlZWdXFoVTZjT09oc0E&single=true&gid=4&output=csv")
datacor<-read.csv(textConnection(data4))

##########################################################################################################################
#Exploratory Data Anlysis
dim(data)
summary(data)


##########################################################################################################################
#EFA
#Determine the number of factors

#Package to be installed to EFa ANalysis
library(nFactors)

#Group of functinos to determine the number os items to be extracted
par(mfrow=c(2,2)) #Command to configure the plot area for the scree plot graph
ev <- eigen(cor(data)) # get eigenvalues - insert the data you want to calculate the scree plot for
ev # Show eigend values
ap <- parallel(subject=nrow(data),var=ncol(data),rep=100,cent=.05) #Calculate the acceleration factor
nS <- nScree(ev$values) #Set up the Scree Plot 
plotnScree(nS) # Plot the Graph

#Functino to calculate the KMO values
kmo = function(data)
  
  {
  
  library(MASS)
  X <- cor(as.matrix(data1))
  iX <- ginv(X)
  S2 <- diag(diag((iX^-1)))
  AIS <- S2%*%iX%*%S2                      
  IS <- X+AIS-2*S2                         
  Dai <- sqrt(diag(diag(AIS)))
  IR <- ginv(Dai)%*%IS%*%ginv(Dai)         
  AIR <- ginv(Dai)%*%AIS%*%ginv(Dai)       
  a <- apply((AIR - diag(diag(AIR)))^2, 2, sum)
  AA <- sum(a)
  b <- apply((X - diag(nrow(X)))^2, 2, sum)
  BB <- sum(b)
  MSA <- b/(b+a)                        
    AIR <- AIR-diag(nrow(AIR))+diag(MSA)  
    kmo <- BB/(AA+BB)                     
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
  
}    # end of the function kmo()

kmo(data) #Run the Kmo function for the data you want to calculate

#Functions to exctract factor loadings for differente factor structures

#Packages needed to run the Data exctraction function
library(psych)
library(GPArotation)

fa(data,3,fm="pa",rotate="promax") #Functino to exctract the factor loadings. Arguments are DATA, Number of factors, rotation method. Look here http://goo.gl/kY3ln for different methods of estimations or rotations

#######################################################################################################
#Using SEM package to rund CFA Models
# CFA model

#Package you will need to run the CFA Model
library (sem)

# Function to specify the model you will run
# You need to type these values that specify the model's relations please see DDD for further instructions
#Summarizing here you will create the paths for the SEM (CFA) model. Use -> for unidirectional relatino and <-> for correlation, covariance or error.
mod.1 <- specifyModel()

#Defining Latent Variables - Set the variable you wnat to create (SST, in the example is the latent variable)
# Put the arrow for direction ->, the arguments are: Latent Variable, observed variable, code, NA (Indicating that you want the program to define this loading number)
SST->q11,var3,NA
SST->q12,var4,NA
SST->q13,var5,NA
SST->q14,var6,NA
SST->q15,var7,NA
SST->q16,var8,NA
SST->q17,var9,NA
SST->q18,var10,NA
SST->q19,var11,NA
SST->q20,var12,NA
SST->q21,var13,NA
SST->q22,var14,NA

#You always have to specify the erros and COvariances, so you will use <-> between the variables
# In the example, we are setting the factor weight for the first error to 1 so we change the order, NA goes int he code argument and 1 to the loading argument (associated with the latent variable is one of doing it)
SST<->SST,NA,1
q11<->q11,err3,NA
q12<->q12,err4,NA
q13<->q13,err5,NA
q14<->q14,err6,NA
q15<->q15,err7,NA
q16<->q16,err8,NA
q17<->q17,err9,NA
q18<->q18,err10,NA
q19<->q19,err11,NA
q20<->q20,err12,NA
q21<->q21,err13,NA
q22<->q22,err14,NA

# Insert de covariance matrix - CFA (or SEM) is always calculated in relation to a covariance or correlation matrix, here we will create the covariance matrix
cov <- cov(data, y = NULL, use = "everything", method = c("pearson", "kendall", "spearman"))
cov #Diplat covariance matrix
cor(data) #Display correlation matrix (Only if you need to check the correlation matrix)

#Run SEM functino to get the SEM model
sem.1 <- sem(mod.1, cov, N=100) #mod.1 is the model you ahve crated before, cov is the covriance matrix and N argument is the number of observation
summary(sem.1) # Show summary resutls for SEM
modIndices(sem.1) #Show modification indexes to enhance goodness of fit (suggested to change the model for modification indexes higher than 40, but it is your call)
standardizedCoefficients(sem.1) #Get the standardize coefficientes for each path
standardizedResiduals(sem.1) #Get the residuals for the regressions models

#################################################################################################

#Alpha de Cronbach by psych package - Gives lots if indicatores, but doesn't give the CI (must have a way that I don't know yet)
alpha(data, keys=NULL,title=NULL,na.rm = TRUE)

#Alpha de Cronbach by ltm package - GIves CI 
install.packages("ltm")
library(ltm)
cronbach.alpha(data, standardized = TRUE, CI = TRUE, 
               probs = c(0.025, 0.975), B = 1000, na.rm = FALSE)

#################################################################################################

# ICC by irr package
attach(dataicc)
summary(dataicc)
#Creating data frames for ICC analysis
PAICC<-data.frame(PA,PA2)
ADLICC<-data.frame(ADL,ADL2)
CRICC<-data.frame(CR,CR2)
VGICC<-data.frame(PA,PA2)

install.packages("irr")
library(irr)

#Command to give ICC
icc(PAICC, model="twoway", type="agreement")
icc(ADLICC, model="twoway", type="agreement")
icc(CRICC, model="twoway", type="agreement")
icc(VGICC, model="twoway", type="agreement")

###########################################################
#Correlation with SF-36

#Clean workspace
summary(datacor)
detach(dataicc)
detach(datacor)
attach(datacor)       
#Creating correlations vectors
sst<-data.frame(PA,ADL,CR,VG)
sf36<-data.frame(CF,Limit,Dor,SG,VIT,AS,Emo,SM)

cor(sst,sf36)
cor<-as.table(cor(sst,sf36)) #Give correlation matrix
cor.test(PA,CF,method=c("spearman")) #Give correlation significance test
cor.test(PA,Limit,method=c("spearman"))
cor.test(PA,Dor,method=c("spearman"))
cor.test(PA,SG,method=c("spearman"))
cor.test(PA,VIT,method=c("spearman"))
cor.test(PA,AS,method=c("spearman"))
cor.test(PA,Emo,method=c("spearman"))
cor.test(PA,SM,method=c("spearman"))

#### IRT Analysis - Still need to work on it. Have tried only once, but will dominate this when we need to use it

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
