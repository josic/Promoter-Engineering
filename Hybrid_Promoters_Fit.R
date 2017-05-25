library(qpcR) # for PRESS
library(DAAG) # for cv.lm
library(GenSA)
library(DEoptim)
library(optimx)
library(xlsx)

#Make sure the file Hybridpromoterdata.xlsx is in the working directory
#Check the files to make sure the correct sheets are loaded
#Check the Excel file for structure of data. Note that the 
#data is organized somewhat differently with the column Condition
#determining the state 

# Condition.  Note that there are four conditions. We denoted the
# four conditions by ++, +−, −+, and −− depending on whether 
# the first sign denoted inducer absent/present, and the second sign 
# denotes repressed state (ie repressed = -, unrepressed = +).  We assumed
# expression was not fully repressed.  Therefore we used a mixed model 
# which assumes a fraction of promoters sites are not repressed.  
# The four conditions in the data

# 1 - this is the -/- state in the notes, repressed and no activator
# 2 - this is the +/- state in the notes, repressed with activator
# 3 - this is the -/+ state in the notes, uninduced
# 4 - this is the +/+ state in the notes, activator, but no repression

# The other conditions correspond to the following 
# alphabetic notation in the manuscript, and corresponding sequences

# A - baseline or	1	- 	CCAGTC
# B	- 2	-	TATGTT
# C -	6 -	TACTGT
# D	- 5	-	TAAATT
# E	- 4	-	GATACT
# F	- 7	-	GATAAT
# G	- 3	-	TATAGT
# H	- 8	-	TATAAT
a

fname <- "Hybridpromoterdata.xlsx"

#Set the data to be analyzed.  Set only one of these to 1.
Ara <- 0
Las <- 1

if(Ara == 1) {sheet <- 1}
if(Las == 1) {sheet <- 2}

#Check the files to make sure the correct sheets are loaded
data <- read.xlsx(fname,sheetIndex = sheet)

# Undo the logarithmic tranformation performed on the data
data$z <- 2^data$Strength

# Normalize the data
y.up <- max(data$z, na.rm=TRUE) * 1.005

# Normalize data
data$z.norm <- data$z / y.up

#Transform the data for logistic fit. 
data$nz <- log(data$z.norm / (1 - data$z.norm))


# Define the different factors
data$Condition <- factor(data$Condition)
data$Activator <- factor(data$Activator)
data$Minus35 <- factor(data$Minus35)
data$Minus10 <- factor(data$Minus10)

# This creates a reduced set where the three replications are collapsed to their mean
# Thus newdata has a third of the data points
attach(data)
newdata <- data.frame()
for (i in levels(Condition)){
  
  for (j in levels(Minus35)){
    
    for (k in levels(Minus10)){
      temp_nz <- data[Condition==i & Minus35==j & Minus10==k, 9]
      temp_z <- data[Condition==i & Minus35==j & Minus10==k, 7]
      temp_znorm <- data[Condition==i & Minus35==j & Minus10==k, 8]
      newdata <- rbind(newdata, data.frame(Condition=i, Minus35=j, Minus10=k, z = mean(temp_z,na.rm=TRUE ), z.norm=mean(temp_znorm,na.rm=TRUE), nz=mean(temp_nz, na.rm=TRUE)))
    }
  }
}
detach(data)

#Fit linear model to the average of the triplicates for comparison
model.linear.averages <- lm(nz~ Minus35+Minus10+Condition, data = newdata)

#Before fitting the modles specify conditions which we will be fitting
#As before we will first fit two conditions



#A nonlinear fit ----- fit epsilon and amplitude
# The idea is exactly the same as in the single transcription factor case.
# The coefficients b[1] - b[4] represent the terms gamma + Delta G_10,1 + Delta G_35,1
# in the four conditions --,+-,-+, and ++ (See Methods). 

# We allow for the possibility that repression may not be full.  Thus some sites may be 
# unrepressed. The fraction (or probability) of unrepressed sites is given by the parameter p.

# The baseline and amplitude, epsilon, and amp, have the same meaning as in the 
# single transcription factor case.

model.func_ampeps_mixed <- function(in.val, bool.predict=FALSE){
  
  b <- c(1:4)
  
  b[1] <- in.val[1]
  bm352 <- in.val[2]
  bm353 <- in.val[3]
  bm354 <- in.val[4]
  bm355 <- in.val[5]
  bm356 <- in.val[6]
  bm102 <- in.val[7]
  bm103 <- in.val[8]
  bm104 <- in.val[9]
  bm105 <- in.val[10]
  bm106 <- in.val[11]
  bm107 <- in.val[12]
  bm108 <- in.val[13]
  b[2] <- in.val[14]
  b[3] <- in.val[15]
  b[4] <- in.val[16]
  
  epsilon <- in.val[17]
  amp <- in.val[18]
  p <- in.val[19]
  
  bm35 <- c(bm352, bm353, bm354, bm355, bm356)
  vm35 <- rep(0,5)
  bm10 <- c(bm102, bm103, bm104, bm105, bm106, bm107, bm108)
  vm10 <- rep(0,7)
  error <- 0
  
  
  
  for(i.row in 1:nrow(newdata)){
    condition <- as.numeric(newdata[i.row,1])
    
    vm35i <- as.integer(newdata[i.row, 'Minus35'])
    vm10i <- as.integer(newdata[i.row, 'Minus10'])
    if(vm35i != 1) {vm35[vm35i-1] <- 1}
    if(vm10i != 1) {vm10[vm10i-1] <- 1}
    exp_factor<- b[1] + vm10 %*% bm10 + vm35 %*% bm35
    
    #If we are in a different condition, add different coefficient
    #Here the parameter b depends on the conditions as defined above.
    if(condition != 1){
      exp_factor<- exp_factor + b[condition]
    }
    
    #In the repressed case, we allow for the possibility that some copies of the
    #promoter are not repressed. We therefore save the energies in the unrepressed
    #conditions (conditions 3 and 4) to define the mixed models
    if (condition == 1 || condition ==2) {
      exp_factor_norepress <- b[1] + vm10 %*% bm10 + vm35 %*% bm35 + b[condition+2]
      
    }
    
    #Predicted values from model - 
    # in conditions 1 and 2, ie in the repressed
    # state these are a mixture of two states, with probability p of being repressed.
    if(condition == 3 || condition == 4){
      pz <- amp * (exp(exp_factor)/(1+exp(exp_factor))) + epsilon}
    else {
      pz <- amp * (p * (exp(exp_factor)/(1+exp(exp_factor))) + 
                    (1-p) * (exp(exp_factor_norepress)/(1+exp(exp_factor_norepress)))) + epsilon
    }
    
    
    #Minimize difference in logs
    if(!is.na(newdata[i.row,]$z)) {
      error <- error +  (log(pz) - log(newdata$z[i.row]))^2
    }
    
    
    # reset vm35 and vm10
    vm35 <- rep(0, 5)
    vm10 <- rep(0, 7)
    i.row <- i.row +1
  }
  
  return (error)  
  
}





model.fit3 <- DEoptim(model.func_ampeps_mixed, 
          c(-40,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.9), c(-10,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,100,200000,1),  
          control=list(NP=190, itermax=5000))


#Create a new parameter vector
#Use this one for DEoptim, ie model 3
par<- model.fit3$optim$bestmem


#Some housekeeping due to the definition of the offsets b[1] - b[4]
b<-c(1:4)
b[1] <- par[1]
b[2] <- par[14] + b[1]
b[3] <- par[15] + b[1]
b[4] <- par[16] + b[1]



index1 <- (newdata$Condition==1)
index2 <- (newdata$Condition==2)
index3 <- (newdata$Condition==3)
index4 <- (newdata$Condition==4)



#Create vector of scores for MIXED MODEL
for (k in 4:1){
  
  for (j in 1:8){
    for (i in 1:6) {
      if(i==1 & j == 1) {eng <- b[k]  } 
      else if (i==1)  {eng <- b[k] +  par[j+5] }
      else if (j==1)  {eng <- b[k] + par[i]} 
      else { eng <- b[k] + par[i] + par[j+5] }
      
      newdata$energy[newdata$Minus10 == j & newdata$Minus35 == i & newdata$Condition == k] <- eng
      newdata$Score[newdata$Minus10 == j & newdata$Minus35 == i & newdata$Condition == k] <- 
        newdata$energy[newdata$Minus10 == j & newdata$Minus35 == i & newdata$Condition == 4]
    }
  }
}

#Plot log(z) ---- MIXED MODEL
epsilon <- par[17]
amp <- par[18]
p <- par[19]

#Compute predictions
pred3 <- amp *exp(newdata[index3,]$energy)/(1 + exp(newdata[index3,]$energy))+epsilon
pred4 <- amp *exp(newdata[index4,]$energy)/(1 + exp(newdata[index4,]$energy))+epsilon
pred1 <- amp * (p * (exp(newdata[index1,]$energy)/(1+exp(newdata[index1,]$energy))) + 
                  (1-p) * (exp(newdata[index3,]$energy)/(1+exp(newdata[index3,]$energy)))) + epsilon
pred2 <- amp * (p * (exp(newdata[index2,]$energy)/(1+exp(newdata[index2,]$energy))) + 
                  (1-p) * (exp(newdata[index4,]$energy)/(1+exp(newdata[index4,]$energy)))) + epsilon


#Plot the results
plot(newdata[index1,]$Score, log(newdata[index1,]$z), pch=20, ylim=c(min(log(newdata$z),na.rm=TRUE), max(log(newdata$z),na.rm=TRUE)), xlim=c(min(newdata$Score), max(newdata$Score)), xlab='Combined Score', ylab='log z')
points(newdata[index2,]$Score, log(newdata[index2,]$z), col="red", pch=20)
points(newdata[index3,]$Score, log(newdata[index3,]$z), col="blue", pch=20)
points(newdata[index4,]$Score, log(newdata[index4,]$z), col="green", pch=20)
points(newdata[index1,]$Score, log(pred1), pch=4)
points(newdata[index2,]$Score, log(pred2), pch=4,col="red")
points(newdata[index3,]$Score, log(pred3), pch=4,col="blue")
points(newdata[index4,]$Score, log(pred4), pch=4,col="green")

newdata$sd <- NaN
attach(data)
for (i in levels(Condition)){
  
  for (j in levels(Minus35)){
    
    for (k in levels(Minus10)){
      temp_logz <- log(data[Condition==i & Minus35==j & Minus10==k, 7])
      if (is.na(sd(temp_logz, na.rm=TRUE)) == FALSE){
        newdata[newdata$Condition==i & newdata$Minus35==j & newdata$Minus10==k,]$sd <- sd(temp_logz, na.rm=TRUE)
      }
    }
  }
}
detach(data)


#Output to file
logz <- log(newdata$z)
score <- newdata$Score  
sdev <- newdata$sd
condition <- newdata$Condition
out <- cbind(condition,score,logz,sdev)

if(Ara == 1) {fileout = "coordsHybrid_Ara_mixed"} else {fileout = "coordsHybrid_Las_mixed"}
write.csv(out,fileout,row.names=FALSE)

if(Ara == 1) {fileout = "parametersHybrid_Ara_mixed"} else {fileout = "parametersHybrid_Las_mixed"}
write.csv(par,fileout)
