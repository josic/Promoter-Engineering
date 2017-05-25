library(qpcR) # for PRESS
library(DAAG) # for cv.lm
library(DEoptim)
library(optimx)
library(xlsx)

fname <- "Inducible_Promoters.xlsx"

#Make sure the file Inducible_Promoters.xlsx is in the working directory
#Check the files to make sure the correct sheets are loaded
#Check the Excel file for structure of data - in this case third and fourth column
#are the ones that determine -35 and -10 condition (in that order), and the fith column
#gives the measured promoter strength.

#Set Ara = 1 or Las = 1, depending on which data set needs to be analyzed.
Ara <- 1
Las <- 0

if (Ara == 1) {sheetIndexBasal <- 1
               sheetIndexInduced <- 2}
if (Las == 1) {sheetIndexBasal <- 3
               sheetIndexInduced <- 4}



#SheetIndex refers to the sheets in Inducible_Promoters.xlsx
#1 - Ara basal
#2 - Ara induced
#3 - Las Basal
#4 - Las induced

# The other conditions correspond to the following 
# alphabetic notation in the manuscript, and corresponding sequences

# a	- baseline or 1 -	CCCGGG
# b	- 5 -	CTGACA
# c	- 4	-	TTGTGA
# d	- 2	-	TTTACA
# e	- 3	-	TAGACA
# f	- 6	- TTGACA



df1 <- read.xlsx(fname,sheetIndex = sheetIndexBasal)
df2 <- read.xlsx(fname,sheetIndex = sheetIndexInduced)
data <- rbind(df1,df2)


# Undo the logarithmic tranformation performed on the data
data$z <- 2^data$Strength

# Normalize the data
y.up <- max(data$z, na.rm=TRUE) * 1.005


# Normalize data
data$z.norm <- data$z / y.up


#Transform the data for logistic fit. 
data$nz <- log(data$z.norm / (1 - data$z.norm))



# Define the different factors for the analysis
# Here and below in both the Las and Ara case
#     Condition 1 stands for uninduced
#     Condition 2 stands for inducded 

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


#DEFINE THE COST FUNCTION
#A nonlinear fit ------- fit amplitude AND epsilon (see model)
#Note that we are fitting to the 
model.func_Ampe <- function(in.val, bool.predict=FALSE){
  
  #These are the values of Delta G_ij in the notes. Here we set -b0 = Delta G_10,1 + Delta G_35,1, bind = gamma_+ (See model description).
  #The other terms are measured relative to b0, so that -Delta G12 = b0 + bm352 
  #and -DeltaG22 = b0 + bm352 + bm102, for example 
  #Note the change in sign as well.  
  b0 <- in.val[1]
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
  bind <- in.val[14]
  
  #This is an amplitude
  Amp <- in.val[15]
  
  #This is an additional offest
  epsilon <- in.val[16]
  
  bm35 <- c(bm352, bm353, bm354, bm355, bm356)
  vm35 <- rep(0,5)
  bm10 <- c(bm102, bm103, bm104, bm105, bm106, bm107, bm108)
  vm10 <- rep(0,7)
  error <- 0
  
  
  #Loop over first condition onlyfi
  cond <- 1 
  i.row <- 1
  
  #first row for condition 2
  shift <- 48
  while(newdata$Condition[i.row] == cond){
    condition <- newdata[i.row,2]
    vm35i <- as.integer(newdata[i.row, 'Minus35'])
    vm10i <- as.integer(newdata[i.row, 'Minus10'])
    if(vm35i != 1) {vm35[vm35i-1] <- 1}
    if(vm10i != 1) {vm10[vm10i-1] <- 1}
    exp_factor_0 <- b0 + vm10 %*% bm10 + vm35 %*% bm35
    exp_factor_ind <- b0+ bind + vm10 %*% bm10 + vm35 %*% bm35
    
    #Predicted values from model
    pz_0 <- Amp*(exp(exp_factor_0)/(1+exp(exp_factor_0))) + epsilon
    pz_ind <-  Amp*(exp(exp_factor_ind)/(1+exp(exp_factor_ind))) + epsilon
    
    
    #This can be used as a sanity check that the minimization recovers linear model
    #error <- error + (newdata$nz[i.row]-exp_factor_0)^2 + (newdata$nz[i.row+48]-exp_factor_ind)^2
    
    #Minimize difference in observed and predicted values
    #     if(!is.na(newdata[i.row,]$z)) {
    #           error <- error +  (pz_0 - newdata$z[i.row])^2
    #     }
    #     if(!is.na(newdata[i.row+shift,]$z)) {
    #           error <- error + (pz_ind - newdata$z[i.row+shift])^2
    #     }
    
    
    #Minimize difference in logs squared
    if(!is.na(newdata[i.row,]$z)) {
      error <- error +  (log(pz_0) - log(newdata$z[i.row]))^2
    }
    if(!is.na(newdata[i.row+shift,]$z)) {
      error <- error + (log(pz_ind) - log(newdata$z[i.row+shift]))^2
    }
    
    # reset vm35 and vm10
    vm35 <- rep(0, 5)
    vm10 <- rep(0, 7)
    i.row <- i.row +1
  }
  
  return (error)  
  
}

# This is the call to DEoptim, the function that minimizes the cost defined in the model above
# This takes a bit to run, but gives the best result among the optimization results we tried.
# itermax can be increased for a slight improvement in fit, but 1000 is fairly good.
model.fit3 <- DEoptim(model.func_Ampe, c(-40,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0), c(20,20,20,20,20,20,20,20,20,20,20,20,20,20,1000000,10000),  control=list(NP=160, itermax=5000))

# These vectors will be used below
newdata$Score <- numeric(length(newdata$z))
newdata$energy <- numeric(length(newdata$z))


# Save the parameters obtained from the fit in a vector called par 
par <- model.fit3$optim$bestmem


# Create vector of scores
# The scores are computed using the same indexing explained when introducing the bm* parameters 
# in the cost function above.

for (i in 1:6){
  
  for (j in 1:8){
    if(i==1 & j == 1) {eng0 <- par[1]
                      eng1 <- eng0 + par[14]} 
      else if (i==1)  {eng0 <- par[1] +  par[j+5]
                       eng1 <- eng0 + par[14]}
      else if (j==1)  {eng0 <- par[1] + par[i]
                       eng1 <- eng0+ par[14]} 
      else { eng0 <- par[1] + par[i] + par[j+5] 
             eng1 <- eng0 + par[14]} 
    newdata$energy[newdata$Minus10 == j & newdata$Minus35 == i & newdata$Condition == 1] <- eng0
    newdata$energy[newdata$Minus10 == j & newdata$Minus35 == i & newdata$Condition == 2] <- eng1
    newdata$Score[newdata$Minus10 == j & newdata$Minus35 == i & newdata$Condition == 1] <- eng0
    newdata$Score[newdata$Minus10 == j & newdata$Minus35 == i & newdata$Condition == 2] <- eng0
    }
}



index1 <- (newdata$Condition==1)



#Plot log(z) ---- epsilon, and amplitude
amp <- par[15]
epsilon <- par[16]


plot(newdata[index1,]$Score, log(newdata[index1,]$z), pch=20, ylim=c(min(log(newdata$z),na.rm=TRUE), max(log(newdata$z),na.rm=TRUE)), xlim=c(min(newdata$Score), max(newdata$Score)), xlab='Combined Score', ylab='log z')
points(newdata[!index1,]$Score, log(newdata[!index1,]$z), col="red", pch=20)
points(newdata[index1,]$Score, log(amp*(exp(newdata[index1,]$energy)/(1 + exp(newdata[index1,]$energy)))+epsilon), pch=4, col="green")
points(newdata[!index1,]$Score, log(amp*(exp(newdata[!index1,]$energy)/(1 + exp(newdata[!index1,]$energy)))+epsilon), pch=4, col="blue")




# Compute the stdv of the three measurements for plotting
# Note that the fit was to the mean

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

# Output to file

logz <- log(newdata[index1,]$z)
prediction <-  log(amp*(exp(newdata[index1,]$energy)/(1 + exp(newdata[index1,]$energy)))+epsilon)
score <- newdata[index1,]$Score
sd <- newdata[index1,]$sd
out_uninduced <- cbind(score,logz,prediction,sd)

if (Ara == 1) {fileout = 'coordsAra_uninduced'}
if (Las == 1) {fileout = 'coordsLas_uninduced'}
write.csv(out_uninduced,file=fileout,row.names=FALSE)


# Prepare uninduced fit data
logz <- log(newdata[!index1,]$z)
prediction <- log(amp*(exp(newdata[!index1,]$energy)/(1 + exp(newdata[!index1,]$energy)))+epsilon)
score <- newdata[!index1,]$Score   
sd <- newdata[!index1,]$sd
out_induced <- cbind(score,logz,prediction,sd)

# Write induced fit data
if (Ara == 1) {fileout = 'coordsAra_induced'}
if (Las == 1) {fileout = 'coordsLas_induced'}
write.csv(out_induced,file=fileout,row.names=FALSE)


#Write parameters
if (Ara == 1) {fileout = 'parameters_Ara'}
if (Las == 1) {fileout = 'parameters_Las'}
write.csv(par,file = fileout)
