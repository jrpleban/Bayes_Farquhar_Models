## Bayesian Estimate of PS paramters from Brassica ACi curves ###
setwd("~/Documents/ACiP/Scripts")   ### Change 4 server
library("rjags", lib.loc="/Library/Frameworks/R.framework/Versions/3.1/Resources/library") ### Change 4 server
#source("farQ_TPU.R")  ### Use output pars to predict net assimilation 
source("CiCc_Jf_model.R") ### Bayesian Modelw/ Temp Dependency & farQ limitation

#### Set up for rjags #########
parameters = c("Vcmax", "Rd","gm",
               "gammaS" ,"Kc", "Ko","tau") ### pars to be monitored
adaptSteps = 2000             # Number of steps to "tune" the samplers.
burnInSteps = 3000            # Number of steps to "burn-in" the samplers.
nChains = 4                   # Number of chains to run.
DICsteps= 20000                # Number of steps of sample DIC
numSavedSteps= 10000        # Total number of steps in chains to save.
thinSteps=50                   # Number of steps to "thin" (1=keep every step).
nPerChain = ceiling( ( numSavedSteps * thinSteps ) / nChains ) # Steps per chain.

######################
#####DATA SET UP######
######################
ACdat<-read.delim("~/Documents/ACiP/Data/ACdata_reduced1_clean.txt")  ### Change 4 server
useID<-read.delim("~/Documents/ACiP/Data/IDs_to_use-2-19-15.txt")     ### Change 4 server
gID<-useID$ID
IDgeno<-as.character(useID$geno)
#row.has.na <- apply(ACdat, 1, function(x){any(is.na(x))})
#AC.filtered <- ACdat[!row.has.na,]
#names(ACi)
#Bad<-c(182, 2312, 2319, 1937, 1872, 2653) -- based on visual inspectin
ACi<-ACdat[ACdat$AC.ID %in% gID, ]
#ACi<-AC.filtered2[order(AC.filtered2$Photo),]
######Pull out geno###

####  DATA SET  UP  ############
#### PULL OUT Data for each Curve   ###
ID<-unique(ACi$AC.ID)

#### Photo Data

A = matrix(ncol=18, nrow=12)
for(i in 1:18){
  A[,i]<-ACi[ACi$AC.ID==ID[i],]$Photo
}


## CP = CO2 concentraion intercelluar space (in partial pressure)
ACi$CP<-ACi$Ci*ACi$Press/1000   ## calcs CO2 concentraion intercelluar space (in partial pressure)
CP = matrix(ncol=18, nrow=12)   ### pulls for each ID
for(i in 1:18){
  CP[,i]<-ACi[ACi$AC.ID==ID[i],]$CP
}


### J using Ch Flour estimate Elelcton tranport Rate ###
Jf = matrix(ncol=18, nrow=12)   ### pulls for each ID
for(i in 1:18){
  Jf[,i]<-ACi[ACi$AC.ID==ID[i],]$ETR
}

###  O2 in pp from atmospheric pressure (kPa)
### PP(Pa) =  MolFact * pressure of the gas mixture
ACi$O<-(.21)*(ACi$Press*1000)
O = matrix(ncol=18, nrow=12)   ### pulls for each ID
for(i in 1:18){
  O[,i]<-ACi[ACi$AC.ID==ID[i],]$O
}

N<-12

datalist1<-list(N=N, An=A[,1], CiP=CP[,1], Jf=Jf[,1], O=O[,1])
datalist2<-list(N=N, An=A[,2], CiP=CP[,2], Jf=Jf[,2], O=O[,2])
datalist3<-list(N=N, An=A[,3], CiP=CP[,3], Jf=Jf[,3], O=O[,3])
datalist4<-list(N=N, An=A[,4], CiP=CP[,4], Jf=Jf[,4], O=O[,4])
datalist5<-list(N=N, An=A[,5], CiP=CP[,5], Jf=Jf[,5], O=O[,5])
datalist6<-list(N=N, An=A[,6], CiP=CP[,6], Jf=Jf[,6], O=O[,6])
datalist7<-list(N=N, An=A[,7], CiP=CP[,7], Jf=Jf[,7], O=O[,7])
datalist8<-list(N=N, An=A[,8], CiP=CP[,8], Jf=Jf[,8], O=O[,8])
datalist9<-list(N=N, An=A[,9], CiP=CP[,9], Jf=Jf[,9], O=O[,9])
datalist10<-list(N=N, An=A[,10], CiP=CP[,10], Jf=Jf[,10], O=O[,10])
datalist11<-list(N=N, An=A[,11], CiP=CP[,11], Jf=Jf[,11], O=O[,11])
datalist12<-list(N=N, An=A[,12], CiP=CP[,12], Jf=Jf[,12], O=O[,12])
datalist13<-list(N=N, An=A[,13], CiP=CP[,13], Jf=Jf[,13], O=O[,13])
datalist14<-list(N=N, An=A[,14], CiP=CP[,14], Jf=Jf[,14], O=O[,14])
datalist15<-list(N=N, An=A[,15], CiP=CP[,15], Jf=Jf[,15], O=O[,15])
datalist16<-list(N=N, An=A[,16], CiP=CP[,16], Jf=Jf[,16], O=O[,16])
datalist17<-list(N=N, An=A[,17], CiP=CP[,17], Jf=Jf[,17], O=O[,17])
datalist18<-list(N=N, An=A[,18], CiP=CP[,18], Jf=Jf[,18], O=O[,18])
#################################
##################################
##### Impliment model in JAGS ####
#################################
### running each curve
print("initialize models")
model1 <- jags.model(textConnection(CiCc_Jf), 
                     data = datalist1, n.chains=nChains , n.adapt=adaptSteps)
model2 <- jags.model(textConnection(CiCc_Jf), 
                     data = datalist2, n.chains=nChains , n.adapt=adaptSteps)
model3 <- jags.model(textConnection(CiCc_Jf), 
                     data = datalist3, n.chains=nChains , n.adapt=adaptSteps)
model4 <- jags.model(textConnection(CiCc_Jf), 
                     data = datalist4, n.chains=nChains , n.adapt=adaptSteps)
model5 <- jags.model(textConnection(CiCc_Jf), 
                     data = datalist5, n.chains=nChains , n.adapt=adaptSteps)
model6 <- jags.model(textConnection(CiCc_Jf), 
                     data = datalist6, n.chains=nChains , n.adapt=adaptSteps)
model7 <- jags.model(textConnection(CiCc_Jf), 
                     data = datalist7, n.chains=nChains , n.adapt=adaptSteps)
model8 <- jags.model(textConnection(CiCc_Jf), 
                     data = datalist8, n.chains=nChains , n.adapt=adaptSteps)
model9 <- jags.model(textConnection(CiCc_Jf), 
                     data = datalist9, n.chains=nChains , n.adapt=adaptSteps)
model10 <- jags.model(textConnection(CiCc_Jf), 
                      data = datalist10, n.chains=nChains , n.adapt=adaptSteps)
model11 <- jags.model(textConnection(CiCc_Jf), 
                      data = datalist11, n.chains=nChains , n.adapt=adaptSteps)
model12 <- jags.model(textConnection(CiCc_Jf), 
                      data = datalist12, n.chains=nChains , n.adapt=adaptSteps)
model13 <- jags.model(textConnection(CiCc_Jf), 
                      data = datalist13, n.chains=nChains , n.adapt=adaptSteps)
model14 <- jags.model(textConnection(CiCc_Jf), 
                      data = datalist14, n.chains=nChains , n.adapt=adaptSteps)
model15 <- jags.model(textConnection(CiCc_Jf), 
                      data = datalist15, n.chains=nChains , n.adapt=adaptSteps)
model16 <- jags.model(textConnection(CiCc_Jf), 
                      data = datalist16, n.chains=nChains , n.adapt=adaptSteps)
model17 <- jags.model(textConnection(CiCc_Jf), 
                      data = datalist17, n.chains=nChains , n.adapt=adaptSteps)
model18 <- jags.model(textConnection(CiCc_Jf), 
                      data = datalist18, n.chains=nChains , n.adapt=adaptSteps)



#################################
print("updating")
update(model1, burnInSteps) # Burnin for burnInSteps samples
update(model2, burnInSteps)
update(model3, burnInSteps)
update(model4, burnInSteps)
update(model5, burnInSteps)
update(model6, burnInSteps)
update(model7, burnInSteps) # Burnin for burnInSteps samples
update(model8, burnInSteps)
update(model9, burnInSteps)
update(model10, burnInSteps)
update(model11, burnInSteps)
update(model12, burnInSteps)
update(model13, burnInSteps) # Burnin for burnInSteps samples
update(model14, burnInSteps)
update(model15, burnInSteps)
update(model16, burnInSteps)
update(model17, burnInSteps)
update(model18, burnInSteps)

##########################################

###### SAMPLE all 18 individual curves ###

##########################################


##########################################
print("sampling chains")
##### mcmc_samples  model 1 IDs  #####
#####   ID order   1845 1894 1902 1937 2157 2208 2304 2319 2320 2349 2603 2629 2655 2696 2712 2723 2737 2760
mcmc_samples1<- coda.samples(model1,
variable.names=parameters,
n.iter=nPerChain , thin=thinSteps )
####### Plot results #####
#plot(mcmc_samples1)
mcmcChain1845_1 = as.matrix( mcmc_samples1)
chainLength = NROW(mcmcChain1845_1)
# Convert precision (tau) to SD###
sigma =1  / sqrt( mcmcChain1845_1[, "tau" ] )
mcmcChain1845_1 = as.data.frame(cbind( mcmcChain1845_1, sigma ))
g1<-gelman.diag(mcmc_samples1)


####
mcmc_samples2<- coda.samples(model2,
variable.names=parameters,
n.iter=nPerChain , thin=thinSteps )
mcmcChain1894_1 = as.matrix( mcmc_samples2)
# Convert precision (tau) to SD###
sigma =1  / sqrt( mcmcChain1845_1[, "tau" ] )
mcmcChain1894_1 = as.data.frame(cbind( mcmcChain1894_1, sigma ))
g2<-gelman.diag(mcmc_samples2)

######
mcmc_samples3 <- coda.samples(model3,
variable.names=parameters,
n.iter=nPerChain , thin=thinSteps )
mcmcChain1902_1 = as.matrix( mcmc_samples3)
# Convert precision (tau) to SD###
sigma = 1  / sqrt( mcmcChain1902_1[, "tau" ] )
mcmcChain1902_1 = as.data.frame(cbind( mcmcChain1902_1, sigma ))
g3<-gelman.diag(mcmc_samples3)

######
mcmc_samples4 <- coda.samples(model4,
variable.names=parameters,
n.iter=nPerChain , thin=thinSteps )
mcmcChain1937_1 = as.matrix( mcmc_samples4)
# Convert precision (tau) to SD###
sigma = 1  / sqrt( mcmcChain1937_1[, "tau" ] )
mcmcChain1937_1 = as.data.frame(cbind( mcmcChain1937_1, sigma ))
g4<-gelman.diag(mcmc_samples4)

######
mcmc_samples5 <- coda.samples(model5,
variable.names=parameters,
n.iter=nPerChain , thin=thinSteps )
mcmcChain2157_1 = as.matrix( mcmc_samples5)
# Convert precision (tau) to SD###
sigma = 1  / sqrt( mcmcChain2157_1[, "tau" ] )
mcmcChain2157_1 = as.data.frame(cbind( mcmcChain2157_1, sigma ))
g5<-gelman.diag(mcmc_samples5)

######
mcmc_samples6 <- coda.samples(model6,
variable.names=parameters,
n.iter=nPerChain , thin=thinSteps )
mcmcChain2208_1 = as.matrix( mcmc_samples6)
# Convert precision (tau) to SD###
sigma = 1  / sqrt( mcmcChain2208_1[, "tau" ] )
mcmcChain2208_1 = as.data.frame(cbind( mcmcChain2208_1, sigma ))
g6<-gelman.diag(mcmc_samples6)


######
mcmc_samples7 <- coda.samples(model7,
variable.names=parameters,
n.iter=nPerChain , thin=thinSteps )
mcmcChain2304_1 = as.matrix( mcmc_samples7)
# Convert precision (tau) to SD###
sigma = 1  / sqrt( mcmcChain2304_1[, "tau" ] )
mcmcChain2304_1 = as.data.frame(cbind( mcmcChain2304_1, sigma ))
g7<-gelman.diag(mcmc_samples7)
######
mcmc_samples8 <- coda.samples(model8,
variable.names=parameters,
n.iter=nPerChain , thin=thinSteps )
mcmcChain2319_1 = as.matrix( mcmc_samples8)
# Convert precision (tau) to SD###
sigma = 1  / sqrt( mcmcChain2319_1[, "tau" ] )
mcmcChain2319_1 = as.data.frame(cbind( mcmcChain2319_1, sigma ))
g8<-gelman.diag(mcmc_samples8)
######
mcmc_samples9 <- coda.samples(model9,
variable.names=parameters,
n.iter=nPerChain , thin=thinSteps )
mcmcChain2320_1 = as.matrix( mcmc_samples9)
# Convert precision (tau) to SD###
sigma = 1  / sqrt( mcmcChain2320_1[, "tau" ] )
mcmcChain2320_1 = as.data.frame(cbind( mcmcChain2320_1, sigma ))
g9<-gelman.diag(mcmc_samples9)
######
mcmc_samples10 <- coda.samples(model10,
variable.names=parameters,
n.iter=nPerChain , thin=thinSteps )
mcmcChain2349_1 = as.matrix( mcmc_samples10)
# Convert precision (tau) to SD###
sigma = 1  / sqrt( mcmcChain2349_1[, "tau" ] )
mcmcChain2349_1 = as.data.frame(cbind( mcmcChain2349_1, sigma ))
g10<-gelman.diag(mcmc_samples10)
######
mcmc_samples11 <- coda.samples(model11,
variable.names=parameters,
n.iter=nPerChain , thin=thinSteps )
mcmcChain2603_1 = as.matrix( mcmc_samples11)
# Convert precision (tau) to SD###
sigma = 1  / sqrt( mcmcChain2603_1[, "tau" ] )
mcmcChain2603_1 = as.data.frame(cbind( mcmcChain2603_1, sigma ))
g11<-gelman.diag(mcmc_samples11)
######
mcmc_samples12 <- coda.samples(model12,
variable.names=parameters,
n.iter=nPerChain , thin=thinSteps )
mcmcChain2629_1 = as.matrix( mcmc_samples12)
# Convert precision (tau) to SD###
sigma = 1  / sqrt( mcmcChain2629_1[, "tau" ] )
mcmcChain2629_1 = as.data.frame(cbind( mcmcChain2629_1, sigma ))
g12<-gelman.diag(mcmc_samples12)
######
mcmc_samples13 <- coda.samples(model13,
variable.names=parameters,
n.iter=nPerChain , thin=thinSteps )
mcmcChain2655_1 = as.matrix( mcmc_samples13)
# Convert precision (tau) to SD###
sigma = 1  / sqrt( mcmcChain2655_1[, "tau" ] )
mcmcChain2655_1 = as.data.frame(cbind( mcmcChain2655_1, sigma ))
g13<-gelman.diag(mcmc_samples13)
######
mcmc_samples14 <- coda.samples(model14,
variable.names=parameters,
n.iter=nPerChain , thin=thinSteps )
mcmcChain2696_1 = as.matrix( mcmc_samples14)
# Convert precision (tau) to SD###
sigma = 1  / sqrt( mcmcChain2696_1[, "tau" ] )
mcmcChain2696_1 = as.data.frame(cbind( mcmcChain2696_1, sigma ))
g14<-gelman.diag(mcmc_samples14)
######
mcmc_samples15 <- coda.samples(model15,
variable.names=parameters,
n.iter=nPerChain , thin=thinSteps )
mcmcChain2712_1 = as.matrix( mcmc_samples15)
# Convert precision (tau) to SD###
sigma = 1  / sqrt( mcmcChain2712_1[, "tau" ] )
mcmcChain2712_1 = as.data.frame(cbind( mcmcChain2712_1, sigma ))
g15<-gelman.diag(mcmc_samples15)
######
mcmc_samples16 <- coda.samples(model16,
variable.names=parameters,
n.iter=nPerChain , thin=thinSteps )
mcmcChain2723_1 = as.matrix( mcmc_samples16)
# Convert precision (tau) to SD###
sigma = 1  / sqrt( mcmcChain2723_1[, "tau" ] )
mcmcChain2723_1 = as.data.frame(cbind( mcmcChain2723_1, sigma ))
g16<-gelman.diag(mcmc_samples16)
######
mcmc_samples17 <- coda.samples(model17,
variable.names=parameters,
n.iter=nPerChain , thin=thinSteps )
mcmcChain2737_1 = as.matrix( mcmc_samples17)
# Convert precision (tau) to SD###
sigma = 1  / sqrt( mcmcChain2737_1[, "tau" ] )
mcmcChain2737_1 = as.data.frame(cbind( mcmcChain2737_1, sigma ))
g17<-gelman.diag(mcmc_samples17)
######
mcmc_samples18 <- coda.samples(model18,
variable.names=parameters,
n.iter=nPerChain , thin=thinSteps )
mcmcChain2760_1 = as.matrix( mcmc_samples18)
# Convert precision (tau) to SD###
sigma = 1  / sqrt( mcmcChain2760_1[, "tau" ] )
mcmcChain2760_1 = as.data.frame(cbind( mcmcChain2760_1, sigma ))
g18<-gelman.diag(mcmc_samples18)


##########################################

##########################################

######### Gelman potential scale reduction factor summary for all par and multivariate #############
gelmin<-c(min(g1$psrf[,1]), min(g2$psrf[,1]),min(g3$psrf[,1]), min(g4$psrf[,1]),
min(g5$psrf[,1]), min(g6$psrf[,1]),min(g7$psrf[,1]), min(g8$psrf[,1]),
min(g9$psrf[,1]), min(g10$psrf[,1]),min(g11$psrf[,1]), min(g12$psrf[,1]),
min(g13$psrf[,1]), min(g14$psrf[,1]),min(g15$psrf[,1]), min(g16$psrf[,1]),
min(g17$psrf[,1]), min(g18$psrf[,1]))


gelmax<-c(max(g1$psrf[,1]), max(g2$psrf[,1]),max(g3$psrf[,1]), max(g4$psrf[,1]),
max(g5$psrf[,1]), max(g6$psrf[,1]),max(g7$psrf[,1]), max(g8$psrf[,1]),
max(g9$psrf[,1]), max(g10$psrf[,1]),max(g11$psrf[,1]), max(g12$psrf[,1]),
max(g13$psrf[,1]), max(g14$psrf[,1]),max(g15$psrf[,1]), max(g16$psrf[,1]),
max(g17$psrf[,1]), max(g18$psrf[,1]))

gelmulit<-c(g1$mpsrf, g2$mpsrf, g3$mpsrf, g4$mpsrf,g5$mpsrf, g6$mpsrf,g7$mpsrf, g8$mpsrf,
g9$mpsrf, g10$mpsrf, g11$mpsrf, g12$mpsrf,g13$mpsrf, g14$mpsrf,g15$mpsrf, g16$mpsrf,
g17$mpsrf, g18$mpsrf)

gelsum<-as.data.frame(cbind(gelmin, gelmax,gelmulit))
colnames(gelsum)<-c("gel_min", "gel_max","gel_multi")



##########################################
#########  SAVE Samples ####################
# ID 1845 1894 1902 1937 2157 2208 2304 2319 2320 2349 2603 2629 2655 2696 2712 2723 2737 2760
###
print("writing samples")


setwd("~/Documents/ACiP/Post_Data")
write.table(mcmcChain1845_1,"mcmcChain1845_1", sep="\t", col.name=TRUE)
write.table(mcmcChain1894_1,"mcmcChain1894_1", sep="\t", col.name=TRUE)
write.table(mcmcChain1902_1,"mcmcChain1902_1", sep="\t", col.name=TRUE)
write.table(mcmcChain1937_1,"mcmcChain1937_1", sep="\t", col.name=TRUE)
write.table(mcmcChain2157_1,"mcmcChain2157_1", sep="\t", col.name=TRUE)
write.table(mcmcChain2208_1,"mcmcChain2208_1", sep="\t", col.name=TRUE)
write.table(mcmcChain2304_1,"mcmcChain2304_1", sep="\t", col.name=TRUE)
write.table(mcmcChain2319_1,"mcmcChain2319_1", sep="\t", col.name=TRUE)
write.table(mcmcChain2320_1,"mcmcChain2320_1", sep="\t", col.name=TRUE)
write.table(mcmcChain2349_1,"mcmcChain2349_1", sep="\t", col.name=TRUE)
write.table(mcmcChain2603_1,"mcmcChain2603_1", sep="\t", col.name=TRUE)
write.table(mcmcChain2629_1,"mcmcChain2629_1", sep="\t", col.name=TRUE)
write.table(mcmcChain2655_1,"mcmcChain2655_1", sep="\t", col.name=TRUE)
write.table(mcmcChain2696_1,"mcmcChain2696_1", sep="\t", col.name=TRUE)
write.table(mcmcChain2712_1,"mcmcChain2712_1", sep="\t", col.name=TRUE)
write.table(mcmcChain2723_1,"mcmcChain2723_1", sep="\t", col.name=TRUE)
write.table(mcmcChain2737_1,"mcmcChain2737_1", sep="\t", col.name=TRUE)
write.table(mcmcChain2760_1,"mcmcChain2760_1", sep="\t", col.name=TRUE)

#########################################################################

#################    PARAMETER ESTIMATES      ###########################

#########################################################################

#########################################################################

EstPars1<-apply(mcmcChain1845_1, 2, median)
EstPars2<-apply(mcmcChain1894_1, 2, median)
EstPars3<-apply(mcmcChain1902_1, 2, median)
EstPars4<-apply(mcmcChain1937_1, 2, median)
EstPars5<-apply(mcmcChain2157_1, 2, median)
EstPars6<-apply(mcmcChain2208_1, 2, median)
EstPars7<-apply(mcmcChain2304_1, 2, median)
EstPars8<-apply(mcmcChain2319_1, 2, median)
EstPars9<-apply(mcmcChain2320_1, 2, median)
EstPars10<-apply(mcmcChain2349_1, 2, median)
EstPars11<-apply(mcmcChain2603_1, 2, median)
EstPars12<-apply(mcmcChain2629_1, 2, median)
EstPars13<-apply(mcmcChain2655_1, 2, median)
EstPars14<-apply(mcmcChain2696_1, 2, median)
EstPars15<-apply(mcmcChain2712_1, 2, median)
EstPars16<-apply(mcmcChain2723_1, 2, median)
EstPars17<-apply(mcmcChain2737_1, 2, median)
EstPars18<-apply(mcmcChain2760_1, 2, median)

ParEst<-rbind(EstPars1,EstPars2,EstPars3,EstPars4,EstPars5,EstPars6,EstPars7, 
              EstPars8,EstPars9,EstPars10,EstPars11,EstPars12,EstPars13,EstPars14, 
              EstPars15,EstPars16,EstPars17,EstPars18 )

###### quantiles for all Pars

qf <- function(d) {
    (quantile(d, c(.05,.5,.95)))
}

qEstPars1<-t(apply(mcmcChain1845_1, 2, qf))
qEstPars2<-t(apply(mcmcChain1894_1, 2, qf))
qEstPars3<-t(apply(mcmcChain1902_1, 2, qf))
qEstPars4<-t(apply(mcmcChain1937_1, 2, qf))
qEstPars5<-t(apply(mcmcChain2157_1, 2, qf))
qEstPars6<-t(apply(mcmcChain2208_1, 2, qf))
qEstPars7<-t(apply(mcmcChain2304_1, 2, qf))
qEstPars8<-t(apply(mcmcChain2319_1, 2, qf))
qEstPars9<-t(apply(mcmcChain2320_1, 2, qf))
qEstPars10<-t(apply(mcmcChain2349_1, 2, qf))
qEstPars11<-t(apply(mcmcChain2603_1, 2, qf))
qEstPars12<-t(apply(mcmcChain2629_1, 2, qf))
qEstPars13<-t(apply(mcmcChain2655_1, 2, qf))
qEstPars14<-t(apply(mcmcChain2696_1, 2, qf))
qEstPars15<-t(apply(mcmcChain2712_1, 2, qf))
qEstPars16<-t(apply(mcmcChain2723_1, 2, qf))
qEstPars17<-t(apply(mcmcChain2737_1, 2, qf))
qEstPars18<-t(apply(mcmcChain2760_1, 2, qf))

qEst<-cbind(qEstPars1,qEstPars2,qEstPars3,qEstPars4,qEstPars5,qEstPars6,qEstPars7,
qEstPars8,qEstPars9,qEstPars10,qEstPars11,qEstPars12,qEstPars13,qEstPars14,
qEstPars15,qEstPars16,qEstPars17,qEstPars18 )

rID<-rep(ID, rep(3,18))
qs<-rep(c("_05","_50","_95"), 18)
clnm<-paste("ID",rID,qs, sep = "" )
colnames(qEst)<-clnm


print("writing quantiles")
write.table(qEst,"Par_quantiles_1", sep="\t", col.name=TRUE)




#########################################################################

#########################################################################

#########################################################################

#########################################################################

#################   Calc DIC with pD criteria  ###########################

#########################################################################

#########################################################################
# ID 1845 1894 1902 1937 2157 2208 2304 2319 2320 2349 2603 2629 2655 2696 2712 2723 2737 2760

print("sampling for DIC")
dic1 <- dic.samples(model1, DICsteps, "pD")
dic1a <- dic.samples(model1, DICsteps, "popt")
dic1845_1<-c(sum(dic1$deviance), sum(dic1$penalty),sum(dic1$deviance, dic1$penalty),
             sum(dic1a$deviance), sum(dic1a$penalty), sum(dic1a$deviance, dic1a$penalty))

dic1 <- dic.samples(model2, DICsteps, "pD")
dic1a <- dic.samples(model2, DICsteps, "popt")
dic1894_1<-c(sum(dic1$deviance), sum(dic1$penalty),sum(dic1$deviance, dic1$penalty),
             sum(dic1a$deviance), sum(dic1a$penalty), sum(dic1a$deviance, dic1a$penalty))

dic1 <- dic.samples(model3, DICsteps, "pD")
dic1a <- dic.samples(model3, DICsteps, "popt")
dic1902_1<-c(sum(dic1$deviance), sum(dic1$penalty),sum(dic1$deviance, dic1$penalty),
             sum(dic1a$deviance), sum(dic1a$penalty), sum(dic1a$deviance, dic1a$penalty))

dic1 <- dic.samples(model4, DICsteps, "pD")
dic1a <- dic.samples(model4, DICsteps, "popt")
dic1937_1<-c(sum(dic1$deviance), sum(dic1$penalty),sum(dic1$deviance, dic1$penalty),
             sum(dic1a$deviance), sum(dic1a$penalty), sum(dic1a$deviance, dic1a$penalty))

dic1 <- dic.samples(model5, DICsteps, "pD")
dic1a <- dic.samples(model5, DICsteps, "popt")
dic2157_1<-c(sum(dic1$deviance), sum(dic1$penalty),sum(dic1$deviance, dic1$penalty),
             sum(dic1a$deviance), sum(dic1a$penalty), sum(dic1a$deviance, dic1a$penalty))

dic1 <- dic.samples(model6, DICsteps, "pD")
dic1a <- dic.samples(model6, DICsteps, "popt")
dic2208_1<-c(sum(dic1$deviance), sum(dic1$penalty),sum(dic1$deviance, dic1$penalty),
             sum(dic1a$deviance), sum(dic1a$penalty), sum(dic1a$deviance, dic1a$penalty))

dic1 <- dic.samples(model7, DICsteps, "pD")
dic1a <- dic.samples(model7, DICsteps, "popt")
dic2304_1<-c(sum(dic1$deviance), sum(dic1$penalty),sum(dic1$deviance, dic1$penalty),
             sum(dic1a$deviance), sum(dic1a$penalty), sum(dic1a$deviance, dic1a$penalty))

dic1 <- dic.samples(model8, DICsteps, "pD")
dic1a <- dic.samples(model8, DICsteps, "popt")
dic2319_1<-c(sum(dic1$deviance), sum(dic1$penalty),sum(dic1$deviance, dic1$penalty),
             sum(dic1a$deviance), sum(dic1a$penalty), sum(dic1a$deviance, dic1a$penalty))

dic1 <- dic.samples(model9, DICsteps, "pD")
dic1a <- dic.samples(model9, DICsteps, "popt")
dic2320_1<-c(sum(dic1$deviance), sum(dic1$penalty),sum(dic1$deviance, dic1$penalty),
             sum(dic1a$deviance), sum(dic1a$penalty), sum(dic1a$deviance, dic1a$penalty))

dic1 <- dic.samples(model10, DICsteps, "pD")
dic1a <- dic.samples(model10, DICsteps, "popt")
dic2349_1<-c(sum(dic1$deviance), sum(dic1$penalty),sum(dic1$deviance, dic1$penalty),
             sum(dic1a$deviance), sum(dic1a$penalty), sum(dic1a$deviance, dic1a$penalty))

dic1 <- dic.samples(model11, DICsteps, "pD")
dic1a <- dic.samples(model11, DICsteps, "popt")
dic2603_1<-c(sum(dic1$deviance), sum(dic1$penalty),sum(dic1$deviance, dic1$penalty),
             sum(dic1a$deviance), sum(dic1a$penalty), sum(dic1a$deviance, dic1a$penalty))

dic1 <- dic.samples(model12, DICsteps, "pD")
dic1a <- dic.samples(model12, DICsteps, "popt")
dic2629_1<-c(sum(dic1$deviance), sum(dic1$penalty),sum(dic1$deviance, dic1$penalty),
             sum(dic1a$deviance), sum(dic1a$penalty), sum(dic1a$deviance, dic1a$penalty))

dic1 <- dic.samples(model13, DICsteps, "pD")
dic1a <- dic.samples(model13, DICsteps, "popt")
dic2655_1<-c(sum(dic1$deviance), sum(dic1$penalty),sum(dic1$deviance, dic1$penalty),
             sum(dic1a$deviance), sum(dic1a$penalty), sum(dic1a$deviance, dic1a$penalty))

dic1 <- dic.samples(model14, DICsteps, "pD")
dic1a <- dic.samples(model14, DICsteps, "popt")
dic2696_1<-c(sum(dic1$deviance), sum(dic1$penalty),sum(dic1$deviance, dic1$penalty),
             sum(dic1a$deviance), sum(dic1a$penalty), sum(dic1a$deviance, dic1a$penalty))

dic1 <- dic.samples(model15, DICsteps, "pD")
dic1a <- dic.samples(model15, DICsteps, "popt")
dic2712_1<-c(sum(dic1$deviance), sum(dic1$penalty),sum(dic1$deviance, dic1$penalty),
             sum(dic1a$deviance), sum(dic1a$penalty), sum(dic1a$deviance, dic1a$penalty))

dic1 <- dic.samples(model16, DICsteps, "pD")
dic1a <- dic.samples(model16, DICsteps, "popt")
dic2723_1<-c(sum(dic1$deviance), sum(dic1$penalty),sum(dic1$deviance, dic1$penalty),
             sum(dic1a$deviance), sum(dic1a$penalty), sum(dic1a$deviance, dic1a$penalty))

dic1 <- dic.samples(model17, DICsteps, "pD")
dic1a <- dic.samples(model17, DICsteps, "popt")
dic2737_1<-c(sum(dic1$deviance), sum(dic1$penalty),sum(dic1$deviance, dic1$penalty),
             sum(dic1a$deviance), sum(dic1a$penalty), sum(dic1a$deviance, dic1a$penalty))

dic1 <- dic.samples(model18, DICsteps, "pD")
dic1a <- dic.samples(model18, DICsteps, "popt")
dic2760_1<-c(sum(dic1$deviance), sum(dic1$penalty),sum(dic1$deviance, dic1$penalty),
             sum(dic1a$deviance), sum(dic1a$penalty), sum(dic1a$deviance, dic1a$penalty))

dics<-rbind(dic1845_1, dic1894_1, dic1902_1, dic1937_1, dic2157_1, 
            dic2208_1, dic2304_1, dic2319_1, dic2320_1, dic2349_1, 
            dic2603_1, dic2629_1, dic2655_1, dic2696_1, dic2712_1, 
            dic2723_1, dic2737_1, dic2760_1)
colnames(dics)<-c("pD_dev", "pD_pen", "pD_DIC", "popt_dev", "popt_pen", "popt_DIC")



###################################################

###################################################

#   FINAL TABLE W/ median estimates and DICS ########

###################################################
modelname<-rep("CiCc_Jf",18)

print("writing PAR Estimates")
ParEstFull<-as.data.frame(cbind(modelname, ID, IDgeno, ParEst,dics, gelsum))
write.table(ParEstFull,"ParEstFull_model_1", sep="\t", col.name=TRUE)




#############   SAVE TRACE PLOTS
setwd("~/Documents/ACiP/Trace_plots")
require(ggmcmc)
pdf("density_1845_1", colormodel='cmyk',width=6.25, height=(10))
ggs_density(ggs(mcmc_samples1))
dev.off()
pdf("density_1894_1", colormodel='cmyk',width=6.25, height=(10))
ggs_density(ggs(mcmc_samples2))
dev.off()
pdf("density_1902_1", colormodel='cmyk',width=6.25, height=(10))
ggs_density(ggs(mcmc_samples3))
dev.off()
pdf("density_1937_1", colormodel='cmyk',width=6.25, height=(10))
ggs_density(ggs(mcmc_samples4))
dev.off()
pdf("density_2157_1", colormodel='cmyk',width=6.25, height=(10))
ggs_density(ggs(mcmc_samples5))
dev.off()
pdf("density_2208_1", colormodel='cmyk',width=6.25, height=(10))
ggs_density(ggs(mcmc_samples6))
dev.off()
pdf("density_2304_1", colormodel='cmyk',width=6.25, height=(10))
ggs_density(ggs(mcmc_samples7))
dev.off()
pdf("density_2319_1", colormodel='cmyk',width=6.25, height=(10))
ggs_density(ggs(mcmc_samples8))
dev.off()
pdf("density_2320_1", colormodel='cmyk',width=6.25, height=(10))
ggs_density(ggs(mcmc_samples9))
dev.off()
pdf("density_2349_1", colormodel='cmyk',width=6.25, height=(10))
ggs_density(ggs(mcmc_samples10))
dev.off()
pdf("density_2603_1", colormodel='cmyk',width=6.25, height=(10))
ggs_density(ggs(mcmc_samples11))
dev.off()
pdf("density_2629_1", colormodel='cmyk',width=6.25, height=(10))
ggs_density(ggs(mcmc_samples12))
dev.off()
pdf("density_2655_1", colormodel='cmyk',width=6.25, height=(10))
ggs_density(ggs(mcmc_samples13))
dev.off()
pdf("density_2696_1", colormodel='cmyk',width=6.25, height=(10))
ggs_density(ggs(mcmc_samples14))
dev.off()
pdf("density_2712_1", colormodel='cmyk',width=6.25, height=(10))
ggs_density(ggs(mcmc_samples15))
dev.off()
pdf("density_2723_1", colormodel='cmyk',width=6.25, height=(10))
ggs_density(ggs(mcmc_samples16))
dev.off()
pdf("density_2737_1", colormodel='cmyk',width=6.25, height=(10))
ggs_density(ggs(mcmc_samples17))
dev.off()
pdf("density_2760_1", colormodel='cmyk',width=6.25, height=(10))
ggs_density(ggs(mcmc_samples18))
dev.off()



pdf("traceplot_1845_1", colormodel='cmyk',width=6.25, height=(10))
ggs_traceplot(ggs(mcmc_samples1))
dev.off()
pdf("traceplot_1894_1", colormodel='cmyk',width=6.25, height=(10))
ggs_traceplot(ggs(mcmc_samples2))
dev.off()
pdf("traceplot_1902_1", colormodel='cmyk',width=6.25, height=(10))
ggs_traceplot(ggs(mcmc_samples3))
dev.off()
pdf("traceplot_1937_1", colormodel='cmyk',width=6.25, height=(10))
ggs_traceplot(ggs(mcmc_samples4))
dev.off()
pdf("traceplot_2157_1", colormodel='cmyk',width=6.25, height=(10))
ggs_traceplot(ggs(mcmc_samples5))
dev.off()
pdf("traceplot_2208_1", colormodel='cmyk',width=6.25, height=(10))
ggs_traceplot(ggs(mcmc_samples6))
dev.off()
pdf("traceplot_2304_1", colormodel='cmyk',width=6.25, height=(10))
ggs_traceplot(ggs(mcmc_samples7))
dev.off()
pdf("traceplot_2319_1", colormodel='cmyk',width=6.25, height=(10))
ggs_traceplot(ggs(mcmc_samples8))
dev.off()
pdf("traceplot_2320_1", colormodel='cmyk',width=6.25, height=(10))
ggs_traceplot(ggs(mcmc_samples9))
dev.off()
pdf("traceplot_2349_1", colormodel='cmyk',width=6.25, height=(10))
ggs_traceplot(ggs(mcmc_samples10))
dev.off()
pdf("traceplot_2603_1", colormodel='cmyk',width=6.25, height=(10))
ggs_traceplot(ggs(mcmc_samples11))
dev.off()
pdf("traceplot_2629_1", colormodel='cmyk',width=6.25, height=(10))
ggs_traceplot(ggs(mcmc_samples12))
dev.off()
pdf("traceplot_2655_1", colormodel='cmyk',width=6.25, height=(10))
ggs_traceplot(ggs(mcmc_samples13))
dev.off()
pdf("traceplot_2696_1", colormodel='cmyk',width=6.25, height=(10))
ggs_traceplot(ggs(mcmc_samples14))
dev.off()
pdf("traceplot_2712_1", colormodel='cmyk',width=6.25, height=(10))
ggs_traceplot(ggs(mcmc_samples15))
dev.off()
pdf("traceplot_2723_1", colormodel='cmyk',width=6.25, height=(10))
ggs_traceplot(ggs(mcmc_samples16))
dev.off()
pdf("traceplot_2737_1", colormodel='cmyk',width=6.25, height=(10))
ggs_traceplot(ggs(mcmc_samples17))
dev.off()
pdf("traceplot_2760_1", colormodel='cmyk',width=6.25, height=(10))
ggs_traceplot(ggs(mcmc_samples18))
dev.off()



pdf("crosscorrelation_1845_1", colormodel='cmyk',width=6.25, height=(10))
ggs_crosscorrelation(ggs(mcmc_samples1))
dev.off()
pdf("crosscorrelation_1894_1", colormodel='cmyk',width=6.25, height=(10))
ggs_crosscorrelation(ggs(mcmc_samples2))
dev.off()
pdf("crosscorrelation_1902_1", colormodel='cmyk',width=6.25, height=(10))
ggs_crosscorrelation(ggs(mcmc_samples3))
dev.off()
pdf("crosscorrelation_1937_1", colormodel='cmyk',width=6.25, height=(10))
ggs_crosscorrelation(ggs(mcmc_samples4))
dev.off()
pdf("crosscorrelation_2157_1", colormodel='cmyk',width=6.25, height=(10))
ggs_crosscorrelation(ggs(mcmc_samples5))
dev.off()
pdf("crosscorrelation_2208_1", colormodel='cmyk',width=6.25, height=(10))
ggs_crosscorrelation(ggs(mcmc_samples6))
dev.off()
pdf("crosscorrelation_2304_1", colormodel='cmyk',width=6.25, height=(10))
ggs_crosscorrelation(ggs(mcmc_samples7))
dev.off()
pdf("crosscorrelation_2319_1", colormodel='cmyk',width=6.25, height=(10))
ggs_crosscorrelation(ggs(mcmc_samples8))
dev.off()
pdf("crosscorrelation_2320_1", colormodel='cmyk',width=6.25, height=(10))
ggs_crosscorrelation(ggs(mcmc_samples9))
dev.off()
pdf("crosscorrelation_2349_1", colormodel='cmyk',width=6.25, height=(10))
ggs_crosscorrelation(ggs(mcmc_samples10))
dev.off()
pdf("crosscorrelation_2603_1", colormodel='cmyk',width=6.25, height=(10))
ggs_crosscorrelation(ggs(mcmc_samples11))
dev.off()
pdf("crosscorrelation_2629_1", colormodel='cmyk',width=6.25, height=(10))
ggs_crosscorrelation(ggs(mcmc_samples12))
dev.off()
pdf("crosscorrelation_2655_1", colormodel='cmyk',width=6.25, height=(10))
ggs_crosscorrelation(ggs(mcmc_samples13))
dev.off()
pdf("crosscorrelation_2696_1", colormodel='cmyk',width=6.25, height=(10))
ggs_crosscorrelation(ggs(mcmc_samples14))
dev.off()
pdf("crosscorrelation_2712_1", colormodel='cmyk',width=6.25, height=(10))
ggs_crosscorrelation(ggs(mcmc_samples15))
dev.off()
pdf("crosscorrelation_2723_1", colormodel='cmyk',width=6.25, height=(10))
ggs_crosscorrelation(ggs(mcmc_samples16))
dev.off()
pdf("crosscorrelation_2737_1", colormodel='cmyk',width=6.25, height=(10))
ggs_crosscorrelation(ggs(mcmc_samples17))
dev.off()
pdf("crosscorrelation_2760_1", colormodel='cmyk',width=6.25, height=(10))
ggs_crosscorrelation(ggs(mcmc_samples18))
dev.off()


pdf("Rhat_1845_1", colormodel='cmyk',width=6.25, height=(10))
ggs_Rhat(ggs(mcmc_samples1)) + xlab("R_hat")
dev.off()
pdf("Rhat_1894_1", colormodel='cmyk',width=6.25, height=(10))
ggs_Rhat(ggs(mcmc_samples2))+ xlab("R_hat")
dev.off()
pdf("Rhat_1902_1", colormodel='cmyk',width=6.25, height=(10))
ggs_Rhat(ggs(mcmc_samples3))+ xlab("R_hat")
dev.off()
pdf("Rhat_1937_1", colormodel='cmyk',width=6.25, height=(10))
ggs_Rhat(ggs(mcmc_samples4))
dev.off()
pdf("Rhat_2157_1", colormodel='cmyk',width=6.25, height=(10))
ggs_Rhat(ggs(mcmc_samples5))+ xlab("R_hat")
dev.off()
pdf("Rhat_2208_1", colormodel='cmyk',width=6.25, height=(10))
ggs_Rhat(ggs(mcmc_samples6))+ xlab("R_hat")
dev.off()
pdf("Rhat_2304_1", colormodel='cmyk',width=6.25, height=(10))
ggs_Rhat(ggs(mcmc_samples7))+ xlab("R_hat")
dev.off()
pdf("Rhat_2319_1", colormodel='cmyk',width=6.25, height=(10))
ggs_Rhat(ggs(mcmc_samples8))+ xlab("R_hat")
dev.off()
pdf("Rhat_2320_1", colormodel='cmyk',width=6.25, height=(10))
ggs_Rhat(ggs(mcmc_samples9))+ xlab("R_hat")
dev.off()
pdf("Rhat_2349_1", colormodel='cmyk',width=6.25, height=(10))
ggs_Rhat(ggs(mcmc_samples10))+ xlab("R_hat")
dev.off()
pdf("Rhat_2603_1", colormodel='cmyk',width=6.25, height=(10))
ggs_Rhat(ggs(mcmc_samples11))+ xlab("R_hat")
dev.off()
pdf("Rhat_2629_1", colormodel='cmyk',width=6.25, height=(10))
ggs_Rhat(ggs(mcmc_samples12))+ xlab("R_hat")
dev.off()
pdf("Rhat_2655_1", colormodel='cmyk',width=6.25, height=(10))
ggs_Rhat(ggs(mcmc_samples13))+ xlab("R_hat")
dev.off()
pdf("Rhat_2696_1", colormodel='cmyk',width=6.25, height=(10))
ggs_Rhat(ggs(mcmc_samples14))+ xlab("R_hat")
dev.off()
pdf("Rhat_2712_1", colormodel='cmyk',width=6.25, height=(10))
ggs_Rhat(ggs(mcmc_samples15))+ xlab("R_hat")
dev.off()
pdf("Rhat_2723_1", colormodel='cmyk',width=6.25, height=(10))
ggs_Rhat(ggs(mcmc_samples16))+ xlab("R_hat")
dev.off()
pdf("Rhat_2737_1", colormodel='cmyk',width=6.25, height=(10))
ggs_Rhat(ggs(mcmc_samples17))
dev.off()
pdf("Rhat_2760_1", colormodel='cmyk',width=6.25, height=(10))
ggs_Rhat(ggs(mcmc_samples18))+ xlab("R_hat")
dev.off()














print("Finished")

















