###############################
###      Bayesian Model    ####
### Following Farq PS model ####
#### Testing Mesophyll, ETR from Chl Flr
######  & no temp dependency ###
###############################

CiCc_Jf <- "model {
for (i in 1:N){
An[i] ~ dnorm( mu.A[i] , tau )
mu.A[i] <- min(Ac[i], Aj[i])

#quadratic solution for net A if limited by Rubisco
a1[i]<-(-1/gm)
b1[i]<-((Vcmax-Rd)/gm)+CiP[i]+(Kc*((1+O[i])/Ko ))
c1[i]<-Rd*(CiP[i]+(Kc*((1+O[i])/Ko )))-Vcmax*(CiP[i]-gammaS)
bac1[i]<-(b1[i]^2)-(4*a1[i]*c1[i])
Ac[i]<- (-b1[i]+sqrt(bac1[i]))/(2*a1[i])

# quadratic solution for net A if limited by light (RuBP regeneration)
a2[i]<-(-1.0/gm)
b2[i]<-(((Jf[i]/4)-Rd)/gm) + CiP[i]+2*gammaS
c2[i]<-(Rd*(CiP[i]+2*gammaS)) -((Jf[i]/4)*(CiP[i]-gammaS))
bac2[i]<-(b2[i]^2)-(4*a2[i]*c2[i])
Aj[i]<- (-b2[i]+sqrt(bac2[i]))/(2*a2[i])


}
### nromal prior from Agricultural Species WULLSCHLEGER
Vcmax ~ dnorm(90, 0.000625)T(0,100000)
### normal prior Dark Respiration from LR data (umol m-2 s-1)
Rd ~ dunif(0, 25)
### normal prior centered on 0  for mesophyl cond (umol m-2 s-1 Pa-1)
gm ~ dnorm(2.5, 0.0025)T(0,142)
### prior on gamma star (Pa)
gammaS ~ dnorm(3.86,0.0625)T(0.24,7.24)
Kc ~ dnorm(27.24,0.01)T(0,100)
Ko ~ dnorm(30400, 0.00000004)T(0,65400)
### flat prior on precision
tau ~ dgamma(.001 , .001)
}
" 