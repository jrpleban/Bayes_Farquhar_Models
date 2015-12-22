###############################
###      Bayesian Model    ####
### Following Farq PS model
#### with TPU, no temp dependency ###
###############################

CaCc_Jm <- "model {
for (i in 1:N){
An[i] ~ dnorm( mu.A[i] , tau )
mu.A[i] <- min(Ac[i], Aj[i])

#### electron transport
Ji[i] <- (Q[i] * phiJ*.85)
bb[i] <- -Ji[i] -Jmax
cc[i] <- Ji[i]*Jmax
bac[i]<-(bb[i]^2)-(4*thetaJ*cc[i])
Jm[i]<- (-bb[i]-sqrt(bac[i]))/(2*thetaJ)

#quadratic solution for net A if limited by Rubisco
a1[i]<-(-1/g[i])
b1[i]<-((Vcmax-Rd)/g[i])+(Ca[i]+(Kc*((1+O[i])/Ko )))
c1[i]<-Rd*(Ca[i]+(Kc*((1+O[i])/Ko )))-Vcmax*(Ca[i]-gammaS)
bac1[i]<-(b1[i]^2)-(4*a1[i]*c1[i])
Ac[i]<- (-b1[i]+sqrt(bac1[i]))/(2*a1[i])

# quadratic solution for net A if limited by light (RuBP regeneration)
a2[i]<-(-1.0/g[i])
b2[i]<-(((Jm[i]/4)-Rd)/g[i]) + Ca[i]+2*gammaS
c2[i]<-(Rd*(Ca[i]+2*gammaS)) -((Jm[i]/4)*(Ca[i]-gammaS))
bac2[i]<-(b2[i]^2)-(4*a2[i]*c2[i])
Aj[i]<- (-b2[i]+sqrt(bac2[i]))/(2*a2[i])


}
### nromal prior from Agricultural Species WULLSCHLEGER
Vcmax ~ dnorm(90, 0.000625)T(0,100000)
### normal prior Dark Respiration from LR data (umol m-2 s-1)
Rd ~ dunif(0, 25)
### prior on gamma star (Pa)
gammaS ~ dnorm(3.86,0.0625)T(0,142)
# curvature parameter on LR curve (thetaJ)
thetaJ ~ dunif(0.1, 0.9999)
### nromal prior from Agricultural Species WULLSCHLEGER
Jmax ~ dnorm(171,0.000308)T(0,100000)
## quantum yield phiJ from theoretical limits
phiJ ~ dunif(0,.5)
Kc ~ dnorm(27.24,0.01)T(0, 100)
Ko ~ dnorm(30400, 0.00000004)T(0,65400)
### flat prior on precision
tau ~ dgamma(.001 , .001)
}
" 
