###############################
###      Bayesian Model    ####
### Following Farq PS model#### 
###############################
CiCc_Temp_Jf <- "model {
for (i in 1:N)
{
An[i] ~ dnorm( mu.A[i] , tau )
mu.A[i] <- min(Ac[i], Aj[i])

### Temp Depedencies on Pars (Arrhenius Temp Functions)###
#  Arrhenius Temp function 
#                                          ( Ee (Tobs - 298))
##     Y = f(Y25, Ee, Tobs) = Y25 * exp (-----------------------)
#                                          ( 298 * R * Tobs)

gammaS[i] <- gammaS25*exp((T[i]-Tref)*EgammaS/(Kref*R*K[i]))
Rd[i] <- Rd25* exp((T[i]-Tref)*ERd/(Kref*R*K[i]))
Kc[i] <- Kc25*exp((T[i]-Tref)*EKc/(Kref*R*K[i]))
Ko[i] <- Ko25*exp((T[i]-Tref)*EKo/(Kref*R*K[i]))
Vcmax[i] <- Vcmax25*exp((T[i]-Tref)*EVcmax/(Kref*R*K[i]))
gm[i] <- gm25*exp((T[i]-Tref)*Egm/(Kref*R*K[i]))

#quadratic solution for net A if limited by Rubisco
a1[i]<-(-1/gm[i])
b1[i]<-((Vcmax[i]-Rd[i])/gm[i])+(CiP[i]+(Kc[i]*((1+O[i])/Ko[i] )))
c1[i]<-Rd[i]*(CiP[i]+(Kc[i]*((1+O[i])/Ko[i] )))-Vcmax[i]*(CiP[i]-gammaS[i])
bac1[i]<-(b1[i]^2)-(4*a1[i]*c1[i])
Ac[i]<- (-b1[i]+sqrt(bac1[i]))/(2*a1[i])

# quadratic solution for net A if limited by light (RuBP regeneration)
a2[i]<-(-1.0/gm[i])
b2[i]<-(((Jf[i]/4)-Rd[i])/gm[i]) + CiP[i]+2*gammaS[i]
c2[i]<-(Rd[i]*(CiP[i]+2*gammaS[i])) -((Jf[i]/4)*(CiP[i]-gammaS[i]))
bac2[i]<-(b2[i]^2)-(4*a2[i]*c2[i])
Aj[i]<- (-b2[i]+sqrt(bac2[i]))/(2*a2[i])

}
### Activation Energy Prior for Arrhenius Temp Functions
#### dnorm ( mean, precision)  precision = 1/(sd^2)
EgammaS ~ dnorm(26.8, 0.01)
ERd ~ dnorm(63.9, 0.0025)
EKc ~ dnorm(70.4,0.0025)
EKo ~ dnorm(36.0,0.0025)
EVcmax ~ dnorm(65.4, 0.0025)
Egm ~ dnorm(49.6, 0.001)


### normal prior from Agricultural Species WULLSCHLEGER (1993)
Vcmax25 ~ dnorm(90, 0.000625)T(0,100000)  ## Truncated at extremes to ease computation
### normal prior Dark Respiration (umol m-2 s-1)
Rd25 ~ dnorm (1.17, 2.5)
### normal prior centered on 0  for mesophyl cond (umol m-2 s-1 Pa-1)
gm25 ~ dnorm(2.5,0.0025)T(0,142)
### prior on gamma star (Pa)
gammaS25 ~ dnorm(3.74, 4)T(.24,7.23)
### normal prior for Kc (Pa)
Kc25 ~ dnorm(27.24,0.01)T(0,100)
###  prior for Ko(Pa)
Ko25 ~ dnorm(30400, 0.00000004)T(0,65400)

### prior on precision
tau ~ dgamma(.001 , .001)}
" 

