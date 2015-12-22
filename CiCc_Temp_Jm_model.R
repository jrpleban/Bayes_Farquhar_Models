###############################
###      Bayesian Model    ####
### Following Farq PS model#### 
###############################
CiCc_Temp_Jm <- "model {
for (i in 1:N)
{
An[i] ~ dnorm( mu.A[i] , tau )
mu.A[i] <- min(Ac[i], Aj[i])

### Temp Depedencies on Pars (Arrhenius Temp Functions)###
#  Arrhenius Temp function 
#                                          ( Ee (Tobs - 298))
##     Y = f(Y25, Ee, Tobs) = Y25 * exp (-----------------------)
#                                          ( 298 * R * Tobs)

farQ_CiCc_Jm_temp.R(cbind(),)

}
### Activation Energy Prior for Arrhenius Temp Functions
#### dnorm ( mean, precision)  precision = 1/(sd^2)
EgammaS ~ dnorm(26.8, 0.01)
ERd ~ dnorm(63.9, 0.0025)
EKc ~ dnorm(70.4,0.0025)
EKo ~ dnorm(36.0,0.0025)
EVcmax ~ dnorm(65.4, 0.0025)
EJmax ~ dnorm(46.1, 00025)
Egm ~ dnorm(49.6, 0.001)
### normal prior Dark Respiration from Brassica LR data (umol m-2 s-1)
Rd25 ~ dunif(0, 25)
### normal prior centered on 0  for mesophyl cond (umol m-2 s-1 Pa-1)
gm25 ~ dnorm(2.5, 0.0025)T(0.1, 142)
### flat prior on gamma star (Pa)
gammaS25 ~ dnorm(3.74,0.86)T(0.24,7.24)
### normal prior for0.01Kc (Pa)
Kc25 ~ dnorm(27.24,0.01)T(0,100)
### flat prior for Ko(Pa)
Ko25 ~ dnorm (30400, 0.00000004)T(0,65400)
### normal prior from Agricultural Species WULLSCHLEGER  (umol m-2 s-1) (1993)
Vcmax25 ~ dnorm(90, 0.000625)T(0,100000)
# curvature parameter on LR curve (thetaJ)
thetaJ ~ dunif(0.1, 0.99999)
##normal prior from Agricultural Species WULLSCHLEGER  (umol m-2 s-1) (1993)
Jmax25 ~ dnorm(171, 0.000308)T(0,100000)
## quantum yield phiJ from Chl Fl data
phiJ ~ dunif(0,.5)
### prior on precision
tau ~ dgamma(.001 , .001)
}
" 

