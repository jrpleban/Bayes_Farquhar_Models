# Bayes_Farquhar_Models

8 Photosynthesis Models based on 3 all combinations of 3 assumtions:

CaCc_Jm - Assumes infinite mesophyll conductance,
           electron transport based on empirical formula using a maximum rate (Jmax), quantum yield and a curvature parameters 
           there is no temperature dependency on parameters.
           
CaCc_Jf - Assumes infinite mesophyll conductance,
           electron transport based on chlorophyll fluoresces
           there is no temperature dependency on parameters.
           
CaCc_Jm_Temp - Assumes infinite mesophyll conductance,
           electron transport based on empirical formula using a maximum rate (Jmax), quantum yield and a curvature parameters 
           there is an Arrhenius style temperature dependency on parameters.
           
CaCc_Jf_Temp - Assumes infinite mesophyll conductance,
           electron transport based on empirical formula using a maximum rate (Jmax), quantum yield and a curvature parameters 
           there is an Arrhenius style temperature dependency on parameters.
           
CiCc_Jm - Assumes a mesophyll conductance limitaion,
           electron transport based on empirical formula using a maximum rate (Jmax), quantum yield and a curvature parameters 
           there is no temperature dependency on parameters.
           
CiCc_Jf - Assumes a mesophyll conductance limitaion,
           electron transport based on chlorophyll fluoresces
           there is no temperature dependency on parameters.
           
CiCc_Jm_Temp - Assumes a mesophyll conductance limitaion,
           electron transport based on empirical formula using a maximum rate (Jmax), quantum yield and a curvature parameters 
           there is an Arrhenius style temperature dependency on parameters.
           
CiCc_Jf_Temp - Assumes a mesophyll conductance limitaion,
           electron transport based on empirical formula using a maximum rate (Jmax), quantum yield and a curvature parameters 
           there is an Arrhenius style temperature dependency on parameters.
           

_Model provides a model text for implementation in rjags

_Script provides model implementatation on A/Ci data and some diagnostics.


Key References

Original Model:

Farquhar, G.D., Caemmerer, S. V. and Berry, J. A. (1980) A biochemical model of photosynthetic CO2 assimilation in leaves of C3 species. Planta, 149, 78-90.

On incorporation of mesophyll conductance:

Ethier G.J. & Livingston N.J. (2004) On the need to incorporate sensitivity to CO2 transfer conductance into the Farquhar–von Caemmerer–Berry leaf photosynthesis model. Plant, Cell & Environment, 27, 137-153.

Similar Bayesian Implimentation:

Patrick L.D., Ogle K. & Tissue D.T. (2009) A hierarchical Bayesian approach for estimation of photosynthetic parameters of C3 plants. Plant, Cell & Environment, 32, 1695-1709.


It you use these models please cite my forthcoming publication, anticipated submission to PCE

Identification of genotype specific photosynthetic traits in Brassica rapa using a multi-model Bayesian evaluation.
Pleban, J. R.1; Mackay, D. S.1; Aston, T.2; Ewers, B. E.2; Wienig, C.2
1-Geography, University at Buffalo, Buffalo, NY, USA 
2-Botany, University of Wyoming, Laramie, WY, USA







