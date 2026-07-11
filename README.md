# COVID-TSDPD
Networked Regional Resilience: Diffusion of Risk and Adaptation across Interconnected Regions

# Abstact
Regions face not only direct exposure to economic and physical shocks but also indirect consequences of actions in interconnected regions. This paper examines how socio-spatial networks shape regional resilience to contagious shocks through dynamic feedback between risk diffusion and adaptation diffusion. We build a **Threshold Spatial Dynamic Panel Data (TSDPD)** framework that incorporates both temporal and spatial transmission, along with observed and latent factors. In our framework, peer effects generate positive spillovers in risk‐mitigating adaptation across social networks, while evolving cross‐regional physical connections mediate the spatial spread of risks. Calibrated to U.S. regional data, we simulate and compare network‐aware intervention strategies that prioritize regions based on social centrality, spatial centrality, and population scale. The results show that the cost‐effectiveness of interventions depends critically on the population scale of socially or spatially central regions. Strategies that maximize aggregate risk reduction do not always align with those that are the most cost-effective, revealing fundamental trade‐offs that call for hybrid, model‐informed intervention design. Our findings underscore the importance of interregional linkages and network-aware policy design in shaping regional resilience to systemic shocks, with broader implications for economic contagion, climate adaptation, and public health.

# Notes
******************************************************************************************************************
This version: 2026/7/10

Written by: Han, Xiaoyi (xiaoyihan@xmu.edu.cn),  Xiamen University
******************************************************************************************************************
In this simulation code for empirical study, we produce parameter estimates for the Threshold SDPD model with formation equations for vaccination, and cross state flow, which correspond to Columns (II) and (III) in Table 1 of the main text.

To implement the Bayesian 95% confidence interval, we include a function 'hpdin' , Largely based on R code by Kruschke, J. K. (2015). Doing Bayesian Data Analysis,  Second Edition: A Tutorial with R, JAGS, and Stan. Academic Press / Elsevier.  See http://www.indiana.edu/~kruschke/BEST/ for R code. The function is written by Nils Winter (nils.winter1@gmail.com) and  Johann-Wolfgang-Goethe University, Frankfurt


We also include the data, as well as the cross state travel flow in the year 2019 and 2021 used in the empirical study.

state0108_0415_v2: The excel file for data used in the empirical study 
Stateflow2019_weights: The M file for cross state flow in the year 2019 starting from Feb
StateflowV_weights: The M file for cross state flow in the year 2021 starting from Feb
fb_staten: Sci social network in the vaccination equation
migration_raw: migration network in the vaccination equation

main_program_empirical_TSDPD_social_network: The code to generate parameter estimates of the TSDPD model with friendshipnetwork as social network in the vaccination equation. It reads the data and spatial weights, and implements mcmctsdpdbaselinei, which is a function to conduct MCMC estimation for the TSDPD model with formation equations

main_program_empirical_TSDPD_migration_network: The code to generate parameter estimates of the TSDPD model with friendshipnetwork as social network in the vaccination equation. It reads the data and spatial weights, and implements mcmctsdpdbaselinei, which is a function to conduct MCMC estimation for the TSDPD model with formation equations

jplv7: files for functions used in the econometric toolbox for matlab, which includes functions to normalize the spatal weights matrices in the model (written by James P. LeSage)
