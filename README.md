# COVID-TSDPD
Seeding efficient large-scale public health interventions in diverse spatial-social networks

# Significance Statement
Implementation of public health interventions at the regional level to achieve a national goal should consider spillovers through spatial-social networks, regional heterogeneity, and the feedback effects on health outcomes and network structures. We build a **threshold spatial dynamic panel data (TSDPD)** model with these features to identify intervention effect channels, compare the clinical outcomes and cost-effectiveness of different seeding strategies, demonstrate the spatial distributions of benefits and costs, and inform the optimal termination of campaign efforts. This analytical framework can be applied to address other challenges such as climate change adaptation, antimicrobial resistance, and biodiversity loss where spillover effects are critical. Analysis based on COVID-19 suggests that leveraging social network influence leads to the most cost-effective intervention.

# Notes
******************************************************************************************************************
This version: 2022/12/2

Written by:
Han, Xiaoyi (xiaoyihan@xmu.edu.cn)

******************************************************************************************************************
In this simulation code, we produce parameter estimates for the Threshold SDPD model with formation equations for vaccination,
 as well as within and cross state flow.

To implement the Bayesian 95% confidence interval, we include a function 'hpdi' provied by James P. LeSage in the file named 'jplv7', 
which can be downloaded  from http://www.spatial-econometrics.com/

We also include the data, as well as the cross state travel flow in the year 2019 and 2021 used in the empirical study.

state0205_0415_v13:        The excel file for data used in the empirical study
Stateflow2019_weights:     The M file for cross state flow in the year 2019
StateflowV_weights:        The M file for cross state flow in the year 2021

Main_program_empirical_withinstateflow49t5b_sci:  The code to generate parameter estimates of the TSDPD model. It reads the data and spatial weights, and implements mcmcwfbd8b, which is a function to conduct MCMC estimation
for the TSDPD model


jplv7: file for functions used in the econometric toolbox for matlab, which includes the function to generate Bayesian 95% confidence interval (written by James P. LeSage)



