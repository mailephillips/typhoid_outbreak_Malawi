# typhoid_outbreak_Malawi


## Overview
This repository has the associated code for simulating typhoid fever outbreaks for the manuscript "Cost-effectiveness analysis of typhoid conjugate vaccines in an outbreak setting: a modeling study" in XXXXX. 

The functions simulate a typhoid fever outbreak under specified conditions. The function carries out one iteration of a stochastic transmission dynamic model for typhoid fever (using one setting as an example).


## How To Use

### Running one iteration

* **Load specified data** (`data` folder)
    * `data_timestep_pt1_randT.Rdata` : the associated data for the randomized outbreak timing simulations (Scenario 1 in manuscript)
    * `data_timestep_pt1.Rdata` : the associated data for the fixed outbreak timing and non-outbreak simulations 
    * Example: `load('data_timestep_pt1_randT.Rdata')`
   
* **Load the specified stochastic simulation function** (`simulate_outbreak` folder)
    * `stochtyph_outbreak.R`: the function with non-randomized timing, used for the non-vaccination (base case) strategy 
    * `stochtyph_outbreak_rand.R`: the function with randomized timing, used for the non-vaccination (base case), reactive vaccination, or preventative strategies in Scenario 1
    * `stochtyph_outbreak_vaxb4.R`: the function with non-randomized timing, used for the preventative strategies 
    * `stochtyph_outbreak_vaxb4_rand.R`: the function with randomized timing, used for the preventative strategies
    * Example: `source('stochtyph_outbreak_rand.R')`
  
* **Run function with specified parameters**
    * `tmax` : how many timesteps (in weeks) to run (we used 4173)
    * `dt` : *delta t*, the fraction of a timestep to use for the stochastic runs (we used 0.1)
    * `vax` : logical -- do you want to introduction vaccination? (TRUE or FALSE)
    * `catchup` : logical -- if vaccinating, do you want to include a one-time catchup campaign up to 15 years of age? (TRUE or FALSE)
    * `threshIR` : if using incidence rate as a threshold to identify an outbreak, specify the numeric incidence rate above which would be considered an outbreak. If not using incidence rate, use a large number (i.e., 1000)
    * `threshSD` : if using standard deviations above the monthly mean as a threshold to identify an outbreak, specify the number of standard deviations above which would be considered an outbreak (in main analysis, we used 15)
    * `t_to_vax` : time delay in implementing vaccination (after identifying the outbreak). Must be a whole number, and must be in the same time measurement specified (for example, if you want to have a 1-month delay, this parameter would be `round(4.345*1)/dt`)
    * Example for vaccination with catchup campaign, 15 SD, 1-month delay: `stochtyph_identify(tmax=4173, dt=.1, vax=T, catchup=T, threshIR=1000, threshSD=15, t_to_vax=round(4.345*1)/dt )`


## Author

Maile Thayer Phillips, PhD, MS

---

> GitHub [@mailephillips](https://github.com/mailephillips) &nbsp;&middot;&nbsp;

