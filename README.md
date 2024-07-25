# S.GA.NOBO.IPM
Data, R and NIMBLE code, and table of parameter estimates for an integrated population model estimating population dynamics and demographic rates for a population of northern bobwhites (Colinus virginianus) in southern Georgia, USA, 1998-2022.
---
Authors: William B. Lewis, Chloé R. Nater, Justin A. Rectenwald, D. Clay Sisson, and James A. Martin
<br />
Manuscript Title: Use of integrated population models for assessing density-dependence and juvenile survival in Northern Bobwhites (Colinus virginanus)

---

# Metadata

# Code_NOBO_IPM_S_GA.R

Code for running the integrated population model in R and NIMBLE is contained in the 'Code_NOBO_IPM_S_GA' R file. The state population model is structured by sex (males and females), age (juveniles, subadults, and adults), and season (breeding season: April - September, non-breeding season: October - March). The biological year starts on April 1st. Subadults surviving to the start of the breeding season transition to adults, and all adults (regardless of sex) survive to each month of the breeding season through a single breeding survival parameter. Adults surviving to the start of each month June-September produce a number of chicks in that month based on the month-and-sex-specific per-capita productivity rate. The sex ratio of chicks is assumed to be 50/50. Juveniles from each of the four monthly cohorts (June-September) survive to the start of the non-breeding season based on a cohort-specific juvenile daily survival, an average of 30 days/month, and the number of months between the cohort month and October. Juveniles transition to subadults at the start of the non-breeding season in October. During the non-breeding season, male and female adults and subadults survive each month based on the sex-and-age-specific survival rates. Adult breeding and sex-and-age-specific non-breeding survival is directly estimated from radiotelemetry data, monthly sex-specific productivity is directly estimated from nest monitoring data, and November population size is directly estimated from covey count surveys. Juvenile survival for each monthly cohort (June-September) is inferred using a combination of monthly breeding productivity, November abundance surveys, ratios of adults to subadults and ratios of subadults from the monthly cohorts in the November/December traping and harvest data, and informative priors. Negative density-dependent effects are directly incorporated on estimates of breeding survival, non-breeding survival, and productivity. The model is run using parallel processing with 3 cores and can take several days. <br />
 

<br />
<br />

# S.GA.NOBO.IPM.data.gzip

The data for the northern bobwhite (Colinus virginianus, NOBO) project are stored in the 'S.GA.NOBO.IPM.data' gzip file. These data were collected from 1998 - 2022 at a managed property in southern Geogria, USA. The data
comes from 7 main databases: radiotelemetry survival data, monthly productivity data from the breeding season (June-September), November covey count data, trapping (capture) ratio data from November, and harvest ratio data from November and December.
## Radiotelemetry survival data
Males and females were captured annually using baited traps during the late fall and spring. A subset of captured birds were fit with a radiotransmitter and tracked at least twice weekly until either mortality or radito failure. Capture histories were condensed into bi-weekly periods.
## Monthly productivity data
Nests were found during the breeding season primarily through radiotelemetry, though some were found incidentally. The majority of nests were incubated by females, though a percentage of nests in each year were incubated by males. Nests were monitored until either hatch or fail, and fate, clutch size, number hatched, and identity of attending adult were recorded. Productivity data was aggregated to the number of chicks produced from successful nests attended by males and females in each month of the breeding season. In some cases, the clutch size or number hatched was not recorded for a nest. For these nests, we performed a preliminary analysis to estimate the missing data based on the clutch sizes and/or hatch rates observed in other nests for that month/year. This dataset also incorporated the number of breeding adults of each sex/month/year. This was calculated as the number of radiotracked adults of each sex alive at the start of the month, plus the number of birds for which nests were found incidentally.
## November covey count data
Coveys were surveyed annually in mid-October - mid-November using a 4-person quadrat-sampling method. About 3/4 of the 12 sampling grids were surveyed in a given year. Bird dogs were used at a subset of grids each year to flush coveys and estimate covey size. Covey counts are related to population abundance through observed average covey sizes, the proportion of grids surveyed, covey calling availability, and detection probability. Calling availability (probability of a covey calling during the survey) is estimated based on the observed neighbor density and parameters reported by Wellendorf et al. (2004). Detection (probability of an available covey being detected by at least one observer) is estimated based on data from Howell et al. (2021). 
## November trapping data
Birds were captured annually using baited traps during mid-October - mid-November. Captured birds were banded with a uniquely-number aluminum leg band and classified by sex and age (adults or subadults) based on plumage. Subadults were aged to the specifically monthly cohort in which they were hatched based on the pattern of primary feather molt (Rosene 1969).
## November and December harvest data
Data from bobwhite harvested on the property was subset to November (mid-November - end of November) and December (end of November - mid-December). Harvested subadults were aged to monthly cohorts as with the trapping data.<br />
<br />
## The 43 objects in the gzip data file are described below
### years
The years of data collection
### nyears
The number of years of the study
### n.b.months
The number of months during the breeding season (April - September)
### n.breeding.months
The number of months during which chicks enter the population (June - September). A few nests hatched at the end of May, these were lumped in with the June cohort.
### CH.sex
The sex of birds tracked via radiotelemetry (1=male, 2=female).
### CH.age
The age/season of birds tracked via radiotelemetry for each bi-weekly period. This data is constrained to the deployment period, from deployment until either mortality or censure. 0 represents mortality, 1 represents alive non-breeding subadult, 2 represents alive non-breeding adult, and 3 represents alive breeding adult.
### CH.state
The alive/dead state of birds tracked via radiotelemetry for each bi-weekly period (0=dead, 1=alive). This data is constrained to the deployment period, from the deployment until either mortality or censure.
### CH.year
The year of the study for each bi-weekly tracking period.
### nCH
The number of birds tracked via radiotelemetry.
### trackdates
The bi-weekly tracking periods, beyond the deployment period, for which each bird was tracked via radiotelemetry.
### trackperiods
The number of bi-weekly tracking periods, beyond the depolyment period, for which each bird was tracked via radiotelemetry.
### chicks.nest.mean and chicks.nest.sd
4x2x25 arrays giving the mean and standard deviation, respectively, of the number of chicks produced in each month of the breeding season (June - September, x-dim), for each sex (1=males, 2=females, y-dim), and in each year of the study (z-dim). Some nests were missing data on the clutch size or number of chicks hatched. We estimated the chick production for these nests in a preliminary analysis using the observed clutch sizes and hatch rates for observed nests in the same month/year, then summed the estimated and observed number of chicks produced from nests hatching in each month, year, and from each sex. These datsets contain the mean and standard deviation of the posterior samples from the preliminary analysis. A value of 0 for chicks.nest.sd represents that no nests were missing information on number hatched for that month, sex, and year.
### N.tracked
A 4x2x25 array giving the number of adult birds from each month of the breeding season (June - September, x-dim), for each sex (1=males, 2=females, y-dim), and in each year of the study (z-dim) from which productivity data was collected.
### count
The number of coveys detected on annual November count surveys.
### csize.mean and csize.sd
The mean and standard deviation, respectively, of the number of birds/covey/year estimated from coveys flushed by bird dogs.
### csize.min and csize.max
The 2.5% and 97.5% quantiles, respectively, of observed covey sizes across the entire study period.
### effort
The proportion of the 12 grid cells surveyed for bobwhite covey abundance in each year.
### avail.mean and avail.sd
The annual mean and standard deviation, respectively, of covey calling availability. The estimates are on the logit scale and were calculated based on observed neighbor density on covey counts (number of coveys detected - 1) and parameters relating neighbor density and calling availability reported by Wellendorf et al. (2004).
### avail.min and avail.max
The minimum and maximum, respectively, of covey calling availability estimates across the entire study period. The estimates are on the logit scale.
### pdet_mean and pdet_sd
The mean and standard deviation, respectively, of the probability that an available coveys will be detected by at least one observer during surveys. The estimates are on the logit scale and were calculated from data provided in Howell et al. (2021).
### trap.am.af.sub
A 25x3 matrix giving the number of adult males (1st column), adult females (2nd column), and subadults of either sex (3rd column) captured and banded during the fall trapping season in each year of the study (rows). Subadults are not grouped by sex because many could not be reliably differentiated.
### trap.sub.cohort
A 25x4 matrix giving the number of subadults from the June (1st column), July (2nd column), August (3rd column), and October (4th column) breeding cohorts captured and banded during the fall trapping season in each year of the study (rows). Subadults captured in the fall were aged and backdated to hatch month based on primary feather molt as in Rosene (1969).
### trap.N
The number of birds captured each year during the fall trapping period.
### trap.cohort.N
The number of subadult birds captured each year during the fall trapping period.
### nyears.trap
The number of years for which trapping data is available.
### harv.cohort.Nov
A 24x4 matrix giving the number of subadults from the June (1st column), July (2nd column), August (3rd column), and October (4th column) breeding cohorts harvested from mid-late November in each year of the study (rows). Subadults captured in the fall were aged and backdated to hatch month based on primary feather molt as in Rosene (1969).
### harv.cohort.Dec
A 23x3 matrix giving the number of subadults from the June/July (1st column), August (2nd column), and October (3rd column) breeding cohorts harvested from late-November to mid-December in each year of the study (rows). Subadults from the June and July cohorts could not be reliably separated based on primary molt during this period, so they were lumped for analysis. Subadults captured in the fall were aged and backdated to hatch month based on primary feather molt as in Rosene (1969).
### harv.ratio.years
The years of the study for which November harvest data is available.
### harv.ratio.years.n
The number of years for which November harvest data is available.
### harv.ratio.years.Dec
The years of the study for which December harvest data is available.
### harv.ratio.years.Dec.n
The number of years for which December harvest data is available.
### Area
The total area (ha) of the survey grids.
### phi.j.prior.mean
The mean daily survival rate of juveniles from the June (1st), July (2nd), August (3rd), and September (4th) monthly breeding cohorts. Estimates are on the logit scale and were taken from Terhune et al. (2019).
### phi.j.prior.sd
The standard deviation for daily survival rates of juveniles. The estimate is on the logit scale and was taken from Terhune et al. (2019). <br />

<br />
<br />

# Param_Estimates_NOBO_IPM_S_GA.pdf

Yearly estimates of population parameters from the integrated population model is stored in the 'Param_Estimates_NOBO_IPM_S_GA' pdf. The mean, standard deviation, median, and 95% credible intervals of posterior samples are provided for parameters. Included parameters are: population abundance, density, and growth rate for April and November; daily and overall (hatch month - start of October) juvenile survival rates for each monthly cohort (June - September); global mean monthly, realized monthly, and overall (April - September) breeding survival rates; global mean monthly, realized monthly, and overall (October - March) sex-and-age-specific non-breeding survival rates; global mean and realized sex-specific per-capita productivity rates for each month (June - September); total realized sex-specific per-capita productivity rates for the entire breeding season; ratios of adult males, adult females, and subadult birds n November; ratios of subadults from each monthly cohort (June - September) in November; average bobwhite covey size in November; covey calling availability during November surveys; covey detection probability (probability of at least one observer detecting a covey, conditional on availability) during November surveys; effect of April population density on breeding survival; effect of October population density on non-breeding survival; and effect of monthly population density (June - September) on monthly per-capita productivity rates.
