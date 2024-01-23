# S.GA.NOBO.IPM
Data, R and NIMBLE code, and table of parameter estimates for Integrate Population Model estimating population dynamics and demographic rates for a population of northern bobwhites (Colinus virginianus) in southern Georgia, USA, 1998-2022.
---
authors: William B. Lewis, Justin A. Rectenwald, D. Clay Sisson, and James A. Martin


---

# Metadata

# S.GA.NOBO.IPM.data.gzip

The data for the northern bobwhite (Colinus virginianus, NOBO) project are stored in the 'S.GA.NOBO.IPM.data' gzip file. These data were collected from 1998 - 2022 at a managed property in southern Geogria, USA. The biological year for the population model was specified as starting on April 1st, with the breeding season extending April - September and the non-breeding season extending from October - March. In the population model, juveniles enter the population in June-September, transition to subadults at the start of the non-breeding season in October, and transition to adults at the start of the subsequent breeding season in April. The data
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
A 4x2x25 array giving the mean and standard deviation, respectively, of the number of chicks produced in each month of the breeding season (June - September, x-dim), for each sex (1=males, 2=females, y-dim), and in each year of the study (z-dim). Some nests were missing data on the clutch size or number of chicks hatched. We estimated the chick production for these nests in a preliminary analysis using the observed clutch sizes and hatch rates for observed nests in the same month/year, then summed the estimated and observed number of chicks produced from nests hatching in each month, year, and from each sex. These datsets contain the mean and standard deviation of the posterior samples from the preliminary analysis. A value of 0 for chicks.nest.sd represents that no nests were missing information on number hatched for that month, sex, and year.
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
The area (ha) of the survey grids.
### phi.j.prior.mean
The mean daily survival rate of juveniles from the June (1st), July (2nd), August (3rd), and September (4th) monthly breeding cohorts. Estimates are on the logit scale and were taken from Terhune et al. (2019).
### phi.j.prior.sd
The standard deviation for daily survival rates of juveniles. The estimate is on the logit scale and was taken from Terhune et al. (2019). <br />

<br />


# Code_NOBO_IPM_S_GA.R

Code for running the integrated population model in R and NIMBLE is contained in the 'Code_NOBO_IPM_SE_GA' R file. 

