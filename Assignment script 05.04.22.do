// Further Stats - Assignment

//Set directory and load in data
clear all
cd "U:\1. MODULES\Further Stats\Assignment"
capture log close
log using "Assignment_SFB.smcl", replace
use "BMIData.dta"

//////////////////////////////////////////////////////////////////////////////// Q1 - Describe the variables using summary statistics and appropriate graphs

describe // dataset contains PID, the dependant variable (pbmichange), plus binary independant variables (treat, male & stress), and continuous independant variables (age, bmi0) 
summarize // look at mean and the spread of the data for each variable 
mdesc // check for missing data (none in this dataset)

//Look at distribution of continuous variables:
hist pbmichange // Relatively normally distributed.  Important as the dependant variable.
graph export pbmichange_dist.png, replace
hist age // Not normally distributed, peaks at 18 and 55.  However, OK as conditioning variable.
graph export age_dist.png, replace
hist bmi0 // Not normally distributed, peaks around 26 and 28. However, OK as conditioning variable.
graph export bmi0_dist.png, replace

// compare variables in each group to see if they are similar in baseline characteristics (including t-test for continuous variables and proportions test for binary variables)
by treat, sort : sum 

by treat, sort: sum if male==1 // Similar proportion of males and females in the treatment group, but more females than males in the control group. 
prtest male, by(treat) // More males in the treatment group (p=0.0000)

by treat, sort: sum if stress==1 // More people with chronic stress in the control group (33.3% versus 20.1%).
prtest stress, by(treat) // More people with chronic stress in the control group (p=0.000)

by treat, sort: sum age
ttest age, by(treat) // People in control group significantly older on average (difference of 8.27, p=0.0000).  Similar SDs.  Same age range (18-55 years)

by treat, sort: sum bmi0
ttest bmi0, by(treat) // People in control group significantly heavier on average (difference of 0.19 points, p=0.0000).  Similar SDs.

// Look at relationship between the independant and dependant variables:
twoway (scatter pbmichange bmi0)
graph export BMI0_scatter.png, replace
pwcorr bmi0 pbmichange, sig // correlation coefficeint and statistical significance

twoway (scatter pbmichange age)
graph export age_scatter.png, replace
pwcorr age pbmichange, sig // correlation coefficeint and statistical significance

graph box pbmichange, over(male, relabel(1 "Female" 2 "Male"))
graph export sex_boxplot.png, replace
twoway (kdensity pbmichange if male==0) (kdensity pbmichange if male==1), xtitle("% BMI change") ytitle("Density") legend(order(1 "Female" 2 "Male"))
graph export sex_kernel.png, replace

graph box pbmichange, over(stress, relabel(1 "No chronic stress" 2 "Chronic stress"))
graph export stress_boxplot.png, replace
twoway (kdensity pbmichange if stress==0) (kdensity pbmichange if stress==1), xtitle("% BMI change") ytitle("Density") legend(order(1 "No chronic stress" 2 "Chronic stress"))
graph export stress_kernel.png, replace

graph box pbmichange, over(treat, relabel(1 "Control" 2 "Treatment"))
graph export treat_boxplot.png, replace

//////////////////////////////////////////////////////////////////////////////// Q2. NAIVE ANALYSIS
// Calculate the (unadjusted) mean pBMIchange for the groups and calculate the significance in differences

// unadjusted mean pBMIchange by treatment group
mean pbmichange, over(treat) // Greater average % change with treatment (7.42% compared to 1.6%). Also, greater variability in outcome in the control group, with some people putting on weight (range -7.84-10.51%)

// test difference in mean cost
ttest pbmichange, by(treat) // interpretation: the difference in unadjusted mean % BMI change between the groups is 5.827 and this is statistically significant (p=0.0000 at 5% significance).

// keep mean pBMI change by treatment group
global a = r(mu_1) // stores the control group mean in the global macro a
global b = r(mu_2) // stores the treatment group mean in the global macro b
global diff= $a -$b // stores the difference in mean costs in the global macro diff
display "Difference in unadjusted means: $diff" // displays what is stored in the global macro diff

//by gender
mean pbmichange, over(male)
ttest pbmichange, by(male) // significant

// by age group
sum age, d // median age is 39
replace agegroup=1 if age <39
replace agegroup=2 if age >=39

mean pbmichange, over(agegroup)
ttest pbmichange, by(agegroup) 

//////////////////////////////////////////////////////////////////////////////// INVERSE PROBABILITY WEIGHTING (IPW)
/* This method estimates the propoensity score for treatment (i.e. the likelihood of being treated) based on the covariates, and uses this to reweight the predicted outcome.  When there are few observations in one group but a lot in the other, the observations that are sparse are given a higher weight.  This allows you to improve (and also check) the comparability of groups. Can use probit or logit*/

// When thinking about which variables to include looking at statistical significance, direction of change and whether needing to make any potential changes to the variable types
// Age - Continuous.  Appears to be a negative correlation between age & pbmichange (linear) from the scatter.  Could test for polynomial.
// Male - Dummy variable.  No changes to be made.  
// BMI - Continuous.  Appears no relationship from the scatter but could test for polynomial.
// Stress - Dummy variable.  No changes can be made.
// Treat - Dummy variable. No changes can be made.

/////////////////////////////////////////////////////////////////////////////// Step 1 - Ensuring correct specification of probit model to be used in the IPW:

//Comparing models using linear and quadratic continuous variables:

// IPW model 1 - all coefficients
probit treat age male bmi0 stress  // Male and baseline bmi not-sgnificant.  Other variables in right direction.
estimates store model1
// Pseduo R2 0.130, Prob > chi2   = 0.0000
estat gof // Simple goodness of fit test (Paerson's chi-squared) - p-value 0.8509 - cannot reject the model in terms of goodness of fit

// IPW model 2 - using both age squared and BMI squared (most complex model)
probit treat c.age c.age#c.age male c.bmi c.bmi0#c.bmi0 stress // Most variables non-significant
estimates store model2
// Pseudo R2: 0.131, Prob > chi2   = 0.0000
estat gof // p-value 0.7663 - cannot reject the model in terms of goodness of fit
lrtest model2 model1 // p>0.05, simpler model (1) may be adequate

// IPW model 3 - using age squared only 
probit treat c.age c.age#c.age male bmi0 stress // Age-squared is significant, but age and baseline BMI now not. 
estimates store model3
// Pseudo R2: 0.131, Prob > chi2   = 0.0000
estat gof //simple goodness of fit test (Paerson's chi-squared): p-value 0.7874 - cannot reject the model in terms of goodness of fit
lrtest model3 model1 // p<0.05, simpler model (1) may not be adequate i.e. may need age-squared.
lrtest model3 model2 // p>0.05, simpler model (3) may be adequate 

// IPW model 4 - using BMI squared only
probit treat age male c.bmi c.bmi0#c.bmi0 stress // BMI and its squared co-efficients still not significant
estimates store model4
// Pseudo R2: 0.130, Prob > chi2   = 0.0000
estat gof // simple goodness of fit test (Paerson's chi-squared): p-value 0.8367 - cannot reject the model in terms of goodness of fit
lrtest model4 model1 // p>0.05, simpler model (1) may be adequate i.e. don't appear to need bmi.squared

//////// Conclusion age.squared may improve fit but bmi.squared does not.  Taking model 3 as best model so far:

// Trialling removing variables that do not appear to be confounders (male and baseline BMI).

//IPW model 5 - Remove male only
probit treat c.age c.age#c.age bmi0 stress // all except age significant
estimates store model5
estat gof // p-value 0.8114 - cannot reject the model in terms of goodness of fit
fitstat // McFadden's R2: 0.131
lrtest model3 model5 // p>0.05 - cannot reject the hypothesis that the simpler model (model 5) is adequate

// IPW model 6 - Remove baseline BMI only
probit treat c.age c.age#c.age male stress // all except age significant
estimates store model6
// Pseudo R-squared: 0.130, Prob > chi2   = 0.0000
estat gof // simple goodness of fit test: p-value 0.0316 - reject the model in terms of goodness of fit

// IPW model 7 - Remove both male and BMI
probit treat c.age c.age#c.age stress // all significant (except age)
estimates store model7
// Pseudo R2: 0.128,  Prob > chi2   = 0.0000
estat gof //simple goodness of fit test: p-value 0.068 - can't reject the model in terms of goodness of fit
lrtest model5 model7 // p<0.05 - can reject the hypothesis that the simpler model (model 7) is adequate

estimates table _all, star stats(aic bic) // Both AIC and BIC favour model 5 

// Based on likelihood ratio tests, Perason's chi-squared GOF test, AIC and BIC scores, model 5 appears best.
// Mc-Fadden's r-squared slightly favours a more complex model, however the difference in r-squared is negligable (0.130 vs 0.131) 
// SELECTED MODEL - MODEL 5: probit treat c.age c.age#c.age bmi0 stress (i.e. removed sex and included a quadratic for age)

//////////////////////////////////////////////////////////////////////////////// Step 2: Estimating ATE:
teffects ipw (pbmichange) (treat c.age c.age#c.age bmi0 stress, probit), pomeans
teffects ipw (pbmichange) (treat c.age c.age#c.age bmi0 stress, probit), aeq 
/* ATE estimate based on the difference in the weighted means = 5.06361  
- Co-efficients provide direction of the effect and significance but not the size of the effect:
- age and age.squared have a negative effect, which is significant for age.squared only
- baseline bmi has a significant negative effect
- stress has a significant negative effect*/

//////////////////////////////////////////////////////////////////////////////// Step 3: Looking at balance & overlap:

// Balance:
tebalance overid // H0: Covariates are balanced.  P>0.05 so we can't reject that the covariates are balanced (overall).  
tebalance summarize c.age c.age#c.age bmi0 stress, baseline // looking at baseline (mean and variance)
tebalance summarize c.age c.age#c.age bmi0 stress // looking at after re-weighting (mean and variance)
// Standardized differences: Close to zero achieved after weighting.  Gap closed less for stress than for other variables
// Variance: Close to 1 achieved after weighting.  Gap closed less for stress than for other variables.
tabstat ps, by(treat) stats(N mean median min max)

// explore balance in more detail for each covariate
tebalance density stress 
graph export balance_stress.png, replace
tebalance density age
graph export balance_age.png, replace
tebalance density bmi0
graph export balance_bmi0.png, replace

// Look at weightings assigned to the data points to achieve the more balanced groups:
predict double pscore, ps tlevel(1)
gen double wei = 1/pscore if treat==1
replace wei =1/(1-pscore) if treat==0
summ wei // Provides the mean, min and max weghting given to an individual data point. 

// Explore how well overlap has been achieved:
teffects overlap, ptlevel(1)
graph export balance_overall.png, replace

// Plot of graph that shows the weightings given:
 // Stress:
twoway (scatter pbmichange stress if treat==0 [w=wei], sort mcolor(blue) msymbol(circle_hollow) legend(label(1 "Control")) xtitle("Chronic stress") ytitle("pbmichange"))  ///
  (scatter pbmichange stress if treat==1 [w=wei], mcolor(red) msymbol(circle_hollow) legend(label(2 "Treated")) )
  graph export overlap_stress.png, replace
   // Age:
 twoway (scatter pbmichange age if treat==0 [w=wei], sort mcolor(blue) msymbol(circle_hollow) legend(label(1 "Control")) xtitle("Age") ytitle("pbmichange"))  ///
  (scatter pbmichange age if treat==1 [w=wei], mcolor(red) msymbol(circle_hollow) legend(label(2 "Treated")) )
    graph export overlap_age.png, replace
   // Baseline BMI: 
 twoway (scatter pbmichange bmi0 if treat==0 [w=wei], sort mcolor(blue) msymbol(circle_hollow) legend(label(1 "Control")) xtitle("Baseline BMI") ytitle("pbmichange"))  ///
  (scatter pbmichange bmi0 if treat==1 [w=wei], mcolor(red) msymbol(circle_hollow) legend(label(2 "Treated")) )
    graph export overlap_bmi0.png, replace
  
// Seen that some individuals have been given heavy weights. One data point was given the equivelant weighting of 36 observations. Also, this is shown graphically, denoted by some very large circles for individuals in the treatment group who were older, had stress, or had particularly high baseline BMI. 
// This can therefore introduce bias as we are heavily relying on a few individual observation and we don't know how usual their observations are.  However, this has resulted in good weighting of the groups.
