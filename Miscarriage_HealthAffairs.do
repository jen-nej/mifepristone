* Created by Sungsik Hwang 
* last updated May 8th 2024 
* shwang97@wisc.edu
* last checked by Jenna Nobles (May 8th 2024)
* jnobles@ssc.wisc.edu

clear
set more off
macro define DTA "Miscarriage"
graph set window fontface helvetica

* abortion law by states
import excel using $DTA/cross_re, first
save $DTA/cross_re, replace

* Mortality rates from Li et al. (2002) 
* converting mortality rates to probabilities
import excel using $DTA/Lifetable_Li, first clear
gen qx = mmid/(1+(1/2)*mmid) 
gen qx_lb = mlow/(1+(1/2)*mlow) 
gen qx_ub = mhigh/(1+(1/2)*mhigh) 

rename week gest 
keep gest qx qx_lb qx_ub
save $DTA/Li_clean, replace

* appending NCHS for 21-45 weeks 

use $DTA/NCHS_fet_all, replace 
rename oegest_comb gest
keep if gest >=21
gen qx = prob 
gen qx_lb =.
gen qx_ub =.
forvalue i = 1/`=_N' {
	local n = enter[`i']
	local d = numFD[`i']
	cii means `n' `d', poisson level(95)
	replace qx_ub = r(ub) in `i'
	replace qx_lb = r(lb) in `i'
}

* generate live birth probabilities
gen qx_birth= numBD/enter
keep gest qx qx_lb qx_ub qx_birth
save $DTA/NCHS_Li, replace

use $DTA/Li_clean, replace
append using $DTA/NCHS_Li
save $DTA/MDLifetable, replace

* calculating miscarriage rates
use $DTA/MDLifetable, replace
* assigned the live birth probability 0 before week 21
replace qx_birth = 0 if qx_birth==. 

* constructing life table
gen lx =100000 if _n==1
gen dx= lx*qx if _n==1
gen dx_birth = lx*qx_birth if _n==1

forval x = 2/46 {
replace lx = lx[_n-1]-dx[_n-1]-dx_birth[_n-1] if _n==`x'
replace dx = lx*qx if _n==`x'
replace dx_birth = lx*qx_birth if _n==`x'
}

* rep. upper and lower bounds
foreach v in lb ub {
gen lx_`v' = 100000 if _n==1 
gen dx_`v' = lx_`v'*qx_`v' if _n==1 
gen dx_birth_`v' = lx_`v'*qx_birth if _n==1
}

foreach v in lb ub {
forval x = 2/46 {
replace lx_`v' = lx_`v'[_n-1]-dx_`v'[_n-1]-dx_birth_`v'[_n-1] if _n==`x'
replace dx_`v' = lx_`v'*qx_`v' if _n==`x'
replace dx_birth_`v' = lx_`v'*qx_birth if _n==`x'
}
}

* calculating the probability of transitioning to live birth
egen prob_surv = total(dx_birth)
scalar prob_surv = prob_surv

* radix cancels out in later calculation
foreach v in lb ub {
egen prob_surv_`v'= total(dx_birth_`v')
scalar prob_surv_`v' = prob_surv_`v'
}

* first trimester miscarriages
egen first = total(dx) if gest<13
scalar first = first
foreach v in lb ub {
egen first_`v' = total(dx_`v') if gest<13
scalar first_`v' = first_`v'
}

* second trimester miscarriages
egen second = total(dx) if gest<26
scalar second = second
foreach v in lb ub {
egen second_`v' = total(dx_`v') if gest<26
scalar second_`v' = second_`v'
}


* calculating scalar for the number of miscarriages
* first trimester 
scalar m_rate = first / prob_surv
foreach v in lb ub {
scalar m_rate_`v' = first_`v' / prob_surv_`v'
}

* first + second trimester
scalar m_rate1 = second / prob_surv
foreach v in lb ub {
scalar m_rate1_`v' = second_`v' / prob_surv_`v'
}

* calculating number of total live births for all 4 years
use $DTA/Natality_2018-2021, replace
collapse (sum) birth=birth, by(state)

* average number of births per year
replace birth = birth/4
gen miscarriage = birth * m_rate
gen miscarriage_1 = birth *m_rate1

foreach x in lb ub  {
gen miscarriage_`x' = birth * m_rate_`x'
gen miscarriage_1_`x' = birth * m_rate1_`x'
}

* merging abortion law by states
merge 1:1 state using $DTA/cross_re
gen state_re  = state
labmask state_re, val(state_ab)
save $DTA/Miscarriage, replace

* Graph 
* Figure 1 
use $DTA/Miscarriage, replace
gen id = 1
collapse (sum) miscarriage miscarriage_1 miscarriage_1_lb miscarriage_1_ub, by(id)

gen mis_1 = (miscarriage/1000)*0.102
gen mis_2 = (miscarriage/1000)
gen mis_3 = (miscarriage_1/1000)

gen c_1 = mis_3 - mis_2 
gen c_2 = mis_2 - mis_1 
gen c_3 = mis_1 

keep c_*
gen idt = 1 
reshape long c_, i(idt)
drop idt 
lab define _j 1 "Second trimester miscarriages" 2 "Would not benefit from access" 3 "Could benefit from access"
lab values _j _j

export excel using $DTA/exhibit_1, replace 


* Figure 2 
* wihtout abortion ban 
use $DTA/Miscarriage, replace
gsort -miscarriage
gen n = _n

lab define n 1 " "
forval x =2/50 { 
lab define n `x' " ", add	
}
lab values n n 

labmask n if abortion>=1 , val(state_ab)

replace miscarriage =. if abortion==0
replace miscarriage_1 =. if abortion==0
foreach x in lb ub {
replace miscarriage_1_`x'=. if abortion==0	
}

gen mis_1 = (miscarriage/1000) *0.102
gen mis_2 = (miscarriage/1000)
gen mis_3 = (miscarriage_1/1000)

foreach x in lb ub {
replace miscarriage_1_`x' = miscarriage_1_`x'/1000 
}

keep if miscarriage~=. 
keep state_re mis_1 mis_2 mis_3 miscarriage_1_lb miscarriage_1_ub

gen First_trimester_benefit = mis_1
gen First_trimester_nobenefit = mis_2 - mis_1 
gen Second_trimester = mis_3 - mis_2

gen LowerCIvalue = miscarriage_1_lb
gen UpperCIvalue = miscarriage_1_ub

rename state_re State 
keep State First_trimester_benefit First_trimester_nobenefit Second_trimester LowerCIvalue UpperCIvalue
export excel using $DTA/exhibit_2_2, replace first(var)

* abortion ban
use $DTA/Miscarriage, replace
gsort -miscarriage
gen n = _n

lab define n 1 " "
forval x =2/50 { 
lab define n `x' " ", add	
}
lab values n n 
labmask n if abortion==0 , val(state_ab)

replace miscarriage =. if abortion==1
replace miscarriage_1 =. if abortion==1
foreach x in lb ub {
replace miscarriage_1_`x'=. if abortion==1	
}

gen mis_1 = (miscarriage/1000) *0.102
gen mis_2 = (miscarriage/1000)
gen mis_3 = (miscarriage_1/1000)

foreach x in lb ub {
replace miscarriage_1_`x' = miscarriage_1_`x'/1000 
}

keep if miscarriage~=. 
keep state_re mis_1 mis_2 mis_3 miscarriage_1_lb miscarriage_1_ub

gen First_trimester_benefit = mis_1
gen First_trimester_nobenefit = mis_2 - mis_1 
gen Second_trimester = mis_3 - mis_2

gen LowerCIvalue = miscarriage_1_lb
gen UpperCIvalue = miscarriage_1_ub

rename state_re State 
keep State First_trimester_benefit First_trimester_nobenefit Second_trimester LowerCIvalue UpperCIvalue
export excel using $DTA/exhibit_2_3, replace first(var)

* All states
foreach x in 2 3 {
import excel using $DTA/exhibit_2_`x' , first clear 
save $DTA/exhibit_2_`x', replace 
}

use $DTA/exhibit_2_2, replace
append using  $DTA/exhibit_2_3, gen(Status)
gen status = Status
lab define status 0 "Current abortion ban" 1 "All other states"
lab values status status
drop Status
order status, after(State)
export excel using $DTA/exhibit_2_1, replace first(var)




