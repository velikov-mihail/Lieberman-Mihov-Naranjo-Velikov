* Import data for tables
import delimited "Data\stata_data.csv", clear

* Create dates
gen year = floor(dates/100)
gen month = dates - year*100
drop dates
gen dates = ym(year, month)
format dates %tm
drop year month

* Reformating ROE as it is read as a string
destring roe, replace force

* Winsorize
winsor2 zpaydex lpc whitedwu distress r ivol logme gp roe pmf marketshare logbm, suffix(_w) cuts(1 99) by(dates)

* Rescale (for presentational purposes)
replace car3 = car3*100

* Create lagged variables
sort permno dates
by permno: gen zpaydex_w_l1 = zpaydex_w[_n-1] if dates==dates[_n-1]+1
by permno: gen zpaydex_w_l2 = zpaydex_w[_n-2] if dates==dates[_n-2]+2
by permno: gen zpaydex_w_l3 = zpaydex_w[_n-3] if dates==dates[_n-3]+3
by permno: gen zpaydex_w_l4 = zpaydex_w[_n-4] if dates==dates[_n-4]+4
by permno: gen zpaydex_w_l5 = zpaydex_w[_n-5] if dates==dates[_n-5]+5
by permno: gen zpaydex_w_l6 = zpaydex_w[_n-6] if dates==dates[_n-6]+6

* Label variables
label variable whitedwu_w "FC"
label variable distress_w "Distress"
label variable r_w "Past Returns"
label variable ivol_w "IVol"
label variable logme_w "log(ME)"
label variable gp_w "GP/A"
label variable zpaydex_w_l1 "Z_{paydex}"
label variable pmf_w "PMF"
label variable marketshare_w "Market Share"
label variable roe_w "ROE"

********************************************************************************
* Table 4, Panel A: Zpaydex panel regressions

* Month fixed effects
local fixed_effects "dates"

reghdfe zpaydex_w whitedwu_w, absorb(`fixed_effects') vce(cluster permno dates)
est store reg1, title("")

reghdfe zpaydex_w distress_w, absorb(`fixed_effects') vce(cluster permno dates)
est store reg2, title("")

reghdfe zpaydex_w r_w, absorb(`fixed_effects') vce(cluster permno dates)
est store reg3, title("")

reghdfe zpaydex_w ivol_w, absorb(`fixed_effects') vce(cluster permno dates)
est store reg4, title("")

reghdfe zpaydex_w logme_w, absorb(`fixed_effects') vce(cluster permno dates)
est store reg5, title("")

reghdfe zpaydex_w roe_w, absorb(`fixed_effects') vce(cluster permno dates)
est store reg6, title("")

reghdfe zpaydex_w whitedwu_w distress_w r_w ivol_w logme_w roe_w, absorb(`fixed_effects') vce(cluster permno dates)
est store reg7, title("")

* Add fixed effects in the table output
estfe . reg*, labels(dates "Year-month Fixed Effects" permno "Firm Fixed Effects")

esttab reg1 reg2 reg3 reg4 reg5 reg6 reg7 ///
	using "Results\table4_panelA.tex" , wide cells(b(fmt(2)) t(fmt(2) par([ ]) abs)) label nomtitles ///
	indicate(`r(indicate_fe)') ///
	keep(whitedwu_w distress_w r_w ivol_w logme_w roe_w) ///
	order(whitedwu_w distress_w r_w ivol_w logme_w roe_w) ///
	modelwidth(8) collabels( ) substitute(_ \_) ar2(4) booktabs alignment(D{.}{.}{-1}) legend style(tex) nostar replace addnotes ("")


********************************************************************************
* Table 4, Panel B: Zpaydex as an independent variable

* Month fixed effects
local fixed_effects "dates"

gen zpaydex_temp = zpaydex_w_l1
label variable zpaydex_temp "Z_{Paydex}"
reghdfe distress_w zpaydex_temp, absorb(`fixed_effects') vce(cluster permno dates)
est store reg1, title("")

replace zpaydex_temp = zpaydex_w_l3
reghdfe ivol_w zpaydex_temp, absorb(`fixed_effects') vce(cluster permno dates)
est store reg2, title("")

* Firm and month fixed effects
local fixed_effects "permno dates"

replace zpaydex_temp = zpaydex_w_l1
reghdfe distress_w zpaydex_temp, absorb(`fixed_effects') vce(cluster permno dates)
est store reg3, title("")

replace zpaydex_temp = zpaydex_w_l3
reghdfe ivol_w zpaydex_temp, absorb(`fixed_effects') vce(cluster permno dates)
est store reg4, title("")

* Add fixed effects in the table output
estfe . reg*, labels(dates "Year-month Fixed Effects" permno "Firm Fixed Effects")

esttab reg1 reg2 reg3 reg4 ///
	using "Results\table4_panelB.tex" , wide cells(b(fmt(2)) t(fmt(2) par([ ]) abs)) label nomtitles ///
	indicate(`r(indicate_fe)') ///
	keep(zpaydex_temp) ///
	order(zpaydex_temp) ///
	modelwidth(8) collabels( ) substitute(_ \_) ar2(4) booktabs alignment(D{.}{.}{-1}) legend style(tex) nostar replace addnotes ("")

drop zpaydex_temp


************************************************;

* Table 6: LPC panel regressions
local fixed_effects "dates"

reghdfe lpc_w whitedwu_w, absorb(`fixed_effects') vce(cluster permno dates)
est store reg1, title("")

reghdfe lpc_w distress_w, absorb(`fixed_effects') vce(cluster permno dates)
est store reg2, title("")

reghdfe lpc_w pmf_w, absorb(`fixed_effects') vce(cluster permno dates)
est store reg3, title("")

reghdfe lpc_w marketshare_w, absorb(`fixed_effects') vce(cluster permno dates)
est store reg4, title("")

reghdfe lpc_w logme_w, absorb(`fixed_effects') vce(cluster permno dates)
est store reg5, title("")

reghdfe lpc_w gp_w, absorb(`fixed_effects') vce(cluster permno dates)
est store reg6, title("")

reghdfe lpc_w whitedwu_w distress_w pmf_w marketshare_w logme_w gp_w, absorb(`fixed_effects') vce(cluster permno dates)
est store reg7, title("")

* Add fixed effects in the table output
estfe . reg*, labels(dates "Year-month Fixed Effects" permno "Firm Fixed Effects")
return list

esttab reg1 reg2 reg3 reg4 reg5 reg6 reg7 ///
	using "Results\table6.tex" , wide cells(b(fmt(2)) t(fmt(2) par([ ]) abs)) label nomtitles ///
	indicate(`r(indicate_fe)') ///
	keep(whitedwu_w distress_w pmf_w marketshare_w logme_w gp_w) ///
	order(whitedwu_w distress_w pmf_w marketshare_w logme_w gp_w) ///
	modelwidth(8) collabels( ) substitute(_ \_) ar2(4) booktabs alignment(D{.}{.}{-1}) legend style(tex) nostar replace addnotes ("")


************************************************;

* Figure 2: Earnings predictability with Z PAYDEX
local fixed_effects "dates"
matrix parameters = J(1, 4, .)

forval i = 1(1)6{
	reghdfe car3 zpaydex_w_l`i', absorb(`fixed_effects') vce(cluster permno dates)

	matrix rtable = r(table)
	matrix parameters`i' = J(1, 4, .)
	matrix parameters`i'[1,1]=`i'
	matrix parameters`i'[1,2]=rtable[1,1]
	matrix parameters`i'[1,3]=rtable[5,1]
	matrix parameters`i'[1,4]=rtable[6,1]
	
	matrix parameters = parameters \ parameters`i'
}

matrix list parameters
svmat parameters, name(var)
keep var*
keep if missing(var1)==0
rename var1 lag
rename var2 beta
rename var3 ll
rename var4 ul

*Graph
graph twoway ///
(scatter beta lag, mcolor(black) lcolor(black)) ///
(rcap ll ul lag, mcolor(blue) connect(direct) lpattern(dash) lcolor(black)), ///
graphregion(color(white)) plotregion(color(white)) bgcolor(white) ///
xtitle("Lag in months ({&tau})", size(medium) margin(medium)) ///
ytitle("{&beta}{sub:{&tau}}", size(medium) margin(medium)) ///
xlabel("1 2 3 4 5 6", nogrid labsize(medsmall)) ///
ylabel(, nogrid angle(0)) ///
legend(off) ///
yline(0, lcolor(black) lpattern(solid))
graph export "Figures\earningsPredictabilityFig.pdf", replace

