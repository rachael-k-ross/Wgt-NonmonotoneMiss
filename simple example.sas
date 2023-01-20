
*********************************************;
*
* SAS code for simple illustrative example
* from Ross et al. 
* Accounting for Nonmonotone Missing Data using Inverse Probability Weights
*
*********************************************;

OPTIONS FORMCHAR="|----|+|---+=|-/\<>*"  symbolgen notes source source2 msglevel=I mcompilenote=all mautosource mprint mlogic;

*************************;
* Generate full data	*;
*************************;
%let n = 20000;

*Full data;
data full;
call streaminit(7);
do id = 1 to &n; 
	z = rand("bernoulli",0.5); *binary confounder;
	x = rand("bernoulli",0.2 + 0.2*z); *binary exposure;
	py = 1/(1 + exp(-1*(log(0.1/0.9) + log(2)*z))); *probability of y under no exposure;
	y0 = rand("bernoulli", py); *binary potential outcome under no exposure;
	y1 = rand("bernoulli", py + 0.1); *binary potential outcome under exposure;
	y = x*y1 + (1-x)*y0; *observed binary outcome;
	output;
end;
run;

proc means data=full;run;

***************************;
* Generate missingness	  *;
***************************;

*Set parameters of missingness models;
%let g20 = -.945; *selected to produce 25% pattern 2;
%let g21 = log(.5);
%let g22 = log(1.7);
%let g23 = 0;
%let g30 = -1.815; *selected to produce 15% pattern 3;
%let g31 = 0;
%let g32 = log(1.3);
%let g33 = 0;
%let g40 = -2.41; *selected to produce 10% pattern 4;
%let g41 = 0;
%let g42 = log(2);
%let g43 = log(.8);

*Data with missingness;
data withmiss;
set full;
call streaminit(23);
pR2 = 1/(1+exp(-1*(&g20 + &g21*z + &g22*x + &g23*y))); *probability of being in pattern 2;
pR3 = 1/(1+exp(-1*(&g30 + &g31*z + &g32*x + &g33*y))); *probability of being in pattern 3;
pR4 = 1/(1+exp(-1*(&g40 + &g41*z + &g42*x + &g43*y))); *probability of being in pattern 4;
pR1 = 1 - pR2 - pR3 - pR4; *probability of being in pattern 1;

R = rand("table",pR1,pR2,pR3,pR4); *multinomial draw using the probabilities for each pattern;
array Rs R1 R2 R3 R4;
do i=1 to 4;
	if R = i then Rs{i} = 1; else Rs{i} = 0; *create indicators for each pattern;
end;

if R in (3,4) then z = .; *set missingness for each pattern;
if R in (2,3) then y = .;
keep id z x y R R1-R4;
run;

proc freq data=withmiss; tables r; run;


*************************************;
* Implement analysis of full data   *;
*************************************;

*Fit propensity score model;
ods exclude all;
proc genmod data=full;
model x(ref='0') = z/ dist=binomial link=logit;
output out=full_p p=pscore;
run;
ods exclude none;

*Create IPTW and combined weights;
data full_pw;
set full_p;
ipw = x/pscore + (1-x)/(1-pscore);
run;

*Fit weighted outcome model;
proc genmod data=full_pw;
class id;
model y(ref='0') = x/ dist=binomial link=identity;
weight ipw;
repeated sub=id/type=ind;
estimate 'RD' intercept 0 x 1;
ods select estimates;
ods output estimates=estrd_full;
run;


*************************************;
* Implement Complete case analysis  *;
*************************************;

*Fit propensity score model;
ods exclude all;
proc genmod data=withmiss;
model x(ref='0') = z/ dist=binomial link=logit;
output out=cc p=pscore;
where R=1;
run;
ods exclude none;

*Create IPTW and combined weights;
data cc_w;
set cc;
ipw = x/pscore + (1-x)/(1-pscore);
run;

*Fit weighted outcome model;
proc genmod data=cc_w;
class id;
model y(ref='0') = x/ dist=binomial link=identity;
weight ipw;
repeated sub=id/type=ind;
estimate 'RD' intercept 0 x 1;
ods select estimates;
ods output estimates=estrd_cc;
run;
title;


***************************;
* Implement ST IPW UMLE   *;
***************************;

*Fit the joint likelihood of missingness models;
ods exclude all;
proc nlmixed data=withmiss maxiter=1000; 
parms g20 g30 g40 -2
	  g21 g22 g32 g42 g43 0;

	p2 = 1/(1+exp(-1*(g20 + g21*z + g22*x        ))); *logistic models for each pattern R>1;
	p3 = 1/(1+exp(-1*(g30         + g32*x        )));
	p4 = 1/(1+exp(-1*(g40         + g42*x + g43*y)));
	sump = p2 + p3 + p4;
	
	if R=1 then loglik=log(1-sump); 
		else if R=2 then loglik=log(p2);
		else if R=3 then loglik=log(p3);
		else if R=4 then loglik=log(p4);

	model R~general(loglik); 

	ods output parameterestimates=umlegams(keep=parameter estimate);
run;
ods exclude none;

*Transpose estimates for merging;
proc transpose data=umlegams out=umlegamsw(keep=g:) ; 
id parameter; 
run; 

*Obtain marginal Pr(R=1);
proc means data=withmiss noprint;
var R1;
output out=mnum(keep=mnum) mean=mnum;
run;

*Merge into complete cases, obtain probability of complete case, create IPMW;
data cc_umle;
set withmiss; 
where R=1;
	if _N_=1 then set umlegamsw; 
	if _N_=1 then set mnum;
	p2 = 1/(1+exp(-1*(g20 + g21*z + g22*x        )));
	p3 = 1/(1+exp(-1*(g30         + g32*x        )));
	p4 = 1/(1+exp(-1*(g40         + g42*x + g43*y)));
	p1 = 1 - p2 - p3 - p4;
	ipmw=mnum/p1; 
run;

proc means data=cc_umle min mean max sum;
	var ipmw;
	run;

*Fit weighted propensity score model;
ods exclude all;
proc genmod data=cc_umle;
model x(ref='0') = z/ dist=binomial link=logit;
weight ipmw;
output out=cc_umlew p=pscore;
run;
ods exclude none;

*Create IPTW and combined weights;
data cc_umlewts;
set cc_umlew;
iptw = x/pscore + (1-x)/(1-pscore);
ipw = ipmw*iptw;
run;

*Fit weighted outcome model;
proc genmod data=cc_umlewts;
class id;
model y(ref='0') = x/ dist=binomial link=identity;
weight ipw;
repeated sub=id/type=ind;
estimate 'RD' intercept 0 x 1;
ods select estimates;
ods output estimates=estrd_umle;
run;

***************************;
* Implement ST IPW CBE    *;
***************************;

*Standardize data;
proc standard data=withmiss mean=0 std=1 out=scaled;
var x y z;
run;

*MH algorithm in proc IML;
proc iml; 
	q=8; *# of parameters;
	c = 1e-8; *constraint recommended by ST;
	
	*Define required functions;
	start log0(x); *log that returns -infinity for 0s;
		if x <= 0 then log0 = - 1e6; *constant("big");
			else log0 = log(x);
		return(log0);
	finish log0;

	start plogis(x); *inverse logit function;
		p = 1/(1 + exp(-1*x));
		return(p);
	finish plogis; 

	start loglike(g) global(c, data, r, n); *likelihood function;
		p2 = plogis(g[1] + g[2]*data[,1] + g[3]*data[,2]                );
		p3 = plogis(g[4]                 + g[5]*data[,2]                );
		p4 = plogis(g[6]                 + g[7]*data[,2] + g[8]*data[,3]);
		sump = p2 + p3 + p4;
		if sump < 1 - c then cons = 1; *1 if constraint met, otherwise 0;
			else cons = 0;
		logl = j(n, 1, 0); * n by 1 vector to hold the individual log-likelihood contributions;
		do i = 1 to n;
			if r[i] = 1 then logl[i] = log0( (1 - sump[i]) * cons ); *if constraint not met, this will be -infinity;
				else if r[i] = 2 then logl[i] = log(p2[i]);
				else if r[i] = 3 then logl[i] = log(p3[i]);
				else if r[i] = 4 then logl[i] = log(p4[i]);
		end;
		loglike = sum(logl); *sum of individual log-likelihoods;
		return(loglike);
	finish loglike;

	start priorf(parm); *function for prior distributions;
		prior = sum(log(pdf("normal", parm, 0, sqrt(100)))); *diffuse prior recommended by ST;
		return(prior);
	finish priorf;

	*Read in data;
	use scaled ;
		read all var {r} into r;
		read all var {z x y} into data;
	close scaled;
	n = nrow(r);

	*Metropolis;
	call randseed(13);
	m = 10000; *number of iterations in chain;
	b = 5000; *burn-in;
	sigma = 0.015; *sd of proposal distribution;

	*to store acceptance indicator and parameter ests (and create starting values);
	accept = j(m, 1, 0); accept[1,] = 1; *set first acceptence to 1;
	chain = j(m, q, 0); chain[1,] = {-2 0 0 -2 0 -2 0 0}; *set first iter in chain to starting values;

	do k = 2 to m;
		oldp = chain[k-1,]; *get current parameter values (from previous iteration);
		old = loglike(oldp);  *log-likelihood value at current parameter values;
	 
		prop = randfun(q, "normal", 0, sigma)`; *draw from proposal distribution for each parameter;
		newp = oldp + prop; *proposed parameter values;
		new = loglike(newp); *log-likelihood at proposed parameter values

		*compare current and proposed log-likelihoods, with priors;
		num = new + priorf(newp);
		den = old + priorf(oldp); 
			p = exp(num - den);
			acc =  p > rand("uniform"); *flip coin to accept/reject proposed values;
			accept[k] = acc; *store whether accepted/rejected;
			if acc = 1 then chain[k,] = newp; *if accepted, store proposed values;
			if acc = 0 then chain[k,] = oldp; *if rejected, store current values;
	end;

	avgacc = accept[+] / m; print avgacc; *acceptence rate;

	post = chain[(b+1):m,]; *remove burn in; 
	ests = median(post); print ests; *save medians as gamma estimates;
	varNames = ("g1":"g8"); *create variable names; 

	create cbegamsw from ests[colname=varNames]; *save estimates as a dataset;
		append from ests; 
	close cbegamsw;

quit;

*Obtain marginal Pr(R=1);
proc means data=scaled noprint;
var R1;
output out=mnum(keep=mnum) mean=mnum;
run;

*Merge into complete cases, obtain probability of complete case, create IPMW;
data cc_scaled;
set scaled; 
where R=1;
	if _N_=1 then set cbegamsw; 
	if _N_=1 then set mnum;
	p2 = 1/(1+exp(-1*(g1 + g2*z + g3*x        )));
	p3 = 1/(1+exp(-1*(g4        + g5*x        )));
	p4 = 1/(1+exp(-1*(g6        + g7*x + g8*y)));
	p1 = 1 - p2 - p3 - p4;
	ipmw=mnum/p1; 
run;

*Merge IPMW into unscaled data;
proc sql;
create table cc_cbe as
select a.*, b.ipmw from withmiss as a inner join cc_scaled as b on a.id=b.id;
quit;

proc means data=cc_cbe min mean max sum;
	var ipmw;
	run;

*Fit weighted propensity score model;
ods exclude all;
proc genmod data=cc_cbe;
model x(ref='0') = z/ dist=binomial link=logit;
weight ipmw;
output out=cc_cbew p=pscore;
run;
ods exclude none;

*Create IPTW and combined weights;
data cc_cbewts;
set cc_cbew;
iptw = x/pscore + (1-x)/(1-pscore);
ipw = ipmw*iptw;
run;

*Fit weighted outcome model;
proc genmod data=cc_cbewts;
class id;
model y(ref='0') = x/ dist=binomial link=identity;
weight ipw;
repeated sub=id/type=ind;
estimate 'RD' intercept 0 x 1;
ods select estimates;
ods output estimates=estrd_cbe;
run;


******************************;
* Comparison of results      *;
******************************;

title "truth in full data";
proc sql;
select mean(y1) as y1, mean(y0) as y0, calculated y1 - calculated y0 as rd
from full; quit;

title "analysis of full data";
proc print data=estrd_full noobs;
run;

title "complete case analysis";
proc print data=estrd_cc noobs;
run;

title "IPW UMLE";
proc print data=estrd_umle noobs;
run;

title "IPW CBE";
proc print data=estrd_cbe noobs;
run;


