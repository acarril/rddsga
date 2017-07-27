{smcl}
{* *! version 0.8.1 Jul 2017}{...}
{title:Title}

{pstd}
{hi:rddsga} {hline 2} Subgroup analysis with propensity score weighting in RDD settings


{title:Syntax}

{p 8 16 2}
{cmd:rddsga} {depvar} {it:assignvar} [{indepvars}] {ifin}
{cmd:,} {it:options}
{p_end}

{phang}
{it:assignvar} is the assignment variable for which there is a known cutoff at which the conditional mean of the treatment variable changes abruptly.{p_end}

{synoptset 22 tabbed}{...}
{synopthdr}
{synoptline}
{syntab :RD design}
{p2coldent:* {opth sg:roup(varname)}}subgroup indicator variable{p_end}
{synopt :{opth t:reatment(varname)}}indicator for the assignment variable above the cutoff; if not specified, a sharp RDD is assumed{p_end}
{synopt :{opt c:utoff(real)}}specifies the cutoff value in {it:assignvar}; default is 0{p_end}
{p2coldent:* {opt bw:idth(real)}}specifies the bandwidth around the cutoff{p_end}

{syntab :Balance}
{p2coldent:+ {opth bal:ance(varlist)}}variables for which the propensity score weighting is calculated; default is {indepvars}{p_end}
{synopt :{opt probit}}predict propensity score using a {manhelp probit R:probit} model; default is {manhelp logit R:logit}{p_end}
{synopt :{opt nocom:sup}}do not restrict sample to area of common support{p_end}

{syntab :Model}
{synopt :{opt first:stage}}estimate the first stage regression model{p_end}
{synopt :{opt reduced:form}}estimate the reduced form regression model{p_end}
{synopt :{opt iv:reg}}estimate the instrumental variables regression model{p_end}
{synopt :{opth vce(vcetype)}}{it:vcetype} may be {opt un:adjusted},
   {opt r:obust}, {opt cl:uster} {it:clustvar}, {opt boot:strap},
   {opt jack:knife}, or {opt hac} {help ivregress##kernel:{it:kernel}}{p_end}
{synopt :{opt quad:ratic}}use quadratic spline; default is linear{p_end}

{syntab :Reporting/Output}
{synopt :{opt dibal:ance}}display original balance and propensity score weighting balance tables and statistics{p_end}
{synopt :{opth psw:eight(newvar)}}name of new variable with propensity score weighting; if not specified, no variable will be generated{p_end}
{synopt :{opth com:sup(newvar)}}name of new binary variable indicating common support; if not specified, no variable will be generated{p_end}
{synopt :{opth psc:ore(newvar)}}name of new variable with propensity score; if not specified, no variable will be generated{p_end}
{synoptline}
{p2colreset}{...}
{p 4 6 2}* These options must be specified.{p_end}
{p 4 6 2}+ This option must be specified if {it:indepvars} is empty.{p_end}
{p 4 6 2}{it:indepvars} may contain factor variables; see {help fvvarlist}.{p_end}


{marker description}{...}
{title:Description}

{pstd}
{cmd:rddsga} allows to conduct a binary subgroup analysis in RDD settings based on propensity score weighting (PSW).
Observations in each subgroup are weighted by the inverse of their conditional probabilities to belong to that subgroup, given a set of covariates.
Performing RDD analysis separately within each weighted subgroup eliminates potential confounding differences due to other observable factors that may vary systematically across (uneweighted) subgroups.

{pstd}
The program computes the PSW vector for the covariates in {it:indepvars}, which are also used as control variables in the model.
A separate set of variables for which the PSW balance is computed may be optionally specified with {opt balance(varlist)}.
The computed vector may be stored as a new variable using {opt psweight(newvar)}.
In order to assess the statistical significance of the difference in means for each covariate, {cmd:rddsga} uses a t-test on the equality of means and reports the resulting p-value, as well as the (weighted) standardized mean difference.
The resulting balance tables are stored as matrices (see {help rddsga##results:stored results} below), and may also be displayed using {opt dibalance}.

{pstd}
Options {opt firststage}, {opt reducedform} and {opt ivreg} may be used to estimate the first stage, reduced form and instrumental variables regression models.
Standard variance estimator options may be passed onto the models with {opt vce(vcetype)}.
The estimated coefficients for the interaction of indicator variables for each subgroup and a treatment indicator are reported, along with (robust) standard errors.
Although {cmd:rddsga} can output results without any additional packages, the presence of {browse "http://repec.sowi.unibe.ch/stata/estout/":estout} is automatically detected and used to produce better-looking output.
In any case, the full estimation results can be retrieved with {help estimates dir}.

{pstd}
Additional details regarding the methodology implemented by {cmd: rddsga} can be found in the project's {browse "https://gitlab.com/acarril/rddsga/wikis/home":repository}.
An application can be found in {help rddsga##mainpaper:Gerardino, Litschig, Olken and Pomeranz (2017)}.


{marker options}{...}
{title:Options}

{dlgtab:RD design}

{phang}
{opt sgroup(varname)} specifies a subgroup indicator variable.
This variable must be a {it:dummy} (values 0 or 1).
This option must be specified.

{phang}
{opt treatment(varname)} specifies an indicator variable for the assignment variable ({it:assignvar}) above the cutoff. If not specified, a sharp RDD is assumed.

{phang}
{opt cutoff(real)} specifies the cutoff value in {it:assignvar}; default is 0 (assuming normalized {it:assignvar}).

{phang}
{opt bwidth(real)} specifies a symmetrical the bandwidth around the cutoff.
This option must be specified.

{dlgtab:Balance}

{phang}
{opt balance(varlist)} specifies variables for which the propensity score weighting is calculated.
If not specified, variables in {it:indepvars} are used.
This option is useful if one wants to balance for a different set of covariates that the ones used as controls in the model.
This option must be specified if {it:indepvars} is empty.

{phang}
{opt probit} indicates that the propensity score is computed after a fitting a  {manhelp probit R:probit} model; default is {manhelp logit R:logit}.

{phang}
{opt nocomsup} indicates that the sample is not to be restricted to the area of common support.

{marker options_model}{...}
{dlgtab:Model}

{phang}
{opt firststage} estimates the first stage regression model.

{phang}
{opt reducedform} estimates the reduced form regression model.

{phang}
{opt ivreg} estimate the instrumental variables regression model.

{phang}
{opt vce(vcetype)} specifies the variance estimators options. See {help  vce_option}.

{phang}
{opt quadratic} indicates that a quadratic spline is to be used for full interaction with subgroup indicators. If not specified, a linear spline is used.

{dlgtab:Reporting}

{phang}
{opt dibalance} display original balance and propensity score weighting balance tables and statistics.
This balance is computed for each covariate in {it:indepvars}, unless {opt balance(varlist)} is specified.

{phang}
{opt psweight(newvar)} specifies a name for a new variable with the propensity score weighting. If not specified, no variable will be generated.

{phang}
{opt comsup(newvar)} specifies a name for a new binary variable indicating common support. If not specified, no variable will be generated.

{phang}
{opt pscore(newvar)} specifies a name for a new variable with the propensity score. If not specified, no variable will be generated.


{marker examples}{...}
{title:Examples}

{pstd}Setup (click {browse "https://github.com/acarril/rddsga#installation":here} for details on getting ancillary files){p_end}
{phang2}{cmd:. use rddsga_synth}{p_end}

{pstd}Assess covariate imbalance using one covariate{p_end}
{phang2}{cmd:. rddsga Y runvar, balance(X1) sgroup(G) bwidth(10) dibal}{p_end}

{pstd}Store computed PSW vector balancing for all covariates{p_end}
{phang2}{cmd:. rddsga Y runvar, balance(X1 X2) sgroup(G) bwidth(10) psweight()}{p_end}

{pstd}Fit reduced form model{p_end}
{phang2}{cmd:. rddsga Y runvar, balance(X1 X2) sgroup(G) bwidth(10) reduced}{p_end}

{pstd}Assess balance and fit reduced form model using minimal-MSE bandwidth ({help rddsga##imbens2012:Imbens and Kalyanaraman, 2012}){p_end}
{phang2}{cmd:. rddsga Y runvar, balance(X1 X2) sgroup(G) bwidth(6.095) reduced dibal}{p_end}


{marker results}{...}
{title:Stored results}

{pstd}
{cmd:rddsga} stores the following in {cmd:r()}:

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Scalars}{p_end}
{synopt:{cmd:r(oribal_N_G0)}}number of observations in subgroup 0 (original balance){p_end}
{synopt:{cmd:r(oribal_N_G1)}}number of observations in subgroup 1 (original balance){p_end}
{synopt:{cmd:r(oribal_Fstat)}}F-statistic (original balance){p_end}
{synopt:{cmd:r(oribal_pvalue)}}F-statistic p-value (original balance){p_end}
{synopt:{cmd:r(oribal_avgdiff)}}Average of absolute values of standardized differences (original balance){p_end}

{synopt:{cmd:r(pswbal_N_G0)}}number of observations in subgroup 0 (PSW balance){p_end}
{synopt:{cmd:r(pswbal_N_G1)}}number of observations in subgroup 1 (PSW balance){p_end}
{synopt:{cmd:r(pswbal_Fstat)}}F-statistic (PSW balance){p_end}
{synopt:{cmd:r(pswbal_pvalue)}}F-statistic p-value (PSW balance){p_end}
{synopt:{cmd:r(pswbal_avgdiff)}}Average of absolute values of standardized differences (PSW balance){p_end}

{p2col 5 15 19 2: Matrices}{p_end}
{synopt:{cmd:r(pswbal)}}balance table matrix (original balance){p_end}
{synopt:{cmd:r(oribal)}}balance table matrix (PSW balance){p_end}
{p2colreset}{...}

{pstd}
Additionally, {cmd:rddsga} stores all estimation results for the specified models ({opt firststage}, {opt firststage} and/or {opt ivregress}; see {help rddsga##options_model:Model options} above).
Both the unweighted and PSW models are stored using {help estimates store}.
The list of stored models can be retrieved using {help estimates dir}.
The full list of estimates that may be stored is described below.

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Estimates}{p_end}
{synopt:{cmd:unw_first}}Unweighted first stage{p_end}
{synopt:{cmd:psw_first}}PSW first stage{p_end}
{synopt:{cmd:unw_reduced}}Unweighted reduced form{p_end}
{synopt:{cmd:psw_reduced}}PSW reduced form{p_end}
{synopt:{cmd:unw_ivreg}}Unweighted instrumental variables{p_end}
{synopt:{cmd:psw_ivreg}}PSW instrumental variables{p_end}


{marker authors}{...}
{title:Authors}

{pstd}
Alvaro Carril (maintainer){break}
J-PAL LAC{break}
acarril@fen.uchile.cl

{pstd}
Andre Cazor{break}
J-PAL LAC{break}
ajcazor@uc.cl

{pstd}
Maria Paula Gerardino{break}
Inter-American Development Bank{break}
mariage@iadb.org

{pstd}
Stephan Litschig{break}
National Graduate Institute for Policy Studies{break}
s-litschig@grips.ac.jp

{pstd}
Dina Pomeranz{break}
University of Zurich and NBER{break}
dina.pomeranz@econ.uzh.ch


{marker disclaimer}{...}
{title:Disclaimer}

{pstd}
This software is provided "as is", without warranty of any kind.
If you have suggestions or want to report problems, please create a new issue in the {browse "https://github.com/acarril/rddsga/issues":project repository} or contact the project maintainer.
All remaining errors are our own.


{marker references}{...}
{title:References}

{marker imbens2012}{...}
{phang}Imbens, Guido, and Karthik Kalyanaraman. "Optimal Bandwidth Choice for the Regression Discontinuity Estimator." Review of Economic Studies 79, no. 3 (2012): 933–59.

{marker mainpaper}{...}
{phang}Gerardino, Maria Paula, Stephan Litschig, Benjamin Olken, and Dina Pomeranz. 2017.
"Can Audits Backfire? Evidence from Public Procurement in Chile".
{it:Working Paper}.
