{smcl}
{* *! version 1.0 Jul 2017}{...}
{title:Title}

{pstd}
{hi:rddsga} {hline 2} Subgroup analysis with propensity score weighting in RDD settings


{title:Syntax}

{p 8 16 2}
{cmd:rddsga} {depvar} {it:assignvar} {indepvars} {ifin}
{cmd:,} {it:options}
{p_end}

{phang}
{it:assignvar} is the assignment variable for which there is a known cutoff at which the conditional mean of the treatment variable changes abruptly.{p_end}

{synoptset 22 tabbed}{...}
{synopthdr}
{synoptline}
{p2col 3 4 4 2:* RD design}{p_end}
{synopt :{opth sg:roup(varname)}}subgroup indicator variable{p_end}
{synopt :{opth t:reatment(varname)}}indicator for the assignment variable above the cutoff; if not specified, a sharp RDD is assumed{p_end}
{synopt :{opt c:utoff(real)}}specifies the cutoff value in {it:assignvar}{p_end}
{synopt :{opt bw:idth(real)}}specifies the bandwidth around the cutoff{p_end}

{syntab :Balance}
{synopt :{opth bal:ance(varlist)}}variables for which the propensity score weighting is calculated; default is {indepvars}{p_end}
{synopt :{opt probit}}predict propensity score using a {manhelp probit R:probit} model; default is {manhelp logit R:logit}{p_end}
{synopt :{opt nocom:sup}}do not restrict sample to area of common support{p_end}

{syntab :Model}
{synopt :{opt first:stage}}estimate the first stage regression model{p_end}
{synopt :{opt reduced:form}}estimate the reduced form regression model{p_end}
{synopt :{opt iv:reg}}estimate the instrumental variables regression model{p_end}
{synopt :{opth vce(vcetype)}}{it:vcetype} may be {opt un:adjusted},
   {opt r:obust}, {opt cl:uster} {it:clustvar}, {opt boot:strap},
   {opt jack:knife}, or {opt hac} {help ivregress##kernel:{it:kernel}}{p_end}

{syntab :Reporting/Output}
{synopt :{opt dibal:ance}}display original balance and propensity score weighting balance tables and statistics{p_end}
{synopt :{opth psw:eight(newvar)}}name of new variable with propensity score weighting; if not specified, no variable will be generated{p_end}
{synopt :{opth com:sup(newvar)}}name of new binary variable indicating common support; if not specified, no variable will be generated{p_end}
{synopt :{opth psc:ore(newvar)}}name of new variable with propensity score; if not specified, no variable will be generated{p_end}
{synoptline}
{p2colreset}{...}
{p 4 6 2}* These options must be specified.{p_end}
{p 4 6 2}
{it:indepvars} may contain factor variables; see {help fvvarlist}.{p_end}


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


{marker options}{...}
{title:Options}

{dlgtab:RD design}

{dlgtab:Balance}

{dlgtab:Model}

{dlgtab:Reporting}


{marker examples}{...}
{title:Examples}


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
Additionally, {cmd:rddsga} stores all results stored in {cmd:e()} by {manhelp ivregress R:ivregress} after fitting the weighted model.


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
If you have suggestions or want to report problems, please create a new issue in the {browse "https://github.com/acarril/rddsga/issues":project repository}.
All remaining errors are our own.


{marker references}{...}
{title:References}

{marker mainpaper}{...}
{phang}Gerardino, Maria Paula, Stephan Litschig, Benjamin Olken, and Dina Pomeranz. 2017.
"Can Audits Backfire? Evidence from Public Procurement in Chile".
{it:Working Paper}.

{marker imbens_rubin_2015}{...}
{phang}Imbens, Guido W. and Donald B. Rubin. 2015.
{it: Causal Inference in Statistics, Social, and Biomedical Sciences}.
New York: Cambridge University Press.

{marker imbens_2015}{...}
{phang}Imbens, Guido W. 2015.
"Matching Methods in Practice: Three Examples".
{it:Journal of Human Resources} 50(2): 373-419.
{p_end}
