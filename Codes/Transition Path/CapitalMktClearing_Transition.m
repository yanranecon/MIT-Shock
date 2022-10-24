function [devV, Seq] = CapitalMktClearing_Transition(K, vNext, muNow, cS, Ts)
%% Documentation:
%{
This function does forward calculation.
It solves for a sequence of distributions and a sequence of K*s 
that clears capital market in each period, given policy function and 
distribution in t=0.
Policy function and distribution in t=0 is the policy function and 
distribution in initial steady state. Hence I can compute the distribution 
in t=1 using law of motion for distribution.
Given distribution in t=1 and the value function in period 2 (sloved 
in BackwardInduc_VF), solve for policy function k' and aggregate K* that 
clears the capital market.

% INPUTS:
% (1). K:      calibration targets, a scalar
% (2). vNext:  next period value function, nk*ns matrix
% (3). muNow:  joint distribution of state variables. In this model, it is
               the joint distribution between asset k and productivity s
               Given today's lambda and today's policy function, I can
               compute tomorrow's lambda. Hence the lamda I'm computing in
               this function is the distribution in next period.
               And the lambda I'm using as an input, is the current
               distribution, which was computed in last period already.


% OUTPUTS:
% (1). devV:  deviation between aggregate KSupply computed by policy function and
              distribution and calibration targets K
% (2). Seq:   other variables of interest

%}

%% Main
% Compute the prices faced by HHs
[Seq.r, Seq.w]        = Prices_Firm_FOC(K, cS);

% Solve household problem (with exogenous prices)
[Seq.kPolM, Seq.IndM] = ForwardCalc_Pol(Seq.r, Seq.w, vNext, cS, Ts);

% Next period distribution
% Compute using law of motion with policy function
Seq.muPrime           = LawMotion_mu(Seq.IndM,muNow,cS);

% Aggregate capital
Seq.wealth            = sum(Ts.kGridV'*muNow);


%%   Deviations  
devV                  = Seq.wealth - K;


end