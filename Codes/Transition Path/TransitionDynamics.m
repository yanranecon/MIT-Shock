function Seq = TransitionDynamics(cS, Ts)
%% Documentation
%{
1. Function Description:
This function solves for the sequences of policy function, value function
and distribution between the initial steady state (t=0) and the end steady 
state (t=T) by finding the value for K (agrregate capital) that clears 
capital markets in each transition period.

Notice that there are T-1 transition periods between the two steady states.

The algorithm is
(1). Have an initial guess for a sequence of K (i.e. T Ks). Each K is the 
     guess for aggregate capital in each transition period.
(2). Given the guessed K in each period, calculate the prices (r and w) in
     each period using firm FOCs
(3). Backward induction: solve the value function backwards from 
     t=T, ... 1, setting V^T = V^SS_end
(4). Forward calculation: given policy function and distribution in t=0,
     which is the policy function and distribution in initial steady state, 
     compute the distribution in t=1 using law of motion for distribution.
     Given distribution in t=1 and the value function in period 2 (sloved 
     in Step-3), solve for policy function k' and aggregate K* that clears
     the capital market
     In summary, in this step, solve for distribution and K* for each period
(5). If max|K* - K_guess| > tolerance level, update guess and go to Step-2

The idea used to find K is similar as the one used in Aiyagari_Calib
Given K, solve the HH problem and get policy function for saving
Compute the aggregate saving S --> a func of prices --> a func of K 
Find K which makes capital market to clear 


2. Input:
(1). K_initial:  aggregate K in initial steady state
                 Use it as my initial guess for K
(2). valueM_end: period t=T+1 value function, i.e. value function in new steady state
                 Use it to do backward induction to solve for value 
                 function in each transition period
(3). mu_initial: period t=0 distribution, i.e. distribution in initial steady state
                 Use it to do forward calculation


3. Output:
(1). r: a T*1 vector, sequence of interest rates during transition
(2). w: a T*1 vector, sequence of wages during transition
(3). K: a T*1 vector, sequence of aggregate capital during transition
%}


%%  Preparation
mu_ts        = zeros(cS.nk, cS.ns, Ts.T);
mu_ts(:,:,1) = Ts.mu_initial;

tol          = 10^(-5);
dif          = tol + 1;

guessK       = Ts.K_initial * ones(Ts.T-1,1);      % Guess a sequence of K


%% Main
% Run loop to find the sequence of K*
while dif > tol
    
    [r,w]              = Prices_Firm_FOC(guessK, cS);
    
    % Backward induction to get value function
    valueM_Backward    = BackwardInduc_VF(Ts.valueM_end, r, w, cS,Ts);

    % Forward calculation to get K* that clears capital market
    K_star             = zeros(Ts.T-1,1);

    for t = 1:Ts.T-1
        
        optS           = optimoptions('fsolve', 'Display', 'none', 'TolFun', 1E-3);
        wrapper        = @(K) CapitalMktClearing_Transition(K, valueM_Backward(:,:,t+1), mu_ts(:,:,t), cS, Ts);
        [K_star(t,1), ~, exitFlag] ...
                       = fsolve(wrapper, guessK(t,1), optS);

        [~, Seq]       = CapitalMktClearing_Transition(K_star(t,1),valueM_Backward(:,:,t+1), mu_ts(:,:,t), cS, Ts);
        mu_ts(:,:,t+1) = Seq.muPrime;

        % Check if fsolve converges
        if exitFlag <= 0
           warning('No convergence');
           keyboard;
        end

    end

    dif    = norm(K_star - guessK);
    guessK = 0.8*guessK + 0.2*K_star;

end


%% Save
% Save other steady state features
for t = 1 : Ts.T-1

    [~, Seq(t,1)] = CapitalMktClearing_Transition(K_star(t,1),valueM_Backward(:,:,t+1), mu_ts(:,:,t), cS, Ts);

end


end