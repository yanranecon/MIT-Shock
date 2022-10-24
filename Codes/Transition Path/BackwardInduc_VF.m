function valueM_ts = BackwardInduc_VF(valueM_LastT,r,w,cS,Ts)
%% Documentation
%{
1. Function Description:
This function solves for the sequences of value functions given the last
period value function, valueM_LastT, using backward induction

The basic idea is
(1). Given the the prices (r and w) in each period and given the last
     period (t=T) value function, I can find the policy function and
     corresponding value function in period t=T-1 using Bellman equation
     V(k,s) = max u(c) + sum_(s->s') Prob(s'|s)V(k',s')
(2). Now I know the prices in each period and the value function in period
     t=T-1. Apply the same method, I can solve for policy function and
     corresponding value function in period t=T-2
     ...
     Keep doing this, I can solve all the value functions backwards from 
     t=T-1, ... 1, setting V^T = valueM_LastT


2. Inputs:
(1). cS.nk:      the number of possible realizations of k
(2). cS.kGridV:  vector of k grid
(3). cS.s:       vector of productivity state
(4). cS.kMin:    lower bound of k
(5). cS.kMax:    Upper bound of k


3. Output:
valueM_ts: a nk*ns*T matrix, each nk*ns is a matrix of value function in
           that transition period, the last valueM (t=T) is the value
           function in new steady state
%}


%%  Compute the utility part in Bellman equation
% ns matrices, whose size is nk*nk
% Each matrix is given productivity s, the consumption for each k (row), and the all the possible k' (column)
% I have T-1 such matrices, one for each transition period
cPolM    = ones(cS.nk*cS.ns, cS.nk, Ts.T-1); 

% nk*nk matrix of k', given each k (row), the possible choices of k' (colum)
kChoice  = ones(cS.nk,1)*Ts.kGridV';
% I need ns kChoice matrices, each one corresponds to one state of s
kChoiceM = kron(ones(cS.ns,1),kChoice);

% The BC is c+k'=ws+(1+r)k, hence c=ws+[(1+r)k-k']
% Then compute the matrix of (1+r)k-k' for each given k and all possible k'
% This is a (ns*nk)*nk matrix called kM, 
% e.g. ns=2, then the first half of matrix kM (nk*nk) corresponds to s=1
% the second half of matrix kM (nk*nk) corresponds to s=2
% e.g. the first row of first half means when s=1
% kM(1,1): when k=kGridV(1) and k'=kGridV(1), the value of (1+r)k-k'
% kM(1,2): when k=kGridV(1) and k'=kGridV(2), the value of (1+r)k-k'
% And I need to get T kM matrices, each one is the kM matrix in each
% transition period. Hence kM is a nk*nk*T matrix
kM        = zeros(size(kChoiceM,1),size(kChoiceM,2),Ts.T-1);
for t = 1:Ts.T-1
kM(:,:,t) = (1+r(t,1))*kron(ones(cS.ns,1),Ts.kGridV)*ones(1,cS.nk) - kChoiceM;
end

% Since c=ws+[(1+r)k-k']
% cPolM contains ns matrices. Each matrix corresponds to one realization of 
% s with kM
for t = 1:Ts.T-1
    cPolM(:,:,t) = w(t,1)*kron(cS.s,ones(cS.nk,1)) + kM(:,:,t);
end

% Pick up the elements in cPolM that are negative
cPolM(cPolM<0) = 0;
% This is a logic value. 
% It means that if c is a negative number, I assign logic value 'false'  (i.e. 0)

% Due to the last step, some c have numeric value (for those >= 0), some c have logic value
% The value of utility which is calculated base on logic value is 'inf'
uM             = CRRA_Utility(cPolM,cS);
% Assign a very large negtive value (i.e. -realmax) to utility when it is inf
uM(isinf(uM))  = -realmax;


%% Backward Induction
valueM_ts           = zeros(size(Ts.valueM_end,1),size(Ts.valueM_end,2),Ts.T);
valueM_ts(:,:,Ts.T) = valueM_LastT;
Ind_ts              = zeros(cS.nk,cS.ns,Ts.T-1);

% Preparation
walue  = zeros(size(uM));
valueM = zeros(cS.nk*cS.ns,1,Ts.T-1);
Ind    = zeros(cS.nk*cS.ns,1,Ts.T-1);
vNext  = zeros(cS.ns, cS.nk, Ts.T-1);

% Run backward induction routine
for t = Ts.T-1 : -1 : 1
    vNext(:,:,t)                = valueM_ts(:,:,t+1)';
    walue(:,:,t)                = uM(:,:,t) + cS.beta*kron(cS.P*vNext(:,:,t),ones(cS.nk,1));    
    [valueM(:,:,t), Ind(:,:,t)] = max(walue(:,:,t),[],2);
    valueM_ts(:,:,t)            = reshape(valueM(:,:,t),cS.nk,cS.ns);
    Ind_ts(:,:,t)               = reshape(Ind(:,:,t),cS.nk,cS.ns);  
end


end