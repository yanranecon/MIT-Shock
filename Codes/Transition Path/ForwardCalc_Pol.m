function [kPolM, IndM] = ForwardCalc_Pol(r, w, vNext, cS, Ts)
%% Documentation
%{
1. Function Description:
This function solves for the sequences of policy functions given the next
period value function, vNext, using Bellman equation.


2. Inputs:
(1). r:     interest rate, a scalar
(2). w:     wage, a scalar
(3). vNext: next period value function, a nk*ns matrix


3. Output:
(1). kPolM: policy function, a nk*ns matrix
(2). IndM:  indicator, showing the optimal decision for k' is the #th
            element in k grid vector
%}


%%  Compute the utility part in Bellman equation
% ns matrices, whose size is nk*nk
% Each matrix is given productivity s, the consumption for each k (row), and the all the possible k' (column)
cPolM     = ones(cS.nk, cS.nk, cS.ns); 

% nk*nk matrix of k', given each k (row), the possible choices of k' (colum)
kChoiceM  = ones(cS.nk,1)*Ts.kGridV';

% The BC is c+k'=ws+(1+r)k, hence c=ws+[(1+r)k-k']
% Then compute the matrix of (1+r)k-k' for each given k and all possible k'
% e.g. the first row of first half means when s=1
% kM(1,1): when k=kGridV(1) and k'=kGridV(1), the value of (1+r)k-k'
% kM(1,2): when k=kGridV(1) and k'=kGridV(2), the value of (1+r)k-k'
% And I need to get T kM matrices, each one is the kM matrix in each
% transition period. Hence kM is a nk*nk*T matrix
kM        = (1+r)*Ts.kGridV*ones(1,cS.nk) - kChoiceM;

% Since c=ws+[(1+r)k-k']
% cPolM contains ns matrices. Each matrix corresponds to one realization of 
% s with kM
for i = 1 : cS.ns
    cPolM(:,:,i) = w*cS.s(i,1) + kM;
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

% uM has ns matrices. Each one is a nk*nk matrix corresponds to one 
% realization of productivity. Combine all these matrices together, and
% make uM to be one big (ns*nk)*nk matrix
uM             = reshape(permute(uM,[1,3,2]),[],size(uM,2));


%% Solve for Policy Function Using Bellman Equation
vNext          = vNext';
walue          = uM + cS.beta*kron(cS.P*vNext,ones(cS.nk,1));
[~, IndM]      = max(walue,[],2);


%% Outcomes
IndM           = reshape(IndM,cS.nk,cS.ns);

kPolM = zeros(cS.nk, cS.ns);
for i = 1 : cS.ns
kPolM(:,i)     = cS.kGridV(IndM(:,i));
end


end