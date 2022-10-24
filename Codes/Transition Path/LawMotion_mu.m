function muprime = LawMotion_mu(Ind, muNow,cS)
%% Documentation
%{
1. Function Description:
   Given prices and the policy function k'=g(k,s), this program find next 
period distribution using law of motion of distribution mu with policy function
and the current distribution muNow

Law of motion for mu(k,s):
          mu(k',s')=sum_k sum_s 1{k'=g(k,s)} prob(s'|s) mu(k,s)


2. Input
(1). Ind:     position indicator that shows the position of k' in k state vector
(2). muNow:   distribution in period t (current period), nk*ns matrix
(3). cS.P:    transition probability (ns*ns matrix)

3. Output:
muprime:      distribution in period t+1 (next period)
%}


%% Law of Motion for mu
muprime = zeros(cS.nk,cS.ns);

for iprime = 1:cS.nk
   for jprime = 1:cS.ns
      for i=1:cS.nk
         for j=1:cS.ns
            if Ind(i,j) == iprime
               muprime(iprime,jprime) = muprime(iprime, jprime) + muNow(i,j)*cS.P(j,jprime); 
            end
         end
      end
   end
end


end