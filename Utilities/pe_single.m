function [pe,pd] = pe_single(TS,m,t)

%==========================================================================
% Calcualtes the permutation entropy
% Input:  TS = time series;
%         m  = order of permuation entropy (pattern length, typically between 3 and 7)
%         t  = delay time between pattern points
% Output: pe = permuation entropy
%         pd = probability distribution
%==========================================================================
permlist = perms(1:m);              % list of possible permutations of order m
c = zeros(1,length(permlist));      % initialise array of permutation counts
    
 for a = 1:length(TS)-t*(m-1)
     [~,V] = sort(TS(a:t:a+t*(m-1)));
     for b = 1:length(permlist)
         if (abs(permlist(b,:)-V')) == 0
             c(b) = c(b) + 1 ;
             break
         end
     end
 end

pd = c/sum(c);              % probability distribution
c = c(c~=0);                % remove zeroes for PE calculation
p = c/sum(c);               % probability dist with zeroes removed
sp = -sum(p .* log(p));     % Shannon entropy of p
pe = sp/log(factorial(m));  % Normalised Permutation entropy