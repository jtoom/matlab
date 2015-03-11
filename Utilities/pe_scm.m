function [PE,SCM,tau] = pe_scm(y,m,td)

%==========================================================================
% DESCRIPTION
% Computes the permutation entropy and statistical complexity measure as a
% function of tau (delay between ordinal pattern points)

% INPUTS
% y = timeseries to be analysed
% m = length of ordinal pattern to use
% td = maximum time delay between points used in ordinal pattern. The
%      function computes from 1 to tau.

%OUTPUTS
% PE = Permutation entropy
% SCM = Statistical complexity measure
% tau = delay used in calculation of PE and SCM
%==========================================================================

ly = length(y);
permlist = perms(1:m);      % all possible permutations of length m
h = waitbar(0,'Calculating...');

% Initialise output arrays
PE = zeros(td,1);
JSD = zeros(td,1);
Qnorm = zeros(td,1);
DISEQ = zeros(td,1);
SCM = zeros(td,1);

for t = 1:td
    c = zeros(1,length(permlist));
     for j = 1:ly-t*(m-1)
         [~,iv] = sort(y(j:t:j+t*(m-1)));
         
         for jj = 1:length(permlist)
             if (abs(permlist(jj,:)-iv')) == 0
                 c(jj) = c(jj) + 1 ;
                 break
             end
         end
     end
     
    c_nz = c(c~=0);                   % remove zeroes for PE calculation
    p_nz = c_nz/sum(c_nz);            % probability dist with zeroes removed
  
    pe = -sum(p_nz .* log(p_nz));     % permutation entropy
    PE(t) = pe/log(factorial(m));     % normalised permutation entropy
    
    p = c/sum(c);                     % probability distribution for SCM

    % Construct uniform probability distribution
    pu(1:length(p)) = 1/factorial(m);

    % Construct probability distribution for max JS divergence
    pmax(1:length(pu)) = zeros;
    pmax(1) = 1;
    
    %======================================================================
    % Jensen-Shannon Divergence
    JSD(t) = (-sum(((p+pu)/2) .*log((p+pu)/2))) - (-sum(p_nz .*log(p_nz)))/2 - (-sum(pu .*log(pu)))/2;
    % Normalisation factor
    Qnorm(t) = (-sum(((pmax+pu)/2) .*log((pmax+pu)/2))) - 0 - (log(factorial(m)))/2;
%     N = factorial(m);
%     Qnorm2(t) = -2*(((N+1)/N)*log(N+1)-2*log(2*N)+log(N))^-1;   % from Masoller's paper
    % Disequilibrium
    DISEQ(t) = (1/Qnorm(t))*JSD(t);
%     DISEQ2(t) = (1/Qnorm2(t))*JSD(t);
    % Statistical complexity measure
    SCM(t) = DISEQ(t)*PE(t);
%     SCM2(t) = DISEQ2(t)*PE(t);
    %======================================================================

    waitbar(t/td,h,'Calculating...');
    
end
tau = 1:td;
close(h)
 