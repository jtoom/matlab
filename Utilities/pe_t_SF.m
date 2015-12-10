function [PE,tau] = pe_t_SF(y,m,start_td,finish_td)

%==========================================================================
% DESCRIPTION
% Calculates permutation entropy of a time series for a given set of tau

% INPUTS
% y = timeseries to be analysed
% m = length of ordinal pattern to use
% start_td = initial time delay between points used in ordinal pattern
% finish_td = maximum time delay between points used in ordinal pattern

% OUTPUTS
% PE = permutation entropy
% tau = list of tau values
%==========================================================================

ly = length(y);
permlist = perms(1:m);      % all possible permutations of length m
h = waitbar(0,'Please wait...');

% Initialise output arrays
PE = zeros(finish_td-start_td+1,1);

% Initialise waitbar counter
count = 0;

for t = start_td:finish_td
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
    PE(count+1) = pe/log(factorial(m));     % normalised permutation entropy

    waitbar(count/(length(start_td:finish_td)),h,'Calculating...');
    count = count+1;
    
end
tau = start_td:finish_td;
close(h)
 