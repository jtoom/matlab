function [PE_norm, PE, tau, N] = pe_t_SF_v2(y,m,start_td,finish_td,int_td,w)

%==========================================================================
% DESCRIPTION
% Calculates permutation entropy of a time series for a given set of tau
%
% Version 2.1.

% INPUTS
% y = timeseries to be analysed
% m = length of ordinal pattern to use
% start_td = initial time delay between points used in ordinal pattern
% finish_td = maximum time delay between points used in ordinal pattern
% int_td = Interval between time delays.
% w = Enter 0 to turn the waitbar off. Useful when running this code
% within scripts.

% OUTPUTS
% PE_norm = Normalised permutation entropy
% PE = Permutation entropy 
% tau = list of tau values.
% N = Counts for each ordinal pattern.
%==========================================================================

%Preallocate matrices;

M = factorial(m);
N = zeros(M,1);
tau = start_td:int_td:finish_td;
PE = zeros(1,length(tau));
PERM = perms(1:m);
Pairs = zeros(M*(m-1),2);
A = zeros(1,sum(1:m-1));
B = zeros(1,sum(1:m-1));
Ineq = zeros(M*(m-1),1);

if w ~= 0
    h = waitbar(0,'Tigers dont play footy real good');
end

%Encode ordinal patterns into sets of inequalities;

for i = 1:m-1                       
    Pairs((i-1)*M+1:i*M,:) = PERM(:,i:i+1);
end

if m/2 == ceil(m/2)                         %Check if m is even.
    for j = 1:(m/2-1)                       %For even m.
        A(1+m*(j-1):m*j) = 1:m;
        B(1+m*(j-1):m*j) = mod((1:m)+j-1,m)+1;
    end
    A(m*(m/2-1)+1:end) = 1:m/2;
    B(m*(m/2-1)+1:end) = (m/2+1):m;
else
    for j = 1:(m-1)/2                       %For odd m.
        A(1+m*(j-1):m*j) = mod(j*(1:m),m)+1;
        B(1+m*(j-1):m*j) = mod(j*((1:m)+1),m)+1;
    end
end

for k = 1:sum(1:m-1)
    Hit = (A(k) == Pairs(:,1)) & (B(k) == Pairs(:,2));
    Hit_Inv = (B(k) == Pairs(:,1)) & (A(k) == Pairs(:,2));
    Ineq(Hit) = k;                          %Encode matching pairs as k.
    Ineq(Hit_Inv) = -1*k;                   %Encode a match in the inverse order as -k.
end

RIneq = reshape(Ineq,M,[]);

%Calculate PE;

for k = 1:length(tau)

    %Compile time series into a single matrix;
    
    Y = zeros(m,length(y)+(1-m)*tau(k));
    Logic_Y = false(sum(1:m-1),size(Y,2));      %Use logical data type for faster speeds.
    
    for i = 1:m
        Y(i,:) = y((1+(i-1)*tau(k)):end-(m-i)*tau(k))';
    end
    
%     Y = Y + 0.1*randn(size(Y));                 %Add random perturbations.
    
    %Take all sum(1:m-1) logical operations;
    
    for ii = 1:sum(1:m-1)
        Logic_Y(ii,:) = Y(A(ii),:) > Y(B(ii),:);
    end
    
    Logic_Ynot = not(Logic_Y);
    
    %Evaluate population of each ordinal pattern using the required logical
    %operations;
    
    for jj = 1:M
        Match = 1;
        for kk = 1:m-1
            if RIneq(jj,kk) > 0
                Match = and(Match,Logic_Y(RIneq(jj,kk),:));
            else
                Match = and(Match,Logic_Ynot(abs(RIneq(jj,kk)),:));
            end
        end
        N(jj) = sum(Match);
    end

    PE(k) = -1*sum(N/sum(N).*log(N/sum(N)));    %If there are any N = 0, then the code probably isn't working as it should;
    
    if w ~= 0
        waitbar(k/length(tau),h,strcat('Calculating (',num2str(k/length(tau)*100),' %)'));
    end
    
end

PE_norm = PE/log(M);

if w ~= 0
    close(h)
end

end
 