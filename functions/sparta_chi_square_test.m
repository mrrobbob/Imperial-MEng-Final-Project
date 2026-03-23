function [chi2_stat,df,p_value] = sparta_chi_square_test(n_success,n_traps,p0)
%
% INPUT:
%   n_successes (column vector):    each element corresponds to number of successful traps of a given group
%   n_traps (column vector):        each element corresponds to total number of traps attempts in the sparta run
%   p0 (scalar):                    probability of success under H0
%
% OUTPUT:
%   chi2_stat:                      chi-square statistics
%   df:                             degrees of freedom
%   p_value:                        probability of H0 being true 

% Christopher Li (christopher.li19@imperial.ac.uk) & Carl Emil Eskildsen
% (c.eskildsen@imperial.ac.uk)
% Imperial College London, 2023

% Reference
% Brockhoff, et al., Introduction to Statistics at DTU, Chapter 7,2018, https://02402.compute.dtu.dk/enotes/book-IntroStatistics

I = size(n_success,1);          % number of groups
df = I-1;                       % degrees of freedom

nhat_success = nan(I,1);        % expected successfull traps under H0
nhat_unsuccess = nan(I,1);      % expected unsuccessfull traps under H0
n_unsuccess = nan(I,1);         % observed unsuccessfull traps
for i = 1:I
    nhat_success(i) = n_traps(i)*p0;
    nhat_unsuccess(i) = n_traps(i)-nhat_success(i);
    n_unsuccess(i) = n_traps(i)-n_success(i);
end

obs = [n_success; n_unsuccess];             % observed values
expected = [nhat_success; nhat_unsuccess];   % expected values

chi2_stat = sum(power(obs-expected,2)./expected);
p_value = 1-chi2cdf(chi2_stat,df);
end