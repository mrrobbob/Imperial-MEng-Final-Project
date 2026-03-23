function [phat,CI] = sparta_CI_trap_prob(n,x)
%
% INPUT
%   n (column vector):  total number of trap attemps
%   x (column vector):  total number of successful traps
%
% OUTPUT
%   phat:               estimated probability of successful trap
%   CI:                 confidence interval

% Christopher Li (christopher.li19@imperial.ac.uk) & Carl Emil Eskildsen
% (c.eskildsen@imperial.ac.uk)
% Imperial College London, 2023

% Reference
% Brockhoff, et al., Introduction to Statistics at DTU, Chapter 7,2018, https://02402.compute.dtu.dk/enotes/book-IntroStatistics

I = size(x,1);

% z = 1.645;  % 90% confidence
z = 1.96;   % 95% confidence
% z = 2.576;  % 99% confedence

phat = nan(I,1);
CI = nan(I,2);

for i = 1:I
    phat(i) = x(i)/n(i);
    sphat = sqrt(phat(i)*(1-phat(i))/n(i));
    CI(i,1) = phat(i)-z*sphat;
    CI(i,2) = phat(i)+z*sphat;
end

end