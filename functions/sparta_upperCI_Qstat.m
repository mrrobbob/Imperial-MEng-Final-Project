function CI_UpperLim = sparta_upperCI_Qstat(E)
%
% INPUT
%   E:              model residuals
%
% OUTPUT
%   CI_UpperLim:    upper confidence limit

% Christopher Li (christopher.li19@imperial.ac.uk) & Carl Emil Eskildsen
% (c.eskildsen@imperial.ac.uk)
% Imperial College London, 2023

% References
% Box, Some Theorems on Quadratic Forms in Applied in the Study of Analysis of Variance Problems: Effects of Inequality of Variance in One-Way Classification, The Annals of Mathematical Statistics, 1954, 25, 290-302
% Momikos & MacGregor, Multivariate SPC Charts for Monitoring Batch Processes, Technometrics, 1995, 37:1, 41-59
% Jackson & Mudholkar, Control Procedures for residuals Associated With
% Principal Component Analysis, Technometrics, 1979, 21, 341-349

[I,J] = size(E);

lampda = eig(cov(E));
lampda = sort(lampda,'descend');
% lampda = lampda(1:min([I,J])-1);

% z = 1.645;  % 90% confidence
% z = 1.96;   % 95% confidence
z = 2.576;  % 99% confedence
% alpha = 0.05;
% z = norminv(1-alpha/2);

theta = nan(3,1);
for i = 1:3
    theta(i) = sum(power(lampda,i));
end

h0 = 1-(2*theta(1)*theta(3))/(3*power(theta(2),2));
if h0<0
    z = -z;
end
beta = nan(2,1);
beta(1,1) = 1-theta(2)*h0*(1-h0)/power(theta(1),2); 
beta(2,1) = z*power(2*theta(2)*power(h0,2),0.5)/theta(1);
CI_UpperLim = theta(1)*power(beta(1)+beta(2),1/h0);

end