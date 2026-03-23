function Xsavgol = sparta_SavGol_filter(X,pol,win,der)
% Savitzky-Golay filtering for SPARTA
% Xsavgol = SavitzkyGolay(spectra, poly.order, window size, derivative order);

% Carl Emil Eskildsen, 2023
% c.eskildsen@imperial.ac.uk
% Department of Materials, Faculty of Engineering, Imperial College London

% References
% A. Savitzky, M.J.E. Golay, Smoothing and differentiation of data by simplified least squares procedure, Anal. Chem. 33 (1964) 1627e1639.
% J. Steiner, Y. Termonia, J. Deltour, Smoothing and differentiation of data by simplified least squares procedure, Anal. Chem. 44 (1972) 1906e1909

[I,J] = size(X);
F = (win-1)/2;
Xtemp = nan(I,J+2*F);
[~,Jtemp] = size(Xtemp);
idx = false(Jtemp,5);
idx(1:F,1) = true;
idx(F+1:Jtemp-F,2) = true;
idx(Jtemp-F+1:Jtemp,3) = true;
idx(F+1:F+win,4) = true; % first window
idx(F+J-win+1:F+J,5) = true; % last window
xx = 1:Jtemp;
Xtemp(:,idx(:,2)) = X;

% extrapolate spectra 
for i = 1:I
    p1 = polyfit(xx(idx(:,4)),Xtemp(i,idx(:,4)),pol);
    p2 = polyfit(xx(idx(:,5)),Xtemp(i,idx(:,5)),pol);
    Xtemp(i,idx(:,1)) = polyval(p1,xx(idx(:,1)));
    Xtemp(i,idx(:,3)) = polyval(p2,xx(idx(:,3)));
end

% SavGol coefficients
W = eye(win);
s = fliplr(vander(-(win-1)/2:(win-1)/2));
S = s(:,1:pol+1);
[~,R] = qr(sqrt(W)*S,0);
G = S/(R)*inv(R)';

Xsavgol = zeros(I,Jtemp);
F1 = (win+1)/2;
F2 = -F1+1:F1-1;

for j = F1:Jtemp-(win-1)/2
   if der == 0
       z = Xtemp(:,j + F2)*G(:,1);
   else
       z = Xtemp(:,j + F2)*(der*G(:,der+1));
   end
   Xsavgol(:,j) = z;
end 
Xsavgol = Xsavgol(:,idx(:,2));
end

