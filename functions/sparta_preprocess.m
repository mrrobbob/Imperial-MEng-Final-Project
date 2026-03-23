function Xp = sparta_preprocess(X,wn)

% Despike
% [X,~]=AutoSpike_Matrix(X,wn,1);

% remove baseline
lambda_x = 1e4;     % smoothing parameter (generally 1e5 to 1e8) 
p_x = 0.001;        % asymmetry parameter (generally 0.001)
d_x = 2;            % order of differences in penalty (generally 2)
Xbck = zeros(size(X));
for i = 1:size(X,1)
    Xbck(i,:) = [asysm(X(i,:)', lambda_x, p_x, d_x)]'; % Asymetric least squares to estimate baseline
end
X = X-Xbck;

% noise filtering
pol = 2;
win = 11;
der = 0;
X = sparta_SavGol_filter(X,pol,win,der);


Xp = X;
end