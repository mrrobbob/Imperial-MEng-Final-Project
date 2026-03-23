function [X,ID, I, it] = sparta_import_shortacq(pat)

current_pat = cd;

cd(pat)                     % set path
fdir = dir('*.csv');        % get all csv files in the folder

I = size(fdir,1);                       % number of traps (i.e., number of csv files in folder)
Mtemp = readmatrix(fdir(1).name);       % read first csv file to extract number of short acquisitions per trap
it = size(Mtemp,2)-1;                   % number of short acquisition measurements per trap (i.e., number of measurements per csv file)

J = 2000;           % number of variables (Raman shift, wavenumbers)
N = I*it;           % total number of measurements
X = nan(N,J);       % preallocation of matrix for spectral data
ID = nan(I*it,1);   % preallocatre vector for ID

X(1:it,:) = Mtemp(2:end,2:end)';                            % read first file
ID(1:it) = repmat(str2num(fdir(1).name(1:end-4)),it,1);     % ID of first file
for i = 2:I
    Mtemp = readmatrix(fdir(i).name);
    X(it*(i-1)+1:it*i,:) = Mtemp(2:end,2:end)';
    ID(it*(i-1)+1:it*i) = repmat(str2num(fdir(i).name(1:end-4)),it,1);
end
[~,idx] = sort(ID); % idx for sorting
ID = ID(idx);
X = X(idx,:);

cd(current_pat);
end