function mycolor(Xax,X,c,cbar,cmax)
%
% function mycolor(Xax,X,c,cbar,cmax)
%
% Input:
%   X           =   data matrix - rows colored by c
%   c           =   column vector (color scale)
%   Xax         =   scale for x-axis
%   cbar        =   optional input, if present colorbar will appear
%
%
% Output:
%   color coded plot

% Carl Emil Eskildsen, 2011

% Make sure dimensions of the data match
if size(X,1) ~= size(c,1);
    error ('dimentions of X and Y do not match')
end

if size(c,2) ~= 1
    error('Y has to be a vector');
end

% delete NaN
% if isdataset(c)
%     c = c.data;
% end
% if isdataset(X)
%     X = X.data;
% end

a = find(isnan(c));

if ~isempty(a)
    X = delsamps(X,a);
    c = delsamps(c,a);
end
% perform color coding

% cmap = colormap;
% cmap = colormap('gray');
% cmap = 1-cmap;
% bcol = [0, 0.25, 0.45]; % drak blue
% graycol = [0.49, 0.60, 0.67]; % gray blue
% gcol = [0.48, 0.72, 0]; % green
% rcol = [0.73, 0.07, 0.24]; % red
% ocol = [0.93, 0.48, 0.03]; % orange
% % 
% C = [graycol; gcol; rcol];
% x = [0; 0.5; 1];
% cmap = interp1(x,C,linspace(0,1,64));
% cmap = colormap('jet');
cmap = colormap('cool');
colormap(cmap);
maxCol = size(cmap,1)-1;


if nargin == 4 && isnumeric(cbar)
    cnew = c-min(c);        % The minimum value of y will be zero
    cnew = cnew/cbar;        % The maximum value will be one
    cnew = round(cnew*maxCol)+1;   % The values of y will be between 1 and 64 (i.e. the number of colors of the colormap)
elseif (nargin == 5 && isnumeric(cmax))
    cnew = c-min(c);        % The minimum value of y will be zero
    cnew = cnew/cmax;        % The maximum value will be one
    cnew = round(cnew*maxCol)+1;   % The values of y will be between 1 and 64 (i.e. the number of colors of the colormap)
else
    cnew = c-min(c);        % The minimum value of y will be zero
    cnew = cnew/max(cnew);        % The maximum value will be one
    cnew = round(cnew*maxCol)+1;   % The values of y will be between 1 and 64 (i.e. the number of colors of the colormap)
    
end



% cnew = c-min(c);        % The minimum value of y will be zero
% cnew = cnew/max(cnew);        % The maximum value will be one
% cnew = round(cnew*63)+1;   % The values of y will be between 1 and 64 (i.e. the number of colors of the colormap)

%make the plotting
for i = 1:size(X,1);
    hold on;
    %     plot(Xax(i,:),X(i,:),'color',cnew(y(i),:));
    plot(Xax,X(i,:),'color',cmap(cnew(i),:),'linewidth',1);
end
not_so_tight('Y')


% colorbar
if nargin == 4 && ~isnumeric(cbar)
    colorbar
    caxis([min(c) max(c)]);
end

grid on;
box on
hold off;
shg