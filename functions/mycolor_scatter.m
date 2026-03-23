function mycolor_scatter(x,y,c,cbar,cmax)
%
% 

% Carl Emil Eskildsen, 2011

% % perform color coding
% cmap = colormap('cool');
% bcol = [0, 0.25, 0.45]; % drak blue
% graycol = [0.49, 0.60, 0.67]; % gray blue
% gcol = [0.48, 0.72, 0]; % green
% rcol = [0.73, 0.07, 0.24]; % red
% ocol = [0.93, 0.48, 0.03]; % orange
% 
% C = [224/255,236/255,244/255; 140/255,150/255,198/255; 94/255, 60/255, 153/255];
% C = [184/255,189/255,220/255; 140/255,150/255,198/255; 94/255, 60/255, 153/255];

% xx = [0; 0.5; 1];
% cmap = interp1(xx,C,linspace(0,1,64));

% Make sure dimensions of the data match
if size(x,1) ~= size(c,1) || size(x,1) ~= size(y,1)
    error ('dimentions of x, y and c do not match')
end

cmap = colormap('cool');
% cmap = colormap;
% cmap = colormap('gray');
% cmap = 1-cmap;
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
for i = 1:size(x,1)
    %     plot(Xax(i,:),X(i,:),'color',cnew(y(i),:));
    plot(x(i),y(i),'.','color',cmap(cnew(i),:),'MarkerSize',15);
    hold on;
end
% not_so_tight
% axis square; 
grid on;

% colorbar
if nargin == 4 && ~isnumeric(cbar)
    pos = get(gca,'Position');
    colorbar
    caxis([min(c) max(c)]);
    set(gca,'Position',pos);
end

box on
hold off;
shg