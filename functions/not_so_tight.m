function not_so_tight(ax_tight)
% the function not_so_tight sets the axis limits such
% that there is a little space around the plot. The function can be called
% with an optional input argument, which specifies whether the x-axis or
% the y-axis should be not_so_tight. When called with no input argument no 
% axis will be tight.
%
% Input:
%       ax_tight: X | Y

% Carl Emil Eskildsen

if nargin == 0;
    ax_tight = 'none';
end
ax_tight = upper(ax_tight);

axis tight              % set tight axes, in order to get the data limits
a = axis;

dx = a(2)-a(1);         % the range of x data
dy = a(4)-a(3);         % ditto y axis


w = 0.05;               % percentage of axis range to be added

% adjust the ranges
if ax_tight ~= 'Y'
    a(1) = a(1) - w*dx;
    a(2) = a(2) + w*dx;
end

if ax_tight ~= 'X'
    a(3) = a(3) - w*dy;
    a(4) = a(4) + w*dy;
end

axis(a);                % set the new ranges
end