function mysetsize(h,hight,width)

set(h,'Units','centimeters')
pos = get(h,'Position');
pos(1,3:4) = [width,hight];
set(gcf,'Position',pos);
end
