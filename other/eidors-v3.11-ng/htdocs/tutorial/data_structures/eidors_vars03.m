% Show EIDORS colours $Id: eidors_vars03.m 2757 2011-07-14 15:56:06Z bgrychtol $

clf;
% Set square figure and make figure fill the axis
axes('position',[0 0 1 1]);
pp= get(gcf,'paperposition');
set(gcf,'paperposition',[pp(1:3),pp(3)]);

calc_colours('greylev',.001); % black background level
show_slices(img);
print_convert eidors_vars03a.png '-density 20'

calc_colours('greylev',.2); % grey background level
show_slices(img);
print_convert eidors_vars03b.png '-density 20'

calc_colours('greylev',-.2); %light grey background level
show_slices(img);
print_convert eidors_vars03c.png '-density 20'

calc_colours('greylev',-.001); %white background level (default)
show_slices(img);
print_convert eidors_vars03d.png '-density 20'


calc_colours('backgnd',[0.2,0.1,0.15]);
show_slices(img);
calc_colours('backgnd',[0.5,0.5,0.15]); %default value
print_convert eidors_vars03e.png '-density 20'
set(gcf,'paperposition',pp(1:4));
