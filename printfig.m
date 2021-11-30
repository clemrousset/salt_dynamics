function printfig(dir,name)
%set(findobj('tag','m_grid_color'),'facecolor','none')  
eval(['print -depsc ',dir,name]);
eval(['print -dpng -r400 ',dir,name]);
%eval(['export_fig ',dir,name,' -eps']);

