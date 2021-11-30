function initfig(orient)

  if findstr(orient,'land'), orient='landscape';
  elseif findstr(orient,'por'), orient='portrait';end
  
  %figure('paperorientation',orient,'paperposition',[.02 .02 .96 .96],...
  %        'paperunits','normalized','papertype','A4');
  %figure('paperorientation',orient);
  figure('paperorientation',orient,'papertype','A4');
    
