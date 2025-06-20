function hfig = stdplot(i,config)

hfig=figure(i);
clf(hfig);

if (config == 1)
    fontsize  = 6;%14;
    linewidth = 1;%2;
    dim = [7.5 5.5 6.25 2.5];
    set(hfig,'Units','Inches','Position',dim);
    set(hfig,'PaperUnits',get(gcf,'Units'));
    pos = get(hfig,'Position');
    set(hfig,'PaperPosition',[0 0 pos(3) pos(4)]);
    set(gca,'Units','Inches');
elseif (config == 2)
    fontsize  = 14;
    linewidth = 2;
else
    error('config must be 1 or 2')
end

set(hfig,'DefaultAxesFontName','Helvetica');
set(hfig,'DefaultTextFontName','Helvetica'); 
set(hfig,'DefaultAxesFontSize',fontsize);
set(hfig,'DefaultTextFontSize',fontsize);
set(hfig,'DefaultLineLineWidth',linewidth)
set(hfig,'DefaultLineLineWidth',linewidth)

end