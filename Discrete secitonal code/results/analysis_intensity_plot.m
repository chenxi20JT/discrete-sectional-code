% This plot shows how to do a contour plot of the cluster size 
%  distribution

clearvars;
clc;

% import data
data = load('test0');
t = linspace(0, data.MP.totalTime, data.MP.nStep) / 3600; % in h
t = data.MP.outputTime/3600;
N = data.resultsNum; % in #/cm3
N(N<0) = 0;
dp = data.MP.dpBins *1e9;   % in nm
I = data.MP.ND; % number of discrete bins
n = N*0;
sLim = data.MP.sLim;    % right limit of each bin, in molecule number
for ii = 1:length(t)
    n(ii,:) = N(ii,:) .* data.MP.convertNumToLog; % convert to dn/dlogdp
end

% use pcolor to plot the particle size distribution at simulation output time
n_plot = n(:,2:end)';
[T,DP]=meshgrid(t,dp(2:end));
p1 = pcolor(T,DP,log10(n_plot));

load('MyColormaps','mycmap');
set(gcf,'Colormap',mycmap,'Units','centimeters','Position',[1 2 30 15]);
grid off;
shading flat;
cb = colorbar;
axis([t(1), t(end), 1, 1e3]);
ax = get(p1,'Parent');
set(ax,'LineWidth',1.5);
xlabel('Time (hour of day)');
ylabel('Particle diameter (nm)');
caxis([0,10]);
ylabel(cb,'d\itN\rm/dlog\itd_p \rm (#/cm{\bf\fontsize{18}^3})','FontName','Arial','FontSize',20);
set(gca,'FontSize',18,'FontName','Arial','YScale','log','box','off','ticklength',[0.007 0.007], 'TickDir', 'out');
set(get(gca,'XLabel'),'FontSize',20,'FontName','Arial');
set(get(gca,'YLabel'),'FontSize',20,'FontName','Arial');

