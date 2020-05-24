%%# This script is an example of post analysis
%   It shows the dimensionless particle size distribuiton at times
%   specified by the variable time_list

load('test0');

dpBins = MP.dpBins; 
sAvg    = MP.sAvg;
convertNumToLog   = MP.convertNumToLog;
totalTime  = MP.totalTime;
nStep   = MP.nStep;
NT      = MP.NT;

% convert the distribution to dN/dlogdp 
dN_dlogdp = zeros(nStep,NT);
for i = 1:nStep
    dN_dlogdp(i,:) = resultsMass(i,1:NT)./sAvg.*convertNumToLog;
end

%specify time to view distribution
time_list = [500 1000 2500];


% plotting
lineStyleCell = {'-','--','-.','--','-','--','-','--','-','--'};
lineColorCell = {[253 51 51],...
    [51 153 51],...
    [0 51 204],...
    [255 153 51],...
    [101 101 0],...
    [0 230 230],...
    [0 0 0],...
    [152 153 102],...
    [102 0 51]};

fig = figure();
%plot the distribution
for i = 1:length(time_list)
    timeStep = round(time_list(i)/(totalTime)*(nStep-1))+1;
    h = loglog(dpBins*1e9, dN_dlogdp(timeStep,:),...
        'linestyle',lineStyleCell{i},...
        'color',lineColorCell{i}/255,...
        'linewidth',2);
    hold on;
end

ax = gca();
font_size = 22;
set(ax,'fontsize',font_size)
xrange = [0.5 2E2];
yrange = [1e3 1e8];
xlim(xrange);
ylim(yrange);
x_ticks = [1.0 10.0 100.0 1000.0];
set(ax,'XTick',x_ticks);
axis square;
set(ax,'TickLength',[0.02, 0.01])
xlabel('$d_p \  [nm]$ ','fontsize',18,'interpreter','latex')
ylabel('$dN/dlogd_{p} \  [cm^{-3}]$ ','fontsize',18, 'interpreter','latex')

%print(fig,'distribuiton','-djpeg','-r600');
