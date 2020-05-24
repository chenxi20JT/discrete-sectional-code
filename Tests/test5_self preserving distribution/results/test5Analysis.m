%% This part is used to compare to the 'real' self-preserving distribution

load('test5.mat');

sLim      = MP.sLim;
totalTime = MP.totalTime;
sAvg      = MP.sAvg;
sLen      = MP.sLen;
nStep     = MP.nStep;

time_list = [1000,100,20]; % in seconds
markerCell = {'s','o','<'};
colorCell  = {[0, 0.4470, 0.7410],[0.4660, 0.6740, 0.1880],[0.8500, 0.3250, 0.0980]};
spacing = [5,2,1];

fig = figure();
filename = 'self_preserving_GR.xlsx';
selfGR  = xlsread(filename); % load self preserving distribution calculated by Graham and Robinson https://doi.org/10.1016/0021-8502(76)90041-0
etaSelf = selfGR(:,1);
psiSelf = selfGR(:,2);
semilogx(etaSelf, psiSelf,'linewidth',2,'color',[0.6 0.6 0.6]);
hold on;

for i = 1:length(time_list)
    time = time_list(i);
    timeStep = round(time/(totalTime)*(nStep-1))+1;    
    numDist = resultsNum(timeStep,:);
    volDist = resultsMass(timeStep,:); % volume distribution is equavalent to mass distribution
    Ntot    = sum(numDist);
    Vtot   = sum(volDist);
    eta = sAvg*Ntot./Vtot;
    psi = numDist./sLen./Ntot^2*Vtot;
    semilogx(eta(1:spacing(i):end), psi(1:spacing(i):end),'linestyle','none','linewidth',2, 'marker',markerCell{i},'color',colorCell{i});
    hold on;
end

xlim([1e-4 1e1]);
ylim([0 1.2]);
legend({ 'self-preserving','1000 s', '100 s', '20 s'},'fontsize',12,'location','northwest');
ax = gca();
set(ax, 'fontsize',16);
% xlabel(ax,'\eta');
% ylabel(ax,'\psi (\eta)');
axis(ax,'square');
set(ax, 'xtick',[1e-3 1e-2 1e-1 1 1e1]);

print(fig,'fig_test5_self_preserving','-dtiff','-r600');
