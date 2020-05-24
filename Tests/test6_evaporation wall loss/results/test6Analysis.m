%% This is the analysis file for test case 5

dataNames = {'test6a.mat','test6b.mat','test6c.mat','test6d.mat'};
fig = figure();
colorCell  = {[0, 0.4470, 0.7410],[0.4660, 0.6740, 0.1880],[0.8500, 0.3250, 0.0980],[0.4940, 0.1840, 0.5560]};
styleCell  = {'-','--','-','--'};

% plot simulation results
for i = 1:length(dataNames)
    load(dataNames{i});
    loglog(MP.dpBins*1e9, resultsMass(end,:)./MP.sAvg.*MP.convertNumToLog, 'linestyle',styleCell{i},'linewidth',1.5, 'color',colorCell{i});
    hold on;
end
xlim([0.5,100]);
ylim([1e4,1e10]);
legend({'Kinetically controlled','Wall loss only', 'Evaporation only','Wall loss and evaporation'} ,'fontsize',14, 'position',[0.48,0.7,0.2,0.2],'box','off')
ax = gca();
set(ax, 'fontsize',16);
set(ax, 'ticklength',[0.02 0.01])
axis square;
% xlabel(ax,'Paricle diameter [nm]');
% ylabel(ax,'Size distribution, d\itN\rmdlog\itd\rm_p [cm^{-3}]');
print(fig,'fig_test6_evaporation_wallLoss','-dtiff','-r600');
