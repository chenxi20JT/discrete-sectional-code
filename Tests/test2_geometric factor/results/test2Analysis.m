%% This is the analysis file for test case 2

dataNames = {'test2z.mat','test2a.mat','test2b.mat','test2c.mat','test2d.mat'};
fig = figure();

% plot simulation results
for i = 1:length(dataNames)
    load(dataNames{i});
    if i == length(dataNames)
        plot(MP.dpBins*1e9, resultsMass(end,:)./MP.sAvg.*MP.convertNumToLog, 'linestyle','--','linewidth',1.5);
    else
        plot(MP.dpBins*1e9, resultsMass(end,:)./MP.sAvg.*MP.convertNumToLog, 'linestyle','-','linewidth',1.5);
    end
    hold on;
end
xlim([3,8]);
ylim([3e6,7e6]);
legend({'\kappa = 2','\kappa = 1.148','\kappa = 1.072', '\kappa = 1.035', '\kappa = 1.018'} ,'fontsize',16, 'location','south')
ax = gca();
axis square;
set(ax, 'fontsize',16);
% xlabel(ax,'Paricle diameter [nm]');
% ylabel(ax,'Size distribution, d\itN\rmdlog\itd\rm_p [cm^{-3}]');
print(fig,'fig_test2_geoFactor','-dtiff','-r600');
