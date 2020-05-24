%%# This script analyze mass conservation of test case 1

fig = figure();
load('test1a');
plot(MP.outputTime, sum(resultsMass,2)','linewidth',2);
hold on;
load('test1b');
plot(MP.outputTime, sum(resultsMass,2)','linewidth',2);

ax = gca();
axis(ax,'square')
xlim([0 1000]);
ylim([0,0.6e9]);
set(ax,'fontsize',16);
legend({'$R_{1}=5\times10^{5}$','$R_{1}=0$'},'box','off','location','northwest','fontsize',16,'interpreter','latex');

print(fig,'test1_mass_conservation','-dtiff','-r600');
 