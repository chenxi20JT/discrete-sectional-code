%%# This script is the analysis file for Test 2 in the tutorial
load('test3');

beta0 = 1e-9;
RR = MP.RR; 
N1_sim = resultsNum(:,1); 
Ntot_sim = sum(resultsNum, 2);
outputTime = MP.outputTime; 

% compare the the asymptotic solution for N1
fig = figure();
N1_asy = 0.5*sqrt(2).*sqrt(RR/beta0).*tanh(outputTime*sqrt(RR*beta0)/sqrt(2))+ outputTime.*sqrt(RR*beta0)./2./(cosh(outputTime*sqrt(RR*beta0)/sqrt(2)).^2);

plot(outputTime, N1_sim', 'linestyle', '-','linewidth', 2);
hold on;
plot(outputTime, N1_asy, 'linestyle', '--','linewidth', 2);
hold on;

% compare the the asymptotic solution for Ntot
Ntot_asy = sqrt(2).*sqrt(RR/beta0).*tanh(outputTime*sqrt(RR*beta0)/sqrt(2));
plot(outputTime, Ntot_sim', 'linestyle', '-', 'linewidth', 2);
hold on;
plot(outputTime, Ntot_asy, 'linestyle', '--','linewidth', 2)

ax = gca();
legend(ax,{ 'Simulated $n_{1}$','Asymptotic $n_{1}$','Simulated $n_{tot}$','Asymptotic $n_{tot}$',},'fontsize',16,'location','southeast','interpreter','latex')
% xlabel('Time [s]'); 
% ylabel('N_{tot} [cm^{-3}');
set(ax,'fontsize',16);
axis(ax, 'square');
print(fig,'fig_test3','-dtiff','-r600');
