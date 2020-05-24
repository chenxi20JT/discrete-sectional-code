%% This is the analysis file for test case 4
load('test4.mat');

ND        = MP.ND;
totalTime = MP.totalTime;
nStep     = MP.nStep;
N1        = resultsNum(1,1); %the monomer concentration
N3        = resultsNum(1,3); %the initial trimer concentration

time_list = [1000,2000,3000]; % in seconds
fig = figure();
markerCell = {'s','o','<'};
colorCell  = {[0, 0.4470, 0.7410],[0.4660, 0.6740, 0.1880],[0.8500, 0.3250, 0.0980]};


% plot simulation results
for i = 1:length(time_list)
    time = time_list(i);
    timeStep = round(time/(totalTime)*(nStep-1))+1;    
    plot(1:ND, resultsNum(timeStep,:),'linestyle','-','linewidth',1.5, 'color',colorCell{i});
    hold on;
end

%overlapping with Poisson distribution
lambda = time_list*1e-9*N1; %expected number of added molecules,1e-9 is the collision frequency when the collision kernel is chosen to be 'uniform' 
x = 1:100;
% geometric standard deviation of a poission process
for i = 1:length(lambda)
    y = poisspdf(x,lambda(i))*N3; 
    plot(x+3,y,'linestyle','--','linewidth',1.5, 'color',colorCell{i},'marker',markerCell{i},'markersize',5)
    hold on;
end
legend({'t = 1000 s','t = 2000 s', 't = 3000 s', 'Poisson, \lambda = 10','Poisson, \lambda = 20','Poisson, \lambda = 30'} ,'fontsize',14, 'location','northeast')
ax = gca();
set(ax, 'fontsize',16);
% xlabel(ax,'Number of molecules');
% ylabel(ax,'Concentration [cm^{-3}]');
ylim([0, 600]);
xlim([4, 50]);
axis(ax,'square');
print(fig,'fig_test4_discrete_growth','-dtiff','-r600');
