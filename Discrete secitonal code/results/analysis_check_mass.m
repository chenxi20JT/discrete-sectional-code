%%# This script is an example of post analysis
%   It plots total system mass as a function of time

load('test0');

% convert dp to dimensional diameter
RR = MP.RR; 
outputTime = MP.outputTime;

plot(outputTime, sum(resultsMass,2)');
xlabel('times [s]');
ylabel('total mass');
%print(fig,'distribuiton','-djpeg','-r600');
 