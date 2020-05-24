%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%  This script is the entry point for dimensional simulation of a homogenous unary system

%  INPUT:
%    All the input parameters are summarized in a Matlab struct named MP (short for model parameters). 
%    The fields of MP that require user input are listed below:
%       resultsFileName : filename in which the output is saved
%       filePrefix      : a string prefix that is added to filename for further distinction between simulations
%       ND              : the number of discrete bins
%       NS              : number of sectional bins
%       DIS_SECT        : 1 -> discrete-sectional model, 2 -> pure sectional model, 3 -> pure discrete model
%       sLim0           : the left limit of sectional bin. The input value is used only for the pure sectional
%                         model; otherwise, sLim0 is overwritten by ND+0.5
%       binsPer2        : (2)**(1/binsPer2) is the geometric factor between neighbouring sectional bins
%       totalTime       : total simulation time, in sec
%       nStep           : total number of output time step
%       wLoss           : the constant part of wall loss rate,in s^-1
%       scgLoss         : loss constant of scavenging by pre-existing particles
%       dilution        : dilution constant, in s^-1
%       satPressure     : saturation vapor pressure in units of #/cm^3
%       wLoss, scgLoss, dilution, satPressure correspond to the second column of Table S1 in the tutorial. 
%       KERNEL          : 'fuchs': Fuchs transition regime limiting sphere collision theory
%                         'free': Free molecular limit
%                         'hogan':https://www.tandfonline.com/doi/abs/10.1080/02786826.2011.601775
%                         'sceats': https://www.sciencedirect.com/science/article/pii/S0021850200000811 
%                         'uniform': produces a size-independent rate constant this 1e-9 cm^3 s-1
%       CAL_COEFFS      : true-> recalculate coefficients between bins, false-> load previous coefficients
%       CAL_SINK_PARAMS : true-> recalculate coefficients for sink processes, false-> load previous coefficients
%       CONST_MONOMER   : true-> the monomer concentration will remain constant throught the simulation
%       COAG_OFF        : true-> coagualtion between particles will be turned off    
%       GROW_BEY_BOUND  : true -> particles are allowed to grow out of the upper boundary; 
%                         false -> no mass loss of the whole system due to particle growth 
%       odeOptions      : solver options for ode15s
%       initConc        ; initial mass concentration
%       RR              : monomer production rate, in #/(cm^3*s)
%       T               : Temperature in K
%       P               : Pressure in Pa
%       rho             : bulk density in kg/m^3
%       mMono           : monomer mass in kg
%       surfaceTension  : suface tension of nulceation species, in N/m

%  OUTPUT:
%       MP : a Matlab struct that contains all the model parameters
%       resultsMass    : particle mass in each bin
%       resultsNum     : particle number in each bin
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% USER INPUT 
MP = struct;

% specify the output filename
MP.resultsFileName = 'test6d';
MP.filePrefix      = 'test6_';

% Setting up the simulation bins and section equaitons
MP.ND       = 200;       
MP.binsPer2 = 20;        
MP.NS       = 20*MP.binsPer2;
MP.NT       = MP.ND + MP.NS; % total number of bins
MP.sLim0    = 100;  % will be overwitten by ND+0.5 when ND > 0

% specify simulation time
MP.totalTime = 3000; 
MP.nStep     = 301; 

% specify particle sink parameters
MP.wLoss  = 1e-2;   
MP.svgLoss = 0;   
MP.dilution = 0;  
MP.satPressure  = 2e4;   

MP.CAL_COEFFS      = false;       
MP.CAL_SINK_PARAMS = false; 
MP.CONST_MONOMER   = false; 
MP.COAG_OFF        = false;  
MP.GROW_BEY_BOUND  = true;

% choose a collision kernel for particles
MP.KERNEL  = 'fuchs';  

%set convergence criteria for the GDE solover 
MP.odeOptions = odeset('RelTol',1e-6,'AbsTol',1e-6,'OutputFcn',@odephas2); 

%set the initial number distribution of particles
MP.initConc    = zeros(MP.NT,1);

%set the monomer generation term
MP.RR = 2e6;   
                          
%Physical parameter of the system
MP.T       = 293.15;   
MP.P       = 1.01e5;   
MP.rho     = 1.47e3;   
MP.mMono   = 2.4e-25;  
MP.surfaceTension = 67.5e-3; 

%Calculate other model parameters based on the user input 
MP = refinemodelparams(MP);
  
%% start simulation
[resultsMass,resultsNum] = main(MP);
                                  