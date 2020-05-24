%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Based on user input, this script calculates other parameters used in the 
% code 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function MP = refinemodelparams(MP) 
    if MP.ND == 0
        error('There have to be one discrete bin for monomers!');
    elseif (MP.ND > 1) && (MP.NS > 0)
        DIS_SECT = 1;   % discrete-sectional model
    elseif (MP.ND == 1) && (MP.NS > 0)
        DIS_SECT = 2;   % pure sectional model
    elseif (MP.ND > 0) && (MP.NS == 0)
        DIS_SECT = 3;   % pure discrete model
    else
        error('There have to be at least one size bin!');
    end
    MP.DIS_SECT = DIS_SECT;
    
    MP.vMono  = MP.mMono/MP.rho; % monomer volume
    MP.dMono  = (MP.vMono*6/pi)^(1/3);    % monomer diameter
    MP.A      = pi*MP.dMono^2*MP.surfaceTension/1.5/1.38e-23/MP.T; %calculate the dimensionless surface tension
    
    MP.outputTime = linspace(0, MP.totalTime,MP.nStep);
   
    MP.geoFactor = 2^(1/MP.binsPer2);
    sLim  = zeros(1, MP.NT);  %right limit of each bin, in number of moleucles
    sLen  = zeros(1, MP.NT);  %length of each bin, in number of molecules
    sAvg  = zeros(1, MP.NT);  %average particle mass in each bin, in number of molecules
    convertNumToLog = zeros(1, MP.NT);  %conversion factor from number per bin to dN/dlog10(Dp)
    
    if (DIS_SECT == 1) || (DIS_SECT == 3)
        for i = 1:MP.ND
            sLim(i) = i;
            sAvg(i) = i;
            convertNumToLog(i) = 3/log10((i+0.5)/(i-0.5));
        end
        sLen(1:MP.ND) = 1;
        MP.sLim0 = sLim(MP.ND) + 0.5;  % overwrite sLim0
    end
    if MP.NS > 0
        sLim(MP.ND+1)  = MP.sLim0*MP.geoFactor;
        sLen(MP.ND+1)  = sLim(MP.ND+1)-MP.sLim0;
        sAvg(MP.ND+1)  = sLen(MP.ND+1)/log(MP.geoFactor);
        convertNumToLog(MP.ND+1) = 3/log10(MP.geoFactor);
    end
    if MP.NS > 1
        for i = MP.ND+2 : MP.NT
            sLim(i)  = sLim(i-1)*MP.geoFactor;
            sLen(i)  = sLim(i)-sLim(i-1);
            sAvg(i)  = sLen(i)/log(MP.geoFactor);
            convertNumToLog(i) = 3/log10(MP.geoFactor);
        end
    end
    dpBins = (sAvg*MP.vMono).^(1/3).*(6/pi)^(1/3);

    MP.sLim    = sLim;
    MP.sLen    = sLen;
    MP.sAvg    = sAvg;
    MP.convertNumToLog   = convertNumToLog;
    MP.dpBins  = dpBins;

end