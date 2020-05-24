%This is the main program for GDE simulation


function [resultsMass,resultsNum] = main(MP)
    
    %get current working directory
    currentDir = pwd;
    % unpack the fiels from MP
    v2struct(MP);
    % construct the collision kernel function
    collisionKernel = makecollisionkernel(KERNEL,T,P,rho,mMono);
    % if CAL_COEFFS== true,calculate the coagulation coefficients in the current simulation
    if CAL_COEFFS 
        [coagSnkCoefs, coagSrcIds, coagSrcCoefs] = calcoagcoeffs(ND,DIS_SECT,sLim,sLim0,sLen,geoFactor,binsPer2,KERNEL,collisionKernel,GROW_BEY_BOUND,filePrefix);
    else
        filename1 = [currentDir, filesep,'coefficients',filesep, filePrefix,'coag_',KERNEL,'_',num2str(ND),'_',num2str(NT),'_',num2str(binsPer2), '.mat'];
        if isfile(filename1)
           load(filename1);
           disp('Coagulation coefficient loaded');
        else
            warning('Coagulation coefficient not found. Recalculating...');
            [coagSnkCoefs, coagSrcIds, coagSrcCoefs] = calcoagcoeffs(ND,DIS_SECT,sLim,sLim0,sLen,geoFactor,binsPer2,KERNEL,collisionKernel,GROW_BEY_BOUND,filePrefix);
        end
    end
    % if CAL_SINK_PARAMS == true,calculate the sink parameters in the current simulation
    if CAL_SINK_PARAMS
        [wall_coeffs,l0_coeffs,e0_coeffs] = calsink(NT,ND,sLim,sLim0,sLen,A,KERNEL,collisionKernel,binsPer2,filePrefix);
    else
        filename2 = [currentDir,filesep,'coefficients',filesep, filePrefix,'sinkParams_',KERNEL,'_',num2str(ND),'_',num2str(NT),'_',num2str(binsPer2),'.mat'];
        if isfile(filename2)
           load(filename2);
        else
            warning('Sink coefficients not found. Recalculating...');
            [wall_coeffs,l0_coeffs,e0_coeffs] = calsink(NT,ND,sLim,sLim0,sLen,A,KERNEL,collisionKernel,binsPer2,filePrefix);
        end
    end

    %-------------------------COMPUTE DISTRIBUTION----------------------------%
    % create empty matrices to store output, each row is the distribution at a given output time
    resultsNum  = zeros(nStep,NT);

    %use ode15s to simulate distribution evolution
    odeFunc = @(t,initConc) GDEqns(t,initConc,RR,geoFactor,ND,NT,coagSnkCoefs, coagSrcIds, coagSrcCoefs,...
        wall_coeffs,l0_coeffs,e0_coeffs,wLoss,svgLoss,satPressure,dilution,CONST_MONOMER,COAG_OFF);
    [~,resultsMass] = ode15s(odeFunc,outputTime,initConc, odeOptions);

    %convert to nubmer of particle per bin
    for i = 1:nStep
        resultsNum(i,:) = resultsMass(i,1:NT)./sAvg;
    end

    %save the output into a folder called  results
    if exist('results') 
        save(['results' filesep resultsFileName],'MP','resultsMass','resultsNum');
    else
        mkdir([currentDir filesep 'results']);
        save(['results' filesep resultsFileName],'MP','resultsMass','resultsNum');
    end
end




