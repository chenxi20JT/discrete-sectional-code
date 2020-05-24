%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function calculate the mass based coagulation coefficient for every
% two size bin
% Discrete and Sectional bins are calculated separately but not
% distinguished in the output matrix (because ODE do not distinguish them)
% Output: coagSnkCoefs, coagSrcIds, coagSrcCoefs
% coagSnkCoefs:stores aij described in section 2.2.2 of the tutorial
% coagScrIds: stores u described in section 2.2.2 of the tutorial 
% coagSrcCoefs: stores bij->u described in section 2.2.2 of the tutorial
% Note that coagSnkCoefs also incoporates condensation described in 
% section 2.2.3 of the tutorial to make the code compact
% Last update: 2020-04-28
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [coagSnkCoefs, coagSrcIds, coagSrcCoefs] = calcoagcoeffs(ND, DIS_SECT, sLim, sLim0, sLen, kappa, binsPer2, KERNEL, coag_kernel, GROW_BEY_BOUND, filePrefix)
    NT = length(sLim);
    % initialize the output variables
    [coagSnkCoefs, coagSrcIds, coagSrcCoefs] = deal(zeros(NT, NT));
    
    % D-D coagulation
    if (DIS_SECT == 1) || (DIS_SECT == 3)   % there are discrete bins
        fprintf('Calculating D-D coagulation...');
        for ii = 1:ND
            for jj = 1:ii
                beta = coag_kernel(ii,jj);
                coagSnkCoefs(ii,jj) = beta/jj; % loss for ii bin        
                coagSnkCoefs(jj,ii) = beta/ii; % loss for jj bin

                nk = ii + jj;
                % determine the corresponding bin of the new particle
                idk = find(nk <= sLim, 1);   % the right limit belong to a given sectional bin, but not the left limit
                if ~isempty(idk) % within the whole size range?
                    coagSrcIds(ii,jj) = idk;
                    coagSrcCoefs(ii,jj) = coagSnkCoefs(ii,jj) + coagSnkCoefs(jj,ii);
                else
                    if GROW_BEY_BOUND == false % remain in the largest bin?
                        coagSrcIds(ii,jj) = NT;
                        coagSrcCoefs(ii,jj) = coagSnkCoefs(ii,jj) + coagSnkCoefs(jj,ii);
                    end % else coagSrcIds(ii,jj) remains 0
                end
                coagSrcIds(jj,ii) = coagSrcIds(ii,jj);
                coagSrcCoefs(jj,ii) = coagSrcCoefs(ii,jj);
            end % end of jj loop
        end % end of ii loop
        fprintf('\t\tcompleted\n');
    end % end of discrete bin check

    % D-S coagulation
    if (DIS_SECT == 1) || (DIS_SECT == 2)   % there are sectional bins
        fprintf('Calculating D-S coagulation...');
        % condensation: ii = 1
        infoDisSec = [num2str(1),'/',num2str(ND)];
        fprintf(infoDisSec);
        for jj = ND+1:NT
            if jj-1 == ND
                coagSnkCoefs(1,jj) = integral(@(nx) coag_kernel(1,nx)./nx, sLim0, sLim(jj), 'AbsTol',0,'RelTol',1e-6) / sLen(jj);
                coagSnkCoefs(jj,1) = integral(@(nx) coag_kernel(1,nx)/1, sLim0, sLim(jj), 'AbsTol',0,'RelTol',1e-6) / sLen(jj);
            else
                coagSnkCoefs(1,jj) = integral(@(nx) coag_kernel(1,nx)./nx, sLim(jj-1), sLim(jj), 'AbsTol',0,'RelTol',1e-6) / sLen(jj);
                coagSnkCoefs(jj,1) = integral(@(nx) coag_kernel(1,nx)/1, sLim(jj-1), sLim(jj), 'AbsTol',0,'RelTol',1e-6) / sLen(jj);
            end
            % condensation source coefficient
            coagSrcIds(1,jj) = jj;
            if (jj == NT) && (GROW_BEY_BOUND == false)
                coagSrcCoefs(1,jj) = coagSnkCoefs(1,jj) + coagSnkCoefs(jj,1);
            else
                coagSrcCoefs(1,jj) = coagSnkCoefs(1,jj) + coagSnkCoefs(jj,1) - coagSnkCoefs(1,jj)/(1-1/kappa);
            end
            coagSrcIds(jj,1) = coagSrcIds(1,jj);
            coagSrcCoefs(jj,1) = coagSrcCoefs(1,jj);
        end
        if DIS_SECT == 1    % there are both discrete and sectional bins
            % coagulation
            for ii = 2:ND
                for nn = 1:length(infoDisSec)
                    fprintf('\b');
                end
                infoDisSec = [num2str(ii),'/',num2str(ND)];
                fprintf(infoDisSec);
                for jj = ND+1:NT
                    if jj-1 == ND
                        coagSnkCoefs(ii,jj) = integral(@(nx) coag_kernel(ii,nx)./nx, sLim0, sLim(jj), 'AbsTol',0,'RelTol',1e-6) / sLen(jj);
                        coagSnkCoefs(jj,ii) = integral(@(nx) coag_kernel(ii,nx)/ii, sLim0, sLim(jj), 'AbsTol',0,'RelTol',1e-6) / sLen(jj);
                    else
                        coagSnkCoefs(ii,jj) = integral(@(nx) coag_kernel(ii,nx)./nx, sLim(jj-1), sLim(jj), 'AbsTol',0,'RelTol',1e-6) / sLen(jj);
                        coagSnkCoefs(jj,ii) = integral(@(nx) coag_kernel(ii,nx)/ii, sLim(jj-1), sLim(jj), 'AbsTol',0,'RelTol',1e-6) / sLen(jj);
                    end
                    % determine the size range of new particles, (nMin, nMax], and find the new bin(s)
                    if jj-1 == ND
                        nMin = ii + sLim0;
                    else
                        nMin = ii + sLim(jj-1);
                    end
                    nMax = ii + sLim(jj);
                    idk = find(nMin < sLim, 1);
                    % coagulation source coefficient
                    if ~isempty(idk)    % within the whole size range?
                        coagSrcIds(ii,jj) = idk;
                        if (nMax <= sLim(idk)) || ...   % all the new particles locates in one size bin
                                (idk == NT) && (GROW_BEY_BOUND == false) % or they are force to be in one size bin
                            coagSrcCoefs(ii,jj) = coagSnkCoefs(ii,jj) + coagSnkCoefs(jj,ii);
                        else    % new particles are distributed over two bins
                            if jj-1 == ND
                                coagSrcCoefs(ii,jj) = integral(@(nx) coag_kernel(ii,nx).*(ii+nx)/ii./nx, sLim0, sLim(idk)-ii, 'AbsTol',0,'RelTol',1e-6) / sLen(jj);
                            else
                                coagSrcCoefs(ii,jj) = integral(@(nx) coag_kernel(ii,nx).*(ii+nx)/ii./nx, sLim(jj-1), sLim(idk)-ii, 'AbsTol',0,'RelTol',1e-6) / sLen(jj);
                            end
                        end
                    else
                        if GROW_BEY_BOUND == false
                            coagSrcIds(ii,jj) = NT;
                            coagSrcCoefs(ii,jj) = coagSnkCoefs(ii,jj) + coagSnkCoefs(jj,ii);
                        end
                    end % if idk == nothing and GROW_BEY_BOUND = true, only sink, no source
                    coagSrcIds(jj,ii) = coagSrcIds(ii,jj);
                    coagSrcCoefs(jj,ii) = coagSrcCoefs(ii,jj);
                end % end of jj loop
            end % end of ii loop
            for nn = 1 : length('Calculating D-S coagulation...')+length(infoDisSec)
                fprintf('\b');
            end
        end % end of discrete and sectional bins check
        fprintf('Calculating D-S coagulation...\t\tcompleted\n');
    end % end of sectional bins check
    
    % S-S coagulation
    if (DIS_SECT == 1) || (DIS_SECT == 2)   % there are sectional bins
        fprintf('Calculating S-S coagulation...');
        for ii = ND+1:NT % upper triangle
            if ii > ND+1
                for nn = 1:length(infoSecSec)
                    fprintf('\b');
                end
            end
            infoSecSec = [num2str(ii-ND),'/',num2str(NT-ND)];
            fprintf(infoSecSec);
            for jj = ND+1:ii
                % in this function, m: mass, d: diameter, n: molecule number, k: new particles, id: bin number, beta: coagulation coefficient
                if ii-1 == ND % then jj-1 == ND
                    coagSnkCoefs(ii,jj) = TwoD(@(nx,ny) coag_kernel(nx,ny)./ny,...
                        sLim0,sLim(ii),sLim0,sLim(jj), 'AbsTol',0,'RelTol',1e-6) /sLen(ii)/sLen(jj);
                elseif (ii-1 > ND) && (jj-1 == ND)
                    coagSnkCoefs(ii,jj) = TwoD(@(nx,ny) coag_kernel(nx,ny)./ny,...
                        sLim(ii-1),sLim(ii),sLim0,sLim(jj), 'AbsTol',0,'RelTol',1e-6) /sLen(ii)/sLen(jj);
                    coagSnkCoefs(jj,ii) = TwoD(@(nx,ny) coag_kernel(nx,ny)./nx,...
                        sLim(ii-1),sLim(ii),sLim0,sLim(jj), 'AbsTol',0,'RelTol',1e-6) /sLen(ii)/sLen(jj);
                else
                    coagSnkCoefs(ii,jj) = TwoD(@(nx,ny) coag_kernel(nx,ny)./ny,...
                        sLim(ii-1),sLim(ii),sLim(jj-1),sLim(jj), 'AbsTol',0,'RelTol',1e-6) /sLen(ii)/sLen(jj);
                    coagSnkCoefs(jj,ii) = TwoD(@(nx,ny) coag_kernel(nx,ny)./nx,...
                        sLim(ii-1),sLim(ii),sLim(jj-1),sLim(jj), 'AbsTol',0,'RelTol',1e-6) /sLen(ii)/sLen(jj);
                end
                % determine the molecule number range of new particles, (nMin, nMax], and find the new bins
                if ii-1 == ND
                    nMin = 2*sLim0;
                elseif (ii-1 > ND) && (jj-1 == ND)
                    nMin = sLim(ii-1) + sLim0;
                else
                    nMin = sLim(ii-1) + sLim(jj-1);
                end
                nMax = sLim(ii) + sLim(jj);
                idk = find(nMin < sLim, 1);
                % coagulation source coefficient
                if ~isempty(idk)    % within the whole size range?
                    coagSrcIds(ii,jj) = idk;
                    if (nMax <= sLim(idk)) || ((idk == NT)&&(GROW_BEY_BOUND == false))
                        coagSrcCoefs(ii,jj) = coagSnkCoefs(ii,jj) + coagSnkCoefs(jj,ii);
                    else    % new particles are distributed over two bins
                        if ii-1 == ND
                            coagSrcCoefs(ii,jj) = TwoD(@(nx,ny) coag_kernel(nx,ny).*(nx+ny)./nx./ny,...
                                sLim0, min(sLim(idk)-sLim0,sLim(ii)), sLim0, min(sLim(idk)-sLim0,sLim(jj)), 'AbsTol',0,'RelTol',1e-6) /sLen(ii)/sLen(jj);
                        elseif (ii-1 > ND) && (jj-1 == ND)
                            coagSrcCoefs(ii,jj) = TwoD(@(nx,ny) coag_kernel(nx,ny).*(nx+ny)./nx./ny,...
                                sLim(ii-1), min(sLim(idk)-sLim0,sLim(ii)), sLim0, min(sLim(idk)-sLim(ii-1),sLim(jj)), 'AbsTol',0,'RelTol',1e-6) /sLen(ii)/sLen(jj);
                        else
                            coagSrcCoefs(ii,jj) = TwoD(@(nx,ny) coag_kernel(nx,ny).*(nx+ny)./nx./ny,...
                                sLim(ii-1), min(sLim(idk)-sLim(jj-1),sLim(ii)), sLim(jj-1), min(sLim(idk)-sLim(ii-1),sLim(jj)), 'AbsTol',0,'RelTol',1e-6) /sLen(ii)/sLen(jj);
                        end
                    end
                else
                    if GROW_BEY_BOUND == false
                        coagSrcIds(ii,jj) = NT;
                        coagSrcCoefs(ii,jj) = coagSnkCoefs(ii,jj) + coagSnkCoefs(jj,ii);
                    end
                end % if idk == nothing and GROW_BEY_BOUND = true, only sink, no source
                coagSrcIds(jj,ii) = coagSrcIds(ii,jj);
                coagSrcCoefs(jj,ii) = coagSrcCoefs(ii,jj);
            end % end of jj loop
        end % end of ii loop
        for nn = 1 : length('Calculating S-S coagulation...')+length(infoSecSec)
            fprintf('\b');
        end
        fprintf('Calculating S-S coagulation...\t\tcompleted\n');
    end % end of sectional bin check
    
    filename = [pwd,filesep,'coefficients',filesep,filePrefix, 'coag_',KERNEL,'_',num2str(ND),'_',num2str(NT),'_',num2str(binsPer2)];
    if exist([pwd,filesep,'coefficients'])      
        save(filename, 'coagSnkCoefs', 'coagSrcIds','coagSrcCoefs');
    else
        mkdir([pwd,filesep,'coefficients']);
        save(filename, 'coagSnkCoefs', 'coagSrcIds','coagSrcCoefs');
    end
    
end