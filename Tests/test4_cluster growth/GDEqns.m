
function dc_dt = GDEqns(t,conc,R,kappa,ND,NT,coagSnkCoefs, coagSrcIds, coagSrcCoefs,wall_coeffs,l0_coeffs,e0_coeffs,w,l0,e0,m0,CONST_MONOMER,COAG_OFF)
    dc_dt = zeros(NT,1);
    
    if COAG_OFF
        coagSnkCoefs(2:end,2:end) = 0;
        coagSrcCoefs(2:end,2:end) = 0;
    end
    % only for test 4
    coagSnkCoefs(1,1) = 0;
    coagSrcCoefs(1,1) = 0;
    
    % coagulation & concensation sink
    for i = 1:NT
        dc_dt(i) = dc_dt(i) - sum(coagSnkCoefs(i,1:end)'.*conc(1:end)) * conc(i);
    end
    % coagulation & concensation source
    for i = 1:NT
        for j = 1:i
            if i == j
                preFactor = 0.5;
            else
                preFactor = 1;
            end
            if coagSrcIds(i,j) ~= 0
                dc_dt(coagSrcIds(i,j)) = dc_dt(coagSrcIds(i,j)) + preFactor*coagSrcCoefs(i,j)*conc(i)*conc(j);
                if coagSrcIds(i,j) < NT
                    dc_dt(coagSrcIds(i,j)+1) = dc_dt(coagSrcIds(i,j)+1) + preFactor * (coagSnkCoefs(i,j)+coagSnkCoefs(j,i)-coagSrcCoefs(i,j)) * conc(i)*conc(j);
                end
            end
        end
    end
    
    % incorporate monomer production, wall loss, scavenging by pre-existing particle,evaporation
    dc_dt(1) = dc_dt(1) + R - w*wall_coeffs(1)*conc(1)-l0*l0_coeffs(1)*conc(1)- m0*conc(1)+ sum(conc.*e0_coeffs*e0);
    
    dc_dt(2) = dc_dt(2) - w*wall_coeffs(2)*conc(2)-l0*l0_coeffs(2)*conc(2)- m0*conc(2) - conc(2)*e0*e0_coeffs(2) + conc(3)*e0*e0_coeffs(3)*2;    
    for i = 3:ND-1
        dc_dt(i) = dc_dt(i) -w*wall_coeffs(i)*conc(i)-l0*l0_coeffs(i)*conc(i)- m0*conc(i)- conc(i)*e0*e0_coeffs(i)*i+conc(i+1)*e0*e0_coeffs(i+1)*i;
    end
    
    if ND < NT  %if there are sectional bins
        dc_dt(ND) = dc_dt(ND)-w*wall_coeffs(ND)*conc(ND)-l0*l0_coeffs(ND)*conc(ND)- m0*conc(ND)- conc(ND)*e0*e0_coeffs(ND)*ND+(e0*e0_coeffs(ND+1)*conc(ND+1))/(kappa-1);
        for i = ND+1:NT-1
            dc_dt(i) = dc_dt(i) -w*wall_coeffs(i)*conc(i)-l0*l0_coeffs(i)*conc(i)- m0*conc(i)-conc(i)*e0*e0_coeffs(i)*(1+1/(kappa-1))+conc(i+1)*e0*e0_coeffs(i+1)/(kappa-1);
        end
        dc_dt(NT) = dc_dt(NT) -w*wall_coeffs(NT)*conc(NT)-l0*l0_coeffs(NT)*conc(NT)- m0*conc(NT)-conc(NT)*e0*e0_coeffs(NT)*(1+1/(kappa-1));
    else % if there are only discrete bins
        dc_dt(ND) = dc_dt(ND)-w*wall_coeffs(ND)*conc(ND)-l0*l0_coeffs(ND)*conc(ND)- m0*conc(ND)- conc(ND)*e0*e0_coeffs(ND)*ND;
    end
    
    if CONST_MONOMER
        dc_dt(1) = 0;
    end
end

