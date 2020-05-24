function [wall_coeffs,l0_coeffs,e0_coeffs] = calsink(NT,ND,slim,slim0,slen,A,kernel,c_kernel,binsPer2,filePrefix)

currentDir = pwd;

wall_coeffs = zeros(NT,1);
l0_coeffs = zeros(NT,1);
e0_coeffs = zeros(NT,1);

%set wall loss coefficients
if ND > 0
    wall_coeffs(1:ND) = 1./((1:ND).^(1/3));
end
if ND+1 <= NT
    wall_coeffs(ND+1) = 1.5/slen(ND+1) * (slim(ND+1)^(2/3) - slim0^(2/3));
end
if ND+2 <= NT
    wall_coeffs(ND+2:NT) = 1.5./slen(ND+2:NT).*...
                            (slim(ND+2:NT).^(2/3) - slim(ND+1:NT-1).^(2/3));
end
disp('complete wall loss parameters');

%set loss coefficients to preexisting particles
if ND > 0
    l0_coeffs(1:ND) = 1./((1:ND).^(1/2));
end
if ND+1 <= NT
    l0_coeffs(ND+1) = 2/slen(ND+1) * (slim(ND+1).^(1/2) - slim0.^(1/2));
end
if ND+2 <= NT
    l0_coeffs(ND+2:NT) = 2./slen(ND+2:NT).*...
                            (slim(ND+2:NT).^(1/2) - slim(ND+1:NT-1).^(1/2));
end
disp('complete loss-to-preexisting-particles parameters');

% %set monomer mass gain parameters due to evaporation
e0_coeffs(1) = 0;
%discrete
if ND > 1
    e0_coeffs(2) = 2/2*calevap(c_kernel,A,2);
end
if ND > 2
    for i = 3:ND
        e0_coeffs(i) = 1/i*calevap(c_kernel,A,i);
    end
end
%sectional
if ND+1 <= NT
    e0_coeffs(ND+1) = integral(@(x) 1*calevap(c_kernel,A,x)./x,...
                                slim0,slim(ND+1),'AbsTol',0,'RelTol',1e-6)/slen(i); 
end
if ND+2 <= NT
    for i = ND+2:NT
        e0_coeffs(i) = integral(@(x) 1*calevap(c_kernel,A,x)./x,...
                                    slim(i-1),slim(i),'AbsTol',0,'RelTol',1e-6)/slen(i); 
    end
end
disp('complete evaporation');
filename = [currentDir,'\','coefficients\',filePrefix,'sinkParams_',kernel,'_',num2str(ND),'_',num2str(NT),'_',num2str(binsPer2)];

save(filename,'wall_coeffs','l0_coeffs','e0_coeffs');
end
