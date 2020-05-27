%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function calculate beta_(i,j)*exp(1.5*A*(k**(2/3)-(k-1)**(2/3)))
% which is the scaling factor of the evaporation rate,see Table S1 in tutorial
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function E_k = calevap(c_kernel,A,k)
E_k = c_kernel(1,k-1).*exp(1.5*A*(k.^(2/3)-(k-1).^(2/3)));
end