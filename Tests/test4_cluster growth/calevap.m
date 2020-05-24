function E_k = calevap(c_kernel,A,k)
% this function calculate beta_(i,j)*exp(1.5*A*(k**(2/3)-(k-1)**(2/3)))
E_k = c_kernel(1,k).*exp(1.5*A*(k.^(2/3)-(k-1).^(2/3)));
end