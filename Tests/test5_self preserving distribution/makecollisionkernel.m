function h = makecollisionkernel(kerneltype,T,P,rho,mMono)

boltz = 1.38e-23;
v_mono = mMono/rho;

switch kerneltype
    case 'free'
        h = @c_kernel_free;
    case 'sceats'
        h = @c_kernel_sceats;
    case 'fuchs'
        h = @c_kernel_fuchs;
    case 'hogan'
        h = @c_kernel_hogan;
    case 'uniform'
        h = @c_kernel_uniform;
    otherwise
        error('Collision frequency function not valid.');
end

    function beta = c_kernel_free(i,j)
        %This function calcuates the collision kernel in the free molecular regime
        %  two iputs i & j -- the number of monomers in the two colliding particles
        beta = (3/4/pi)^(1/6)*(6*boltz*T/rho)^0.5*(v_mono)^(1/6)*(1./i+1./j).^0.5.*(i.^(1/3)+...
            j.^(1/3)).^2.*1E6;
    end

    function beta=c_kernel_fuchs(i,j)
        
        %reference temperature, mean free path, pressure and viscosity
        T0 = 300;
        mfp0 = 67E-9;
        P0 = 1.01E5;
        mu0 = 1.82E-5;
        
        mu = mu0*(T0+110.4)/(T+110.4)*(T/T0)^(1.5);  %corrected viscosity
        mfp = mfp0*(T/T0)*(P0/P)*(1+110.4/T0)/(1+110.4/T); %corrected mean free path
        
        
        rpi = ((i.*mMono)/rho*6/pi).^(1/3)/2;
        rpj = ((j.*mMono)/rho*6/pi).^(1/3)/2; %size in m
        
        kni = mfp./rpi;
        knj = mfp./rpj;
        
        a = 1.257; %Cunningham correction coefficients
        b = 0.4;
        c = 1.1;
        
        Cc_i = 1 + kni.*(a + b*exp(-c./kni));
        Cc_j = 1 + knj.*(a + b*exp(-c./knj));
        
        Bi = Cc_i./(6*pi*mu.*rpi); %s/kg
        Bj = Cc_j./(6*pi*mu.*rpj);
        Di = Bi.*boltz*T; %m^2/s, close to Table 2.1 in Friedlander's book
        Dj = Bj.*boltz*T;
        
        mi = mMono*i;  %mass
        mj = mMono*j;  %mass
        
        ci = (8*boltz*T./(pi*mi)).^0.5; % m/s, mean thermal velocity of the particles
        cj = (8*boltz*T./(pi*mj)).^0.5;
        
        %mean traveling time
        tau_i = (i.*mMono.*Cc_i)./6./pi./mu./rpi;
        tau_j = (j.*mMono.*Cc_j)./6./pi./mu./rpj;
        
        %'mean free path' of the particles
        lbi = tau_i.*ci;
        lbj = tau_j.*cj;
        
        %detla values
        deltai = 1./6./rpi./lbi.*((2*rpi + lbi).^3 -(4*rpi.^2+lbi.^2).^1.5)-2.*rpi;
        deltaj = 1./6./rpj./lbj.*((2*rpj + lbj).^3 -(4*rpj.^2+lbj.^2).^1.5)-2.*rpj;
        
        r_mean = (rpi+rpj)./2;
        beta_corr = r_mean./(r_mean+ sqrt(deltai.^2 + deltaj.^2)/2) + 4*(Di+Dj)./2./sqrt(ci.^2+cj.^2)./r_mean;
        beta = 2.*8.*pi.*(rpi+rpj)./2.*(Di+Dj)./2./beta_corr*1e6;
    end

    function beta = c_kernel_hogan( i,j )
        %calculate the collision kernel for clusters containing i and j monomers
        %   return the value in cm3/s
        %   to accomodate for the integration, i and j are vectors
        %   using equation from Gopalakrishnan & Hogan (2011)
                
        %reference temperature, mean free path, pressure and viscosity
        T0 = 300;
        mfp0 = 67E-9;
        P0 = 1.01E5;
        mu0 = 1.82E-5;
        
        mu = mu0*(T0+110.4)/(T+110.4)*(T/T0)^(1.5);  %corrected viscosity
        mfp = mfp0*(T/T0)*(P0/P)*(1+110.4/T0)/(1+110.4/T); %corrected mean free path
        
        mi = i*mMono; %mass of the clusters
        mj = j*mMono;
        mij = mi.*mj./(mi+mj);
        ri = ((i*v_mono*6)/pi).^(1/3)/2;   %radius in m
        rj = ((j*v_mono*6)/pi).^(1/3)/2;
        
        kni = mfp./ri; %Knudsen number
        knj = mfp./rj;
        
        a = 1.257; %Cunningham correction coefficients
        b = 0.4;
        c = 1.1;
        Cci = 1+kni.*(a+b*exp(-c./kni)); %Cunningham correction
        Ccj = 1+knj.*(a+b*exp(-c./knj));
        
        fi = 3*pi*mu*2.*ri./Cci; %friction coeffcient
        fj = 3*pi*mu*2.*rj./Ccj;
        fij = fi.*fj./(fi+fj);
        
        knd = sqrt(boltz*T*mij)./fij./(ri+rj);
        
        c1 = 25.84;
        c2 = 11.21;
        c3 = 3.50;
        c4 = 7.21;
        fitting = (4*pi*knd.^2+c1*knd.^3+sqrt(8*pi)*c2*knd.^4)./(1+c3*knd+c4*knd.^2+c2*knd.^3);
        beta = (ri+rj).^3.*fij./mij.*fitting*1e6;
    end

    function beta = c_kernel_sceats( i,j )
                
        T0 = 300;
        mfp0 = 67E-9;
        P0 = 1.01E5;
        mu0 = 1.82E-5;
        
        mu = mu0*(T0+110.4)/(T+110.4)*(T/T0)^(1.5);  %corrected viscosity
        mfp = mfp0*(T/T0)*(P0/P)*(1+110.4/T0)/(1+110.4/T); %corrected mean free path
        
        mi = i*mMono; %mass of the clusters
        mj = j*mMono;
        mij = mi.*mj./(mi+mj);
        ri = ((i*v_mono*6)/pi).^(1/3)/2;   %radius in m
        rj = ((j*v_mono*6)/pi).^(1/3)/2;
        A = 6.4e-20; %Hamaker constant
        
        
        a = 1.257; %Cunningham correction coefficients
        b = 0.4;
        c = 1.1;
        kni = mfp./ri; %Knudsen number
        knj = mfp./rj;
        Cci = 1+kni.*(a+b*exp(-c./kni)); %Cunningham correction
        Ccj = 1+knj.*(a+b*exp(-c./knj));
        fi = 3*pi*mu*2.*ri./Cci; %friction coeffcient
        fj = 3*pi*mu*2.*rj./Ccj;
        Di = boltz*T./fi;
        Dj = boltz*T./fj;
        
        A_prime =  A/boltz/T*4.*ri.*rj/(ri+rj).^2;
        x = log(1+A_prime);
        
        a1 = 0.0757;
        a3 = 0.0015;
        b0 = 0.0151;
        b1 = -0.186;
        b3 = -0.0163;
        E_zero = 1 + a1.*x + a3.*x.^3;
        E_inf = 1 + sqrt(A_prime/3)./(1+b0*sqrt(A_prime))+b1*x + b3*x.^3;
        
        cij = sqrt(8*boltz*T./pi./mij);
        
        c_K = pi/2*(ri+rj).^2.*cij.*E_inf;
        c_D = 2*pi*(ri+rj).*(Di+Dj).*E_zero;
        tij = c_K./2./c_D;
        %a factor of two is added comapared to the paper
        beta = 2*c_K.*(sqrt(1+tij.^2)-tij)*1e6;        
    end

    function beta = c_kernel_uniform(i,j )
        beta = i.*j./i./j*1e-9;     
    end
end
