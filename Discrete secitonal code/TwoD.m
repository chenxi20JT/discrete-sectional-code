function Q = TwoD(fun,a,b,c,d,varargin)
%TwoD Numerically evaluate integral of f(x,y) over a plane region.
%   Q = TWOD(FUN,A,B,C,D) approximates the integral
%
%                         /b    /d
%                   I =  |     |   f(x,y) dy dx
%                       /a    /c
%
%   FUN evaluates f(x,y).  It must accept arrays X and Y and return an array
%   Z = f(X,Y) of corresponding values. These arrays are all 14x14. Two
%   kinds of finite regions are allowed.  The default is a generalized
%   rectangle which is described in Cartesian coordinates (x,y) by
%      a <= x <= b -- A and B must be constants
%      c <= y <= d -- C and D can be constants or functions of x that 
%   describe the lower and upper boundaries. Functions must be vectorized.
%   Q is an approximation to the integral with an estimated bound ERRBND on 
%   the error that satisfies  ERRBND <= 1e-5.
%
%   Q = TWOD(FUN,A,B,C,D,PARAM1,VAL1,PARAM2,VAL2,...) performs the integration
%   as above with specified values of optional parameters: 
%
%   'Sector' With value TRUE, the region is a generalized sector that is 
%        described in polar coordinates (r,theta) by
%           0 <= a <= theta <= b -- A and B must be constants
%           c <= r <= d -- C and D can be constants or functions of theta
%        that describe the lower and upper boundaries. Functions must be
%        vectorized.  NOTE Polar coordinates  are used only to describe 
%        the region--the integrand is f(x,y) for both kinds of regions.
% 
%   'AbsTol', absolute error tolerance
%   'RelTol', relative error tolerance
%       TwoD attempts to satisfy ERRBND <= max(AbsTol,RelTol*|Q|). This
%       is absolute error control when |Q| is sufficiently small and
%       relative error control when |Q| is larger. A default tolerance
%       value is used when a tolerance is not specified. The default value
%       of 'AbsTol' is 1e-5. The default value of 'RelTol' is 0.
%
%   'Singular' TwoD can be applied to f(x,y) that are singular on a 
%       boundary. With value TRUE, this option causes TwoD to use
%       transformations to weaken singularities for better performance.
%
%   'Vectorized' The efficiency of TwoD depends STRONGLY on vectorization
%       of FUN, but if this option has value FALSE, TwoD calls FUN only
%       with scalar X and Y.  If f(x,y) is not expensive to evaluate, the 
%       cost may be quite acceptable.
%
%   EXAMPLES USING DEFAULT PURE ABSOLUTE ERROR TOLERANCE OF 1e-5.
%
%   #1 Integrate y*sin(x)+x*cos(y) over pi <= x <= 2*pi, 0 <= y <= pi. 
%   The true value of the integral is -pi^2.
%      Q = TwoD(@(x,y) y.*sin(x)+x.*cos(y),pi,2*pi,0,pi) 
%
%   #2 Integrate 1/( sqrt(x + y)*(1 + x + y)^2 ) over the triangle
%   0 <= x <= 1, 0 <= y <= 1 - x.  The integrand is infinite at (0,0).  
%   The true value of the integral is pi/4 - 1/2. The region can
%   be specified in more than one way. Either way the integrand is
%       fun = @(x,y) 1./( sqrt(x + y) .* (1 + x + y).^2 );
%
%   #2a The region can be described as a generalized rectangle with 
%   one side that degenerates to a point:
%       ymax = @(x) 1 - x;
%       Q = TwoD(fun,0,1,0,ymax);
% 
%   #2b It is more natural to describe this region as a sector. The
%   computation is more efficient with the 'Singular' option:
%       rmax = @(theta) 1./(sin(theta) + cos(theta));
%       Q = TwoD(fun,0,pi/2,0,rmax,'Sector',true,'Singular',true);



%   Copyright 2008 Lawrence F. Shampine

  if isa(c,'function_handle')
    phiBvar = c;
  elseif isnumeric(c)
    phiBvar = @(x) c*ones(size(x));
  else
    error('MATLAB:TwoD:invalidC','C must be a constant or a function handle')
  end
  if isa(d,'function_handle')
    phiTvar = d;
  elseif isnumeric(d)
    phiTvar = @(x) d*ones(size(x));
  else
    error('MATLAB:TwoD:invalidD','D must be a constant or a function handle')
  end

  % Set defaults and then process optional values.    
  ATOL = 1e-5;
  RTOL = 0; 
  SECTOR = false;  
  SINGULAR = false;
  VECTORIZED = true;
  parseOptions(varargin{:});
  if SECTOR
      % x = r*cos(theta), y = r*sin(theta);
      FUN = @(theta,r) fun(r.*cos(theta),r.*sin(theta)).*r;
  else
      FUN = fun;
  end 
  if SINGULAR
      thetaL = 0; 
      thetaR = pi;
      phiB = 0;   
      phiT = pi; 
  else
      thetaL = a; 
      thetaR = b;
      phiB = 0;   
      phiT = 1;
  end
  AREA = (thetaR - thetaL)*(phiT - phiB);
 
  % Gauss-Kronrod (3,7) pair with degrees of precision 5 and 11.
  NODES = [ -0.9604912687080202, -0.7745966692414834, -0.4342437493468026, ...
             0, 0.4342437493468026,  0.7745966692414834,  0.9604912687080202 ];  
  NNODES = length(NODES);    
  ONEVEC = ones(2*NNODES,1);
  NARRAY = 0.25*[ NODES(:), NODES(:) ];
  WT3 = [ 0, 5/9, 0, 8/9, 0, 5/9, 0];
  WT7 = [ 0.1046562260264672, 0.2684880898683334, 0.4013974147759622, ...
          0.4509165386584744, 0.4013974147759622, 0.2684880898683334, ...
          0.1046562260264672 ];
      
  % Compute initial approximations on four subrectangles.  Initialize LIST
  % of information about subrectangles for which the approximations are
  % not sufficiently accurate.  NLIST is the number of subrectangles that
  % remain to be processed.  ERRBND is a bound on the error.
  [Qsub,esub] = tensor(thetaL,thetaR,phiB,phiT);
  Q = sum(Qsub);
  % Use an artificial value of TOL to force the program to refine.
  TOL = 100*eps*abs(Q);  
  ERR_OK = 0;
  ADJUST = 1;
  LIST = zeros(200,7);
  NLIST = 0;
  ERRBND = Inf; % Will be updated in Save2LIST()
  Save2LIST(Qsub,esub,thetaL,thetaR,phiB,phiT);
  if  NLIST == 0 || ERRBND <= TOL
      return;
  end
  
  while true
    % Get entries from LIST corresponding to the biggest (adjusted) error.
    [q,e,thetaL,thetaR,phiB,phiT] = NextEntry;

    % Approximate integral over four subrectangles.
    [Qsub,esub] = tensor(thetaL,thetaR,phiB,phiT);
    
    % Saved in LIST is "e", a conservative estimate of the error in the 
    % approximation "q" of the integral over a rectangle.  Newq = sum(Qsub)
    % is a much better approximation to the integral.  It is used here to 
    % estimate the error in "q" and thereby determine that the estimator is
    % conservative by a factor of "ADJUST".  This factor is applied to the 
    % estimate of the error in "Newq" to get a more realistic estimate.
    % This scheme loses the sign of the error, so a conservative local test 
    % is used to decide convergence.
    Newq = sum(Qsub);
    ADJUST = min(1,abs(q - Newq)/e);
    Q = Q + (Newq - q);
    TOL = max(ATOL,RTOL*abs(Q))/8;   
    TOL = max(TOL,100*eps*abs(Q));
    Save2LIST(Qsub,esub,thetaL,thetaR,phiB,phiT);
    
    % Test for convergence and failures.
    if  NLIST == 0 || ERRBND <= TOL
        break;
    elseif NLIST >= 2000
        if ERRBND > max(ATOL,max(100*eps,RTOL)*abs(Q))
            error('MATLAB:TwoD:maxsubs',...
                  'Maximum number of subintervals: Q does NOT pass error test.')
        else
            warning('MATLAB:TwoD:maxsubs',...
                    'Maximum number of subintervals: Q appears to pass error test.')
        end
        break;
    end

  end % while

  
  %==Nested functions======================================================
  
  function [Qsub,esub] = tensor(thetaL,thetaR,phiB,phiT)
  % Compute the integral with respect to theta from thetaL to thetaR of the 
  % integral with respect to phi from phiB to phiT of F in four blocks.
    Qsub = zeros(4,1);
    esub = zeros(4,1);
    
    dtheta = thetaR - thetaL;
    theta = bsxfun(@plus,dtheta*NARRAY,thetaL+dtheta*[0.25,0.75]);
    theta = reshape(theta,1,[]);                      
    if SINGULAR
        x = 0.5*(b + a) + 0.5*(b - a)*cos(theta);     
    else
        x = theta;
    end
    X = x(ONEVEC,:);   
    dphi = phiT - phiB;
    phi = bsxfun(@plus,dphi*NARRAY,phiB+dphi*[0.25,0.75]);
    phi = reshape(phi,1,[])';                    
    top = phiTvar(x); top = top(:)';             
    bottom = phiBvar(x); bottom = bottom(:)';    
    dydt = top - bottom;                         
    if SINGULAR
        t = 0.5 + 0.5*cos(phi);                   
    else
        t = phi;
    end
    Y = bottom(ONEVEC,:) + t * dydt;    
    if VECTORIZED
        Z = FUN(X,Y);
    else
        Z = arrayfun(FUN,X,Y);
    end
    if SINGULAR
        % Full matrix formed as outer product:
        temp = 0.25*(b - a)*sin(phi)*(dydt .* sin(theta));  
    else
        temp = dydt(ONEVEC,:);
    end
    Z = Z .* temp;
    
    % Tensor product Gauss 3 point formula:
    esub(1) = WT3 * (WT3 * Z(1:NNODES,1:NNODES))';
    esub(2) = WT3 * (WT3 * Z(1:NNODES,NNODES+1:end))';   
    esub(3) = WT3 * (WT3 * Z(NNODES+1:end,1:NNODES))';    
    esub(4) = WT3 * (WT3 * Z(NNODES+1:end,NNODES+1:end))';
    esub = (esub/4)*(dtheta/2)*(dphi/2);  
    
    % Tensor product Kronrod 7 point formula:    
    Qsub(1) = WT7 * (WT7 * Z(1:NNODES,1:NNODES))';
    Qsub(2) = WT7 * (WT7 * Z(1:NNODES,NNODES+1:end))';
    Qsub(3) = WT7 * (WT7 * Z(NNODES+1:end,1:NNODES))';
    Qsub(4) = WT7 * (WT7 * Z(NNODES+1:end,NNODES+1:end))';   
    Qsub = (Qsub/4)*(dtheta/2)*(dphi/2);  
    
    esub = abs(esub - Qsub);

  end % tensor

%--------------------------------------------------------------------------

  function [q,e,L,R,B,T] = NextEntry
    % Take entries from LIST corresponding to the biggest adjusted error:
    [~,index] = max(abs(LIST(1:NLIST,7)));
    temp = LIST(index,:);
    q = temp(1); e = temp(2); L = temp(3); 
    R = temp(4); B = temp(5); T = temp(6); 
    LIST(index,:) = [];
    NLIST = NLIST - 1;
  end % NextEntry

%--------------------------------------------------------------------------

  function Save2LIST(Qsub,esub,thetaL,thetaR,phiB,phiT)
    % Save to the list information about subrectangles for which the
    % integral is not sufficiently accurate. NLIST is the number of
    % subrectangles to be processed.
    dtheta = thetaR - thetaL;
    dphi = phiT - phiB;
    localtol = TOL*(dtheta/2)*(dphi/2)/AREA;
    localtol = max(localtol,100*eps*abs(sum(Qsub)));

    adjerr = ADJUST*esub;
    if NLIST+4 > size(LIST,1)
        LIST = [LIST; zeros(100,7)];
    end
    
    if adjerr(1) > localtol 
      NLIST = NLIST + 1;
      LIST(NLIST,:) = [Qsub(1),esub(1),thetaL,thetaL+dtheta/2,...
                                       phiB,phiB+dphi/2,adjerr(1)];
    else
      ERR_OK = ERR_OK + adjerr(1);
    end
    if adjerr(2) > localtol 
      NLIST = NLIST + 1;
      LIST(NLIST,:) = [Qsub(2),esub(2),thetaL+dtheta/2,thetaR,...
                                       phiB,phiB+dphi/2,adjerr(2)];
    else
      ERR_OK = ERR_OK + adjerr(2);      
    end
    if adjerr(3) > localtol
      NLIST = NLIST + 1;
      LIST(NLIST,:) = [Qsub(3),esub(3),thetaL,thetaL+dtheta/2,...
                                       phiB+dphi/2,phiT,adjerr(3)];
    else
      ERR_OK = ERR_OK + adjerr(3);     
    end
    if adjerr(4) > localtol
      NLIST = NLIST + 1;
      LIST(NLIST,:) = [Qsub(4),esub(4),thetaL+dtheta/2,thetaR,...
                                       phiB+dphi/2,phiT,adjerr(4)];
    else
      ERR_OK = ERR_OK + adjerr(4);
    end   
    
    ERRBND = ERR_OK + sum(LIST(:,7));
  end % Save2LIST

%--------------------------------------------------------------------------

  function parseOptions(varargin)
    % Parse optional input arguments.
    k = 1;
    while k < nargin
      propname = lower(varargin{k});
      k = k + 1;
      propvalue = varargin{k};
      switch propname
        case 'abstol'
          ATOL = propvalue;  
          if ~(isfloat(ATOL) && isscalar(ATOL) && ...
               isreal(ATOL) && ATOL >= 0)
            error('MATLAB:TwoD:invalidReltol','Invalid ABSTOL');
          end
        case 'reltol'
          RTOL = propvalue;
          if ~(isfloat(RTOL) && isscalar(RTOL) && ...
               isreal(RTOL) && RTOL >= 0 && RTOL < 1)
            error('MATLAB:TwoD:invalidReltol','Invalid RELTOL');
          end
        case 'sector'
          SECTOR = propvalue;
          if ~(isscalar(SECTOR) && islogical(SECTOR))
            error('MATLAB:TwoD:invalidSector','Invalid SECTOR');
          end          
        case 'singular'
          SINGULAR = propvalue;
          if ~(isscalar(SINGULAR) && islogical(SINGULAR))
            error('MATLAB:TwoD:invalidSingular','Invalid SINGULAR');
          end
        case 'vectorized'
          VECTORIZED = propvalue;
          if ~(isscalar(VECTORIZED) && islogical(VECTORIZED))
            error('MATLAB:TwoD:invalidVectorized','Invalid VECTORIZED');
          end           
        otherwise
          error('MATLAB:TwoD:unknownOption', 'Unrecognized option.');
      end
      k = k + 1;
    end
    if k == nargin
      error('MATLAB:TwoD:PVPairError', ...
            'Invalid option.  Expected option name followed by option value.');
    end
  end % parseOptions

%==========================================================================
  
end % TwoD
