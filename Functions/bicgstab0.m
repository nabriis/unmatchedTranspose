function [x, flag, relres, iter, resvec] = bicgstab0(A, b, tol, maxit, M1, M2, x)

%BICGSTAB0   BiConjugate Gradients Stabilized Method
%
% $Revision: 1.19.4.8 $ $Date: 2011/05/17 02:33:02 $
% Adapted by Michiel Hochstenbach, March 20, 2014

% Assign default values to unspecified parameters
if nargin < 3 || isempty(tol)
  tol = 1e-6;
end
if nargin < 4 || isempty(maxit)
  maxit = 10;
end
if nargin < 5
  M1 = [];
end
if nargin < 6
  M2 = [];
end
if nargin < 7
  x = [];
end
n = length(b);

n2b = sqrt(b'*b);                  % Norm of rhs vector, b

% Set up for the method
flag = 1;
imin = 0;                          % Iteration at which xmin was computed
tolb = tol * n2b;                  % Relative tolerance
if isempty(x)
  r = b;
  normr = n2b;
  x = zeros(n,1);
else
  r = b - mv(A,x);
  normr = sqrt(r'*r);              % Norm of residual
end
xmin = x;                          % Iterate which has minimal residual so far
normr_act = normr;

rt = r;                            % Shadow residual
resvec = zeros(2*maxit+1,1);       % Preallocate vector for norm of residuals
resvec(1) = normr;                 % resvec(1) = norm(b-A*x0)
normrmin = normr;                  % Norm of residual from xmin
rho = 1;
omega = 1;
stag = 0;                          % stagnation of the method
alpha = [];                        % overshadow any functions named alpha
moresteps = 0;
maxmsteps = min([floor(n/50),10,n-maxit]);
maxstagsteps = 3;
% loop over maxit iterations (unless convergence or failure)

for ii = 1 : maxit
  rho1 = rho;
  rho = rt' * r;
  if (rho == 0.0) || isinf(rho)
    flag = 4;
    resvec = resvec(1:2*ii-1);
    break
  end
  if ii == 1
    p = r;
  else
    beta = (rho/rho1)*(alpha/omega);
    if (beta == 0) || ~isfinite(beta)
      flag = 4;
      break
    end
    p = r + beta * (p - omega * v);
  end
  if ~isempty(M1)
    ph = minvv(M1,p);
  else
    ph = p;
  end
  if ~isempty(M2)
    ph = minvv(M2,ph);
  end
  v = mv(A,ph);
  rtv = rt' * v;
  if (rtv == 0) || isinf(rtv)
    flag = 4;
    resvec = resvec(1:2*ii-1);
    break
  end
  alpha = rho / rtv;
  if isinf(alpha)
    flag = 4;
    resvec = resvec(1:2*ii-1);
    break
  end

  if abs(alpha)*sqrt(ph'*ph) < eps*sqrt(x'*x)
    stag = stag + 1;
  else
    stag = 0;
  end

  xhalf = x + alpha * ph;        % form the "half" iterate
  s = r - alpha * v;             % residual associated with xhalf
  normr = sqrt(s'*s);
  normr_act = normr;
  resvec(2*ii) = normr;

  % check for convergence
  if (normr <= tolb || stag >= maxstagsteps || moresteps)
    s = b - mv(A,xhalf);
    normr_act = sqrt(s'*s);
    resvec(2*ii) = normr_act;
    if normr_act <= tolb
      x = xhalf;
      flag = 0;
      iter = ii - 0.5;
      resvec = resvec(1:2*ii);
      break
    else
      if stag >= maxstagsteps && moresteps == 0
        stag = 0;
      end
      moresteps = moresteps + 1;
      if moresteps >= maxmsteps
        if 1 % ~warned
          warning(message('MATLAB:bicgstab:tooSmallTolerance'));
        end
        flag = 3;
        x = xhalf;
        resvec = resvec(1:2*ii);
        break;
      end
    end
  end

  if stag >= maxstagsteps
    flag = 3;
    resvec = resvec(1:2*ii);
    break
  end

  if normr_act < normrmin        % update minimal norm quantities
    normrmin = normr_act;
    xmin = xhalf;
    imin = ii - 0.5;
  end

  if ~isempty(M1)
    sh = minvv(M1,s);
    if ~all(isfinite(sh))
      flag = 2;
      resvec = resvec(1:2*ii);
      break
    end
  else
    sh = s;
  end
  if ~isempty(M2)
    sh = minvv(M2,sh);
  end
  t = mv(A,sh);
  tt = t' * t;
  if (tt == 0) || isinf(tt)
    flag = 4;
    resvec = resvec(1:2*ii);
    break
  end
  omega = (t' * s) / tt;
  if isinf(omega)
    flag = 4;
    resvec = resvec(1:2*ii);
    break
  end

  if abs(omega)*sqrt(sh'*sh) < eps*sqrt(xhalf'*xhalf)
    stag = stag + 1;
  else
    stag = 0;
  end

  x = xhalf + omega * sh;        % x = (x + alpha * ph) + omega * sh
  r = s - omega * t;
  normr = sqrt(r'*r);
  normr_act = normr;
  resvec(2*ii+1) = normr;

  % check for convergence        
  if (normr <= tolb || stag >= maxstagsteps || moresteps)
    r = b - mv(A,x);
    normr_act = sqrt(r'*r);
    resvec(2*ii+1) = normr_act;
    if normr_act <= tolb
      flag = 0;
      iter = ii;
      resvec = resvec(1:2*ii+1);
      break
    else
      if stag >= maxstagsteps && moresteps == 0
        stag = 0;
      end
      moresteps = moresteps + 1;
      if moresteps >= maxmsteps
        if 1 %~warned
          warning(message('MATLAB:bicgstab:tooSmallTolerance'));
        end
        flag = 3;
        resvec = resvec(1:2*ii+1);
        break;
      end
    end        
  end

  if normr_act < normrmin        % update minimal norm quantities
    normrmin = normr_act;
    xmin = x;
    imin = ii;
  end

  if stag >= maxstagsteps
    flag = 3;
    resvec = resvec(1:2*ii+1);
    break
  end
    
end                                % for ii = 1 : maxit

% Returned solution is first with minimal residual
if flag == 0
  relres = normr_act / n2b;
else
  r = b - mv(A,xmin);
  normr = sqrt(r'*r);
  if normr <= normr_act
    x = xmin;
    iter = imin;
    relres = normr / n2b;
  else
    iter = ii;
    relres = normr_act / n2b;
  end
end
