function [x, relres, k, hist, Ax] = gmres_fast(A, b, maxit, tol, M1, M2, Mtype, reortho)

%GMRES_FAST  GMRES version that is (much) faster than Matlab's version
% function [x, relres, k, hist, Ax] = gmres_fast(A, b, maxit, tol, M1, M2, Mtype, reortho)
% Assumes x0 = 0
% In:
%   maxit : maximum number of steps (no restarts)
%   Mtype : left | right | implicit | sym1 | sym2 (preconditioning)
%           left:     solve inv(M)Ax = inv(M)b
%           right:    solve A inv(M)y = b
%           implicit: solve A inv(M)y = b, and store inv(M)V
%           sym1:     solve inv(M1)A inv(M2)y = inv(M1)b
%           sym2:     solve inv(M1)A inv(M2)y = inv(M1)b, and store inv(M1)V
% Out:
%   relres: obtained residual reduction
%   k     : number of steps performed
%
% See also FOM
%
% Adapted from code by Gerard Sleijpen
%
% Revision date: January 15, 2015
% (C) Michiel Hochstenbach 2015

% Optional removal of warning
warning('off', 'MATLAB:rankDeficientMatrix');

% NB: if 1 step of gmres then x is a multiple of b
if nargin < 3 || isempty(maxit)
  maxit = 10;
end
if nargin < 4 || isempty(tol)
  tol = 1e-6;
end
if nargin < 5
  M1 = [];
end
if nargin < 6
  M2 = [];
end
if nargin < 7 || isempty(Mtype)
  Mtype = 'left';
end
if nargin < 8 || isempty(reortho)
  reortho = 1;
end

if maxit < 1 || tol >= 1
  relres = 1;
  x = b;
  return
end

left_prec_M   = ~isempty(M1) && strcmp(Mtype, 'left');
left_prec_M1  = ~isempty(M1) && (strcmp(Mtype, 'sym1')  || strcmp(Mtype, 'sym2'));
right_prec_M  = ~isempty(M1) && (strcmp(Mtype, 'right') || strcmp(Mtype, 'implicit'));
right_prec_M2 = left_prec_M1;
store         = ~isempty(M1) && (strcmp(Mtype, 'implicit') || strcmp(Mtype, 'sym2'));
extra_prec    = ~isempty(M1) && (strcmp(Mtype, 'right') || strcmp(Mtype, 'sym1'));
H = zeros(maxit+1, maxit);

if left_prec_M       % Precondition right-hand side
  y = minvv(M2, minvv(M1, b, 0), 0);
elseif left_prec_M1
  y = minvv(M1, b, 0);
end
rho0 = sqrt(b'*b);
v    = b * (1 / rho0);
Gamma = 1;

rho = 1;
W = [];              % W = inv(M) V
if tol == 0
  tol0 = Inf;
else
  tol0 = 1/(tol*tol); 
end
V = v;

for k = 1:maxit
  if right_prec_M || right_prec_M2
    if right_prec_M
      v = minvv(M2, minvv(M1, v, 0), 0);
    else
      v = minvv(M2, v, 0);
    end
    if store
      W(:,k) = v;    % Store W = inv(M) V, to avoid one extra action with inv(M)
    end
  end
  v = mv(A, v);
  if left_prec_M
    v = minvv(M2, minvv(M1, v, 0), 0);
  elseif left_prec_M1
    v = minvv(M1, v, 0);
  end
  if reortho
    [v, h] = rgs(v, V);
  else
    [v, h] =  gs(v, V);    
  end
  H(1:k+1,k) = h;
  gamma = h(k+1);
  
  if gamma == 0      % Lucky breakdown
    break
  end
  V(:,k+1) = v;
  gamma = -Gamma*h(1:k) / gamma; 
  Gamma(1,k+1) = gamma;
  rho = rho + gamma'*gamma;
  if nargout > 3
    hist(1,k) = 1 / sqrt(rho);
  end
  if rho >= tol0
    break
  end
end

if gamma == 0        % Lucky breakdown
  relres = 0;
  e1 = zeros(k,1); e1(1) = 1;
  c = H(1:k,1:k) \ e1;
  if store
    x = W(:,1:k)*c;
  else
    x = V(:,1:k)*c;
  end
  if nargout > 4
    Ax = V*(H*c);
  end
else                 % Solve in least square sense 
  relres = 1 / sqrt(rho);
  e1 = zeros(k+1,1); e1(1) = 1;
  c = H(1:k+1,1:k) \ e1;
  if store
    x = W*c;
  else
    x = V(:,1:k)*c;
  end
  if nargout > 4
    Ax = V*(H(1:k+1,1:k)*c);
  end
end

if extra_prec
  if right_prec_M
    x = minvv(M2, minvv(M1, x, 0), 0);
  else
    x = minvv(M2, x, 0);
  end
end
x = rho0 * x;
if nargout > 4
  Ax = rho0 * Ax;
end
