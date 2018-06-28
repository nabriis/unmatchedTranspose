function [H, V, lambda, hist, nr_mv, flag] = krylov_schur(A, opts)

%KRYLOV_SCHUR  Krylov-Schur method
% function [H, V, lambda, hist, nr_mv, flag] = krylov_schur(A, opts)
%
% Opts can have the following fields:
%   A            square matrix, or function
%                  if function, v1 should be given, or opts.n should be size
%   nr           number of desired eigenpairs                        1   
%   v1           initial space (may be more-dim.)                    randn(n,1)
%   tol          tolerance of the outer iteration                    1e-6
%   absrel       absolute or relative tolerance outer iteration      'rel'
%                  relative tolerance: ||AU-UL|| < tol * ||A||_1
%   mindim       minimum dimension of subspaces                      10
%   maxdim       maximum dimension of subspaces                      30
%   maxit        maximum number of outer iterations                  100 
%   target       complex nr|'-inf'|'inf'|'infreal'|'-infreal'|'imag'|'abs'     'abs'
%   extraction   standard | harmonic                                 'standard'       
%   verbosity    output desired                                      0
%
% References:
%   A Krylov--Schur algorithm for large eigenproblems
%   GW Stewart - SIAM Journal on Matrix Analysis and Applications, 2002
%
%   Addendum to "A Krylov--Schur Algorithm for Large Eigenproblems"
%   GW Stewart - SIAM journal on matrix analysis and applications, 2002
%
% Revision date: January 19, 2015
% (C) Michiel Hochstenbach 2015

if nargin < 2
  opts = [];
end
lambda = [];
if isfield(opts, 'n'),           n           = opts.n;           else n           =  size(A,2); end
if isfield(opts, 'target'),      target      = opts.target;      else target      =      'abs'; end
if isfield(opts, 'extraction'),  extraction  = opts.extraction;  else extraction  = 'standard'; end
if isfield(opts, 'v1'),          v1          = opts.v1;          else v1          = randn(n,1); end
if isfield(opts, 'nr'),          nr          = opts.nr;          else nr          =          1; end
if isfield(opts, 'mindim'),      mindim      = opts.mindim;      else mindim      =         10; end
if isfield(opts, 'maxdim'),      maxdim      = opts.maxdim;      else maxdim      =         20; end
if isfield(opts, 'maxit'),       maxit       = opts.maxit;       else maxit       =       1000; end
if isfield(opts, 'tol'),         tol         = opts.tol;         else tol         =       1e-6; end
if isfield(opts, 'absrel'),      absrel      = opts.absrel;      else absrel      =      'rel'; end
if isfield(opts, 'verbosity'),   verbosity   = opts.verbosity;   else verbosity   =          0; end
if strcmp(absrel, 'rel') && isnumeric(A)
  tol = tol * norm(A,1);
end

if mindim < nr
  mindim = nr;
end
if maxdim < 2*mindim
  maxdim = 2*mindim;
end
if strcmp(extraction, 'harmonic') && ~isnumeric(target)
  histor('Target needs to be complex number for harmonic extraction')
end

flag = 0;

if verbosity
  fprintf('\n*** Krylov-Schur ***\n\n');
  fprintf('  Size of problem:               %d\n', n);
  fprintf('  Target:                        '); disp_target(target);
  fprintf('  Number of eigenvalues:         %d\n', nr);
  fprintf('  Extraction:                    %s\n', extraction);
  fprintf('  Max number of iterations:      %d\n', maxit);
  fprintf('  Tolerance for ||H(k+1,1:nr)||: %g\n', tol);
  fprintf('  Dim search spaces:             min %d, max %d\n', mindim, maxdim);
  disp(' ');
end
if verbosity > 1
  fprintf(' Iter ||H(k+1,1:nr)||    theta\n');
  fprintf('-----------------------------------\n')
end

[V, H] = krylov(A, v1, mindim);
mindim1 = mindim;

for k = 1:maxit
  [V, H] = krylov_expand(A, V, H, maxdim-mindim1);
  if nargout > 4 || verbosity
    if k > 1
      nr_mv(k) = nr_mv(k-1) + maxdim-mindim1;
    else
      nr_mv = maxdim;
    end
  end
  if strcmp(extraction, 'standard')
    % Ritz values
    [Q,S] = schur0(H(1:maxdim,1:maxdim), target);
    % shifts = element(diag(S), mindim+1:maxdim);
    if S(mindim+1, mindim)
      mindim1 = mindim+1;
    else
      mindim1 = mindim;
    end
    V = [element(V(:,1:maxdim)*Q, [], 1:mindim1) V(:,maxdim+1)];
    H = [S(1:mindim1,1:mindim1); H(maxdim+1,:)*Q(:,1:mindim1)];
  else
    % Harmonic Ritz values
    g = (H(1:maxdim,1:maxdim)-target*eye(maxdim))' \ (H(maxdim+1,:)');
    vhat = V(:,maxdim+1)-V(:,1:maxdim)*g;
    [Q,S] = schur0(H(1:maxdim,1:maxdim) + g*H(maxdim+1,:), target);
    V = element(V(:,1:maxdim)*Q, 1:n, 1:mindim);
    [V(:,mindim+1), h] = rgs(vhat, V);
    H = [S(1:mindim,1:mindim); element(H(maxdim+1,:)*Q, 1:mindim)];
    H(1:mindim,1:mindim) = H(1:mindim,1:mindim)+h(1:mindim)*H(mindim+1,:);
    H(mindim+1,:) = H(mindim+1,:)*h(mindim+1);
  end
  hist(k) = norm(H(mindim1+1, 1:nr));
  if (verbosity == 2) || (verbosity > 2 && ~mod(k, verbosity))
    fprintf('%4d  %6.2e  | ', k, hist(k));
    if H(nr+1,nr)
      lambda = ordeig(H(1:nr+1,1:nr+1));
    else
      lambda = ordeig(H(1:nr,1:nr));
    end
    fprintf('%6.3g + i* %.3g ', real(lambda(1:nr)), imag(lambda(1:nr)));
    fprintf('\n');
  end
  if hist(k) < tol
    if H(nr+1,nr) == 0
      lambda = ordeig(H(1:nr,1:nr));
    else
      lambda = ordeig(H(1:nr+1,1:nr+1));
    end
    if verbosity
      fprintf('**** Found after %d restarts, %d mvs, with norm residual = %8.3g\n', k, nr_mv(k), hist(k));
      for i = 1:nr
        fprintf('  lambda = %9g %+9gi\n', real(lambda(i)), imag(lambda(i)));
      end
    end
    return
  end
end

% No convergence
flag = 1;
warning('Reached maxit')
