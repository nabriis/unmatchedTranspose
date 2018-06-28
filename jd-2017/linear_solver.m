function [x, relres, mvs, precs] = linear_solver(A, b, solver, maxit, tol, M1, M2, Mtype)

%LINEAR_SOLVER  Call linear solver with additional checks
% function [x, relres, mvs, precs] = linear_solver(A, b, solver, maxit, tol, M1, M2, Mtype)
%   
%   solver: gmres | minres | bicgstab | qmr | cg | bicg | idrs | fom | symmlq | tfqmr
%     Solver is not restarted, maxit = max number of steps
%   Mtype : left | right | implicit | sym1 | sym2 (preconditioning)
%           left:     solve inv(M)Ax = inv(M)b
%           right:    solve A inv(M)y = b
%           implicit: solve A inv(M)y = b, and store inv(M)V
%                       (saves 1 action with preconditioner at the cost of more storage)
%           sym1:     solve inv(M1)A inv(M2)y = inv(M1)b
%           sym2:     solve inv(M1)A inv(M2)y = inv(M1)b, and store inv(M1)V
%   Assumes x0 = 0
%
% See also FOM, GMRES_FAST, MINRES0, SYMMLQ0
%
% Revision date: December 2, 2014
% (C) Michiel Hochstenbach 2014

if nargin < 3 || isempty(solver)
  solver = 'gmres';
end
if nargin < 4 || isempty(maxit)
  maxit = 10;
end
if nargin < 5 || isempty(tol)
  tol = 0;
end
if nargin < 6
  M1 = [];
end
if nargin < 7
  M2 = [];
end
if nargin < 8
  Mtype = [];
end

precs = 0;
mvs = 0;
if maxit == 0
  if isempty(M1)
    x = b;
%     warning('0 steps of linear solver without preconditioning');
  else % Left-precondition right-hand side
    x = prec(M1, M2, b);
    precs = 1;
  end
  relres = NaN;
elseif maxit == size(A,2)
  % Solve with backslash
  x = minvv(A,b);
  relres = 0;
  mvs = 0;
else
  % Own routines, (slightly) faster than Matlab version
  if strcmp(solver, 'gmres')
    [x, relres, mvs] = gmres_fast(A, b, maxit, tol, M1, M2, Mtype, 1);
  elseif strcmp(solver, 'minres')
%     warning('off', 'MATLAB:minres:tooSmallTolerance') % needed for Matlab's minres
    [x, ~, relres, mvs] = minres0(A, b, tol, maxit, M1, M2);
%     [x, relres, mvs] = minres_ger(A, b, maxit, tol); % does not have preconditioning at the moment
  elseif strcmp(solver, 'bicgstab')
    [x, ~, relres, mvs] = bicgstab(A, b, tol, maxit, M1, M2);
  elseif strcmp(solver, 'qmr')
    [x, ~, relres, mvs] = qmr0(A, b, tol, maxit, M1, M2);
    mvs = 2*mvs;
  elseif strcmp(solver, 'tfqmr')
    [x, ~, relres, mvs] = tfqmr0(A, b, tol, maxit, M1, M2);
    mvs = 2*mvs;
  elseif strcmp(solver, 'cg')
    [x, ~, relres, mvs] = cg(A, b, tol, maxit, M1, M2);
  elseif strcmp(solver, 'bicg')
    [x, ~, relres, mvs] = bicg0(A, b, tol, maxit, M1, M2);
%     [x, relres, mvs] = bicg_ger(A, b, tol, maxit, M1, M2);
    mvs = 2*mvs;
  elseif strcmp(solver, 'cr')
    [x, relres, mvs] = cr(A, b, tol, maxit, M1, M2);
  elseif strcmp(solver, 'lsqr')
    [x, ~, relres, mvs] = lsqr0(A, b, tol, maxit, M1, M2);

  elseif strcmp(solver, 'fom')
    [x, relres, mvs] = fom(A, b, maxit, tol, M1, M2, Mtype);

  % Matlab routines    
  elseif strcmp(solver, 'gmres_ml')
    [x, ~, relres, mvs] = gmres(A, b, maxit, tol, 1, M1, M2);

  elseif strcmp(solver, 'idrs')
    [x, ~, relres, mvs] = idrs(A, b, [], tol, maxit, M1, M2, []);

  elseif strcmp(solver, 'minres_jdsvd')
    [x, relres, mvs] = minres_jdsvd(A, b, maxit, tol, M1, M2, Mtype);

  elseif strcmp(solver, 'symmlq')
    [x, relres, mvs] = symmlq0(A, b, maxit, tol, M1, M2, Mtype);
  else
    fprintf('Unknown solver\n')
  end
  
  if relres >= 1
    x = b;
  end
  if nargout > 3
    if ~isempty(M1)
      if strcmp(Mtype, 'implicit')
        precs = mvs;
      else
        precs = mvs+1;
      end
    end
  end
end
