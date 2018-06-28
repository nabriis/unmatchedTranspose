function disp_info(str, m, n, target, nr, extraction, expansion, ...
  maxit, tol, fix, mindim, maxdim, maxit_inner, inner_tol, solver, ...
  thick, M1, Mtype, krylov, p)

%DISP_INFO  Display information for Jacobi-Davidson methods
% function disp_info(str, m, n, target, nr, extraction, expansion, ...
%   maxit, tol, fix, mindim, maxdim, maxit_inner, inner_tol, solver, ...
%   thick, M1, Mtype, krylov, p)
%
% Revision date: August 10, 2017
% (C) Michiel Hochstenbach 2017

if ~isempty(str)
  fprintf('\n*** %s ***\n\n', str);
end
for j = 1:length(m)
  fprintf('  Size of problem:          %d x %d\n', m(j), n(j));
end
fprintf('  Target:                   '); disp_target(target);
fprintf('  Number of values:         %d\n', nr);
if ~isa(extraction, 'function_handle')
  fprintf('  Extraction:               %s\n', extraction);
else
  fprintf('  Extraction:               special selection\n');
end
fprintf('  Expansion:                %s\n', expansion);
fprintf('  Max number of iterations: %d\n', maxit);
fprintf('  Absolute tolerance:       %g\n', tol);
% fprintf('  fix target until ||r|| <  '); yes_no(~ischar(fix), num2str(fix), fix);
if ischar(fix)
  fprintf('  Fix target until ||r|| <  %s\n', fix);
else
  fprintf('  Fix target until ||r|| <  %g\n', fix);
end
fprintf('  Dim search spaces:        min %d, max %d\n', mindim, maxdim);
fprintf('  Inner iteration:          max %d steps, tolerance %g\n', maxit_inner, inner_tol);
fprintf('  Inner solver:             %s\n', upper(solver));
fprintf('  Thick restart with:       %d vector(s)\n', thick);
fprintf('  Preconditioning:          '); yes_no(~isempty(M1));
if ~isempty(M1)
  fprintf('  Type of preconditioning:  %s\n', Mtype);
end
if nargin > 19 && ~isempty(p) && p > 1
  fprintf('  Block size:               %d\n', p);
end
fprintf('  Start with Krylov spaces: '); yes_no(krylov);
disp(' ');
