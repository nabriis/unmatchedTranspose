function [Q, S, Z, T] = qz0(A, B, target)

%QZ0  Schur or QZ factorization with some extra features
% function [Q, S, Z, T] = qz0(A, B, target)
% 
% Z*A*Q = S,  Z*B*Q = T
%
% Revision date: October 5, 2010
% (C) Michiel Hochstenbach 2010

if nargin < 3 || isempty(target)
  target = 'abs';
end
if nargin < 2 || isempty(B)
  if issparse(A)
    [Q, S] = schur0(full(A), target);
  else
    [Q, S] = schur0(A, target);
  end
else
  if issparse(A) || issparse(B)
    [S, T, Z, Q] = qz(full(A), full(B));
  else
    [S, T, Z, Q] = qz(A, B);
  end
  [S, T, Z, Q] = ordqz(S, T, Z, Q, select_target(diag(S)./diag(T), target));
  if nargout < 3
    S = diag(S)./diag(T);
  end
end
