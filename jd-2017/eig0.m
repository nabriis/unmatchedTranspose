function [X, D, Y] = eig0(A, B, target)

%EIG0  EIG with some extras
% function [X, D, Y] = eig0(A, B, target)
%   - return diag(D) instead of D
%   - sort D if desired
%   - normalize eigenvectors for generalized eigenvalue problem
%   - take full(A), full(B) if A, B are sparse
%
% See also EIG
%
% Revision date: August 8, 2015
% (C) Michiel Hochstenbach 2015

if nargin < 2
  B = [];
end
if nargin < 3 || isempty(target)
  target = 'abs';
end

if issparse(A)
  A = full(A);
end
if issparse(B)
  B = full(B);
end
n = size(A,2);

if isempty(B)
  if nargout == 1
    D = eig(A);
  elseif nargout == 2
    [X,D] = eig(A);
    D = diag(D);
  else
    [X,D,Y] = eig(A);
    D = diag(D);
  end
else
  if nargout == 1
    D = eig(A,B);
  elseif nargout == 2
    [X,D] = eig(A,B);
    D = diag(D);
  else
    [X,D,Y] = eig(A,B);
    D = diag(D);
  end
  if nargout > 1
    for j = 1:n
      x = X(:,j);
      X(:,j) = x / sqrt(x'*x);
    end
  end
  if nargout > 2
    for j = 1:n
      x = Y(:,j);
      Y(:,j) = x / sqrt(x'*x);
    end
  end
end

[D, index] = sort_target(D, target);
if nargout > 1
  X = X(:,index);
end
if nargout > 2
  Y = Y(:,index);
end
 
if nargout < 2
  X = D;
end
