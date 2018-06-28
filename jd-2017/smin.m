function [s, y, x] = smin(A, k, sparse)

%SMIN  Smallest singular triple(s)
% function [s, y, x] = smin(A, k, sparse)
% Minimal singular value (or triple) of A
%   s = sigma_min
%   y = smallest right singular vector
%   x = smallest left  singular vector
%
% See also SMAX, LMIN, LMAX
%
% Revision date: January 23, 2011
% (C) Michiel Hochstenbach 2011

if nargin < 2
  k = 1;
end
if nargin < 3
  sparse = 0;
end

if ~sparse && issparse(A)
  A = full(A);
end
if ~issparse(A)
  if nargout > 1
    [U,S,V] = svd(A,0);
    S = diag(S);
  else
    S = svd(A);
  end
  [S, index] = sort(S);
  if nargout > 1
    U = U(:, index);
    V = V(:, index);
  end
else
  [U,S,V] = svds(A,k+1,0);
  S = diag(S);
end
if min(size(A)) > 1 && abs(S(k) - S(k+1)) < 1e-15
  fprintf('Warning: sigma_min close to multiple:  %8.2e  %8.2e\n', S(k), S(k+1));
%     S(k)-S(k+1), relerr(S(k),S(k+1)));
end;
s = S(1:k);
if nargout > 1
	x = U(:,1:k);
	y = V(:,1:k);
end
