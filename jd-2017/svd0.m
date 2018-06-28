function [U, S, V] = svd0(A, fullsize, target, tall)

%SVD0  SVD with some extras
% function [U, S, V] = svd0(A, fullsize, target, tall)
% - Return diag(S) instead of S
% - Sort S if desired
%
% See also SVD
%
% Revision date: March 11, 2014
% (C) Michiel Hochstenbach 2014

if nargin < 2 || isempty(fullsize)
  fullsize = 1;
end
if nargin < 3
  target = [];
end
if nargin < 4 || isempty(tall)
  tall = 0;
end

if ~tall
  if fullsize
    if issparse(A)
      [U,S,V] = svd(full(A));
    else
      [U,S,V] = svd(A);
    end
  else
    if issparse(A)
      [U,S,V] = svd(full(A),0);
    else
      [U,S,V] = svd(A,0);
    end
  end
else
  % First determine QR decomposition: was faster for tall matrices in the past
  %   not any longer
  [Q,R] = qr(A,0);
  [U,S,V] = svd(R,0);
  if nargout > 2
    U = Q*U;
  end
end

if min(size(S)) == 1
  S = S(1);
else
  S = diag(S);
end

if ~isempty(target) && ~ischar(target)
  if target == 0   % Smallest singular value
    l = length(S);
    S = S(l:-1:1);
    U(:,1:l) = U(:,l:-1:1);
    V(:,1:l) = V(:,l:-1:1);
  else             % Other singular value: should be rare
    [S, index] = sort_target(S, target);
    U = U(:,index);
    V = V(:,index);
  end
end

if nargout < 2
  U = S;
elseif nargout == 2
  U = S;
  S = V;
end
