function [V, H] = krylov(A, b, k, full)

%KRYLOV  Construct orthonormal basis for Krylov subspace K_k(A,b)
% function [V, H] = krylov(A, b, k, full)
%
% Input:
%   A : n x n matrix
%   b : n vector
%   k : positive integer
% Output:
%   V : orthonormal basis of size k+1 (full) or k (~full)
%   H : Hessenberg size (k+1,k) (full) or (k,k) (~full)
%   AV_k = V_{k+1} H_{k+1,k},  H_{k,k} = V_k'AV_k
%
% See also RGS, KRYLOV_SYM, KRYLOV_EXPAND
%
% Revision date: June 9, 2016
% (C) Michiel Hochstenbach 2016


if isnumeric(A)
  n = size(A,2);
else % A is function
  if isempty(b)
    error('Need a starting vector if A is a function handle')
  else
    n = length(b);
  end
end
if nargin < 2 || isempty(b)
  b = randn1(n);
end
if nargin < 3 || isempty(k)
  k = 20;
end
if nargin < 4 || isempty(full)
  full = 1;
end

% Preallocation
if nargout > 1
  if full
    V = zeros(n,k+1);
    H = zeros(k+1,k);
  else
    V = zeros(n,k);
    H = zeros(k,k);
  end
end

normb = sqrt(b'*b);
V(:,1) = b * (1 / normb);
for j = 1:k
  if j < k || full
    if nargout > 1
      [V(:,j+1), H(1:j+1,j)] = rgs(mv(A,V(:,j),0), V(:,1:j));
    else
      V(:,j+1) = rgs(mv(A,V(:,j),0), V(:,1:j));
    end
  else
    if nargout > 1
      [~, h] = rgs(mv(A,V(:,j),0), V(:,1:j));
      H(1:j,j) = h(1:j);
    end
  end
end
