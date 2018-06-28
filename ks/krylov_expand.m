function [V, H] = krylov_expand(A, V, H, k)

%KRYLOV_EXPAND  Expand orthonormal basis for Krylov subspace K_m(A,b) to K_{k+m}(A,b)
% function [V, H] = krylov_expand(A, V, H, k)
%
% Input:
%   V  n x (m+1) orthonormal columns, nonempty (if empty use krylov)
%   H  (m+1) x m Hessenberg
% Output:
%   V  n x (k+m+1) orthonormal columns, nonempty (if empty use krylov)
%   H  (k+m+1) x m Hessenberg
%
% See also KRYLOV
%
% Revision date: June 9, 2016
% (C) Michiel Hochstenbach 2016

if nargin < 4 || isempty(k)
  k = 10;
end

% Preallocation
m = size(V,2);
V = [V zeros(size(V,1), k)];

for j = m:(k+m-1)
  [V(:,j+1), H(1:j+1,j)] = rgs(mv(A, V(:,j)), V(:,1:j));
end
