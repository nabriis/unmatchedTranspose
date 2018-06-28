function [u, normu] = normalize_vec(u)

%NORMALIZE_VEC  Normalize vector(s)
% function [u, normu] = normalize_vec(u)
% perform u = u / norm(u) on every column of u
%
% Revision date: February 9, 2006
% (C) Michiel Hochstenbach 2002-2006

for j = 1:size(u,2)
  normu(j) = norm(u(:,j));
  u(:,j)   = u(:,j) / normu(j);
end
