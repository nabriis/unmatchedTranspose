function e = unv(j, n, a)

%UNV  Unit vector(s) in R^n or C^n
% function e = unv(j, n, a)
%
% Revision date: February 27, 2010
% (C) Michiel Hochstenbach 2016

if nargin < 3 || isempty(a)
  a = 1;
end

e = zeros(n,length(j));
for k = 1:length(j)
  e(j(k),k) = a;
end
