function n = norm1(A, m)

%NORM1  One-norm
% function n = norm1(A, m)
% Input: m, necessary if A is function
%
% See also NORM, FNORM, INFNORM
%
% Revision date: July 14, 2009
% (C) Michiel Hochstenbach 2002-2009

if ~isnumeric(A)
  n = max(abs(A(ones(m,1), 1)));
else
  n = norm(A,1);
end
