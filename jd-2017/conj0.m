function x = transpose0(x)

%CONJ0  Complex conjugate with extra feature
% function x = transpose0(x)
%
% See also CONJ
%
% Revision date: November 10, 2009
% (C) Michiel Hochstenbach 2002-2009

if isnumeric(x)
  % do nothing for strings such as 'abs'
  x = conj(x);
end
