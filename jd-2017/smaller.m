function b = smaller(x, y)

%SMALLER  Check if x<y, including case 'inf's
% function b = smaller(x, y)
% In:  x : real number
%      y : real number | 'inf'
% Out: b = true iff y = 'inf' | x < y
%
% Revision date: December 30, 2003
% (C) Michiel Hochstenbach 2002-2004

if ischar(y)
  b = true;
else
  b = (x < y);
end
