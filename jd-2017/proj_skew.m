function y = proj_skew(x, a, b, ba)

%PROJ_SKEW  One-dimensional oblique projection
% function y = proj_skew(x, a, b, ba)
%         ab'
% y = (I- ---) x
%         b'a
%
% See also PROJ
%
% Revision date: March 5, 2014
% (C) Michiel Hochstenbach 2014

if nargin < 4 || isempty(ba)
  % Assume b'a = 1
  y = x - a * (b'*x);
else
  y = x - a * ((b'*x) / ba);
end
