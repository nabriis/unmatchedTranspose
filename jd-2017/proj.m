function y = proj(x, a)

%PROJ  One-dimensional orthogonal projection
% function y = proj(x, a)
% In:  a'*a = 1
% Out: y = (I-aa')x
%
% See also PROJSS, PROJ_SKEW
%
% Revision date: May 27, 2004
% (C) Michiel Hochstenbach 2002-2008

y = x-a*(a'*x);
