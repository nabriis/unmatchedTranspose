function s = select_target(x, target, rel)

%SELECT_TARGET  Select on target
% function s = select_target(x, target, rel)
% Give an integer number indicating closeness, but do not sort
% The higher, the closer
% target: complex number, '-inf', 'inf', 'x-axis', 'abs'
%
% See also SORT_TARGET, SCHUR0
%
% Revision date: December 11, 2007
% (C) Michiel Hochstenbach 2014

if nargin < 2 || isempty(target)
  target = 'abs';
end
if nargin < 3 || isempty(rel)
  rel = 0; % Sort on absolute, not on relative distance
end

[~, index] = sort_target(x, target, rel);
s(index) = length(x):-1:1;
