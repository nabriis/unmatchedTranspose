function [index, flag] = sort_compl_conj(x, y, index)

% function [index, flag] = sort_compl_conj(x, y, index)
% In: x == conj(y(index)) except some possible places
%   where there occur complex conjugates
%   e.g. x(1) = y(1) = x(2)' = y(2)'
% Out: index is changed such that x == conj(y)
%
% See also SORT, SORT_TARGET
%
% Revision date: November 10, 2009
% (C) Michiel Hochstenbach 2002-2009

if nargin < 3 || isempty(index)
  index = 1:length(x);
end
flag = 0;

for k = 1:length(x)-1
  if abs(x(k)-y(k)') > 1e-2 * abs(x(k))
    % abs(x(k)-x(k+1)') > 1e-2 * abs(x(k)) || abs(y(k)-y(k+1)') > 1e-2 * abs(y(k))
    if abs(x(k)-y(k+1)') > 1e-2 * abs(x(k))
      disp('mismatch detected between computation of right and left eigenvalues in sort_compl_conj');
      flag = 1;
      return
    end
    [index(k), index(k+1)] = change_order(index(k), index(k+1));
    [y(k), y(k+1)] = change_order(y(k), y(k+1));
  end
end
