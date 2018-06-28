function [y, index] = sort_target(x, target, rel)

%SORT_TARGET  Sort on target, y(1) = x(index(1)) is closest to target
% function [y, index] = sort_target(x, target, rel)
% target: complex number, '-inf', 'inf', '-infreal', 'infreal',
%   'real', 'imag', 'abs'
%
% See also SORT, SORT_COMPL_CONJ, SELECT_TARGET
%
% Revision date: December 17, 2009
% (C) Michiel Hochstenbach 2013

if nargin < 2 || isempty(target)
  target = 'abs';
end
if nargin < 3 || isempty(rel)
  rel = 0; % sort on absolute, not on relative distance
end

if isnumeric(target)
  if ~rel
    [~, index] = sort(abs(x-target));
  else
    [~, index] = sort(abs((x-target)./x));
  end
elseif isa(target, 'function_handle')
  % maximize function
  [~, index] = sort(target(x), 'descend');
else % string
  if strcmp(target, 'abs')
    [~, index] = sort(-abs(x));
  elseif strcmp(target, 'inf')
    [~, index] = sort(-real(x));
  elseif strcmp(target, 'infreal')
    if size(x,1) == 1
      x = x.';
    end
    [~, index] = sort(-real(x));
    xx = x(index);
    i1 = find(abs(imag(xx)) <= 1e-12);
    i2 = find(abs(imag(xx)) > 1e-12);
    index = index([i1; i2]);
  elseif strcmp(target, '-inf')
    [~, index] = sort(real(x));
  elseif strcmp(target, '-infreal')
    if size(x,1) == 1
      x = x.';
    end
    [~, index] = sort(real(x));
    xx = x(index);
    i1 = find(abs(imag(xx)) <= 1e-12);
    i2 = find(abs(imag(xx)) > 1e-12);
    index = index([i1; i2]);
  elseif strcmp(target, 'real')     % most real = smallest imaginary part
    [~, index] = sort(abs(imag(x)));
  elseif strcmp(target, 'imag')     % most imaginary = smallest real part
    [~, index] = sort(abs(real(x)));
  elseif strcmp(target, 'imagpart') % largest imaginary part
    [~, index] = sort(-imag(x));
  else
    error('Unknown target in sort_target');
  end
end
y = x(index);
