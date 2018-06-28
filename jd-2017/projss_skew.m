function x = projss_skew(x, U, V, VU)

%PROJSS_SKEW  Oblique projection onto a subspace
% function x = projss_skew(x, U, V, VU)
%
% Out: y = (I-UV')x
%
% See also PROJSS, PROJ, PROJ_SKEW, RGS_TWO
%
% Revision date: February 14, 2006
% (C) Michiel Hochstenbach 2014

if ~isempty(U) && ~isempty(x)
  if nargin < 4 || isempty(VU)
    % Assume V'U = I
    x = x - U*(V'*x);
  else
    x = x - U*(VU \ (V'*x));
  end
end
