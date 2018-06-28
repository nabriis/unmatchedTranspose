function f = getfield0(s, fieldname, default)

%GETFIELD0  Get field from a structure with extra feature
% function f = getfield0(s, fieldname, default)
%
% See also GETFIELD
%
% Revision date: January 22, 2014
% (C) Michiel Hochstenbach 2014

if nargin < 3
  default = [];
end

try
  f = s.(fieldname);
catch
  f = default;
end
