function b = strcmp0(s1, s2)

%STRCMP0  Compare strings with extra features
% function b = strcmp0(s1, s2)
%
% See also STRCMP
%
% Revision date: June 10, 2007
% (C) Michiel Hochstenbach 2014

if ~ischar(s1) || ~ischar(s2)
  b = false;
else
  b = strcmp(s1,s2);
end
