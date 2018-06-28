function opts = write_options(opts, varargin)

%WRITE_OPTIONS  Write options to a structure
%function opts = write_options(opts, varargin)
%
% See also READ_OPTIONS
%
% Revision date: May 28, 2007
% (C) Michiel Hochstenbach 2002-2007

for j = 2:nargin
  opts.(inputname(j)) = varargin{j-1};
% MIQ old code  opts = setfield(opts, inputname(j), varargin{j-1});
end