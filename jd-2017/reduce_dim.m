function varargout = reduce_dim(X, varargin)

%REDUCE_DIM  Tool to restart with part of subspace
% varargout = reduce_dim(X, varargin)
%
% See also REDUCE_DIM_LEFT, REDUCE_DIM_TWO
%
% Revision date: January 22, 2014
% (C) Michiel Hochstenbach 2014

varargout = cell(1, nargin-1);
for j = 1:nargin-1
  if ~isempty(varargin{j})
    varargout{j} = varargin{j}*X;
  end
end
