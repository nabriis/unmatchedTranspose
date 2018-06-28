function varargout = reduce_dim_two(X, Y, varargin)

%REDUCE_DIM_TWO  Tool to restart with part of subspace in two-sided process
% varargout = reduce_dim_two(X, Y, varargin)
%
% See also REDUCE_DIM, REDUCE_DIM_LEFT
%
% Revision date: January 22, 2014
% (C) Michiel Hochstenbach 2014

varargout = cell(1, nargin-2);
for j = 1:nargin-2
  if ~isempty(varargin{j})
    varargout{j} = Y'*varargin{j}*X;
  end
end
