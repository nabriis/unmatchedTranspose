function [Q,S] = schur0(A, target, real)

%SCHUR0  Schur form with some extra features
% function [Q,S] = schur0(A, target, real)
%
% When called with <= 1 output argument, S is given
%
% Revision date: December 29, 2017
% (C) Michiel Hochstenbach 2017

if nargin < 2 || isempty(target)
  target = 'abs';
end
if nargin < 3 || isempty(real)
  real = 0;
end

if issparse(A)
  A = full(A);
end

if ~real
  [Q,S] = schur(A, 'complex');
  [~, index] = sort_target(diag(S), target);
  index2(index) = length(index):-1:1;
  [Q,S] = ordschur(Q, S, index2);
else
  [Q,S] = schur(A);
  [~, index] = sort_target(ordeig(S), target);
  % Treat complex conjugates together
  index2(index) = length(x):-1:1;
  for i = 1:size(S)-1
    if S(i+1,i)
      index2(i) = index2(i+1);
    end
  end
  [Q,S] = ordschur(Q, S, index2);
end
if nargout < 2
  Q = S;
end
