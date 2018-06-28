function [v, h] = rgs(v, U)

%RGS  Repeated classical Gram-Schmidt
% function [v, h] = rgs(v, U)
%
% The projection is performed at least once, possibly twice ("twice is enough")
%
% In:  U'U = I
% Out: v = (I-UU')v and normalization computed in a stable way
%      h = coefficients u_j'*v, for instance for use in Krylov methods
%
% See also PROJSS
%
% Revision date: June 9, 2016
% (C) Michiel Hochstenbach 2016

% Defining constant for repetition
c = 0.70; % ~= 1/sqrt(2);

norm_input = sqrt(v'*v);
k = size(U,2);

if k == 0
  h = sqrt(v'*v);
  v = v / h;
  return
end

for i = 1:2
  if i == 2 && (norm_output > c * norm_input)
    return
  end
  a = U'*v;
  v = v-U*a;
  norm_output = sqrt(v'*v);   % Faster than norm
  %keyboard
  v = v * (1 / norm_output);  % Faster than division
  if nargout > 1
    if i == 1
      h = [a; norm_output];
    else
      h(1:k,1) = h(1:k,1) + h(k+1,1) * a;
      h(k+1,1) = h(k+1,1) * norm_output;
    end
  end
end
