function [theta, D, x, X] = eig_real(theta, D, x, X)

%EIG_REAL  Ensure real eigenvalues and eigenvectors
% function [theta, D, x, X] = eig_real(theta, D, x, X)
%
% Revision date: April 18, 2009
% (C) Michiel Hochstenbach 2002-2009

i1 = 1;
while i1 <= length(D)
  if isreal(D(i1))
    if ~isreal(X(:,i1))
      X(:,i1) = real(sign(X(1,i1))'*X(:,i1));
    end
  else
    X(:,i1:i1+1) = orth([real(X(:,i1)) imag(X(:,i1))]);
    % X(:,i1:i1+1) = [real(X(:,i1)) imag(X(:,i1))];
    i1 = i1+1;
  end
  i1 = i1+1;
end
x = X(:,1);
if ~isreal(theta)
  theta = real(theta);
end
