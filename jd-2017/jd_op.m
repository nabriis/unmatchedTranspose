function y = jd_op(x, A, B, theta, u1, v1, v1u1, u2, v2, v2u2, U1, V1, U2, V2, transp_flag)

%JD_OP  Perform action of JD operator in correction equation
% function y = jd_op(x, A, B, theta, u1, v1, v1u1, u2, v2, v2u2, U1, V1, U2, V2, transp_flag)
% In:  V1'*U1 = I
%      V2'*U2 = I
% Out: y = (I-u1*v1')(I-U1*V1')(A-theta*B)(I-U2*V2')(I-u2*v2') x
% If transp_flag then transpose is performed
%
% See also JD, JD_OP_SYM, JD_OP_PREC
%
% Revision date: March 9, 2014
% (C) Michiel Hochstenbach 2016

if nargin < 15 || isempty(transp_flag)
  transp_flag = 0;
end

if ~transp_flag
  if isempty(B)
    if theta ~= 0
      x = mv(A,x) - theta*x;
    else
      x = mv(A,x);
    end
  else
    if theta ~= 0
      x = mv(A,x) - mv(B,x)*theta;
    else
      x = mv(A,x);
    end
  end
  x = projss_skew(x, U1, V1);
  y = proj_skew(x, u1, v1, v1u1);
else
  % Transpose
  if isempty(B)
    if theta ~= 0
      x = mv(A,x,1) - theta'*x;
    else
      x = mv(A,x,1);
    end
  else
    if theta ~= 0
      x = mv(A,x,1) - mv(B,x,1)*theta';
    else
      x = mv(A,x,1);
    end
  end
  x = projss_skew(x, V2, U2);
  y = proj_skew(x, v2, u2, v2u2');
end
