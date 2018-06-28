function printnr(x, form, re, newline)

%PRINTNR  Display (complex) number
% function printnr(x, form, re, newline)
%
% Revision date: December 27, 2007
% (C) Michiel Hochstenbach 2002-2007

if nargin < 2 || isempty(form)
  form = '7.4g';
end
if nargin < 3 || isempty(re)
  re = 0;
end
if nargin < 4 || isempty(newline)
  newline = 0;
end

for j = 1:length(x)
  if re
    s = remove_zeros_in_exp(sprintf(strcat(['%' form]), real(x(j))));
  else
    if imag(x(j)) > 0
      s = strcat([remove_zeros_in_exp(sprintf(strcat(['%' form]), real(x(j)))) ...
      remove_zeros_in_exp(sprintf(strcat([' + %' num2str(str2num(form(1))-1) form(2:end) 'i']), ...
      imag(x(j))))]);
    else
      s = strcat([remove_zeros_in_exp(sprintf(strcat(['%' form]), real(x(j)))) ...
      remove_zeros_in_exp(sprintf(strcat([' - %' num2str(str2num(form(1))-1) form(2:end) 'i']), ...
      abs(imag(x(j)))))]);
    end
  end
  fprintf('%s ', s);
  if newline
    fprintf('\n');
  else
    if j < length(x)
      fprintf(' & ');
    end
  end
end
