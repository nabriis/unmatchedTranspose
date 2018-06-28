function disp_target(target)

%DISP_TARGET  Displays target
% function disp_target(target)
%
% target  = real or complex number
%           | Inf (largest)| -Inf (smallest)
%           | 'abs' (max magnitude)
%           | 'inf' (largest real part) | '-inf' (smallest real part)
%           | 'real' (smallest abs imag part) | 'imag' (smallest abs real part)
%           | 'imagpart' (largest imaginary part)
%
% Revision date: December 23, 2016
% (C) Michiel Hochstenbach 2016

if isnumeric(target)
  if isreal(target)
    fprintf('%g\n', target);
  else
    fprintf('%g + %gi\n', real(target), imag(target));
  end
elseif ischar(target) % Special target such as 'abs' etc
  if strcmp(target, 'abs')
    disp('maximal magnitude');
  elseif strcmp(target, 'inf')
    disp('largest, or maximal real part');
  elseif strcmp(target, '-inf')
    disp('smallest, or minimal real');
  elseif strcmp(target, '-infreal')
    disp('real with minimal real part');
  elseif strcmp(target, 'infreal')
    disp('real with maximal real part');
  elseif strcmp(target, 'real')
    disp('minimal absolute imaginary part ("most real")');
  elseif strcmp(target, 'imag')
    disp('minimal absolute real part ("most purely imaginary")');
  elseif strcmp(target, 'imagpart')
    disp('largest imaginary part');
  end
elseif isa(target, 'function_handle')
  disp('special function'); 
end
