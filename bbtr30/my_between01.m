function [Flag]=my_between01(a)
%
% [Flag]=my_between01(a)
%
% This function checks if 0 <= a <= 1 considering the tolerance
% setted in global_parameter file.
%
% Output of this function is Flag, that may assume following values:
%   1 if 0 <= a <= 1, within the tolerance
%   0 otherwise, within the tolerance
%


global tolerance

if a > tolerance
   if (1-a) > tolerance
       Flag = 1;
   else
       Flag = 0;
   end
else
   Flag = 0;
end

return