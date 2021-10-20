function [Flag]=my_equal(a,b)
%
% [Flag]=my_equal(a,b)
%
% This function checks if a and b are equal considering the tolerance
% setted in global_parameter file.
%
% Output of this function is Flag, that may assume following values:
%   1 if a = b, within the tolerance
%   0 in a <> b, within the tolerance
%


global tolerance

if abs(a-b)< tolerance
   Flag = 1;
else
   Flag = 0;
end

return