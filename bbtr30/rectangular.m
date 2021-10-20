function [Flag] = rectangular (Descriptors, P)
%
% [Flag] = rectangular (Descriptors,P)
%
% -----------------------------------------------------------------------
% METHOD FOR BTR
% -----------------------------------------------------------------------
%
% "rectangular" is a method for the triangler, and, like all methods,
% describes a specific region in which a specific refining must be applied.
% 
% "rectangular" defines a rectangular region described by four values:
% x_max, x_min, y_max, y_min. These are stored in the descriptors.
%
% Like every method "rectangular" must evaluate a point whose coordinates
% are P = [xP,yP] and find if the point belongs to the specified region or not.
% If P belongs to the region Flag must be put = true, otherwise must be put = false
%
% Inputs & Output of this function are:
%
% Descriptors: it specify the limits of the ractangle and it's in the form:
%    Descriptors = [x_max, x_min, y_max, y_min]
%
% P: the generic point (P = [xP,yP]) that must control if belongs or not no
% the region
%
% Flag: the boolean that is true if P belongs to the region

if P(1) > Descriptors (2)
    if P(1) < Descriptors (1)
        if P(2) > Descriptors (4)
            if P(2) < Descriptors (3)
                Flag = true;
            else
                Flag = false;
            end
        else
            Flag = false;
        end
    else
        Flag = false;
    end
else
    Flag = false;
end

    
return