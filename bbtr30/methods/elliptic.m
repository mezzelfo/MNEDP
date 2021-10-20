function [Flag] = elliptic (Descriptors, P)
%
% [Flag] = elliptic (Descriptors,P)
%
% -----------------------------------------------------------------------
% METHOD FOR BTR
% -----------------------------------------------------------------------
%
% "elliptic" is a method for the triangler, and, like all methods,
% describes a specific region in which a specific refining must be applied.
% 
% "elliptic" defines an elliptic region described by four values:
% xC, yC, a, b. These are stored in the descriptors.
%
% Like every method "elliptic" must evaluate a point whose coordinates
% are P = [xP,yP] and find if the point belongs to the specified region or not.
% If P belongs to the region Flag must be put = true, otherwise must be put = false
%
% Inputs & Output of this function are:
%
% Descriptors: it specify the limits of the ractangle and it's in the form:
%    Descriptors = [xC, yC, a, b]
%       c = [xC yC] is the center of the ellipse
%       a is the semiaxis in x-direction
%       b is the semiaxis in y-direction
%
%    ( The ellipse equation used is (xP-xC)^2/a^2 + (yP-yC)^2/b^2 = 1 )
%
% P: the generic point (P = [xP,yP]) that must control if belongs or not no
% the region
%
% Flag: the boolean that is true if P belongs to the region

if ( (P(1) - Descriptors(1)) / Descriptors(3) )^2 + ( (P(2) - Descriptors(2)) / Descriptors(4) )^2 <= 1
    Flag = true;
else
    Flag = false;
end
    
return