function [] = find_B_of_memorized_triangle(T);
%
% [] = find_B_of_memorized_triangle(T);
%
% This function find the grid quality B of triangle T. B is calculated as the
% minimum between every possible r/l, where r is the circumradius of the
% triangle and l is the length of one of the triangles.
%
% The function doesn't consider borders that would delimit domain
% boundary (* see attachments for further informations)
%
% Input & output of this function are:
%
% T: the considered triangle
% (output is already written in global variable Tinfo)
% 

global TT BInfo TInfo

% Initialize MinimumLength
MinimumLength = Inf;

% Controls borders
if BInfo ( TT(T,4) , 1 ) == 0
    
    l = find_length_of_memorized_border ( TT(T,5) );
    if l < MinimumLength
        MinimumLength = l;
    end
    
    l = find_length_of_memorized_border ( TT(T,6) );
    if l < MinimumLength
        MinimumLength = l;
    end
    
    if BInfo ( TT(T,5) , 1 ) == 0
        
        l = find_length_of_memorized_border ( TT(T,4) );
        if l < MinimumLength
            MinimumLength = l;
        end
        
    elseif BInfo ( TT(T,6) , 1 ) == 0
        
        l = find_length_of_memorized_border ( TT(T,4) );
        if l < MinimumLength
            MinimumLength = l;
        end
        
    end
    
else
    
    if BInfo ( TT(T,5) , 1 ) == 0
        
        l = find_length_of_memorized_border ( TT(T,4) );
        if l < MinimumLength
            MinimumLength = l;
        end

        l = find_length_of_memorized_border ( TT(T,6) );
        if l < MinimumLength
            MinimumLength = l;
        end
        
        if BInfo ( TT(T,6) , 1 ) == 0
            
            l = find_length_of_memorized_border ( TT(T,5) );
            if l < MinimumLength
                MinimumLength = l;
            end
            
        end
        
    elseif BInfo ( TT(T,6) , 1 ) == 0
        
        l = find_length_of_memorized_border ( TT(T,4) );
        if l < MinimumLength
            MinimumLength = l;
        end
        
        l = find_length_of_memorized_border ( TT(T,5) );
        if l < MinimumLength
            MinimumLength = l;
        end
        
    end
    
end

TInfo(T).B = sqrt(TInfo(T).Circumradius) / MinimumLength;

return

       



    
        

