function [Flag] = check_encroaching_border (Vx, Vy);
%
% [Flag] = check_encroaching_border (Vx, Vy);
%
% This function considers a vertex V(Vx,Vy) and discovers if any one of the
% EncroachableBorders would be encroached by the vertex. This function
% finds randomly only one enchroaching border even if there are more than
% one.
%
% Output of this funcion if Flag, that assumes the value -1 if none of the
% border is encroached or the reference of an encroached border.

global EB BCircle

% Find a single Subdomain that contains (Vx,Vy). If it belongs to more than
% one subregion only one is chosen
% -------------------------------------------------------------------------

CurrentDomain = 1;

Found = false;
while Found == false
    
    if EB.Br ( CurrentDomain ).Unsplitted == true
        
        % A subdomain is found!
        Found = true;
        
    else
        
        % Search in subdomains
        I = EB.Br ( CurrentDomain ).KindOfSplitting;
        Vertex = [Vx Vy];
        
        if Vertex(I) <= EB.Br(CurrentDomain).SplittingThreshold;
            CurrentDomain = EB.Tree(CurrentDomain,1);
        else
            CurrentDomain = EB.Tree(CurrentDomain,2);
        end
        
    end
    
end


% Check the effective encroaching
% --------------------------------

Flag = -1;

for iB = EB.Br(CurrentDomain).Leaves
    PositionInCircle = check_if_in_circle (Vx,Vy, BCircle(iB).Circumcenter , BCircle(iB).r2);
    
    if PositionInCircle <= 0
        if PositionInCircle == 0
            
            % The points belongs if the point belongs to the boundary of the circle 
            global SafeMeshing
            if SafeMeshing == true
                
                global V B
                if (V(B(iB,1),1) == Vx & V(B(iB,1),2) == Vy)
                elseif (V(B(iB,2),1) == Vx & V(B(iB,2),2) == Vy)
                else
                    Flag = iB;
                    break
                end
                
            else
                Flag = iB;
                break
            end
            
        else
            % The points lyes within the circle
            Flag = iB;
            break
        end
    end
end 



return