function [] = checks_for_recursive_encroaching (Vertex,EncroachingLevel);
%
% [] = checks_for_recursive_encroaching (Vertex,EncroachingLevel);
%
% This function checks if a recently inserted vertex Vertex encroaches a
% border. In this case a new vertex is inserted; then it is controlled if
% this new vertex encroaches a border or don't. 
% 
% This function may also stop inserting if it understands that a small angle
% is recursively splitted. In this case it stops the splitting, inserting
% two particular verteces (Reference to FIG 7)
%
% Input of this function are:
%
% Vertex: the vertex that may encroach a border
% EncroachingLevel: a counter that alert about the number of recursive call
% of this function


global V nV
global RecursiveEncroachingLimit

if EncroachingLevel == RecursiveEncroachingLimit
    
    % There is a small angle! A particular inserction must be done!
    Border1 = check_encroaching_border  ( V(Vertex,1) , V(Vertex,2) );
    
    global B TT
    
    % Search Border1
    if B(Border1,3) == -1
        T = B(Border1,4);
    else
        T = B(Border1,3);
    end
    
    % Search Border2
    for I = 1:3
        if TT(T,I) == -1
            if TT(T,I+3) ~= Border1
                Border2 = TT(T,I+3);
                break
            end
        end
    end
    
    % Search Center & VP2
    if B (Border1,1) == B (Border2,1)
        Center = B (Border1,1);
        VP2 = B (Border2,2);
        V3 = B (Border1,2);
    elseif B (Border1,1) == B (Border2,2)
        Center = B (Border1,1);
        VP2 = B (Border2,1);
        V3 = B (Border1,2);
    else
        Center = B (Border1,2);
        V3 = B (Border1,1);
        if B (Border2,1) == Center
            VP2 = B (Border2,2);
        else
            VP2 = B (Border2,1);
        end
    end
    
    % Find VP1
    Length2 = find_length_of_memorized_border (Border2);
    Length1 = find_length_of_memorized_border (Border1);
    Lambda1 = Length2/Length1;
    Lambda0 = 1 - Lambda1;
    xP1 = V(Center,1)*Lambda0 + V(V3,1)*Lambda1;
    yP1 = V(Center,2)*Lambda0 + V(V3,2)*Lambda1;

    % Insert VP1
    insert_vertex_given_border (xP1,yP1,Border1);
        
    
else
    
    % There hasn't been found a small angle.
    global B
    
    nVertexToControl = 0;
    VertexToControl = [];
    
    % Discovers which borders are encroached
    EncroachedBorder = check_encroaching_border  ( V(Vertex,1) , V(Vertex,2) ); %gives reference, -1 for none
    while EncroachedBorder ~= -1
        
        % insert the border
        xV = (V (B(EncroachedBorder,1) , 1) + V (B(EncroachedBorder,2) , 1))/2;
        yV = (V (B(EncroachedBorder,1) , 2) + V (B(EncroachedBorder,2) , 2))/2;
        insert_vertex_given_border (xV,yV,EncroachedBorder);
        
        nVertexToControl = nVertexToControl + 1;
        VertexToControl (nVertexToControl) = nV;
        
        EncroachedBorder = check_encroaching_border  ( V(Vertex,1) , V(Vertex,2) ); %gives reference, -1 for none
        
    end
    
    % Disencroach other borders
    for iV = VertexToControl
        checks_for_recursive_encroaching (iV,EncroachingLevel+1);
    end

end

return
    
    

