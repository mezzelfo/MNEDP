function [] = try_rescuing_kidnapped_borders ();
%
% [] = try_rescuing_kidnapped_borders ();
%
% This function must be called only if sure that there esists at least one
% kidnapped border (*see attachments for further informations). The
% function tries to force a kidnapped border into the triangulation by
% inserting its middle point. If this point encroaches a border then 
% the encroached border middle point is inserted


global nKidnappedBorders KidnappedBorders
global B V
global nV

% Find center of last kidnapped border
xV = ( V ( KidnappedBorders (nKidnappedBorders,1) , 1) + V ( KidnappedBorders (nKidnappedBorders,2) , 1) ) / 2;
yV = ( V ( KidnappedBorders (nKidnappedBorders,1) , 2) + V ( KidnappedBorders (nKidnappedBorders,2) , 2) ) / 2;

% check if some borders are encroached
EncroachedBorder = check_encroaching_border ( xV , yV ); %gives reference, -1 for none

if EncroachedBorder ~= -1
    
    % Some border is encroached
    xV = (V (B(EncroachedBorder,1) , 1) + V (B(EncroachedBorder,2) , 1))/2;
    yV = (V (B(EncroachedBorder,1) , 2) + V (B(EncroachedBorder,2) , 2))/2;
    insert_vertex_given_border (xV,yV,EncroachedBorder);

else
    
    % No border is encroached - Find which triangle contains the vertex nV+1
    [T RelativePosition] = prioritary_triangle_research (xV , yV , KidnappedBorders (nKidnappedBorders,3) );

    % Find information about border that we're trying to rescue
    V1 = KidnappedBorders (nKidnappedBorders,1);
    V2 = KidnappedBorders (nKidnappedBorders,2);
    BorderInfo = KidnappedBorders (nKidnappedBorders,4:5);

    % Update kidnapping list
    KidnappedBorders (nKidnappedBorders , 1:5) = [V1 nV+1 T BorderInfo];
    KidnappedBorders (nKidnappedBorders+1 , 1:5) = [V2 nV+1 T BorderInfo];
    nKidnappedBorders = nKidnappedBorders + 1;
    
    % Insert new vertex
    if RelativePosition == 0
        % The vertex lyes inside T
        insert_vertex_given_triangle ( xV , yV , T );
    else
         % The vertex lyes in a border
        insert_vertex_given_border ( xV , yV , RelativePosition );
    end
    
    % Check if some kidanpped border already belong to the grid
    if nKidnappedBorders ~= 0
        check_if_some_kidnapped_is_already_rescued (nV);
    end

end

return