function [] = load_domain_in_kidnapping_list (Domain,BC,InitialVertexReference)
%
% [] = load_domain_in_kidnapping_list (Domain,BC,InitialVertexReference)
%
% This function tries to find boundary borders in the existing grid. 
% If these borders are not present in the triangulations the funcion
% updates Kidnapping list (*see attachments).
%
% Inputs of this function are:
%
% Domain:
%   Domain.Boundary.Values: references to boundary verteces
%   Domain.Hole.Values: references to verteces belonging to a hole
%   Domain.Segment.Values: references to verteces belonging to a segment
%   Domain.InputVertex coordinates of input verteces
%   (*see attachments to find usage informations on Domain variable)
%
% BC:
%   BC.Boundary.Values: references to boundary BC
%   BC.Hole.Values: references to holes BC
%   BC.Segment.Values: references to segments BC
%   BC.Values: references to BC values
%   BC.InputVertexValues: references to BC to be imposed in InputVertex
%   (*see attachments to find usage informations on BC variable)
%
% InitialVertexReference: the i-esim component of this vector says that
%   InputVertex(i) is inserted in the triangulation as the
%   InitialVertexReference(i).esim element of V array

% declaration of global 
global VB BInfo B
global KidnappedBorders nKidnappedBorders



% -----------------------
% Insert domain boundary
% -----------------------

CurrentList = Domain.Boundary.Values;
CurrentBC = BC.Boundary.Values;
BorderInfoDomain = 1;

n = length(CurrentList);
NextVertex = InitialVertexReference ( CurrentList(1) );
for i = 1 : n
    
    % Find the verteces that are involved in this inserction
    CurrentVertex = NextVertex;
    if i ~= n
        NextVertex = InitialVertexReference ( CurrentList(i+1) );
    else
        NextVertex = InitialVertexReference ( CurrentList(1) );
    end
    BCToSet = abs (CurrentBC(i));
    
    % Check if the border already exists in current triangulation
    Found = false;
    for j = 1 : VB(CurrentVertex).n
        if VB(CurrentVertex).V(j) == NextVertex
            Found = true;
            BInfo(VB(CurrentVertex).B(j),1:3) = [BorderInfoDomain 0 BCToSet]; %as defined in BInfo
            break
        end
    end
    
    % If the border doesn't already exist state that the considererd border
    % is a kidnapped one!
    
    if Found == false
        V1 = CurrentVertex;
        V2 = NextVertex;
        %%% Choose a random triangle near CurrentVertex as the eyewitness (*see attachments)
        T = B(VB( CurrentVertex ).B(1),3); % 1 e 3 are random
        nKidnappedBorders = nKidnappedBorders + 1;
        KidnappedBorders(nKidnappedBorders,:) = [V1 V2 T BorderInfoDomain BCToSet];
    end
end





% -----------------------
% Insert holes boundary
% -----------------------

for CurrentElement = 1: length (Domain.Holes.Hole)

    CurrentList = Domain.Holes.Hole(CurrentElement).Values;
    CurrentBC = BC.Holes.Hole(CurrentElement).Values;
    BorderInfoDomain = 2;
    
    n = length(CurrentList);
    NextVertex = InitialVertexReference ( CurrentList(1) );
    for i = 1 : n
        
        % Find the verteces that are involved in this inserction
        CurrentVertex = NextVertex;
        if i ~= n
            NextVertex = InitialVertexReference ( CurrentList(i+1) );
        else
            NextVertex = InitialVertexReference ( CurrentList(1) );
        end
        BCToSet = abs (CurrentBC(i));
        
        % Check if the border already exists in current triangulation
        Found = false;
        for j = 1 : VB(CurrentVertex).n
            if VB(CurrentVertex).V(j) == NextVertex
                Found = true;
                BInfo(VB(CurrentVertex).B(j),1:3) = [BorderInfoDomain 0 BCToSet]; %as defined in BInfo
                break
            end
        end
        
        % If the border doesn't already exist state that the considererd border
        % is a kidnapped one!
        
        if Found == false
            V1 = CurrentVertex;
            V2 = NextVertex;
            %%% Choose a random triangle near CurrentVertex as the eyewitness (*see attachments)
            T = B(VB( CurrentVertex ).B(1),3); % 1 e 3 are random
            nKidnappedBorders = nKidnappedBorders + 1;
            KidnappedBorders(nKidnappedBorders,:) = [V1 V2 T BorderInfoDomain BCToSet];
        end
    end

end


% -------------------------
% Insert segment boundary
% -------------------------

for CurrentElement = 1: length (Domain.Segments.Segment)

    CurrentList = Domain.Segments.Segment(CurrentElement).Values;
    CurrentBC = BC.Segments.Segment(CurrentElement).Values;
    BorderInfoDomain = 3;
    
    n = length(CurrentList);
    NextVertex = InitialVertexReference ( CurrentList(1) );
    for i = 1 : n - 1
        
        % Find the verteces that are involved in this inserction
        CurrentVertex = NextVertex;
        if i ~= n
            NextVertex = InitialVertexReference ( CurrentList(i+1) );
        else
            NextVertex = InitialVertexReference ( CurrentList(1) );
        end
        BCToSet = abs (CurrentBC(i));
        
        % Check if the border already exists in current triangulation
        Found = false;
        for j = 1 : VB(CurrentVertex).n
            if VB(CurrentVertex).V(j) == NextVertex
                Found = true;
                BInfo(VB(CurrentVertex).B(j),1:3) = [BorderInfoDomain 0 BCToSet]; %as defined in BInfo
                break
            end
        end
        
        % If the border doesn't already exist state that the considererd border
        % is a kidnapped one!
        
        if Found == false
            V1 = CurrentVertex;
            V2 = NextVertex;
            %%% Choose a random triangle near CurrentVertex as the eyewitness (*see attachments)
            T = B(VB( CurrentVertex ).B(1),3); % 1 e 3 are random
            nKidnappedBorders = nKidnappedBorders + 1;
            KidnappedBorders(nKidnappedBorders,:) = [V1 V2 T BorderInfoDomain BCToSet];
        end
    end

end



return