function [Di,Ne,Nodelist] = evaluate_boundary_conditions (Domain,BC);
%
% [Di Ne Nodelist] = evaluate_boundary_conditions (Domain,BC)
%
% This function find the correct boundary conditions to assign to a domain.
%
% Input & Output of this function are:
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
% Nodelist: an array that keep information of BCs to impose in every vertex
%   Nodelist(i) = 0 : the i-vertex has no BC related to it
%   Nodelist(i) = j : the i-vertex has BC = j (where j is a reference to 
%       BC.Values). This condition must be a Dirichlet one
%
% Ne: this is a list of borders whose assigned BC are Neumann ones.
%   The i-esim row of Ne is [i,j], where i is the i-esim border and j is
%   the reference to the condition related to it.
%
% Di: this is a list of verteces whose assigned BC are Dirichlet ones.
%   The i-esim row of De is [i,j], where i is the i-esim vertex and j is
%   the reference to the condition related to it.
%

global nB nV
global B V BInfo

% Define the searched variables
Nodelist = zeros (nV,1);
nNe = 0;
Ne = zeros (nV,2);
nDi = 0;
Di = zeros (nB,2);

% Check every border
% -------------------

for iB = 1:nB
    
    if BInfo (iB,1) ~= 0
        % The border has some BC setted on it
        
        if mod ( BInfo(iB,3) , 2 ) == 1
            
            % The condition is a Dirichlet one (reference is an even number)
            
            if Nodelist ( B(iB,1) ) == 0
                % The vertex hasn't still been recognized as a Dirichlet one
                Nodelist ( B(iB,1) ) = BInfo(iB,3);
                nDi = nDi + 1;
                Di (nDi,:) = [B(iB,1) BInfo(iB,3)];
            end

            if Nodelist ( B(iB,2) ) == 0
                % The vertex hasn't still been recognized as a Dirichlet one
                Nodelist ( B(iB,2) ) = BInfo(iB,3);
                nDi = nDi + 1;
                Di (nDi,:) = [B(iB,2) BInfo(iB,3)];
            end

        else
            
            % The condition is a Neumann one (reference is an odd number)
            nNe = nNe + 1;
            Ne (nNe,:) = [iB BInfo(iB,3)];
            
        end
        
    end
            
end

% Impose remaining InputVertex
% -----------------------------

nIV = size(Domain.InputVertex,1);

for iIV = 1 : nIV
    
    for jV = 1: nV
        
        if Domain.InputVertex(iIV,1) == V (jV,1)
            if Domain.InputVertex(iIV,2) == V (jV,2)
                
                % a InputVertex has been found
                
                if Nodelist(jV) == 0
                    
                    if BC.InputVertexValues(iIV) ~= 0
                        
                        Nodelist ( jV ) = BC.InputVertexValues(iIV);
                        nDi = nDi + 1;
                        Di (nDi,:) = [jV BC.InputVertexValues(iIV)];
                        
                    end
                    
                else
                    
                    % Nodelist(jV) ~= 0. Find jDi, index that says where jV
                    % is memorized in Di
                    
                    for jDi = 1:nDi
                        if Di(jDi,1) == jV
                            break
                        end
                    end
                    
                    if BC.InputVertexValues(iIV) ~= 0
                        
                        % Di is to update
                        Di (jDi,2) = BC.InputVertexValues(iIV);
                        Nodelist ( jV ) = BC.InputVertexValues(iIV);
                        
                    else
                        
                        % A component of Di must be removed
                        Di ( jDi:nDi-1 , : ) = Di ( jDi+1:nDi , : );
                        nDi = nDi - 1;
                        Nodelist ( jV ) = 0;
                        
                    end
                    
                end
                
                break
                
            end
        end
        
    end
    
end


% clean variables
Temp = Di(1:nDi,:);
clear Di
Di = Temp;
clear Temp;
Temp = Ne(1:nNe,:);
clear Ne
Ne = Temp;

return