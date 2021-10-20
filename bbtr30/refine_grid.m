function refine_grid (RefiningOptions)
%
% refine_grid (RefiningOptions)
%
% This function refines existing grid inserting verteces in particular
% positions and locally updating the grid
%
% Input of this function is:
%
% RefiningOptions: a structure that explain how to impose the
%   condition that will identify the grid (i.e. the minumum angle checked
%   and the maximum area allowed). It is possible to refine only a
%   particular subregion of the domain. (Particular usages of this variable
%   may be found in the attachments

global TInfo B V nKidnappedBorders
global nT nV nB

global SafeMeshing

% Set variables to start the refining
RefiningIsFinished = false;
RefiningRoof = 0; %%% RefiningRoof sets the maximum number of triangles studied in a single sorting (*see attachments for further informations)
                  %%% 0 means that that Roof must be extimated.
PreviousRoof = 0; %%% says where we stopped in previous sorting (must be <= RefiningRoof)
WhatHasBeenRefined = 0; %%% It will be imposed
                         %%% First of all area is sorted
                         %%% WhatHasBeenRefined assumes following values
                         %%% -1: working on area of domain - still some elements missing |
                         %%% 1: working on area of domain - not any more element missing |

% Find which are current element limits:
nTLimit = nT * 2;
nVLimit = nV * 2;
nBLimit = nB * 2;
reinitialize_variables (nTLimit,nVLimit,nBLimit);

WhereWeMustArriveInRefining = find_end_of_refining (RefiningOptions);

if WhereWeMustArriveInRefining ~= 0
    
    while RefiningIsFinished == false
        
        % Controls if some KidnappedBorder is to insert
        while nKidnappedBorders ~= 0
            try_rescuing_kidnapped_borders
        end
        
        % Create the inserction list
        [SortedInsertionList RefiningRoof WhatHasBeenRefined] = create_inserction_list (RefiningOptions,RefiningRoof,PreviousRoof,WhatHasBeenRefined,WhereWeMustArriveInRefining);

        % Update element lists
        if ( nBLimit < (nB + RefiningRoof) ) | ( nTLimit < (nT + RefiningRoof) ) | ( nVLimit < (nV + RefiningRoof) )
            nTLimit = nT * 2;
            nVLimit = nV * 2;
            nBLimit = nB * 2;
            reinitialize_variables (nTLimit,nVLimit,nBLimit);
        end
        
        disp (['currently inserted ',num2str(nT),' triangles']);
        
        % Continue inserting point until there's area coincidance & roof is
        % not reached
        % -------------------------------------------------------------------
        i = 1;
        AreaCoincidance = true;
        while (i <= RefiningRoof) & AreaCoincidance == true
            
            iT = SortedInsertionList(i,1);
            
            %Check if the refining is made on quality or on Area
            if mod(WhatHasBeenRefined,2) == 1
                ComparisonValue = TInfo(iT).Area;
            else
                ComparisonValue = TInfo(iT).B;
            end
            
            if ComparisonValue == SortedInsertionList(i,2)
                
                
                % Insert a new vertex
                % --------------------
            
                i= i + 1;
            
                xV = TInfo(iT).Circumcenter(1);
                yV = TInfo(iT).Circumcenter(2);
                
                % Check encroaching
                EncroachedBorder = check_encroaching_border  (xV,yV); %gives reference, -1 for none
                
                if EncroachedBorder == -1
                   
                    % In this case no border is encroached and InputVertex is accepted
                    [T PositionVertexTriangle] = prioritary_triangle_research (xV,yV,iT);
                    
                    if PositionVertexTriangle == 0 
                        
                        % In this case the point belongs to the interior of triangle T
                        insert_vertex_given_triangle (xV,yV,T);
                        
                    else
                        
                        % In this case the point belongs to the border PositionVertexTriangle
                        insert_vertex_given_border (xV,yV,PositionVertexTriangle);
                        
                    end
                    
                else
                    
                    % In this case a border is encroached and InputVertex is rejected
                    xV = (V (B(EncroachedBorder,1) , 1) + V (B(EncroachedBorder,2) , 1))/2;
                    yV = (V (B(EncroachedBorder,1) , 2) + V (B(EncroachedBorder,2) , 2))/2;
                    insert_vertex_given_border (xV,yV,EncroachedBorder);
                    
                    % Eventually checks for recursive encroaching
                    if SafeMeshing == true
                        checks_for_recursive_encroaching (nV,1);
                    end

                    
                end
                
            else
                AreaCoincidance = false;
            end
            
        end
        
        % Control limits for further inserction
        % --------------------------------------
        
        % Check if a Roof is really reached or it is necessary to insert other
        % points
        if WhatHasBeenRefined > 0
            if i == (RefiningRoof + 1)
                WhatHasBeenRefined = check_roof_reaching (WhatHasBeenRefined,RefiningOptions);
            else
                WhatHasBeenRefined = - WhatHasBeenRefined;
            end
        end
        
        if WhatHasBeenRefined == WhereWeMustArriveInRefining
            RefiningIsFinished = true;
        end
        
        % Update Roof
        PreviousRoof = i - 1;
        
    end
    
end

return