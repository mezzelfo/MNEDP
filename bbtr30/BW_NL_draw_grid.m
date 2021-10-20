function [] = BW_NL_draw_grid(geom,Fig)
%
% [] = draw_grid(geom,fig,parameters)
%
% This functions draws in figure (fig) the grid passed ad geom.
%
% Input of this function are:
%
% geom: a structure containing every information on the grid
%       (*see attachments for further informations)
% Fig: a reference to the number of the figure
%



% Assign easy-usable variables
V = geom.elements.coordinates;
B = geom.elements.borders;
BInfo = geom.support.BInfo;
Nodelist = geom.pivot.nodelist;
nB = geom.nelements.nBorders;
nV = geom.nelements.nVertexes;


% Set the figure
figure (Fig)
hold on
grid on
axis equal

% Draw lines
for iB = 1 : nB
    
    
    if BInfo(iB,1) ~= 0
    LineStyle = plot ( [ V(B(iB,1),1) V(B(iB,2),1) ] , [ V(B(iB,1),2) V(B(iB,2),2) ] );
        % iB is a Dirichlet border
        set (LineStyle,'Color',[0 0 0],'linewidth',2);
    end
    
end

%Mark points
% for iV = 1: nV
%     
%     PointStyle = plot ( V(iV,1) , V(iV,2));
%     
%     if Nodelist (iV) == 0
%         set (PointStyle,'Color',[0 0 0]);
%     else
%         %iV is a Dirichlet node
%         set (PointStyle,'Color',[0 0 0],'Marker','o','MarkerSize',8,'linewidth',2);
%     end
%     
% end



return

        
