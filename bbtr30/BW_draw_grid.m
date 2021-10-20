function [] = draw_grid(geom,Fig)
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
    
    LineStyle = plot ( [ V(B(iB,1),1) V(B(iB,2),1) ] , [ V(B(iB,1),2) V(B(iB,2),2) ] );
    
    if BInfo(iB,1) == 0
        % iB is a normal border
        set (LineStyle,'Color',[0 0 0],'linewidth',1);
    elseif mod ( BInfo(iB,3) , 2 ) == 1
        % iB is a Dirichlet border
        set (LineStyle,'Color',[0 0 0],'linewidth',2);
    else
        % iB is a Neumann border
        set (LineStyle,'Color',[0 0 0],'linewidth',2,'LineStyle','--');
    end
    
end

%Mark points
for iV = 1: nV
    
    PointStyle = plot ( V(iV,1) , V(iV,2));
    
    if Nodelist (iV) == 0
        set (PointStyle,'Color',[0 0 0]);
    else
        %iV is a Dirichlet node
        set (PointStyle,'Color',[0 0 0],'Marker','o','MarkerSize',8,'linewidth',2);
    end
    
end


% if input > 10 & mod (input,2) == 1
%     
%     global BCircle TInfo nT nV
%     
%     for iT = 1: nT
%         textTria = text (TInfo(iT).CG(1),TInfo(iT).CG(2),num2str(iT));
%         set (textTria,'Color',[0 0 0],'Fontsize',12,'Fontweight','Bold');
%     end
%     
%     for iB = 1: nB
%         textTria = text (BCircle(iB).Circumcenter(1),BCircle(iB).Circumcenter(2),num2str(iB));
%         set (textTria,'Color',[1 0 0],'Fontsize',12,'Fontweight','Bold');
%     end
%     
%     for iV = 1: nV
%         textTria = text (V(iV,1),V(iV,2),num2str(iV));
%         set (textTria,'Color',[0 0 1],'Fontsize',12,'Fontweight','Bold');
%     end
%     
% end
% 
% if input > 20
%     
%     global EB
%     
%     for iS = 1: EB.nBr
%         
%         if EB.Br(iS).Unsplitted == true
%             
%             Limits = EB.Br(iS).SplittingLimits;
%             
%             limit = plot ( [Limits(1),Limits(1),Limits(2),Limits(2),Limits(1)],[Limits(4),Limits(3),Limits(3),Limits(4),Limits(4)],'Color',[1 0 1],'linewidth',1.5);
%             textLimit = text (( Limits(1) + Limits(2) )/2 , ( Limits(3) + Limits(4) )/2,num2str(iS));
%             set (textLimit,'Color',[1 0 1],'Fontsize',12,'Fontweight','Bold');
%             
%         end
%         
%     end
%     
% end
% 
% return
% 

return

        
