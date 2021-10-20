function draw_domain

% Draw all

global B V BInfo
global nB nV

figure (1)
hold on
axis equal

for iB = 1 : nB
    if BInfo(iB,1) == 1 | BInfo(iB,1) == 2 | BInfo(iB,1) == 3
        line = plot ( [ V(B(iB,1),1) V(B(iB,2),1) ] , [ V(B(iB,1),2) V(B(iB,2),2) ] );
        set (line,'Color',[0 0 0],'linewidth',2);
    end
end

% for iV = 1: nV
%     textTria = text (V(iV,1),V(iV,2),num2str(iV));
%     set (textTria,'Color',[0 0 0],'Fontsize',12,'Fontweight','Bold');
% end


return

        
