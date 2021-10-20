function draw_triangles(input)

% Draw all

global B V BInfo
global nB

global foo %%%will be deleted
foo = foo + 1; %%%will be deleted

figure (foo)
hold on
grid on
axis equal

for iB = 1 : nB
    
    line = plot ( [ V(B(iB,1),1) V(B(iB,2),1) ] , [ V(B(iB,1),2) V(B(iB,2),2) ] );
    
    if BInfo(iB,1) == 0
        set (line,'Color',[0 0 0],'linewidth',1);
    else
        set (line,'Color',[1 0 0],'linewidth',1);
    end
    
    if BInfo(iB,2) == 1
        set (line,'linewidth',2.5);
    end
    
end

if input > 10 & mod (input,2) == 1
    
    global BCircle TInfo nT nV
    
    for iT = 1: nT
        textTria = text (TInfo(iT).CG(1),TInfo(iT).CG(2),num2str(iT));
        set (textTria,'Color',[0 0 0],'Fontsize',12,'Fontweight','Bold');
    end
    
    for iB = 1: nB
        textTria = text (BCircle(iB).Circumcenter(1),BCircle(iB).Circumcenter(2),num2str(iB));
        set (textTria,'Color',[1 0 0],'Fontsize',12,'Fontweight','Bold');
    end
    
    for iV = 1: nV
        textTria = text (V(iV,1),V(iV,2),num2str(iV));
        set (textTria,'Color',[0 0 1],'Fontsize',12,'Fontweight','Bold');
    end
    
end

if input > 20
    
    global EB
    
    for iS = 1: EB.nBr
        
        if EB.Br(iS).Unsplitted == true
            
            Limits = EB.Br(iS).SplittingLimits;
            
            limit = plot ( [Limits(1),Limits(1),Limits(2),Limits(2),Limits(1)],[Limits(4),Limits(3),Limits(3),Limits(4),Limits(4)],'Color',[1 0 1],'linewidth',1.5);
            textLimit = text (( Limits(1) + Limits(2) )/2 , ( Limits(3) + Limits(4) )/2,num2str(iS));
            set (textLimit,'Color',[1 0 1],'Fontsize',12,'Fontweight','Bold');
            
        end
        
    end
    
end

return

        
