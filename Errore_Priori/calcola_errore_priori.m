function [tot_err] = calcola_errore_priori(true_sol_handle,u)
    global geom
    [zita,csi,eta,omega] = int_nodes_weights(5); %TODO: da rendere globale
    % Calcoliamo una volta per tutte le funzioni di base su tutti i punti
    % della quadratura sul triangolo di riferimento
    
    P = [csi; eta; zita]'; %Per P1
    
    tot_err = 0;
    for e=1:geom.nelements.nTriangles
        points_idx = geom.elements.triangles(e,:);
        coords = geom.elements.coordinates(points_idx,:);
        d = (coords([3,1,2],:)-coords([2,3,1],:)).*[1,-1];
        area = geom.support.TInfo(e).Area;
        
        jjs = geom.pivot.pivot(points_idx);
        contrib_approx = P(:, jjs > 0 )*u(jjs(jjs > 0));
        
        pts = coords(end,:)'+[d(2,1) -d(1,1); -d(2,2) d(1,2)]*[csi; eta];
        contrib_exact = true_sol_handle(pts);
        
        tot_err = tot_err + area*omega*((contrib_exact'-contrib_approx).^2);
        
        %external_pts = coords(end,:)'+[d(2,1) -d(1,1); -d(2,2) d(1,2)]*[0 1 0;0 0 1]
        %plot(external_pts(1,[1,2,3,1]),external_pts(2,[1,2,3,1]),'o--r')
        %plot(pts(1,:),pts(2,:),'*g')
    end
    tot_err = sqrt(tot_err);
end