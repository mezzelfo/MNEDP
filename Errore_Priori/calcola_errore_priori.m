function [tot_err_inf, tot_err_L2, tot_err_H1] = calcola_errore_priori(true_sol_handle, grad_true_sol_handle, sol)
    global geom
    [zita,csi,eta,omega] = int_nodes_weights(5); %TODO: da rendere globale
    % Calcoliamo una volta per tutte le funzioni di base su tutti i punti
    % della quadratura sul triangolo di riferimento
    
    P = [csi; eta; zita]'; %Per P1
    grad_P = [[1 0];[0 1];[-1 -1]]';
    grad_P = repmat(grad_P,1,1,length(omega));
    
    
    
    tot_err_L2 = 0;
    tot_err_H1 = 0;
    tot_err_inf = 0;
    for e=1:geom.nelements.nTriangles
        points_idx = geom.elements.triangles(e,:);
        coords = geom.elements.coordinates(points_idx,:);
        d = (coords([3,1,2],:)-coords([2,3,1],:)).*[1,-1];
        B = [d(2,1) -d(1,1); -d(2,2) d(1,2)];
        pts = coords(end,:)'+B*[csi; eta];
        area = geom.support.TInfo(e).Area;
        
        %% L2
        contrib_exact_L2 = true_sol_handle(pts);
        contrib_approx_L2 = P*sol(points_idx);
        tot_err_L2 = tot_err_L2 + 2*area*omega*((contrib_exact_L2'-contrib_approx_L2).^2);
        
        %% H1
        contrib_exact_H1 = grad_true_sol_handle(pts);
        contrib_approx_H1 = inv(B')*squeeze(pagemtimes(grad_P, sol(points_idx)));
        tot_err_H1 = tot_err_H1 + 2*area*omega*(vecnorm(contrib_exact_H1-contrib_approx_H1,2).^2)';
        
        %% inf
        tot_err_inf = max(tot_err_inf, max(abs(sol(points_idx)' - true_sol_handle(coords'))));
        
    end
    tot_err_L2 = sqrt(tot_err_L2);
    tot_err_H1 = sqrt(tot_err_H1);

end