function tot_err_L2 = calcolaErrorePrioriParabolico(true_sol_handle, sol, Pk, delta_t)
global geom
[zita,csi,eta,omega] = int_nodes_weights(5); %TODO: da rendere globale
% Calcoliamo una volta per tutte le funzioni di base su tutti i punti
% della quadratura sul triangolo di riferimento

switch Pk
    case 'P1'
        P = [csi; eta; zita]';
    case 'P2'
        P = [2.*csi.*(csi-0.5);2.*eta.*(eta-0.5);2.*zita.*(zita-0.5);4*csi.*eta;4.*eta.*zita;4.*csi.*zita]';
    otherwise
        error('Pk only implemented with P1 or P2');
end


% tot_err_L2 = 0;
% for t = 1:size(sol,2)
%     tot_err_L2_t = 0;
%     for e=1:geom.nelements.nTriangles
%         points_idx = geom.elements.triangles(e,:);
%         coords = geom.elements.coordinates(points_idx,:)';
%         d = [1;-1].*(coords(:,[3,1,2])-coords(:,[2,3,1]));
%         B = [d(1,2) -d(1,1); -d(2,2) d(2,1)];
%
%         pts = coords(:,3)+B*[csi; eta];
%         area = geom.support.TInfo(e).Area;
%
%         contrib_exact_L2 = true_sol_handle(pts,(t-1)*delta_t);
%         contrib_approx_L2 = P*sol(points_idx,t);
%         tot_err_L2_t = tot_err_L2_t + 2*area*omega*((contrib_exact_L2'-contrib_approx_L2).^2);
%     end
%     tot_err_L2 = tot_err_L2 + tot_err_L2_t;
% end
% tot_err_L2 = sqrt(delta_t*tot_err_L2);

tot_err_L2 = 0;

for e=1:geom.nelements.nTriangles
    points_idx = geom.elements.triangles(e,:);
    coords = geom.elements.coordinates(points_idx,:)';
    d = [1;-1].*(coords(:,[3,1,2])-coords(:,[2,3,1]));
    B = [d(1,2) -d(1,1); -d(2,2) d(2,1)];
    
    pts = coords(:,3)+B*[csi; eta];
    area = geom.support.TInfo(e).Area;
    
    for t = 1:size(sol,2)
        contrib_exact_L2 = true_sol_handle(pts,(t-1)*delta_t);
        contrib_approx_L2 = P*sol(points_idx,t);
        tot_err_L2 = tot_err_L2 + 2*area*omega*((contrib_exact_L2'-contrib_approx_L2).^2);
    end
end
tot_err_L2 = sqrt(delta_t*tot_err_L2);
end