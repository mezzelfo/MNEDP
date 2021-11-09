[zita,csi,eta,omega] = int_nodes_weights(5);
P = [csi; eta; zita]'; %Per P1

%profile on
totjjs = geom.pivot.pivot(geom.elements.triangles);
totu = (u(max(totjjs,1)).*(totjjs > 0))';
totcontrib_approx = P*totu;

totcoords = reshape(...
    geom.elements.coordinates(geom.elements.triangles(:),:),...
    [geom.nelements.nTriangles,3,2]);
totcoords = permute(totcoords,[3,2,1]);
totpts = totcoords(:,3,:) + pagemtimes(totcoords, [1 0; 0 1; -1 -1] * [csi; eta]);
totcontrib_exact = squeeze(true_sol_handle(totpts));

totarea = [geom.support.TInfo(:).Area];
toterr = omega*((totcontrib_exact-totcontrib_approx).^2)*totarea';
toterr = sqrt(2*toterr);
%profile viewer


% tic
% tot_err = 0;
% for e=1:geom.nelements.nTriangles
%     points_idx = geom.elements.triangles(e,:);
%     coords = geom.elements.coordinates(points_idx,:);
%     d = (coords([3,1,2],:)-coords([2,3,1],:)).*[1,-1];
%     area = geom.support.TInfo(e).Area;
% 
%     jjs = geom.pivot.pivot(points_idx);
%     contrib_approx = P(:, jjs > 0 )*u(jjs(jjs > 0));
% 
%     pts = coords(end,:)'+[d(2,1) -d(1,1); -d(2,2) d(1,2)]*[csi; eta];
%     contrib_exact = true_sol_handle(pts);
% 
%     tot_err = tot_err + 2*area*omega*((contrib_exact'-contrib_approx).^2);
% 
%     %external_pts = coords(end,:)'+[d(2,1) -d(1,1); -d(2,2) d(1,2)]*[0 1 0;0 0 1]
%     %plot(external_pts(1,[1,2,3,1]),external_pts(2,[1,2,3,1]),'o--r')
%     %plot(pts(1,:),pts(2,:),'*g')
% end
% tot_err = sqrt(tot_err);
% toc