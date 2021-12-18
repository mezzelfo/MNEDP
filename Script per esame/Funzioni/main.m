clc
close all
clear all

global problem

true_sol_handle = @(x) sin(3*x(1,:)).*cos(4*x(2,:));
grad_true_sol_handle = @(x) [3*cos(3*x(1,:)).*cos(4*x(2,:)); -4*sin(3*x(1,:)).*sin(4*x(2,:))];

problem.epsilon = @(x) x(1,:)*0+1;
problem.beta = @(x) [x(2,:);-x(1,:)];
problem.sigma = @(x) x(1,:);
problem.f = @(x) (25+x(1,:)).*true_sol_handle(x)+dot(problem.beta(x), grad_true_sol_handle(x));
problem.bordo_dirichlet = @(x, marker) true_sol_handle(x);
problem.bordo_neumann = @(x, marker) 3*cos(3*x(1,:)).*cos(4*x(2,:));


%% PLOT?
% model = fitlm(log(h_ax), log(errors_ax_P1(:,2)));
% h = model.plot;
% h(4).Color = 'green';
% h

%% P1 vs P2
[h_ax,errors_ax_P1] = convergenzaErrorePriori(...
    'P1',false,false,...
    true_sol_handle,grad_true_sol_handle,...
    'QuadratoMisto'...
    );

[h_ax,errors_ax_P2] = convergenzaErrorePriori(...
    'P2',false,false,...
    true_sol_handle,grad_true_sol_handle,...
    'QuadratoMisto'...
    );

loglog(h_ax,errors_ax_P1(:,2),'-o',h_ax,errors_ax_P2(:,2),'-o')
legend('P1','P2');

polyfit(log(h_ax), log(errors_ax_P1(:,1)),1)
polyfit(log(h_ax), log(errors_ax_P1(:,2)),1)
polyfit(log(h_ax), log(errors_ax_P1(:,3)),1)

polyfit(log(h_ax), log(errors_ax_P2(:,1)),1)
polyfit(log(h_ax), log(errors_ax_P2(:,2)),1)
polyfit(log(h_ax), log(errors_ax_P2(:,3)),1)

%% P1 - SUPG vs no SUPG ============ Fatto meglio sotto
% [h_ax,errors_ax_P1_noSUPG] = convergenzaErrorePriori(...
%     'P1',false,false,...
%     true_sol_handle,grad_true_sol_handle,...
%     'QuadratoMisto'...
%     );
% 
% [h_ax,errors_ax_P1_yesSUPG] = convergenzaErrorePriori(...
%     'P1',true,false,...
%     true_sol_handle,grad_true_sol_handle,...
%     'QuadratoMisto'...
%     );
% 
% loglog(h_ax,errors_ax_P1_noSUPG(:,2),'-o',h_ax,errors_ax_P1_yesSUPG(:,2),'-o')
% legend('without SUPG','with SUPG');
% 
% polyfit(log(h_ax), log(errors_ax_P1_noSUPG(:,1)),1)
% polyfit(log(h_ax), log(errors_ax_P1_noSUPG(:,2)),1)
% polyfit(log(h_ax), log(errors_ax_P1_noSUPG(:,3)),1)
% 
% 
% polyfit(log(h_ax), log(errors_ax_P1_yesSUPG(:,1)),1)
% polyfit(log(h_ax), log(errors_ax_P1_yesSUPG(:,2)),1)
% polyfit(log(h_ax), log(errors_ax_P1_yesSUPG(:,3)),1)

%% P1 - massLumping vs no massLumping
[h_ax,errors_ax_P1_noMassLumping] = convergenzaErrorePriori(...
    'P1',false,false,...
    true_sol_handle,grad_true_sol_handle,...
    'QuadratoMisto'...
    );

[h_ax,errors_ax_P1_yesMassLumping] = convergenzaErrorePriori(...
    'P1',false,true,...
    true_sol_handle,grad_true_sol_handle,...
    'QuadratoMisto'...
    );

loglog(h_ax,errors_ax_P1_noMassLumping(:,2),'-o',h_ax,errors_ax_P1_yesMassLumping(:,2),'-o')
legend('without mass lumping','with mass lumping');

polyfit(log(h_ax), log(errors_ax_P1_noMassLumping(:,1)),1)
polyfit(log(h_ax), log(errors_ax_P1_noMassLumping(:,2)),1)
polyfit(log(h_ax), log(errors_ax_P1_noMassLumping(:,3)),1)


polyfit(log(h_ax), log(errors_ax_P1_yesMassLumping(:,1)),1)
polyfit(log(h_ax), log(errors_ax_P1_yesMassLumping(:,2)),1)
polyfit(log(h_ax), log(errors_ax_P1_yesMassLumping(:,3)),1)

%% P1 - SUPG vs no SUPG, variando epsilon
clc
close all
clear all

global problem

true_sol_handle = @(x) 5+16*x(1,:).*(1-x(1,:)).*x(2,:).*(1-x(2,:));
grad_true_sol_handle = @(x) [   16*x(2,:).*(1-x(2,:)).*(1-2*x(1,:));
    16*x(1,:).*(1-x(1,:)).*(1-2*x(2,:))];

problem.beta = @(x) x*0+[1;1];
problem.sigma = @(x) x(1,:)*0;
problem.bordo_dirichlet = @(P, marker) 5;

epsilon_ax = logspace(1,-6,10);
noSUPG_res = [];
yesSUPG_res = [];

for epsilon = epsilon_ax
    problem.epsilon = @(x) x(1,:)*0+epsilon;
    problem.f = @(x) epsilon*32*(x(1,:).*(1-x(1,:))+x(2,:).*(1-x(2,:))) + dot(problem.beta(x), grad_true_sol_handle(x));
    problem.bordo_neumann = @(x, marker) epsilon*(-16*(x(2,:)-x(2,:).^2));

    [h_ax,errors_ax_P1_noSUPG] = convergenzaErrorePriori(...
        'P1',false,false,...
        true_sol_handle,grad_true_sol_handle,...
        'QuadratoMisto'...
        );

    noSUPG_res(end+1,:,:) = errors_ax_P1_noSUPG;

    [h_ax,errors_ax_P1_yesSUPG] = convergenzaErrorePriori(...
        'P1',true,false,...
        true_sol_handle,grad_true_sol_handle,...
        'QuadratoMisto'...
        );

    yesSUPG_res(end+1,:,:) = errors_ax_P1_yesSUPG;

end

polyfit(log(h_ax), log(yesSUPG_res(end,:,2)),1)

%%
subplot(1,2,1)
plot(h_ax,noSUPG_res(:,:,2)) %Can also look at H1 or Linf error
title('without SUPG')
xlabel('h')
ylabel('L2 Error')

subplot(1,2,2)
plot(h_ax,yesSUPG_res(:,:,2))
title('with SUPG')
xlabel('h')
ylabel('L2 Error')


sgtitle('L2 error of P1 FEM: SUPG vs no SUPG')

%%
orders = [];
for i = 1:size(noSUPG_res,1)
    fit = polyfit(log(h_ax), log(noSUPG_res(i,:,2)),1); %Can also look at H1 or Linf error
    orders(i,1) = fit(1);
    fit = polyfit(log(h_ax), log(yesSUPG_res(i,:,2)),1);
    orders(i,2) = fit(1);
end
semilogx(epsilon_ax,orders(:,1),epsilon_ax,orders(:,2))
legend('without SUPG','with SUPG');
xlabel('epsilon')
ylabel('convergence order of L2 error')