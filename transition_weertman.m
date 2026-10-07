% Plots of the transition region for a Weertman sliding law (Figures 3 and 6)

% Load default parameter values
[scales, params, funcs] = default_values();
Hd = scales.Hd;   xd_km = scales.xd_km; Nd_MPa = scales.Nd_MPa; ud_myr = scales.ud_myr;
params.Sigma = 0;
delta = params.delta;   eps = params.eps;   E = params.E;
lambda = params.lambda;     n = params.n;   mW = params.mW;
C = params.C;   

% Functions for ODEs and nullcline
HNf = @(h, n, kap) [ h.^(-mW-1) ; h.^(-mW-1) - n./(1+ kap.*n) ];
hnull = @(n, kap) (n./(1+kap*n)).^(-1/(mW+1));

% Phase plane for kappa = 0 (Figure 3)
figure(1); clf; set(gcf(),'Position',[1 1 1000 400])
aK0 = make_nice_axes(subplot('Position',[0.07 0.15 0.6 0.8]),'$\hat H$','$\hat N$','',5);
ylim(aK0,[0 1])

% Quiver plot of phase trajectories
hs = 0.5:0.5:5;
ns = 0.1:0.1:1;
[hh, nn] = meshgrid(hs, ns);
HxiNxi = HNf(hh,nn,0);
Hxi  = HxiNxi(1:length(ns),:);
Nxi  = HxiNxi(length(ns)+1:end,:);
quiver(aK0,hh,nn,Hxi,Nxi,'k-');

% plot nullcline(s)
ns2 = 0.05:0.01:1;
pnull = plot(aK0,hnull(ns2,0),ns2,'k--','LineWidth',2);


% example trajectories
starts = [0.001 0.001; 1 0.001; 3 1; 3 0.001; 2 0.001];
for i = 1:size(starts,1)
    [~, HN] = ode45( @(t,hn) HNf(hn(1), hn(2), 0) ,[0 20], [starts(i,1) starts(i,2)]);
    ptraj = plot(HN(:,1),HN(:,2),'k-','LineWidth',1);
end

mrates = [.5 1 2 4];
clrs = colormap(cmocean('haline',length(mrates)));
for mrate = mrates
    
    % Plot the scaled version of the actual solution
    [H_, N_, ~] = get_scaled_sol(mrate, 0, mW, E, delta, Hd);
    plot(H_,N_,'color',clrs(mrate==mrates,:),'LineWidth',2);
    
    % % Phase trajectory for the same initial condition
    % [~, HN] = ode45( @(t,hn) HNf(hn(1), hn(2), 0) ,[0 20], [H_(end) 0.001]);
    % plot(HN(:,1),HN(:,2),'r-','LineWidth',1)
end

legend(aK0,aK0.Children([end-1 5:-1:1]),{'$\hat N$ nullcline','Example trajectories','$m=0.5$ mm yr$^{-1}$','$m=1$ mm yr$^{-1}$','$m=2$ mm yr$^{-1}$','$m=4$ mm yr$^{-1}$'},"fontsize", 14, 'interpreter', 'latex');

% Panel (b) (asymptotic behaviour)
H0s = [0.001 0.1:0.1:8];
Nmax = zeros(size(H0s));
for H0 = H0s
    [~, HN] = ode45( @(t,hn) HNf(hn(1), hn(2), 0) ,[0 20], [H0 0.001]);
    Nmax(H0==H0s) = max(HN(:,2));
end

aK0b = make_nice_axes(subplot('Position',[0.75 0.15 0.23 0.8]),'$\hat{H}_g$','$\hat{N}_{max}$','',8);
plot(aK0b,H0s,Nmax,'k-d');
plot(aK0b,H0s,H0s.^(-mW-1),'r-');
legend(aK0b,{'${\hat N}_{max}(\hat{H}_g)$','$\hat{H}_g^{-m_W-1}$'},'FontSize',14,'Interpreter','latex')
ylim(aK0b,[0 1])

text(aK0,0.25,0.9,'(a)','Interpreter','latex','BackgroundColor','w','FontSize',14);
text(aK0b,0.5,0.9,'(b)','Interpreter','latex','BackgroundColor','w','FontSize',14);
exportgraphics(figure(1),'figures/transition_weertman_K0.png','Resolution',300)

% Additional figure: Phase plane for nonzero kappa
figure(2); clf;
a2 = make_nice_axes(gca(),'$\hat H$','$\hat N$','',6);
ylim([0 2])

[H_, N_, kappa] = get_scaled_sol(1, 1, mW, E, delta, Hd);

% Quiver plot
hs = 0.5:0.5:6;
ns = 0.1:0.1:2;
[hh, nn] = meshgrid(hs, ns);
HxiNxi = HNf(hh,nn,kappa);
Hxi  = HxiNxi(1:length(ns),:);
Nxi  = HxiNxi(length(ns)+1:end,:);
quiver(a2,hh,nn,Hxi,Nxi,'k-');

title(sprintf('$\\hat \\kappa= %g$',kappa),'Interpreter','latex')

% plot nullcline(s)
ns2 = 0.05:0.01:2;
pnull2 = plot(hnull(ns2,kappa),ns2,'k--','LineWidth',2);

% some trajectories
starts = [0.001 0.001; 1 0.001; 3 0.001 ; 2 0.001 ; 5 0.001];
for i = 1:size(starts,1)
    [~, HN] = ode45( @(t,hn) HNf(hn(1), hn(2), kappa) ,[0 30], [starts(i,1) starts(i,2)]);
    ptraj2 = plot(HN(:,1),HN(:,2),'k-','LineWidth',1);
end

% Scaled solution
psol2 = plot(H_,N_,'b');
legend([pnull2 ptraj2 psol2],{'N nullcline','Example trajectories','Scaled solution (m=1, K=1)'} )

% Figure 6 (commented out sections plot extra values of kappa): 
figure(3); clf; set(gcf(),'Position',[1 1 840 620])
a3 = make_nice_axes(subplot('Position',[0.08 0.6 0.9 0.38]),'$\hat H$','$\hat N$','',8); 
ylim([0 1]);

% for kappa = 0:3
%     plot(hnull(ns2,kappa),ns2,'k--','LineWidth',2);
% end

hs = 0.5:0.5:8;
ns = 0.1:0.1:1;
[hh, nn] = meshgrid(hs, ns);

clrs = colormap(cmocean('dense',4));
 
[H_, N_, kappa1] = get_scaled_sol(8, 0, mW, E, delta, Hd);
% plot(a3,hnull(ns2,kappa1),ns2,'--','LineWidth',2,'Color',clrs(2,:));
% psol0  = plot(a3,H_,N_,'Color',clrs(2,:),'LineWidth',1.5);
% 
% [~, HN] = ode45( @(t,hn) HNf(hn(1), hn(2), kappa1) ,[0 50], [1 0.001]);
% plot(a3,HN(:,1),HN(:,2),'-','Color',clrs(2,:),'LineWidth',.75);
% [~, HN] = ode45( @(t,hn) HNf(hn(1), hn(2), kappa1) ,[0 50], [4 0.001]);
% plot(a3,HN(:,1),HN(:,2),'-','Color',clrs(2,:),'LineWidth',.75);
% 
% HxiNxi = HNf(hh,nn,kappa1);
% Hxi  = HxiNxi(1:length(ns),:);
% Nxi  = HxiNxi(length(ns)+1:end,:);
% %quiver(a3,hh,nn,Hxi,Nxi,'Color',clrs(2,:));

[H_, N_, kappa2] = get_scaled_sol(8, 2, mW, E, delta, Hd);
pnull = plot(a3,hnull(ns2,kappa2),ns2,'--','LineWidth',2,'Color',clrs(3,:));
psol2 = plot(a3,H_(160:end),N_(160:end),'Color',clrs(3,:),'LineWidth',1.5);

[~, HN] = ode45( @(t,hn) HNf(hn(1), hn(2), kappa2) ,[0 50], [1 0.001]);
plot(a3,HN(:,1),HN(:,2),'-','Color',clrs(3,:),'LineWidth',.75);
[~, HN] = ode45( @(t,hn) HNf(hn(1), hn(2), kappa2) ,[0 50], [4 0.001]);
ptraj = plot(a3,HN(:,1),HN(:,2),'-','Color',clrs(3,:),'LineWidth',.75);

HxiNxi = HNf(hh,nn,kappa2);
Hxi  = HxiNxi(1:length(ns),:);
Nxi  = HxiNxi(length(ns)+1:end,:);
pq = quiver(a3,hh,nn,Hxi,Nxi,'Color',clrs(3,:));

[H_, N_, kappa3] = get_scaled_sol(8, 4, mW, E, delta, Hd);
% pnull = plot(a3,hnull(ns2,kappa3),ns2,'--','LineWidth',2,'Color',clrs(4,:));
% psol2 = plot(a3,H_,N_,'Color',clrs(4,:),'LineWidth',1.5);
% 
% [~, HN] = ode45( @(t,hn) HNf(hn(1), hn(2), kappa3) ,[0 50], [1 0.001]);
% plot(a3,HN(:,1),HN(:,2),'-','Color',clrs(4,:),'LineWidth',.75);
% [~, HN] = ode45( @(t,hn) HNf(hn(1), hn(2), kappa3) ,[0 50], [4 0.001]);
% ptraj = plot(a3,HN(:,1),HN(:,2),'-','Color',clrs(4,:),'LineWidth',.75);
% 
% HxiNxi = HNf(hh,nn,kappa3);
% Hxi  = HxiNxi(1:length(ns),:);
% Nxi  = HxiNxi(length(ns)+1:end,:);
% %pq = quiver(a3,hh,nn,Hxi,Nxi,'Color',clrs(4,:));

% legend([pq pnull ptraj psol2],{'Trajectory directions','N nullclines','Example trajectories','$m=8$ mm yr$^{-1}$'} ,...
%     'FontSize',14,'Interpreter','latex','location','northeast')
legend([pnull ptraj psol2],{'N nullclines','Example trajectories','$m=8$ mm yr$^{-1}$'} ,...
     'FontSize',14,'Interpreter','latex','location','northeast')

text(a3,0.2,0.9,'(a)','Interpreter','latex','BackgroundColor','w','FontSize',14);

% Panel b
H0s = [0.001 0.1:0.1:8];
Nmax1 = zeros(size(H0s));
Nmax2 = zeros(size(H0s));
Nmax3 = zeros(size(H0s));
for H0 = H0s
    [~, HN] = ode45( @(t,hn) HNf(hn(1), hn(2), kappa1) ,[0 20], [H0 0.001]);
    Nmax1(H0==H0s) = max(HN(:,2));
    [~, HN] = ode45( @(t,hn) HNf(hn(1), hn(2), kappa2) ,[0 20], [H0 0.001]);
    Nmax2(H0==H0s) = max(HN(:,2));
    [~, HN] = ode45( @(t,hn) HNf(hn(1), hn(2), kappa3) ,[0 20], [H0 0.001]);
    Nmax3(H0==H0s) = max(HN(:,2));
end

%a3b = make_nice_axes(subplot('Position',[0.75 0.15 0.23 0.8]),'$\hat{H}_g$','$\hat{N}_{max}$','',8);
a3b = make_nice_axes(subplot('Position',[0.08 0.1 0.42 0.38]),'$\hat{H}_g$','$\hat{N}_{max}$','',8);
plot(a3b,H0s,Nmax1,'-d','Color',clrs(2,:));
plot(a3b,H0s,Nmax2,'-d','Color',clrs(3,:));
plot(a3b,H0s,Nmax3,'-d','Color',clrs(4,:));
ylim(a3b,[0 2])

legend(a3b,{'$\hat \kappa = 0$ ($K=0$)', ['$\hat \kappa = ' num2str(kappa2,3) '$ ($K=2$)'], ['$\hat \kappa = ' num2str(kappa3,3) '$ ($K=4$)']}, 'FontSize',14,'Interpreter','latex','location','northeast')
text(a3b,0.2,1.8,'(b)','Interpreter','latex','BackgroundColor','w','FontSize',14);

% Panel c
a3c = make_nice_axes(subplot('Position',[0.56 0.1 0.42 0.38]),'$\hat{\kappa}$','$\hat{N}_{max}$','');
kappas = 0:0.25:10;
Nmax1b = zeros(size(kappas));
Nmax2b = zeros(size(kappas));
Nmax3b = zeros(size(kappas));
for kappa = kappas
    [~, HN] = ode45( @(t,hn) HNf(hn(1), hn(2), kappa) ,[0 100], [0.001 0.001]);
    Nmax1b(kappa==kappas) = max(HN(:,2));
    [~, HN] = ode45( @(t,hn) HNf(hn(1), hn(2), kappa) ,[0 100], [1 0.001]);
    Nmax2b(kappa==kappas) = max(HN(:,2));
    [~, HN] = ode45( @(t,hn) HNf(hn(1), hn(2), kappa) ,[0 100], [2 0.001]);
    Nmax3b(kappa==kappas) = max(HN(:,2));
end

plot(a3c,kappas,Nmax1b,'-d','Color','k');
plot(a3c,kappas,Nmax2b,':d','Color','k');
%plot(a3c,kappas,Nmax3b,'-s','Color','k');

plot(a3c,kappas(9:end),((mW+1)/(mW+2))*kappas(9:end).^(1/(mW+1)),'-','Color','r');
plot(a3c,kappas(9:end),((mW+1)/(mW+2))*kappas(9:end).^(1/(mW+1))-1+1.^(mW+2)./((mW+2)*kappas(9:end)),':','Color','r');
%plot(a3c,kappas,((mW+1)/(mW+2))*kappas.^(1/(mW+1))-2+2.^(mW+2)./((mW+2)*kappas),'-s','Color','b');

%legend(a3c,{'Numerical, $\hat{H}_g=0$','Numerical, $\hat{H}_g=1$','Numerical, $\hat{H}_g=2$','Asymptotic, $\hat{H}_g=0$','Asymptotic, $\hat{H}_g=1$','Asymptotic, $\hat{H}_g=2$'} , 'FontSize',14,'Interpreter','latex')
legend(a3c,{'${\hat N}_{max}$, $\hat{H}_g=0$','${\hat N}_{max}$, $\hat{H}_g=1$','Asymp., $\hat{H}_g=0$','Asymp., $\hat{H}_g=1$'},...
    'FontSize',14,'Interpreter','latex','location','northwest','NumColumns',2)
ylim(a3c,[0 5])
text(a3c,0.25,3,'(c)','Interpreter','latex','BackgroundColor','w','FontSize',14);
exportgraphics(figure(3),'figures/transition_weertman_m8K.png','Resolution',300)

% Function to obtain a scaled solution, along with values of kappa for a
% particular scenario
function [H_, N_, kappa] = get_scaled_sol(mrate, K, mW, E, delta, Hd)
    load(sprintf('steady_states/SS_m%g_K%g.mat',mrate,K),'HUxN');
    J = 300;
    
    H = HUxN(1:J+1);
    xg = HUxN(2*J+4);
    N = HUxN(2*J+5:3*J+5);
    
    Hscale = (xg^(mW-1)*E/mrate)^(1/(mW+2));
    Nscale = ((1-delta)/E)*(xg^(mW-1)*E/mrate)^(1/(mW+2));
    H_ = H/Hscale;
    N_ = N/Nscale;
    kappa = Nscale*K*2000/Hd;

    beta = 0.1;
    xscale = E / (mrate*xg) ;
    beta_ = beta * xscale / Hscale;
end