% Plots results of parameter sweep (additional panels for Figures 2 and 5)
[scales, params, funcs] = default_values();
delta = params.delta;
mW= params.mW;
bedf = funcs.bedf;
mrates  = 0:0.25:20;
Ks = 0:0.25:8;

 C = params.C;
load("variable_grid_300.mat","XH","XU");
J = 300;

% commented out Figure showing where the parameter sweep has obtained
% solutions 
% bedf = funcs.bedf;
% [mm, KK] = meshgrid(mrates,Ks);
% existence = false(size(mm));
% 
% for l = 1:length(mrates)
%     for ll = 1:length(Ks) 
%         existence(ll,l) = isfile(sprintf('steady_states/SS_m%g_K%g.mat', mrates(l), Ks(ll)));
%     end
% end
% 
% pcolor(mm,KK,double(existence))

% Additional figure plotting the steady state ice sheet thickness for
% various m (K=0), and comparing this to the approximate solutions C/mx and (1-delta)C/mx for large m 
figure(3); clf;
aH = subplot(2,1,1); hold on;
aH2 = subplot(2,1,2); hold on;

clrs1 = cmocean('-thermal',length(mrates));
for l = 1:length(mrates)
    load(sprintf('steady_states/SS_m%g_K%g.mat', mrates(l), 0),'HUxN');
    plot(aH,HUxN(2*J+4)*XH,HUxN(1:J+1) + bedf(HUxN(2*J+4)*XH),'color',clrs1(l,:));
    plot(aH2,HUxN(2*J+4)*XH,HUxN(1:J+1)*mrates(l),'color',clrs1(l,:));
end

plot(aH2,(0.3:0.01:1.8),C./(0.3:0.01:1.8),'k');
plot(aH2,(0.3:0.01:1.8),(1-delta)*C./(0.3:0.01:1.8),'k--');

% Additional panel for figure 2
% Generating mesh
[Xmesh, mmesh] = meshgrid(XH,mrates);
xmesh = zeros(size(Xmesh));
Hmesh = zeros(size(Xmesh));
Umesh = zeros(size(Xmesh));
Nmesh = zeros(size(Xmesh));

% Populating mesh
for l = 1:length(mrates)
    load(sprintf('steady_states/SS_m%g_K%g.mat', mrates(l), 0),'HUxN');
    xmesh(l,:) = HUxN(2*J+4)*XH;
    Hmesh(l,:) = HUxN(1:J+1);
    Umesh(l,:) = interp1(XU,HUxN(J+2:2*J+3),XH);
    Nmesh(l,:) = HUxN(2*J+5:3*J+5);
end

R = (C*Nmesh- Umesh.^mW )./(Umesh.^mW + C*Nmesh);

% Figure creation
figure(5); clf;  set(gcf(),'Position',[1   1  1260   430]);
aslm = make_nice_axes(gca(),'$x$ / km', '$m$ / mm yr$^{-1}$', '', 903);
ylim(aslm,[0 20]);
pcolor(aslm,scales.xd_km*xmesh,scales.md_mmyr*mmesh,R,'EdgeColor','none','facecolor','interp');
colormap(cmocean('balance'));

plot(aslm,[0 scales.xd_km*1.8],[1 1],'k--','linewidth',1);
plot(aslm,[0 scales.xd_km*1.8],[8 8],'k--','linewidth',1);

text(aslm,10,10,'Weertman-like','BackgroundColor','w','interpreter','latex','fontsize',14);
text(aslm,400,2,'Weertman-like','BackgroundColor','w','interpreter','latex','fontsize',14);
text(aslm,400,14,'Coulomb-like','BackgroundColor','w','interpreter','latex','fontsize',14);
text(aslm,800,10,{'Extensional', 'stress BL'},'color','k','interpreter','latex','fontsize',14);
text(aslm,750,14,{'Transition', 'region'},'color','k','interpreter','latex','fontsize',14);

plot(aslm,[820 780],[5 12.5],'k-','linewidth',.5);
plot(aslm,[875 850],[3 8.5],'k-','linewidth',.5);

text(aslm,10,1,'$m=1$','Interpreter','latex','BackgroundColor','w','FontSize',14);
text(aslm,10,8,'$m=8$','Interpreter','latex','BackgroundColor','w','FontSize',14);

text(aslm,10,19,'(c)','Interpreter','latex','BackgroundColor','w','FontSize',14);

cb = colorbar(aslm,'eastoutside','FontSize',14,'TickLabelInterpreter','latex');
cb.Label.String = '$R$';
cb.Label.FontSize = 14;
cb.Label.Interpreter = 'latex';

exportgraphics(aslm,'figures/param_sweep_m.png','Resolution',300)

% Additional panel for figure 2
% Generating mesh
[Xmesh2, Kmesh] = meshgrid(XH,Ks);
xmesh2 = zeros(size(Xmesh2));
Hmesh2 = zeros(size(Xmesh2));
Umesh2 = zeros(size(Xmesh2));
Nmesh2 = zeros(size(Xmesh2));

% Populating mesh
for l = 1:length(Ks)
    load(sprintf('steady_states/SS_m%g_K%g.mat', 8, Ks(l)),'HUxN');
    xmesh2(l,:) = HUxN(2*J+4)*XH;
    Hmesh2(l,:) = HUxN(1:J+1);
    Umesh2(l,:) = interp1(XU,HUxN(J+2:2*J+3),XH);
    Nmesh2(l,:) = HUxN(2*J+5:3*J+5);
end

R2 = (C*Nmesh2- Umesh2.^mW )./(Umesh2.^mW + C*Nmesh2);

% Figure creation
figure(6); clf; set(gcf(),'Position',[1   1  1260   430]);
aslm2 = make_nice_axes(gca(),'$x$ / km', '$K$', '', 903);
ylim(aslm2,[0 8]);
set(aslm2,"FontSize",14,"TickLabelInterpreter",'latex')
pcolor(aslm2,scales.xd_km*xmesh2,Kmesh,R2,'EdgeColor','none','facecolor','interp');
colormap(cmocean('balance'));

text(aslm2,10,5,'Weertman-like','BackgroundColor','w','interpreter','latex','fontsize',14);
text(aslm2,225,3,'Weertman-like','BackgroundColor','w','interpreter','latex','fontsize',14);
text(aslm2,625,7.5,'Weertman-like','BackgroundColor','w','interpreter','latex','fontsize',14);
text(aslm2,450,5,'Coulomb-like','BackgroundColor','w','interpreter','latex','fontsize',14);
text(aslm2,810,.75,{'Extensional', 'stress BL'},'color','k','interpreter','latex','fontsize',14);
text(aslm2,825,2.75,{'Transition', 'region'},'color','k','interpreter','latex','fontsize',14);
plot(aslm2,[830 810],[1 1.5],'k-','linewidth',.5);
plot(aslm2,[825 790],[2.5 1.5],'k-','linewidth',.5);

plot(aslm2,[0 scales.xd_km*1.8],0*[1 1],'k--','linewidth',1);
plot(aslm2,[0 scales.xd_km*1.8],2*[1 1],'k--','linewidth',1);
plot(aslm2,[0 scales.xd_km*1.8],4*[1 1],'k--','linewidth',1);
plot(aslm2,[0 scales.xd_km*1.8],6*[1 1],'k--','linewidth',1);

text(aslm2,10,7.5,'(e)','Interpreter','latex','BackgroundColor','w','FontSize',14);

text(aslm2,10,0.5,'$K=0$','Interpreter','latex','BackgroundColor','w','FontSize',14);
text(aslm2,10,2,'$K=2$','Interpreter','latex','BackgroundColor','w','FontSize',14);
text(aslm2,10,4,'$K=4$','Interpreter','latex','BackgroundColor','w','FontSize',14);
text(aslm2,10,6,'$K=6$','Interpreter','latex','BackgroundColor','w','FontSize',14);

cb = colorbar(aslm2,'eastoutside','FontSize',14,'TickLabelInterpreter','latex');
cb.Label.String = '$R$';
cb.Label.FontSize = 14;
cb.Label.Interpreter = 'latex';
exportgraphics(aslm2,'figures/param_sweep_m8K.png','Resolution',300)