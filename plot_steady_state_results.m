% Finds a selection of steady state solutions
% Generates Figures 2, 3, 4 and 5 using these

% Default parameters and functions
[scales, params, funcs] = default_values();
Hd = scales.Hd;   xd_km = scales.xd_km; Nd_MPa = scales.Nd_MPa; ud_myr = scales.ud_myr;
md_mmyr = scales.md_mmyr;
delta = params.delta;  E = params.E;  
params.Sigma = 0;
mW = params.mW;

params.lambda = 1e-3;

bedf = funcs.bedf;
bedxf = funcs.bedxf;
Hsbf = funcs.Hsbf;

% Alternative sedimentary basin geometry (Fig. 5)
Hsbf2 = @(x) 750/Hd -  (500/Hd)*cos(x*2*pi/.5); 

% Mesh grid used for numerical solutions
load("variable_grid_300.mat","XA","XH","XU");
J = 300;

% Folders to store generated output
if ~isfolder('steady_states');   mkdir('steady_states');  end
if ~isfolder('figures');   mkdir('figures');  end

% Colormaps for plots
clrs = flip(cmocean('thermal',7));
clrs2 = cmocean('haline',8);

%% Finding steady states

% First: find solution for m=K=0
funcs.meltf = @(x) zeros(size(x));
funcs.intmeltf = @(x) zeros(size(x));
params.K = 0;

if ~isfile('steady_states/SS_m0_K0_rough.mat')     
    make_rough_steady_state(params,funcs,XH,XU,'steady_states/SS_m0_K0')
end
find_steady_state('steady_states/SS_m0_K0.mat','steady_states/SS_m0_K0_rough.mat',XA,params,funcs);

% Next: find solution for various m with K = 0
 %mrates = [0:0.5:6 6.1:.1:8] ;
mrates  = 0:0.25:8;
for l = 2:length(mrates)
    funcs.meltf = @(x) mrates(l)*ones(size(x));
    find_steady_state(sprintf('steady_states/SS_m%g_K0.mat', mrates(l)), sprintf('steady_states/SS_m%g_K0.mat', mrates(l-1)),XA,params,funcs);
end

funcs.meltf = @(x) 1*ones(size(x));

funcs8 = funcs;
funcs8.meltf = @(x) 8*ones(size(x));

% funcs6 = funcs;
% funcs6.meltf = @(x) 8*ones(size(x));

% Next, find solutions for various K for the cases m=1 and m=8 respectively
Ks =0:0.5:6 ;
for l = 2:length(Ks)
    params.K = Ks(l);  

    % Solution for m=1
    initfilename= sprintf('steady_states/SS_m1_K%g.mat', Ks(l-1));
    if isfile(initfilename)
        find_steady_state(sprintf('steady_states/SS_m1_K%g.mat', Ks(l)), initfilename, XA,params,funcs);
    end

    % Solution for m=8
    initfilename= sprintf('steady_states/SS_m8_K%g.mat', Ks(l-1));
    if isfile(initfilename)
    find_steady_state(sprintf('steady_states/SS_m8_K%g.mat', Ks(l)), initfilename, XA,params,funcs8);
    end    
    % initfilename= sprintf('steady_states/SS_m6_K%g.mat', Ks(l-1));
    % if isfile(initfilename)
    %      find_steady_state(sprintf('steady_states/SS_m6_K%g.mat', Ks(l)), initfilename, XA,params,funcs6);
    % end
end


% Initial solution for alternative geometry when K=0 can be copied over
if ~isfile('steady_states/SS_nonu_m8_K0.mat')
    load('steady_states/SS_m8_K0.mat','HUxN');
    save('steady_states/SS_nonu_m8_K0.mat','HUxN');
end
% if ~isfile('steady_states/SS_nonu_m6_K0.mat')
%      load('steady_states/SS_m6_K0.mat','HUxN');
%      save('steady_states/SS_nonu_m6_K0.mat','HUxN');
% end

% Find solutions for various K, m=8, in the alternative geometry
Ks = 0:0.5:3;
funcs8.Hsbf = Hsbf2;
for l = 2:length(Ks)
    params.K = Ks(l);    
    initfilename= sprintf('steady_states/SS_nonu_m8_K%g.mat', Ks(l-1));
    if isfile(initfilename)
        find_steady_state(sprintf('steady_states/SS_nonu_m8_K%g.mat', Ks(l)), initfilename, XA,params,funcs8);
    end
end
% funcs6.Hsbf = Hsbf2;
% for l = 2:length(Ks)
%     params.K = Ks(l);    
%     initfilename= sprintf('steady_states/SS_nonu_m6_K%g.mat', Ks(l-1));
%     if isfile(initfilename)
%         find_steady_state(sprintf('steady_states/SS_nonu_m6_K%g.mat', Ks(l)), initfilename, XA,params,funcs6);
%     end
% end
%% Plotting results

% Figure 2: Solutions when K=0 for m=1 and 8 respectively 

HUxN1 = load('steady_states/SS_m1_K0.mat','HUxN').HUxN;
HUxN2 = load('steady_states/SS_m8_K0.mat','HUxN').HUxN; 
%HUxN2 = load('steady_states/SS_m6_K0.mat','HUxN').HUxN; 

% Approximate solutions for N
% P1 = HUxN1(1:J+1)+ bedf(HUxN1(2*J+4)*XH)/(1-delta);
% G1 = HUxN1(2*J+4)^(-1)*(P1(2:J+1)-P1(1:J))./(XH(2:J+1)-XH(1:J));
% Nouter1 = -G1./(HUxN1(2*J+4)*1*XU(2:J+1));
% Ninner1 = P1/E;
% 
% P2 = HUxN2(1:J+1)+ bedf(HUxN2(2*J+4)*XH)/(1-delta);
% G2 = HUxN2(2*J+4)^(-1)*(P2(2:J+1)-P2(1:J))./(XH(2:J+1)-XH(1:J));
% Nouter2 = -G2./(HUxN2(2*J+4)*8*XU(2:J+1));
% Ninner2 = P2/E;

P1 = (1-delta)*HUxN1(1:J+1)+ bedf(HUxN1(2*J+4)*XH);
G1 = HUxN1(2*J+4)^(-1)*(P1(2:J+1)-P1(1:J))./(XH(2:J+1)-XH(1:J));
Nouter1 = -G1./(HUxN1(2*J+4)*1*XU(2:J+1));
Ninner1 = P1/E;
%Nnear01 = (1/(E*HUxN1(1)^(mW+1)))*exp(1*(HUxN1(2*J+4)*XH).^2/(2*E)).*cumtrapz(XH*HUxN1(2*J+4),HUxN1(2*J+4)*XH).^mW.*exp(-1*(HUxN1(2*J+4)*XH).^2/(2*E));

% N_exact = (-1/E)* exp(1*(HUxN1(2*J+4)*XU(2:J+1)).^2/(2*E)).*cumsum( G1.*exp(-1*(HUxN1(2*J+4)*XU(2:J+1)).^2/(2*E)).*(XH(2:J+1)-XH(1:J))*HUxN1(2*J+4)  ,'reverse');
% N_int1 = (1/E)* exp(1*(HUxN1(2*J+4)*XU(2:J+1)).^2/(2*E)).*cumsum( (HUxN1(2*J+4)*(XU(2:J+1))).^mW/HUxN1(1).^(mW+1).*exp(-1*(HUxN1(2*J+4)*XU(2:J+1)).^2/(2*E)).*(XH(2:J+1)-XH(1:J))*HUxN1(2*J+4)  ,'reverse');
% N_int1 = HUxN1(1).^(-mW-1)*(1/E)* exp(1*(HUxN1(2*J+4)*XU(2:J+1)).^2/(2*E)).*cumsum( (HUxN1(2*J+4)*(XU(2:J+1))).^mW.*exp(-1*(HUxN1(2*J+4)*XU(2:J+1)).^2/(2*E)).*(XH(2:J+1)-XH(1:J))*HUxN1(2*J+4)  ,'reverse');
% 

Nnear01 = ((2*E)^((mW-1)/2)/(1^((mW+1)/2)*HUxN1(1)^(mW+1)))*exp(1*(HUxN1(2*J+4)*XH).^2/(2*E)).*igamma((mW+1)/2,1*(HUxN1(2*J+4)*XH).^2/(2*E));

P2 = (1-delta)*HUxN2(1:J+1)+ bedf(HUxN2(2*J+4)*XH);
G2 = HUxN2(2*J+4)^(-1)*(P2(2:J+1)-P2(1:J))./(XH(2:J+1)-XH(1:J));
Nouter2 = -G2./(HUxN2(2*J+4)*8*XU(2:J+1));
Ninner2 = P2/E;
Nnear02 = ((2*E)^((mW-1)/2)/(8^((mW+1)/2)*HUxN2(1)^(mW+1)))*exp(8*(HUxN2(2*J+4)*XH).^2/(2*E)).*igamma((mW+1)/2,8*(HUxN2(2*J+4)*XH).^2/(2*E));

% Initialisation of axes
figure(1); clf; 
set(gcf(),'Position',[170   700   1110   640])
a2H1 = make_nice_axes(subplot('Position', [0.08  0.41  0.41  0.515]),'','$z$ / m','$m=1$ mm yr$^{-1}$',xd_km*2);
a2H2 = make_nice_axes(subplot('Position', [0.58  0.41  0.41  0.515]),'','$z$ / m','$m=8$ mm yr$^{-1}$',xd_km*1.8);
a2N1 = make_nice_axes(subplot('Position', [0.08  0.11  0.41  0.235]),'$x$ / km','$N$ / MPa','',xd_km*2);
a2N2 = make_nice_axes(subplot('Position', [0.58  0.11  0.41  0.235]),'$x$ / km','$N$ / MPa','',xd_km*1.8);
set(a2H1,'XTickLabel','');
set(a2H2,'XTickLabel','');

a2H1I = axes('Position', [0.39  0.785  0.09 0.12]); hold on;
a2H2I = axes('Position', [0.89  0.785  0.09  0.12]); hold on;
% a2N1I = axes('Position', [0.39  0.22  0.09 0.1]); hold on;
% a2N2I = axes('Position', [0.89  0.22  0.09  0.1]); hold on;
set(a2H1I,'TickLabelInterpreter','latex')
set(a2H2I,'TickLabelInterpreter','latex')
% set(a2N1I,'TickLabelInterpreter','latex')
% set(a2N2I,'TickLabelInterpreter','latex')

% Plot sea level
plot(a2H1,xd_km*[0 2], [0 0], 'LineWidth',1, 'Color', [0.5 0.5 0.5]);
plot(a2H2,xd_km*[0 2], [0 0], 'LineWidth',1, 'Color', [0.5 0.5 0.5]);

% Bed, sedimentary basin, ice sheet, grounding line
plot(a2H1, xd_km*HUxN1(2*J+4)*[XH; 2], Hd*bedf(HUxN1(2*J+4)*[XH; 2]),'LineWidth',2,'Color','k');
p2b = plot(a2H2, xd_km*HUxN2(2*J+4)*[XH; 2], Hd*bedf(HUxN2(2*J+4)*[XH; 2]),'LineWidth',2,'Color','k');
plot(a2H1, xd_km*HUxN1(2*J+4)*[XH; 2], Hd*(bedf(HUxN1(2*J+4)*[XH; 2])-Hsbf(HUxN1(2*J+4)*[XH; 2])),'LineWidth',2,'Color','k','Linestyle','--');
p2sb = plot(a2H2, xd_km*HUxN2(2*J+4)*[XH; 2], Hd*(bedf(HUxN2(2*J+4)*[XH; 2])-Hsbf(HUxN1(2*J+4)*[XH; 2])),'LineWidth',2,'Color','k','Linestyle','--');
plot(a2H1, xd_km*HUxN1(2*J+4)*XH, Hd*(HUxN1(1:J+1)+ bedf(HUxN1(2*J+4)*XH)) ,'LineWidth',2,'Color',clrs(6,:));
p2h = plot(a2H2, xd_km*HUxN2(2*J+4)*XH, Hd*(HUxN2(1:J+1)+ bedf(HUxN2(2*J+4)*XH)) ,'LineWidth',2,'Color',clrs(6,:));
plot(a2H1, xd_km*HUxN1(2*J+4)*[1 1], Hd*bedf(HUxN1(2*J+4))*[1 -delta/(1-delta)] ,':','LineWidth',2,'Color',clrs(6,:));
plot(a2H2, xd_km*HUxN2(2*J+4)*[1 1], Hd*bedf(HUxN2(2*J+4))*[1 -delta/(1-delta)] ,':','LineWidth',2,'Color',clrs(6,:));

plot(a2H1I, xd_km*HUxN1(2*J+4)*XH, Hd*(HUxN1(1:J+1)+ bedf(HUxN1(2*J+4)*XH)) ,'LineWidth',2,'Color',clrs(6,:));
plot(a2H2I, xd_km*HUxN2(2*J+4)*XH, Hd*(HUxN2(1:J+1)+ bedf(HUxN2(2*J+4)*XH)) ,'LineWidth',2,'Color',clrs(6,:));
plot(a2H1, xd_km*HUxN1(2*J+4)*[1 1], Hd*bedf(HUxN1(2*J+4))*[1 0] ,':','LineWidth',2,'Color',clrs(6,:));
plot(a2H2, xd_km*HUxN2(2*J+4)*[1 1], Hd*bedf(HUxN2(2*J+4))*[1 0] ,':','LineWidth',2,'Color',clrs(6,:));

xlim(a2H1I,[870 895])
xlim(a2H2I,[775 800])

% Effective pressure
plot(a2N1, xd_km*HUxN1(2*J+4)*XH, Nd_MPa*HUxN1(2*J+5:3*J+5),'LineWidth',2,'Color',clrs2(6,:));
plot(a2N2, xd_km*HUxN2(2*J+4)*XH, Nd_MPa*HUxN2(2*J+5:3*J+5),'LineWidth',2,'Color',clrs2(6,:));

% Approximations to effective pressure
plot(a2N1,xd_km*HUxN1(2*J+4)*XU(2:J+1), Nd_MPa*Nouter1, '--' ,'Color', clrs2(4,:) , 'LineWidth',2);
plot(a2N1,xd_km*HUxN1(2*J+4)*XH, Nd_MPa*Ninner1, '-.' ,'Color', clrs2(2,:) , 'LineWidth',1.5);
plot(a2N2,xd_km*HUxN2(2*J+4)*XU(2:J+1), Nd_MPa*Nouter2, '--' ,'Color', clrs2(4,:) , 'LineWidth',2);
plot(a2N2,xd_km*HUxN2(2*J+4)*XH, Nd_MPa*Ninner2, '-.' ,'Color', clrs2(2,:) , 'LineWidth',1.5);

plot(a2N1,xd_km*HUxN1(2*J+4)*XH(1:40), Nd_MPa*Nnear01(1:40), '--' ,'Color', clrs2(1,:) , 'LineWidth',1.5);
plot(a2N2,xd_km*HUxN2(2*J+4)*XH(1:40), Nd_MPa*Nnear02(1:40), '--' ,'Color', clrs2(1,:) , 'LineWidth',1.5);

% plot(a2N1I, xd_km*HUxN1(2*J+4)*XH, Nd_MPa*HUxN1(2*J+5:3*J+5),'LineWidth',2,'Color',clrs2(6,:));
% plot(a2N2I, xd_km*HUxN2(2*J+4)*XH, Nd_MPa*HUxN2(2*J+5:3*J+5),'LineWidth',2,'Color',clrs2(6,:));
% plot(a2N1I,xd_km*HUxN1(2*J+4)*XU(2:J+1), Nd_MPa*Nouter1, '--' ,'Color', clrs2(4,:) , 'LineWidth',2);
% plot(a2N1I,xd_km*HUxN1(2*J+4)*XH, Nd_MPa*Ninner1, '-.' ,'Color', clrs2(2,:) , 'LineWidth',1.5);
% plot(a2N2I,xd_km*HUxN2(2*J+4)*XU(2:J+1), Nd_MPa*Nouter2, '--' ,'Color', clrs2(4,:) , 'LineWidth',2);
% plot(a2N2I,xd_km*HUxN2(2*J+4)*XH, Nd_MPa*Ninner2, '-.' ,'Color', clrs2(2,:) , 'LineWidth',1.5);
% 
% xlim(a2N1I,[870 895])
% xlim(a2N2I,[775 800])

% Axis limits
ylim(a2H1, [-3000 4000]);
ylim(a2H2, [-3000 4000]);
ylim(a2N1, [0 2.25]);
ylim(a2N2, [0 0.75]);

% Generating streamfunction for streamlines
z = (-3000:10:0)/Hd;
[Psi1,xx1,zz1] = make_SF(HUxN1(1:J+1),HUxN1(2*J+5:3*J+5),bedf(HUxN1(2*J+4)*XH),Hsbf(HUxN1(2*J+4)*XH),HUxN1(2*J+4)*XH,z,delta,E);
[Psi2,xx2,zz2] = make_SF(HUxN2(1:J+1),HUxN2(2*J+5:3*J+5),bedf(HUxN2(2*J+4)*XH),Hsbf(HUxN2(2*J+4)*XH),HUxN2(2*J+4)*XH,z,delta,E);

% Streamlines
contour(a2H1,xd_km*xx1', Hd*zz1', Psi1, (linspace(0.015^.2,max(Psi1,[],'all')^.2,10)).^5, 'Color','r','LineWidth',.5);
contour(a2H2, xd_km*xx2', Hd*zz2', Psi2,linspace(0.015^.2,(max(Psi2,[],'all')-0.01)^.2,10).^5, 'Color','r','LineWidth',.5);
p2sl = a2H2.Children(1);

% Legends
legend(a2H2,[p2b,p2sb,p2h,p2sl],{'$b$','$b-H_{sb}$','$H+b$','Streamlines'},'location','northwest','FontSize',14,'Interpreter','latex','NumColumns',2);
% legend(a2N2,{'$N$','$E=0$','$\psi =0$'},'location','northeast','FontSize',14,'Interpreter','latex','NumColumns',2);
legend(a2N2,{'$N$','$E=0$','$\psi =0$','Small $x$'},'location','northeast','FontSize',14,'Interpreter','latex','NumColumns',2);

% Subplot labels
text(a2H1,20,3000,'(a)(i)','Interpreter','latex','BackgroundColor','w','FontSize',14);
text(a2N1,20,2,'(a)(ii)','Interpreter','latex','BackgroundColor','w','FontSize',14);
text(a2H2,20,1900,'(b)(i)','Interpreter','latex','BackgroundColor','w','FontSize',14);
text(a2N2,50,0.65,'(b)(ii)','Interpreter','latex','BackgroundColor','w','FontSize',14);

% Export
exportgraphics(figure(1),'figures/SL_K0.pdf','ContentType','vector')

% Figure 4: Solution for m=8 and varying K.

% Initialising axes
figure(2); clf; %set(gcf(),'Position',[ 250   400   630   735]) 
set(gcf(),'Position',[250   424   800  700])
a4q = make_nice_axes(subplot('Position', [.13 .42 .8 .15]),'','$q_E$ / mm yr$^{-1}$','',xd_km*1.8);
a4H = make_nice_axes(subplot('Position', [.13 .61 .8 .25]),'','$z$ / m','',xd_km*1.8);
a4U = make_nice_axes(subplot('Position', [.13 .08 .8 .15]),'$x$ / km','$u$ / m yr$^{-1}$','',xd_km*1.8);
a4N = make_nice_axes(subplot('Position', [.13 .25 .8 .15]),'','$N$ / MPa','',xd_km*1.8);
set(a4H,'XTickLabel','');
set(a4N,'XTickLabel','','Ylim',[0 3]);
set(a4q,'XTickLabel','','YLim',[-100 100]);
set(a4U,'YLim',[0 450]);

load('steady_states/SS_m8_K0.mat','HUxN');
%load('steady_states/SS_m6_K0.mat','HUxN');
% Sea level, bed, sedimentary basin
plot(a4H,xd_km*[0 2], [0 0], 'LineWidth',1, 'Color', [0.5 0.5 0.5]);
plot(a4H, xd_km*HUxN(2*J+4)*[XH; 2], Hd*bedf(HUxN(2*J+4)*[XH; 2]),'k','LineWidth',2);
plot(a4H, xd_km*HUxN(2*J+4)*[XH; 2], Hd*(bedf(HUxN(2*J+4)*[XH; 2])-Hsbf(HUxN(2*J+4)*[XH; 2])),'k--','LineWidth',2);
% Initial streamlines
contour(a4H, xd_km*xx2', Hd*zz2', Psi2,linspace(0.015^.2,max(Psi2,[],'all')^.2,12).^5, 'Color',clrs2(1,:),'LineWidth',1);

% Plotting each variable for various K
for K = 0:2:6
    load(sprintf('steady_states/SS_m8_K%g.mat',K),'HUxN');
    %load(sprintf('steady_states/SS_m6_K%g.mat',K),'HUxN');
    qE = get_qE(HUxN(1:J+1),HUxN(2*J+5:3*J+5),bedf(HUxN(2*J+4)*XH),Hsbf(HUxN(2*J+4)*XU(2:J+1)),HUxN(2*J+4)*XH,K,delta,E);
    plot(a4q, xd_km*HUxN(2*J+4)*XH(2:J), md_mmyr*qE,'linewidth',2,'Color',clrs2(K+1,:));
    plot(a4H, xd_km*HUxN(2*J+4)*XH, Hd*(HUxN(1:J+1)+ bedf(HUxN(2*J+4)*XH)),'linewidth',2,'Color',clrs2(K+1,:));
    plot(a4H, xd_km*HUxN(2*J+4)*[1 1], Hd*bedf(HUxN(2*J+4))*[1 -delta/(1-delta)] ,':','LineWidth',2,'Color',clrs2(K+1,:));
    plot(a4U, xd_km*HUxN(2*J+4)*XU, ud_myr*HUxN(J+2:2*J+3),'LineWidth',2,'Color',clrs2(K+1,:));
    plot(a4N, xd_km*HUxN(2*J+4)*XH, Nd_MPa*HUxN(2*J+5:3*J+5),'LineWidth',2,'Color',clrs2(K+1,:));
end

% Streamlines for largest value of K
[Psi3,xx3,zz3] = make_SF(HUxN(1:J+1),HUxN(2*J+5:3*J+5),bedf(HUxN(2*J+4)*XH),Hsbf(HUxN(2*J+4)*XH),HUxN(2*J+4)*XH,z,delta,E);
contour(a4H, xd_km*xx3', Hd*zz3', Psi3,linspace(0.005^.2,max(Psi3,[],'all')^.2,12).^5, 'Color',clrs2(7,:),'LineWidth',1);

% Legend
% legend(a4q,{'$k_{sb} = 0$ m$^2$ $(K \ll 1)$','$k_{sb} = 3.3 \times 10^{-12}$ m$^2$',...
%     '$k_{sb} = 6.6 \times 10^{-12}$ m$^2$','$k_{sb} = 9.8 \times 10^{-12}$ m$^2$'},...
%     'Position', [.13 .88 .8 .1],'FontSize',14,'Interpreter','latex','Numcolumns',2);
legend(a4q,['$k_{sb} = 0$ m$^2$ $(K \ll 1)$'; strcat(legend_text('k_{sb}', scales.K_to_ksb*(2:2:6),'m$^2$'),' (',legend_text('K', 2:2:6,''),')')],...
     'Position', [.13 .88 .8 .1],'FontSize',14,'Interpreter','latex','Numcolumns',2);

% Subplot labels
text(a4H,20,1800,'(a)','Interpreter','latex','BackgroundColor','w','FontSize',14);
text(a4q,20,70,'(b)','Interpreter','latex','BackgroundColor','w','FontSize',14);
text(a4N,20,2.4,'(c)','Interpreter','latex','BackgroundColor','w','FontSize',14);
text(a4U,20,375,'(d)','Interpreter','latex','BackgroundColor','w','FontSize',14);

% Aligning axis labels
a4H.YLabel.Position(1) = -80;
a4N.YLabel.Position(1) = -80;
a4q.YLabel.Position(1) = -80;
a4U.YLabel.Position(1) = -80;

% Exporting figure
exportgraphics(figure(2),'figures/SS_varyK_m8.pdf','ContentType','vector');

% Figure 3: Solution for m=1 and varying K.

% Initialising axes
figure(3); clf;
set(gcf(),'Position',[ 250   400   600   600])
a3H = make_nice_axes(subplot('Position',[0.13 0.45 0.85 0.4]),'','$z$ / m','',xd_km*1.9);
a3N = make_nice_axes(subplot('Position',[0.13 0.1 0.85 0.3]),'$x$ / km','$N$ / MPa','',xd_km*1.9);
set(a3H,'XTickLabel','','YLim',[-3000 4000]);
set(a3N,'YLim',[0 30]);

% Sea level, bed, sedimentary basin
load('steady_states/SS_m1_K0.mat','HUxN');   
plot(a3H,xd_km*[0 2], [0 0], 'LineWidth',1, 'Color', [0.5 0.5 0.5]);
plot(a3H, xd_km*HUxN(2*J+4)*[XH; 2], Hd*bedf(HUxN(2*J+4)*[XH; 2]),'k','LineWidth',2);
plot(a3H, xd_km*HUxN(2*J+4)*[XH; 2], Hd*(bedf(HUxN(2*J+4)*[XH; 2])-Hsbf(HUxN(2*J+4)*[XH; 2])),'k--','LineWidth',2);
% Initial streamlines
contour(a3H, xd_km*xx1', Hd*zz1', Psi1,linspace(0.015^.2,max(Psi1,[],'all')^.2,12).^5, 'Color',clrs2(1,:),'LineWidth',1);
% Plotting variables for each K
for K = 0:2:6
    load(sprintf('steady_states/SS_m1_K%g.mat',K),'HUxN');
    plot(a3H, xd_km*HUxN(2*J+4)*XH, Hd*(HUxN(1:J+1)+ bedf(HUxN(2*J+4)*XH)),'linewidth',2,'Color',clrs2(K+1,:));
    plot(a3H, xd_km*HUxN(2*J+4)*[1 1], Hd*bedf(HUxN(2*J+4))*[1 -delta/(1-delta)] ,':','LineWidth',2,'Color',clrs2(K+1,:));
    plot(a3N, xd_km*HUxN(2*J+4)*XH, Nd_MPa*HUxN(2*J+5:3*J+5),'LineWidth',2,'Color',clrs2(K+1,:));
end

% Streamlines for largest value of K
[Psi4,xx4,zz4] = make_SF(HUxN(1:J+1),HUxN(2*J+5:3*J+5),bedf(HUxN(2*J+4)*XH),Hsbf(HUxN(2*J+4)*XH),HUxN(2*J+4)*XH,z,delta,E);
contour(a3H, xd_km*xx4', Hd*zz4', Psi4,linspace(0.005^.2,max(Psi4,[],'all')^.2,12).^5, 'Color',clrs2(7,:),'LineWidth',1);

% Legend
% legend(a3H,flip(a3N.Children),{'$k_{sb} = 0$ m$^2$ $(K\ll 1)$','$k_{sb} = 3.3 \times 10^{-12}$ m$^2$',...
%     '$k_{sb} = 6.6 \times 10^{-12}$ m$^2$','$k_{sb} = 9.8 \times 10^{-12}$ m$^2$'},...
%     'position',[0.13 0.9 0.85 0.08],'FontSize',14,'Interpreter','latex','Numcolumns',2);
legend(a3H,flip(a3N.Children),['$k_{sb} = 0$ m$^2$ $(K \ll 1)$'; legend_text('k_{sb}', scales.K_to_ksb*(2:2:6),'m$^2$')],...
     'position',[0.13 0.9 0.85 0.08],'FontSize',14,'Interpreter','latex','Numcolumns',2); 

% Subplot labels
text(a3H,875,3500,'(a)','Interpreter','latex','BackgroundColor','w','FontSize',14);
text(a3N,875,25,'(b)','Interpreter','latex','BackgroundColor','w','FontSize',14);

% Export
exportgraphics(figure(3),'figures/SS_varyK_m1.pdf','ContentType','vector');

% Figure 5: alternative sedimentary basin geometry
% Initialise axes
figure(4); clf; 
set(gcf(),'Position',[ 250   400   600   500])
a5H = make_nice_axes(subplot('Position',[0.13 0.45 0.85 0.4]),'','$z$ / m','',xd_km*1.7);
a5N = make_nice_axes(subplot('Position',[0.13 0.1 0.85 0.3]),'$x$ / km','$N$ / MPa','',xd_km*1.7);
set(a5H,'XTickLabel','','YLim',[-2000 2750]);
set(a5N,'YLim',[0 1]);

% Sea level, bed, sedimentary basin
load('steady_states/SS_nonu_m8_K0.mat','HUxN'); 
%load('steady_states/SS_nonu_m6_K0.mat','HUxN'); 
plot(a5H,xd_km*[0 2], [0 0], 'LineWidth',1, 'Color', [0.5 0.5 0.5]);
plot(a5H, xd_km*HUxN(2*J+4)*[XH; 2], Hd*bedf(HUxN(2*J+4)*[XH; 2]),'k','LineWidth',2);
plot(a5H, xd_km*(0:0.01:2), Hd*(bedf((0:0.01:2))-Hsbf2(0:0.01:2)),'k--','LineWidth',2);
% Initial streamlines
[Psi6,xx6,zz6] = make_SF(HUxN(1:J+1),HUxN(2*J+5:3*J+5),bedf(HUxN(2*J+4)*XH),Hsbf2(HUxN(2*J+4)*XH),HUxN(2*J+4)*XH,z,delta,E);
contour(a5H, xd_km*xx6', Hd*zz6', Psi6,linspace(0.015^.2,max(Psi6,[],'all')^.2,12).^5, 'Color',clrs2(1,:),'LineWidth',1);
% Plot variables for various K
for K = 0:3
    load(sprintf('steady_states/SS_nonu_m8_K%g.mat',K),'HUxN');
    %load(sprintf('steady_states/SS_nonu_m6_K%g.mat',K),'HUxN');
    plot(a5H, xd_km*HUxN(2*J+4)*XH, Hd*(HUxN(1:J+1)+ bedf(HUxN(2*J+4)*XH)),'linewidth',2,'Color',clrs2(2*K+1,:));
    plot(a5H, xd_km*HUxN(2*J+4)*[1 1], Hd*bedf(HUxN(2*J+4))*[1 -delta/(1-delta)] ,':','LineWidth',2,'Color',clrs2(K+1,:));
    plot(a5N, xd_km*HUxN(2*J+4)*XH, Nd_MPa*HUxN(2*J+5:3*J+5),'LineWidth',2,'Color',clrs2(2*K+1,:));
end
% Streamlines for final K
[Psi5,xx5,zz5] = make_SF(HUxN(1:J+1),HUxN(2*J+5:3*J+5),bedf(HUxN(2*J+4)*XH),Hsbf2(HUxN(2*J+4)*XH),HUxN(2*J+4)*XH,z,delta,E);
contour(a5H, xd_km*xx5', Hd*zz5', Psi5,linspace(0.015^.2,max(Psi5,[],'all')^.2,12).^5, 'Color',clrs2(7,:),'LineWidth',1);

% Legend
% legend(a5H,flip(a5N.Children),{'$k_{sb} = 0$ m$^2$ $(K\ll 1)$','$k_{sb} = 3.3 \times 10^{-12}$ m$^2$',...
%     '$k_{sb} = 6.6 \times 10^{-12}$ m$^2$','$k_{sb} = 9.8 \times 10^{-12}$ m$^2$'},...
%     'position',[0.13 0.9 0.85 0.08],'FontSize',13,'Interpreter','latex','Numcolumns',2);
legend(a5H,flip(a5N.Children),['$k_{sb} = 0$ m$^2$ $(K \ll 1)$'; legend_text('k_{sb}', scales.K_to_ksb*(1:3),'m$^2$')],...
     'position',[0.13 0.9 0.85 0.08],'FontSize',13,'Interpreter','latex','Numcolumns',2);


% Subplot labels
text(a5H,25,2500,'(a)','Interpreter','latex','BackgroundColor','w','FontSize',14);
text(a5N,25,0.8,'(b)','Interpreter','latex','BackgroundColor','w','FontSize',14);

% Export
exportgraphics(figure(4),'figures/SS_nonu.pdf','ContentType','vector');

% A function to find the streamfunction of groundwater flow
function [Psi,xx,zz] = make_SF(H,N,b,Hsb,x,z,delta,E)  
    [xx,zz] = meshgrid(x,z);    
    
    %P = H - E*N + b/(1-delta);
    P = (1-delta)*H - E*N + b;
    Usb = - gradient(P,x);
   
    Psi = (z - b + Hsb).*Usb;
   
    % Mask values outside of sedimntary basin
    Psi(zz'>b) = nan;
    Psi(zz'<b-Hsb) = nan;
end

% A function to calculate the exfiltration (in steady state)
function qE = get_qE(H,N,b,Hsb,x,K,delta,E)
    %Psi = H - E*N + b/(1-delta) ;
    Psi = (1-delta)*H - E*N + b ;
    flux =  K*Hsb.* (Psi(2:end)-Psi(1:end-1))./(x(2:end)-x(1:end-1));

    qE = (flux(2:end)-flux(1:end-1))./(x(2:end-1)-x(1:end-2));
end
