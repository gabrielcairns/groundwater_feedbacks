% Plots the results of experiments involving grounding line advance and
% retreat (Figures 6, A1 and A2)

% Scales, parameters and functions
[scales, params, funcs] = default_values();
Hd = scales.Hd;     xd_km = scales.xd_km;     md_mmyr = scales.md_mmyr;
td_kyr = scales.td_kyr;     Nd_MPa = scales.Nd_MPa;
delta = params.delta;   E = params.E; xi = params.xi;
ud_myr = scales.ud_myr;

bedf = funcs.bedf;
bedxf = funcs.bedxf;
Hsbf = funcs.Hsbf;

% Meshgrid
load("variable_grid_300.mat","XA","XH","XU");
J = 300;

% Colormaps for plots
clrsr = cmocean('thermal',8);
clrsb = cmocean('haline',7);

% Folder to store output figures
if ~isfolder('figures');   mkdir('figures');  end

%% Running experiments
SigmasAR = [0 60 120 240 480];

retreat_or_advance('r',1,1,0,500);
retreat_or_advance('a',1,1,0,1000);

for Sigma = SigmasAR(2:end)
    retreat_or_advance('r',1,1,Sigma);
    retreat_or_advance('a',1,1,Sigma);
end

SigmaP = 30;

% Retreat for Sigma = 0,30 with m=1 and various K
Ks = 0.5:0.5:2.5;
for K = Ks
    retreat_or_advance('r',1,K,0);
    retreat_or_advance('r',1,K,SigmaP);
end

% Retreat for Sigma = 0,30 with K=1 and various constant m
    retreat_or_advance('r',0,1,0,800);
    retreat_or_advance('r',0,1,SigmaP);
ms = 0:0.5:2;
for mrate = ms(2:end)
    retreat_or_advance('r',mrate,1,0);
    retreat_or_advance('r',mrate,1,SigmaP);
end

%% Plotting results 

% Additional figure showing ice flux as a function of position in both
% cases
figure(4); clf;
set(gcf(),'Position',[ 1 1 915 580])
aqI = make_nice_axes(subplot(2,1,1), '$x_g$ / km', '$q_g$ / m kyr$^{-1}$', 'Retreat');
set ( aqI, 'XDir', 'reverse' )
aqIa = make_nice_axes(subplot(2,1,2), '$x_g$ / km', '$q_g$ / m kyr$^{-1}$', 'Advance');

% Figure 6: retreat for m=K=1 and various Sigma

% Extracting the data
[H1, N1, xg1, ~, t, ks,qi1, U1] = get_data('r', 1,1,SigmasAR(1),  Hd,XH,J,E,xi,delta);
[~, ~, xg2, ~, ~, ~,qi2] = get_data('r', 1,1,SigmasAR(2),  Hd,XH,J,E,xi,delta);
[~, ~, xg3, ~, ~, ~,qi3] = get_data('r', 1,1,SigmasAR(3),  Hd,XH,J,E,xi,delta);
[~, ~, xg4, ~, ~, ~,qi4] = get_data('r', 1,1,SigmasAR(4),  Hd,XH,J,E,xi,delta);
[H5, N5, xg5, qE5, ~, ~,qi5, U5] = get_data('r', 1,1,SigmasAR(5),  Hd,XH,J,E,xi,delta);

% Initialising axes for xg against t
figure(1); clf;
set(gcf(),'Position',[400 400 1000 1130],'Color','w');
axg = make_nice_axes(subplot('Position',[0.15 0.77 0.8 0.22]), '$t$ / kyr', '$x_g$ / km', '');
xlim(axg,[0 31]); ylim(axg,[700 900])

% Plotting xg against t for different cases in Sigma
plot(axg, td_kyr*t, xd_km*xg1, 'color', clrsr(6,:), 'LineWidth', 1.5);
plot(axg, td_kyr*t, xd_km*xg2, 'color', clrsr(5,:) ,'LineWidth', 1.5);
plot(axg, td_kyr*t, xd_km*xg3, 'color', clrsr(4,:) ,'LineWidth', 1.5);
plot(axg, td_kyr*t, xd_km*xg4, 'color', clrsr(3,:) ,'LineWidth', 1.5);
plot(axg, td_kyr*t, xd_km*xg5, 'color', clrsr(2,:), 'LineWidth', 1.5);

% Ice flux (additional figure)
plot(aqI, xd_km*xg1, Hd*ud_myr*1e-3*qi1, 'color', clrsr(6,:), 'LineWidth', 1.5);
plot(aqI, xd_km*xg2, Hd*ud_myr*1e-3*qi2, 'color', clrsr(5,:) ,'LineWidth', 1.5);
plot(aqI, xd_km*xg3, Hd*ud_myr*1e-3*qi3, 'color', clrsr(4,:) ,'LineWidth', 1.5);
plot(aqI, xd_km*xg4, Hd*ud_myr*1e-3*qi4, 'color', clrsr(3,:) ,'LineWidth', 1.5);
plot(aqI, xd_km*xg5, Hd*ud_myr*1e-3*qi5, 'color', clrsr(2,:), 'LineWidth', 1.5);

% Initalising the rest of the axes (snapshots of H, N, qE)
ah = make_nice_axes(subplot('Position',[0.15 0.55 0.8 0.17]), '', '$z$ / m', '',1000);
an = make_nice_axes(subplot('Position',[0.15 0.21 0.8 0.17]), '', '$N$ / MPa', '',1000);
aq = make_nice_axes(subplot('Position',[0.15 0.42 0.8 0.09]), '', '$q_{Ec}$ / mm yr$^{-1}$', '',1000);
aflux = make_nice_axes(subplot('Position',[0.15 0.08 0.8 0.09]), '$x$ / km', '$H u $ / km$^2$ yr$^{-1}$', '',1000);
set(ah, 'YLim', [-1000 4000], 'XTickLabel', '');
set(an, 'YLim', [0 6], 'XTickLabel', '');
set(aflux,  'YLim', [0 .3]);
set(aq, 'YLim', [-120 10], 'XTickLabel', '');

% Bed, sea level
plot(ah, xd_km*2*XH, Hd*bedf(2*XH), 'k-', 'LineWidth', 2)
plot(ah, xd_km*[0 2], [0 0], 'LineWidth',1, 'Color', [0.5 0.5 0.5]);

% The entries at which we plot snapshots
ksp = [142 200 284 399];

% Plotting initial condition
load('retreat_advance/max_xg_m1_K1.mat','HUxN');
plot(ah, xd_km*HUxN(2*J+4)*XH, Hd*(bedf(HUxN(2*J+4)*XH)+ HUxN(1:J+1)), 'color', clrsb(1,:),'linestyle', '-', 'LineWidth', 2);
plot(an, xd_km*HUxN(2*J+4)*XH, Nd_MPa*HUxN(2*J+5:3*J+5), 'color',  clrsb(1,:),'linestyle','-', 'LineWidth', 2);
plot(ah, xd_km*HUxN(2*J+4)*[1 1], Hd*bedf(HUxN(2*J+4))*[1 -delta/(1-delta)] ,':','LineWidth',2,'Color',clrsb(1,:));

% Initial states
plot(aq, xd_km*xg1(1)*XH(1:J), md_mmyr*qE5(:,1),  'Color', clrsb(1,:)  , 'LineWidth', 2);  
plot(aflux, xd_km*xg1(1)*XH, 1e-6*Hd*ud_myr*H1(:,1).*U1(2:end,1),  'Color', clrsb(1,:)  , 'LineWidth', 1); 
plot(aflux, xd_km*xg5(1)*XH, 1e-6*Hd*ud_myr*H5(:,1).*U5(2:end,1),  'Color', clrsb(1,:)  , 'LineWidth', 2); 

ains = axes('Position',[0.86 0.3 0.09 0.08]); hold on;
set(ains,'TickLabelInterpreter','latex');
xlabel(ains,'$x-x_g$','FontSize',12,'Interpreter','latex')
xlim(ains,[-100 0]);

% Plotting snapshots of H, N, qE at each desired time
L = length(ksp); 
for l = 1:L
    XG1 = xg1(ks(ksp(l)));
    XG5 = xg5(ks(ksp(l)));
    plot(ah, xd_km*XG1*XH, Hd*(bedf(XG1*XH)+ H1(:,ksp(l))), 'Color', clrsb(l+1,:) , 'LineWidth', 1, 'Linestyle','-' );
    plot(an, xd_km*XG1*XH, Nd_MPa*N1(:,ksp(l)), 'Color', clrsb(l+1,:)  , 'LineWidth', 1,'Linestyle','-');   
    plot(ah, xd_km*XG5*XH, Hd*(bedf(XG5*XH)+ H5(:,ksp(l))), 'Color', clrsb(l+1,:)   , 'LineWidth', 2);
    plot(an, xd_km*XG5*XH, Nd_MPa*N5(:,ksp(l)), 'Color', clrsb(l+1,:)  , 'LineWidth', 2);
    plot(aq, xd_km*XG5*XH(1:J), md_mmyr*qE5(:,ksp(l)-1),  'Color', clrsb(l+1,:)  , 'LineWidth', 2);  

    plot(aflux, xd_km*XG1*XH, 1e-6*Hd*ud_myr*H1(:,ksp(l)-1).*U1(2:end,ksp(l)-1),  'Color', clrsb(l+1,:)  , 'LineWidth', 1); 
    plot(aflux, xd_km*XG5*XH, 1e-6*Hd*ud_myr*H5(:,ksp(l)-1).*U5(2:end,ksp(l)-1),  'Color', clrsb(l+1,:)  , 'LineWidth', 2); 

    plot(ah, xd_km*XG1*[1 1], Hd*bedf(XG1)*[1 -delta/(1-delta)] ,':','LineWidth',1,'Color',clrsb(l+1,:));
    plot(ah, xd_km*XG5*[1 1], Hd*bedf(XG5)*[1 -delta/(1-delta)] ,':','LineWidth',2,'Color',clrsb(l+1,:));

    plot(ains, xd_km*XG1*(XH-1), Nd_MPa*N1(:,ksp(l)), 'Color', clrsb(l+1,:)  , 'LineWidth', .5,'Linestyle','-');   
    plot(ains, xd_km*XG5*(XH-1), Nd_MPa*N5(:,ksp(l)), 'Color', clrsb(l+1,:)  , 'LineWidth', 2);

end

% Plotting new steady state
load('retreat_advance/min_xg_m1_K1.mat','HUxN');
plot(axg, td_kyr*t([1 end]), xd_km*HUxN(2*J+4)*[1 1], 'color', [0.5 0.5 0.5], 'linestyle', '--', 'LineWidth', 2);
plot(ah, xd_km*HUxN(2*J+4)*XH, Hd*(bedf(HUxN(2*J+4)*XH)+ HUxN(1:J+1)), 'color', [0.5 0.5 0.5],'linestyle','--', 'LineWidth', 2);
plot(an, xd_km*HUxN(2*J+4)*XH, Nd_MPa*HUxN(2*J+5:3*J+5), 'color', [0.5 0.5 0.5], 'linestyle', '--', 'LineWidth', 2);
plot(ah, xd_km*HUxN(2*J+4)*[1 1], Hd*bedf(HUxN(2*J+4))*[1 -delta/(1-delta)] ,'--','LineWidth',2,'Color',[0.5 0.5 0.5]);

% Legends

legend(axg,[legend_text('S_s', SigmasAR*scales.Sigma_to_Ss, ' m$^{-1}$' ) ; 'New steady state'],...
     'FontSize',13,'Interpreter','latex','Position',[0.64 0.84 0.31 0.15])

legend(aqI,legend_text('S_s', SigmasAR*scales.Sigma_to_Ss, ' m$^{-1}$' ),...
     'FontSize',13,'Interpreter','latex','Location','northeast')

legend(ah,[ah.Children(2), ah.Children(4*L+2), ah.Children(4*L+1)],...
    ['New steady state' ; legend_text('S_s', [0 max(SigmasAR)]*scales.Sigma_to_Ss, ' m$^{-1}$' ) ],...
    'FontSize',13,'Interpreter','latex','Position',[0.16 0.555 0.31 .09]);

% Subplot labels
text(axg,0.4,750,'(a)','Interpreter','latex','BackgroundColor','w','FontSize',14);
text(ah,950,3500,'(b)','Interpreter','latex','BackgroundColor','w','FontSize',14);
text(an,950,1.5,'(d)','Interpreter','latex','BackgroundColor','w','FontSize',14);
text(aq,950,0,'(c)','Interpreter','latex','BackgroundColor','w','FontSize',14);
text(aflux,950,.25,'(e)','Interpreter','latex','BackgroundColor','w','FontSize',14);

% Aligning axis labels
ah.YLabel.Position(1) = -96;
aq.YLabel.Position(1) = -96;
an.YLabel.Position(1) = -96;

% Export
exportgraphics(figure(1),'figures/retreat.pdf','ContentType','vector')

% Figure A2: Effect of varying m and K on the difference between Sigma = 0
% and 2

% Initalising axes
figure(2); clf; set(gcf(),'Position',[400 400 950 400],'Color','w');
am = make_nice_axes(subplot('Position',[0.1 0.425 0.375 0.55]), '$t$ / kyr', '$\Delta x_g$ / km', '');
aK = make_nice_axes(subplot('Position',[0.55 0.425 0.375 0.55]), '$t$ / kyr', '', '');
ylim(am,[-15 15]);       ylim(aK,[-6 6]);

% For each m, extract and plot the data
for l=1:length(ms)
    [~, ~, xg1, ~, t, ~] = get_data('r', ms(l),1,0,  Hd,XH,J,E,xi,delta);
    [~, ~, xg2, ~, ~, ~] = get_data('r', ms(l),1,SigmaP,  Hd,XH,J,E,xi,delta);
    plot(am, td_kyr*t, xd_km*(xg2-xg1), 'LineWidth',2, 'Color', clrsr(l,:))
end

% For each K, extract and plot the data
for l=1:length(Ks)
    [~, ~, xg1, ~, ~, ~] = get_data('r', 1, Ks(l), 0,  Hd,XH,J,E,xi,delta);
    [~, ~, xg2, ~, ~, ~] = get_data('r', 1, Ks(l), SigmaP,  Hd,XH,J,E,xi,delta);
    plot(aK, td_kyr*t, xd_km*(xg2-xg1), 'LineWidth',2, 'Color', clrsb(7-l,:))
end

% Legends
legend(am,legend_text('m', md_mmyr*ms,'mm yr$^{-1}$'),...
    'FontSize',13,'Interpreter','latex','Position',[0.1 0.08 0.375 0.2],'numcolumns',2)
legend(aK,legend_text('k_{sb}', scales.K_to_ksb*Ks,'m$^2$'),...
     'FontSize',13,'Interpreter','latex','Position',[0.5 0.08 0.425 0.2],'numcolumns',2)

% Subplot labels
text(am,28,10,'(a)','Interpreter','latex','BackgroundColor','w','FontSize',14);
text(aK,28,5,'(b)','Interpreter','latex','BackgroundColor','w','FontSize',14);

text(am,15,12,'Accelerated retreat','Interpreter','latex','BackgroundColor','w','FontSize',14);
text(am,15,-12,'Decelerated retreat','Interpreter','latex','BackgroundColor','w','FontSize',14);

% Export
exportgraphics(figure(2),'figures/retreat_mK.pdf','ContentType','vector');

% Figure A1: retreat for m=K=1 and various Sigma

% Extracting data to plot
[H1, N1, xg1, qE1, t,~,qi1, U1] = get_data('a',1,1,SigmasAR(1),  Hd,XH,J,E,xi,delta);
[~, ~, xg2, ~, ~,~,qi2] = get_data('a',1,1,SigmasAR(2),  Hd,XH,J,E,xi,delta);
[~, ~, xg3, ~, ~,~,qi3] = get_data('a',1,1,SigmasAR(3),  Hd,XH,J,E,xi,delta);
[~, ~, xg4, ~, ~,~,qi4] = get_data('a',1,1,SigmasAR(4),  Hd,XH,J,E,xi,delta);
[H5, N5, xg5, qE5, ~,~,qi5, U5] = get_data('a',1,1,SigmasAR(5),  Hd,XH,J,E,xi,delta);

% Initialise axes for xg against t
figure(3); clf;
set(gcf(),'Position',[400 400 1000 1130],'Color','w');
axga  = make_nice_axes(subplot('Position',[0.15 0.75 0.8 0.225]), '$t$ / kyr', '$x_g$ / km', '');
xlim(axga,[0 31]); ylim(axga,[700 900]);

% Plot xg against t for different Sigma
plot(axga, td_kyr*t, xd_km*xg1, 'color', clrsr(6,:),  'LineWidth', 1.5);
plot(axga, td_kyr*t, xd_km*xg2, 'color', clrsr(5,:), 'LineWidth', 1.5);
plot(axga, td_kyr*t, xd_km*xg3, 'color', clrsr(4,:)  ,'LineWidth', 1.5);
plot(axga, td_kyr*t, xd_km*xg4, 'color', clrsr(3,:)  ,'LineWidth', 1.5);
plot(axga, td_kyr*t, xd_km*xg5, 'color', clrsr(2,:) , 'LineWidth', 1.5);

% Ice flux (additional figure)
plot(aqIa, xd_km*xg1, Hd*ud_myr*1e-3*qi1, 'color', clrsr(6,:), 'LineWidth', 1.5);
plot(aqIa, xd_km*xg2, Hd*ud_myr*1e-3*qi2, 'color', clrsr(5,:) ,'LineWidth', 1.5);
plot(aqIa, xd_km*xg3, Hd*ud_myr*1e-3*qi3, 'color', clrsr(4,:) ,'LineWidth', 1.5);
plot(aqIa, xd_km*xg4, Hd*ud_myr*1e-3*qi4, 'color', clrsr(3,:) ,'LineWidth', 1.5);
plot(aqIa, xd_km*xg5, Hd*ud_myr*1e-3*qi5, 'color', clrsr(2,:), 'LineWidth', 1.5);

% Initialise other axes for snapshots at given times
aha = make_nice_axes(subplot('Position',[0.15 0.55 0.8 0.17]), '', '$z$ / m', '',1000); 
ana = make_nice_axes(subplot('Position',[0.15 0.21 0.8 0.17]), '', '$N$ / MPa', '',1000);
aqa = make_nice_axes(subplot('Position',[0.15 0.42 0.8 0.09]), '', '$q_{Ec}$ / mm yr$^{-1}$', '',1000); 
afluxa = make_nice_axes(subplot('Position',[0.15 0.08 0.8 0.09]), '$x$ / km', '$H u $ / km$^2$ yr$^{-1}$', '',1000);
set(aha, 'YLim', [-1000 4000], 'XTickLabel', '');
set(ana, 'YLim', [0 8], 'XTickLabel', '');
set(afluxa, 'YLim', [0 .3]);
set(aqa, 'YLim', [-10 100], 'XTickLabel', '');

% Bed and sea level
plot(aha, xd_km*2*XH, Hd*bedf(2*XH), 'k-', 'LineWidth', 2)
plot(aha, xd_km*[0 2], [0 0], 'LineWidth',1, 'Color', [0.5 0.5 0.5]);

% Plot final state
load('retreat_advance/max_xg_m1_K1.mat','HUxN')
plot(axga, td_kyr*t([1 end]), xd_km*HUxN(2*J+4)*[1 1], 'color', [0.5 0.5 0.5],'linestyle','--', 'LineWidth', 2);
plot(aha, xd_km*HUxN(2*J+4)*XH, Hd*(bedf(HUxN(2*J+4)*XH)+ HUxN(1:J+1)), 'color', [0.5 0.5 0.5],'linestyle','--', 'LineWidth', 2);
plot(ana, xd_km*HUxN(2*J+4)*XH, Nd_MPa*HUxN(2*J+5:3*J+5), 'color', [0.5 0.5 0.5],'linestyle','--', 'LineWidth', 2);
plot(aha, xd_km*HUxN(2*J+4)*[1 1], Hd*bedf(HUxN(2*J+4))*[1 -delta/(1-delta)] ,'--','LineWidth',2,'Color',[0.5 0.5 0.5]);

% Legend for xg against t
legend(axga,[legend_text('S_s', SigmasAR*scales.Sigma_to_Ss, ' m$^{-1}$' ) ; 'New steady state'],...
     'FontSize',13,'Interpreter','latex','Position',[0.64 0.76 0.31 0.15])

% Initial states
plot(aqa, xd_km*xg1(1)*XH(1:J), md_mmyr*qE5(:,1),  'Color', clrsb(1,:)  , 'LineWidth', 2);  
plot(afluxa, xd_km*xg1(1)*XH, 1e-6*Hd*ud_myr*H1(:,1).*U1(2:end,1),  'Color', clrsb(1,:)  , 'LineWidth', 1); 
plot(afluxa, xd_km*xg5(1)*XH, 1e-6*Hd*ud_myr*H5(:,1).*U5(2:end,1),  'Color', clrsb(1,:)  , 'LineWidth', 2); 

ains2 = axes('Position',[0.86 0.3 0.09 0.08]); hold on;
set(ains2,'TickLabelInterpreter','latex');
xlabel(ains2,'$x-x_g$','FontSize',12,'Interpreter','latex')
xlim(ains2,[-100 0]);

% Go through and plot snapshots of H, N, qE at selected times
L = length(ksp);
for l = 1:L
    XG1 = xg1(ks(ksp(l)));
    XG5 = xg5(ks(ksp(l)));
    plot(aha, xd_km*XG1*XH, Hd*(bedf(XG1*XH)+ H1(:,ksp(l))), 'Color', clrsb(l+1,:)  , 'LineWidth', 1);
    plot(ana, xd_km*XG1*XH, Nd_MPa*N1(:,ksp(l)), 'Color', clrsb(l+1,:)  , 'LineWidth', 1);  
    plot(aha, xd_km*XG5*XH, Hd*(bedf(XG5*XH)+ H5(:,ksp(l))), 'Color', clrsb(l+1,:)  , 'LineWidth', 2);
    plot(ana, xd_km*XG5*XH, Nd_MPa*N5(:,ksp(l)), 'Color', clrsb(l+1,:)  , 'LineWidth', 2);
    plot(aqa, xd_km*XG5*XH(1:J), md_mmyr*qE5(:,ksp(l)-1),  'Color', clrsb(l+1,:)  , 'LineWidth', 2);  
    
    plot(afluxa, xd_km*XG1*XH, 1e-6*Hd*ud_myr*H1(:,ksp(l)-1).*U1(2:end,ksp(l)-1),  'Color', clrsb(l+1,:)  , 'LineWidth', 1); 
    plot(afluxa, xd_km*XG5*XH, 1e-6*Hd*ud_myr*H5(:,ksp(l)-1).*U5(2:end,ksp(l)-1),  'Color', clrsb(l+1,:)  , 'LineWidth', 2); 

    plot(aha, xd_km*XG1*[1 1], Hd*bedf(XG1)*[1 -delta/(1-delta)] ,':','LineWidth',1,'Color',clrsb(l+1,:));
    plot(aha, xd_km*XG5*[1 1], Hd*bedf(XG5)*[1 -delta/(1-delta)] ,':','LineWidth',2,'Color',clrsb(l+1,:));

    plot(ains2, xd_km*XG1*(XH-1), Nd_MPa*N1(:,ksp(l)), 'Color', clrsb(l+1,:)  , 'LineWidth', .5,'Linestyle','-');   
    plot(ains2, xd_km*XG5*(XH-1), Nd_MPa*N5(:,ksp(l)), 'Color', clrsb(l+1,:)  , 'LineWidth', 2);

end

% Plot initial state
load('retreat_advance/min_xg_m1_K1.mat','HUxN');
plot(aha, xd_km*HUxN(2*J+4)*XH, Hd*(bedf(HUxN(2*J+4)*XH)+ HUxN(1:J+1)), 'color', clrsb(1,:),'linestyle','-', 'LineWidth', 2);
plot(ana, xd_km*HUxN(2*J+4)*XH, Nd_MPa*HUxN(2*J+5:3*J+5), 'color', clrsb(1,:),'linestyle','-', 'LineWidth', 2);
plot(aha, xd_km*HUxN(2*J+4)*[1 1], Hd*bedf(HUxN(2*J+4))*[1 -delta/(1-delta)] ,':','LineWidth',2,'Color',clrsb(1,:));

% Legend for snapshots
legend(aha,[aha.Children(4*L+3), aha.Children(4*L+2), aha.Children(4*L+1)],...
    ['New steady state' ; legend_text('S_s', [0 max(SigmasAR)]*scales.Sigma_to_Ss, ' m$^{-1}$' ) ],...
    'FontSize',13,'Interpreter','latex','Position',[0.16 0.555 0.31 .09]);

% Subplot labels
text(axga,0.4,860,'(a)','Interpreter','latex','BackgroundColor','w','FontSize',14);
text(aha,950,3500,'(b)','Interpreter','latex','BackgroundColor','w','FontSize',14);
text(ana,950,2,'(d)','Interpreter','latex','BackgroundColor','w','FontSize',14);
text(aqa,950,90,'(c)','Interpreter','latex','BackgroundColor','w','FontSize',14);
text(afluxa,950,.25,'(e)','Interpreter','latex','BackgroundColor','w','FontSize',14);

% Align axis text
aha.YLabel.Position(1) = -96;
aqa.YLabel.Position(1) = -96;
ana.YLabel.Position(1) = -96;

% Export
exportgraphics(figure(3),'figures/advance.pdf','ContentType','vector')

text(aqI,880,275,'(a)','Interpreter','latex','BackgroundColor','w','FontSize',14);
text(aqIa,720,275,'(b)','Interpreter','latex','BackgroundColor','w','FontSize',14);
exportgraphics(figure(4),'figures/qg_ret_adv.pdf','ContentType','vector')

% Function to extract the required data from a saved file
function [H, N, xg, qE, t, ks, qI, U] = get_data(mode,mrate,K,Sig, Hd,XH,J,E,xi,delta)

    % Loading the relevant file
    if mode == 'r'
        load(sprintf("retreat_advance/retreat_m%g_K%g_Sig%g.mat", mrate, K, Sig),'Hs','Ns','ks','t','xgs','Us')
    elseif mode == 'a'
        load(sprintf("retreat_advance/advance_m%g_K%g_Sig%g.mat", mrate, K, Sig),'Hs','Ns','ks','t','xgs','Us')
    end

    % Resample to a consistent grid size
    if length(t) ~= 400
        tN = linspace(0,2,400).^2;
        xg = interp1(t,xgs,tN);
        H = interp1(t,Hs',tN)'; 
        N = interp1(t,Ns',tN)';  
        U = interp1(t,Us',tN)'; 
        Uend = interp1(t,Us(end,:),tN) ;
        ks = 1:400;
        t = tN;
    else 
        H = Hs; N = Ns;  xg = xgs; U = Us;
        Uend = Us(end,:);
    end
    

    % Calculating the compaction-driven exfiltration
    % qE = -Sig*(2000/Hd)*(...
    %     (H(1:J,2:end) -  H(1:J,1:end-1) - E*N(1:J,2:end)  + E*N(1:J,1:end-1) )./(t(ks(2:end)) - t(ks(1:end-1)))-...
    %      (((1./xg(ks(2:end))).*(xg(ks(2:end)) - xg(ks(1:end-1)) )./(t(ks(2:end)) - t(ks(1:end-1))))...
    %      .*XH(1:J).*( (H(2:J+1,2:end)-H(1:J,2:end))  -E*(N(2:J+1,2:end)-N(1:J,2:end))  )./(XH(2:J+1)-XH(1:J)) ));
    qE = -Sig*(2000/Hd)*(...
         ( (1-xi)*(1-delta)*(H(1:J,2:end) -  H(1:J,1:end-1)) - E*N(1:J,2:end)  + E*N(1:J,1:end-1) )./(t(ks(2:end)) - t(ks(1:end-1)))-...
          (((1./xg(ks(2:end))).*(xg(ks(2:end)) - xg(ks(1:end-1)) )./(t(ks(2:end)) - t(ks(1:end-1))))...
          .*XH(1:J).*( (1-xi)*(1-delta)*((H(2:J+1,2:end)-H(1:J,2:end)))  -E*(N(2:J+1,2:end)-N(1:J,2:end))  )./(XH(2:J+1)-XH(1:J)) ));
    
    qI = H(end,:).*Uend;
end