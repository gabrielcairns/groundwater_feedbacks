% Plots the results of experiments involving grounding line advance and
% retreat (Figures 6, A1 and A2)

% Scales, parameters and functions
[scales, params, funcs] = default_values();
Hd = scales.Hd;     xd_km = scales.xd_km;     md_mmyr = scales.md_mmyr;
td_kyr = scales.td_kyr;     Nd_MPa = scales.Nd_MPa;
delta = params.delta;   E = params.E;

bedf = funcs.bedf;
bedxf = funcs.bedxf;
Hsbf = funcs.Hsbf;

% Meshgrid
load("variable_grid_300.mat","XA","XH","XU");
J = 300;

% Colormaps for plots
clrsr = cmocean('thermal',7);
clrsb = cmocean('haline',7);

% Folder to store output figures
if ~isfolder('figures');   mkdir('figures');  end

%% Running experiments
% Run retreat and advance experiments for m=1, K=1, Sigma=0,4,8,12
for Sigma = 0:4:12
    retreat_or_advance('r',1,1,Sigma);
    retreat_or_advance('a',1,1,Sigma);
end

% Retreat for Sigma = 0,2 with m=1 and various K
for K = 0:0.5:2.5
    retreat_or_advance('r',1,K,0);
    retreat_or_advance('r',1,K,2);
end

% Retreat for Sigma = 0,2 with K=1 and various constant m
for mrate = 0:0.5:2.5
    retreat_or_advance('r',mrate,1,0);
    retreat_or_advance('r',mrate,1,2);
end

%% Plotting results 

% Figure 6: retreat for m=K=1 and various Sigma

% Extracting the data
[H1, N1, xg1, ~, t, ks] = get_data('r', 1,1,0,  Hd,XH,J,E);
[~, ~, xg2, ~, ~, ~] = get_data('r', 1,1,4,  Hd,XH,J,E);
[~, ~, xg3, ~, ~, ~] = get_data('r', 1,1,8,  Hd,XH,J,E);
[H4, N4, xg4, qE4, ~, ~] = get_data('r', 1,1,12,  Hd,XH,J,E);

% Initialising axes for xg against t
figure(1); clf;
set(gcf(),'Position',[400 400 750 900],'Color','w');
axg  = subplot('Position',[0.15 0.75 0.8 0.225]);   
make_nice_axes(axg, '$t$ / kyr', '$x_g$ / km', '')
xlim(axg,[0 15.5]); ylim(axg,[725 910])

% Plotting xg against t for different cases in Sigma
plot(axg, td_kyr*t, xd_km*xg1, 'color', clrsr(6,:), 'LineWidth', 2);
plot(axg, td_kyr*t, xd_km*xg2, 'color', clrsr(5,:) ,'LineWidth', 2);
plot(axg, td_kyr*t, xd_km*xg3, 'color', clrsr(4,:) ,'LineWidth', 2);
plot(axg, td_kyr*t, xd_km*xg4, 'color', clrsr(3,:), 'LineWidth', 2);

% Initalising the rest of the axes (snapshots of H, N, qE)
ah = subplot('Position',[0.15 0.475 0.8 0.2]);  
an = subplot('Position',[0.15 0.075 0.8 0.2]); 
aq = subplot('Position',[0.15 0.325 0.8 0.1]); 
make_nice_axes(ah, '', '$H$ / m', '')
make_nice_axes(an, '$x$ / km', '$N$ / MPa', '')
make_nice_axes(aq, '', '$q_{Ec}$ / mm yr$^{-1}$', '')
set(ah, 'XLim', [0 1000], 'YLim', [-1000 4000], 'XTickLabel', '');
set(an, 'XLim', [0 1000], 'YLim', [0 9]);
set(aq, 'XLim', [0 1000], 'YLim', [0 220], 'XTickLabel', '');

% Bed, sea level
plot(ah, xd_km*2*XH, Hd*bedf(2*XH), 'k-', 'LineWidth', 2)
plot(ah, xd_km*[0 2], [0 0], 'LineWidth',1, 'Color', [0.5 0.5 0.5]);

% The entries at which we plot snapshots
ksp = [45 71 100 142 199];

% Plotting initial condition
load('retreat_advance/max_xg_m1_K1.mat','HUxN');
plot(ah, xd_km*HUxN(2*J+4)*XH, Hd*(bedf(HUxN(2*J+4)*XH)+ HUxN(1:J+1)), 'color', clrsb(1,:),'linestyle', '-', 'LineWidth', 2);
plot(an, xd_km*HUxN(2*J+4)*XH, Nd_MPa*HUxN(2*J+5:3*J+5), 'color',  clrsb(1,:),'linestyle','-', 'LineWidth', 2);
plot(ah, xd_km*HUxN(2*J+4)*[1 1], Hd*bedf(HUxN(2*J+4))*[1 -delta/(1-delta)] ,':','LineWidth',2,'Color',clrsb(1,:));

% Plotting snapshots of H, N, qE at each desired time
L = length(ksp); 
for l = 1:L
    XG1 = xg1(ks(ksp(l)));
    XG4 = xg4(ks(ksp(l)));
    plot(ah, xd_km*XG1*XH, Hd*(bedf(XG1*XH)+ H1(:,ksp(l))), 'Color', clrsb(l+1,:) , 'LineWidth', 1, 'Linestyle','-' );
    plot(an, xd_km*XG1*XH, Nd_MPa*N1(:,ksp(l)), 'Color', clrsb(l+1,:)  , 'LineWidth', 1,'Linestyle','-');   
    plot(ah, xd_km*XG4*XH, Hd*(bedf(XG4*XH)+ H4(:,ksp(l))), 'Color', clrsb(l+1,:)   , 'LineWidth', 2);
    plot(an, xd_km*XG4*XH, Nd_MPa*N4(:,ksp(l)), 'Color', clrsb(l+1,:)  , 'LineWidth', 2);
    plot(aq, xd_km*XG4*XH(1:J), md_mmyr*qE4(:,ksp(l)-1),  'Color', clrsb(l+1,:)  , 'LineWidth', 2);  

    plot(ah, xd_km*XG1*[1 1], Hd*bedf(XG1)*[1 -delta/(1-delta)] ,':','LineWidth',1,'Color',clrsb(l+1,:));
    plot(ah, xd_km*XG4*[1 1], Hd*bedf(XG4)*[1 -delta/(1-delta)] ,':','LineWidth',2,'Color',clrsb(l+1,:));
end

% Plotting new steady state
load('retreat_advance/min_xg_m1_K1.mat','HUxN');
plot(axg, td_kyr*t([1 end]), xd_km*HUxN(2*J+4)*[1 1], 'color', [0.5 0.5 0.5], 'linestyle', '--', 'LineWidth', 2);
plot(ah, xd_km*HUxN(2*J+4)*XH, Hd*(bedf(HUxN(2*J+4)*XH)+ HUxN(1:J+1)), 'color', [0.5 0.5 0.5],'linestyle','--', 'LineWidth', 2);
plot(an, xd_km*HUxN(2*J+4)*XH, Nd_MPa*HUxN(2*J+5:3*J+5), 'color', [0.5 0.5 0.5], 'linestyle', '--', 'LineWidth', 2);
plot(ah, xd_km*HUxN(2*J+4)*[1 1], Hd*bedf(HUxN(2*J+4))*[1 -delta/(1-delta)] ,'--','LineWidth',2,'Color',[0.5 0.5 0.5]);

% Legends
legend(axg,{'$S_s=0$ m$^{-1}$', '$S_s=7.9 \times 10^{-5}$ m$^{-1}$', '$S_s=1.58 \times 10^{-4}$ m$^{-1}$', '$S_s=2.37 \times 10^{-4}$ m$^{-1}$', 'New steady state'},...
    'FontSize',14,'Interpreter','latex','Position',[0.64 0.81 0.31 0.15])
legend(ah,[ah.Children(2), ah.Children(4*L+2), ah.Children(4*L+1)],...
    {'New steady state', '$S_s=0$ m$^{-1}$', '$S_s=2.37 \times 10^{-4}$ m$^{-1}$'},...
    'FontSize',14,'Interpreter','latex','Position',[0.16 0.48 0.31 0.1]);

% Subplot labels
text(axg,0.4,750,'(a)','Interpreter','latex','BackgroundColor','w','FontSize',14);
text(ah,950,3500,'(b)','Interpreter','latex','BackgroundColor','w','FontSize',14);
text(an,950,8,'(d)','Interpreter','latex','BackgroundColor','w','FontSize',14);
text(aq,950,190,'(c)','Interpreter','latex','BackgroundColor','w','FontSize',14);

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
am = subplot('Position',[0.1 0.425 0.375 0.55]);     hold(am,'on') ;
aK = subplot('Position',[0.55 0.425 0.375 0.55]);     hold(aK,'on') ;
make_nice_axes(am, '$t$ / kyr', '$\Delta x_g$ / km', '');
make_nice_axes(aK, '$t$ / kyr', '', '');
ylim(am,[-0.1 4.5]);       ylim(aK,[-0.1 4.5]);

% Range of m and K to work with
ms = 0:0.5:2.5;
Ks = 0:0.5:2.5;

% For each m, extract and plot the data
for l=1:length(ms)
    [~, ~, xg1, ~, ~, ~] = get_data('r', ms(l),1,0,  Hd,XH,J,E);
    [~, ~, xg2, ~, ~, ~] = get_data('r', ms(l),1,2,  Hd,XH,J,E);
    plot(am, td_kyr*t, xd_km*(xg1-xg2), 'LineWidth',2, 'Color', clrsr(l,:))
end

% For each K, extract and plot the data
for l=1:length(Ks)
    [~, ~, xg1, ~, ~, ~] = get_data('r', 1, Ks(l), 0,  Hd,XH,J,E);
    [~, ~, xg2, ~, ~, ~] = get_data('r', 1, Ks(l), 2,  Hd,XH,J,E);
    plot(aK, td_kyr*t, xd_km*(xg1-xg2), 'LineWidth',2, 'Color', clrsb(7-l,:))
end

% Legends
legend(am,{'$m=0$ mm yr$^{-1}$', '$m=5$ mm yr$^{-1}$','$m=10$ mm yr$^{-1}$',...
    '$m=15$ mm yr$^{-1}$','$m=20$ mm yr$^{-1}$','$m=25$ mm yr$^{-1}$'},...
    'FontSize',13,'Interpreter','latex','Position',[0.1 0.15 0.375 0.2],'numcolumns',2)
legend(aK,{'$k_{sb} = 0$ m$^2$', '$k_{sb} = 0.8 \times 10^{-12}$ m$^2$','$k_{sb} = 1.6 \times 10^{-12}$ m$^2$',...
    '$k_{sb} = 2.4 \times 10^{-12}$ m$^2$','$k_{sb} = 3.3 \times 10^{-12}$ m$^2$','$k_{sb} = 4.1 \times 10^{-12}$ m$^2$'},...
     'FontSize',13,'Interpreter','latex','Position',[0.5 0.15 0.425 0.2],'numcolumns',2)

% Subplot labels
text(am,14.5,4,'(a)','Interpreter','latex','BackgroundColor','w','FontSize',14);
text(aK,14.5,4,'(b)','Interpreter','latex','BackgroundColor','w','FontSize',14);

% Export
exportgraphics(figure(2),'figures/retreat_mK.pdf','ContentType','vector');

% Figure A1: retreat for m=K=1 and various Sigma

% Extracting data to plot
[H1, N1, xg1, qE1] = get_data('a',1,1,0,  Hd,XH,J,E);
[~, ~, xg2, ~] = get_data('a',1,1,4,  Hd,XH,J,E);
[~, ~, xg3, ~] = get_data('a',1,1,8,  Hd,XH,J,E);
[H4, N4, xg4, qE4] = get_data('a',1,1,12,  Hd,XH,J,E);

% Initialise axes for xg against t
figure(3); clf;
set(gcf(),'Position',[400 400 750 900],'Color','w');
axga  = subplot('Position',[0.15 0.75 0.8 0.225]); 
make_nice_axes(axga, '$t$ / kyr', '$x$ / km', '');
xlim(axga,[0 15.5]); ylim(axga,[735 920]);

% Plot xg against t for different Sigma
plot(axga, td_kyr*t, xd_km*xg1, 'color', clrsr(6,:),  'LineWidth', 2);
plot(axga, td_kyr*t, xd_km*xg2, 'color', clrsr(5,:), 'LineWidth', 2);
plot(axga, td_kyr*t, xd_km*xg3, 'color', clrsr(4,:)  ,'LineWidth', 2);
plot(axga, td_kyr*t, xd_km*xg4, 'color', clrsr(3,:) , 'LineWidth', 2);
 
% Initialise other axes for snapshots at given times
aha = subplot('Position',[0.15 0.475 0.8 0.2]);  
ana = subplot('Position',[0.15 0.075 0.8 0.2]); 
aqa = subplot('Position',[0.15 0.325 0.8 0.1]); 
make_nice_axes(aha, '', '$H$ / m', '')
make_nice_axes(ana, '$x$ / km', '$N$ / MPa', '')
make_nice_axes(aqa, '', '$q_{Ec}$ / mm yr$^{-1}$', '')
set(aha, 'XLim', [0 1000], 'YLim', [-1000 4000], 'XTickLabel', '');
set(ana, 'XLim', [0 1000], 'YLim', [0 12]);
set(aqa, 'XLim', [0 1000], 'YLim', [-40 0], 'XTickLabel', '');

% Bed and sea level
plot(aha, xd_km*2*XH, Hd*bedf(2*XH), 'k-', 'LineWidth', 2)
plot(aha, xd_km*[0 2], [0 0], 'LineWidth',1, 'Color', [0.5 0.5 0.5]);

% Plot initial condition
load('retreat_advance/max_xg_m1_K1.mat','HUxN')
plot(axga, td_kyr*t([1 end]), xd_km*HUxN(2*J+4)*[1 1], 'color', [0.5 0.5 0.5],'linestyle','--', 'LineWidth', 2);
plot(aha, xd_km*HUxN(2*J+4)*XH, Hd*(bedf(HUxN(2*J+4)*XH)+ HUxN(1:J+1)), 'color', [0.5 0.5 0.5],'linestyle','--', 'LineWidth', 2);
plot(ana, xd_km*HUxN(2*J+4)*XH, Nd_MPa*HUxN(2*J+5:3*J+5), 'color', [0.5 0.5 0.5],'linestyle','--', 'LineWidth', 2);
plot(aha, xd_km*HUxN(2*J+4)*[1 1], Hd*bedf(HUxN(2*J+4))*[1 -delta/(1-delta)] ,'--','LineWidth',2,'Color',[0.5 0.5 0.5]);

% Legend for xg against t
legend(axga,{'$S_s=0$ m$^{-1}$', '$S_s=7.9 \times 10^{-5}$ m$^{-1}$', '$S_s=1.58 \times 10^{-4}$ m$^{-1}$', '$S_s=2.37 \times 10^{-4}$ m$^{-1}$', 'New steady state'},...
    'FontSize',14,'Interpreter','latex','Position',[0.64 0.76 0.31 0.15])

% Go through and plot snapshots of H, N, qE at selected times
L = length(ksp);
for l = 1:L
    XG1 = xg1(ks(ksp(l)));
    XG4 = xg4(ks(ksp(l)));
    plot(aha, xd_km*XG1*XH, Hd*(bedf(XG1*XH)+ H1(:,ksp(l))), 'Color', clrsb(l+1,:)  , 'LineWidth', 1);
    plot(ana, xd_km*XG1*XH, Nd_MPa*N1(:,ksp(l)), 'Color', clrsb(l+1,:)  , 'LineWidth', 1);  
    plot(aha, xd_km*XG4*XH, Hd*(bedf(XG4*XH)+ H4(:,ksp(l))), 'Color', clrsb(l+1,:)  , 'LineWidth', 2);
    plot(ana, xd_km*XG4*XH, Nd_MPa*N4(:,ksp(l)), 'Color', clrsb(l+1,:)  , 'LineWidth', 2);
    plot(aqa, xd_km*XG4*XH(1:J), md_mmyr*qE4(:,ksp(l)-1),  'Color', clrsb(l+1,:)  , 'LineWidth', 2);  

    plot(aha, xd_km*XG1*[1 1], Hd*bedf(XG1)*[1 -delta/(1-delta)] ,':','LineWidth',1,'Color',clrsb(l+1,:));
    plot(aha, xd_km*XG4*[1 1], Hd*bedf(XG4)*[1 -delta/(1-delta)] ,':','LineWidth',2,'Color',clrsb(l+1,:));
end

% Plot new steady state
load('retreat_advance/min_xg_m1_K1.mat','HUxN');
plot(aha, xd_km*HUxN(2*J+4)*XH, Hd*(bedf(HUxN(2*J+4)*XH)+ HUxN(1:J+1)), 'color', clrsb(1,:),'linestyle','-', 'LineWidth', 2);
plot(ana, xd_km*HUxN(2*J+4)*XH, Nd_MPa*HUxN(2*J+5:3*J+5), 'color', clrsb(1,:),'linestyle','-', 'LineWidth', 2);
plot(aha, xd_km*HUxN(2*J+4)*[1 1], Hd*bedf(HUxN(2*J+4))*[1 -delta/(1-delta)] ,':','LineWidth',2,'Color',clrsb(1,:));

% Legend for snapshots
legend(aha,[aha.Children(2), aha.Children(4*L+2), aha.Children(4*L+1)],...
    {'New steady state', '$S_s=0$ m$^{-1}$', '$S_s=2.37 \times 10^{-4}$ m$^{-1}$'},...
    'FontSize',14,'Interpreter','latex','Position',[0.16 0.48 0.31 0.1])

% Subplot labels
text(axga,0.4,880,'(a)','Interpreter','latex','BackgroundColor','w','FontSize',14);
text(aha,950,3500,'(b)','Interpreter','latex','BackgroundColor','w','FontSize',14);
text(ana,950,10,'(d)','Interpreter','latex','BackgroundColor','w','FontSize',14);
text(aqa,950,-5,'(c)','Interpreter','latex','BackgroundColor','w','FontSize',14);

% Align axis text
aha.YLabel.Position(1) = -96;
aqa.YLabel.Position(1) = -96;
ana.YLabel.Position(1) = -96;

% Export
exportgraphics(figure(3),'figures/advance.pdf','ContentType','vector')


% Function to extract the required data from a saved file
function [H, N, xg, qE, t, ks] = get_data(mode,mrate,K,Sig, Hd,XH,J,E)

    % Loading the relevant file
    if mode == 'r'
        load(sprintf("retreat_advance/retreat_m%g_K%g_Sig%g.mat", mrate, K, Sig),'Hs','Ns','ks','t','xgs')
    elseif mode == 'a'
        load(sprintf("retreat_advance/advance_m%g_K%g_Sig%g.mat", mrate, K, Sig),'Hs','Ns','ks','t','xgs')
    end
    H = Hs; N = Ns;  xg = xgs; 

    % Calculating the compaction-driven exfiltration
    qE = -Sig*(2000/Hd)*(...
        (H(1:J,2:end) -  H(1:J,1:end-1) - E*N(1:J,2:end)  + E*N(1:J,1:end-1) )./(t(ks(2:end)) - t(ks(1:end-1)))-...
         (((1./xg(ks(2:end))).*(xg(ks(2:end)) - xg(ks(1:end-1)) )./(t(ks(2:end)) - t(ks(1:end-1))))...
         .*XH(1:J).*( (H(2:J+1,2:end)-H(1:J,2:end))  -E*(N(2:J+1,2:end)-N(1:J,2:end))  )./(XH(2:J+1)-XH(1:J)) ));
end

% Function to help initialise axes
function make_nice_axes(a, xlab, ylab, titl)
    hold(a, 'on');
    set(a, 'TickLabelInterpreter', 'latex', 'FontSize', 14);
    xlabel(a, xlab, 'FontSize', 14, 'Interpreter', 'latex');
    ylabel(a, ylab, 'FontSize', 14, 'Interpreter', 'latex');
    title(a, titl, 'FontSize', 14, 'Interpreter', 'latex');
end
