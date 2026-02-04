[scales, params, funcs] = default_values();
Hd = scales.Hd;     rhoi = scales.rhoi;     g = scales.g;
xdkm = scales.xdkm; md_mmyr = scales.md_mmyr; td = scales.td;
td_kyr = scales.td_kyr; Nd_MPa = scales.Nd_MPa;
mW = params.mW; n = params.n;  delta = params.delta;   C = params.C;
E = params.E;

load("variable_grid_300.mat","XA","XH","XU");
J = 300;

clrsr = cmocean('thermal',7);
clrsb = cmocean('haline',7);

if ~isfolder('figures');   mkdir('figures');  end

% arate = 1;

% bedf = @(x) -100/Hd - .1*x;  
% bedxf = @(x) - .1 + 0*x;  
% Hdf = @(x) 2000/Hd + 0*x;
bedf = funcs.bedf;
bedxf = funcs.bedxf;
Hsbf = funcs.Hsbf;

for Sigma = 0:4:12
    retreat_or_advance('r',1,1,Sigma);
    retreat_or_advance('a',1,1,Sigma);
end

for K = 0:0.5:2.5
    retreat_or_advance('r',1,K,0);
    retreat_or_advance('r',1,K,2);
end

for mrate = 0:0.5:2.5
    retreat_or_advance('r',mrate,1,0);
    retreat_or_advance('r',mrate,1,2);
end

% Cases to plot
%mrate = .25;  Sig1 = 0;   Sig2 = 8;   Sig3 = 16;
[H1, N1, xg1, ~, t, kgif] = get_data('r', 1,1,0,  Hd,XH,J,E);
[~, ~, xg2, ~, ~, ~] = get_data('r', 1,1,4,  Hd,XH,J,E);
[~, ~, xg3, ~, ~, ~] = get_data('r', 1,1,8,  Hd,XH,J,E);
[H4, N4, xg4, qE4, ~, ~] = get_data('r', 1,1,12,  Hd,XH,J,E);

% FIRST FIGURE - the standard tale of grounding line retreat for different
% sigma

% * The highest Sigma-case looks a little different to what we have
% elsewhere. The final steady state is in the wrong position
figure(1); clf;
set(gcf(),'Position',[400 400 750 900],'Color','w');
axg  = subplot('Position',[0.15 0.75 0.8 0.225]);   
make_nice_axes(axg, '$t$ / kyr', '$x_g$ / km', '')

plot(axg, td_kyr*t, xdkm*xg1, 'color', clrsr(6,:), 'LineWidth', 2);
plot(axg, td_kyr*t, xdkm*xg2, 'color', clrsr(5,:) ,'LineWidth', 2);
plot(axg, td_kyr*t, xdkm*xg3, 'color', clrsr(4,:) ,'LineWidth', 2);
plot(axg, td_kyr*t, xdkm*xg4, 'color', clrsr(3,:), 'LineWidth', 2);

%load(sprintf('retreat_advance/final_SS_m%g_K%g.mat',mrate,K),'HUxN')
%load('retreat_advance/min_xg_m1_K1.mat','HUxN')

xlim(axg,[0 15.5]); ylim(axg,[725 910])

ah = subplot('Position',[0.15 0.475 0.8 0.2]);  
an = subplot('Position',[0.15 0.075 0.8 0.2]); 
aq = subplot('Position',[0.15 0.325 0.8 0.1]); 

make_nice_axes(ah, '', '$H$ / m', '')
make_nice_axes(an, '$x$ / km', '$N$ / MPa', '')
make_nice_axes(aq, '', '$q_{Ec}$ / mm yr$^{-1}$', '')

set(ah, 'XLim', [0 1000], 'YLim', [-1000 4000], 'XTickLabel', '');
set(an, 'XLim', [0 1000], 'YLim', [0 9]);
set(aq, 'XLim', [0 1000], 'YLim', [0 210], 'XTickLabel', '');

plot(ah, xdkm*2*XH, Hd*bedf(2*XH), 'k-', 'LineWidth', 2)

%kgtp = [11 26 51 101 201];
kgtp = [11 32 57 87 125];
% Initial and final state
%load(sprintf('retreat_advance/IC_m%g_K%g',mrate,K),'HUxN')
load('retreat_advance/max_xg_m1_K1.mat','HUxN');
plot(ah, xdkm*HUxN(2*J+4)*XH, Hd*(bedf(HUxN(2*J+4)*XH)+ HUxN(1:J+1)), 'color', clrsb(1,:),'linestyle', '-', 'LineWidth', 2);
plot(an, xdkm*HUxN(2*J+4)*XH, Nd_MPa*HUxN(2*J+5:3*J+5), 'color',  clrsb(1,:),'linestyle','-', 'LineWidth', 2);

%load(sprintf('retreat_advance/final_SS_m%g_K%g.mat',mrate,K),'HUxN')
load('retreat_advance/min_xg_m1_K1.mat','HUxN');
plot(axg, td_kyr*t([1 end]), xdkm*HUxN(2*J+4)*[1 1], 'color', [0.5 0.5 0.5], 'linestyle', '--', 'LineWidth', 2);
plot(ah, xdkm*HUxN(2*J+4)*XH, Hd*(bedf(HUxN(2*J+4)*XH)+ HUxN(1:J+1)), 'color', [0.5 0.5 0.5],'linestyle','--', 'LineWidth', 2);
plot(an, xdkm*HUxN(2*J+4)*XH, Nd_MPa*HUxN(2*J+5:3*J+5), 'color', [0.5 0.5 0.5], 'linestyle', '--', 'LineWidth', 2);

L = length(kgtp);
for l = 1:L
    XG1 = xg1(kgif(kgtp(l)));
    XG4 = xg4(kgif(kgtp(l)));
    plot(ah, xdkm*XG1*XH, Hd*(bedf(XG1*XH)+ H1(:,kgtp(l))), 'Color', clrsb(l+1,:) , 'LineWidth', 1, 'Linestyle','-' );
    plot(an, xdkm*XG1*XH, Nd_MPa*N1(:,kgtp(l)), 'Color', clrsb(l+1,:)  , 'LineWidth', 1,'Linestyle','-');   
    plot(ah, xdkm*XG4*XH, Hd*(bedf(XG4*XH)+ H4(:,kgtp(l))), 'Color', clrsb(l+1,:)   , 'LineWidth', 2);
    plot(an, xdkm*XG4*XH, Nd_MPa*N4(:,kgtp(l)), 'Color', clrsb(l+1,:)  , 'LineWidth', 2);
    plot(aq, xdkm*XG4*XH(1:J), md_mmyr*qE4(:,kgtp(l)-1),  'Color', clrsb(l+1,:)  , 'LineWidth', 2);  
end

legend(axg,{'$S_s=0$ m$^{-1}$', '$S_s=7.9 \times 10^{-5}$ m$^{-1}$', '$S_s=1.58 \times 10^{-4}$ m$^{-1}$', '$S_s=2.37 \times 10^{-4}$ m$^{-1}$', 'Steady state'},...
    'FontSize',14,'Interpreter','latex','Position',[0.64 0.81 0.31 0.15])
legend(ah,[ah.Children(2*L+1), ah.Children(2*L), ah.Children(2*L-1)],...
    {'Steady state', '$S_s=0$ m$^{-1}$', '$S_s=2.37 \times 10^{-4}$ m$^{-1}$'},...
    'FontSize',14,'Interpreter','latex','Position',[0.16 0.48 0.31 0.1]);

text(axg,0.4,750,'(a)','Interpreter','latex','BackgroundColor','w','FontSize',14);
text(ah,950,3500,'(b)','Interpreter','latex','BackgroundColor','w','FontSize',14);
text(an,950,8,'(d)','Interpreter','latex','BackgroundColor','w','FontSize',14);
text(aq,950,190,'(c)','Interpreter','latex','BackgroundColor','w','FontSize',14);


ah.YLabel.Position(1) = -96;
aq.YLabel.Position(1) = -96;
an.YLabel.Position(1) = -96;
exportgraphics(figure(1),'figures/retreat.pdf','ContentType','vector')

% THIRD: effect of m,k
figure(3); clf; set(gcf(),'Position',[400 400 950 400],'Color','w');
am = subplot('Position',[0.1 0.425 0.375 0.55]);     hold(am,'on') ;
aK = subplot('Position',[0.55 0.425 0.375 0.55]);     hold(aK,'on') ;
make_nice_axes(am, '$t$ / kyr', '$\Delta x_g$ / km', '');
make_nice_axes(aK, '$t$ / kyr', '', '');

ms = 0:0.5:2.5;
Ks = 0:0.5:2.5;
for l=1:length(ms)
    [~, ~, xg1, ~, ~, ~] = get_data('r', ms(l),1,0,  Hd,XH,J,E);
    [~, ~, xg2, ~, ~, ~] = get_data('r', ms(l),1,2,  Hd,XH,J,E);
    plot(am, td_kyr*t, xdkm*(xg1-xg2), 'LineWidth',2, 'Color', clrsr(l,:))
end

for l=1:length(Ks)
    [~, ~, xg1, ~, ~, ~] = get_data('r', 1, Ks(l), 0,  Hd,XH,J,E);
    [~, ~, xg2, ~, ~, ~] = get_data('r', 1, Ks(l), 2,  Hd,XH,J,E);
    plot(aK, td_kyr*t, xdkm*(xg1-xg2), 'LineWidth',2, 'Color', clrsb(7-l,:))
end

ylim(am,[-0.1 4.5]);       ylim(aK,[-0.1 4.5]);
legend(am,{'$m=0$ mm yr$^{-1}$', '$m=5$ mm yr$^{-1}$','$m=10$ mm yr$^{-1}$',...
    '$m=15$ mm yr$^{-1}$','$m=20$ mm yr$^{-1}$','$m=25$ mm yr$^{-1}$'},...
    'FontSize',13,'Interpreter','latex','Position',[0.1 0.15 0.375 0.2],'numcolumns',2)
legend(aK,{'$k_{sb} = 0$ m$^2$', '$k_{sb} = 0.8 \times 10^{-12}$ m$^2$','$k_{sb} = 1.6 \times 10^{-12}$ m$^2$',...
    '$k_{sb} = 2.4 \times 10^{-12}$ m$^2$','$k_{sb} = 3.3 \times 10^{-12}$ m$^2$','$k_{sb} = 4.1 \times 10^{-12}$ m$^2$'},...
     'FontSize',13,'Interpreter','latex','Position',[0.5 0.15 0.425 0.2],'numcolumns',2)

text(am,14.5,4,'(a)','Interpreter','latex','BackgroundColor','w','FontSize',14);
text(aK,14.5,4,'(b)','Interpreter','latex','BackgroundColor','w','FontSize',14);

exportgraphics(figure(3),'figures/retreat_mK.pdf','ContentType','vector');


% Fourth: advance.
[H1, N1, xg1, qE1] = get_data('a',1,1,0,  Hd,XH,J,E);
[~, ~, xg2, ~] = get_data('a',1,1,4,  Hd,XH,J,E);
[~, ~, xg3, ~] = get_data('a',1,1,8,  Hd,XH,J,E);
[H4, N4, xg4, qE4] = get_data('a',1,1,12,  Hd,XH,J,E);

figure(4); clf;
set(gcf(),'Position',[400 400 750 900],'Color','w');
axga  = subplot('Position',[0.15 0.75 0.8 0.225]); 
make_nice_axes(axga, '$t$ / kyr', '$x$ / km', '')
plot(axga, td_kyr*t, xdkm*xg1, 'color', clrsr(6,:),  'LineWidth', 2);
plot(axga, td_kyr*t, xdkm*xg2, 'color', clrsr(5,:), 'LineWidth', 2);
plot(axga, td_kyr*t, xdkm*xg3, 'color', clrsr(4,:)  ,'LineWidth', 2);
plot(axga, td_kyr*t, xdkm*xg4, 'color', clrsr(3,:) , 'LineWidth', 2);

%load(sprintf('retreat_advance/IC_m%g_K%g.mat',mrate,K),'HUxN')
%load('retreat_advance/max_xg_m1_K1.mat','HUxN')

xlim(axga,[0 15.5]); ylim(axga,[735 915])
legend(axga,{'$S_s=0$ m$^{-1}$', '$S_s=7.9 \times 10^{-5}$ m$^{-1}$', '$S_s=1.58 \times 10^{-4}$ m$^{-1}$', '$S_s=2.37 \times 10^{-4}$ m$^{-1}$', 'Steady state'},...
    'FontSize',14,'Interpreter','latex','Position',[0.64 0.76 0.31 0.15])

% SECOND: SNAPSHOTS of H, N, qE during advance for different cases
aha = subplot('Position',[0.15 0.475 0.8 0.2]);  
ana = subplot('Position',[0.15 0.075 0.8 0.2]); 
aqa = subplot('Position',[0.15 0.325 0.8 0.1]); 

make_nice_axes(aha, '', '$H$ / m', '')
make_nice_axes(ana, '$x$ / km', '$N$ / MPa', '')
make_nice_axes(aqa, '', '$q_{Ec}$ / mm yr$^{-1}$', '')

set(aha, 'XLim', [0 1000], 'YLim', [-1000 4000], 'XTickLabel', '');
set(ana, 'XLim', [0 1000], 'YLim', [0 12]);
set(aqa, 'XLim', [0 1000], 'YLim', [-40 0], 'XTickLabel', '');

plot(aha, xdkm*2*XH, Hd*bedf(2*XH), 'k-', 'LineWidth', 2)

kgtp = [11 26 51 101 201];
%load(sprintf('retreat_advance/IC_m%g_K%g',mrate,K),'HUxN')
load('retreat_advance/max_xg_m1_K1.mat','HUxN')
plot(axga, td_kyr*t([1 end]), xdkm*HUxN(2*J+4)*[1 1], 'color', [0.5 0.5 0.5],'linestyle','--', 'LineWidth', 2);
plot(aha, xdkm*HUxN(2*J+4)*XH, Hd*(bedf(HUxN(2*J+4)*XH)+ HUxN(1:J+1)), 'color', [0.5 0.5 0.5],'linestyle','--', 'LineWidth', 2);
plot(ana, xdkm*HUxN(2*J+4)*XH, Nd_MPa*HUxN(2*J+5:3*J+5), 'color', [0.5 0.5 0.5],'linestyle','--', 'LineWidth', 2);

%load(sprintf('retreat_advance/final_SS_m%g_K%g.mat',mrate,K),'HUxN')
load('retreat_advance/max_xg_m1_K1.mat','HUxN');
plot(aha, xdkm*HUxN(2*J+4)*XH, Hd*(bedf(HUxN(2*J+4)*XH)+ HUxN(1:J+1)), 'color', clrsb(1,:),'linestyle','-', 'LineWidth', 2);
plot(ana, xdkm*HUxN(2*J+4)*XH, Nd_MPa*HUxN(2*J+5:3*J+5), 'color', clrsb(1,:),'linestyle','-', 'LineWidth', 2);

L = length(kgtp);
for l = 1:L
    XG1 = xg1(kgif(kgtp(l)));
    XG4 = xg4(kgif(kgtp(l)));
    plot(aha, xdkm*XG1*XH, Hd*(bedf(XG1*XH)+ H1(:,kgtp(l))), 'Color', clrsb(l+1,:)  , 'LineWidth', 1);
    plot(ana, xdkm*XG1*XH, Nd_MPa*N1(:,kgtp(l)), 'Color', clrsb(l+1,:)  , 'LineWidth', 1);  
    plot(aha, xdkm*XG4*XH, Hd*(bedf(XG4*XH)+ H4(:,kgtp(l))), 'Color', clrsb(l+1,:)  , 'LineWidth', 2);
    plot(ana, xdkm*XG4*XH, Nd_MPa*N4(:,kgtp(l)), 'Color', clrsb(l+1,:)  , 'LineWidth', 2);
    plot(aqa, xdkm*XG4*XH(1:J), 10*qE4(:,kgtp(l)-1),  'Color', clrsb(l+1,:)  , 'LineWidth', 2);  
end

legend(aha,[aha.Children(2*L+2), aha.Children(2*L), aha.Children(2*L-1)],...
    {'Steady state', '$S_s=0$ m$^{-1}$', '$S_s=2.37 \times 10^{-4}$ m$^{-1}$'},...
    'FontSize',14,'Interpreter','latex','Position',[0.16 0.48 0.31 0.1])

text(axga,0.4,880,'(a)','Interpreter','latex','BackgroundColor','w','FontSize',14);
text(aha,950,3500,'(b)','Interpreter','latex','BackgroundColor','w','FontSize',14);
text(ana,950,10,'(d)','Interpreter','latex','BackgroundColor','w','FontSize',14);
text(aqa,950,-5,'(c)','Interpreter','latex','BackgroundColor','w','FontSize',14);

aha.YLabel.Position(1) = -96;
aqa.YLabel.Position(1) = -96;
ana.YLabel.Position(1) = -96;
exportgraphics(figure(4),'figures/advance.pdf','ContentType','vector')


function [H, N, xg, qE, t, kgif] = get_data(mode,mrate,K,Sig, Hd,XH,J,E)
    if mode == 'r'
        load(sprintf("retreat_advance/retreat_m%g_K%g_Sig%g.mat", mrate, K, Sig),'Hgif','Ngif','kgif','t','xggif')
    elseif mode == 'a'
        load(sprintf("retreat_advance/advance_m%g_K%g_Sig%g.mat", mrate, K, Sig),'Hgif','Ngif','kgif','t','xggif')
    end
    H = Hgif; N = Ngif;  xg = xggif; 
    qE = -Sig*(2000/Hd)*(...
        (H(1:J,2:end) -  H(1:J,1:end-1) - E*N(1:J,2:end)  + E*N(1:J,1:end-1) )./(t(kgif(2:end)) - t(kgif(1:end-1)))-...
         (((1./xg(kgif(2:end))).*(xg(kgif(2:end)) - xg(kgif(1:end-1)) )./(t(kgif(2:end)) - t(kgif(1:end-1))))...
         .*XH(1:J).*( (H(2:J+1,2:end)-H(1:J,2:end))  -E*(N(2:J+1,2:end)-N(1:J,2:end))  )./(XH(2:J+1)-XH(1:J)) ));
end


function make_nice_axes(a, xlab, ylab, titl)
    hold(a, 'on');
    set(a, 'TickLabelInterpreter', 'latex', 'FontSize', 14);
    xlabel(a, xlab, 'FontSize', 14, 'Interpreter', 'latex');
    ylabel(a, ylab, 'FontSize', 14, 'Interpreter', 'latex');
    title(a, titl, 'FontSize', 14, 'Interpreter', 'latex');
end
