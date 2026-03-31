clc; clear all; close all;
format long;

% fem
fem1 = load('./NonlinearFEM2D_T6_BendingBlock_2x4.mat');
fem2 = load('./NonlinearFEM2D_T6_BendingBlock_2x8.mat');
fem3 = load('./NonlinearFEM2D_T6_BendingBlock_2x12.mat');
fem4 = load('./NonlinearFEM2D_T6_BendingBlock_2x16.mat');
fem5 = load('./NonlinearFEM2D_T6_BendingBlock_2x20.mat');
fem6 = load('./NonlinearFEM2D_T6_BendingBlock_2x24.mat');
fem7 = load('./NonlinearFEM2D_T6_BendingBlock_2x28.mat');
fem8 = load('./NonlinearFEM2D_T6_BendingBlock_2x32.mat');
fem1_relerrdisp = [fem1.result.Relerrdisp, fem2.result.Relerrdisp,...
    fem3.result.Relerrdisp, fem4.result.Relerrdisp, fem5.result.Relerrdisp,...
    fem6.result.Relerrdisp, fem7.result.Relerrdisp, fem8.result.Relerrdisp];
% 
fem1 = load('./NonlinearFEM2D_T6_BendingBlock_4x4.mat');
fem2 = load('./NonlinearFEM2D_T6_BendingBlock_4x8.mat');
fem3 = load('./NonlinearFEM2D_T6_BendingBlock_4x12.mat');
fem4 = load('./NonlinearFEM2D_T6_BendingBlock_4x16.mat');
fem5 = load('./NonlinearFEM2D_T6_BendingBlock_4x20.mat');
fem6 = load('./NonlinearFEM2D_T6_BendingBlock_4x24.mat');
fem7 = load('./NonlinearFEM2D_T6_BendingBlock_4x28.mat');
fem8 = load('./NonlinearFEM2D_T6_BendingBlock_4x32.mat');
fem2_relerrdisp = [fem1.result.Relerrdisp, fem2.result.Relerrdisp,...
    fem3.result.Relerrdisp, fem4.result.Relerrdisp, fem5.result.Relerrdisp,...
    fem6.result.Relerrdisp, fem7.result.Relerrdisp, fem8.result.Relerrdisp];


% csfem
csfem1 = load('./NonlinearSFEM2D_cell_T6_BendingBlock_2x4.mat');
csfem2 = load('./NonlinearSFEM2D_cell_T6_BendingBlock_2x8.mat');
csfem3 = load('./NonlinearSFEM2D_cell_T6_BendingBlock_2x12.mat');
csfem4 = load('./NonlinearSFEM2D_cell_T6_BendingBlock_2x16.mat');
csfem5 = load('./NonlinearSFEM2D_cell_T6_BendingBlock_2x20.mat');
csfem6 = load('./NonlinearSFEM2D_cell_T6_BendingBlock_2x24.mat');
csfem7 = load('./NonlinearSFEM2D_cell_T6_BendingBlock_2x28.mat');
csfem8 = load('./NonlinearSFEM2D_cell_T6_BendingBlock_2x32.mat');
csfem1_relerrdisp = [csfem1.result.Relerrdisp, csfem2.result.Relerrdisp,...
    csfem3.result.Relerrdisp, csfem4.result.Relerrdisp, csfem5.result.Relerrdisp,...
    csfem6.result.Relerrdisp, csfem7.result.Relerrdisp, csfem8.result.Relerrdisp];
% 
csfem1 = load('./NonlinearSFEM2D_cell_T6_BendingBlock_4x4.mat');
csfem2 = load('./NonlinearSFEM2D_cell_T6_BendingBlock_4x8.mat');
csfem3 = load('./NonlinearSFEM2D_cell_T6_BendingBlock_4x12.mat');
csfem4 = load('./NonlinearSFEM2D_cell_T6_BendingBlock_4x16.mat');
csfem5 = load('./NonlinearSFEM2D_cell_T6_BendingBlock_4x20.mat');
csfem6 = load('./NonlinearSFEM2D_cell_T6_BendingBlock_4x24.mat');
csfem7 = load('./NonlinearSFEM2D_cell_T6_BendingBlock_4x28.mat');
csfem8 = load('./NonlinearSFEM2D_cell_T6_BendingBlock_4x32.mat');
csfem2_relerrdisp = [csfem1.result.Relerrdisp, csfem2.result.Relerrdisp,...
    csfem3.result.Relerrdisp, csfem4.result.Relerrdisp, csfem5.result.Relerrdisp,...
    csfem6.result.Relerrdisp, csfem7.result.Relerrdisp, csfem8.result.Relerrdisp];
% 

% h = 2*[size(fem1.model.Nodes,1), size(fem2.model.Nodes,1), size(fem3.model.Nodes,1),...
%     size(fem4.model.Nodes,1), size(fem5.model.Nodes,1), size(fem6.model.Nodes,1),...
%     size(fem7.model.Nodes,1), size(fem8.model.Nodes,1)];
h = 4./[4 8 12 16 20 24 28 32];
%----- settings
set(0,'DefaultAxesFontName','Times New Roman');
set(0,'DefaultTextFontName','Times New Roman');
set(0,'DefaultAxesFontSize',17.5);
set(0,'DefaultAxesLineWidth',1.5);
set(0,'DefaultAxesFontWeight','bold');

set(0,'DefaultLineLineWidth',1.5);
set(0,'DefaultTextFontSize',15);
set(0,'DefaultTextFontWeight','bold');

% Displacements
figure
    hold on
    plot(abs(log(h)),log10(fem1_relerrdisp),'ko-','MarkerSize',10,'MarkerFaceColor','k');
    plot(abs(log(h)),log10(fem2_relerrdisp),'bs-','MarkerSize',10,'MarkerFaceColor','b');
    plot(abs(log(h)),log10(csfem1_relerrdisp),'r^-','MarkerSize',10,'MarkerFaceColor','r');
    plot(abs(log(h)),log10(csfem2_relerrdisp),'Color',CustomColors('Magenta'),'LineStyle',...
        '-','Marker','*','MarkerSize',10,'MarkerFaceColor',CustomColors('Magenta'));
    % xlab = get(gca,'xlabel');   
    % set(xlab,'string','|log(1/h)|','fontsize',15);
    % ylab = get(gca,'ylabel');
    % set(ylab,'string','log(U_{relative error})','fontsize',17);
    xlabel('$\left|\mathrm{log}\left(1/\mathrm{h}\right)\right|$','interpreter','latex','FontSize',14);
    ylabel('$\mathrm{log}\left(\mathbf{U}_\mathrm{relative\; error}\right)$','interpreter','latex','FontSize',14);
    fig = legend('FEM T6 (case 1)','FEM T6 (case 2)','CS-FEM T6 (case 1)','CS-FEM T6 (case 2)');
    set(fig,'Location','southeast');
    set(fig,'Box','on');
    set(fig,'FontSize',15);
    % set(gcf,'PaperPosition',[0 0 7.5 5.5]);
    % set(gcf,'PaperSize',[7.5 5.5]);
    % set(gcf,'PaperPosition',[0 0 10 7.5]);
    % set(gcf,'PaperSize',[10 7.5]);
    set(gcf,'PaperPosition',[0 0 20 15]);
    set(gcf,'PaperSize',[20 15]);
    % saveas(gcf,['Convergence_disp_BendingBlock'],'epsc2');
    fig = gcf; 
    print('Convergence_disp_BendingBlock','-dpdf');
hold off

% fem
fem1 = load('./NonlinearFEM2D_T6_BendingBlock_2x4.mat');
fem2 = load('./NonlinearFEM2D_T6_BendingBlock_2x8.mat');
fem3 = load('./NonlinearFEM2D_T6_BendingBlock_2x12.mat');
fem4 = load('./NonlinearFEM2D_T6_BendingBlock_2x16.mat');
fem5 = load('./NonlinearFEM2D_T6_BendingBlock_2x20.mat');
fem6 = load('./NonlinearFEM2D_T6_BendingBlock_2x24.mat');
fem7 = load('./NonlinearFEM2D_T6_BendingBlock_2x28.mat');
fem8 = load('./NonlinearFEM2D_T6_BendingBlock_2x32.mat');
fem1_relerrstne = [fem1.result.stnE(end), fem2.result.stnE(end),...
    fem3.result.stnE(end), fem4.result.stnE(end), fem5.result.stnE(end),...
    fem6.result.stnE(end), fem7.result.stnE(end), fem8.result.stnE(end)];
% 
fem1 = load('./NonlinearFEM2D_T6_BendingBlock_4x4.mat');
fem2 = load('./NonlinearFEM2D_T6_BendingBlock_4x8.mat');
fem3 = load('./NonlinearFEM2D_T6_BendingBlock_4x12.mat');
fem4 = load('./NonlinearFEM2D_T6_BendingBlock_4x16.mat');
fem5 = load('./NonlinearFEM2D_T6_BendingBlock_4x20.mat');
fem6 = load('./NonlinearFEM2D_T6_BendingBlock_4x24.mat');
fem7 = load('./NonlinearFEM2D_T6_BendingBlock_4x28.mat');
fem8 = load('./NonlinearFEM2D_T6_BendingBlock_4x32.mat');
fem2_relerrstne = [fem1.result.stnE(end), fem2.result.stnE(end),...
    fem3.result.stnE(end), fem4.result.stnE(end), fem5.result.stnE(end),...
    fem6.result.stnE(end), fem7.result.stnE(end), fem8.result.stnE(end)];


% csfem
csfem1 = load('./NonlinearSFEM2D_cell_T6_BendingBlock_2x4.mat');
csfem2 = load('./NonlinearSFEM2D_cell_T6_BendingBlock_2x8.mat');
csfem3 = load('./NonlinearSFEM2D_cell_T6_BendingBlock_2x12.mat');
csfem4 = load('./NonlinearSFEM2D_cell_T6_BendingBlock_2x16.mat');
csfem5 = load('./NonlinearSFEM2D_cell_T6_BendingBlock_2x20.mat');
csfem6 = load('./NonlinearSFEM2D_cell_T6_BendingBlock_2x24.mat');
csfem7 = load('./NonlinearSFEM2D_cell_T6_BendingBlock_2x28.mat');
csfem8 = load('./NonlinearSFEM2D_cell_T6_BendingBlock_2x32.mat');
csfem1_relerrstne = [csfem1.result.stnE(end), csfem2.result.stnE(end),...
    csfem3.result.stnE(end), csfem4.result.stnE(end), csfem5.result.stnE(end),...
    csfem6.result.stnE(end), csfem7.result.stnE(end), csfem8.result.stnE(end)];
% 
csfem1 = load('./NonlinearSFEM2D_cell_T6_BendingBlock_4x4.mat');
csfem2 = load('./NonlinearSFEM2D_cell_T6_BendingBlock_4x8.mat');
csfem3 = load('./NonlinearSFEM2D_cell_T6_BendingBlock_4x12.mat');
csfem4 = load('./NonlinearSFEM2D_cell_T6_BendingBlock_4x16.mat');
csfem5 = load('./NonlinearSFEM2D_cell_T6_BendingBlock_4x20.mat');
csfem6 = load('./NonlinearSFEM2D_cell_T6_BendingBlock_4x24.mat');
csfem7 = load('./NonlinearSFEM2D_cell_T6_BendingBlock_4x28.mat');
csfem8 = load('./NonlinearSFEM2D_cell_T6_BendingBlock_4x32.mat');
csfem2_relerrstne = [csfem1.result.stnE(end), csfem2.result.stnE(end),...
    csfem3.result.stnE(end), csfem4.result.stnE(end), csfem5.result.stnE(end),...
    csfem6.result.stnE(end), csfem7.result.stnE(end), csfem8.result.stnE(end)];


%----- settings
set(0,'DefaultAxesFontName','Times New Roman');
set(0,'DefaultTextFontName','Times New Roman');
set(0,'DefaultAxesFontSize',17.5);
set(0,'DefaultAxesLineWidth',1.5);
set(0,'DefaultAxesFontWeight','bold');

set(0,'DefaultLineLineWidth',1.5);
set(0,'DefaultTextFontSize',15);
set(0,'DefaultTextFontWeight','bold');

analytic = 4.485618*ones(1,length(h));
% Displacements
figure
    hold on
    plot(abs(log(h)),log(analytic./analytic),'k--');
    plot(abs(log(h)),log10(fem1_relerrstne./analytic),'ko-','MarkerSize',10,'MarkerFaceColor','k');
    plot(abs(log(h)),log10(fem2_relerrstne./analytic),'bs-','MarkerSize',10,'MarkerFaceColor','b');
    plot(abs(log(h)),log10(csfem1_relerrstne./analytic),'r^-','MarkerSize',10,...
        'MarkerFaceColor','r');
    plot(abs(log(h)),log10(csfem2_relerrstne./analytic),'Color',CustomColors('Magenta'),...
        'LineStyle','-','Marker','*','MarkerSize',10,'MarkerFaceColor',CustomColors('Magenta'));
    % xlab = get(gca,'xlabel');   
    % set(xlab,'string','|log(1/h)|','fontsize',15);
    % ylab = get(gca,'ylabel');
    % set(ylab,'string','log(W^{numerical}/W^{exact})','fontsize',17);
    xlabel('$\left|\mathrm{log}\left(1/\mathrm{h}\right)\right|$','interpreter','latex','FontSize',14);
    ylabel('$\mathrm{log}\left(\frac{\mathrm{W}^\mathrm{numerical}}{\mathrm{W}^\mathrm{exact}}\right)$','interpreter','latex','FontSize',14);
    fig = legend('Exact','FEM T6 (case 1)','FEM T6 (case 2)','CS-FEM T6 (case 1)','CS-FEM T6 (case 2)');
    set(fig,'Location','east');
    set(fig,'Box','on');
    set(fig,'FontSize',15);
    % set(gcf,'PaperPosition',[0 0 7.5 5.5]);
    % set(gcf,'PaperSize',[7.5 5.5]);
    % set(gcf,'PaperPosition',[0 0 10 7.5]);
    % set(gcf,'PaperSize',[10 7.5]);
    set(gcf,'PaperPosition',[0 0 20 15]);
    set(gcf,'PaperSize',[20 15]);
    % saveas(gcf,['Convergence_stne_BendingBlock'],'epsc2');
    fig = gcf;
    print('Convergence_stne_BendingBlock','-dpdf');
hold off