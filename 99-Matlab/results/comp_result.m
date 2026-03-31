clc; clear all; close all;
format long;

% fem
fem1 = load('./Comp_DBC/NonlinearFEM2D_T6_Comp_DBC_2x2.mat');
fem2 = load('./Comp_DBC/NonlinearFEM2D_T6_Comp_DBC_4x4.mat');
fem3 = load('./Comp_DBC/NonlinearFEM2D_T6_Comp_DBC_6x6.mat');
fem4 = load('./Comp_DBC/NonlinearFEM2D_T6_Comp_DBC_8x8.mat');
fem5 = load('./Comp_DBC/NonlinearFEM2D_T6_Comp_DBC_10x10.mat');
fem_dbc_relerrdisp = [fem1.result.Relerrdisp, fem2.result.Relerrdisp,...
    fem3.result.Relerrdisp, fem4.result.Relerrdisp, fem5.result.Relerrdisp];

% csfem
csfem1 = load('./Comp_DBC/NonlinearSFEM2D_cell_T6_Comp_DBC_2x2.mat');
csfem2 = load('./Comp_DBC/NonlinearSFEM2D_cell_T6_Comp_DBC_4x4.mat');
csfem3 = load('./Comp_DBC/NonlinearSFEM2D_cell_T6_Comp_DBC_6x6.mat');
csfem4 = load('./Comp_DBC/NonlinearSFEM2D_cell_T6_Comp_DBC_8x8.mat');
csfem5 = load('./Comp_DBC/NonlinearSFEM2D_cell_T6_Comp_DBC_10x10.mat');
csfem_dbc_relerrdisp = [csfem1.result.Relerrdisp, csfem2.result.Relerrdisp,...
    csfem3.result.Relerrdisp, csfem4.result.Relerrdisp, csfem5.result.Relerrdisp];

% fem
fem1 = load('./Comp_NBC/NonlinearFEM2D_T6_Comp_NBC_2x2.mat');
fem2 = load('./Comp_NBC/NonlinearFEM2D_T6_Comp_NBC_4x4.mat');
fem3 = load('./Comp_NBC/NonlinearFEM2D_T6_Comp_NBC_6x6.mat');
fem4 = load('./Comp_NBC/NonlinearFEM2D_T6_Comp_NBC_8x8.mat');
fem5 = load('./Comp_NBC/NonlinearFEM2D_T6_Comp_NBC_10x10.mat');
fem_nbc_relerrdisp = [fem1.result.Relerrdisp, fem2.result.Relerrdisp,...
    fem3.result.Relerrdisp, fem4.result.Relerrdisp, fem5.result.Relerrdisp];

% csfem
csfem1 = load('./Comp_NBC/NonlinearSFEM2D_cell_T6_Comp_NBC_2x2.mat');
csfem2 = load('./Comp_NBC/NonlinearSFEM2D_cell_T6_Comp_NBC_4x4.mat');
csfem3 = load('./Comp_NBC/NonlinearSFEM2D_cell_T6_Comp_NBC_6x6.mat');
csfem4 = load('./Comp_NBC/NonlinearSFEM2D_cell_T6_Comp_NBC_8x8.mat');
csfem5 = load('./Comp_NBC/NonlinearSFEM2D_cell_T6_Comp_NBC_10x10.mat');
csfem_nbc_relerrdisp = [csfem1.result.Relerrdisp, csfem2.result.Relerrdisp,...
    csfem3.result.Relerrdisp, csfem4.result.Relerrdisp, csfem5.result.Relerrdisp];

h = 1./[2 4 6 8 10];

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
    plot(abs(log10(h)),log10(fem_dbc_relerrdisp),'ro-','MarkerSize',10);
    plot(abs(log10(h)),log10(fem_nbc_relerrdisp),'ro-.','MarkerSize',10,'MarkerFaceColor','r')
    plot(abs(log10(h)),log10(csfem_dbc_relerrdisp),'b^-','MarkerSize',10);
    plot(abs(log10(h)),log10(csfem_nbc_relerrdisp),'b^-.','MarkerSize',10,'MarkerFaceColor','b');
%     xlab = get(gca,'xlabel');   
%     set(xlab,'string','|log(1/h)|','fontsize',15);
    xlabel('$\left|\mathrm{log}\left(1/\mathrm{h}\right)\right|$','interpreter','latex','FontSize',14);
%     ylab = get(gca,'ylabel');
%     set(ylab,'string','log(U_{relative error})','fontsize',17);
    ylabel('$\mathrm{log}\left(\mathbf{U}_\mathrm{relative\; error}\right)$','interpreter','latex','FontSize',14);
    fig = legend('FEM T6 (Dirichlet)','FEM T6 (Dirichlet & Neumann)','CS-FEM T6 (Dirichlet)',...
    	'CS-FEM T6 (Dirichlet & Neumann)');
    set(fig,'Location','east');
    set(fig,'Box','on');
    set(fig,'FontSize',15);
    % set(gcf,'PaperPosition',[0 0 7.5 5.5]);
    % set(gcf,'PaperSize',[7.5 5.5]);
%     set(gcf,'PaperPosition',[0 0 10 7.5]);
%     set(gcf,'PaperSize',[10 7.5]);
    set(gcf,'PaperPosition',[0 0 20 15]);
    set(gcf,'PaperSize',[20 15]);
%     saveas(gcf,['Comp'],'epsc2');
    fig=gcf;                                     % your figure
%     fig.PaperPositionMode='auto';
    print('Comp','-dpdf');
hold off