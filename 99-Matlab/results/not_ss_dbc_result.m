clc; clear all; close all;
format long;
set(groot,'DefaultFigureVisible','on');
set(groot,'DefaultFigureWindowStyle','normal');

baseDir = fileparts(mfilename('fullpath'));
if isempty(baseDir)
    % Fallback when running selected lines/sections from editor
    baseDir = fileparts(which('not_ss_dbc_result'));
    if isempty(baseDir)
        baseDir = pwd;
    end
end
dataDir = fullfile(baseDir,'NOT_SS_DBC');
levels = [2 4 6 8 10];
h = 1./levels;

% FEM
fem_relerrdisp = nan(1,length(levels));
fem_relerrstne = nan(1,length(levels));
for i = 1:length(levels)
    f = fullfile(dataDir,sprintf('NonlinearFEM2D_T6_NOT_SS_DBC_%dx%d.mat',levels(i),levels(i)));
    if exist(f,'file')
        s = load(f);
        fem_relerrdisp(i) = s.result.Relerrdisp;
        fem_relerrstne(i) = s.result.stnE(end);
    end
end

% CS-FEM
csfem_relerrdisp = nan(1,length(levels));
csfem_relerrstne = nan(1,length(levels));
for i = 1:length(levels)
    f = fullfile(dataDir,sprintf('NonlinearSFEM2D_cell_T6_NOT_SS_DBC_%dx%d.mat',levels(i),levels(i)));
    if exist(f,'file')
        s = load(f);
        csfem_relerrdisp(i) = s.result.Relerrdisp;
        csfem_relerrstne(i) = s.result.stnE(end);
    end
end

% ES-FEM (may be partially available)
esfem_relerrdisp = nan(1,length(levels));
esfem_relerrstne = nan(1,length(levels));
for i = 1:length(levels)
    f = fullfile(dataDir,sprintf('NonlinearSFEM2D_edge_T6_NOT_SS_DBC_%dx%d.mat',levels(i),levels(i)));
    if exist(f,'file')
        s = load(f);
        esfem_relerrdisp(i) = s.result.Relerrdisp;
        esfem_relerrstne(i) = s.result.stnE(end);
    end
end

%----- settings
set(0,'DefaultAxesFontName','Times New Roman');
set(0,'DefaultTextFontName','Times New Roman');
set(0,'DefaultAxesFontSize',17.5);
set(0,'DefaultAxesLineWidth',1.5);
set(0,'DefaultAxesFontWeight','bold');
set(0,'DefaultLineLineWidth',1.5);
set(0,'DefaultTextFontSize',15);
set(0,'DefaultTextFontWeight','bold');

% Displacement relative error
f1 = figure('Visible','on','WindowStyle','normal','Name','NOT_SS_DBC: Displacement Error');
hold on
x = abs(log10(h));
idx_fem = isfinite(fem_relerrdisp) & (fem_relerrdisp > 0);
idx_csfem = isfinite(csfem_relerrdisp) & (csfem_relerrdisp > 0);
idx_esfem = isfinite(esfem_relerrdisp) & (esfem_relerrdisp > 0);
if any(idx_fem)
    plot(x(idx_fem),log10(fem_relerrdisp(idx_fem)),'ro-','MarkerSize',10,'MarkerFaceColor','r');
end
if any(idx_csfem)
    plot(x(idx_csfem),log10(csfem_relerrdisp(idx_csfem)),'b^-','MarkerSize',10,'MarkerFaceColor','b');
end
if any(idx_esfem)
    plot(x(idx_esfem),log10(esfem_relerrdisp(idx_esfem)),'ks-','MarkerSize',10,'MarkerFaceColor','g');
end
xlabel('$\left|\mathrm{log}\left(1/\mathrm{h}\right)\right|$','interpreter','latex','FontSize',14);
ylabel('$\mathrm{log}\left(\mathbf{U}_\mathrm{relative\; error}\right)$','interpreter','latex','FontSize',14);
fig = legend('FEM T6','CS-FEM T6','ES-FEM T6');
set(fig,'Location','north');
set(fig,'Box','on');
set(fig,'FontSize',15);
set(gcf,'PaperPosition',[0 0 20 15]);
set(gcf,'PaperSize',[20 15]);
print(f1,fullfile(dataDir,'disp_NOT_SS_DBC'),'-dpdf');
drawnow expose;
hold off

% Strain energy ratio
analytic = 1.6*ones(1,length(h));
f2 = figure('Visible','on','WindowStyle','normal','Name','NOT_SS_DBC: Strain Energy Ratio');
hold on
plot(x,log10(analytic./analytic),'k-');
idx_fem_e = isfinite(fem_relerrstne) & (fem_relerrstne > 0);
idx_csfem_e = isfinite(csfem_relerrstne) & (csfem_relerrstne > 0);
idx_esfem_e = isfinite(esfem_relerrstne) & (esfem_relerrstne > 0);
if any(idx_fem_e)
    plot(x(idx_fem_e),log10(fem_relerrstne(idx_fem_e)./analytic(idx_fem_e)),'ro-','MarkerSize',10,'MarkerFaceColor','r');
end
if any(idx_csfem_e)
    plot(x(idx_csfem_e),log10(csfem_relerrstne(idx_csfem_e)./analytic(idx_csfem_e)),'b^-','MarkerSize',10,'MarkerFaceColor','b');
end
if any(idx_esfem_e)
    plot(x(idx_esfem_e),log10(esfem_relerrstne(idx_esfem_e)./analytic(idx_esfem_e)),'gs-','MarkerSize',10,'MarkerFaceColor','g');
end
xlabel('$\left|\mathrm{log}\left(1/\mathrm{h}\right)\right|$','interpreter','latex','FontSize',14);
ylabel('$\mathrm{log}\left(\frac{\mathrm{W}^\mathrm{numerical}}{\mathrm{W}^\mathrm{exact}}\right)$','interpreter','latex','FontSize',14);
fig = legend('Exact','FEM T6','CS-FEM T6','ES-FEM T6');
set(fig,'Location','east');
set(fig,'Box','on');
set(fig,'FontSize',15);
set(gcf,'PaperPosition',[0 0 20 15]);
set(gcf,'PaperSize',[20 15]);
print(f2,fullfile(dataDir,'stne_NOT_SS_DBC'),'-dpdf');
drawnow expose;
hold off

if ~any(idx_fem) && ~any(idx_csfem) && ~any(idx_esfem)
    warning('No valid displacement-error data found. Check MAT files in: %s', dataDir);
end

% Print strain-energy relative error values: |W_num/W_exact - 1|
fem_relerrstne_ratio = abs(fem_relerrstne./analytic - 1);
csfem_relerrstne_ratio = abs(csfem_relerrstne./analytic - 1);
esfem_relerrstne_ratio = abs(esfem_relerrstne./analytic - 1);

T = table(levels(:),fem_relerrstne_ratio(:),csfem_relerrstne_ratio(:),esfem_relerrstne_ratio(:), ...
    'VariableNames',{'Mesh','FEM_T6','CSFEM_T6','ESFEM_T6'});
disp('Strain-energy relative error |W_num/W_exact - 1|');
disp(T);
