% nonlinear cell-based smoothed finite element method
% patch test
% neo-Hookeam model
close all; clear all; clc;

format long;
rootDir = fileparts(mfilename('fullpath'));
addpath(fullfile(rootDir,'sources'));
addpath(fullfile(rootDir,'models'));

% Model parameters
bcType = 'NOT_SS_DBC';
% bcCase = '5';
% numEls = [12,4];
numEls = 8*ones(1,2);
param = str2func(['model_' bcType]);
% model = param('T6',[8 460],numEls);
model = param('T6',[0.6 100],numEls);
model.bcType = bcType;
model.flag = 0;
model.runMode = 'low_memory'; % 'low_memory' or 'fast'
model.useParfor = strcmp(model.runMode,'fast');

% Create result folder
resultpath = ['./results/' model.bcType];
if (~exist(resultpath,'dir')); mkdir(resultpath); end

% plot meshes
if model.flag == 1
	clf; axis equal; axis on; hold on;
    plotMesh(model.Nodes,model.Elements,model.elemType,'k');
	for in = 1:size(model.Nodes,1)
    	xc = model.Nodes(in,1) + 0.001;
    	yc = model.Nodes(in,2) - 0.003;
    	text(xc,yc,num2str(in),'color','blue');
	end
	for iel = 1:size(model.Elements,1)
		econ = model.Elements(iel,:);
		nod = model.Nodes(econ,:);
		text(mean(nod(:,1)),mean(nod(:,2)),num2str(iel),'color','r');
	end
end

% Load increment
uu = sparse(2*size(model.Nodes,1),1);
duu = sparse(2*size(model.Nodes,1),1);
nstep = 100;
Tolerance = 1e-9;
maxiter = 80;
stnE = zeros(nstep,1);

if isfield(model,'useParfor') && model.useParfor && license('test','Distrib_Computing_Toolbox')
    if isempty(gcp('nocreate'))
        parpool('threads');
    end
end

for istp = 1:nstep

    LoadFactor = istp/nstep;
    niter = 0.0;
    condition = 1.0;

    fprintf(1,'\n Step %d\t Scale %g%%\n',istp,LoadFactor*100);

    % Newton-Raphson
    while (condition>Tolerance) && (niter<maxiter)
        niter = niter + 1;

        % Global tangent stiffness & internal force
        model = nonlinearSFEM2D_cell_stiffness(model,uu);
        
        % Solve the system
        [Kmod,Bmod] = feaplyc2_nl(model,uu,LoadFactor);
        duu = Kmod\Bmod;
        uu = uu + duu;

        condition = norm(duu)/max(norm(uu),eps);
        resi_cond = norm(Bmod)/(2*size(model.Nodes,1));
             
        fprintf(1,'Iteration Number %d Condition %f Residual %f Tolerance %f\n',...
            niter,condition,resi_cond,Tolerance);
    end
    stnE(istp) = model.EE;
end


% Solve the system
result.uu = uu;
result.stnE = stnE;

result.U = [result.uu(1:2:end) result.uu(2:2:end)];

% Exact solution
if strcmp(model.bcType,'SS_DBC')  | strcmp(model.bcType,'SS_NBC')
	exactU = [model.Nodes(:,2), zeros(size(model.Nodes,1),1)];
elseif strcmp(model.bcType,'Comp_DBC') | strcmp(model.bcType,'Comp_NBC')
    exactU = [0.15*model.Nodes(:,1), -0.130434782608696*model.Nodes(:,2)];
elseif strcmp(model.bcType,'NOT_SS_DBC')
	exactU = [0.5*model.Nodes(:,2).^2, zeros(size(model.Nodes,1),1)];
elseif strcmp(model.bcType,'BendingBlock')
    alpha = 0.9;
    exactU = [sqrt(2*alpha*model.Nodes(:,1)).*cos(model.Nodes(:,2)/alpha)-model.Nodes(:,1),...
        sqrt(2*alpha*model.Nodes(:,1)).*sin(model.Nodes(:,2)/alpha)-model.Nodes(:,2)];
end
% plot deformed shape
if model.flag == 1
    textFlag = 0; 
	newNodes = model.Nodes + result.U;
	figure(2)
	hold on
    plotMesh(model.Nodes,model.Elements,model.elemType,'k--','LineWidth',0.8);
	plotMesh(newNodes,model.Elements,model.elemType,'r','LineWidth',1.5);
    if textFlag == 1
        for in = 1:size(newNodes,1)
            xc = newNodes(in,1) + 0.001;
            yc = newNodes(in,2) - 0.003;
            text(xc,yc,num2str(in),'color','blue');
        end 
        for iel = 1:size(model.Elements,1)
            econ = model.Elements(iel,:);
            nod = newNodes(econ,:);
            text(mean(nod(:,1)),mean(nod(:,2)),num2str(iel),'color','r');
        end 
    end
    if strcmp(model.bcType,'SS_DBC')
        bcTitle = 'Simple Shear with DBC';
    elseif strcmp(model.bcType,'SS_NBC')
        bcTitle = 'Simple shear with mixed DBC & NBC';
    elseif strcmp(model.bcType,'Comp_DBC')
        bcTitle = 'Compression with lateral extension with DBC';
    elseif strcmp(model.bcType,'Comp_NBC')
        bcTitle = 'Compression with lateral extension with mixed DBC & NBC';
    elseif strcmp(model.bcType,'NOT_SS_DBC')
        bcTitle = 'Not-So-Simple shear with DBC';
    elseif strcmp(model.bcType,'BendingBlock')
        bcTitle = 'Bending of a Block';
    elseif strcmp(model.bcType,'Cantilever')
        bcTitle = 'Cantilever';
    elseif strcmp(model.bcType,'PlateHole')
        bcTitle = 'Plate with a hole';
    elseif strcmp(model.bcType,'Indentation')
        bcTitle = 'Indentation';
    end
	title(bcTitle);
    saveas(gcf,[resultpath '/NonlinearSFEM2D_cell_' model.elemType '_' model.bcType '_' ...
        num2str(model.numEls(1)) 'x' num2str(model.numEls(2))],'epsc2');
	hold off
end

% Compute L2 norm relative error
if ~strcmp(model.bcType,'Cantilever') & ~strcmp(model.bcType,'PlateHole') &...
    ~strcmp(model.bcType,'Cook') & ~strcmp(model.bcType,'Indentation')
    diffU = exactU - result.U;
    Err = sum(diffU(:,1).^2 + diffU(:,2).^2);
    De = sum(exactU(:,1).^2 + exactU(:,2).^2);
    result.Relerrdisp = sqrt(Err/De)

    % [U(:,1) exactU(:,1) U(:,2) exactU(:,2)]

    % paraview
    exact.uu(1:2:2*size(model.Nodes,1),1) = exactU(:,1);
    exact.uu(2:2:2*size(model.Nodes,1),1) = exactU(:,2);
    VTKPostProcess([resultpath '/exact_' num2str(model.numEls(1)) 'x' num2str(model.numEls(2))...
        '.vtu'],model,exact);
end
if strcmp(model.bcType,'PlateHole')
    vtk = [resultpath '/NonlinearSFEM2D_cell_' model.elemType '_' model.bcType '_' bcCase '.vtu'];
    save([resultpath '/NonlinearSFEM2D_cell_' model.elemType '_' model.bcType '_' bcCase '.mat'],'model','result');
else
    vtk = [resultpath '/NonlinearSFEM2D_cell_' model.elemType '_' model.bcType '_' ...
        num2str(model.numEls(1)) 'x' num2str(model.numEls(2)) '.vtu'];
    save([resultpath '/NonlinearSFEM2D_cell_' model.elemType '_' model.bcType '_' ...
        num2str(model.numEls(1)) 'x' num2str(model.numEls(2)) '.mat'],'model','result');
end
VTKPostProcess(vtk,model,result);
