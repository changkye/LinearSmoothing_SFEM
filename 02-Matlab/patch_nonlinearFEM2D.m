% nonlinear finite element method
% patch test
% neo-Hookeam model
close all; clear all; clc;

format long;
restoredefaultpath;
path(path,'./sources');
path(path,'./models');

% Model parameters
bcType = 'NOT_SS_DBC';
% bcCase = '5';
numEls = 10*ones(1,2);
param = str2func(['model_' bcType]);
% model = param('T6',[8 460],numEls);
model = param('T6',[0.6 100],numEls);
model.bcType = bcType;
model.flag = 0;

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

for istp = 1:nstep

    LoadFactor = istp/nstep;
    niter = 0.0;
    condition = 1.0;

    fprintf(1,'\n Step %d\t Scale %g%%\n',istp,LoadFactor*100);
    
%     if istp == 24 keyboard; format short; else format long; end

    % Newton-Raphson
    while (condition>Tolerance) && (niter<maxiter)
        niter = niter + 1;

        % Global tangent stiffness & internal force
        model = nonlinearFEM2D_stiffness(model,uu);
        
        % Solve the system
        [Kmod,Bmod] = feaplyc2_nl(model,uu,LoadFactor);
        duu = Kmod\Bmod;
        uu = uu + duu;

        condition = sqrt(dot(full(duu),full(duu))/dot(uu,uu));
        resi_cond = sqrt(dot(full(Bmod),full(Bmod)))/(2*size(model.Nodes,1));
             
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
if strcmp(model.bcType,'SS_DBC') | strcmp(model.bcType,'SS_NBC')
	for i = 1:size(model.Nodes,1)
		exactU(i,1) = 1.0*model.Nodes(i,2);
		exactU(i,2) = 0.0;
	end
elseif strcmp(model.bcType,'Comp_DBC') | strcmp(model.bcType,'Comp_NBC')
    for i = 1:size(model.Nodes,1)
        exactU(i,1) = 0.15*model.Nodes(i,1);
        exactU(i,2) = -0.130434782608696*model.Nodes(i,2);
    end
elseif strcmp(model.bcType,'NOT_SS_DBC')
	for i = 1:size(model.Nodes,1)
		exactU(i,1) = 0.5*model.Nodes(i,2)^2;
		exactU(i,2) = 0.0;
	end
elseif strcmp(model.bcType,'BendingBlock')
    alpha = 0.9;
    for i = 1:size(model.Nodes,1)
        exactU(i,1) = sqrt(2*alpha*model.Nodes(i,1))*cos(model.Nodes(i,2)/alpha)...
            - model.Nodes(i,1);
        exactU(i,2) = sqrt(2*alpha*model.Nodes(i,1))*sin(model.Nodes(i,2)/alpha)...
            - model.Nodes(i,2);
    end
end
% plot deformed shapes
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
        bcTitle = 'Plate with hole';
    elseif strcmp(model.bcType,'Cook')
        bcTitle = 'Cook''s membrane';
    elseif strcmp(model.bcType,'Indentation')
        bcTitle = 'Indentation';
    end
	title(bcTitle);
    saveas(gcf,[resultpath '/NonlinearFEM2D_' model.elemType '_' model.bcType '_' ...
        num2str(model.numEls(1)) 'x' num2str(model.numEls(2))],'epsc2');
	hold off
end

% Compute L2 norm relative error
if ~strcmp(model.bcType,'Cantilever') & ~strcmp(model.bcType,'PlateHole') &...
    ~strcmp(model.bcType,'Cook') & ~strcmp(model.bcType,'Indentation')
    Err = 0; De = 0;
    for i = 1:size(model.Nodes,1)
	   Err = Err + ((exactU(i,1)-result.U(i,1))^2 + (exactU(i,2)-result.U(i,2))^2);
	   De = De + (exactU(i,1)^2 + exactU(i,2)^2);
    end
    result.Relerrdisp = sqrt(Err/De)

    % [U(:,1) exactU(:,1) U(:,2) exactU(:,2)]

    % paraview
    exact.uu(1:2:2*size(model.Nodes,1),1) = exactU(:,1);
    exact.uu(2:2:2*size(model.Nodes,1),1) = exactU(:,2);
    VTKPostProcess([resultpath '/exact_' num2str(model.numEls(1)) 'x' ...
        num2str(model.numEls(2)) '.vtu'],model,exact);
end
if strcmp(model.bcType,'PlateHole') | strcmp(model.bcType,'Cook')
    vtk = [resultpath '/NonlinearFEM2D_' model.elemType '_' model.bcType '_' bcCase '.vtu'];
    save([resultpath '/NonlinearFEM2D_' model.elemType '_' model.bcType '_' bcCase '.mat'],...
        'model','result');
else
    vtk = [resultpath '/NonlinearFEM2D_' model.elemType '_' model.bcType '_' ...
        num2str(model.numEls(1)) 'x' num2str(model.numEls(2)) '.vtu'];
    save([resultpath '/NonlinearFEM2D_' model.elemType '_' model.bcType '_' ...
        num2str(model.numEls(1)) 'x' num2str(model.numEls(2)) '.mat'],'model','result');
end
VTKPostProcess(vtk,model,result);