%
close all; clear all; clc;

format long;
restoredefaultpath;
path(path,'./sources');
path(path,'./models');

% Model parameters
% model.flag = 1;
% model.txtFlag = 0;
% model.elemType = 'T6';
% model.stressType = 'plane_stn';
% model.bcType = 'lin_patch';
% model.L = [0 1; 0 1];
% model.param = [1, 0];		% Young's modulus, Poisson's ratio
% model.P = 0.0;
% model.numEls = 1*ones(1,2);
bcType = 'quad_patch';
numEls = 10;
param = str2func(['model_' bcType]);
model = param('T6',[1, 0],numEls);
model.bcType = bcType;
model.stressType = 'plane_stn';
model.flag = 0; model.txtFlag = 0;


% Create result folder
resultpath = ['./results/' model.bcType];
if (~exist(resultpath,'dir')); mkdir(resultpath); end

% plot meshes
if model.flag == 1
    NPE = size(model.Elements,2);
	clf; axis equal; axis on; hold on;
	plot_mesh(model.Nodes,model.Elements,model.elemType,'k');
    if model.txtFlag == 1
        for in = 1:size(model.Nodes,1)
            xc = model.Nodes(in,1) + 0.001;
            yc = model.Nodes(in,2) - 0.003;
            text(xc,yc,num2str(in),'color','b');
        end
        for iel = 1:size(model.Elements,1)
            nod = model.Nodes(model.Elements(iel,:),:);
            text(mean(nod(:,1)),mean(nod(:,2)),num2str(iel),'color','r');
        end
    end
end

% % Dirichlet BCs
% if strcmp(model.bcType,'lin_patch')
% 	% get boundary nodes
%     botNodes = find(model.Nodes(:,2) == min(model.Nodes(:,2)));
%     topNodes = find(model.Nodes(:,2) == max(model.Nodes(:,2)));
%     rightNodes = find(model.Nodes(:,1) == max(model.Nodes(:,1)));
%     leftNodes = find(model.Nodes(:,1) == min(model.Nodes(:,1)));
%     fixNode = unique([botNodes; topNodes; rightNodes; leftNodes]);
% 
%     % fixed boundary nodes' coordinates
%     fixXY = model.Nodes(fixNode,:);
%     % 
%     bcdof = zeros(1,2*length(fixNode));
%     bcdof(1:2:2*length(fixNode)) = 2*fixNode - 1;
%     bcdof(2:2:2*length(fixNode)) = 2*fixNode;
%     % 
%     bcval = zeros(1,2*length(fixNode));
%     bcval(1:2:2*length(fixNode)) = 0.1 + 0.1*fixXY(:,1) + 0.2*fixXY(:,2);
%     bcval(2:2:2*length(fixNode)) = 0.05 + 0.15*fixXY(:,1) + 0.1*fixXY(:,2);
% elseif strcmp(model.bcType,'quad_patch')
% 	% get boundary nodes
%     botNodes = find(model.Nodes(:,2) == min(model.Nodes(:,2)));
%     topNodes = find(model.Nodes(:,2) == max(model.Nodes(:,2)));
%     rightNodes = find(model.Nodes(:,1) == max(model.Nodes(:,1)));
%     leftNodes = find(model.Nodes(:,1) == min(model.Nodes(:,1)));
%     fixNode = unique([botNodes; topNodes; rightNodes; leftNodes]);
% 
%     % fixed boundary nodes' coordinates
%     fixXY = model.Nodes(fixNode,:);
%     % 
%     bcdof = zeros(1,2*length(fixNode));
%     bcdof(1:2:2*length(fixNode)) = 2*fixNode - 1;
%     bcdof(2:2:2*length(fixNode)) = 2*fixNode;
%     % 
%     bcval = zeros(1,2*length(fixNode));
%     bcval(1:2:2*length(fixNode)) = 0.1*fixXY(:,1).^2 + ...
%         0.1*fixXY(:,1).*fixXY(:,2) + 0.2*fixXY(:,2).^2;
%     bcval(2:2:2*length(fixNode)) = 0.05*fixXY(:,1).^2 + ...
%         0.15*fixXY(:,1).*fixXY(:,2) + 0.1*fixXY(:,2).^2;	
% elseif strcmp(model.bcType,'SS_DBC')
%     % get boundary nodes
%     botNodes = find(model.Nodes(:,2) == min(model.Nodes(:,2)));
%     topNodes = find(model.Nodes(:,2) == max(model.Nodes(:,2)));
%     rightNodes = find(model.Nodes(:,1) == max(model.Nodes(:,1)));
%     leftNodes = find(model.Nodes(:,1) == min(model.Nodes(:,1)));
%     fixNode = unique([botNodes; topNodes; rightNodes; leftNodes]);
% 
%     % fixed boundary nodes' coordinates
%     fixXY = model.Nodes(fixNode,:);
%     % 
%     bcdof = zeros(1,2*length(fixNode));
%     bcdof(1:2:2*length(fixNode)) = 2*fixNode - 1;
%     bcdof(2:2:2*length(fixNode)) = 2*fixNode;
%     % 
%     bcval = zeros(1,2*length(fixNode));
%     bcval(1:2:2*length(fixNode)) = 1.0*fixXY(:,2);
%     bcval(2:2:2*length(fixNode)) = 0.0;
% elseif strcmp(model.bcType,'NOT_SS_DBC')
% 	% get boundary nodes
%     botNodes = find(model.Nodes(:,2) == min(model.Nodes(:,2)));
%     topNodes = find(model.Nodes(:,2) == max(model.Nodes(:,2)));
%     rightNodes = find(model.Nodes(:,1) == max(model.Nodes(:,1)));
%     leftNodes = find(model.Nodes(:,1) == min(model.Nodes(:,1)));
%     fixNode = unique([botNodes; topNodes; rightNodes; leftNodes]);
% 
%     % fixed boundary nodes' coordinates
%     fixXY = model.Nodes(fixNode,:);
%     % 
%     bcdof = zeros(1,2*length(fixNode));
%     bcdof(1:2:2*length(fixNode)) = 2*fixNode - 1;
%     bcdof(2:2:2*length(fixNode)) = 2*fixNode;
%     % 
%     bcval = zeros(1,2*length(fixNode));
%     bcval(1:2:2*length(fixNode)) = 0.5*fixXY(:,2).^2;
%     bcval(2:2:2*length(fixNode)) = 0.0;
% end
% 
% model.bcdof = bcdof; model.bcval = bcval;

% Linear Constitutive matrix
model = linearConstitutive(model);

% Linear SFEM stiffness matrix
if strcmp(model.elemType,'T3') | strcmp(model.elemType,'Q4')
	model = linearSFEM2D_cell_stiffness_T3_1(model);
elseif strcmp(model.elemType,'T6')
	model = linearSFEM2D_cell_stiffness_T6(model);
elseif strcmp(model.elemType,'Q8')
    model = linearSFEM2D_cell_stiffness_Q8(model);
elseif strcmp(model.elemType,'Q9')
    model = linearSFEM2D_cell_stiffness_Q9(model);
end

% External force vector
if strcmp(model.bcType,'quad_patch')
    model = quad_patch_force_cell(model);
else
    model.F = zeros(2*size(model.Nodes,1),1);
end

% Solve the system
[Kmod,Fmod] = feaplyc2(model);
uu = Kmod\Fmod;
result.uu = uu;

U = [result.uu(1:2:end) result.uu(2:2:end)];

% Strain energy
stnE = 0.5*(uu'*model.K*uu);

% Exact solution
if strcmp(model.bcType,'lin_patch')
	for i = 1:size(model.Nodes,1)
		exactU(i,1) = 0.1 + 0.1*model.Nodes(i,1) + 0.2*model.Nodes(i,2);
		exactU(i,2) = 0.05 + 0.15*model.Nodes(i,1) + 0.1*model.Nodes(i,2);
	end	
elseif strcmp(model.bcType,'quad_patch')
	for i = 1:size(model.Nodes,1)
		exactU(i,1) = 0.1*model.Nodes(i,1)^2 + ...
            0.1*model.Nodes(i,1)*model.Nodes(i,2) + 0.2*model.Nodes(i,2)^2;
		exactU(i,2) = 0.05*model.Nodes(i,1)^2 + ...
            0.15*model.Nodes(i,1)*model.Nodes(i,2) + 0.1*model.Nodes(i,2)^2;
	end
elseif strcmp(model.bcType,'SS_DBC')
    for i = 1:size(model.Nodes,1)
        exactU(i,1) = 1.0*model.Nodes(i,2);
        exactU(i,2) = 0.0;
    end
elseif strcmp(model.bcType,'NOT_SS_DBC')
	for i = 1:size(model.Nodes,1)
		exactU(i,1) = 0.5*model.Nodes(i,2)^2;
		exactU(i,2) = 0.0;
	end
end
% plot deformed shapes
if model.flag == 1
	newNodes = model.Nodes + U;
    figure(2)
    hold on
    plot_mesh(newNodes,model.Elements,model.elemType,'r');
    if model.txtFlag == 1
        for in = 1:size(newNodes,1)
            xc = newNodes(in,1) + 0.001;
            yc = newNodes(in,2) - 0.003;
            text(xc,yc,num2str(in),'color','b');
        end
        for iel = 1:size(model.Elements,1)
            nod = newNodes(model.Elements(iel,:),:);
            text(mean(nod(:,1)),mean(nod(:,2)),num2str(iel),'color','r');
        end
    end
    if strcmp(model.bcType,'lin_patch')
        bcTitle = 'Linear Patch';
    elseif strcmp(model.bcType,'quad_patch')
        bcTitle = 'Quadratic Patch';
    elseif strcmp(model.bcType,'SS_DBC')
        bcTitle = 'Simple Shear with DBC';
    elseif strcmp(model.bcType,'NOT_SS_DBC')
        bcTitle = 'Not-So-Simple shear with DBC';
    end
    title(bcTitle);
    saveas(gcf,[resultpath '/LinearSFEM2D_cell_' model.elemType '_' model.bcType '_' ...
        num2str(model.numEls(1)) 'x' num2str(model.numEls(2))],'epsc2');
    hold off
end

% Compute L2 norm relative error
Err = 0; De = 0;
for i = 1:size(model.Nodes,1)
	Err = Err + ((exactU(i,1)-U(i,1))^2 + (exactU(i,2)-U(i,2))^2);
	De = De + (exactU(i,1)^2 + exactU(i,2)^2);
end
Relerrdisp = sqrt(Err/De)

[U(:,1) exactU(:,1) U(:,2) exactU(:,2)];

% paraview
vtk = [resultpath '/LinearSFEM2D_cell_' model.elemType '_' model.bcType '_' ...
	num2str(model.numEls(1)) 'x' num2str(model.numEls(2)) '.vtu'];
VTKPostProcess(vtk,model,result);
exact.uu(1:2:2*size(model.Nodes,1),1) = exactU(:,1);
exact.uu(2:2:2*size(model.Nodes,1),1) = exactU(:,2);
VTKPostProcess([resultpath '/exact_' num2str(model.numEls(1)) 'x' ...
    num2str(model.numEls(2)) '.vtu'],model,exact);

% save data
save([resultpath '/LinearSFEM2D_cell_' model.elemType '_' model.bcType '_' ...
    num2str(model.numEls(1)) 'x' num2str(model.numEls(2)) '.mat']);