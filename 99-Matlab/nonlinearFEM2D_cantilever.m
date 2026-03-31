%
close all; clear all; clc;

format long;
restoredefaultpath;
path(path,'./sources');

% Model parameters
model.flag = 1;
model.txtFlag = 0;
model.elemType = 'T6';
model.stressType = 'plane_stn';
model.bcType = 'cantilever';
model.L = [0 10; -1 1];
model.param = [3e7, 0.25];		% Young's modulus, Poisson's ratio
model.P = 150.0;
model.numEls = [10 8];

% Create result folder
resultpath = ['./results/' model.bcType];
if (~exist(resultpath,'dir')); mkdir(resultpath); end

% Mesh generation
model = meshGeneration2D(model);
% plot meshes
if model.flag == 0
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

% Dirichlet BCs
if strcmp(model.bcType,'cantilever')
	% get boundary nodes
	L = model.L(1,2) - model.L(1,1);
	D = model.L(2,2) - model.L(2,1);
	I = D^3/12;
	PEI = model.P/(6*model.param(1)*I);
	leftNodes = find(model.Nodes(:,1) == min(model.Nodes(:,1)));
	fixNode = unique(leftNodes);

	% fixed boundary nodes' coordinates
	fixX = model.Nodes(fixNode,:);
	% 
	bcdof = zeros(1,2*length(fixNode));
	bcdof(1:2:2*length(fixNode)) = 2*fixNode - 1;
	bcdof(2:2:2*length(fixNode)) = 2*fixNode;
	% 
	bcval = [];
	for i = 1:length(fixNode)
		ux = PEI*fixX(i,2)*((6*L-3*fixX(i,1))*fixX(i,1) + ...
			(2+model.param(2))*fixX(i,2)*fixX(i,2) - D*D/4);
		uy = -PEI*(3*model.param(2)*fixX(i,2)*fixX(i,2)*(L-fixX(i,1)) + ...
			(4+5*model.param(2))*D*D*fixX(i,1)/4 + (3*L-fixX(i,1))*fixX(i,1)*fixX(i,1));

		bcval = [bcval ux uy];
	end
end
model.bcdof = bcdof; model.bcval = bcval;

% Linear Constitutive matrix
model = linearConstitutive(model);

% Linear FEM stiffness matrix
model = linearFEM2D_stiffness(model);

% External force vector
model = FEM2D_cantilever_force(model);

% Solve the system
[Kmod,Fmod] = feaplyc2(model);
uu = Kmod\Fmod;
result.uu = uu;

U = [result.uu(1:2:end) result.uu(2:2:end)];

% Strain energy
stnE = 0.5*(uu'*model.K*uu);

% Exact solution
result = FEM2D_cantilever_analytic(model,result);
% plot deformed shapes
if model.flag == 1
	scaleFactor = 2e2;
	newNodes = model.Nodes + scaleFactor*U;
	figure(2)
	hold on
	plot_mesh(model.Nodes,model.Elements,model.elemType,'k--','LineWidth',0.8);
	plot_mesh(newNodes,model.Elements,model.elemType,'r-o','LineWidth',1.5);
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
	title('Cantilever');
	saveas(gcf,[resultpath '/LinearFEM2D_' model.elemType '_' model.bcType '_' ...
		num2str(model.numEls(1)) 'x' num2str(model.numEls(2))].'epsc2');
	hold off
end
result.L2norm
result.H1norm

% paraview
vtk = [resultpath '/LinearFEM2D_' model.elemType '_' model.bcType '_' ...
	num2str(model.numEls(1)) 'x' num2str(model.numEls(2)) '.vtu'];
VTKPostProcess(vtk,model,result);

% save data
save([resultpath '/LinearFEM2D_' model.elemType '_' model.bcType '_' ...
    num2str(model.numEls(1)) 'x' num2str(model.numEls(2)) '.mat']);