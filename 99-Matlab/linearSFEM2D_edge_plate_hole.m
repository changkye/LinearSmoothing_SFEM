%
close all; clear all; clc;

format long;
restoredefaultpath;
path(path,'./sources');
path(path,'./gmsh');

% Model parameters
model.flag = 1;
model.txtFlag = 1;
model.elemType = 'T6';
model.stressType = 'plane_stn';
model.bcType = 'plate_hole';
model.param = [1e3, 0.3];		% Young's modulus, Poisson's ratio
model.P = 0.0;

% Create result folder
resultpath = ['./results/' model.bcType];
if (~exist(resultpath,'dir')); mkdir(resultpath); end

% Mesh generation
probType = '1';
model = msh2mat(model,probType);
% plot meshes
if model.flag == 1
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
tol = 1e-6;
btmNodes = find(abs(model.Nodes(:,2)-min(model.Nodes(:,2)))<tol);
topNodes = find(abs(model.Nodes(:,2)-max(model.Nodes(:,2)))<tol);
ltNodes = find(abs(model.Nodes(:,1)-min(model.Nodes(:,1)))<tol);
rtNodes = find(abs(model.Nodes(:,1)-max(model.Nodes(:,1)))<tol);

bcdof = []; bcval = [];
for i = 1:length(btmNodes)
    bcdof = [bcdof 2*btmNodes(i)-1 2*btmNodes(i)];

    gpt = model.Nodes(btmNodes(i),:);
    [S,U] = exactPlateHole(model,gpt,1);

    bcval = [bcval U(1) U(2)];
end
for i = 1:length(ltNodes)
    bcdof = [bcdof 2*ltNodes(i)-1 2*ltNodes(i)];

    gpt = model.Nodes(ltNodes(i),:);
    [S,U] = exactPlateHole(model,gpt,1);

    bcval = [bcval U(1) U(2)];
end
model.bcdof = bcdof; model.bcval = bcval;

% Linear Constitutive matrix
model = linearConstitutive(model);

% External force vector
model = forcePlateHole(model,rtNodes,topNodes);

% Linear SFEM stiffness matrix
if strcmp(model.elemType,'T3')
	model = linearSFEM2D_edge_stiffness_T3(model);
elseif strcmp(model.elemType,'T6')
	model = linearSFEM2D_edge_stiffness_T6_0(model);
end

% Solve the system
[Kmod,Fmod] = feaplyc2(model);
uu = Kmod\Fmod;
result.uu = uu;

U = [result.uu(1:2:end) result.uu(2:2:end)];

% Strain energy
stnE = 0.5*(uu'*model.K*uu);

% Exact solution
% if strcmp(model.bcType,'lin_patch')
% 	for i = 1:size(model.Nodes,1)
% 		exactU(i,1) = 0.1 + 0.1*model.Nodes(i,1) + 0.2*model.Nodes(i,2);
% 		exactU(i,2) = 0.05 + 0.15*model.Nodes(i,1) + 0.1*model.Nodes(i,2);
% 	end	
% elseif strcmp(model.bcType,'quad_patch')
% 	for i = 1:size(model.Nodes,1)
% 		exactU(i,1) = 0.1*model.Nodes(i,1)^2 + ...
%             0.1*model.Nodes(i,1)*model.Nodes(i,2) + 0.2*model.Nodes(i,2)^2;
% 		exactU(i,2) = 0.05*model.Nodes(i,1)^2 + ...
%             0.15*model.Nodes(i,1)*model.Nodes(i,2) + 0.1*model.Nodes(i,2)^2;
% 	end
% elseif strcmp(model.bcType,'SS_DBC')
%     for i = 1:size(model.Nodes,1)
%         exactU(i,1) = 1.0*model.Nodes(i,2);
%         exactU(i,2) = 0.0;
%     end
% elseif strcmp(model.bcType,'NOT_SS_DBC')
% 	for i = 1:size(model.Nodes,1)
% 		exactU(i,1) = 0.5*model.Nodes(i,2)^2;
% 		exactU(i,2) = 0.0;
% 	end
% end
% plot deformed shapes
if model.flag == 1
    scaleFactor = 1e2;
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
    bcTitle = 'Plate with a hole';
    title(bcTitle);
    saveas(gcf,[resultpath '/LinearSFEM2D_edge_' model.elemType '_' ...
        model.bcType '_' probType],'epsc2');
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
vtk = [resultpath '/LinearSFEM2D_edge_' model.elemType '_' model.bcType '_' ...
	num2str(model.numEls(1)) 'x' num2str(model.numEls(2)) '.vtu'];
VTKPostProcess(vtk,model,result);
exact.uu(1:2:2*size(model.Nodes,1),1) = exactU(:,1);
exact.uu(2:2:2*size(model.Nodes,1),1) = exactU(:,2);
VTKPostProcess([resultpath '/exact_' num2str(model.numEls(1)) 'x' ...
    num2str(model.numEls(2)) '.vtu'],model,exact);

% save data
save([resultpath '/LinearSFEM2D_edge_' model.elemType '_' model.bcType '_' ...
    num2str(model.numEls(1)) 'x' num2str(model.numEls(2)) '.mat']);