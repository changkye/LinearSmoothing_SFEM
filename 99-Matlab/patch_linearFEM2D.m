% linear finite element mehtod
% patch test

close all; clear all; clc;

format long;
restoredefaultpath;
path(path,'./sources');
path(path,'./models');

% Model parameters
model.bcType = 'lin_patch';
model.elemType = 'T6';
model.stressType = 'plane_stn';
model.params = [1 0];
model.numEls = 2;
model.flag = 1; model.txtFlag = 0;
param = str2func(['model_' model.bcType]);
model = param(model);

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
        text(xc,yc,num2str(in),'color','b');
    end
    for iel = 1:size(model.Elements,1)
        nod = model.Nodes(model.Elements(iel,:),:);
        text(mean(nod(:,1)),mean(nod(:,2)),num2str(iel),'color','r');
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

%     % fixed boundary nodes' coordinates
%     fixXY = model.Nodes(fixNode,:);
%     % 
%     bcdof = zeros(1,2*length(fixNode));
%     bcdof(1:2:2*length(fixNode)) = 2*fixNode - 1;
%     bcdof(2:2:2*length(fixNode)) = 2*fixNode;
%     % 
%     bcval = zeros(1,2*length(fixNode));
%     bcval(1:2:2*length(fixNode)) = 0.1*fixXY(:,1).^2 + 0.1*fixXY(:,1).*fixXY(:,2) + 0.2*fixXY(:,2).^2;
%     bcval(2:2:2*length(fixNode)) = 0.05*fixXY(:,1).^2 + 0.15*fixXY(:,1).*fixXY(:,2) + 0.1*fixXY(:,2).^2;	
% elseif strcmp(model.bcType,'SS_DBC')
% 	% get boundary nodes
%     botNodes = find(model.Nodes(:,2) == min(model.Nodes(:,2)));
%     topNodes = find(model.Nodes(:,2) == max(model.Nodes(:,2)));
%     rightNodes = find(model.Nodes(:,1) == max(model.Nodes(:,1)));
%     leftNodes = find(model.Nodes(:,1) == min(model.Nodes(:,1)));
%     fixNode = unique([botNodes; topNodes; rightNodes; leftNodes]);

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
% model.bcdof = bcdof; model.bcval = bcval;

% Linear Constitutive matrix
model = linearConstitutive(model);

if strcmp(model.bcType,'quad_patch')
   model = quad_patch_force_edge(model); 
end

% Linear FEM stiffness matrix
model = linearFEM2D_stiffness(model);

% % External force vector
% if strcmp(model.bcType,'quad_patch')
% 	F = zeros(2*size(model.Nodes,1),1);
% 	% Gauss points
%     if strcmp(model.elemType,'T3') 
%         ng = 1; element = 'TRIANGULAR';
%     elseif strcmp(model.elemType,'T6')
%         ng = 2; element = 'TRIANGULAR';
%     elseif strcmp(model.elemType,'Q4')
%         ng = 1; element = 'GAUSS';
%     elseif strcmp(model.elemType,'Q8')
%         ng = 2; element = 'GAUSS';
%     end
% 	[W,Q] = quadrature(ng,element,2);

% 	% exact body force
% 	bx = -0.2*model.Cmat(1,1) - 0.15*model.Cmat(1,2) - 0.55*model.Cmat(3,3);
% 	by = -0.1*model.Cmat(1,2) - 0.2*model.Cmat(2,2) - 0.2*model.Cmat(3,3);
	
% 	% loop over elements
% 	for iel = 1:size(model.Elements,1)
% 		% current element connectivity
% 		wkInd = model.Elements(iel,:);

% 		% current element coordinates
% 		wkX = model.Nodes(wkInd,:);

% 		% loop over Gauss points
% 		for ig = 1:size(W,1)
% 			% shape functions & their derivs.
% 			[N,dNdxi] = lagrange_basis(model.elemType,Q(ig,:));

% 			% Jacobian
% 			J0 = wkX'*dNdxi;
% 			detJ = det(J0);

% 			% assemble to global force vector 
% 			F(2*wkInd-1,1) = F(2*wkInd-1,1) + N*bx*detJ*W(ig);
% 			F(2*wkInd,1)   = F(2*wkInd,1)   + N*by*detJ*W(ig);
% 		end
% 	end
% 	model.F = F;
% else
% 	model.F = zeros(2*size(model.Nodes,1),1);
% end

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
		exactU(i,1) = 0.1*model.Nodes(i,1)^2 + 0.1*model.Nodes(i,1)*model.Nodes(i,2)...
			+ 0.2*model.Nodes(i,2)^2;
		exactU(i,2) = 0.05*model.Nodes(i,1)^2 + 0.15*model.Nodes(i,1)*model.Nodes(i,2)...
			+ 0.1*model.Nodes(i,2)^2;
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
    plotMesh(model.Nodes,model.Elements,model.elemType,'k--','LineWidth',0.8);
	plotMesh(newNodes,model.Elements,model.elemType,'r','LineWidth',1.5);
    for in = 1:size(newNodes,1)
        xc = newNodes(in,1) + 0.001;
        yc = newNodes(in,2) - 0.003;
        text(xc,yc,num2str(in),'color','b');
    end 
	for iel = 1:size(model.Elements,1)
        nod = newNodes(model.Elements(iel,:),:);
        text(mean(nod(:,1)),mean(nod(:,2)),num2str(iel),'color','r');
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
    saveas(gcf,[resultpath '/LinearFEM2D_' model.elemType '_' model.bcType '_' ...
        num2str(model.numEls(1)) 'x' num2str(model.numEls(2))],'epsc2');
	hold off
end
if model.flag == 1
    newNodes = model.Nodes + exactU;
    figure(3)
    hold on
    plotMesh(model.Nodes,model.Elements,model.elemType,'k--','LineWidth',0.8);
    plotMesh(newNodes,model.Elements,model.elemType,'r','LineWidth',1.5);
    for in = 1:size(newNodes,1)
        xc = newNodes(in,1) + 0.001;
        yc = newNodes(in,2) - 0.003;
        text(xc,yc,num2str(in),'color','b');
    end 
    for iel = 1:size(model.Elements,1)
        nod = newNodes(model.Elements(iel,:),:);
        text(mean(nod(:,1)),mean(nod(:,2)),num2str(iel),'color','r');
    end 
    title('Exact');
    saveas(gcf,[resultpath '/Exact_' model.elemType '_' model.bcType '_' ...
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

U_comp = [U(:,1) exactU(:,1) U(:,2) exactU(:,2)];

% save data
save([resultpath '/LinearFEM2D_' model.elemType '_' model.bcType '_' ...
    num2str(model.numEls(1)) 'x' num2str(model.numEls(2)) '.mat']);

% paraview
vtk = [resultpath '/LinearFEM2D_' model.elemType '_' model.bcType '_' ...
	num2str(model.numEls(1)) 'x' num2str(model.numEls(2)) '.vtu'];
VTKPostProcess(vtk,model,result);
exact.uu(1:2:2*size(model.Nodes,1),1) = exactU(:,1);
exact.uu(2:2:2*size(model.Nodes,1),1) = exactU(:,2);
VTKPostProcess([resultpath '/exact_' num2str(model.numEls(1)) 'x' num2str(model.numEls(2)) '.vtu'],...
    model,exact);