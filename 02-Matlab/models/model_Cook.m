function model = model_Cook(elemType,params,probCase)
	path(path,'./gmsh')
    P = 1/16;

	% -------------------------------------------------------------------------
	%	mesh generation
	% -------------------------------------------------------------------------
	prob = load(['./gmsh/Cook/Cook_' probCase '.mat']);
	Nodes = prob.Nodes;
    Elements = prob.Elements;
	
	% -------------------------------------------------------------------------
	% 	boundary conditions
	% -------------------------------------------------------------------------
	% get boundary nodes
	leftNodes = find(Nodes(:,1) == min(Nodes(:,1)));
	fixNode = unique(leftNodes);

	% fixed boundary nodes' coordinates
	bcdof = zeros(1,2*length(fixNode));
	bcdof(1:2:2*length(fixNode)) = 2*fixNode - 1;
	bcdof(2:2:2*length(fixNode)) = 2*fixNode;
	% 
    bcval = zeros(1,2*length(fixNode));
	
	% -------------------------------------------------------------------------
	%	external force
    % -------------------------------------------------------------------------
    F = zeros(2*size(Nodes,1),1);
	F = computeForces(F,Nodes,P);

	% -------------------------------------------------------------------------
	model.elemType = elemType;
	model.param = params;
	model.Nodes = Nodes;
	model.Elements = Elements;
	model.bcdof = bcdof;
	model.bcval = bcval;
	model.F = F;
end

%%
function F = computeForces(F,Nodes,P)
	NPE = 3; ng = 2; element = 'L3';
	% 
	rightNodes = find(Nodes(:,1) == max(Nodes(:,1)));
	% 
	[~,t] = sort(Nodes(rightNodes,2));
	rightNodes = rightNodes(t);

	%
	numEl = (length(rightNodes)-1)/(NPE-1);
	cnt = 0;
	for ie = 1:numEl
		for in = 1:NPE
			cnt = cnt + 1;
			rightEdge(ie,in) = rightNodes(cnt);
		end
		cnt = cnt - 1;
	end
	rightEdge(:,[2,3]) = rightEdge(:,[3,2]);
    % 
	[W,Q] = quadrature(ng,'GAUSS',1);
	for iel = 1:size(rightEdge,1)
		wkInd = rightEdge(iel,:);
		wkX = Nodes(wkInd,:);

		nn = length(wkInd);
		edof = 2*wkInd;
		
		for ig = 1:ng
			[N,dNdxi] = lagrange_basis(element,Q(ig,:));

			detJ = det(dNdxi'*wkX(:,2));
			
			F(edof,1) = F(edof,1) + N*P*W(ig)*detJ;
		end
	end
end
		
