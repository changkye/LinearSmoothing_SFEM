function model = model_Comp_DBC(elemType,params,numEls)
	% models for SS_DBC
	L = [0 1; 0 1];
	P = 0.0;
	numEls = numEls(1)*ones(1,2);

	% -------------------------------------------------------------------------
	%	mesh generation
	% -------------------------------------------------------------------------
	mesh = str2func(['element' elemType]);
	[Nodes,Elements] = mesh(L,numEls);
	
	% -------------------------------------------------------------------------
	% 	boundary conditions
	% -------------------------------------------------------------------------
	% get boundary nodes
	btmNodes = find(Nodes(:,2) == min(Nodes(:,2)));
	topNodes = find(Nodes(:,2) == max(Nodes(:,2)));
	rightNodes = find(Nodes(:,1) == max(Nodes(:,1)));
	leftNodes = find(Nodes(:,1) == min(Nodes(:,1)));
	fixNode = unique([btmNodes,topNodes,rightNodes,leftNodes]);
	
	% fixed boundary nodes' coordinates
	fixX = Nodes(fixNode,:);
	% 
	bcdof = zeros(1,2*length(fixNode));
	bcdof(1:2:end) = 2*fixNode - 1;
	bcdof(2:2:end) = 2*fixNode;
	% 
	bcval = zeros(1,length(bcdof));
	bcval(1:2:end) = 0.15*fixX(:,1);
	bcval(2:2:end) = -0.130434782608696*fixX(:,2);
	
	% -------------------------------------------------------------------------
	%	external force
	% -------------------------------------------------------------------------
	F = zeros(2*size(Nodes,1),1);

	% -------------------------------------------------------------------------
	model.elemType = elemType;
	model.L = L;
	model.param = params;
	model.numEls = numEls;
	model.Nodes = Nodes;
	model.Elements = Elements;
	model.bcdof = bcdof;
	model.bcval = bcval;
	model.F = F;
end
