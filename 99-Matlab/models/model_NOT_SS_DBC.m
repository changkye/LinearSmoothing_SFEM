function model = model_NOT_SS_DBC(elemType,params,numEls)
	% models for SS_DBC
	L = [0 2; 0 2];
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
	fixNode = unique([btmNodes; topNodes; rightNodes; leftNodes]);
	
	% fixed boundary nodes' coordinates
	fixXY = Nodes(fixNode,:);
	% 
	bcdof = zeros(1,2*length(fixNode));
	bcdof(1:2:end) = 2*fixNode - 1;
	bcdof(2:2:end) = 2*fixNode;
	% 
	bcval = zeros(1,2*length(fixNode));
	bcval(1:2:end) = 0.5*fixXY(:,2).^2;
	bcval(2:2:end) = 0.0;
	
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
