function model = model_lin_patch(model)
	% models for linear patch
	L = [0 1; 0 1];
	P = 0.0;
	numEls = model.numEls*ones(1,2);

	% -------------------------------------------------------------------------
	%	mesh generation
	% -------------------------------------------------------------------------
	mesh = str2func(['element' model.elemType]);
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
    bcdof(1:2:2*length(fixNode)) = 2*fixNode - 1;
    bcdof(2:2:2*length(fixNode)) = 2*fixNode;
    % 
    bcval = zeros(1,2*length(fixNode));
    bcval(1:2:2*length(fixNode)) = 0.1 + 0.1*fixXY(:,1) + 0.2*fixXY(:,2);
    bcval(2:2:2*length(fixNode)) = 0.05 + 0.15*fixXY(:,1) + 0.1*fixXY(:,2);

    % -------------------------------------------------------------------------
	%	external force
	% -------------------------------------------------------------------------
	F = zeros(2*size(Nodes,1),1);

	% -------------------------------------------------------------------------
	model.L = L;
	model.numEls = numEls;
	model.Nodes = Nodes;
	model.Elements = Elements;
	model.bcdof = bcdof;
	model.bcval = bcval;
	model.F = F;

end