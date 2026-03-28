function model = model_BendingBlock(elemType,params,numEls)
	% models for SS_DBC
	L = [2 3; -2 2];
	P = 0.0;
	
	% -------------------------------------------------------------------------
	%	mesh generation
	% -------------------------------------------------------------------------
	mesh = str2func(['element' elemType]);
	[Nodes,Elements] = mesh(L,numEls);
	
	% -------------------------------------------------------------------------
	% 	boundary conditions
	% -------------------------------------------------------------------------
	% get boundary nodes
	alpha = 0.9;
	btmNodes = find(Nodes(:,2) == min(Nodes(:,2)));
	topNodes = find(Nodes(:,2) == max(Nodes(:,2)));
	rightNodes = find(Nodes(:,1) == max(Nodes(:,1)));
	leftNodes = find(Nodes(:,1) == min(Nodes(:,1)));

	% 
	% bcdof1 = zeros(1,2*length(btmNodes));
	% bcdof1(1:2:end) = 2*btmNodes - 1;
	% bcdof1(2:2:end) = 2*btmNodes;

	% bcdof2 = zeros(1,2*length(topNodes));
	% bcdof2(1:2:end) = 2*topNodes - 1;
	% bcdof2(2:2:end) = 2*topNodes;
	
	% bcdof3 = zeros(1,2*length(leftNodes));
	% bcdof3(1:2:end) = 2*leftNodes - 1;
	% bcdof3(2:2:end) = 2*leftNodes;
	
	% bcdof4 = zeros(1,2*length(rightNodes));
	% bcdof4(1:2:end) = 2*rightNodes - 1;
	% bcdof4(2:2:end) = 2*rightNodes;

	% bcdof = [bcdof1 bcdof2 bcdof3 bcdof4];
	% % 
	% bcval1 = zeros(1,2*length(btmNodes));
	% btmX = Nodes(btmNodes,:);
	% bcval1(1:2:end) = sqrt(2*alpha*btmX(:,1)).*cos(btmX(:,2)/alpha) - btmX(:,1);
	% bcval1(2:2:end) = sqrt(2*alpha*btmX(:,1)).*sin(btmX(:,2)/alpha) - btmX(:,2);

	% bcval2 = zeros(1,2*length(topNodes));
	% topX = Nodes(topNodes,:);
 %    bcval2(1:2:end) = sqrt(2*alpha*topX(:,1)).*cos(topX(:,2)/alpha) - topX(:,1);
	% bcval2(2:2:end) = sqrt(2*alpha*topX(:,1)).*sin(topX(:,2)/alpha) - topX(:,2);

	% bcval3 = zeros(1,2*length(leftNodes));
	% leftX = Nodes(leftNodes,:);
 %    bcval3(1:2:end) = sqrt(2*alpha*leftX(:,1)).*cos(leftX(:,2)/alpha) - leftX(:,1);
	% bcval3(2:2:end) = sqrt(2*alpha*leftX(:,1)).*sin(leftX(:,2)/alpha) - leftX(:,2);

	% bcval4 = zeros(1,2*length(rightNodes));
	% rightX = Nodes(rightNodes,:);
 %    bcval4(1:2:end) = sqrt(2*alpha*rightX(:,1)).*cos(rightX(:,2)/alpha) - rightX(:,1);
	% bcval4(2:2:end) = sqrt(2*alpha*rightX(:,1)).*sin(rightX(:,2)/alpha) - rightX(:,2);
	
 %    bcval = [bcval1 bcval2 bcval3 bcval4];

 	% 
	fixNode = unique([btmNodes',topNodes',rightNodes',leftNodes']);
	bcdof = zeros(1,2*length(fixNode));
	bcdof(1:2:end) = 2*fixNode - 1;
	bcdof(2:2:end) = 2*fixNode;

	bcval = zeros(1,2*length(fixNode));
	fixX = Nodes(fixNode,1);
    fixY = Nodes(fixNode,2);
    bcval(1:2:end) = sqrt(2*alpha*fixX).*cos(fixY/alpha) - fixX;
	bcval(2:2:end) = sqrt(2*alpha*fixX).*sin(fixY/alpha) - fixY;

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
