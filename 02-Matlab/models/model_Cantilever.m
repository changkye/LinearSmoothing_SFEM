function model = model_Cantilever(elemType,params,numEls)
	% models for SS_DBC
	Lx = [0 10; -1 1];
	P = 150.0; 	% for linear
    % P = 1/Lx(1,2);
	
	% -------------------------------------------------------------------------
	%	mesh generation
	% -------------------------------------------------------------------------
	mesh = str2func(['element' elemType]);
	[Nodes,Elements] = mesh(Lx,numEls);
	
	% -------------------------------------------------------------------------
	% 	boundary conditions
	% -------------------------------------------------------------------------
	% get boundary nodes
% 	E = (9*params(1)*params(2))/(3*params(2)+params(1));
% 	nu = (3*params(2)-2*params(1))/(2*(3*params(2)+params(1)));
    E = params(1); nu = params(2);
	L = Lx(1,2) - Lx(1,1);
	D = Lx(2,2) - Lx(2,1);
	I = D^3/12;
	PEI = P/(6*E*I);
	leftNodes = find(Nodes(:,1) == min(Nodes(:,1)));
	fixNode = unique(leftNodes);

	% fixed boundary nodes' coordinates
	fixX = Nodes(fixNode,:);
	% 
	bcdof = zeros(1,2*length(fixNode));
	bcdof(1:2:2*length(fixNode)) = 2*fixNode - 1;
	bcdof(2:2:2*length(fixNode)) = 2*fixNode;
	% 
	bcval = [];
	for i = 1:length(fixNode)
		ux = PEI*fixX(i,2)*((6*L-3*fixX(i,1))*fixX(i,1) + (2+nu)*fixX(i,2)*fixX(i,2) - D*D/4);
		uy = -PEI*(3*nu*fixX(i,2)*fixX(i,2)*(L-fixX(i,1)) + ...
			(4+5*nu)*D*D*fixX(i,1)/4 + (3*L-fixX(i,1))*fixX(i,1)*fixX(i,1));

		bcval = [bcval ux uy];
	end

	% -------------------------------------------------------------------------
	%	external force
	% -------------------------------------------------------------------------
	F = zeros(2*size(Nodes,1),1);
	% get bounday nodes
	rightNodes = find(Nodes(:,1) == max(Nodes(:,1)));
	rightEdge = [];
	if strcmp(elemType,'T3') | strcmp(elemType,'Q4')
		for i = 1:length(rightNodes)-1
			rightEdge = [rightEdge; rightNodes(i) rightNodes(i+1)];
		end
		edgeElem = 'L2';
	elseif strcmp(elemType,'T6')
        n = 1;
		for i = 1:numEls(2)
			rightEdge = [rightEdge; rightNodes(n)...
                rightNodes(n+2) rightNodes(n+1)];
            n = n + 2;
		end
		edgeElem = 'L3';
	end

	I = D*D*D/12.0;
	% Gauss points
	[W,Q] = quadrature(2,'GAUSS',1);

	% loop over edge
	for i = 1:size(rightEdge,1)
		% current edge connectivity
		wkInd = rightEdge(i,:);

		% current element coordinates
		wkX = Nodes(wkInd,:);
		edof = 2*wkInd;

		% loop over Gauss points
		for ig = 1:size(W,1)
			% shape functions & their derivs.
			[N,dNdxi] = lagrange_basis(edgeElem,Q(ig,:));

			% Jacobian
			detJ = norm(dNdxi'*wkX(:,2));
			gpt = N'*wkX;

			fy = -P/(2*I)*(D*D/4 - gpt(2)*gpt(2));

			% assemble to global force vector 
			F(edof,1) = F(edof,1) + N*fy*detJ*W(ig);
		end
	end


	% -------------------------------------------------------------------------
	model.elemType = elemType;
	model.L = Lx;
	model.param = params;
	model.numEls = numEls;
	model.Nodes = Nodes;
	model.Elements = Elements;
	model.bcdof = bcdof;
	model.bcval = bcval;
	model.F = F;
end
