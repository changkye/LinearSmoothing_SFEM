function model = model_SS_NBC(elemType,params,numEls)
	% models for SS_DBC
	L = [0 1; 0 1];
	P = 0.6;
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
	fixNode = unique(btmNodes);
	%  
	bcdof = zeros(1,2*length(fixNode));
	bcdof(1:2:end) = 2*fixNode - 1;
	bcdof(2:2:end) = 2*fixNode;
	% 
	bcval = zeros(1,length(bcdof));
	
	% -------------------------------------------------------------------------
	%	external force
	% -------------------------------------------------------------------------
	F = zeros(2*size(Nodes,1),1);
	F = computeForces(F,Nodes,P);
	
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
%%
function F = computeForces(F,Nodes,P)
	NPE = 3; ng = 2; element = 'L3';
	% 
	leftNodes = find(Nodes(:,1) == min(Nodes(:,1)));
	rightNodes = find(Nodes(:,1) == max(Nodes(:,1)));
	topNodes = find(Nodes(:,2) == max(Nodes(:,2)));
	% 
	[~,t] = sort(Nodes(leftNodes,2));
	leftNodes = leftNodes(t);
	
	[~,t] = sort(Nodes(rightNodes,2));
	rightNodes = rightNodes(t);

	[~,t] = sort(Nodes(topNodes,1));
	topNodes = topNodes(t);
	%
	numEl = (length(leftNodes)-1)/(NPE-1);
	cnt = 0;
	for ie = 1:numEl
		for in = 1:NPE
			cnt = cnt + 1;
			leftEdge(ie,in) = leftNodes(cnt);
		end
		cnt = cnt - 1;
	end
	leftEdge(:,[2,3]) = leftEdge(:,[3,2]);

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

	numEl = (length(topNodes)-1)/(NPE-1);
	cnt = 0;
	for ie = 1:numEl
		for in = 1:NPE
			cnt = cnt + 1;
			topEdge(ie,in) = topNodes(cnt);
		end
		cnt = cnt - 1;
	end
	topEdge(:,[2,3]) = topEdge(:,[3,2]);
	% 
	[W,Q] = quadrature(ng,'GAUSS',1);
	for iel = 1:size(leftEdge,1)
		wkInd = leftEdge(iel,:);
		wkX = Nodes(wkInd,:);

		nn = length(wkInd);
		edof = 2*wkInd;
		
		for ig = 1:ng
			[N,dNdxi] = lagrange_basis(element,Q(ig,:));

			detJ = det(dNdxi'*wkX(:,2));
			
			F(edof,1) = F(edof,1) + N*(-P)*W(ig)*detJ;
		end
	end

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

	for iel = 1:size(topEdge,1)
		wkInd = topEdge(iel,:);
		wkX = Nodes(wkInd,:);

		nn = length(wkInd);
		edof = 2*wkInd - 1;
		
		for ig = 1:ng
			[N,dNdxi] = lagrange_basis(element,Q(ig,:));

			detJ = det(dNdxi'*wkX(:,1));
			
			F(edof,1) = F(edof,1) + N*P*W(ig)*detJ;
		end
	end

end
		


