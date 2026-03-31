function model = model_Indentation(elemType,params,numEls)
	path(path,'./gmsh/Indentation/');

	P = -0.5;
    
    % -------------------------------------------------------------------------
	%	mesh generation
	% -------------------------------------------------------------------------
	mesh = load(['Indentation_' num2str(numEls(1)) 'x' num2str(numEls(2)) '.mat']);
	Nodes = mesh.Nodes;
	Elements = mesh.Elements;
	
	% -------------------------------------------------------------------------
	% 	boundary conditions
	% -------------------------------------------------------------------------
	% get boundary nodes
	leftNodes = find(Nodes(:,1) == min(Nodes(:,1)));
    btmNodes = find(Nodes(:,2) == min(Nodes(:,2)));
	
    % roller on left side
    roller_left = zeros(1,length(leftNodes));
    roller_left(1:end) = 2*leftNodes - 1 ;
    roller_btm = zeros(1,length(btmNodes));
    roller_btm(1:end) = 2*btmNodes;    
    % 
	bcdof = unique([roller_left,roller_btm]);
	bcval = zeros(1,length(bcdof));

	% -------------------------------------------------------------------------
	%	external force
    % -------------------------------------------------------------------------
    F = zeros(2*size(Nodes,1),1);
	F = computeForces(F,Nodes,P);

	% -------------------------------------------------------------------------
	model.elemType = elemType;
    model.numEls = numEls;
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
	topNodes = find(Nodes(:,2) == max(Nodes(:,2)) & Nodes(:,1) <= 1.0);
	% 
	[~,t] = sort(Nodes(topNodes,2));
	topNodes = topNodes(t);

	%
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
	for iel = 1:size(topEdge,1)
		wkInd = topEdge(iel,:);
		wkX = Nodes(wkInd,:);

		nn = length(wkInd);
		edof = 2*wkInd;
		
		for ig = 1:ng
			[N,dNdxi] = lagrange_basis(element,Q(ig,:));

			detJ = det(dNdxi'*wkX(:,1));
			
			F(edof,1) = F(edof,1) + N*P*W(ig)*detJ;
		end
	end
end
		
