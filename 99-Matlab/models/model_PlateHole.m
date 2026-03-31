function model = model_PlateHole(elemType,params,probCase)
	% material parameters for this problem
    % referred from Dubuis, Avril, Debayle, Badel, Identification of the material parameters of soft
    % tissues in the compressed leg, doi:10.1080/10255842.2011.50666
    %
    % shear modulus: 1) fat 16~32 kPa, 2) muscles 8~145 kPa
    % bulk modulus: 460~870 kPa
    %
    % for this problem mu=8; kappa=460
    path(path,'./gmsh')
    
	% -------------------------------------------------------------------------
	%	mesh generation
	% -------------------------------------------------------------------------
	prob = load(['./gmsh/PlateHole/PlateHole_' probCase '.mat']);
	Nodes = prob.Nodes;
    Elements = prob.Elements;
	
	% -------------------------------------------------------------------------
	% 	boundary conditions
	% -------------------------------------------------------------------------
	% get boundary nodes
	TOL = 1e-6;
	btmNodes = find(abs(Nodes(:,2) - min(Nodes(:,2)))<TOL); 
	topNodes = find(abs(Nodes(:,2) - max(Nodes(:,2)))<TOL); 
	leftNodes = find(abs(Nodes(:,1) - min(Nodes(:,1)))<TOL);
	rightNodes = find(abs(Nodes(:,1) - max(Nodes(:,1)))<TOL);

	bcdof = []; bcval = [];
	for i = 1:length(btmNodes)
		bcdof = [bcdof 2*btmNodes(i)-1 2*btmNodes(i)];

		gpt = Nodes(btmNodes(i),:);
		[~,U] = exactPlateHole_nl(params,gpt,1);

		bcval = [bcval U(1) U(2)];
	end
	% 
	for i = 1:length(leftNodes)
		bcdof = [bcdof 2*leftNodes(i)-1 2*leftNodes(i)];

		gpt = Nodes(leftNodes(i),:);
		[~,U] = exactPlateHole_nl(params,gpt,1);

		bcval = [bcval U(1) U(2)];
	end
	
	% -------------------------------------------------------------------------
	%	external force
	% -------------------------------------------------------------------------
	F = zeros(2*size(Nodes,1),1);
	[t1,t2] = sort(Nodes(rightNodes,2));
	rightNodes = rightNodes(t2);
	[t1,t2] = sort(Nodes(topNodes,1));
	topNodes = topNodes(t2);
	 
	% right edge
	nume = (length(rightNodes)-1)/2;
	cnt = 0;
	for in = 1:nume
		for jn = 1:3
			cnt = cnt + 1;
			rightEdge(in,jn) = rightNodes(cnt);
		end
		cnt = cnt - 1;
	end
	rightEdge(:,[2,3]) = rightEdge(:,[3,2]);

	% top edge
	cnt = 0;
    nume = (length(topNodes)-1)/2;
	for in = 1:nume
		for jn = 1:3
			cnt = cnt + 1;
			topEdge(in,jn) = topNodes(cnt);
		end
		cnt = cnt - 1;
	end
	topEdge(:,[2,3]) = topEdge(:,[3,2]);

	% 
	ng = 4;
	[W,Q] = quadrature(ng,'GAUSS',1);
	for el = 1:size(rightEdge,1)
		wkInd = rightEdge(el,:);
		wkX = Nodes(wkInd,:);

		nn = length(wkInd); edof = zeros(2,nn);
		edof(1,:) = 2*wkInd - 1; edof(2,:) = 2*wkInd;

		for ig = 1:ng
			[N,dNdxi] = lagrange_basis('L3',Q(ig,:));

			detJ = det(dNdxi'*wkX(:,2));
			xieta_gp = N'*wkX;

			[S,~] = exactPlateHole_nl(params,xieta_gp,1);

			F(edof(1,:),1) = F(edof(1,:),1) + N*S(1)*W(ig)*detJ;
			F(edof(2,:),1) = F(edof(2,:),1) + N*S(2)*W(ig)*detJ;
		end
	end

	for el = 1:size(topEdge,1)
		wkInd = topEdge(el,:);
		wkX = Nodes(wkInd,:);

		nn = length(wkInd); edof = zeros(2,nn);
		edof(1,:) = 2*wkInd - 1; edof(2,:) = 2*wkInd;

		for ig = 1:ng
			[N,dNdxi] = lagrange_basis('L3',Q(ig,:));

			detJ = det(dNdxi'*wkX(:,2));
			xieta_gp = N'*wkX;

			[S,~] = exactPlateHole_nl(params,xieta_gp,1);

			F(edof(1,:),1) = F(edof(1,:),1) + N*S(1)*W(ig)*detJ;
			F(edof(2,:),1) = F(edof(2,:),1) + N*S(2)*W(ig)*detJ;
		end
	end

	% -------------------------------------------------------------------------
	model.elemType = elemType;
	model.param = params;
	model.Nodes = Nodes;
	model.Elements = Elements;
	model.bcdof = bcdof;
	model.bcval = bcval;
	model.F = F;
end
