function model = FEM2D_T3_cantilever_force(model)
	% 
	F = zeros(2*size(model.Nodes,1),1);

	% get bounday nodes
	rightNodes = find(model.Nodes(:,1) == max(model.Nodes(:,1)));
	rightEdge = [];
	for i = 1:length(rightNodes)-1
		rightEdge = [rightEdge; rightNodes(i) rightNodes(i+1)];
	end

	D = model.L(2,2) - model.L(2,1);
	I = D*D*D/12.0;
	% Gauss points
	[W,Q] = quadrature(2,'GAUSS',1);

	% loop over edge
	for i = 1:size(rightEdge,1)
		% current edge connectivity
		wkInd = rightEdge(i,:);

		% current element coordinates
		wkX = model.Nodes(wkInd,:);
		edof = 2*wkInd;

		% loop over Gauss points
		for ig = 1:size(W,1)
			% shape functions & their derivs.
			[N,dNdxi] = lagrange_basis('L2',Q(ig,:));

			% Jacobian
			detJ = norm(dNdxi'*wkX(:,2));
			gpt = N'*wkX;

			fy = -model.P/(2*I)*(D*D/4 - gpt(2)*gpt(2));

			% assemble to global force vector 
			F(edof,1) = F(edof,1) + N*fy*detJ*W(ig);
		end
	end
	model.F = F;
end