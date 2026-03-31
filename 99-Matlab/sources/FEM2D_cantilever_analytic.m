function result = FEM2D_cantilever_analytic(model,result)
	% 
	Nodes = model.Nodes; Elements = model.Elements;
	
	enNum = 0; enDenom = 0;
	edNum = 0; edDenom = 0;

	% Gauss points
	if strcmp(model.elemType,'T3')
		ng = 1;
	else
		ng = 2;
	end
	[W,Q] = quadrature(ng,'TRIANGULAR',2);

	% loop over elements
	for iel = 1:size(Elements,1)
		% current element connectivity
		wkInd = Elements(iel,:);

		% current element coordinates
		wkX = Nodes(wkInd,:);

		edof = zeros(1,2*length(wkInd));
		edof(1:2:end) = 2*wkInd - 1; 
		edof(2:2:end) = 2*wkInd;

		% current element displacements
		wkU = result.uu(edof);

		% loop over Gauss points
		for ig = 1:size(W,1)

			% shape functions
			[N,dNdxi] = lagrange_basis(model.elemType,Q(ig,:));

			% Jacobian
			J0 = dNdxi'*wkX;

			% map natural to physical
			dNdX = dNdxi*inv(J0);

			% strain-displacement matrix
			Bmat = zeros(3,2*length(wkInd));
			Bmat(1,1:2:end) = dNdX(:,1);
			Bmat(2,2:2:end) = dNdX(:,2);
			Bmat(3,1:2:end) = dNdX(:,2);
			Bmat(3,2:2:end) = dNdX(:,1);
			% 
			Uh = zeros(2,length(wkInd));
			Uh = [N'*wkU(1:2:end); N'*wkU(2:2:end)];

			% strain
			Eh = Bmat*wkU;

			% analytical solution
			xygp = N'*wkX;
			[U,E,S] = analytic_solution(model,xygp);

			% error in displacement
			edNum = edNum + ((U(1,1)-Uh(1,1))^2 + (U(2,1)-Uh(2,1))^2)*W(ig);
			edDenom = edDenom + (U(1,1)*U(1,1) + U(2,1)*U(2,1))*W(ig);

			% error in energy
% 			enNum = enNum + (Eh-E)'*model.Cmat*(Eh-E)*W(ig);
            enNum = enNum + (E-Eh)'*model.Cmat*(E-Eh)*W(ig);
			enDenom = enDenom + E'*model.Cmat*E*W(ig);

		end
	end

	result.L2norm = sqrt(edNum/edDenom);
    result.H1norm = sqrt(enNum/enDenom);
end

%%
function [U,E,S] = analytic_solution(model,xy)
	% 
	%-----------------------------------------------------------------------
	%  analytical solution:
	%  ux=Py/(6EI)*[(6L-3x)*x+(2+poisson)*(y^2-D^2/4)];
	%  uy=-P/(6EI)*[3*poisson*y^2*(L-x)+(4+5*poisson)*D^2*x/4+(3*L-x)*x^2];
	%-----------------------------------------------------------------------
	P = 150;
	D = model.L(2,2) - model.L(2,1);
	L = model.L(1,2) - model.L(1,1);
	E0 = model.param(1); nu = model.param(2);
	
	inert = D*D*D/12;
	pei = P/(6*E0*inert);

	uxy = []; E = []; S = [];
	% displacements
	uxy(1,:) = pei*xy(2)*((6*L - 3*xy(1))*xy(1) + (2+nu)*(xy(2)*xy(2) - D*D/4));
	uxy(2,:) = -pei*(3*nu*xy(2)*xy(2)*(L-xy(1)) + (4+5*nu)*D*D*xy(1)/4 + (3*L-xy(1))*xy(1)*xy(1));

	% strains
	E(1,:) = pei*6*(L-xy(1))*xy(2);
	E(2,:) = -pei*6*nu*xy(2)*(L-xy(1));
	E(3,:) = pei*6*(1+nu)*(xy(2)*xy(2)-D*D/4);

	% stresses
	S(1,:) = P*(L-xy(1))*xy(2)/inert;
	S(2,:) = 0;
	S(3,:) = P*(xy(2)*xy(2)-D*D/4)/2/inert;

	U = uxy;
end