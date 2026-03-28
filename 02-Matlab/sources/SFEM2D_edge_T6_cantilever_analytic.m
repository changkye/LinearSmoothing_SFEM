function result = SFEM2D_edge_T6_cantilever_analytic(model,result)
	% 
	enNum = 0; enDenom = 0;
	edNum = 0; edDenom = 0;

	targetEdge = model.targetEdge;
    
    % loop over edges
    for ivo = 1:size(targetEdge,1)
        if targetEdge(ivo,end) == 0
            neighbour = targetEdge(ivo,3);
        else
            neighbour = targetEdge(ivo,3:end-2);
        end
        nc = length(neighbour);
        
        % 1 subcell
        if nc == 1 
        	% Gauss points
            [W,Q] = quadrature(2,'TRIANGULAR',2);
            
            % current subcell connectivity
            wkInd = model.Elements(neighbour,:);
            nndof = length(wkInd);
            edof(1:2:2*nndof) = 2*wkInd - 1;
            edof(2:2:2*nndof) = 2*wkInd;
            
            % current subcell coordinates
            wkX = model.Nodes(wkInd,:);
            
            % current subcell displacements
			wkU = result.uu(edof);
            
            % loop over Gauss points
            for ig = 1:size(W,1)
                
                % shape functions
                [N,dNdxi] = lagrange_basis('T6',Q(ig,:));
                
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
				% enNum = enNum + (Eh-E)'*model.Cmat*(Eh-E)*W(ig);
            	enNum = enNum + (E-Eh)'*model.Cmat*(E-Eh)*W(ig);
				enDenom = enDenom + E'*model.Cmat*E*W(ig);
            end
        else
        	% Gauss points
        	[W,Q] = quadrature(2,'GAUSS',2);

        	% re-numbering
            for ic = 1:nc
                nsf = length(model.supp{neighbour(ic)});
                if ic == 1
                    nodL = model.supp{neighbour(ic)};
                    nn = nsf;
                else
                    i0 = 0;
                    for jj = 1:nsf
                        nod = model.supp{neighbour(ic)}(jj);
                        flag = 0;
                        for j = 1:nn
                            if nodL(j) == nod
                                flag = 1;
                                break;
                            end
                        end
                        if flag == 0
                            i0 = i0 + 1;
                            nodL(nn+i0) = nod;
                        end
                    end
                    nn = nn + i0;
                end
            end
            
            % current subcell connectivity
            if targetEdge(ivo,5) == 1
                node_sc = [1 7 2 3 9 8 5 6 4];
            elseif targetEdge(ivo,5) == 2
                node_sc = [1 2 7 3 4 8 9 6 5];
            else
                node_sc = [1 2 3 7 4 5 8 9 6];
            end
            wkInd = nodL(node_sc); 
            
            edof = zeros(1,2*length(wkInd));
            edof(1:2:end) = 2*wkInd - 1;
            edof(2:2:end) = 2*wkInd;

            % subcell coordinates
            wkX = model.Nodes(wkInd,:);
            
            % current subcell displacements
            wkU = result.uu(edof);
            
            % shape functions at internal Gauss points
            for ig = 1:size(W,1)
                % shape functions
                [N,dNdxi] = lagrange_basis('Q9',Q(ig,:));

                % Jacobian
                J0 = dNdxi'*wkX;
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

                % strains
                Eh = [];
                Eh = Bmat*wkU;

                % compute anlaytical solution for displacement and strains
                xygp = N'*wkX;
                [U,E,S] = analytic_solution(model,xygp);

                % error in displacement
                edNum = edNum + ((U(1,1)-Uh(1,1))^2 + (U(2,1)-Uh(2,1))^2)*W(ig);
                edDenom = edDenom + (U(1,1)^2 + U(2,1)^2)*W(ig);

                % error in energy
                enNum = enNum + (Eh-E)'*model.Cmat*(Eh-E)*W(ig);
                enDenom = enDenom + E'*model.Cmat*E*W(ig);
            end
        end
        clear nodL wkInd wkX wkU edof xygp Q W
    end

    result.L2norm = sqrt(edNum/edDenom);
    result.H1norm = sqrt(enNum/enDenom);
end
%%
function [U,E,S] = analytic_solution(model,xy)
	%-----------------------------------------------------------------------
	%  analytical solution:
	%  ux=Py/(6EI)*[(6L-3x)*x+(2+poisson)*(y^2-D^2/4)];
	%  uy=-P/(6EI)*[3*poisson*y^2*(L-x)+(4+5*poisson)*D^2*x/4+(3*L-x)*x^2];
	%-----------------------------------------------------------------------
	P = model.P;
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