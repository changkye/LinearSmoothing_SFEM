function model = linearSFEM2D_edge_stiffness_T3_2(model)
	% Edge-based smoothed finite element method
	% 	T3 element with linear smoothing funcitons
    %
    % Changkye Lee

    % Initialisation
    K = sparse(2*size(model.Nodes,1),2*size(model.Nodes,1));

    % target edge and sharing elements
    model = getTargetEdge(model);
    targetEdge = model.targetEdge;
    
    % boundary Gauss points
    ng = 2;
    [Wb,Qb] = quadrature(ng,'GAUSS',1);
    
 %    % interior Gauss points
	% [Wi,Qi] = quadrature(2,'TRIANGULAR',2);
    
    % Loop over edges
    for ivo = 1:size(targetEdge,1)
    	% current target edge
        if (targetEdge(ivo,end-2) == 0); neighbour = targetEdge(ivo,3);
        else; neighbour = targetEdge(ivo,3:end-2); end
        nc = length(neighbour);
        
        % loop over sub-cells
        if nc == 1

            [Wi,Qi] = quadrature(2,'TRIANGULAR',2);

            wkInd = model.Elements(neighbour,:);
            wkX = model.Nodes(wkInd,:);
            
            % subcell boundary length
            side = cal_side(wkX(:,1),wkX(:,2));

            % outward normal vectors
            [nx,ny] = cal_nx_ny(wkX(:,1),wkX(:,2),side);

            % initialise
            fx = zeros(3,size(wkX,1)); fy = zeros(3,size(wkX,1));
            detJ0 = zeros(size(Wi,1),1);
            mR0 = zeros(size(Wi,1),2);
            Ni = zeros(3,size(wkX,1));
            Wmat = zeros(3);
            
            bound = [1 2; 2 3; 3 1];
            % loop over subcell boundary
            for is = 1:size(bound,1)

                % location of internagl Gauss points
                [N1,dN1dxi] = lagrange_basis('T3',Qi(is,:));
                detJ0(is) = det(dN1dxi'*wkX);
                mR0(is,:) = N1'*wkX;

                % shape functions & W matrix
                xieta = cal_xieta('T3',mR0(is,:),wkX);
                N = lagrange_basis('T3',xieta);
                Ni(1,:) = Ni(1,:) + (N*Wi(is)*detJ0(is))';

                Wmat(:,is) = Wmat(:,is) + Wi(is)*detJ0(is)*[1; mR0(is,1); mR0(is,2)];

                bxy = zeros(3,size(wkX,1));
                X = wkX(bound(is,:),:);

                % loop over Gauss points on boundary
                for ig = 1:ng
                    % map 1D Gauss point to 2D
                    [Ng,dNdxi] = lagrange_basis('L2',Qb(ig));
                    xieta_gp = Ng'*X;
                    detJ = norm(dNdxi'*X);

                    % shape functions on boundary
                    N_T = zeros(size(wkX,1),1);
                    xieta = cal_xieta('T3',xieta_gp,wkX);
                    N_T = lagrange_basis('T3',xieta);
                    
                    N = N_T'*detJ*Wb(ig);
                    bxy = bxy + [N; N*xieta_gp(1); N*xieta_gp(2)];
                    
                end
                fx = fx + nx(is)*bxy;
                fy = fy + ny(is)*bxy;
            end
            fx(2,:) = fx(2,:) - Ni(1,:);
            fy(3,:) = fy(3,:) - Ni(1,:);
            
            dx = Wmat\fx;
            dy = Wmat\fy;
            
        

            % compute smoothed stiffness matrix
            nndof = 2*length(nodL);
            edof = zeros(1,nndof);
            edof(1:2:end) = 2*nodL - 1;
            edof(2:2:end) = 2*nodL;
            for ig = 1:size(Wi,1)
                % smoothed strain-displacement matrix
                Bmat = zeros(3,nndof);
                Bmat(1,1:2:end) = dx(ig,:);
                Bmat(2,2:2:end) = dy(ig,:);
                Bmat(3,1:2:end) = dy(ig,:);
                Bmat(3,2:2:end) = dx(ig,:);

                % smoothed stiffness matrix
                K(edof,edof) = K(edof,edof) + Bmat'*model.Cmat*Bmat*Wi(ig)*detJ0(ig);
            end
            clear N nodL
    end
    model.K = K;
end