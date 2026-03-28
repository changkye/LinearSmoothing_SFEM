function model = linearSFEM2D_edge_stiffness_T6(model)
	% Edge-based smoothed finite element method
	% 	T6 element
    %   linear smoothing function for T6 
    %   bilinear smoothing funciton for Q9
    %
    % Changkye Lee

    % Initialisation
    K = sparse(2*size(model.Nodes,1),2*size(model.Nodes,1));

    % Gauss points on subcell boundary
    ng = 2;
    [Wb,Qb] = quadrature(ng,'GAUSS',1);

    % target edge and sharing elements
    model = getTargetEdge(model);
    targetEdge = model.targetEdge;
    
    % Loop over edges
    for ivo = 1:size(targetEdge,1)
    	% current target edge
        if targetEdge(ivo,end-2) == 0; neighbour = targetEdge(ivo,3);
        else; neighbour = targetEdge(ivo,3:end-2); end
        nc = length(neighbour);
        
        % loop over sub-cells
        if nc == 1
            % Gauss points on subcell boundary
            ng = 2;
            [Wb,Qb] = quadrature(ng,'GAUSS',1);
            
            % internal Gauss points
            [Wi,Qi] = quadrature(2,'TRIANGULAR',2);

            % current element connectivity
            wkInd = model.Elements(neighbour,:);

            % current element coordinates
            wkX = model.Nodes(wkInd,:);
            edof = zeros(1,2*length(wkInd));
            edof(1:2:end) = 2*wkInd - 1;
            edof(2:2:end) = 2*wkInd;

            % shape functions at internal Gauss points
            Ni = zeros(1,size(wkX,1));
            Wmat = zeros(3);
            
            % loop over internal Gauss points
            detJ0 = zeros(size(Wi,1),1);
            mR = zeros(size(Wi,1),2);
            N = zeros(size(wkX,1),1);
            for ig = 1:size(Wi,1)
                % location of internagl Gauss points
                [N1,dN1dxi] = lagrange_basis('T3',Qi(ig,:));
                detJ0(ig) = det(dN1dxi'*wkX(1:3,:));
                mR(ig,:) = N1'*wkX(1:3,:);

                % shape functions
                N = getSerendipityShapeFunc_lagrange(model.elemType,wkX,mR(ig,:));
                Ni = Ni + (N*Wi(ig)*detJ0(ig))';

                % W matrix
                Wmat(:,ig) = Wi(ig)*detJ0(ig)*[1; mR(ig,1); mR(ig,2)];
            end
            if model.flag==1; plot(mR(:,1),mR(:,2),'r*'); end
            
            % subcell boundary length
            side = cal_side(wkX(1:3,1),wkX(1:3,2));

            % outward normal vectors
            [nx,ny] = cal_nx_ny(wkX(1:3,1),wkX(1:3,2),side);

            % initialise
            fx = zeros(3,size(wkX,1)); fy = zeros(3,size(wkX,1));
            bound = [1 2 4; 2 3 5; 3 1 6];  

            % loop over subcell boundary
            for is = 1:size(bound,1)
                bxy = zeros(3,size(wkX,1));
                X = wkX(bound(is,:),:);

                % loop over Gauss points on boundary
                for ig = 1:ng
                    % map 1D Gauss point to 2D
                    [Ng,dNdxi] = lagrange_basis('L3',Qb(ig));
                    xieta_gp = Ng'*X;
                    if model.flag==1; plot(xieta_gp(1),xieta_gp(2),'b^'); end
                    detJ = norm(dNdxi'*X);

                    % shape functions on boundary
                    N_T = zeros(size(wkX,1),1);
                    N_T = getSerendipityShapeFunc_lagrange(model.elemType,wkX,xieta_gp);
                    
                    bxy = bxy + (N_T'*detJ*Wb(ig)).*[1; xieta_gp(1); xieta_gp(2)];
                end
                fx = fx + nx(is)*bxy;
                fy = fy + ny(is)*bxy;
            end
            fx(2,:) = fx(2,:) - Ni;
            fy(3,:) = fy(3,:) - Ni;

            % get derivatives of basis functions
            dx = Wmat\fx;
            dy = Wmat\fy;

            % compute smoothed stiffness matrix
            for ig = 1:size(Wi,1)
                % smoothed strain-displacement matrix
                Bmat = zeros(3,2*length(wkInd));
                Bmat(1,1:2:end) = dx(ig,:);
                Bmat(2,2:2:end) = dy(ig,:);
                Bmat(3,1:2:end) = dy(ig,:);
                Bmat(3,2:2:end) = dx(ig,:);

                % smoothed stiffness matrix
                K(edof,edof) = K(edof,edof) + Bmat'*model.Cmat*Bmat*Wi(ig)*detJ0(ig);
            end
            clear side nx ny

        else
            % Gauss points on subcell boundary
            ng = 2; 
            [Wb,Qb] = quadrature(ng,'GAUSS',1);
            
            % internal Gauss points
            [Wi,Qi] = quadrature(2,'GAUSS',2);

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

            % current subcell coordinates
            wkX = model.Nodes(wkInd,:);
            
            % shape functions at internal Gauss points
            Ni = zeros(3,size(wkX,1));
            Wmat = zeros(4);
            
            % loop over internal Gauss points
            detJ0 = zeros(size(Wi,1),1);
            mR = zeros(size(Wi,1),2);
            N = zeros(size(wkX,1),1);
            for ig = 1:size(Wi,1)
                % location of internal Gauss points
                [N1,dN1dxi] = lagrange_basis('Q4',Qi(ig,:));
                detJ0(ig) = det(dN1dxi'*wkX(1:4,:));
                mR(ig,:) = N1'*wkX(1:4,:);
                if model.flag == 1
                    plot(mR(ig,1),mR(ig,2),'r*');
                end

                % shape functions
                N = getSerendipityShapeFunc_lagrange('Q9',wkX,mR(ig,:));
                Ni = Ni + (N*Wi(ig)*detJ0(ig))'.*[1; mR(ig,2); mR(ig,1)];
                
                % W matrix
                Wmat(:,ig) = Wi(ig)*detJ0(ig)*[1; mR(ig,1); mR(ig,2); mR(ig,1)*mR(ig,2)];
            end 

            % subcell boundary length
            side = cal_side(wkX(1:4,1),wkX(1:4,2));

            % outward normal vectors
            [nx,ny] = cal_nx_ny(wkX(1:4,1),wkX(1:4,2),side);

            % initialise
            fx = zeros(4,size(wkX,1)); fy = zeros(4,size(wkX,1));
            bound = [1 2 5; 2 3 6; 3 4 7; 4 1 8];

            % loop over subcell boundary
            for is = 1:size(bound,1)
                bxy = zeros(4,size(wkX,1));
                X = wkX(bound(is,:),:);

                % loop over Gauss points on boundary
                for ig = 1:ng
                    % map 1D Gauss point to 2D
                    [Ng,dNdxi] = lagrange_basis('L3',Qb(ig));
                    xieta_gp = Ng'*X;
                    if model.flag==1; plot(xieta_gp(1),xieta_gp(2),'b^'); end
                    detJ = norm(dNdxi'*X);

                    % shape functions on boundary
                    N_T = zeros(size(wkX,1),1);
                    N_T = getSerendipityShapeFunc_lagrange('Q9',wkX,xieta_gp);
                    bxy = bxy + (N_T'*detJ*Wb(ig)).*[1; xieta_gp(1); xieta_gp(2); xieta_gp(1)*xieta_gp(2)];
                end
                fx = fx + nx(is)*bxy;
                fy = fy + ny(is)*bxy;
            end
            fx(2,:) = fx(2,:) - Ni(1,:);
            fx(4,:) = fx(4,:) - Ni(2,:);
            % 
            fy(3,:) = fy(3,:) - Ni(1,:);
            fy(4,:) = fy(4,:) - Ni(3,:);

            % get derivatives of basis functions
            dx = Wmat\fx;
            dy = Wmat\fy;

            % compute smoothed stiffness matrix
            for ig = 1:size(Wi,1)
                % smoothed strain-displacement matrix
                Bmat = zeros(3,2*length(wkInd));
                Bmat(1,1:2:end) = dx(ig,:);
                Bmat(2,2:2:end) = dy(ig,:);
                Bmat(3,1:2:end) = dy(ig,:);
                Bmat(3,2:2:end) = dx(ig,:);

                % smoothed stiffness matrix
                K(edof,edof) = K(edof,edof) + Bmat'*model.Cmat*Bmat*Wi(ig)*detJ0(ig);
            end
            clear side nx ny nodL wkInd wkX X  
        end
    end

    model.K = K;
end