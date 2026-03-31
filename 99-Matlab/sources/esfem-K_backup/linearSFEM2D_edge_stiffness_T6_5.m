function model = linearSFEM2D_edge_stiffness_T6_5(model)
	% Edge-based smoothed finite element method
	% 	T6 element
    %
    % Changkye Lee

    % Initialisation
    K = sparse(2*size(model.Nodes,1),2*size(model.Nodes,1));

    % target edge and sharing elements
    model = getTargetEdge(model);
    targetEdge = model.targetEdge;
    
    % Loop over edges
    for ivo = 1:size(targetEdge,1)
    	% current target edge
        if targetEdge(ivo,end-2) == 0
            neighbour = targetEdge(ivo,3);
        else
            neighbour = targetEdge(ivo,3:end-2); 
        end
        nc = length(neighbour);
        
        % loop over sub-cells
        for ic = 1:nc
            
            ng = 2;
            % boundary Gauss points
            [Wb,Qb] = quadrature(ng,'GAUSS',1);
    
            % interior Gauss points
            [Wi,Qi] = quadrature(ng,'TRIANGULAR',2);
            
            % support domain numbering
            nsf = length(model.supp{neighbour(ic)});
            
            % current element connectivity
            wkInd = model.Elements(neighbour(ic),:);

            % current element coordinates
            wkX = model.Nodes(wkInd,:);
            
            % loop over interior Gauss points
            detJ0 = zeros(size(Wi,1),1);
            mQ0 = zeros(size(Wi,1),2);
            mW0 = zeros(size(Wi,1),1);
            Ni = zeros(3,size(wkX,1));
            for ig = 1:size(Wi,1)
                % location of interior Gauss points
                [N1,dN1dxi] = lagrange_basis('T3',Qi(ig,:));
                detJ0(ig) = det(dN1dxi'*wkX(1:3,:));
                mQ0(ig,:) = N1'*wkX(1:3,:);
                mW0(ig) = Wi(ig)*detJ0(ig);

                if model.flag==1; plot(mQ0(ig,1),mQ0(ig,2),'ko'); end

                % shape functions & W matrix
                N = getSerendipityShapeFunc_lagrange('T6',wkX,mQ0(ig,:));
                if nc == 1
                    Ni = Ni + (N*mW0(ig))'.*[1; 0; 0];
                else
                    Ni = Ni + (N*mW0(ig))'.*[1; mQ0(ig,1); mQ0(ig,2)]; 
                end
            end
            
            % subcell boundary length
            side = cal_side(wkX(1:3,1),wkX(1:3,2));

            % outward normal vectors
            [nx,ny] = cal_nx_ny(wkX(1:3,1),wkX(1:3,2),side);

            % initialise
            if nc == 1
                fx = zeros(3,size(wkX,1)); 
                fy = zeros(3,size(wkX,1));
            else
                fx = zeros(4,size(wkX,1)); 
                fy = zeros(4,size(wkX,1)); 
            end
            bound = [1 2 4; 2 3 5; 3 1 6];

            % loop over subcell boundary
            M1 = 0.0; M2 = 0.0; M12 = 0.0;
            for is = 1:size(bound,1)
                if nc == 1; bxy = zeros(3,size(wkX,1)); else; bxy = zeros(4,size(wkX,1)); end
                X = wkX(bound(is,:),:);

                % loop over Gauss points on boundary
                for ig = 1:ng
                    % map 1D Gauss point to 2D
                    [Ng,dNdxi] = lagrange_basis('L3',Qb(ig));
                    xieta_gp = Ng'*X;
                    detJ = norm(dNdxi'*X);
                    
                    if model.flag==1; plot(xieta_gp(1),xieta_gp(2),'bs'); end

                    % shape functions on boundary
                    N_T = zeros(size(wkX,1),1);
                    N_T = getSerendipityShapeFunc_lagrange('T6',wkX,xieta_gp);
                    
                    N = N_T'*detJ*Wb(ig);
                    if nc == 1
                        bxy = bxy + N.*[1; xieta_gp(1); xieta_gp(2)];
                    else
                        bxy = bxy + N.*[1; xieta_gp(1); xieta_gp(2); xieta_gp(1)*xieta_gp(2)];
                    end
                    
                    M1 = M1 + xieta_gp(1)*detJ*Wb(ig);
                    M2 = M2 + xieta_gp(2)*detJ*Wb(ig);
                    M12 = M12 + xieta_gp(1)*xieta_gp(2)*detJ*Wb(ig);
                end
                fx = fx + nx(is)*bxy;
                fy = fy + ny(is)*bxy;
            end
            
            % re-numbering
            if ic == 1
                nodL = model.supp{neighbour(ic)};
                nn = nsf;
                Fx = fx;
                Fy = fy;
                Nj = Ni;
            else
                i0 = 0;
                for jj = 1:nsf
                    nod = model.supp{neighbour(ic)}(jj);
                    flag = 0;
                    for j = 1:nn
                        if nodL(j) == nod
                            Fx(:,j) = Fx(:,j) + fx(:,jj);
                            Fy(:,j) = Fy(:,j) + fy(:,jj);
                            Nj(:,j) = Nj(:,j) + Ni(:,jj);
                            flag = 1;
                            break;
                        end
                    end
                    if flag == 0
                        i0 = i0 + 1;
                        nodL(nn+i0) = nod;
                        Fx(:,nn+i0) = fx(:,jj);
                        Fy(:,nn+i0) = fy(:,jj);
                        Nj(:,nn+i0) = Ni(:,jj);
                    end
                end
                nn = nn + i0;
            end
        end
            
        if nc == 1
            Fx(2,:) = Fx(2,:) - Nj(1,:);
            %
            Fy(3,:) = Fy(3,:) - Nj(1,:);
            
            % re-numbering
            node_sc = [1 2 3 4 5 6];
            nodL = nodL(node_sc);
            % 
            [W,Q] = quadrature(ng,'TRIANGULAR',2);
            element = 'T3';
            gcoord = model.Nodes(nodL(1:3),:);
            %
            Fx = Fx(:,node_sc);
            Fy = Fy(:,node_sc);
        else
            Fx(2,:) = Fx(2,:) - Nj(1,:);
            Fx(4,:) = Fx(4,:) - Nj(3,:);
            % 
            Fy(3,:) = Fy(3,:) - Nj(1,:);
            Fy(4,:) = Fy(4,:) - Nj(2,:);
            
            % re-numbering
            if (targetEdge(ivo,5) == 1)
                node_sc = [1 7 2 3 9 8 5 6 4];
            elseif (targetEdge(ivo,5) == 2)
                node_sc = [1 2 7 3 4 8 9 6 5];
            else
                node_sc = [1 2 3 7 4 5 8 9 6]; 
            end
            nodL = nodL(node_sc);
            %
            [W,Q] = quadrature(ng,'GAUSS',2);
            element = 'Q4';
            gcoord = model.Nodes(nodL(1:4),:);
            %
            Fx = Fx(:,node_sc);
            Fy = Fy(:,node_sc);
        end
        M0 = polyarea(gcoord(:,1),gcoord(:,2));
        Fx(1,:) = Fx(1,:)/M0; Fy(1,:) = Fy(1,:)/M0;
        Fx(2,:) = Fx(2,:)/M1; Fy(2,:) = Fy(2,:)/M1;
        Fx(3,:) = Fx(3,:)/M2; Fy(3,:) = Fy(3,:)/M2;
        if nc == 2
            Fx(4,:) = Fx(4,:)/M12;
            Fy(4,:) = Fy(4,:)/M12;
        end
        
        % W matrix at interior Gauss points
        detJ = zeros(size(W,1),1);
        mQ = zeros(size(W,1),2);
        mW = zeros(size(W,1),1);
        if nc == 1; Wmat = zeros(3,size(W,1)); else; Wmat = zeros(4,size(W,1)); end 
        for ig = 1:size(W,1)
            % location of internagl Gauss points
            [N,dNdxi] = lagrange_basis(element,Q(ig,:));
            detJ(ig) = det(dNdxi'*gcoord);
            mQ(ig,:) = N'*gcoord;
            mW(ig) = W(ig)*detJ(ig);
            
            if model.flag==1; plot(mQ(ig,1),mQ(ig,2),'r*'); end

            % W matrix
            if nc == 1
                Wmat(:,ig) = mW(ig)*[1; mQ(ig,1); mQ(ig,2)];
            else
                Wmat(:,ig) = mW(ig)*[1; mQ(ig,1); mQ(ig,2); mQ(ig,1)*mQ(ig,2)]; 
            end
        end

        % get derivatives of basis functions
        dx = Wmat\Fx;
        dy = Wmat\Fy;

        % compute smoothed stiffness matrix
        nndof = 2*length(nodL);
        edof = zeros(1,nndof);
        edof(1:2:end) = 2*nodL - 1;
        edof(2:2:end) = 2*nodL;
        for ig = 1:size(W,1)
            % smoothed strain-displacement matrix
            Bmat = zeros(3,nndof);
            Bmat(1,1:2:end) = dx(ig,:);
            Bmat(2,2:2:end) = dy(ig,:);
            Bmat(3,1:2:end) = dy(ig,:);
            Bmat(3,2:2:end) = dx(ig,:);

            % smoothed stiffness matrix
            K(edof,edof) = K(edof,edof) + Bmat'*model.Cmat*Bmat*mW(ig);
        end
        clear N nodL
    end
    model.K = K;
end