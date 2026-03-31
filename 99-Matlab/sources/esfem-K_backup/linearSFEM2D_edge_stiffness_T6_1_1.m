function model = linearSFEM2D_edge_stiffness_T6_1_1(model)
	% Edge-based smoothed finite element method
	% 	T6 element
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
    	
    % Loop over edges
    for ivo = 1:size(targetEdge,1)
    	% current target edge
        if (targetEdge(ivo,end-2) == 0)
            neighbour = targetEdge(ivo,3);
        else 
            neighbour = targetEdge(ivo,3:end-2); 
        end
        nc = length(neighbour);
        
        % loop over sub-cells
        for ic = 1:nc
            
            % support domain numbering
            nsf = length(model.supp{neighbour(ic)});
            
            % current element connectivity
            wkInd = model.Elements(neighbour(ic),:);

            % current element coordinates
            wkX = model.Nodes(wkInd,:);
            
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
            for is = 1:size(bound,1)
                if nc == 1; bxy = zeros(3,size(wkX,1)); else; bxy = zeros(4,size(wkX,1)); end
                X = wkX(bound(is,:),:);

                % loop over Gauss points on boundary
                for ig = 1:ng
                    % map 1D Gauss point to 2D
                    [Ng,dNdxi] = lagrange_basis('L3',Qb(ig));
                    xieta_gp = Ng'*X;
                    detJ = norm(dNdxi'*X);
                    
                    if model.flag == 1; plot(xieta_gp(1),xieta_gp(2),'bs'); end

                    % shape functions on boundary
                    N_T = zeros(size(wkX,1),1);
                    N_T = getSerendipityShapeFunc_lagrange('T6',wkX,xieta_gp);
                    
                    N = N_T'*detJ*Wb(ig);
                    if nc == 1
                        bxy = bxy + N.*[1; xieta_gp(1); xieta_gp(2)];
                    else
                        bxy = bxy + N.*[1; xieta_gp(1); xieta_gp(2); xieta_gp(1)*xieta_gp(2)];
                    end
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
            else
                i0 = 0;
                for jj = 1:nsf
                    nod = model.supp{neighbour(ic)}(jj);
                    flag = 0;
                    for j = 1:nn
                        if nodL(j) == nod
                            Fx(:,j) = Fx(:,j) + fx(:,jj);
                            Fy(:,j) = Fy(:,j) + fy(:,jj);
                            flag = 1;
                            break;
                        end
                    end
                    if flag == 0
                        i0 = i0 + 1;
                        nodL(nn+i0) = nod;
                        Fx(:,nn+i0) = fx(:,jj);
                        Fy(:,nn+i0) = fy(:,jj);
                    end
                end
                nn = nn + i0;
            end
        end
        % re-numbering
        if nc == 1            
            node_sc = [1 2 3 4 5 6];
            % 
            element = 'T6';
            %
            [Wi,Qi] = quadrature(ng,'TRIANGULAR',2);
        else
            
            if (targetEdge(ivo,5) == 1)
                node_sc = [1 7 2 3 9 8 5 6 4];
            elseif (targetEdge(ivo,5) == 2)
                node_sc = [1 2 7 3 4 8 9 6 5];
            else
                node_sc = [1 2 3 7 4 5 8 9 6]; 
            end
            %
            element = 'Q9';
            % 
            [Wi,Qi] = quadrature(ng,'GAUSS',2);
        end
        nodL = nodL(node_sc);
        gcoord = model.Nodes(nodL,:);
        Fx = Fx(:,node_sc);
        Fy = Fy(:,node_sc);
        
        % W matrix at interior Gauss points
        Ni = zeros(3,size(gcoord,1));
        if nc == 1
            Wmat = zeros(3,size(Wi,1)); 
        else
            Wmat = zeros(4,size(Wi,1)); 
        end
        
        % loop over interior Gauss points
        detJ0 = zeros(size(Wi,1),1);
        mQ = zeros(size(Wi,1),2);
        mW = zeros(size(Wi,1),1);
        for ig = 1:size(Wi,1)
            % location of interior Gauss points
            if nc == 1
                [N1,dN1dxi] = lagrange_basis('T3',Qi(ig,:));
                detJ0(ig) = det(dN1dxi'*gcoord(1:3,:));
                mQ(ig,:) = N1'*gcoord(1:3,:);
            else
                [N1,dN1dxi] = lagrange_basis('Q4',Qi(ig,:));
                detJ0(ig) = det(dN1dxi'*gcoord(1:4,:));
                mQ(ig,:) = N1'*gcoord(1:4,:);
            end
            mW(ig) = Wi(ig)*detJ0(ig);
                
            if model.flag == 1; plot(mQ(ig,1),mQ(ig,2),'r*'); end

            % shape functions & W matrix
            N = getSerendipityShapeFunc_lagrange(element,gcoord,mQ(ig,:));
            if nc == 1
                Ni = Ni + (N*mW(ig))'.*[1; 0; 0];
                Wmat(:,ig) = mW(ig)*[1; mQ(ig,1); mQ(ig,2)];
            else
                Ni = Ni + (N*mW(ig))'.*[1; mQ(ig,1); mQ(ig,2)]; 
                Wmat(:,ig) = mW(ig)*[1; mQ(ig,1); mQ(ig,2); mQ(ig,1)*mQ(ig,2)];
            end
        end
        
        if nc == 1
            %
            Fx(2,:) = Fx(2,:) - Ni(1,:);
            %
            Fy(3,:) = Fy(3,:) - Ni(1,:);
        else
            %
            Fx(2,:) = Fx(2,:) - Ni(1,:);
            Fx(4,:) = Fx(4,:) - Ni(3,:);
            % 
            Fy(3,:) = Fy(3,:) - Ni(1,:);
            Fy(4,:) = Fy(4,:) - Ni(2,:);
        end
        
        % get derivatives of basis functions
        dx = Wmat\Fx;
        dy = Wmat\Fy;

        % compute smoothed stiffness matrix
        nndof = 2*length(nodL);
        edof = zeros(1,nndof);
        edof(1:2:end) = 2*nodL - 1;
        edof(2:2:end) = 2*nodL;
        for ig = 1:size(mW,1)
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