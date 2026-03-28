function model = linearSFEM2D_edge_stiffness_T6_1_2(model)
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
            
            % if nc == 1; ng = 2;  else; ng = 2; end
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
            
            % subcell boundary length
            side = cal_side(wkX(1:3,1),wkX(1:3,2));

            % outward normal vectors
            [nx,ny] = cal_nx_ny(wkX(1:3,1),wkX(1:3,2),side);

            % initialise
            if nc == 1
                fx = zeros(3,size(wkX,1)); fy = zeros(3,size(wkX,1));
            else
                fx = zeros(4,size(wkX,1)); fy = zeros(4,size(wkX,1)); 
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
                    
                    if model.flag==1; plot(xieta_gp(1),xieta_gp(2),'bs'); end

                    % shape functions on boundary
                    N_T = zeros(size(wkX,1),1);
%                     N_T = getSerendipityShapeFunc_lagrange('T6',wkX,xieta_gp);
                    N_T = wachspress2D(wkX,xieta_gp);
                    
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
            
        if nc == 1       
            % re-numbering
            node_sc = [1 2 3 4 5 6];
            nodL = nodL(node_sc);
            % 
            [W,Q] = quadrature(ng,'TRIANGULAR',2);
            element = {'T3'; 'T6'};
            gcoord = model.Nodes(nodL,:);
            %
            Fx = Fx(:,node_sc);
            Fy = Fy(:,node_sc);
            %
            NPE = 3;
        else
            % re-numbering
            if targetEdge(ivo,5) == 1
                node_sc = [1 7 2 3 9 8 5 6 4];
            elseif targetEdge(ivo,5) == 2
                node_sc = [1 2 7 3 4 8 9 6 5];
            else 
                node_sc = [1 2 3 7 4 5 8 9 6]; 
            end
            nodL = nodL(node_sc);
            %
            [W,Q] = quadrature(ng,'GAUSS',2);
            element = {'Q4'; 'Q9'};
            gcoord = model.Nodes(nodL,:);
            %
            Fx = Fx(:,node_sc);
            Fy = Fy(:,node_sc);
            %
            NPE = 4;
        end
        
        % W matrix at interior Gauss points
        detJ0 = zeros(size(W,1),1);
        mQ = zeros(size(W,1),2);
        mW = zeros(size(W,1),1);
        Ni = zeros(3,size(gcoord,1));
        if nc == 1; Wmat = zeros(3,size(W,1)); else; Wmat = zeros(4,size(W,1)); end 
        for ig = 1:size(W,1)
            % location of internagl Gauss points
            [N0,dNdxi] = lagrange_basis(element{1},Q(ig,:));
            detJ0(ig) = det(dNdxi'*gcoord(1:NPE,:));
            mQ(ig,:) = N0'*gcoord(1:NPE,:);
            mW(ig) = W(ig)*detJ0(ig);
            
            if model.flag==1; plot(mQ(ig,1),mQ(ig,2),'r*'); end

            % shape functions & W matrix
%             N = getSerendipityShapeFunc_lagrange(element{2},gcoord,mQ(ig,:));
            N = wachspress2D(gcoord,mQ(ig,:));
            if nc == 1
                Ni = Ni + (N*mW(ig))'.*[1; 0; 0];
                Wmat(:,ig) = mW(ig)*[1; mQ(ig,1); mQ(ig,2)];
            else
                Ni = Ni + (N*mW(ig))'.*[1; mQ(ig,1); mQ(ig,2)]; 
                Wmat(:,ig) = mW(ig)*[1; mQ(ig,1); mQ(ig,2); mQ(ig,1)*mQ(ig,2)];
            end
        end

        if nc == 1
            Fx(2,:) = Fx(2,:) - Ni(1,:);
            %
            Fy(3,:) = Fy(3,:) - Ni(1,:);
        else
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