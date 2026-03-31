function model = linearSFEM2D_edge_stiffness_T6_0(model)
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
        if targetEdge(ivo,end-2) == 0
            neighbour = targetEdge(ivo,3);
        else
            neighbour = targetEdge(ivo,3:end-2); 
        end
        nc = length(neighbour);
        
        % loop over sub-cells
        if nc == 1
            
            % internal Gauss points
            [W,Q] = quadrature(2,'TRIANGULAR',2);

            % current element connectivity
            wkInd = model.Elements(neighbour,:);

            % current element coordinates
            wkX = model.Nodes(wkInd,:);
            
            % create subcell
            tmp = [wkX; mean(wkX)];
            gcoord = [tmp; mean(tmp([1,7],:)); mean(tmp([2,7],:)); mean(tmp([3,7],:))];
            clear tmp;

            % subcell coordinates & connectivities
            node_sc = [1 2 7 4 9 8; 2 3 7 5 10 9; 3 1 7 6 8 10];
            if targetEdge(ivo,end-1) == 1
                X = gcoord(node_sc(1,:),:);
            elseif targetEdge(ivo,end-1) == 2
                X = gcoord(node_sc(2,:),:);
            elseif targetEdge(ivo,end-1) == 3
                X = gcoord(node_sc(3,:),:);
            end
            if model.flag==1; trimesh(node_sc,gcoord(:,1),gcoord(:,2)); end

            % side length of the smoothing domain
            side = cal_side(X(1:3,1),X(1:3,2));
            
            % normal vectors
            [nx,ny] = cal_nx_ny(X(1:3,1),X(1:3,2),side);
            
            % initialise
            fx = zeros(3,size(wkX,1)); 
            fy = zeros(3,size(wkX,1));
            bound = [1 2 4; 2 3 5; 3 1 6];

            % loop over subcell boundary
            for is = 1:size(bound,1)
                bxy = zeros(3,size(wkX,1));
                subX = X(bound(is,:),:);

                % loop over Gauss points on boundary
                for ig = 1:ng
                    N_T = zeros(size(gcoord,1),1);
                    % map 1D Gauss point to 2D
                    [Ng,dNdxi] = lagrange_basis('L3',Qb(ig));
                    xieta_gp = Ng'*subX;
                    detJ0 = norm(dNdxi'*subX);
                    
                    if model.flag==1; plot(xieta_gp(1),xieta_gp(2),'bs'); end

                    % shape functions on boundary
                    xieta = cal_xieta('T3',xieta_gp,wkX);        
                    N_T = lagrange_basis('T6',xieta);
%                     N_T = wachspress2D(X,xieta_gp);
                    
                    N = N_T'*detJ0*Wb(ig);
                    bxy = [N; N*xieta_gp(1); N*xieta_gp(2)];
                    fx = fx + nx(is)*bxy;
                    fy = fy + ny(is)*bxy;
                end
            end

            detJ = zeros(size(W,1),1);
            mQ = zeros(size(W,1),2);
            mW = zeros(size(W,1),1);
            Ni = zeros(1,size(wkX,1));
            Wmat = zeros(3,size(W,1));
            
            % loop over interior Gauss points
            for ig = 1:size(W,1)
                % location of internagl Gauss points
                [N,dNdxi] = lagrange_basis('T3',Q(ig,:));
                detJ(ig) = det(dNdxi'*X(1:3,:));
                mQ(ig,:) = N'*X(1:3,:);
                mW(ig) = W(ig)*detJ(ig);
            
                if model.flag==1; plot(mQ(ig,1),mQ(ig,2),'r*'); end

                % interior shape functions
                xieta = cal_xieta('T3',mQ(ig,:),wkX);           
                N_T = lagrange_basis('T6',xieta);
                
                Ni = Ni + N_T'*mW(ig);

                % W matrix
                Wmat(:,ig) = mW(ig)*[1; mQ(ig,1); mQ(ig,2)];
            
            end
            % 
            fx(2,:) = fx(2,:) - Ni;
            fy(3,:) = fy(3,:) - Ni;

            % get derivatives of basis functions
            dx = Wmat\fx;
            dy = Wmat\fy;

            % compute smoothed stiffness matrix
            nodL = wkInd;
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
            
        else
            for ic = 1:nc

                % current element connectivity
                wkInd = model.Elements(neighbour(ic),:);

                 % current element coordinates
                wkX = model.Nodes(wkInd,:);
            
                % create subcell
                tmp = [wkX; mean(wkX)];
                gcoord = [tmp; mean(tmp([1,7],:)); mean(tmp([2,7],:)); mean(tmp([3,7],:))];
                clear tmp;

                % subcell coordinates & connectivities
                node_sc = [1 2 7 4 9 8; 2 3 7 5 10 9; 3 1 7 6 8 10];
                if targetEdge(ivo,ic+4) == 1
                    X = gcoord(node_sc(1,:),:);
                elseif targetEdge(ivo,ic+4) == 2
                    X = gcoord(node_sc(2,:),:);
                elseif targetEdge(ivo,ic+4) == 3
                    X = gcoord(node_sc(3,:),:);
                end
                if model.flag==1; trimesh(node_sc,gcoord(:,1),gcoord(:,2)); end

                % side length of the smoothing domain
                side = cal_side(X(1:3,1),X(1:3,2));
            
                % normal vectors
                [nx,ny] = cal_nx_ny(X(1:3,1),X(1:3,2),side);
            
                % initialise
                fx = zeros(4,size(wkX,1)); 
                fy = zeros(4,size(wkX,1));
                bound = [1 2 4; 2 3 5; 3 1 6];

                % loop over subcell boundary
                for is = 1:size(bound,1)
                    bxy = zeros(4,size(wkX,1));
                    subX = X(bound(is,:),:);

                    % loop over Gauss points on boundary
                    for ig = 1:ng
                        N_T = zeros(size(gcoord,1),1);
                        % map 1D Gauss point to 2D
                        [Ng,dNdxi] = lagrange_basis('L3',Qb(ig));
                        xieta_gp = Ng'*subX;
                        detJ0 = norm(dNdxi'*subX);
                    
                        if model.flag==1; plot(xieta_gp(1),xieta_gp(2),'bs'); end

                        % shape functions on boundary
                        xieta = cal_xieta('T3',xieta_gp,wkX);           
                        N_T = lagrange_basis('T6',xieta);
                    
                        N = N_T'*detJ0*Wb(ig);
                        bxy = N.*[1; xieta_gp(1); xieta_gp(2); xieta_gp(1)*xieta_gp(2)];
                        fx = fx + nx(is)*bxy;
                        fy = fy + ny(is)*bxy;
                    end
                end

                nsf = length(model.supp{neighbour(ic)});

                if ic == 1
                    nodL = wkInd;
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
                                flag = 1;
                                Fx(:,j) = Fx(:,j) + fx(:,jj);
                                Fy(:,j) = Fy(:,j) + fy(:,jj);
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

            % renumbering
            if targetEdge(ivo,5) == 1
                node_sc = [1 7 2 3 9 8 5 6 4];
                nodL = nodL(node_sc);
                tmp1 = model.Nodes(nodL([1,3,4]),:);
                tmp2 = model.Nodes(nodL([1,2,3]),:);
                X = [mean(tmp1); model.Nodes(nodL(1),:); mean(tmp2); model.Nodes(nodL(3),:)];
                X = [X; mean(X([1,2],:)); mean(X([2,3],:)); mean(X([3,4],:)); 
                    mean(X([4,1],:)); mean(model.Nodes(nodL([1,3]),:))];
                if model.flag==1; plot(X(:,1),X(:,2),'ko'); end
            elseif targetEdge(ivo,5) == 2
                node_sc = [1 2 7 3 4 8 9 6 5];
                nodL = nodL(node_sc);
                tmp1 = model.Nodes(nodL([1,2,4]),:);
                tmp2 = model.Nodes(nodL([2,3,4]),:);
                X = [mean(tmp1); model.Nodes(nodL(2),:); mean(tmp2); model.Nodes(nodL(4),:)];
                X = [X; mean(X([1,2],:)); mean(X([2,3],:)); mean(X([3,4],:)); 
                    mean(X([4,1],:)); mean(X([1,2,3,4],:))];
                if model.flag==1; plot(X(:,1),X(:,2),'ko'); end
            else
                node_sc = [1 2 3 7 4 5 8 9 6];
                nodL = nodL(node_sc);
                tmp1 = model.Nodes(nodL([1,3,4]),:);
                tmp2 = model.Nodes(nodL([1,2,3]),:);
                X = [mean(tmp1); model.Nodes(nodL(1),:); mean(tmp2); model.Nodes(nodL(3),:)];
                X = [X; mean(X([1,2],:)); mean(X([2,3],:)); mean(X([3,4],:)); 
                    mean(X([4,1],:)); mean(model.Nodes(nodL([1,3]),:))];
                if model.flag==1; plot(X(:,1),X(:,2),'ko'); end
            end
            % 
            Fx = Fx(:,node_sc); Fy = Fy(:,node_sc);
            
            % W matrix
            [W,Q] = quadrature(ng,'GAUSS',2);
            detJ = zeros(size(W,1),1);
            mQ = zeros(size(W,1),2);
            mW = zeros(size(W,1),1);
            Ni = zeros(3,size(X,1));
            Wmat = zeros(4,size(W,1));
            for ig = 1:size(W,1)
                [N0,dNdxi] = lagrange_basis('Q4',Q(ig,:));
                detJ(ig) = det(dNdxi'*X(1:4,:));
                mQ(ig,:) = N0'*X(1:4,:);
                mW(ig) = W(ig)*detJ(ig);
                                
                if model.flag==1; plot(mQ(ig,1),mQ(ig,2),'r*'); end
                
%                 N = getSerendipityShapeFunc_lagrange('Q9',model.Nodes(nodL,:),mQ(ig,:));
                xieta = cal_xieta('Q4',mQ(ig,:),model.Nodes(nodL,:));  
                N_T = lagrange_basis('Q9',xieta);
                
                Ni = Ni + (N_T*mW(ig))'.*[1; mQ(ig,1); mQ(ig,2)];
                Wmat(:,ig) = mW(ig)*[1; mQ(ig,1); mQ(ig,2); mQ(ig,1)*mQ(ig,2)];                
            end
            %
            Fx(2,:) = Fx(2,:) - Ni(1,:);
            Fx(4,:) = Fx(4,:) - Ni(3,:);
            %
            Fy(3,:) = Fy(3,:) - Ni(1,:);
            Fy(4,:) = Fy(4,:) - Ni(2,:);
            
            % get derivatives of basis funcitons
            dx = Wmat\Fx;
            dy = Wmat\Fy;
            
            % compute smoothed sitffness matrix
            nndof = 2*length(nodL);
            edof = zeros(1,nndof);
            edof(1:2:end) = 2*nodL - 1;
            edof(2:2:end) = 2*nodL;
            for ig = 1:size(W,1)
                % smoothed strain-displacements matrix
                Bmat = zeros(3,nndof);
                Bmat(1,1:2:end) = dx(ig,:);
                Bmat(2,2:2:end) = dy(ig,:);
                Bmat(3,1:2:end) = dy(ig,:);
                Bmat(3,2:2:end) = dx(ig,:);
                
                % smoothed stiffness matrix
                K(edof,edof) = K(edof,edof) + Bmat'*model.Cmat*Bmat*mW(ig);
            end
        end
    end
    model.K = K;
end