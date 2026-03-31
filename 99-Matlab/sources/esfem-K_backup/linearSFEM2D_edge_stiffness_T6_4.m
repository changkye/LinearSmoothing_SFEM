function model = linearSFEM2D_edge_stiffness_T6_4(model)
	% Edge-based smoothed finite element method
	% 	T6 element
    %
    % Changkye Lee

    % Initialisation
    K = sparse(2*size(model.Nodes,1),2*size(model.Nodes,1));

    % Gauss points on subcell boundary
    ng = 2;
    [Wb,Qb] = quadrature(ng,'GAUSS',1);
            
    % internal Gauss points
    [Wi,Qi] = quadrature(2,'TRIANGULAR',2);

    % target edge and sharing elements
    model = getTargetEdge(model);
    targetEdge = model.targetEdge;
    
    % Loop over edges
    for ivo = 1:size(targetEdge,1)
    	% current target edge
        if targetEdge(ivo,end-2) == 0
            neighbour = targetEdge(ivo,3);
            edges = targetEdge(ivo,5);
        else
            neighbour = targetEdge(ivo,3:end-2);
            edges = targetEdge(ivo,5:6);
        end
        nc = length(neighbour);
        
        if nc == 1
            fx = zeros(3,6); fy = zeros(3,6);
        else
            fx = zeros(4,6); fy = zeros(4,6);
        end
        % loop over sub-cells
        for ic = 1:nc
            
            % support domain numbering
            nsf = length(model.supp{neighbour(ic)});

            % current element connectivity
            wkInd = model.Elements(neighbour(ic),:);

            % current element coordinates
            wkX = model.Nodes(wkInd,:);

            % current boundary of subcell
            bound = [1 2 4; 2 3 5; 3 1 6];
            node_sc = bound(edges(ic),:);

            % subcell boundary length & unit outward normal vectors
            side = cal_side(wkX(1:3,1),wkX(1:3,2));
            [nx,ny] = cal_nx_ny(wkX(1:3,1),wkX(1:3,2),side);

            % initialisation
            if nc == 1
                bxy = zeros(3,size(wkX,1));
            else
                bxy = zeros(4,size(wkX,1));
            end

            % current subcell boundary coordinates
            X = wkX(node_sc,:);

            % loop over Gauss points on current boundary
            for ig = 1:ng
                % map 1D Gauss points to 2D
                [Ng,dNdxi] = lagrange_basis('L3',Qb(ig));
                xieta_gp = Ng'*X;
                detJ = norm(dNdxi'*X);

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
            fx = fx + nx(edges(ic))*bxy;
            fy = fy + ny(edges(ic))*bxy;
        
            % re-numbering
            if ic == 1
                nodL = model.supp{neighbour(ic)};
                nn = nsf;
                Fx = fx; Fy = fy;     
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
        % 
        if nc == 1
            tmp = [1 2 3 4 5 6];
            nodL = nodL(tmp);
            % 
            [W,Q] = quadrature(2,'TRIANGULAR',2);
            element = {'T3'; 'T6'};
            gcoord = model.Nodes(nodL(1:3),:);
                
%             Wmat = zeros(3);
        else 
            if edges(1) == 1
                tmp = [1 7 2 3 9 8 5 6 4];
            elseif edges(1) == 2
                tmp = [1 2 7 3 4 8 9 6 5];
            else
                tmp = [1 2 3 7 4 5 8 9 6];
            end
            nodL = nodL(tmp);
            % 
            [W,Q] = quadrature(2,'GAUSS',2);
            element = {'Q4'; 'Q9'};
            gcoord = model.Nodes(nodL(1:4),:);

%             Wmat = zeros(4);
        end
        Fx = Fx(:,tmp); Fy = Fy(:,tmp);
        
        % shape functions
        Ni = zeros(3,length(nodL));
        for ig = 1:size(W,1)
            % location of internal Gauss points in global space
            [N0,dNdxi] = lagrange_basis(element{1},Q(ig,:));
            detJ(ig) = det(dNdxi'*gcoord);
            mR(ig,:) = N0'*gcoord;
            mW(ig,1) = detJ(ig)*W(ig);
        end

        for ig = 1:size(mW,1)
            % serendipity shape functions
            N = getSerendipityShapeFunc_lagrange(element{2},wkX,mR(ig,:));
            if nc == 1
                Ni(1,:) = Ni(1,:) + (N*mW(ig))';
            else
                Ni = Ni + (N*mW(ig))'.*[1; mR(ig,1); mR(ig,2)];
            end
        end

        % W matrix at interior Gauss points
        Wmat = [];
        for ig = 1:size(mW,1)
            if nc == 1
                Wmat(:,ig) = mW(ig)*[1; mR(ig,1); mR(ig,2)];
            else
                Wmat(:,ig) = mW(ig)*[1; mR(ig,1); mR(ig,2); mR(ig,1)*mR(ig,2)];
            end
        end

        if nc == 1
            Fx(2,:) = Fx(2,:) - Ni(1,:);
            Fy(3,:) = Fy(3,:) - Ni(1,:);
        else
            Fx(2,:) = Fx(2,:) - Ni(1,:);
            Fx(4,:) = Fx(4,:) - Ni(3,:);
            Fy(3,:) = Fy(3,:) - Ni(1,:);
            Fy(4,:) = Fy(4,:) - Ni(2,:);
        end
        % get derivation of basis functions
        dx = Wmat\Fx; dy = Wmat\Fy;

        % smoothed stiffness matrix
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
            % stiffness matrix
            K(edof,edof) = K(edof,edof) + Bmat'*model.Cmat*Bmat*mW(ig);
        end
        clear N nodL
    end
    
    model.K = K;
end