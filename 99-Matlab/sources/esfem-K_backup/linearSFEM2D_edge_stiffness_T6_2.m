function model = linearSFEM2D_edge_stiffness_T6_2(model)
	% Edge-based smoothed finite element method
	% 	T6 element
    %
    % Changkye Lee

    % Initialisation
    K = sparse(2*size(model.Nodes,1),2*size(model.Nodes,1));

    % target edge and sharing elements
    model = getTargetEdge(model);
    targetEdge = model.targetEdge;
    
    % Gauss points on subcell boundary
    ng = 2;
    [Wb,Qb] = quadrature(ng,'GAUSS',1);
    
    % internal Gauss points
	[Wi,Qi] = quadrature(2,'TRIANGULAR',2);
    
    % Loop over edges
    for ivo = 1:size(targetEdge,1)
    	% current target edge
        if targetEdge(ivo,end-2) == 0
            neighbour = targetEdge(ivo,3);
        else
            neighbour = targetEdge(ivo,3:end-2);
        end
        nc = length(neighbour);
        
        subNodes = [];
        % loop over sub-cells
        for ic = 1:nc
            
            % current element connectivity
            wkInd = model.Elements(neighbour(ic),:);

            % current element coordinates
            wkX = model.Nodes(wkInd,:);

            % create subcell
            gcoord = [wkX; mean(wkX(1:3,:))];
            gcoord = [gcoord; mean(gcoord([1,7],:)); 
                mean(gcoord([2,7],:)); mean(gcoord([3,7],:))];

            % support domain numbering
            nsf = length(model.supp{neighbour(ic)});
            
            if targetEdge(ivo,4+ic) == 1
                node_sc = [1,2,7,4,9,8];
            elseif targetEdge(ivo,4+ic) == 2
                node_sc = [2,3,7,5,10,9];
            else
                node_sc = [3,1,7,6,8,10];
            end 
%             if sum(targetEdge(ivo,1:2)) == sum(wkInd(1:2))
%                 node_sc = [1,2,7,4,9,8];
%             elseif sum(targetEdge(ivo,1:2)) == sum(wkInd(2:3))
%                 node_sc = [2,3,7,5,10,9];
%             elseif sum(targetEdge(ivo,1:2)) == sum(wkInd([1,3]))
%                 node_sc = [3,1,7,6,8,10];
%             end 
            X = gcoord(node_sc,:);
            
            if model.flag == 1
                tri = [1,2,7,4,9,8; 2,3,7,5,10,9; 3,1,7,6,8,10];
                trimesh(tri,gcoord(:,1),gcoord(:,2)); 
            end
            
            subNodes = [subNodes; X];
            if model.flag==1; plot(subNodes(:,1),subNodes(:,2),'k^'); end
            
            % shape functions at internal Gauss points
            Ni = zeros(3,size(wkX,1));
            if nc == 1
                Wmat = zeros(3);
            else
                Wmat = zeros(4);
            end
            
            % loop over internal Gauss points
            for ig = 1:size(Wi,1)
                % location of internagl Gauss points
                [N1,dN1dxi] = lagrange_basis('T3',Qi(ig,:));
                detJ0(ig) = det(dN1dxi'*X(1:3,:));
                mR(ig,:) = N1'*X(1:3,:);

                % shape functions & W matrix
                N = getSerendipityShapeFunc_lagrange('T6',wkX,mR(ig,:));
                if nc == 1
                    Ni(1,:) = Ni(1,:) + (N*Wi(ig)*detJ0(ig))';
                elseif nc == 2
                    Ni = Ni + [(N*Wi(ig)*detJ0(ig))'; (N*Wi(ig)*detJ0(ig)*mR(ig,1))';
                        (N*Wi(ig)*detJ0(ig)*mR(ig,2))'];
                end
            end
            
            % subcell boundary length
            side = cal_side(X(1:3,1),X(1:3,2));

            % outward normal vectors
            [nx,ny] = cal_nx_ny(X(1:3,1),X(1:3,2),side);

            % initialise
            if nc == 1
                fx = zeros(3,size(wkX,1)); fy = zeros(3,size(wkX,1));
            else
                fx = zeros(4,size(wkX,1)); fy = zeros(4,size(wkX,1));
            end
            bound = [1 2 4; 2 3 5; 3 1 6];

            % loop over subcell boundary
            for is = 1:size(bound,1)
                if nc == 1
                    bxy = zeros(3,size(wkX,1));
                else
                    bxy = zeros(4,size(wkX,1));
                end
                subX = X(bound(is,:),:);

                % loop over Gauss points on boundary
                for ig = 1:ng
                    % map 1D Gauss point to 2D
                    [Ng,dNdxi] = lagrange_basis('L3',Qb(ig));
                    xieta_gp = Ng'*subX;
                    detJ = norm(dNdxi'*subX);

                    % shape functions on boundary
                    N_T = zeros(size(wkX,1),1);
                    N_T = getSerendipityShapeFunc_lagrange('T6',wkX,xieta_gp);
                    
                    N = N_T'*detJ*Wb(ig);
                    if nc == 1
                        bxy = bxy + [N; N*xieta_gp(1); N*xieta_gp(2)];
                    else
                        bxy = bxy + [N; N*xieta_gp(1); N*xieta_gp(2); N*xieta_gp(1)*xieta_gp(2)];
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
        else
            Fx(2,:) = Fx(2,:) - Nj(1,:);
            Fx(4,:) = Fx(4,:) - Nj(3,:);
            % 
            Fy(3,:) = Fy(3,:) - Nj(1,:);
            Fy(4,:) = Fy(4,:) - Nj(2,:);
        end

        % loop over internal Gauss points
        if nc == 1
            [W,Q] = quadrature(2,'TRIANGULAR',2);
            element = 'T3';
            subNodes = model.Nodes(nodL(1:3),:);
            subNode = subNodes(1:3,:);
            node_sc = [1 2 3 4 5 6];
            Fx = Fx(:,node_sc);
            Fy = Fy(:,node_sc);
            %
            subNode = model.Nodes(nodL(node_sc(1:3)),:);
        else
            [W,Q] = quadrature(2,'GAUSS',2);
            element = 'Q4';
            
            if targetEdge(ivo,5) == 1
                node_sc = [1 7 2 3 9 8 5 6 4];
                
            elseif targetEdge(ivo,5) == 2
                node_sc = [1 2 7 3 4 8 9 6 5];
            else
                node_sc = [1 2 3 7 4 5 8 9 6];
            end
            subNode = model.Nodes(nodL(node_sc(1:4)),:);
            Fx = Fx(:,node_sc);
            Fy = Fy(:,node_sc);
        end
        
        for ig = 1:size(W,1)
            % location of internagl Gauss points
            [N,dNdxi] = lagrange_basis(element,Q(ig,:));
            detJ0(ig) = det(dNdxi'*subNode);
            mR(ig,:) = N'*subNode;

            % W matrix
            if nc == 1
                Wmat(:,ig) = W(ig)*detJ0(ig)*[1; mR(ig,1); mR(ig,2)];
            elseif nc == 2
                Wmat(:,ig) = W(ig)*detJ0(ig)*[1; mR(ig,1); mR(ig,2); mR(ig,1)*mR(ig,2)];
            end
        end

        % get derivatives of basis functions
        dx = Wmat\Fx;
        dy = Wmat\Fy;

        % compute smoothed stiffness matrix
        edof = zeros(1,2*length(nodL));
        edof(1:2:end) = 2*nodL(node_sc) - 1;
        edof(2:2:end) = 2*nodL(node_sc);
        for ig = 1:size(W,1)
            % smoothed strain-displacement matrix
            Bmat = zeros(3,2*length(nodL));
            Bmat(1,1:2:end) = dx(ig,:);
            Bmat(2,2:2:end) = dy(ig,:);
            Bmat(3,1:2:end) = dy(ig,:);
            Bmat(3,2:2:end) = dx(ig,:);

            % smoothed stiffness matrix
            K(edof,edof) = K(edof,edof) + Bmat'*model.Cmat*Bmat*W(ig)*detJ0(ig);
        end
        clear N nodL
    end
    model.K = K;
end