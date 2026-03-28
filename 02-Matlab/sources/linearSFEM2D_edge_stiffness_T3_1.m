function model = linearSFEM2D_edge_stiffness_T3_1(model)
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
    
    % interior Gauss points
	[Wi,Qi] = quadrature(2,'TRIANGULAR',2);
    
    % Loop over edges
    for ivo = 1:size(targetEdge,1)
    	% current target edge
        if (targetEdge(ivo,end-2) == 0); neighbour = targetEdge(ivo,3);
        else; neighbour = targetEdge(ivo,3:end-2); end
        nc = length(neighbour);
        
        % loop over sub-cells
        for ic = 1:nc
            
            % support domain numbering
            nsf = length(model.supp{neighbour(ic)});
            
            % current element connectivity
            wkInd = model.Elements(neighbour(ic),:);

            % current element coordinates
            wkX = model.Nodes(wkInd,:);
            
            % loop over internal Gauss points
            detJ0 = zeros(size(Wi,1),1);
            mR0 = zeros(size(Wi,1),2);
            Ni = zeros(3,size(wkX,1));
            for ig = 1:size(Wi,1)
                % location of internagl Gauss points
                [N1,dN1dxi] = lagrange_basis('T3',Qi(ig,:));
                detJ0(ig) = det(dN1dxi'*wkX);
                mR0(ig,:) = N1'*wkX;

                % shape functions & W matrix
                xieta = cal_xieta('T3',mR0(ig,:),wkX);
                N = lagrange_basis('T3',xieta);
                if (nc == 1); Ni(1,:) = Ni(1,:) + (N*Wi(ig)*detJ0(ig))';
                else; Ni = Ni + [(N*Wi(ig)*detJ0(ig))'; (N*Wi(ig)*detJ0(ig)*mR0(ig,1))';
                        (N*Wi(ig)*detJ0(ig)*mR0(ig,2))']; 
                end
            end
            
            % subcell boundary length
            side = cal_side(wkX(:,1),wkX(:,2));

            % outward normal vectors
            [nx,ny] = cal_nx_ny(wkX(:,1),wkX(:,2),side);

            % initialise
            if nc == 1; fx = zeros(3,size(wkX,1)); fy = zeros(3,size(wkX,1));
            else; fx = zeros(4,size(wkX,1)); fy = zeros(4,size(wkX,1)); end
            bound = [1 2; 2 3; 3 1];

            % loop over subcell boundary
            for is = 1:size(bound,1)
                if nc == 1; bxy = zeros(3,size(wkX,1));
                else; bxy = zeros(4,size(wkX,1)); end
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
                    if nc == 1; bxy = bxy + [N; N*xieta_gp(1); N*xieta_gp(2)];
                    else; bxy = bxy + [N; N*xieta_gp(1); N*xieta_gp(2); N*xieta_gp(1)*xieta_gp(2)];
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
            
            % re-numbering
            node_sc = [1 2 3];
            nodL = nodL(node_sc);
            % 
            [W,Q] = quadrature(2,'TRIANGULAR',2);
            element = 'T3';
            gcoord = model.Nodes(nodL,:);
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
            if (targetEdge(ivo,5) == 1); node_sc = [1 4 2 3];
            elseif (targetEdge(ivo,5) == 2); node_sc = [1 2 4 3];
            else; node_sc = [1 2 3 4]; end
            nodL = nodL(node_sc);
            %
            [W,Q] = quadrature(2,'GAUSS',2);
            element = 'Q4';
            gcoord = model.Nodes(nodL,:);
            %
            Fx = Fx(:,node_sc);
            Fy = Fy(:,node_sc);
        end
        
        % W matrix at interior Gauss points
        detJ1 = zeros(size(W,1),1);
        mR1 = zeros(size(W,1),2);
        if nc == 1; Wmat = zeros(3); else; Wmat = zeros(4); end 
        for ig = 1:size(W,1)
            % location of internagl Gauss points
            [N,dNdxi] = lagrange_basis(element,Q(ig,:));
            detJ1(ig) = det(dNdxi'*gcoord);
            mR1(ig,:) = N'*gcoord;

            % W matrix
            if (nc == 1); Wmat(:,ig) = W(ig)*detJ1(ig)*[1; mR1(ig,1); mR1(ig,2)];
            else; Wmat(:,ig) = W(ig)*detJ1(ig)*[1; mR1(ig,1); mR1(ig,2); mR1(ig,1)*mR1(ig,2)]; end
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
            K(edof,edof) = K(edof,edof) + Bmat'*model.Cmat*Bmat*W(ig)*detJ1(ig);
        end
        clear N nodL
    end
    model.K = K;
end