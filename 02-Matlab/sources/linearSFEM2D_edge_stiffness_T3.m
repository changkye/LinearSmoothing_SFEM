function model = linearSFEM2D_edge_stiffness_T3(model)
	% Edge-based smoothed finite element method
	% 	linear T3 element
    %
    % Changkye Lee, School of Engineering,
    % iMAM, Cardiff University,
    % LeeC15@cardiff.ac.uk

    % Initialisation
    K = sparse(2*size(model.Nodes,1),2*size(model.Nodes,1));

    % Gauss points
    if strcmp(model.elemType,'T3') | strcmp(model.elemType,'Q4')
        ng = 1;
        element = 'T3';
    elseif strcmp(model.elemType,'T4')
        ng = 2;
        element = 'T4';
    end
    [W,Q] = quadrature(ng,'GAUSS',1);

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
            % current element connectivity
            wkInd = model.Elements(neighbour(ic),:);

            % current element coordinates
            wkX = model.Nodes(wkInd,:);

            % sub-cell coordinates
            if strcmp(model.elemType,'T4')
                gcoord = wkX;
            else
                gcoord = [wkX; mean(wkX)];
            end

            % sub-cell connectivity
            if sum(targetEdge(ivo,1:2)) == sum(wkInd(1:2))
                X = gcoord([1,2,4],:);
            elseif sum(targetEdge(ivo,1:2)) == sum(wkInd(2:3))
                X = gcoord([2,3,4],:);
            else
                X = gcoord([3,1,4],:);
            end

            % support node numbering
            nsf = length(model.supp{neighbour(ic)});

            % boundary length of sub-cell
            side = cal_side(X(:,1),X(:,2));

            % outward normal vectors
            [nx,ny] = cal_nx_ny(X(:,1),X(:,2),side);

            % loop over boundary
            bx = zeros(nsf,1); by = zeros(nsf,1);
            for is = 1:length(nx)
                % loop over Gauss point
                for ig = 1:ng
                    % 1D line element shape functions
                    N_g = lagrange_basis('L2',Q(ig));
                    [xi1,eta1] = xi_eta4xy('T3',is,N_g);
                    % map 1D Gauss point to 2D triangle
                    N_xy = lagrange_basis('T3',[xi1 eta1]);
                    xy_g = N_xy'*[X(:,1) X(:,2)];
                    xieta = cal_xieta('T3',xy_g,wkX);
                    % 2D shape function s
                    N_T = lagrange_basis(model.elemType,xieta);
                    % Jacobian
                    J = side(is)/2;
                    % smoothed shape funcs.
                    bx = bx + nx(is)*N_T*J*W(ig);
                    by = by + ny(is)*N_T*J*W(ig);
                end   
            end
            bx = bx/model.subA(ivo); by = by/model.subA(ivo);

            % assemble shape functions
            if ic == 1
                nodL = model.supp{neighbour(ic)};
                nn = nsf;
                bxy = [bx by];
            else  
                i0 = 0;
                for jj = 1:nsf
                    nod = model.supp{neighbour(ic)}(jj);
                    flag = 0;
                    for j = 1:nn
                        if nodL(j) == nod
                            bxy(j,:) = bxy(j,:) + [bx(jj) by(jj)];
                            flag = 1;
                            break; 
                        end  
                    end  
                    if flag == 0 
                        i0 = i0 + 1;
                        nodL(nn+i0) = nod;
                        bxy(nn+i0,:) = [bx(jj) by(jj)];
                    end
                end  
                nn = nn + i0;
            end  
        end  
        
        % smoothing domain connectivity
        numnod = length(nodL);
        edof = [];
        edof(1:2:2*numnod) = 2*nodL - 1;
        edof(2:2:2*numnod) = 2*nodL;

        % strain-displacement matrix
        Bmat = zeros(3,2*numnod);
        Bmat(1,1:2:2*numnod) = bxy(:,1);
        Bmat(2,2:2:2*numnod) = bxy(:,2);
        Bmat(3,1:2:2*numnod) = bxy(:,2);
        Bmat(3,2:2:2*numnod) = bxy(:,1);
        
        % assemble global smoothed stiffness matrix
        K(edof,edof) = K(edof,edof) + Bmat'*model.Cmat*Bmat*model.subA(ivo);
        
        clear nodL
    end
    
    model.K = K;
   
end