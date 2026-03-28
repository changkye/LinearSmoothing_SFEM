function model = linearSFEM2D_cell_stiffness_Q8(model)
	% compute linear 2D cell-based smoothing stiffness matrix
	% 	with the linear smoothing function
    %   -> no subcell division
	% 
	% Changkye Lee
	% DSOC National Research Center for Disaster-free & Safe Ocean City,
	% Dong-A University, Korea.
	% changkyelee@gmail.com, December 2018.

	% initialise
    K = sparse(2*size(model.Nodes,1),2*size(model.Nodes,1));

    % Gauss points for internal 
    [Wi,Qi] = quadrature(2,'GAUSS',2);
    % Gauss points for boundary
    ng = 2;
    [Wb,Qb] = quadrature(ng,'GAUSS',1);

    % loop for elements
    for ivo = 1:size(model.Elements,1)
        wkInd = model.Elements(ivo,:);
        wkX = model.Nodes(wkInd,:);
        nndof = length(wkInd);
    
        % get the global index
        edof = zeros(1,2*nndof);
        edof(1:2:2*nndof) = 2*wkInd - 1;
        edof(2:2:2*nndof) = 2*wkInd;
    
        % create subcells
        gcoord = wkX;
        sub_quad = [1 2 3 4 5 6 7 8];

        % compute shape functions at internal Gauss points
        xy = gcoord;
        Ni = zeros(3,size(wkX,1));
        Wmat = zeros(size(Wi,1));
        for ig = 1:size(Wi,1)
            % location of internal Gauss points
            [N1,dN1dxi] = lagrange_basis('Q4',Qi(ig,:));
            detJ0(ig) = det(dN1dxi'*xy(1:4,:));
            mR(ig,:) = N1'*xy(1:4,:);
            % shape functions at internal Gauss points
            N = getSerendipityShapeFunc_lagrange(model.elemType,wkX,mR(ig,:));
            Ni(1,:) = Ni(1,:) + (N*Wi(ig)*detJ0(ig))';
            Ni(2,:) = Ni(2,:) + (N*Wi(ig)*detJ0(ig)*mR(ig,1))';
            Ni(3,:) = Ni(3,:) + (N*Wi(ig)*detJ0(ig)*mR(ig,2))';
            % compute W matrix
            Wmat(:,ig) = Wi(ig)*detJ0(ig)*[1; mR(ig,1); mR(ig,2);
                mR(ig,1)*mR(ig,2)];
        end
        
        % construct smoothed shape functions
        bound = [1 2 5; 2 3 6; 4 7 3; 4 1 8];
           
        % outward normal vectors
        [nx,ny] = cal_nx_ny(xy(1:4,1),xy(1:4,2),cal_side(xy(1:4,1),xy(1:4,2)));
        fx = zeros(4,size(wkX,1)); fy = zeros(4,size(wkX,1));
        for is = 1:size(bound,1)
            bxy = zeros(4,size(wkX,1));
            node_sc = sub_quad(bound(is,:));
            for ig = 1:ng
                % location of Gauss points on the boundary
                [Ng,dNgdxi] = lagrange_basis('L3',Qb(ig));
                xieta_gp = Ng'*xy(bound(is,:),:);
                % map 1D points to the element
                detJ = norm(dNgdxi'*xy(bound(is,:),:));
                % compute shape functions
                N_T = zeros(size(wkX,1),1);
                N_T = getSerendipityShapeFunc_lagrange('Q8',wkX,xieta_gp);
                
                bxy = bxy + [N_T'*detJ*Wb(ig); N_T'*detJ*Wb(ig)*xieta_gp(1); 
                        N_T'*detJ*Wb(ig)*xieta_gp(2); 
                        N_T'*detJ*Wb(ig)*xieta_gp(1)*xieta_gp(2)];
            end
            fx = fx + nx(is)*bxy;
            fy = fy + ny(is)*bxy;
        end
        fx(2,:) = fx(2,:) - Ni(1,:);
        fx(4,:) = fx(4,:) - Ni(3,:);
        %
        fy(3,:) = fy(3,:) - Ni(1,:);
        fy(4,:) = fy(4,:) - Ni(2,:);

        % get derivatives basis functions
        dx = Wmat\fx; dy = Wmat\fy;
                           
        % compute stiffness matrix
        for ig = 1:size(Wi,1)
            % strain-displacement matrix
            Bmat = zeros(3,2*nndof);
            Bmat(1,1:2:2*nndof) = dx(ig,:);
            Bmat(2,2:2:2*nndof) = dy(ig,:);
            Bmat(3,1:2:2*nndof) = dy(ig,:);
            Bmat(3,2:2:2*nndof) = dx(ig,:);
            % smoothed stiffness matrix
            K(edof,edof) = K(edof,edof) + Bmat'*model.Cmat*Bmat*Wi(ig)*detJ0(ig);
        end
    end
    model.K = K;
end