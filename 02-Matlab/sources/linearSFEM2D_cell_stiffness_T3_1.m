function model = linearSFEM2D_cell_stiffness_T3_1(model)
	% compute linear 2D cell-based smoothing stiffness matrix
	% 	with the linear smoothing function
    %   -> no subcell division
	% 
	% Changkye Lee
	% DSOC National Research Center for Disaster-free & Safe Ocean City,
	% Dong-A University, Korea.
	% changkyelee@gmail.com, Novemver 2018.

	% initialise
    K = sparse(2*size(model.Nodes,1),2*size(model.Nodes,1));

    % Gauss points for internal 
    [Wi,Qi] = quadrature(2,'TRIANGULAR',2);
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
        sub_tri = [1 2 3];

        % compute shape functions at internal Gauss points
        xy = gcoord;
        Ni = zeros(1,size(wkX,1));
        Wmat = zeros(size(Wi,1));
        for ig = 1:size(Wi,1)
            % location of internal Gauss points
            [N1,dN1dxi] = lagrange_basis('T3',Qi(ig,:));
            detJ0(ig) = det(dN1dxi'*xy);
            mR(ig,:) = N1'*xy;
            % shape functions at internal Gauss points
            xieta = cal_xieta('T3',mR(ig,:),wkX);
            N = lagrange_basis('T3',xieta);
            
            Ni = Ni + (N*Wi(ig)*detJ0(ig))';
            % compute W matrix
            Wmat(:,ig) = Wi(ig)*detJ0(ig)*[1; mR(ig,1); mR(ig,2)];
        end
        
        % construct smoothed shape functions
        bound = [1 2; 2 3; 3 1];
           
        % outward normal vectors
        [nx,ny] = cal_nx_ny(xy(:,1),xy(:,2),cal_side(xy(:,1),xy(:,2)));
        fx = zeros(3,size(wkX,1)); fy = zeros(3,size(wkX,1));
        for is = 1:size(bound,1)
            bxy = zeros(3,size(wkX,1));
            node_sc = sub_tri(bound(is,:));
            for ig = 1:ng
                % location of Gauss points on the boundary
                [Ng,dNgdxi] = lagrange_basis('L2',Qb(ig));
                xieta_gp = Ng'*xy(bound(is,:),:);
                % map 1D points to the element
                detJ = norm(dNgdxi'*xy(bound(is,:),:));
                % compute shape functions
                N_T = zeros(size(wkX,1),1);
                xieta = cal_xieta('T3',xieta_gp,wkX);
                N_T = lagrange_basis('T3',xieta);
                
                bxy = bxy + [N_T'*detJ*Wb(ig); N_T'*detJ*Wb(ig)*xieta_gp(1); 
                        N_T'*detJ*Wb(ig)*xieta_gp(2)];
            end
            fx = fx + nx(is)*bxy;
            fy = fy + ny(is)*bxy;
        end
        fx(2,:) = fx(2,:) - Ni;
        fy(3,:) = fy(3,:) - Ni;

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