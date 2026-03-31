function model = linearSFEM2D_cell_stiffness(model)
	% Cell-based smoothed finite element method
	% 	linear T3/Q4 elements
    %
    % Changkye Lee, School of Engineering,
    % iMAM, Cardiff University,
    % LeeC15@cardiff.ac.uk

    % Initialisation
    K = sparse(2*size(model.Nodes,1),2*size(model.Nodes,1));

    % Gauss points
    if strcmp(model.elemType,'Q4')
        numSD = 4;
        numSub = 1;
        element = 'Q4';
    else
        numSD = 3;
        numSub = 1;
        element = 'T3';
    end
    ng = 1;
    [W,Q] = quadrature(ng,'GAUSS',1);

    % Loop over elements
    for ivo = 1:size(model.Elements,1)
    	% current element connectivity
        wkInd = model.Elements(ivo,:);

        % current element coordinates 
        wkX = model.Nodes(wkInd,:);
        
        % current element dof
        nndof = length(wkInd);
        edof(1:2:2*nndof) = 2*wkInd - 1;
        edof(2:2:2*nndof) = 2*wkInd;
        
        % loop over subcell
        for ic = 1:numSub
            
            % subdomain coordinates
            subXY = subCellXY(model,ic,numSub,wkInd);
            
            if strcmp(model.elemType,'Q4')
                if numSD == 1
                    % 4-------------3
                    % |             |
                    % |     sc1     |
                    % |             |
                    % 1-------------2
                    gcoord = subXY;
                    node_sc = [1 2 3 4];
                elseif numSD == 2
                    % 4------6------3
                    % |      |      |
                    % | sc1  |  xc2 |
                    % |      |      |
                    % 1------5------2
                    xc1 = mean(subXY(1:2,:));
                    xc2 = mean(subXY(3:4,:));
                    node_sc = [1 5 6 4; 5 2 3 6];
                elseif numSD == 3
                    % 4------8------3
                    % |      |  sc3 |
                    % | sc1  6------7
                    % |      |  sc1 |
                    % 1------5------2
                    xc = mean(subXY);
                    xc1 = mean(subXY(1:2,:));
                    xc2 = mean(subXY(2:3,:));
                    xc3 = mean(subXY(3:4,:));
                    gcoord = [subXY; xc1; xc; xc2; xc3];
                    node_sc = [1 5 8 4; 5 2 7 6; 6 7 3 8];
                elseif numSD == 4
                    % 4------7------3
                    % | sc4  |  sc3 |
                    % 8------9------6
                    % | sc1  |  sc2 |
                    % 1------5------2
                    xc = mean(subXY);
                    xc1 = mean(subXY(1:2,:));
                    xc2 = mean(subXY(2:3,:));
                    xc3 = mean(subXY(3:4,:));
                    xc4 = mean(subXY([1,4],:));
                    gcoord = [subXY; xc1; xc2; xc3; xc4; xc];
                    node_sc = [1 5 9 8; 5 2 6 9;
                       9 6 3 7; 8 9 7 4];
                elseif numSD == 8
                    % 04----11----10----09----03
                    % | sc5 | sc6 | sc7 | sc8 |
                    % 12----13----14----15----08
                    % | sc1 | sc2 | sc3 | sc4 |
                    % 01----05----06----07----02
                    xc = mean(subXY);
                    xc1 = (3*subXY(1,:)+subXY(2,:))/4.0;
                    xc2 = mean(subXY(1:2,:));
                    xc3 = (subXY(1,:)+3*subXY(2,:))/4.0;
                    xc4 = mean(subXY(2:3,:));
                    xc5 = (subXY(4,:)+3*subXY(3,:))/4.0;
                    xc6 = mean(subXY(3:4,:));
                    xc7 = (3*subXY(4,:)+subXY(3,:))/4.0;
                    xc8 = mean(subXY([1,4],:));
                    xc9 = (3*subXY(1,:)+subXY(2,:)+subXY(3,:)+3*subXY(4,:))/8.0;
                    xc10 = (subXY(1,:)+3*subXY(2,:)+3*subXY(3,:)+subXY(4,:))/8.0;
                    gcoord = [subXY; xc1; xc2; xc3; xc4; xc5; xc6; xc7; xc8; xc9; xc; xc10];
                    node_sc = [1 5 13 12; 5 6 14 13; 6 7 15 14; 7 2 8 15; 
                        12 13 11 4; 13 14 10 11; 14 15 9 10; 15 8 3 9];
                end
            else
                if numSD == 1
                    % 3
                    % | -
                    % |    -
                    % |       -
                    % |   sc1    -
                    % |             -
                    % 1----------------2
                    gcoord = subXY;
                    node_sc = [1 2 3];
                elseif numSD == 2
                    % 3
                    % | -
                    % |    -
                    % |  sc2  4
                    % |     -   -
                    % |  -   sc1    -
                    % 1----------------2
                    xc = mean(subXY(2:3,:));
                    gcoord = [subXY; xc];
                    node_sc = [1 2 4; 1 4 3];
                elseif numSD == 3
                    % 3
                    % | -
                    % | .  -
                    % |  .    -
                    % |   .4    -
                    % |  .         -
                    % 1----------------2
                    xc = mean(subXY);
                    gcoord = [subXY; xc];
                    node_sc = [1 2 4; 2 3 4; 3 1 4];
                end
            end
            
            % loop over smoothing domain
            for i = 1:numSD

                bx = zeros(nndof,1); by = zeros(nndof,1);
                
                % current smoothing domain coordinates
                xy = gcoord(node_sc(i,:),:);

                % current smoothing domain area
                subA = polyarea(xy(:,1),xy(:,2));

                % length of each smoothing domain boundary
                side = cal_side(xy(:,1),xy(:,2));

                % outward normal vectors on each boundary
                [nx,ny] = cal_nx_ny(xy(:,1),xy(:,2),side);
                ns = length(xy(:,1));
                
                % loop over each smoothing domain boundary
                for is = 1:ns
                    
                    % loop over Gauss point on the boundary
                    for ig = 1:ng
                    	% shape functions
                        N_g = lagrange_basis('L2',Q(ig));

                        % map 1D Gauss points to 2D T3/Q4 elements
                        [xi1,eta1] = xi_eta4xy(element,is,N_g);
                        N_xy = lagrange_basis(element,[xi1 eta1]);
                        xy_g = N_xy'*[xy(:,1) xy(:,2)];
                        xieta = cal_xieta(element,xy_g,wkX);

                        N_T = lagrange_basis(model.elemType,xieta);
                        
                        % Jacobian
                        J = side(is)/2;
                        
                        % Smoothed Shape funcs.
                        bx = bx + nx(is)*N_T*J*W(ig)/subA;
                        by = by + ny(is)*N_T*J*W(ig)/subA;
                    end 
                end
            
                % strain-displacement matrix
               	Bmat = zeros(3,2*length(wkInd));
                Bmat(1,1:2:end) = bx;
                Bmat(2,2:2:end) = by;
                Bmat(3,1:2:end) = by;
                Bmat(3,2:2:end) = bx;
       
                % global stiffness matrix
                K(edof,edof) = K(edof,edof) + Bmat'*model.Cmat*Bmat*subA;
            end
  
        end
         
    end

    model.K = K;

end

%% 
function subXY = subCellXY(model,ic,numSub,wkInd)
    % 
    Nodes = model.Nodes;

    if strcmp(model.elemType,'Q4')
        if numSub == 1
            subXY = Nodes(wkInd,:);
        elseif numSub == 2
            xc1 = mean(Nodes(wkInd(1:2),:));
            xc2 = mean(Nodes(wkInd(3:4),:));

            if ic == 1
                subXY = [Nodes(wkInd(1),:); xc1; 
                    xc2; Nodes(wkInde(4),:)];
            else
                subXY = [xc1; Nodes(wkInd(2),:);
                    Nodes(wkInd(3),:); xc2];
            end
        elseif numSub == 4
            xc = mean(Nodes(wkInd,:));
            xc1 = mean(Nodes(wkInd(1:2),:));
            xc2 = mean(Nodes(wkInd(2:3),:));
            xc3 = mean(Nodes(wkInd(3:4),:));
            xc4 = mean(Nodes(wkInd([1,4],:)));
            if ic == 1
                subXY = [Nodes(wkInd(1),:); xc1; xc; xc4];
            elseif ic == 2
                subXY = [xc1; Nodes(wkInd(2),:); xc2  xc];
            elseif ic == 3
                subXY = [xc; xc2; Nodes(wkInd(3),:); xc3];
            else
                subXY = [xc4; xc; xc3; Nodes(wkInd(4),:)];
            end
        end
    else
        Nodes = model.Nodes;
        if numSub == 1
            subXY = Nodes(wkInd,:);
        elseif numSub == 3
            xc = mean(Nodes(wkInd,:));
            if ic == 1
                subXY = [Nodes(wkInd(1),:); Nodes(wkInd(2),:); xc];
            elseif ic == 2
                subXY = [Nodes(wkInd(2),:); Nodes(wkInd(3),:); xc];
            else
                subXY = [Nodes(wkInd(3),:); Nodes(wkInd(1),:); xc];
            end
        end
    end

end