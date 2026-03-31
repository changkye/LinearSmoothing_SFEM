function model = linearFEM2D_stiffness(model)
    % LinearFEM2D : compute 2D stiffness matrix
    % 
    % Changkye Lee, Dept. of Mechanical Engineering,
    % CISTIB, The University of Sheffield,
    % changkye.lee@sheffield.ac.uk, November 2015.
    elemType = model.elemType;
    Nodes = model.Nodes;
    Elements = model.Elements;

    % Stiffness matrix
    K = sparse(2*size(Nodes,1),2*size(Nodes,1));

    % Gauss points
    if strcmp(model.elemType,'T3')
        ng = 1;
        element = 'TRIANGULAR';
    elseif strcmp(model.elemType,'T6')
        ng = 2;
        element = 'TRIANGULAR';
    elseif strcmp(model.elemType,'Q4')
        ng = 1;
        element = 'GAUSS';
    elseif strcmp(model.elemType,'Q8')
        ng = 2;
        element = 'GAUSS';
    end
    [W,Q] = quadrature(ng,element,2);

    % Loop over elements
    for el = 1:size(Elements,1)
        % current element connectivity
        wkInd = Elements(el,:);

        % current element coordinates 
        wkX = Nodes(wkInd,:);
        nndof = length(wkInd);

        % current element dof
        edof(1:2:2*nndof) = 2*wkInd - 1;
        edof(2:2:2*nndof) = 2*wkInd;
        
        % Loop over Gauss points
        for ig = 1:size(W,1)
            
            % Shape functions & their derivs
            [N,dNdxi] = lagrange_basis(model.elemType,Q(ig,:));
            
            % Jacobian
            J0 = wkX'*dNdxi;
            dNdX = dNdxi*inv(J0);
            
            % Strain-displacement matrix
            Bmat = zeros(3,2*nndof);
            Bmat(1,1:2:2*nndof) = dNdX(:,1);
            Bmat(2,2:2:2*nndof) = dNdX(:,2);
            Bmat(3,1:2:2*nndof) = dNdX(:,2);
            Bmat(3,2:2:2*nndof) = dNdX(:,1);
            
            % assemble global stiffness matrix
            K(edof,edof) = K(edof,edof) + Bmat'*model.Cmat*Bmat*W(ig)*det(J0);
        end
    end

    model.K = K;
end