function model = nonlinearFEM2D_stiffness(model,uu)
    % NonlinearFEM2D : compute 2D tangent stiffness matrix
    % 
    % Changkye Lee, Dept. of Mechanical Engineering,
    % CISTIB, The University of Sheffield,
    % changkye.lee@sheffield.ac.uk, November 2015.
    elemType = model.elemType;
    Nodes = model.Nodes;
    Elements = model.Elements;

    % Stiffness matrix
    K = sparse(2*size(Nodes,1),2*size(Nodes,1));
    R = sparse(2*size(Nodes,1),1);
    EE = 0.0;

    % Gauss points
    if strcmp(model.elemType,'T3')
        ng = 1;
        element = 'TRIANGULAR';
    elseif strcmp(model.elemType,'T6')
        ng = 3;
        element = 'TRIANGULAR';
    elseif strcmp(model.elemType,'Q4')
        ng = 1;
        element = 'GAUSS';
    elseif strcmp(model.elemType,'Q8')
        ng = 2;
        element = 'GAUSS';
    end
    [W,Q] = quadrature(ng,element,2);
    [We,Qe] = quadrature(1,element,2);


    % Loop over elements
    for el = 1:size(Elements,1)
        % current element connectivity
        wkInd = Elements(el,:);

        % current element coordinates 
        wkX = Nodes(wkInd,:);
        nndof = length(wkInd);

        % current element displacement
        wkU = zeros(nndof,2);
        wkU(:,1) = uu(2*wkInd-1,1);
        wkU(:,2) = uu(2*wkInd,1);

        % current element dof
        edof = zeros(1,2*nndof);
        edof(1:2:2*nndof) = 2*wkInd - 1;
        edof(2:2:2*nndof) = 2*wkInd;
        
        % current element area
        Ae = polyarea(wkX(:,1),wkX(:,2));

        % Loop over Gauss points
        for ig = 1:size(W,1)
    
            % Shape functions & their derivs
            [~,dNdxi] = lagrange_basis(model.elemType,Q(ig,:));
    
            % Jacobian
            J0 = wkX'*dNdxi;
            dNdX = dNdxi*inv(J0);
    
            % Strain-displacement matrix
            [Bmat,Bgeo,Fmat] = getNonlinearBmat2D(dNdX,wkU,nndof);
            [Cmat,Smat,~] = nonlinearConstitutive(model.param,Fmat);
            
            % assemble global tangent stiffness matrix
            K(edof,edof) = K(edof,edof) + (Bmat'*Cmat*Bmat + ...
                Bgeo'*Smat*Bgeo)*W(ig)*det(J0);
            
            % assemble gobal internal force vector
            R(edof,1) = R(edof,1) + Bmat'*[Smat(1,1); Smat(2,2); Smat(1,2)]*W(ig)*det(J0);

        end
        for jg = 1:size(We,1)
            
            % Shape functions & their derivs
            [~,dNdxi] = lagrange_basis(model.elemType,Qe(jg,:));
    
            % Jacobian
            J0 = wkX'*dNdxi;
            dNdX = dNdxi*inv(J0);
    
            % Strain-displacement matrix
            [~,~,Fmat] = getNonlinearBmat2D(dNdX,wkU,nndof);
            [~,~,W0] = nonlinearConstitutive(model.param,Fmat);
            
            % compute strain energy
            EE = EE + W0*Ae;
        end
    end
    model.K = K;
    model.R = R;
    model.EE = EE;
end
