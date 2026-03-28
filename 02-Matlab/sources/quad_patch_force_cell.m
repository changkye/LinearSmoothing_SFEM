function model = quad_patch_force_cell(model)
	% 
	F = zeros(2*size(model.Nodes,1),1);
    
    % Guass points
    if strcmp(model.elemType,'T3') | strcmp(model.elemType,'T6')
        element = 'TRIANGULAR';
    else
        element = 'GAUSS';
    end
    [W,Q] = quadrature(2,element,2);
    
    % exact body force
    bx = -0.2*model.Cmat(1,1) - 0.15*model.Cmat(1,2) - 0.55*model.Cmat(3,3);
    by = -0.1*model.Cmat(1,2) - 0.2*model.Cmat(2,2) - 0.2*model.Cmat(3,3);
    
    % loop over elements
    for iel = 1:size(model.Elements,1)
        
        % current element connectivity
        wkInd = model.Elements(iel,:);
        
        % current element coordinates
        wkX = model.Nodes(wkInd,:);
        
        % loop over Gauss points
        for ig = 1:size(W,1)
            [N,dNdxi] = lagrange_basis(model.elemType,Q(ig,:));
%             J0 = wkX'*dNdxi;
            detJ = det(dNdxi'*wkX);
            F(2*wkInd-1,1) = F(2*wkInd-1,1) + N*bx*detJ*W(ig);
            F(2*wkInd,1)   = F(2*wkInd,1)   + N*by*detJ*W(ig);
        end
    end

    model.F = F;
end