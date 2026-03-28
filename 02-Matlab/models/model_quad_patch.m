function model = model_quad_patch(elemType,params,num_elem)


    
	% models for quadratic patch
	L = [0 1; 0 1];
	P = 0.0;
	numEls = num_elem*ones(1,2);

	% -------------------------------------------------------------------------
	%	mesh generation
	% -------------------------------------------------------------------------
	mesh = str2func(['element' elemType]);
	[Nodes,Elements] = mesh(L,numEls); 

	% -------------------------------------------------------------------------
	% 	boundary conditions
	% -------------------------------------------------------------------------
	% get boundary nodes
	btmNodes = find(Nodes(:,2) == min(Nodes(:,2)));
    topNodes = find(Nodes(:,2) == max(Nodes(:,2)));
    rightNodes = find(Nodes(:,1) == max(Nodes(:,1)));
    leftNodes = find(Nodes(:,1) == min(Nodes(:,1)));
    fixNode = unique([btmNodes; topNodes; rightNodes; leftNodes]);

    % fixed boundary nodes' coordinates
    fixXY = Nodes(fixNode,:);
    % 
    bcdof = zeros(1,2*length(fixNode));
    bcdof(1:2:2*length(fixNode)) = 2*fixNode - 1;
    bcdof(2:2:2*length(fixNode)) = 2*fixNode;
    % 
    bcval = zeros(1,2*length(fixNode));
	bcval(1:2:2*length(fixNode)) = 0.1*fixXY(:,1).^2 + ...
        0.1*fixXY(:,1).*fixXY(:,2) + 0.2*fixXY(:,2).^2;
    bcval(2:2:2*length(fixNode)) = 0.05*fixXY(:,1).^2 + ...
        0.15*fixXY(:,1).*fixXY(:,2) + 0.1*fixXY(:,2).^2;	;

    % -------------------------------------------------------------------------
	%	external force
	% -------------------------------------------------------------------------
	F = zeros(2*size(Nodes,1),1);
    
    stressType = 'plane_stn';
	Cmat = zeros(3);
    E0 = params(1); nu = params(2);
	if strcmp(stressType,'plane_stn')
        C0 = E0/(1+nu)/(1-2*nu);
        Cmat(1,1) = 1 - nu;
        Cmat(1,2) = nu;
        Cmat(2,1) = Cmat(1,2);
        Cmat(2,2) = Cmat(1,1);
        Cmat(3,3) = (1-2*nu)/2.;
        Cmat = C0*Cmat;
	elseif strcmp(stressType,'plane_str')
        C0 = E0/(1-nu*nu);
        Cmat(1,1) = 1.;
        Cmat(1,2) = nu;
        Cmat(2,1) = nu;
        Cmat(2,2) = 1.;
	    Cmat(3,3) = (1-nu)/2.;
        Cmat = C0*Cmat;
    end
    
    
	% Guass points
    if strcmp(elemType,'T3') | strcmp(elemType,'T6')
        element = 'TRIANGULAR';
    else
        element = 'GAUSS';
    end
    [W,Q] = quadrature(2,element,2);
    
    % exact body force
    bx = -0.2*Cmat(1,1) - 0.15*Cmat(1,2) - 0.55*Cmat(3,3);
    by = -0.1*Cmat(1,2) - 0.2*Cmat(2,2) - 0.2*Cmat(3,3);
    
    % loop over elements
    for iel = 1:size(Elements,1)
        
        % current element connectivity
        wkInd = Elements(iel,:);
        
        % current element coordinates
        wkX = Nodes(wkInd,:);
        
        % loop over Gauss points
        for ig = 1:size(W,1)
            [N,dNdxi] = lagrange_basis(elemType,Q(ig,:));
%             J0 = wkX'*dNdxi;
            detJ = det(dNdxi'*wkX);
            F(2*wkInd-1,1) = F(2*wkInd-1,1) + N*bx*detJ*W(ig);
            F(2*wkInd,1)   = F(2*wkInd,1)   + N*by*detJ*W(ig);
        end
    end

	% -------------------------------------------------------------------------
	model.elemType = elemType;
	model.L = L;
	model.param = params;
	model.numEls = numEls;
	model.Nodes = Nodes;
	model.Elements = Elements;
	model.bcdof = bcdof;
	model.bcval = bcval;
	model.F = F;

end