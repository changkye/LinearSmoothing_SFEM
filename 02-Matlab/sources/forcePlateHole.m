function model = forcePlateHole(model,rtNodes,topNodes)
	%
	F = zeros(2*size(model.Nodes,1),1);
    if strcmp(model.elemType,'T3')
        NPE = 2; % node per edge
        ng = 3;
        element = 'L2';
    else
        NPE = 3; 
        ng = 4;
        element = 'L3';
    end
    
	[t1,t2] = sort(model.Nodes(rtNodes,2));
	rtNodes = rtNodes(t2);
	[t1,t2] = sort(model.Nodes(topNodes,1));
	topNodes = topNodes(t2);

	% right edge
    nume = (length(rtNodes)-1)/(NPE-1);
	
    % loop over the number of elements
    cnt = 0;
    for in = 1:nume
        
        % loop over the number of nodes
        for jn = 1:NPE
            cnt = cnt + 1;
            rtEdges(in,jn) = rtNodes(cnt);
        end
        cnt = cnt - 1;
    end
    if strcmp(model.elemType,'T6'); rtEdges(:,[2,3]) = rtEdges(:,[3,2]); end

    % top edge
    nume = (length(topNodes)-1)/(NPE-1);
	
    % loop over the number of elements
    cnt = 0;
    for in = 1:nume
        
        % loop over the number of nodes
        for jn = 1:NPE
            cnt = cnt + 1;
            topEdges(in,jn) = topNodes(cnt);
        end
        cnt = cnt - 1;
    end
    if strcmp(model.elemType,'T6'); topEdges(:,[2,3]) = topEdges(:,[3,2]); end

    % compute force vector
    [W,Q] = quadrature(ng,'GAUSS',1);
    for iel = 1:size(rtEdges,1)
    	% current edge connectivities
    	wkInd = rtEdges(iel,:);
    	
    	% current edge coordinates
    	wkX = model.Nodes(wkInd,:);

    	% global dof
    	nn = length(wkInd);
    	edof = zeros(2,nn);
    	edof(1,:) = 2*wkInd - 1;
    	edof(2,:) = 2*wkInd;

    	% loop over Gauss points
    	for ig = 1:ng
    		[N,dNdxi] = lagrange_basis(element,Q(ig,:));

    		detJ = det(dNdxi'*wkX(:,2));
    		xieta_gp = N'*wkX;

    		plot(xieta_gp(1),xieta_gp(2),'k^');

    		[S,U] = exactPlateHole(model,xieta_gp,1);

    		F(edof(1,:),1) = F(edof(1,:),1) + N*S(1)*W(ig)*detJ;
    		F(edof(2,:),1) = F(edof(2,:),1) + N*S(3)*W(ig)*detJ;

    	end
    end

    for iel = 1:size(topEdges,1)
    	% current edge connectivities
    	wkInd = topEdges(iel,:);
    	
    	% current edge coordinates
    	wkX = model.Nodes(wkInd,:);

    	% global dof
    	nn = length(wkInd);
    	edof = zeros(2,nn);
    	edof(1,:) = 2*wkInd - 1;
    	edof(2,:) = 2*wkInd;

    	% loop over Gauss points
    	for ig = 1:ng
    		[N,dNdxi] = lagrange_basis(element,Q(ig,:));

    		detJ = det(dNdxi'*wkX(:,2));
    		xieta_gp = N'*wkX;

    		plot(xieta_gp(1),xieta_gp(2),'k^');

    		[S,U] = exactPlateHole(model,xieta_gp,1);

    		F(edof(1,:),1) = F(edof(1,:),1) + N*S(3)*W(ig)*detJ;
    		F(edof(2,:),1) = F(edof(2,:),1) + N*S(2)*W(ig)*detJ;

    	end
    end
    
    model.F = F;
end