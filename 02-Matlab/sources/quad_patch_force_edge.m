function model = quad_patch_force_edge(model)
	% 
	F = zeros(2*size(model.Nodes,1),1);

    targetEdge = model.targetEdge;
        
    % exact body force
    bx = -0.2*model.Cmat(1,1) - 0.15*model.Cmat(1,2) - 0.55*model.Cmat(3,3);
    by = -0.1*model.Cmat(1,2) - 0.2*model.Cmat(2,2) - 0.2*model.Cmat(3,3);
    
    % loop over target edges
    for ivo = 1:size(targetEdge,1)
        % current target edge
        if targetEdge(ivo,end-2) == 0
            neighbour = targetEdge(ivo,3);
        else
            neighbour = targetEdge(ivo,3:end-2);
        end
        nc = length(neighbour);

        % loop over sub-cells
        if nc == 1
            % Guass points
            [W,Q] = quadrature(2,'TRIANGULAR',2);
            
            % current subcell coordinates
            wkX = model.Nodes(model.Elements(neighbour,:),:);
            
            % loop over Gauss points
            for ig = 1:size(W,1)
                [N,dNdxi] = lagrange_basis('T6',Q(ig,:));
                detJ = det(dNdxi'*wkX);
                F(2*model.Elements(neighbour,:)-1) = F(2*model.Elements(neighbour,:)-1) + ...
                    N*bx*detJ*W(ig);
                F(2*model.Elements(neighbour,:)) = F(2*model.Elements(neighbour,:)) + ...
                    N*by*detJ*W(ig);
            end
        else
            % Guass points
            [W,Q] = quadrature(2,'GAUSS',2);

            % current subcell connectivity
            for ic = 1:nc
                nsf = length(model.supp{neighbour(ic)});
                if ic == 1
                    nodL = model.supp{neighbour(ic)};
                    nn = nsf;
                else
                    i0 = 0;
                    for jj = 1:nsf
                        nod = model.supp{neighbour(ic)}(jj);
                        flag = 0;
                        for j = 1:nn
                            if nodL(j) == nod
                                flag = 1;
                                break; 
                            end  
                        end  
                        if flag == 0 
                            i0 = i0 + 1;
                            nodL(nn+i0) = nod;
                        end
                    end  
                    nn = nn + i0;
                end  
            end  
            if targetEdge(ivo,5) == 1
                node_sc = [1 7 2 3 9 8 5 6 4];
            elseif targetEdge(ivo,5) == 2
                node_sc = [1 2 7 3 4 8 9 6 5];
            else
                node_sc = [1 2 3 7 4 5 8 9 6];
            end
            wkInd = nodL(node_sc); 

            % current subcell coordinates
            wkX = model.Nodes(wkInd,:);

            for ig = 1:size(W,1)
                [N,dNdxi] = lagrange_basis('Q9',Q(ig,:));
                detJ = det(dNdxi'*wkX);
                F(2*wkInd-1) = F(2*wkInd-1) + N*bx*detJ*W(ig);
                F(2*wkInd) = F(2*wkInd) + N*by*detJ*W(ig);
            end

            clear wkInd wkX N dNdxi
        end
    end

    model.F = F;

end