function model = quad_patch_force_edge_1(model)
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
        
        for ic = 1:nc
            % support domain numbering
            nsf = length(model.supp{neighbour(ic)});

            % Guass points
            [W,Q] = quadrature(2,'TRIANGULAR',2);
            
            % loop over Gauss points
            N = zeros(6,2);
            for ig = 1:size(W,1)
                [N1,dNdxi] = lagrange_basis('T6',Q(ig,:));
                N(:,1) = N(:,1) + N1*bx;
                N(:,2) = N(:,2) + N1*by;
            end

            if ic == 1
                nodL = model.supp{neighbour(ic)};
                nn = nsf;
                N_T = N;
            else
                i0 = 0;
                for jj = 1:nsf
                    nod = model.supp{neighbour(ic)}(jj);
                    flag = 0;
                    for j = 1:nn
                        if nodL(j) == nod
                            N_T(j,:) = N_T(j,:) + N(jj,:);
                            flag = 1;
                            break; 
                        end   
                    end   
                    if flag == 0 
                        i0 = i0 + 1;
                        nodL(nn+i0) = nod;
                        N_T(nn+i0,:) = N(jj,:);
                    end 
                end   
                nn = nn + i0;
            end
        end
           
        % loop over internal Gauss points
        if nc == 1
            [W,Q] = quadrature(2,'TRIANGULAR',2);
            element = 'T6';
            gcoord = model.Nodes(nodL,:);
            node_sc = [1 2 3 4 5 6];
        else
            [W,Q] = quadrature(2,'GAUSS',2);
            element = 'Q9';
            
            if targetEdge(ivo,5) == 1
                node_sc = [1 7 2 3 9 8 5 6 4];
            elseif targetEdge(ivo,5) == 2
                node_sc = [1 2 7 3 4 8 9 6 5];
            else
                node_sc = [1 2 3 7 4 5 8 9 6];
            end
            gcoord = model.Nodes(nodL(node_sc),:);    
        end
        N_T = N_T(node_sc,:);

        for ig = 1:size(W,1)
            % location of internagl Gauss points
            [~,dNdxi] = lagrange_basis(element,Q(ig,:));
            detJ0 = det(dNdxi'*gcoord);
                
            F(2*nodL-1,:) = F(2*nodL-1,:) + N_T(:,1)*detJ0*W(ig);
            F(2*nodL,:) = F(2*nodL,:) + N_T(:,2)*detJ0*W(ig);
        end

	    clear nodL node_sc
        
    end

    model.F = F;

end