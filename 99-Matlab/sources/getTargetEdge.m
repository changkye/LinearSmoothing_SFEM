function model = getTargetEdge(model)
	% 
	% target edge & sharing element index
    for i = 1:size(model.Elements,1)
        supp{i} = model.Elements(i,:);
    end
    subCell = edgeBased2D(model);
    
    if strcmp(model.elemType,'T3')
        % element area
        for i = 1:size(model.Elements,1)
            wkInd = model.Elements(i,:);
            wkX = model.Nodes(wkInd,:);
            triA(i) = abs(polyarea(wkX(:,1),wkX(:,2)));
        end
    
        % subcell area
        subA = zeros(size(subCell,1),1);
        NPE = size(model.Elements,2);
        for i = 1:size(subCell,1)
            edges = find(subCell(i,3:end-2));
            if length(edges) == 1
                subA(i) = triA(subCell(i,3))/NPE;
            elseif length(edges) == 2
                subA(i) = triA(subCell(i,3))/NPE + triA(subCell(i,4))/NPE;
            end
        end
        
    elseif strcmp(model.elemType,'T6')
        % element area
        for i = 1:size(model.Elements,1)
            wkInd = model.Elements(i,:);
            wkX =  model.Nodes(wkInd,:);
            triA(i) = abs(polyarea(wkX(1:3,1),wkX(1:3,2)));
        end

        % subcell area
        subA = zeros(size(subCell,1),1);
        NPE = size(model.Elements(:,1:3),2); % number of edge-based SD per triangle (=3)
        for i = 1:size(subCell,1)
            edges = find(subCell(i,3:end-2));
            if length(edges) == 1
                subA(i) = triA(subCell(i,3))/NPE;
            else
                subA(i) = triA(subCell(i,3))/NPE + triA(subCell(i,4))/NPE;
            end
        end
    end

    model.triA = triA;
    model.subA = subA;
    model.supp = supp;
    model.targetEdge = subCell;
    

end

%%
function subCell = edgeBased2D(model)
    % Define_Edges : Find target edges and elements associated with target
    % edges
    %
    % Changkye Lee, Dept. of Mechancial Engineering,
    % CISTIB, The University of Sheffield,
    % changkye.lee@sheffield.ac.uk, November 2015.
    
    [numEls,NPE] = size(model.Elements(:,1:3));

    subCell = [];
    subCell = [subCell; model.Elements(1,1) model.Elements(1,2) 1 0 1 0;
            model.Elements(1,2) model.Elements(1,3) 1 0 2 0; 
            model.Elements(1,1) model.Elements(1,3) 1 0 3 0];
    
    for i = 2:numEls
        for j = 1:NPE
            n1 = j;
            if n1 == NPE
                n2 = 1;
            else
                n2 = n1 + 1;
            end
            flag = 0;
            for m = 1:size(subCell,1)
                if (model.Elements(i,n1)==subCell(m,1) & model.Elements(i,n2)==subCell(m,2)) | ...
                    (model.Elements(i,n2)==subCell(m,1) & model.Elements(i,n1)==subCell(m,2))
                    flag = 1;
                    subCell(m,4) = i;
                    subCell(m,6) = j;
                    break; 
                end                         
            end
            if flag == 0
                subCell = [subCell; model.Elements(i,n1) model.Elements(i,n2) i 0 j 0];           
            end            
        end 
    end 
    
end
