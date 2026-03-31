function model = meshGeneration2D(model)
    % 
    % number of nodes
    if strcmp(model.elemType,'T3')
        xnode = model.numEls(1) + 1;
        ynode = model.numEls(2) + 1;

        % nodal coordinates
        [x,y] = meshgrid(linspace(model.L(1,1),model.L(1,2),xnode),...
            linspace(model.L(2,1),model.L(2,2),ynode));
        x = x'; y = y';
        coord = [x(:) y(:)];
        
        % element connectivities
        idx = reshape(1:prod([xnode,ynode]),[xnode,ynode]);
        v1 = idx(1:end-1,1:end-1);
        v2 = idx(1:end-1,2:end);
        v3 = idx(2:end,1:end-1);
        v4 = idx(2:end,2:end);
        v1 = v1(:); v2 = v2(:); v3 = v3(:); v4 = v4(:); 
        conn = [v1 v3 v2; v3 v4 v2];
    elseif strcmp(model.elemType,'T4')
        xnode = model.numEls(1) + 1;
        ynode = model.numEls(2) + 1;

        % nodal coordinates
        [x,y] = meshgrid(linspace(model.L(1,1),model.L(1,2),xnode),...
            linspace(model.L(2,1),model.L(2,2),ynode));
        x = x'; y = y';
        coord = [x(:) y(:)];
        
        % element connectivities
        idx = reshape(1:prod([xnode,ynode]),[xnode,ynode]);
        v1 = idx(1:end-1,1:end-1);
        v2 = idx(1:end-1,2:end);
        v3 = idx(2:end,1:end-1);
        v4 = idx(2:end,2:end);
        v1 = v1(:); v2 = v2(:); v3 = v3(:); v4 = v4(:); 
        conn = [v1 v3 v2; v3 v4 v2];

        % addition of nodes at the centre of each element for the bubble
        [numelem,~] = size(conn);
        numnode = size(coord,1);
        
        % additional nodes
        for ie = 1:numelem
            coord(numnode+ie,:) = mean(coord(conn(ie,:),:));
        end
        
        % element connectivities for the added nodes
        for ie = 1:numelem
            conn(ie,4) = numnode + ie;
        end
    elseif strcmp(model.elemType,'T6')
        xnode = 2*model.numEls(1) + 1;
        ynode = 2*model.numEls(2) + 1;

        % nodal coordinates
        [x,y] = meshgrid(linspace(model.L(1,1),model.L(1,2),xnode),...
            linspace(model.L(2,1),model.L(2,2),ynode));
        x = x'; y = y';
        coord = [x(:) y(:)];

        % element connectivity
        ic = 2; jc = 2*xnode;
        idx1 = [1 3 2*xnode+1 2 xnode+2 xnode+1];
        idx2 = [3 2*xnode+3 2*xnode+1 xnode+3 2*xnode+2 xnode+2];
        conn = [makeConnectivity(model,ic,jc,idx1);
            makeConnectivity(model,ic,jc,idx2)];

    elseif strcmp(model.elemType,'Q4')
        xnode = model.numEls(1) + 1;
        ynode = model.numEls(2) + 1;
        
        % nodal coordinates
        [x,y] = meshgrid(linspace(model.L(1,1),model.L(1,2),xnode),...
            linspace(model.L(2,1),model.L(2,2),ynode));
        x = x'; y = y';
        coord = [x(:) y(:)];
        
        % element connectivities
        idx = reshape(1:prod([xnode,ynode]),[xnode,ynode]);
        v1 = idx(1:end-1,1:end-1);
        v2 = idx(1:end-1,2:end);
        v3 = idx(2:end,1:end-1);
        v4 = idx(2:end,2:end);
        v1 = v1(:); v2 = v2(:); v3 = v3(:); v4 = v4(:); 
        conn = [v1 v3 v4 v2];
    elseif strcmp(model.elemType,'Q8')
        xnode = model.numEls(1) + 1;
        ynode = model.numEls(2) + 1;

        % nodal coordinates
        for j = 1:ynode
            for i = 1:xnode
                elem = (3*model.numEls(1)+2)*(j-1)+(2*i-1);
                coord(elem,1) = (model.L(1,2)/model.numEls(1))*(i-1);
                coord(elem,2) = (model.L(2,2)/model.numEls(2))*(j-1);
            end
        end

        % element connectivity
        index = 1;
        for i = 1:model.numEls(2)
            for k = 1:model.numEls(1)
                conn(index,1) = (3*model.numEls(1)+2)*(i-1) + (2*k-1);
                conn(index,2) = (3*model.numEls(1)+2)*(i-1) + (2*k+1);
                conn(index,3) = (3*model.numEls(1)+2)*i+(2*k+1);
                conn(index,4) = (3*model.numEls(1)+2)*i+(2*k-1);
                conn(index,5) = mean(conn(index,[1,2]));
                conn(index,6) = (3*model.numEls(1)+2)*i - model.numEls(1) + k;
                conn(index,7) = mean(conn(index,[3,4]));
                conn(index,8) = (3*model.numEls(1)+2)*i - (model.numEls(1)+1)+ k;
                index = index + 1;
            end
        end
        for i = 1:size(conn,1)
            coord(conn(i,5),:) = mean(coord(conn(i,1:2),:));
            coord(conn(i,6),:) = mean(coord(conn(i,2:3),:));
            coord(conn(i,7),:) = mean(coord(conn(i,3:4),:));
            coord(conn(i,8),:) = mean(coord(conn(i,[4,1]),:));
        end
    elseif strcmp(model.elemType,'Q9')
        xnode = 2*model.numEls(1) + 1;
        ynode = 2*model.numEls(2) + 1;

        % nodal coordinates
        [x,y] = meshgrid(linspace(model.L(1,1),model.L(1,2),xnode),...
            linspace(model.L(2,1),model.L(2,2),ynode));
        x = x'; y = y';
        coord = [x(:) y(:)];

        % element connectivity
        ic = 2; jc = 2*xnode;
        idx = [1 3 2*xnode+3 2*xnode+1 2 xnode+3 2*xnode+2 xnode+1 xnode+2];
        conn = makeConnectivity(model,ic,jc,idx);

    end
   
    model.Nodes = coord;
    model.Elements = conn;

end

%%
function element = makeConnectivity(model,ic,jc,idx)
    % 
    num_u = model.numEls(1); num_v = model.numEls(2);
    
    if nargin < 4
        disp('Not enough parameters specified for make_elem function')
    end

    inc = zeros(1,size(idx,2));
    e = 1;
    element = zeros(num_u*num_v,size(idx,2));

    for row = 1:num_v
        for col = 1:num_u
            element(e,:) = idx + inc;
            inc = inc + ic;
            e = e + 1;
        end
        inc = row*jc;
    end
end
