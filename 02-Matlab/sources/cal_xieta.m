function coord = cal_xieta(elemType,xy_g,nodes1)

    epsilon = 0.00001;
%     if strcmp(elemType,'Q4') | strcmp(elemType,'Q8') | strcmp(elemType,'Q9')
%         nodes = nodes1(1:4,:); 
%         element = 'Q4';
%     elseif strcmp(elemType,'Q8')
%         nodes = nodes1;
%         element = 'Q8';
    % elseif strcmp(elemType,'T3') | strcmp(elemType,'T4') |...
    %         strcmp(elemType,'T6')
    %     nodes = nodes1(1:3,:);
    %     element = 'T3';
%     elseif strcmp(elemType,'T3') | strcmp(elemType,'T4')
%         nodes = nodes1(1:3,:);
%         element = 'T3';
%     elseif strcmp(elemType,'T6')
%         nodes = nodes1;
%         element = 'T6';
%     end 
    
    element = elemType;
    nodes = nodes1;

    coord = zeros(1,2);
    ksi = 0;
    eta = 0;
    iter = 10;
    
    inc = 1;

    while inc < iter
        [N,dNdxi] = lagrange_basis(element,coord);  % compute shape functions

        x = N'*nodes(:,1);
        y = N'*nodes(:,2);
        dx_dxi = dNdxi(:,1)'*nodes(:,1);
        dx_deta = dNdxi(:,2)'*nodes(:,1);
    
        dy_dxi = dNdxi(:,1)'*nodes(:,2);
        dy_deta = dNdxi(:,2)'*nodes(:,2);
    
        deltaX = x - xy_g(1);
        deltaY = y - xy_g(2);
        delta = [deltaX; deltaY];
        F = [dx_dxi dx_deta; dy_dxi dy_deta];
        invF = inv(F);
    
        ksi = ksi - invF(1,:)*delta;
        eta = eta - invF(2,:)*delta;
    
        coord(1) = ksi;
        coord(2) = eta;
%         if (abs(ksi - coord(1)) < epsilon) && ...
%             (abs(eta - coord(2)) < epsilon)
%             inc  = iter + 1;
%             xieta = coord;
%         else
            inc = inc + 1;
%         end
    end

end
