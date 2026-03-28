function [N,dNdxi] = getSerendipityShapeFunc_lagrange(elemType,wkX,pt)
	% 
    if strcmp(elemType,'T6')
%         xieta = cal_xieta('T3',pt,wkX);        
%         [phi,dphi] = lagrange_basis('T3',xieta);
%        
%         % construct the pairwise products
%         mu = [phi(1)*phi(1); phi(2)*phi(2); phi(3)*phi(3);
%         phi(1)*phi(2); phi(2)*phi(3); phi(3)*phi(1)];
%         A = eye(6,6);
%         
%         % linear transformation matrix B
%         bb = [-1 0 -1; -1 -1 0; 0 -1 -1];
%         B = [eye(3,3) bb; zeros(3,3) 4*eye(3,3)];
%         
%         % combinations
%         comb = [1 1; 2 2; 3 3; 1 2; 2 3; 3 1];
%     
%         % get serendipity shape functions
%         xi = A*mu;
%     
%         % get serendipity shape functions with Lagrange property
%         psi = B*xi;
%         N = psi;
%     
%         % derivatives of shape functions
%         for i = 1:size(comb,1)
%             in = comb(i,1); jn = comb(i,2);
%         
%             DmuDx(i,1) = dphi(in,1)*phi(jn) + phi(in)*dphi(jn,1);
%             DmuDy(i,1) = dphi(in,2)*phi(jn) + phi(in)*dphi(jn,2);
%         end 
%         dxix = A*DmuDx; dxiy = A*DmuDx;
%         dpsix = B*dxix; dpsiy = B*dxiy;
%         dNdxi = [dpsix dpsiy];
        xieta = cal_xieta(elemType,pt,wkX);
        [N,dNdxi] = lagrange_basis(elemType,xieta);
    elseif strcmp(elemType,'Q8')
        xieta = cal_xieta(elemType,pt,wkX);
        [phi,dphi] = lagrange_basis(elemType,xieta);
        
        % construct the pairwise products
        mu = [phi(1)*phi(1); phi(2)*phi(2); phi(3)*phi(3); phi(4)*phi(4);
            phi(1)*phi(2); phi(2)*phi(3); phi(3)*phi(4); phi(4)*phi(1);
            phi(1)*phi(3); phi(2)*phi(4)];
        A = getAmat_quad(wkX(1:4,:));
        
        % linear transformation matrix B
        bb = [-1 0 0 -1; -1 -1 0 0; 0 -1 -1 0; 0 0 -1 -1];
        B = [eye(4,4) bb; zeros(4,4) 4*eye(4,4)];
        % combinations
        comb = [1 1; 2 2; 3 3; 4 4; 1 2; 2 3 ; 3 4; 4 1; 1 3; 2 4];
    
        % get serendipity shape functions
        xi = A*mu;
    
        % get serendipity shape functions with Lagrange property
        psi = B*xi;
        N = psi;
    
        % derivatives of shape functions
        for i = 1:size(comb,1)
            in = comb(i,1); jn = comb(i,2);
        
            DmuDx(i,1) = dphi(in,1)*phi(jn) + phi(in)*dphi(jn,1);
            DmuDy(i,1) = dphi(in,2)*phi(jn) + phi(in)*dphi(jn,2);
        end 
        dxix = A*DmuDx; dxiy = A*DmuDx;
        dpsix = B*dxix; dpsiy = B*dxiy;
        dNdxi = [dpsix dpsiy];
        
    elseif strcmp(elemType,'Q9')
        xieta = cal_xieta(elemType,pt,wkX);
        [N,dNdxi] = lagrange_basis(elemType,xieta);
       
    end

end
%%
function [A] = getAmat_quad(nds)
    % filename: [A] = getAmat_quad(nds)
    %   purpose: to compute the A matrix that relates the pairwise products
    %   with the serendipity shape functions.
    %
    % Note: this only works for quadrilateral element. The procedure involves
    % mapping to a unit diamater polygon and the number of diagonals vary.
    % Hence, separate routine for each polygon.
    %
    % Inputs:
    %       nds - coordinates of the element. only vertex coordinates are
    %       required.
    %
    % Output:
    %       A - matrix A that relates the pairwise products with serendipity
    %       functions. Ref: Quadratic serendipity shape functions, A Rand, A
    %       Gillette, C Bajaj.
    %
    % Sundar, 2015
    %--------------------------------------------------------------------------
    
    % step1 - find the length of the diagonals...
    % for the case of a generic quadrilateral, as you know, there will be only
    % two :)
    d13 = sqrt( (nds(1,1)-nds(3,1))^2 + (nds(1,2)-nds(3,2))^2) ;
    d24 = sqrt( (nds(2,1)-nds(4,1))^2 + (nds(2,2)-nds(4,2))^2) ;
    
    if abs(d13-d24) < 1e-06;
        c13 = [-1 0 -1 0 1/2 1/2 1/2 1/2];
        c24 = [0 -1 0 -1 1/2 1/2 1/2 1/2];
        A = [eye(8,8) c13' c24'];
        
        return;
    end
    
    % step 2 - coordinates of the centroid of the region
    [egeom,~,~] = polygeom(nds(:,1),nds(:,2));
    
    % step 3 - scale the region to unit diameter polygon
    alpha = 1./max([d13,d24]);
    smatrix = [alpha 0 0;0 alpha 0;0 0 1];
    T1 = [1 0 0;0 1 0;-egeom(2) -egeom(3) 1];
    T2 = [1 0 0;0 1 0;egeom(2) egeom(3) 1];
    
    newc = [nds(:,1) nds(:,2) ones(4,1)]*T1*smatrix*T2;
    v = newc(:,1:2);
    
    masterV = v ;
    
    % step 4 - define a horizontal line, say [-1,1]. This is with which, we
    % will do all rotations...!!!! Note that we dont have to change the
    % coordinates, as we have already scaled the polygon such that its diamter
    % is 1.
    v1 = [-1,0; 1,0] ;
    
    % step 5 - one diagonal at a time, rotate it such that the y coordinate is
    % zero and then compute the coefficients...
    
    %--------------------------------------------------------------------------
    % step 5a - diagonal 1-3
    l1 = [v(1,1)-v(3,1) v(1,2)-v(3,2)] ;
    l2 = [v1(1,1)-v1(2,1) v1(2,2)-v1(2,2)] ;
    
    angle = acos( dot(l1,l2)/(norm(l1)*norm(l2))) ;
    
    T = [cos(angle) sin(angle);-sin(angle) cos(angle)];
    
    % step 5a.1 - rotate the polygon
    for in = 1:size(v,1)
        x = v(in,:) ;
        new(in,:) = (T*x')';
    end
    
    clear v
    v = new;
    
    % step 5a.2 - shift the region such that y-coordinate is zero...
    v(:,2) = v(:,2) - v(1,2);
    
    % step 5a.3 - before proceeding check if its zero. sometimes, we might
    % have to add the coordinates...
    if v(1,2) ~= 0
        new(:,2) = new(:,2) + new(1,2);
        v = new ;
    end
    
    % step 5a.4 - shift the region horizontally such that the coordinates of
    % the end point of the diagonal are at (-l,+1)
    cl = 0.5*( v(1,1)+v(3,1));
    v(:,1) = v(:,1) - cl ;
    
    % step 5a.5 - all set to compute the coefficients for the diagonal 1-3
    ell = abs(v(1,1));
    
    va = [v(4,:); v(1,:); v(2,:)];
    vb = [v(2,:); v(3,:); v(4,:)];
    
    [caa,cbb,caa1,ca1a,cbb1,cb1b]=getCoeff(va,vb,ell);
    
    % step 5a.6 - save the coefficients
    c13 = [caa 0 cbb 0 caa1 cb1b cbb1 ca1a];
    
    clear caa cbb caa1 ca1a cb1b cbb1
    %--------------------------------------------------------------------------
    clear v; 
    v = masterV ;
    
    % step 5b - diagonal 2-4
    l1 = [v(2,1)-v(4,1) v(2,2)-v(4,2)] ;
    l2 = [v1(1,1)-v1(2,1) v1(2,2)-v1(2,2)] ;
    
    angle = acos( dot(l1,l2)/(norm(l1)*norm(l2))) ;
    
    T = [cos(angle) sin(angle);-sin(angle) cos(angle)];
    
    % step 5b.1 - rotate the polygon
    for in = 1:size(v,1)
        x = v(in,:) ;
        new(in,:) = (T*x')';
    end
    
    clear v
    v = new;
    
    % step 5b.2 - shift the region such that y-coordinate is zero...
    v(:,2) = v(:,2) - v(2,2);
    
    % step 5a.3 - before proceeding check if its zero. sometimes, we might
    % have to add the coordinates...
    if v(2,2) ~= 0
        new(:,2) = new(:,2) + new(2,2);
        v = new ;
    end
    
    % step 5b.4 - shift the region horizontally such that the coordinates of
    % the end point of the diagonal are at (-l,+1)
    cl = 0.5*( v(2,1)+v(4,1));
    v(:,1) = v(:,1) - cl ;
    
    % step 5b.5 - all set to compute the coefficients for the diagonal 1-3
    ell = abs(v(2,1));
    
    va = [v(1,:); v(2,:); v(3,:)];
    vb = [v(3,:); v(4,:); v(1,:)];
    
    [caa,cbb,caa1,ca1a,cbb1,cb1b]=getCoeff(va,vb,ell);
    
    % step 5b.6 - save the coefficients
    c24 = [0 caa 0 cbb ca1a caa1 cb1b cbb1] ;
    
    clear caa cbb caa1 ca1a cb1b cbb1
    
    %--------------------------------------------------------------------------
    
    % step 6 - construct the A matrix...
    A = [eye(8,8) c13' c24'];

end