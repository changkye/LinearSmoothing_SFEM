function [N,dNdxi] = lagrange_basis(type,coord,dim)
    % returns the lagrange interpolant basis and its gradients w.r.t the
    % parent coordinate system.
    %
    %         [N(xi),dNdxi(xi)]=lagrange_basis(type-order,coord,dim)
    %
    %   type is the toplogical class of finite element it is in the general
    %   form 'topology-#of nodes' ie a three node triangel is T3 a four
    %   node quadralateral is Q4 a 4 node tetrahedra is H4 a 27 node brick
    %   is B27 etc
    %
    %   coord is the parent coordinates at which the basis and its
    %   gradients are to be evaluated at.
    %
    %   presently defined are L2, L3, T3, T4(cubic bubble), T6, Q4, Q9,
    %   H4, H10, B8 and B27
    %
    %   If dim is set to 2 then the vector representation of the N
    %   matrix is returned.
    %
    % written by Jack Chessa
    %            j-chessa@northwestern.edu
    % Department of Mechanical Engineering
    % Northwestern University

    if nargin == 2; dim = 1; end

    switch type
        case 'L2'
            %%%%%%%%%%%%%%%%%%%%% L2 TWO NODE LINE ELEMENT %%%%%%%%%%%%%%%%%%%%%
            %
            %    1---------2
            %
            if size(coord,2) < 1
                disp('Error coordinate needed for the L2 element')
            else
                x = coord(1);
                % shape functions
                N = [1-x; 1+x]/2;
                % shape function derivatives
                dNdxi = [-1; 1]/2;
            end

        case 'L3'
            %%%%%%%%%%%%%%%%%%% L3 THREE NODE LINE ELEMENT %%%%%%%%%%%%%%%%%%%%%
            %
            %    1---------3----------2
            %
            if size(coord,2) < 1
                disp('Error two coordinates needed for the L3 element')
            else
                x = coord(1);
                % shape functions
                N = [-(1-x)*x/2; (1+x)*x/2; 1-x^2];
                % shape function derivatives
                dNdxi = [x-.5; x+.5; -2*x];
            end

        case 'T3'
            %%%%%%%%%%%%%%%% T3 THREE NODE TRIANGULAR ELEMENT %%%%%%%%%%%%%%%%%%
            %
            %               3
            %             /  \
            %            /    \
            %           /      \
            %          /        \
            %         /          \
            %        /            \
            %       /              \
            %      /                \
            %     /                  \
            %    1--------------------2
            %
            if size(coord,2) < 2
                disp('Error two coordinates needed for the T3 element')
            else
                x = coord(1); y = coord(2);
                % shape functions
                N = [1-x-y; x; y];
                % shape function derivatives
                dNdxi = [-1 -1; 1  0; 0  1];
            end

        case 'T4'
            %%%%%%%%%% T4 FOUR NODE TRIANGULAR CUBIC BUBBLE ELEMENT %%%%%%%%%%%%
            %
            %               3
            %             /  \
            %            /    \
            %           /      \
            %          /        \
            %         /          \
            %        /      4     \
            %       /              \
            %      /                \
            %     /                  \
            %    1--------------------2
            %
            if size(coord,2) < 2
                disp('Error two coordinates needed for the T4 element')
            else
                x = coord(1); y = coord(2);
                I = 9*x*y*(1-x-y); Ix = 9*y*(2*x+y-1); Iy = 9*x*(x+2*y-1);
                % shape functions
                N = [1-x-y-I; x-I; y-I; 3*I];
                % shape function derivatives
                dNdxi = [Ix-1 Iy-1; Ix+1 Iy; Ix Iy+1; -3*Ix -3*Iy];
            end

        case 'T6'
            %%%%%%%%%%%%%%%%%% T6 SIX NODE TRIANGULAR ELEMENT %%%%%%%%%%%%%%%%%%
            %
            %               3
            %             /  \
            %            /    \
            %           /      \
            %          /        \
            %         6          5
            %        /            \
            %       /              \
            %      /                \
            %     /                  \
            %    1---------4----------2
            %
            if size(coord,2) < 2
                disp('Error two coordinates needed for the T6 element')
            else
                x = coord(1); y = coord(2);
                % shape functions
                N = [(1-x-y)*(1-2*x-2*y); x*(2*x-1); y*(2*y-1); 4*x*(1-x-y); 4*x*y; 4*y*(1-x-y)];
                % shape function derivatives
                dNdxi(:,1) = [4*(x+y)-3; 4*x-1; 0; -4*(2*x+y-1); 4*y; -4*y];
                dNdxi(:,2) = [4*(x+y)-3; 0; 4*y-1; -4*x; 4*x; -4*(x+2*y-1)];
            end

        case 'T10'
            %%%%%%%%%%%%%%%%%% T10 TEN NODE TRIANGULAR ELEMENT %%%%%%%%%%%%%%%%%%
            %
            %               3
            %             /  \
            %            /    \
            %           6      8
            %          /        \
            %         /          \
            %        /     10     \
            %       9              5
            %      /                \
            %     /                  \
            %    1-----4--------7-----2
            %
            if size(coord,2) < 2
                disp('Error two coordinates needed for the T6 element')
            else
                x = coord(1); y = coord(2);
                % shape functions
                N(1,1) = 1/2*(1-3*x-3*y)*(2-3*x-3*y)*(1-x-y);
                N(2,1) = 9/2*x*(x-1/3)*(x-2/3);
                N(3,1) = 9/2*y*(y-1/3)*(y-2/3);
                N(4,1) = 9/2*x*(2-3*x-3*y)*(1-x-y);
                N(5,1) = 9/2*x*y*(3*x-1);
                N(6,1) = 9/2*y*(3*y-1)*(1-x-y);
                N(7,1) = 9/2*x*(3*x-1)*(1-x-y);
                N(8,1) = 9/2*x*y*(3*y-1);
                N(9,1) = 9/2*y*(2-3*x-3*y)*(1-x-y);
                N(10,1) = 27*x*y*(1-x-y);
                % shape function derivatives 
                dNdxi(1,:) = -0.5*[9*x*(3*x+6*y-4)+9*y*(3*y-4)+11, 9*y*(3*y+6*x-4)+9*x*(3*x-4)+11];
                dNdxi(2,:) = [9/2*x*(3*x-2)+1, 0];
                dNdxi(3,:) = [0, 9/2*y*(3*y-2)+1];
                dNdxi(4,:) = [9/2*(9*x^2+2*x*(6*y-5)+3*y^2-5*y+2), 27*x*(x+y-5/6)];
                dNdxi(5,:) = [9/2*y*(6*x-1), 9/2*x*(3*x-1)];
                dNdxi(6,:) = [-9/2*y*(3*y-1), -9/2*(3*y*(3*y+2*x-8/3)-x+1)];
                dNdxi(7,:) = [-9/2*(3*x*(3*x+2*y-8/3)-y+1), -9/2*x*(3*x-1)];
                dNdxi(8,:) = [9/2*y*(3*y-1), 9/2*x*(6*y-1)];
                dNdxi(9,:) = [27*y*(y+x-5/6), 9/2*(2*y*(9/2*y+6*x-5)+x*(3*x-5)+2)];
                dNdxi(10,:) = [-27*y*(2*x+y-1), -27*x*(x+2*y-1)];
            end

        case 'Q4'
            %%%%%%%%%%%%%%% Q4 FOUR NODE QUADRILATERIAL ELEMENT %%%%%%%%%%%%%%%%
            %
            %    4--------------------3
            %    |                    |
            %    |                    |
            %    |                    |
            %    |                    |
            %    |                    |
            %    |                    |
            %    |                    |
            %    |                    |
            %    |                    |
            %    1--------------------2
            %
            if size(coord,2) < 2
                disp('Error two coordinates needed for the Q4 element')
            else
                x = coord(1); y = coord(2);
                % shape functions
                N = 1/4*[(1-x)*(1-y); (1+x)*(1-y); (1+x)*(1+y); (1-x)*(1+y)];
                % shape funciton derivatives
                dNdxi = 1/4*[y-1 x-1; -y+1 -x-1; y+1 x+1; -y-1 -x+1];
            end 

        case 'Q8'
            %%%%%%%%%%%%%%% Q8 EIGHT NODE QUADRILATERIAL ELEMENT %%%%%%%%%%%%%%%%
            %
            %    4---------7----------3
            %    |                    |
            %    |                    |
            %    |                    |
            %    |                    |
            %    8                    6
            %    |                    |
            %    |                    |
            %    |                    |
            %    |                    |
            %    1----------5---------2
            %
            if size(coord,2) < 2
                disp('Error two coordinates needed for the Q8 element')
            else
                x = coord(1); y = coord(2);
                % shape functions
                N(1,1) = -1/4*(1-x)*(1-y)*(1+x+y); 
                N(2,1) = -1/4*(1+x)*(1-y)*(1-x+y); 
                N(3,1) = -1/4*(1+x)*(1+y)*(1-x-y);
                N(4,1) = -1/4*(1-x)*(1+y)*(1+x-y);
                N(5,1) =  1/2*(1-x^2)*(1-y);
                N(6,1) =  1/2*(1+x)*(1-y^2);
                N(7,1) =  1/2*(1-x^2)*(1+y);
                N(8,1) =  1/2*(1-x)*(1-y^2);
                % shape function derivatives
                N(:,1) = -1/4*[(y-1)*(2*x+y); (y-1)*(2*x-y); -(y+1)*(2*x+y); -(y+1)*(2*x-y);
                    -4*x*(y-1); 2*(y^2-1); 4*x*(y+1); -2*(y^2-1)];
                N(:,2) = -1/4*[(x-1)*(x+2*y); (x+1)*(x-2*y); -(x+1)*(x+2*y); -(x-1)*(x-2*y);
                    -2*(x^2-1); 4*y*(x+1); 2*(x^2-1); -4*y*(x-1)];
            end

        case 'Q9'
            %%%%%%%%%%%%%%% Q9 NINE NODE QUADRILATERIAL ELEMENT %%%%%%%%%%%%%%%%
            %
            %    4---------7----------3
            %    |                    |
            %    |                    |
            %    |                    |
            %    |                    |
            %    8          9         6
            %    |                    |
            %    |                    |
            %    |                    |
            %    |                    |
            %    1----------5---------2
            %
            if size(coord,2) < 2
                disp('Error two coordinates needed for the Q9 element')
            else
                x = coord(1); y = coord(2);
                % shape functions
                N(1,1) =  1/4*x*y*(1-x)*(1-y);
                N(2,1) = -1/4*x*y*(1+x)*(1-y);
                N(3,1) =  1/4*x*y*(1+x)*(1+y);
                N(4,1) = -1/4*x*y*(1-x)*(1+y);
                N(5,1) = -1/2*(1-x^2)*(1-y)*y;
                N(6,1) =  1/2*(1-y^2)*(1+x)*x;
                N(7,1) =  1/2*(1-x^2)*(1+y)*y;
                N(8,1) = -1/2*(1-y^2)*(1-x)*x;
                N(9,1) =  (1-x^2)*(1-y^2);
                % shape function derivatives
                dNdxi(:,1) = 1/4*[y*(2*x-1)*(y-1); y*(2*x+1)*(y-1); y*(2*x+1)*(y+1); 
                    y*(2*x-1)*(y+1); -4*x*(y-1)*y; -2*(2*x+1)*(y^2-1); 
                    -4*x*(y+1)*y; -2*(2*x-1)*(y^2-1); 8*x*(y^2-1)];
                dNdxi(:,2) = 1/4*[x*(x-1)*(2*y-1); x*(x+1)*(2*y-1); x*(x+1)*(2*y+1); 
                    x*(x-1)*(2*y+1); -2*(x^2-1)*(2*y-1); -4*x*(x+1)*y; 
                    -2*(x^2-1)*(2*y+1); -4*x*(x-1)*y; 8*(x^2-1)*y];
            end

        case 'Q16'
            %%%%%%%%%%%%%%% Q16 SIXTEEN NODE QUADRILATERIAL ELEMENT %%%%%%%%%%%%%%%%
            %
            %    4-----11----7-----3
            %    |                 | 
            %    |                 | 
            %    8     16    15    10
            %    |                 | 
            %    |                 | 
            %    12    13    14    6 
            %    |                 | 
            %    |                 | 
            %    1-----5-----9-----2
            %
            if size(coord,2) < 2
                disp('Error two coordinates needed for the Q9 element')
            else
                x = coord(1); y = coord(2);
                % shape functions
                N(1,1)  = 81/256*(1-x)*(1-y)*(1/9-x^2)*(1/9-y^2);
                N(2,1)  = 81/256*(1+x)*(1-y)*(1/9-x^2)*(1/9-y^2);
                N(3,1)  = 81/256*(1+x)*(1+y)*(1/9-x^2)*(1/9-y^2);
                N(4,1)  = 81/256*(1-x)*(1+y)*(1/9-x^2)*(1/9-y^2);
                N(5,1)  = 243/256*(1-x^2)*(y^2-1/9)*(1/3-3*x)*(1-y);
                N(9,1)  = 243/256*(1-x^2)*(y^2-1/9)*(1/3+3*x)*(1-y);
                N(7,1)  = 243/256*(1-x^2)*(y^2-1/9)*(1/3+3*x)*(1+y);
                N(11,1) = 243/256*(1-x^2)*(y^2-1/9)*(1/3-3*x)*(1+y);
                N(6,1)  = 243/256*(1-y^2)*(x^2-1/9)*(1/3-3*y)*(1+x);
                N(10,1) = 243/256*(1-y^2)*(x^2-1/9)*(1/3+3*y)*(1+x);
                N(8,1)  = 243/256*(1-y^2)*(x^2-1/9)*(1/3+3*y)*(1-x);
                N(12,1) = 243/256*(1-y^2)*(x^2-1/9)*(1/3-3*y)*(1-x);
                N(13,1) = 729/256*(1-x^2)*(1-y^2)*(1/3-3*x)*(1/3-3*y);
                N(14,1) = 729/256*(1-x^2)*(1-y^2)*(1/3+3*x)*(1/3-3*y);
                N(15,1) = 729/256*(1-x^2)*(1-y^2)*(1/3+3*x)*(1/3+3*y);
                N(16,1) = 729/256*(1-x^2)*(1-y^2)*(1/3-3*x)*(1/3+3*y);
                % shape function derivatives
                dNdxi(1,:) =  1/256*[(27*x^2-18*x-1)*(y-1)*(9*y^2-1) (x-1)*(9*x^2-1)*(27*y^2-18*y-1)];
                dNdxi(2,:) = -1/256*[(27*x^2+18*x-1)*(y-1)*(9*y^2-1) (x+1)*(9*x^2-1)*(27*y^2-18*y-1)];
                dNdxi(3,:) =  1/256*[(27*x^2+18*x-1)*(y+1)*(9*y^2-1) (x+1)*(9*x^2-1)*(27*y^2+18*y-1)];
                dNdxi(4,:) = -1/256*[(27*x^2-18*x-1)*(y+1)*(9*y^2-1) (x-1)*(9*x^2-1)*(27*y^2+18*y-1)];
                dNdxi(5,:) = -9/256*[(27*x^2-2*x-9)*(y-1)*(9*y^2-1) (9*x-1)*(x^2-1)*(27*y^2-18*y-1)];
                dNdxi(9,:) =  9/256*[(27*x^2+2*x-9)*(y-1)*(9*y^2-1) (9*x+1)*(x^2-1)*(27*y^2-18*y-1)];
                dNdxi(7,:) = -9/256*[(27*x^2+2*x-9)*(y+1)*(9*y^2-1) (9*x+1)*(x^2-1)*(27*y^2+18*y-1)];
                dNdxi(11,:) =  9/256*[(27*x^2-2*x-9)*(y+1)*(9*y^2-1) (9*x-1)*(x^2-1)*(27*y^2+18*y-1)];
                dNdxi(6,:)  =  9/256*[(27*x^2+18*x-1)*(9*y-1)*(y^2-1) (x+1)*(9*x^2-1)*(27*y^2-2*y-9)];
                dNdxi(10,:) = -9/256*[(27*x^2+18*x-1)*(9*y+1)*(y^2-1) (x+1)*(9*x^2-1)*(27*y^2+2*y-9)];
                dNdxi(8,:)  =  9/256*[(27*x^2-18*x-1)*(9*y+1)*(y^2-1) (x-1)*(9*x^2-1)*(27*y^2+2*y-9)];
                dNdxi(12,:) = -9/256*[(27*x^2-18*x-1)*(9*y-1)*(y^2-1) (x-1)*(9*x^2-1)*(27*y^2-2*y-9)];
                dNdxi(13,:) =  81/256*[(27*x^2-2*x-9)*(9*y-1)*(y^2-1) (9*x-1)*(x^2-1)*(27*y^2-2*y-9)];
                dNdxi(14,:) = -81/256*[(27*x^2+2*x-9)*(9*y-1)*(y^2-1) (9*x+1)*(x^2-1)*(27*y^2-2*y-9)];
                dNdxi(15,:) =  81/256*[(27*x^2+2*x-9)*(9*y+1)*(y^2-1) (9*x+1)*(x^2-1)*(27*y^2+2*y-9)];
                dNdxi(16,:) = -81/256*[(27*x^2-2*x-9)*(9*y+1)*(y^2-1) (9*x-1)*(x^2-1)*(27*y^2+2*y-9)];
            end

        case 'H4'
            %%%%%%%%%%%%%%%% H4 FOUR NODE TETRAHEDRAL ELEMENT %%%%%%%%%%%%%%%%%%
            %
            %             4
            %           / | \
            %          /  |  \
            %         /   |   \
            %        /    |    \
            %       /     |     \
            %      1 -----|------3
            %         -   2  -
            if size(coord,2) < 3
                disp ('Error three coordinates needed for the H4 element')
            else
                x = coord(1); y = coord(2); z = coord(3);
                % shape functions
                N = [1-x-y-z; x; y; z];
                % shape function derivatives
                dNdxi=[-1 -1 -1; 1 0 0; 0 1 0; 0 0 1];
            end

        case 'H4b'
            %%%%%%%%%%%%%%%% H4 FOUR NODE TETRAHEDRAL ELEMENT WITH BUBBLE %%%%%%%%%%%%%%%%%%
            %
            %             4
            %           / | \
            %          /  |  \
            %         /   |   \
            %        /    |    \
            %       /     |     \
            %      1 -----|------3
            %         -   2  -
            if size(coord,2) < 3
                disp ('Error three coordinates needed for the H4 element')
            else
                % 
                x = coord(1); y = coord(2); z = coord(3);
                I = 64*x*y*z*(1-x-y-z);
                Ix = 64*y*z*(2*x+y+z-1); Iy = 64*x*z*(x+2*y+z-1); Iz = 64*x*y*(x+y+2*z-1);
                % shape functions 
                N = [1-x-y-z-I; x-I; y-I; z-I; 4*I];
                % shape function derivatives
                dNdxi = [Ix-1 Iy-1 Iz-1; Ix+1 Iy Iz; Ix Iy+1 Iz; Ix Iy Iz+1; -4*Ix -4*Iy -4*Iz];
            end

        case 'B8'
            %%%%%%%%%%%%%%%%%%% B8 EIGHT NODE BRICK ELEMENT %%%%%%%%%%%%%%%%%%%%
            %
            %                  8
            %               /    \
            %            /          \
            %         /                \
            %      5                     \
            %      |\                     7
            %      |   \                / |
            %      |     \     4    /     |
            %      |        \    /        |
            %      |           6          |
            %      1           |          |
            %       \          |          3
            %          \       |        /
            %            \     |     /
            %               \  |  /
            %                  2
            %
            if size(coord,2) < 3
                disp ('Error three coordinates needed for the B8 element')
            else
                x = coord(1); y = coord(2); z = coord(3);
                I1 = 1/2 - coord/2; I2 = 1/2 + coord/2;
                % shape functions
                N = [I1(1)*I1(2)*I1(3); I2(1)*I1(2)*I1(3); I2(1)*I2(2)*I1(3); I1(1)*I2(2)*I1(3);
                    I1(1)*I1(2)*I2(3); I2(1)*I1(2)*I2(3); I2(1)*I2(2)*I2(3); I1(1)*I2(2)*I2(3)];
                % shape function derivatives
                dNdxi(:,1) = [-1+y+z-y*z; 1-y-z+y*z; 1+y-z-y*z; -1-y+z+y*z; -1+y-z+y*z; 1-y+z-y*z;
                    1+y+z+y*z; -1-y-z-y*z]/8;
                dNdxi(:,2) = [-1+x+z-x*z; -1-x+z+x*z; 1+x-z-x*z; 1-x-z+x*z; -1+x-z+x*z; -1-x-z-x*z;
                    1+x+z+x*z; 1-x+z-x*z]/8;
                dNdxi(:,3) = [-1+x+y-x*y; -1-x+y+x*y; -1-x-y-x*y; -1+x-y+x*y; 1-x-y+x*y; 1+x-y-x*y;
                    1+x+y+x*y; 1-x+y-x*y]/8;
            end

        case 'B20'
            %%%%%%%%%%%%%%%%%%% B20 TWENTY NODE BRICK ELEMENT %%%%%%%%%%%%%%%%%%%%
            %
            %                  8
            %               /    \
            %            16         15
            %         /                \
            %      5          20         \
            %      |\                     7
            %      |   \                / |
            %      17    13    4    14    |
            %      |        \    /        19
            %      |     12    6     11   |
            %      1           |          |
            %       \          |          3
            %          9       18       /
            %            \     |    10
            %               \  |  /
            %                  2
            %
            if size(coord,2) < 3
                disp ('Error three coordinates needed for the B20 element')
            else
                N = zeros(20,1); dNdxi = zeros(20,3);
                % 
                x = coord(1); y = coord(2); z = coord(3);
                % shape functions
                N = 1/8*[-(1-x)*(1-y)*(1-z)*(x+y+z+2); (1+x)*(1-y)*(1-z)*(x-y-z-2);
                        (1+x)*(1-y)*(1+z)*(x-y+z-2); -(1-x)*(1-y)*(1+z)*(x+y-z+2);
                        -(1-x)*(1+y)*(1-z)*(x-y+z+2); (1+x)*(1+y)*(1-z)*(x+y-z-2);
                        (1+x)*(1+y)*(1+z)*(x+y+z-2); -(1-x)*(1+y)*(1+z)*(x-y-z+2);
                        2*(1-x^2)*(1-y)*(1-z); 2*(1-z^2)*(1+x)*(1-y); 2*(1-x^2)*(1-y)*(1+z);
                        2*(1-z^2)*(1-x)*(1-y); 2*(1-x^2)*(1+y)*(1-z); 2*(1-z^2)*(1+x)*(1+y);
                        2*(1-x^2)*(1+y)*(1+z); 2*(1-z^2)*(1-x)*(1+y); 2*(1-y^2)*(1-x)*(1-z);
                        2*(1-y^2)*(1+x)*(1-z); 2*(1-y^2)*(1+x)*(1+z); 2*(1-y^2)*(1-x)*(1+z)];
                % shape function derivatives
                dNdxi(1,:) = 1/8*[(y-1)*(z-1)*(2*x+y+z+1) (x-1)*(z-1)*(x+2*y+z+1) (x-1)*(y-1)*(x+y+2*z+1)];
                dNdxi(2,:) = 1/8*[-(y-1)*(z-1)*(-2*x+y+z+1) (x+1)*(z-1)*(x-2*y-z-1) (x+1)*(y-1)*(x-y-2*z-1)];
                dNdxi(3,:) = 1/8*[(y-1)*(z+1)*(-2*x+y-z+1) -(x+1)*(z+1)*(x-2*y+z-1) -(x+1)*(y-1)*(x-y+2*z-1)];
                dNdxi(4,:) = 1/8*[-(y-1)*(z+1)*(2*x+y-z+1) -(x-1)*(z+1)*(x+2*y-z+1) -(x-1)*(y-1)*(x+y-2*z+1)];
                dNdxi(5,:) = 1/8*[(y+1)*(z-1)*(-2*x+y-z-1) -(x-1)*(z-1)*(x-2*y+z+1) -(x-1)*(y+1)*(x-y+2*z+1)];
                dNdxi(6,:) = 1/8*[-(y+1)*(z-1)*(2*x+y-z-1) -(x+1)*(z-1)*(x+2*y-z-1) -(x+1)*(y+1)*(x+y-2*z-1)];
                dNdxi(7,:) = 1/8*[(y+1)*(z+1)*(2*x+y+z-1) (x+1)*(z+1)*(x+2*y+z-1) (x+1)*(y+1)*(x+y+2*z-1)];
                dNdxi(8,:) = 1/8*[-(y+1)*(z+1)*(-2*x+y+z-1) (x-1)*(z+1)*(x-2*y-z+1) (x-1)*(y+1)*(x-y-2*z+1)];
                dNdxi(9,:)  = -1/4*[2*x*(y-1)*(z-1) (x^2-1)*(z-1)   (x^2-1)*(y-1)];
                dNdxi(10,:) =  1/4*[(y-1)*(z^2-1)   (x+1)*(z^2-1)   2*(x+1)*(y-1)*z];
                dNdxi(11,:) =  1/4*[2*x*(y-1)*(z+1) (x^2-1)*(z+1)   (x^2-1)*(y-1)];
                dNdxi(12,:) = -1/4*[(y-1)*(z^2-1)   (x-1)*(z^2-1)   2*(x-1)*(y-1)*z];
                dNdxi(13,:) =  1/4*[2*x*(y+1)*(z-1) (x^2-1)*(z-1)   (x^2-1)*(y+1)];
                dNdxi(14,:) = -1/4*[(y+1)*(z^2-1)   (x+1)*(z^2-1)   2*(x+1)*(y+1)*z];
                dNdxi(15,:) = -1/4*[2*x*(y+1)*(z+1) (x^2-1)*(z+1)   (x^2-1)*(y+1)];
                dNdxi(16,:) =  1/4*[(y+1)*(z^2-1)   (x-1)*(z^2-1)   2*(x-1)*(y+1)*z];
                dNdxi(17,:) = -1/4*[(y^2-1)*(z-1)   2*(x-1)*y*(z-1) (x-1)*(y^2-1)];
                dNdxi(18,:) =  1/4*[(y^2-1)*(z-1)   2*(x+1)*y*(z-1) (x+1)*(y^2-1)];
                dNdxi(19,:) = -1/4*[(y^2-1)*(z+1)   2*(x+1)*y*(z+1) (x+1)*(y^2-1)];
                dNdxi(20,:) =  1/4*[(y^2-1)*(z+1)   2*(x-1)*y*(z+1) (x-1)*(y^2-1)];
            end

        otherwise
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            disp (['Element ',type,' not yet supported'])
            N = [];
            dNdxi = [];
    end

    I = eye(dim);
    Nv = [];
    for i = 1:size(N,1)
        Nv = [Nv; I*N(i)];
    end

end