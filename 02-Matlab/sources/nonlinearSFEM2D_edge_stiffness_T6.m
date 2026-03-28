function model = nonlinearSFEM2D_edge_stiffness_T6(model,uu)
    % nonlinear edge-based smoothed tangent stiffness for T6 elements
    % moment-consistent recovery with quadratic basis:
    % p = [1, x, y, x^2, x*y, y^2]^T
    %
    % This version redesigns recovery to solve moment equations
    %   M * a = f
    % with Tikhonov regularization for robustness.
    %
    % Inputs:
    %   model, uu
    % Optional fields in model:
    %   edgeBasisOrder      : currently supports 2 (default)
    %   edgeBoundaryNg      : boundary Gauss points on L3 edges (default 3)
    %   edgeTriQuadOrder    : triangle interior quadrature order (default 4)
    %   edgeQuadQuadOrder   : quad interior Gauss order (default 3)
    %   edgeRegParam        : Tikhonov scale (default 1e-10)

    if nargin < 2
        error('nonlinearSFEM2D_edge_stiffness_T6:NotEnoughInputs', ...
            'Use nonlinearSFEM2D_edge_stiffness_T6(model,uu).');
    end

    % options
    basisOrder = getfield_with_default(model,'edgeBasisOrder',2); %#ok<GFLD>
    if basisOrder ~= 2
        warning('nonlinearSFEM2D_edge_stiffness_T6:BasisOrderReset', ...
            'edgeBasisOrder=%d is not supported here. Using 2.', basisOrder);
        basisOrder = 2;
    end
    ngb = getfield_with_default(model,'edgeBoundaryNg',3); %#ok<GFLD>
    triOrder = getfield_with_default(model,'edgeTriQuadOrder',4); %#ok<GFLD>
    quadOrder = getfield_with_default(model,'edgeQuadQuadOrder',3); %#ok<GFLD>
    regParam = getfield_with_default(model,'edgeRegParam',1e-10); %#ok<GFLD>

    % Initialisation
    ndof = 2*size(model.Nodes,1);
    K = sparse(ndof,ndof);
    R = sparse(ndof,1);
    EE = 0.0;

    % target edge and sharing elements
    model = getTargetEdge(model);
    targetEdge = model.targetEdge;

    % boundary quadrature
    [Wb,Qb] = quadrature(ngb,'GAUSS',1);

    % basis dimension for quadratic basis
    m = 6;
    nEdge = size(targetEdge,1);
    useParfor = isfield(model,'useParfor') && model.useParfor && ...
        license('test','Distrib_Computing_Toolbox');

    if useParfor
        Iall = cell(nEdge,1);
        Jall = cell(nEdge,1);
        Vall = cell(nEdge,1);
        RIall = cell(nEdge,1);
        RVall = cell(nEdge,1);
        EEall = zeros(nEdge,1);

        parfor ivo = 1:nEdge
            [edof,Ke,Re,EEe] = localEdgeContrib(model,uu,targetEdge,ivo,m,ngb,Wb,Qb,triOrder,quadOrder,regParam);
            [ii,jj] = ndgrid(edof,edof);
            Iall{ivo} = ii(:);
            Jall{ivo} = jj(:);
            Vall{ivo} = Ke(:);
            RIall{ivo} = edof(:);
            RVall{ivo} = Re(:);
            EEall(ivo) = EEe;
        end

        K = sparse(vertcat(Iall{:}),vertcat(Jall{:}),vertcat(Vall{:}),ndof,ndof);
        R = sparse(vertcat(RIall{:}),ones(sum(cellfun(@numel,RIall)),1),vertcat(RVall{:}),ndof,1);
        EE = sum(EEall);
    else
        for ivo = 1:nEdge
            [edof,Ke,Re,EEe] = localEdgeContrib(model,uu,targetEdge,ivo,m,ngb,Wb,Qb,triOrder,quadOrder,regParam);
            K(edof,edof) = K(edof,edof) + Ke;
            R(edof,1) = R(edof,1) + Re;
            EE = EE + EEe;
        end
    end

    model.K = K;
    model.R = R;
    model.EE = EE;
end

function p = basisVec(x,y)
    p = [1; x; y; x*x; x*y; y*y];
end

function [edof,Ke,Re,EEe] = localEdgeContrib(model,uu,targetEdge,ivo,m,ngb,Wb,Qb,triOrder,quadOrder,regParam)
    % neighboring elements sharing this target edge
    if targetEdge(ivo,end-2) == 0
        neighbour = targetEdge(ivo,3);
    else
        neighbour = targetEdge(ivo,3:end-2);
    end
    nc = length(neighbour);

    % accumulate boundary moments from each contributing element
    for ic = 1:nc
        nsf = length(model.supp{neighbour(ic)});
        wkInd = model.Elements(neighbour(ic),:);
        wkX = model.Nodes(wkInd,:);

        side = cal_side(wkX(1:3,1),wkX(1:3,2));
        [nx,ny] = cal_nx_ny(wkX(1:3,1),wkX(1:3,2),side);

        fx = zeros(m,size(wkX,1));
        fy = zeros(m,size(wkX,1));
        bound = [1 2 4; 2 3 5; 3 1 6];

        for is = 1:size(bound,1)
            bxy = zeros(m,size(wkX,1));
            X = wkX(bound(is,:),:);
            for ig = 1:ngb
                [Ng,dNdxi] = lagrange_basis('L3',Qb(ig));
                xieta_gp = Ng'*X;
                detJ = norm(dNdxi'*X);

                N_T = getSerendipityShapeFunc_lagrange('T6',wkX,xieta_gp);
                N = N_T'*detJ*Wb(ig);
                p = basisVec(xieta_gp(1),xieta_gp(2));
                bxy = bxy + p*N;
            end
            fx = fx + nx(is)*bxy;
            fy = fy + ny(is)*bxy;
        end

        if ic == 1
            nodL = model.supp{neighbour(ic)};
            nn = nsf;
            Fx = fx;
            Fy = fy;
        else
            i0 = 0;
            for jj = 1:nsf
                nod = model.supp{neighbour(ic)}(jj);
                flag = 0;
                for j = 1:nn
                    if nodL(j) == nod
                        Fx(:,j) = Fx(:,j) + fx(:,jj);
                        Fy(:,j) = Fy(:,j) + fy(:,jj);
                        flag = 1;
                        break;
                    end
                end
                if flag == 0
                    i0 = i0 + 1;
                    nodL(nn+i0) = nod;
                    Fx(:,nn+i0) = fx(:,jj);
                    Fy(:,nn+i0) = fy(:,jj);
                end
            end
            nn = nn + i0; %#ok<NASGU>
        end
    end

    % reorder into local smoothing-domain numbering
    if nc == 1
        node_sc = [1 2 3 4 5 6];
        element = 'T6';
        [Wi,Qi] = quadrature(triOrder,'TRIANGULAR',2);
    else
        if targetEdge(ivo,5) == 1
            node_sc = [1 7 2 3 9 8 5 6 4];
        elseif targetEdge(ivo,5) == 2
            node_sc = [1 2 7 3 4 8 9 6 5];
        else
            node_sc = [1 2 3 7 4 5 8 9 6];
        end
        element = 'Q9';
        [Wi,Qi] = quadrature(quadOrder,'GAUSS',2);
    end
    nodL = nodL(node_sc);
    gcoord = model.Nodes(nodL,:);
    Fx = Fx(:,node_sc);
    Fy = Fy(:,node_sc);

    % interior moment integrals
    nGp = size(Wi,1);
    Ni = zeros(5,size(gcoord,1));
    M = zeros(m,m);
    P = zeros(m,nGp);
    mW = zeros(nGp,1);

    for ig = 1:nGp
        if nc == 1
            [N1,dN1dxi] = lagrange_basis('T6',Qi(ig,:));
            detJ0 = det(dN1dxi'*gcoord(1:6,:));
            mQ = N1'*gcoord(1:6,:);
        else
            [N1,dN1dxi] = lagrange_basis('Q9',Qi(ig,:));
            detJ0 = det(dN1dxi'*gcoord(1:9,:));
            mQ = N1'*gcoord(1:9,:);
        end
        w = Wi(ig)*detJ0;
        mW(ig) = w;

        xg = mQ(1); yg = mQ(2);
        p = basisVec(xg,yg);
        P(:,ig) = p;
        M = M + w*(p*p');

        N = getSerendipityShapeFunc_lagrange(element,gcoord,mQ);
        Ni(1,:) = Ni(1,:) + (N*w)';
        Ni(2,:) = Ni(2,:) + 2*(N*w*xg)';
        Ni(3,:) = Ni(3,:) + 2*(N*w*yg)';
        Ni(4,:) = Ni(4,:) + (N*w*xg)';
        Ni(5,:) = Ni(5,:) + (N*w*yg)';
    end

    % divergence theorem correction terms
    Fx(2,:) = Fx(2,:) - Ni(1,:);
    Fx(4,:) = Fx(4,:) - Ni(2,:);
    Fx(5,:) = Fx(5,:) - Ni(5,:);
    Fy(3,:) = Fy(3,:) - Ni(1,:);
    Fy(5,:) = Fy(5,:) - Ni(4,:);
    Fy(6,:) = Fy(6,:) - Ni(3,:);

    alpha = regParam * max(trace(M)/m,1.0);
    Mreg = M + alpha*eye(m);
    ax = Mreg\Fx;
    ay = Mreg\Fy;

    dx = (P')*ax;
    dy = (P')*ay;

    nndof = length(nodL);
    wkU = zeros(nndof,2);
    wkU(:,1) = uu(2*nodL-1,1);
    wkU(:,2) = uu(2*nodL,1);

    edof = zeros(1,2*nndof);
    edof(1:2:end) = 2*nodL - 1;
    edof(2:2:end) = 2*nodL;
    Ke = zeros(2*nndof,2*nndof);
    Re = zeros(2*nndof,1);
    W0_acc = 0.0;
    w_acc = 0.0;
    for ig = 1:nGp
        [Bmat,Bgeo,Fmat] = getNonlinearBmat2D([dx(ig,:)',dy(ig,:)'],wkU,nndof);
        [Cmat,Smat,W0] = nonlinearConstitutive(model.param,Fmat);

        Ke = Ke + (Bmat'*Cmat*Bmat + Bgeo'*Smat*Bgeo)*mW(ig);
        Re = Re + Bmat'*[Smat(1,1); Smat(2,2); Smat(1,2)]*mW(ig);
        W0_acc = W0_acc + W0*mW(ig);
        w_acc = w_acc + mW(ig);
    end

    EEe = 0.0;
    if w_acc > 0
        EEe = (W0_acc/w_acc)*model.subA(ivo);
    end
end

function v = getfield_with_default(s,field,defaultValue)
    if isfield(s,field)
        v = s.(field);
    else
        v = defaultValue;
    end
end
