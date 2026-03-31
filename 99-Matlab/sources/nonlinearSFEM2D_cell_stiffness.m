function model = nonlinearSFEM2D_cell_stiffness(model,uu)
	% compute nonlinear 2D cell-based smoothing tangent stiffness matrix
	% 	with the linear smoothing function
    %   -> no subcell division
	% 
	% Changkye Lee
	% DSOC National Research Center for Disaster-free & Safe Ocean City,
	% Dong-A University, Korea.
	% changkyelee@gmail.com, April 2019.
	
	% initialisation
    ndof = 2*size(model.Nodes,1);
	K = sparse(ndof,ndof);
    R = sparse(ndof,1);
    EE = 0.0;

    % Gauss points for internal & boundary
    ng = 2;
    [Wi,Qi] = quadrature(ng,'TRIANGULAR',2);
    [Wb,Qb] = quadrature(ng,'GAUSS',1);

    numElem = size(model.Elements,1);
    isLowMemoryMode = isfield(model,'runMode') && strcmpi(model.runMode,'low_memory');
    useParfor = ~isLowMemoryMode && isfield(model,'useParfor') && model.useParfor && ...
        license('test','Distrib_Computing_Toolbox');

    if useParfor
        Iall = cell(numElem,1);
        Jall = cell(numElem,1);
        Vall = cell(numElem,1);
        RIall = cell(numElem,1);
        RVall = cell(numElem,1);
        EEall = zeros(numElem,1);

        parfor ivo = 1:numElem
            [edof,Ke,Re,EEe] = localElementContrib(model,uu,ivo,Wi,Qi,Wb,Qb,ng);
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
        for ivo = 1:numElem
            [edof,Ke,Re,EEe] = localElementContrib(model,uu,ivo,Wi,Qi,Wb,Qb,ng);
            K(edof,edof) = K(edof,edof) + Ke;
            R(edof,1) = R(edof,1) + Re;
            EE = EE + EEe;
        end
    end

    model.K = K;
    model.R = R;
    model.EE = EE/length(Wi);

end

function [edof,Ke,Re,EEe] = localElementContrib(model,uu,ivo,Wi,Qi,Wb,Qb,ng)
    wkInd = model.Elements(ivo,:);
    wkX = model.Nodes(wkInd,:);
    nndof = length(wkInd);

    % current element displacements
    wkU = zeros(nndof,2);
    wkU(:,1) = uu(2*wkInd-1,1);
    wkU(:,2) = uu(2*wkInd,1);

    % global dof indices
    edof = zeros(2*nndof,1);
    edof(1:2:end) = 2*wkInd - 1;
    edof(2:2:end) = 2*wkInd;

    % create subcells
    gcoord = wkX;
    subTri = [1 2 3 4 5 6];
    Asc = polyarea(gcoord(:,1),gcoord(:,2));

    % compute shape functions at internal Gauss points
    xy = gcoord;
    Ni = zeros(1,size(wkX,1));
    Wmat = zeros(3,length(Wi));
    detJ0 = zeros(length(Wi),1);
    mR = zeros(length(Wi),2);
    for ig = 1:length(Wi)
        [N1,dN1dxi] = lagrange_basis('T3',Qi(ig,:));
        detJ0(ig) = det(dN1dxi'*xy(1:3,:));
        mR(ig,:) = N1'*xy(1:3,:);
        N = getSerendipityShapeFunc_lagrange('T6',wkX,mR(ig,:));
        Ni = Ni + (N*Wi(ig)*detJ0(ig))';
        Wmat(:,ig) = Wi(ig)*detJ0(ig)*[1; mR(ig,1); mR(ig,2)];
    end

    % construct smoothed shape functions
    bound = [1 2 4; 2 3 5; 3 1 6];

    % outward normal vectors
    [nx,ny] = cal_nx_ny(xy(1:3,1),xy(1:3,2),cal_side(xy(1:3,1),xy(1:3,2)));
    fx = zeros(3,size(wkX,1));
    fy = zeros(3,size(wkX,1));
    for is = 1:size(bound,1)
        bxy = zeros(3,size(wkX,1));
        node_sc = subTri(bound(is,:)); %#ok<NASGU>
        for ig = 1:ng
            [Ng,dNgdxi] = lagrange_basis('L3',Qb(ig));
            xieta_gp = Ng'*xy(bound(is,:),:);
            detJ = norm(dNgdxi'*xy(bound(is,:),:));
            N_T = getSerendipityShapeFunc_lagrange('T6',wkX,xieta_gp);
            bxy = bxy + [N_T'*detJ*Wb(ig); N_T'*detJ*Wb(ig)*xieta_gp(1); ...
                N_T'*detJ*Wb(ig)*xieta_gp(2)];
        end
        fx = fx + nx(is)*bxy;
        fy = fy + ny(is)*bxy;
    end
    fx(2,:) = fx(2,:) - Ni;
    fy(3,:) = fy(3,:) - Ni;

    % get derivatives basis functions
    dx = Wmat\fx;
    dy = Wmat\fy;

    % local tangent stiffness and internal force
    Ke = zeros(2*nndof,2*nndof);
    Re = zeros(2*nndof,1);
    EEe = 0.0;
    for ig = 1:length(Wi)
        [Bmat,Bgeo,Fmat] = getNonlinearBmat2D([dx(ig,:)',dy(ig,:)'],wkU,nndof);
        [Cmat,Smat,W0] = nonlinearConstitutive(model.param,Fmat);
        w = Wi(ig)*detJ0(ig);

        Ke = Ke + (Bmat'*Cmat*Bmat + Bgeo'*Smat*Bgeo)*w;
        Re = Re + Bmat'*[Smat(1,1);Smat(2,2);Smat(1,2)]*w;
        EEe = EEe + W0*Asc;
    end
end
