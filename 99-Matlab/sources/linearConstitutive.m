function model = linearConstitutive(model)
    % ConputeConstitutive_Linear : compute constitutive matrix for linear problem
    %
    % Changkye Lee, Dept. of Mechanical Engineering,
    % CISTIB, The University of Sheffield,
    % changkye.lee@sheffield.ac.uk, November 2015.

    [~,DIM] = size(model.Nodes);
    if isfield(model,'param')
        E0 = model.param(1); nu = model.param(2);
    elseif isfield(model,'params')
        E0 = model.params(1); nu = model.params(2);
    else
        error('linearConstitutive:MissingMaterialParams', ...
            'model.param or model.params must be defined.');
    end
    if DIM == 2
        stressType = model.stressType;
        Cmat = zeros(3);
        if strcmp(stressType,'plane_stn')
            C0 = E0/(1+nu)/(1-2*nu);
            Cmat(1,1) = 1 - nu;
            Cmat(1,2) = nu;
            Cmat(2,1) = Cmat(1,2);
            Cmat(2,2) = Cmat(1,1);
            Cmat(3,3) = (1-2*nu)/2.;
            Cmat = C0*Cmat;
        elseif strcmp(stressType,'plane_str')
            C0 = E0/(1-nu*nu);
            Cmat(1,1) = 1.;
            Cmat(1,2) = nu;
            Cmat(2,1) = nu;
            Cmat(2,2) = 1.;
            Cmat(3,3) = (1-nu)/2.;
            Cmat = C0*Cmat;
        end
    elseif DIM == 3    % 3D
        Cmat = zeros(6);
        D0 = E0/(1+nu)/(1-2*nu);
        Cmat(1,1) = (1-nu);
        Cmat(1,2) = nu;
        Cmat(1,3) = Cmat(1,2);
        Cmat(2,1) = Cmat(1,2);
        Cmat(2,2) = Cmat(1,1);
        Cmat(2,3) = Cmat(1,2);
        Cmat(3,1) = Cmat(1,3);
        Cmat(3,2) = Cmat(2,3);
        Cmat(3,3) = Cmat(1,1);
        Cmat(4,4) = 0.5*(1-2*nu);
        Cmat(5,5) = Cmat(4,4);
        Cmat(6,6) = Cmat(4,4);
        Cmat = D0*Cmat;
    end

    model.Cmat = Cmat;
    
end
