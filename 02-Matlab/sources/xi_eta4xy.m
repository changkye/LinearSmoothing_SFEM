function [xi,eta] = xi_eta4xy(elemType,side,N_g)

    if strcmp(elemType,'Q4')
        if side == 1
            xi  = [-1 1]*N_g;
            eta = -1;
        elseif side == 2
            xi = 1;
            eta = [-1 1]*N_g;
        elseif side == 3
            xi = [-1 1]*N_g;
            eta = 1;
        else %side==4
            eta = [-1 1]*N_g;
            xi = -1;
        end 
    elseif strcmp(elemType,'T3')
        if side == 1
            xi = [0 1]*N_g;
            eta = 0;
        elseif side == 2
%         xi = 0.5 + 0.5*[-1 1]*N_g;
%         eta= 0.5 - 0.5*[-1 1]*N_g;
            xi = [0 1]*N_g;
            eta = [1 0]*N_g;
        else%if side==3
            eta = [0 1]*N_g;
            xi = 0;
        end 
    elseif strcmp(elemType,'T6')
        if side == 1
            xi = [0 1 0.5]*N_g;
            eta = 0;
        elseif side == 2
            xi = [0 1 0.5]*N_g;
            eta = [1 0 0.5]*N_g;
        else%if side==3
            eta = [0 1 0.5]*N_g;
            xi = 0;
        end 
    end 

end