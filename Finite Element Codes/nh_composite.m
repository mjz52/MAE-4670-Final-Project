% Given the Green-Lagrange strain tensor, E, and material information,
% return the 2nd PK stress tensor. Implementation of compressible 
% NeoHookean constitutive model. Uses a composite material description

function [Se,D,W] = nh_compos(Ee, material, elmID)
    
    % Get index for fiber (1) or matrix (2)
    comp_mat = material.matList(elmID); % check whether fiber or matrix
    % Obtain the material parameters for the fiber or matrix
    mu = material.mu(comp_mat);
    beta = material.beta(comp_mat);
    
    % 2nd Pk
    Ce = 2*Ee+eye(2); % Right stretch tensor
    C_inv = Ce^-1;
    J = sqrt(det(Ce));
    Se = mu*(eye(2)-C_inv) + beta*J*(J-1)*C_inv;
    
    % C1111
    D11 = (mu - beta*J*(J-1))*(C_inv(1,1)*C_inv(1,1) + C_inv(1,1)*C_inv(1,1)) + ...
          beta*J*(2*J-1)*C_inv(1,1)*C_inv(1,1);
    % C2222
    D22 = (mu - beta*J*(J-1))*(C_inv(2,2)*C_inv(2,2) + C_inv(2,2)*C_inv(2,2)) + ...
          beta*J*(2*J-1)*C_inv(2,2)*C_inv(2,2);
    % C1212
    D33 = (mu - beta*J*(J-1))*(C_inv(1,1)*C_inv(2,2) + C_inv(1,2)*C_inv(1,2)) + ...
          beta*J*(2*J-1)*C_inv(1,2)*C_inv(1,2);
    % C2212
    D23 = (mu - beta*J*(J-1))*(C_inv(2,1)*C_inv(2,2) + C_inv(2,2)*C_inv(1,2)) + ...
          beta*J*(2*J-1)*C_inv(2,2)*C_inv(1,2);
    % C1112
    D13 = (mu - beta*J*(J-1))*(C_inv(1,1)*C_inv(2,1) + C_inv(1,2)*C_inv(1,1)) + ...
          beta*J*(2*J-1)*C_inv(1,1)*C_inv(1,2);
    % C1122
    D12 = (mu - beta*J*(J-1))*(C_inv(1,2)*C_inv(2,1) + C_inv(1,2)*C_inv(2,1)) + ...
          beta*J*(2*J-1)*C_inv(1,1)*C_inv(2,2);
    D = [D11, D12, D13; D12, D22, D23; D13, D23, D33];
    
    % Strain energy density
    W = mu/2*(trace(Ce)-2-2*log(J)) + beta/2*(J-1)^2;

end