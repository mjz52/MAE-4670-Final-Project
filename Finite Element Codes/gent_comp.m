% Given the Green-Lagrange strain tensor, E, and material information,
% return the 2nd PK stress tensor. Implementation of compressible Gent 
% constitutive model

function [Se,D,W] = gent_comp(Ee, material, elmID)
    mu = material.mu;
    Jm = material.Jm;
    d = material.d;
    
    Ce = 2*Ee+eye(2); % Right stretch tensor
    C_inv = Ce^-1;
    J = sqrt(det(Ce));
    % A few attempts at the 2nd PK for Gent material ...
%     Se = mu*(Jm/(Jm-trace(Ce)+2)*eye(2) + (d*J^2-(d+1)*J)*C_inv);
%     Se = mu*Jm/(Jm-(trace(Ce)-2))*eye(2) + 1/d*(J^2-1)*C_inv;
    Se = mu*(Jm/(Jm-trace(Ce)+2)*eye(2) - C_inv);
    
    % This technique is apparently slightly faster ...
    C_ = zeros(2,2,2,2);
    for i = 1:2
        for j = 1:2
            for k = 1:2
                for l = 1:2
                    % A few attempts at the constitutive matrix ...
%                     C_(i,j,k,l) = 2*mu*(Jm/(Jm-(trace(Ce)-2))^2 * ((i==j)&&(k==l)&&(i==k)) + (d*J^2-(d+1)/2*J)*C_inv(i,j)*C_inv(k,l) - 1/2*(d*J^2-(d+1)*J)*(C_inv(i,k)*C_inv(l,j) + C_inv(i,l)*C_inv(k,j)));
%                     C_(i,j,k,l) = 2*mu*(Jm/(Jm-(trace(Ce)-2))^2 * ((i==j)&&(k==l)&&(i==k)) + 2/d*J^2*C_inv(i,j)*C_inv(k,l) - 1/d*(J^2-1)*(C_inv(i,k)*C_inv(l,j) + C_inv(i,l)*C_inv(k,j)));
                    C_(i,j,k,l) = 2*mu*(Jm/(Jm-(trace(Ce)-2))^2 * ((i==j)&&(k==l)&&(i==k)) + 1/2*(C_inv(i,k)*C_inv(l,j) + C_inv(i,l)*C_inv(k,j)));
                end
            end
        end
    end
    
    D = [C_(1,1,1,1), C_(1,1,2,2), C_(1,1,1,2);
         C_(2,2,1,1), C_(2,2,2,2), C_(2,2,1,2);
         C_(1,2,1,1), C_(1,2,2,2), C_(1,2,1,2)];
     
%     % C1111
%     D11 = 2*mu*(Jm/(Jm-(trace(Ce)-2))^2 + J*(d*J-(d+1)/2)*C_inv(1,1)*C_inv(1,1) + ...
%                 - 1/2*J*(d*J-(d+1))*(C_inv(1,1)*C_inv(1,1)+C_inv(1,1)*C_inv(1,1)));
%     % C2222
%     D22 = 2*mu*(Jm/(Jm-(trace(Ce)-2))^2 + J*(d*J-(d+1)/2)*C_inv(2,2)*C_inv(2,2) + ...
%                 - 1/2*J*(d*J-(d+1))*(C_inv(2,2)*C_inv(2,2)+C_inv(2,2)*C_inv(2,2)));
%     % C1212
%     D33 = 2*mu*(J*(d*J-(d+1)/2)*C_inv(1,2)*C_inv(1,2) + ...
%                 - 1/2*J*(d*J-(d+1))*(C_inv(1,1)*C_inv(2,2)+C_inv(1,2)*C_inv(1,2)));
%     % C2212
%     D23 = 2*mu*(J*(d*J-(d+1)/2)*C_inv(2,2)*C_inv(1,2) + ...
%                 - 1/2*J*(d*J-(d+1))*(C_inv(2,1)*C_inv(2,2)+C_inv(2,2)*C_inv(1,2)));
%     % C1112
%     D13 = 2*mu*(J*(d*J-(d+1)/2)*C_inv(1,1)*C_inv(1,2) + ...
%                 - 1/2*J*(d*J-(d+1))*(C_inv(1,1)*C_inv(2,1)+C_inv(1,2)*C_inv(1,1)));
%     % C1122
%     D12 = 2*mu*(J*(d*J-(d+1)/2)*C_inv(1,1)*C_inv(2,2) + ...
%                 - 1/2*J*(d*J-(d+1))*(C_inv(1,2)*C_inv(2,1)+C_inv(1,2)*C_inv(2,1)));
%     D = [D11, D12, D13; D12, D22, D23; D13, D23, D33];

    W = 0; % fix
    
end