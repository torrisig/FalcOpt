function [jac_x_structure,jac_u_structure,jac_n_structure] = jacobian_structure()

% structure of jacobian_x matrix
jac_x_structure = [ 1 0 ;...
                    2 3 ];
                
% structure of jacobian_u matrix
jac_u_structure = [ 1 0  ; ...
                    0 2 ];

% structure of  jacobian of nonlinear constraints w.r.t u
jac_n_structure = [ 1   ;...
                    2   ];

end