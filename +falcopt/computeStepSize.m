%% computeStepSize
% this function works only with CasADi and yalmip
%
%
% Copyright (c) 2018, ETH Zurich, Automatic Control Laboratory 
%                    Damian Frick <falcopt@damianfrick.com>
%                    Giampaolo Torrisi <giampaolo.torrisi@gmail.com>
%                    Tommaso Robbiani <tommasro@student.ethz.ch>
%
% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, including without limitation the rights
% to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
% copies of the Software, and to permit persons to whom the Software is
% furnished to do so, subject to the following conditions:
%
% The above copyright notice and this permission notice shall be included in all
% copies or substantial portions of the Software.
%
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
% IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
% AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
% OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
% SOFTWARE.
%

function alpha_opt = computeStepSize(x0, u_ref, o)
    % TODO: allow to also give w_ref!

    % check dimensions inputs
    if any(size(x0) ~= [o.nx, 1]) || any(size(u_ref) ~= [o.nu, o.N])
        throw(MException(['computeStepSize:InvalidParameter'], 'Invalid parameter. Note: computeStepSize is designed only for internal use.'));
    end

    if o.nw > 0
        w_ref = zeros(o.nw,1);
    end

    x = casadi.SX.sym('x',o.nx,1);
    u = casadi.SX.sym('u',o.nu,1);
    u_n = casadi.SX.sym('u_n',o.nu,o.N);
    psi = casadi.SX.sym('psi',o.nx,o.N);
    J = casadi.SX.sym('J',1);

    % define function
    if o.nw > 0
        w = casadi.SX.sym('w',o.nw);
        dynamics = casadi.Function('y_fun',{x,u,w},{o.dynamics(x,u,w)});
    else
        dynamics = casadi.Function('y_fun',{x,u},{o.dynamics(x,u)});
    end

    % define psi as function of x0 and u
    if o.nw > 0
        for k=1:o.N
            if k==1
                psi(:,k) = dynamics(x,u_n(:,k),w);
            else
                psi(:,k) = dynamics(psi(:,k-1),u_n(:,k),w);
            end
        end
    else
        for k=1:o.N
            if k==1
                psi(:,k) = dynamics(x,u_n(:,k));
            else
                psi(:,k) = dynamics(psi(:,k-1),u_n(:,k));
            end
        end
    end

    % define cost J
    if isfield(o.objective,'nonlinear')
        for k = 1:o.N-1
            J = J + o.objective.nonlinear(psi(:,k),u_n(:,k));
        end
        if isfield(o.objective,'nonlinearN')
            J = J + o.objective.nonlinearN(psi(:,end));
        end
    else
        for k = 1:o.N-1
            J = J + 0.5*(psi(:,k)'*o.objective.Q*psi(:,k) + u_n(:,k)'*o.objective.R*u_n(:,k));
        end
        J = J + 0.5*(psi(:,end)'*o.objective.Q*psi(:,end));
    end

    % compute Hessian
    try 
        HJ = hessian(J,u_n);
    catch
        error(['error while computing variable_stepSize.alpha_max, consider to set '...
               'manualy a value for it']);
    end

    % find alpha max
    if o.nw > 0
        HJ_fun = casadi.Function('HJ_fun',{x,u_n,w},{HJ});
        alpha_opt = 1/max(eig(full(HJ_fun(x0,u_ref,w_ref))));
    else
        HJ_fun = casadi.Function('HJ_fun',{x,u_n},{HJ});
        alpha_opt = 1/max(eig(full(HJ_fun(x0,u_ref))));
    end
    
end