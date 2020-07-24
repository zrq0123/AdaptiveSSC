function [ Z ] = ssc_relaxed_adaptive( X, lambda,W)
%% Solves the following
%
% min || X - XZ ||_F^2 + lambda || Z ||_1
%
% via ADM
%
% Created by Stephen Tierney
% stierney@csu.edu.au
%

max_iterations = 200;

func_vals = zeros(max_iterations, 1);

n = size(X, 2);

Z = zeros(n);
J = zeros(n);
Y = zeros(n);

mu = 1;

tol_1 = 1*10^-2;
tol_2 = 1*10^-4;
u=50;

lambda=lambda./abs(W);

lambda=lambda-diag(diag(lambda));

tmp1=X'*X + mu*speye(size(J));
tmp2=X'*X;

for k = 1 : max_iterations
    
    Z_prev = Z;
    J_prev = J;
    
    % Solve for J
    
    J = (tmp1) \ (tmp2 + mu*Z-mu*Y);
    J = J - diag(diag(J));
    J(J<0)=0;
    % Solve for Z
    
    V = J + Y;
    
    Z = solve_l1(V, lambda ./ mu );
    
    Z = Z - diag(diag(Z)); % make the diaginal element to zero
    Z(Z<0)=0;
    % Update Y
    
    Y = Y + mu*(J-Z);
    
    % Check convergence
    
    %func_vals(k) = 0.5*norm(X - X*Z, 'fro')^2 + lambda*norm_l1(Z);
    
    if (errorCoef(J,Z)<tol_2)
        break;
    end
    if (mod(k,10)==0)
    %vary the Penalty Parameter mu
        primal_r = norm(Z-J,'fro');
        dual_r = norm(mu*(Z-Z_prev),'fro');
        if(primal_r>u*dual_r)
            mu=mu*2;
        elseif(dual_r>u*primal_r)
            mu=mu/2;
        end
        %disp(mu);
    end
end
sprintf('%.5f\t%.5f',primal_r,dual_r)
disp(k)
end

