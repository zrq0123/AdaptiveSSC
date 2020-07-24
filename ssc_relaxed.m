function [ Z ] = ssc_relaxed( X, lambda,W)
%% Solves the following
%
% min || X - XZ ||_F^2 + lambda || Z ||_1
%
% via ADM
%
% Created by Stephen Tierney
% stierney@csu.edu.au
%

max_iterations = 500;

func_vals = zeros(max_iterations, 1);

n = size(X, 2);

Z = zeros(n);
J = zeros(n);
Y = zeros(n);

mu = 1;

tol_1 = 1*10^-2;
tol_2 = 1*10^-4;

lambda=lambda./abs(W);

lambda=lambda-diag(diag(lambda));

for k = 1 : max_iterations
    
    Z_prev = Z;
    J_prev = J;
    
    % Solve for J
    
    J = (X'*X + mu*speye(size(J))) \ (X'*X + mu*Z-mu*Y);
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
    
end
disp(k)
end

