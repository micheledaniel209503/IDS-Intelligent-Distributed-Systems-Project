function [x_opt,Qlist,Qprod] = RoundRobinQmatrix(n)
% ROUNDROBINQMATRIX Constructs Qn..Q1 matrices of size nxn with
% double stochasticity constraints and finds a non-negative solution
% with minimal norm.
%
% INPUT:
%   n = number of matrices <--> dimension
% OUTPUT:
%   x_opt = solution vector [x1..x2n]
%   Qlist = cell array containing matrices Qk
%   Qprod = product Qn*Qn-1*...*Q1

    tol = 0.05; % lower bound for the unknown

    % Number of variables: 2 per each Qk
    nVars = 2*n;

    % cost function
    fun = @(x) norm(x-0.5);

    % Nonlinear constraints (double stochasticity)
    nonlcon = @(x) stochastic_constraints_n(x,n);

    %% Bounds: 
    lb = zeros(nVars,1) + tol;
    ub = ones(nVars,1);

    % Initial point: 0.5
    x0 = 0.5*ones(nVars,1);

    % Solve with fmincon
    options = optimoptions('fmincon','Display','iter','Algorithm','sqp');
    [x_opt,~,exitflag] = fmincon(fun,x0,[],[],[],[],lb,ub,nonlcon,options);

    if exitflag <= 0
        warning('Optimization did not converge!');
    end

    % Reconstruct the Qk matrices
    Qlist = build_Q_matrices_n(x_opt,n);

    % Product Q = Qn*...*Q1
    Qprod = eye(n);
    for k = 1:n
        Qprod =  Qlist{k} * Qprod;
    end

end

% --- Function to construct Qk matrices from numeric vector x ---
function Qlist = build_Q_matrices_n(x,n)
    Qlist = cell(1,n);
    for k = 1:n
        Qk = zeros(n);
        Qk(k,k) = 1; % set the k-th row diagonal to 1
        for i = 1:n
            if i ~= k
                Qk(i,k) = x(2*k-1); % k-th column
                Qk(i,i) = x(2*k);   % diagonal
            end
        end
        Qlist{k} = Qk;
    end
end

% --- Double stochasticity constraints for vector x ---
function [c,ceq] = stochastic_constraints_n(x,n)
    % Construct numeric Qk matrices
    Qlist = build_Q_matrices_n(x,n);
    Qprod = eye(n);
    for k = 1:n
        Qprod = Qlist{k} * Qprod;
    end

    % Double stochasticity constraints
    ceq = [];
    for i = 1:n
        ceq = [ceq; sum(Qprod(i,:)) - 1]; % row sums = 1
    end
    for j = 1:n
        ceq = [ceq; sum(Qprod(:,j)) - 1]; % column sums = 1
    end

    c = []; 
end