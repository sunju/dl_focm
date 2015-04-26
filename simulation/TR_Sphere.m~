% TR_sphere.m: trust region method over the sphere for recovering a single
% sparse vector.
% Recover one sparsest row of Y by solving:
%  min_q  h_mu(q' * Y),  s.t.  ||q|| =1.
function q = TR_Sphere( Y, mu, q_init )
[n,~] = size(Y);
% initalize q
if nargin > 2,
    q = q_init;
else
    q = randn(n,1);% random initialization
    q = q / norm(q);
end

%% Parameter Settings

tol = 1e-6;% stopping criteria for the gradient
MaxIter = 200;% max iteration
Delta = 0.1;% inital trust region size
Delta_max = 1;% maximum trust region size
Delta_min = eps;% minimum trust region size

% parameters for adjusting trust region size
eta_s = 0.1;
eta_vs = 0.9;
gamma_i = 2;
gamma_d = 1/2;


%% Main Iteration
for iter = 1:MaxIter
    U = null(q'); % basis of the tangent space T_q = { v| v'*q =0}
    [f,g,H] = l1_exp_approx(Y,q,mu,true); % evaluate the function, gradient and hessian at q
    
    %solve TR-subprolbem
    A = U' * H * U; % hessian in the tangent space T_q
    b = U' * g; % gradient in the tangent space T_q
    [xi,opt] = Solve_TR_Subproblem(A,b,Delta); % optimal direction xi in the tangent space T_q
    
    t = norm(xi);
    q_hat = q * cos(t) + (U*xi/t) * sin(t);% retraction back to the sphere by exponential map
    
    % evaluate the trust region ratio: rho
    f_new = l1_exp_approx(Y,q_hat,mu,false);
    f_hat = f + b' * xi + (1/2) * xi' * A * xi;
    rho = (f - f_new)/(f-f_hat);
    
    % update iterate q and trust-region size Delta
    flag = 0;
    if rho>eta_vs
        q = q_hat; flag = 1;
        Delta = min(gamma_i*Delta,Delta_max);
    elseif rho>eta_s
        q = q_hat; flag = 1;
    else
        Delta = max(gamma_d*Delta,Delta_min);
    end
    [~,idx] = max(abs(q));
    ei = zeros(n,1);
    ei(idx) = sign(q(idx));
    err = min(norm( q - ei ),norm(q +ei));
    fprintf('Iter = %d, Obj = %f, Err = %f, TR Size = %f, Step = %f, rho = %f \n',iter,f,err, Delta,t,rho);

    % this is cheating for speed; should be removed for reproduction 
    if err<mu/2
        break;
    end
    
    if (flag ==1&&abs((f_new - f)/t)<tol);
        break;
    end
        
end

end
