clc;close all;clear all;
% Implementation for the Phase Transition of the Trust Region (TR) Method presented in the paper
% "Complete Dictionary Recovery over the Sphere",
% by Ju Sun, Qing Qu, and John Wright. 
% http://arxiv.org/abs/1504.06785
% 
% 1. l1_exp_approx.m: evaluate the value, gradient and hessian 
% for the function: h_mu(x) = mu * log cosh(x/mu)
%
% 2. TR_sphere.m: trust region method over the sphere for recovering a single
% sparse vector
%
% 3. Solve_TR_Subproblem.m: solving the trust region subproblem in
% TR_sphere.m using cvx toolbox, MUST BE INSTALLED ahead!
% It can be downloaded from http://cvxr.com/cvx/download/. 
%
% Problem formulation: Y = A_0*X_0, 
% A_0: n-by-n orthogonal matrix, fixed to  be A_0 = I;
% X_0: n-by-p matrix, each column with fixed k support and standard Gaussian
% distribution for nonzero entries
% We recover one of the sparsest rows in Y = X_0 by solving
%  min_q  h_mu(q' * Y),  s.t.  ||q|| =1.
%
% 
% Code written by Ju Sun, Qing Qu and John Wright. 
% Last Updated: Tue 02 Jun 2015 02:52:57 PM EDT 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Environment Settings
n = [2:5:150]; % dictionary size n
mu = 1/100;% smoothing parameter mu
p = ceil(5*n.^2 *log(n)); % number of samples p
k = [1:5:130]; % sparsity level
num_rep = 5; % number of simulation repetitions
tol = mu; % error tolerance  
Prob = zeros(length(k),length(n)); % record the probability of success

rng(1,'twister'); % fix the seed for random number generation
%% Main Simulation
for t = 1:num_rep
    for i = 1:length(n)
        for j = 1:length(k)
            if k(j)>n(i)
                continue;
            end
            % generate the data 
            Y = zeros(n(i),p(i));
            for ll = 1:p(i)
                % generate p columns of k-sparse vectors for Y
                Y(randperm(n(i),k(j)),ll) = randn(k(j),1);
            end
            
            % solve the TR-problem for finding one sparse vector
            q = TR_Sphere(Y,mu);
            [~,idx] = max(abs(q));
            ei = zeros(n(i),1);
            ei(idx) = sign(q(idx));
            err = min(norm( q - ei ),norm(q +ei));
            if err<=tol %judging correctness
                Prob(j,i) = Prob(j,i) + 1;
            end
            
            % print intermediate results
            fprintf('Simulation=%d, Dimension = %d, Sparsity = %d\n'...
                ,t,n(i),k(j));
            figure(1);
            imagesc(n, k, Prob/t);
            set(gca,'YDir','normal');
            xlabel('Dictionary Dimension n');
            ylabel('Sparsity Level k');
            title('$$p = 5n^2 \log(n)$$', 'interpreter','latex');
            colormap('gray');
            colorbar;
            pause(.25);
        end
    end
    
end

%% Plot Result
figure(1);
imagesc(n, k, Prob/num_rep);
set(gca,'YDir','normal');
xlabel('Dictionary Dimension n');
ylabel('Sparsity Level k');
title('$$p = 5n^2 \log(n)$$', 'interpreter','latex');
colorbar;
colormap('gray');
