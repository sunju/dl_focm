clc;close all; clear all;
% Reproduce the reparameterized function landscape in Fig. 2 and Fig. 4 of the paper
% "Complete Dictionary Recovery over the Sphere",
% by Ju Sun, Qing Qu, and John Wright.
% 
% Dictionary Learning (DL) Problem formulation: Y = A_0*X_0,
% A_0: n-by-n orthogonal matrix, fixed to be A_0 = I;
% X_0: n-by-p matrix, row independent or column sparse
% Reparameterize q = [w ; sqrt(1-||w||^2)]', ||w||<=1,
% Plot the function landscape h_mu = mu*log cosh(q'*Y/mu) in the w space.
%
% 
% Code written by Ju Sun, Qing Qu and John Wright. 
% Last updated: Sat 25 Apr 2015 10:29:09 PM EDT 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rng(1,'twister'); % fix the seed for random number generation

%% parameter settings
tVec = 0:.02:1;
thetaVec = 0:.02:2*pi;
numPts = length(tVec) * length(thetaVec);

mu = 1/50;% smoothing parameters
theta = 0.1;% sparsity
p = 100000;% number of samples

% data model selection:
% 1 = planted sparse vector model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Model 1 is assumed in the paper:   
% Finding a Sparse Vector in a Subspace: Linear Sparsity Using Alternating Direction. 
% By  Qing Qu, Ju Sun, John Wright. NIPS'14. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2 = DL model, independent row, i.i.d. Bernoulli-Gaussian
% 3 = DL model, independent row, i.i.d. Bernoulli-Uniform
% 4 = DL model, correlated row, Bernoulli-Gaussian
% 5 = DL model, correlated row, Bernoulli-Uniform
choice = 2;


%% generate random data
%% B is square root of the covariance matrix 
A = 0.1*randn(3, 3);
A = A - diag(diag(A)) + diag(ones(3, 1)); 
B = 0.5*(A' + A);
% generate data
switch choice
    case 1
        % PSV
%         X = (2*double(rand(3,p) < .5) - 1) .* [randn(2,p); double( rand(1,p) < theta ) / sqrt(theta)];
        X = [randn(2,p); randn(1,p) .* double( rand(1,p) < theta ) / sqrt(theta)];
    case 2
        X = randn(3,p) .* double( rand(3,p) < theta );
    case 3
        X = (rand(3,p)-0.5) .* double( rand(3,p) < theta );
    case 4
        X = (B *randn(3,p)).* double( rand(3,p) < theta );
    case 5
        X = (B *(rand(3,p)-0.5)).* double( rand(3,p) < theta );
end

fnVals = zeros(numPts,1);
uVals = zeros(numPts,1);
vVals = zeros(numPts,1);

%% evaluate function values
k = 1;
for i = 1:length(tVec)
    if mod(i,50) == 1,
        disp(i);
    end
    
    for j = 1:length(thetaVec),
        t = tVec(i);
        theta = thetaVec(j);
        
        u = t*cos(theta);
        v = t*sin(theta);
        
        q = [ u; v; sqrt(1-u^2-v^2) ];
        
        % calculate the gradient of the objective wrt w = [u v]
        z = q'*X;
        
        fnVals(k) = sum( mu * log( cosh(z/mu) ) ) / p;
        
        uVals(k) = u;
        vVals(k) = v;
        
        k = k + 1;
    end
end

%% plot function landscape
figure(1); clf; plot3( uVals, vVals, fnVals, 'b.' );

hold on;

% plot a circle
thetaVec = 0:.001:(2*pi);

cx = zeros(length(thetaVec),1);
cy = zeros(length(thetaVec),1);
cz = zeros(length(thetaVec),1);

for i = 1:length(thetaVec),
    theta = thetaVec(i);
    
    cx(i) = cos(theta);
    cy(i) = sin(theta);
end

cz = min( fnVals ) * ones(size(cz));

[x, y, z] = ellipsoid(0,0,min(real(fnVals)),1,1,0,30);
surf(x, y, z, 25 * ones(size(z)), 'EdgeColor', 'none', 'CDataMapping', 'direct');

plot3(cx,cy,cz,'k-','LineWidth',2);
