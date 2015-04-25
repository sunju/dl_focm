% evaluate the value f, gradient g and hessian H 
% for the function h_mu(q'*Y) 
function [f,g,H] = l1_exp_approx(Y,q,mu,flag)
f=[];g=[];H =[];
z = q'*Y;
t = z/mu;
ind = abs(t)>50;%for the purpose of preventing overflow
n = size(Y,1);

f_each_1 = mu * log( cosh(t(~ind)) );
f_each_2 = abs(z(ind)) - mu*log(2);
f = sum(f_each_1)+sum(f_each_2); %function value

if(flag == true)    
    g_each =  tanh(t);
    H_each = (1/mu) * (1-g_each.^2);    
    g = Y * g_each';% gradient 
    H = (Y.* repmat(H_each, n, 1))*Y'; %hessian H = Y * diag(H_each) * Y
end

end
