function X = prox_L1(X,d)

X = sign(X) .* max( abs(X) - d, 0 ); 


