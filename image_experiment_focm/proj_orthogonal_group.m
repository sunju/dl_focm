function Q = proj_orthogonal( M )

% project a square matrix M onto the orthgonal group, in the Frobenius
% sense

[U,S,V] = svd(M);

Q = U*V';
