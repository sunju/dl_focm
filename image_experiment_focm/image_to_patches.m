function Y = image_to_patches(I,w)

[m1,m2] = size(I);

n1 = floor(m1/w);
n2 = floor(m2/w);

Y = zeros(w*w,n1*n2);

k = 1;

for i = 1:n1,
    for j = 1:n2,
        Y(:,k) = vec( I( ((i-1)*w+1):(i*w), ((j-1)*w+1):(j*w) ) ); 
        k = k+1;
    end
end