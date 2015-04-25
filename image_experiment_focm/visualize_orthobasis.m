function visualize_orthobasis( A )

m = size(A,1);

sp_handle = tight_subplot(sqrt(m), sqrt(m), [0.01, 0], [0.01, 0.01], [0.01, 0.01]); 

for i = 1:m,
	axes(sp_handle(i)); 
    imagesc(reshape(A(:,i),sqrt(m),sqrt(m)));
    colormap('gray');
    axis off;
    axis image;
end


