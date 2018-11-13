function [laplacian_eigvecs, laplacian_basis_coefs] = draw_pca_vs_laplacian_blob1d()
L = 100;
N_BASIS_VECTORS = 25;
BLOB_WIDTH = 4;

data = gen_moving_blob_1d(L, BLOB_WIDTH);
figure;
subplot(3,3,1)

imagesc(data);
title('Input data (n by p)');


subplot(3,3,4);
[coeff,score,latent,tsquared,explained,mu] = pca(data);
assert(size(coeff,1) == L);
assert(size(coeff,2) >= N_BASIS_VECTORS);

fprintf('Percentage of variance explained using %d principal components: %f\n', N_BASIS_VECTORS, sum(explained(1:N_BASIS_VECTORS)));
pca_basis = coeff(:,1:N_BASIS_VECTORS); 
assert(all(size(pca_basis) == [L N_BASIS_VECTORS]));
imagesc(pca_basis');
title('PCA')
caxis([-0.2 0.2]);
colorbar;

subplot(3,3,5);
pca_basis_coefs = pca_basis \ data;
%assert(all(size(pca_basis_coefs) == [N_BASIS_VECTORS L]));
imagesc(pca_basis_coefs);
title('PCA coefficients');
caxis([-0.2 0.2]);
colorbar;

subplot(3,3,6)
imagesc(pca_basis*pca_basis_coefs);
title('PCA reconstruction');

W = graph_knn(data, 5);
laplacian_matrix = laplacian(W, 'normalized');
[laplacian_eigvecs,D] = eigs(laplacian_matrix,N_BASIS_VECTORS,'smallestabs'); 
assert(all(size(laplacian_eigvecs) == [L N_BASIS_VECTORS]));
subplot(3,3,7);
imagesc(laplacian_eigvecs')
title('Laplacian Eigenvectors')
caxis([-0.2 0.2]);
colorbar;

laplacian_basis_coefs = laplacian_eigvecs \ (data');
assert(all(size(laplacian_basis_coefs) == [N_BASIS_VECTORS L]));
subplot(3,3,8);
imagesc(laplacian_basis_coefs);
title('Diffusion arrays');
caxis([-0.2 0.2]);
colorbar;

subplot(3,3,9)
imagesc(laplacian_eigvecs*laplacian_basis_coefs);
title('Laplacian reconstruction');
end

