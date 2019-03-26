function [laplacian_eigvecs, laplacian_basis_coefs] = draw_pca_vs_laplacian_blob1d(noise)
L = 100;
N_BASIS_VECTORS = 10;
BLOB_WIDTH = 4;

data = gen_moving_blob_1d(L, BLOB_WIDTH, noise);
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
pca_reconstruction = pca_basis*pca_basis_coefs; 
imagesc(pca_reconstruction);
title('PCA reconstruction');
fprintf('PCA reconstruction error: %f', norm(pca_reconstruction-data,'fro'));

W = graph_gaussian_kernel(data, 0.1);
laplacian_matrix = laplacian(W, 'normalized');
[laplacian_eigvecs,D] = eigs(laplacian_matrix,N_BASIS_VECTORS+1,'smallestabs');
assert(all(size(laplacian_eigvecs) == [L N_BASIS_VECTORS+1]));
laplacian_eigvecs = laplacian_eigvecs(:,2:end); % Get rid of DC
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
laplacian_reconstruction = laplacian_eigvecs*laplacian_basis_coefs;
imagesc(laplacian_reconstruction);
title('Laplacian reconstruction');
fprintf('Laplacian reconstruction error: %f', norm(laplacian_reconstruction-data,'fro'));

end