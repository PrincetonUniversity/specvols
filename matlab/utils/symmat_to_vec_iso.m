% SYMMAT_TO_VEC_ISO Isometrically maps a symmetric matrix to a packed vector
%
% Usage
%    vec = symmat_to_vec_iso(mat);
%
% Input
%    mat: An array of size N-by-N-by-... where the first two dimensions
%       constitute symmetric or Hermitian matrices.
%
% Output
%    vec: A vector of size N*(N+1)/2-by-... consisting of the lower triangular
%       part of each matrix, reweighted so that the Frobenius inner product
%       is mapped to the Euclidean inner product.
%
% See also
%    symmat_to_vec

% Author
%    Joakim Anden <janden@flatironinstitute.org>

function vec = symmat_to_vec_iso(mat)
    [mat, sz_roll] = unroll_dim(mat, 3);

    N = size(mat, 1);

    mat = mat_to_vec(mat);

    mat(1:N+1:N^2,:) = 1/sqrt(2)*mat(1:N+1:N^2,:);
    mat = sqrt(2)*mat;

    mat = vec_to_mat(mat);

    mat = roll_dim(mat, sz_roll);

    vec = symmat_to_vec(mat);
end
