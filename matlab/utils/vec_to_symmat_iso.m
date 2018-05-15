% VEC_TO_SYMMAT_ISO Isometrically map packed vector to symmetric matrix
%
% Usage
%    mat = vec_to_symmat_iso(vec);
%
% Input
%    vec: A vector of size N*(N+1)/2-by-... describing a symmetric (or
%       Hermitian) matrix.
%
% Output
%    mat: An array of size N-by-N-by-... which indexes symmetric/Hermitian
%       matrices that occupy the first two dimensions. The lower triangular
%       parts of these matrices consists of the corresponding vectors in vec,
%       reweighted so that the Euclidean inner product maps to the Frobenius
%       inner product.
%
% See also
%    vec_to_symmat

% Author
%    Joakim Anden <janden@flatironinstitute.org>

function mat = vec_to_symmat_iso(vec)
    mat = vec_to_symmat(vec);

    [mat, sz_roll] = unroll_dim(mat, 3);

    N = size(mat, 1);

    mat = mat_to_vec(mat);

    mat(1:N+1:N^2,:) = sqrt(2)*mat(1:N+1:N^2,:);
    mat = 1/sqrt(2)*mat;

    mat = vec_to_mat(mat);

    mat = roll_dim(mat, sz_roll);
end
