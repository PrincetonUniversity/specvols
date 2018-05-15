% SYMMAT_TO_VEC Packs a symmetric matrix into a lower triangular vector
%
% Usage
%    vec = symmat_to_vec(mat);
%
% Input
%    mat: An array of size N-by-N-by-... where the first two dimensions
%       constitute symmetric or Hermitian matrices.
%
% Output
%    vec: A vector of size N*(N+1)/2-by-... consisting of the lower triangular
%       part of each matrix.

% Author
%    Joakim Anden <janden@flatironinstitute.org>

function vec = symmat_to_vec(mat)
    N = size(mat, 1);

    if size(mat, 2) ~= N
        error('matrix must be square');
    end

    [mat, sz_roll] = unroll_dim(mat, 3);

    Jl = tril(true(N*ones(1, 2)));
    Jl = find(Jl(:));

    mat = reshape(mat, [N^2 size(mat, 3)]);

    vec = mat(Jl,:);

    vec = roll_dim(vec, sz_roll);
end
