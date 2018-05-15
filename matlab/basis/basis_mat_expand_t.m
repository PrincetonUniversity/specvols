% BASIS_MAT_EXPAND_T Expand array matrix in dual basis
%
% Usage
%    X = basis_mat_expand_t(basis, V);
%
% Input
%    basis: A basis object (for example, obtained from `dirac_basis`) in whose
%       dual basis the array is to be expanded.
%    V: A matrix of size `V` to be expanded in the dual basis.
%
% Output
%    X: The coefficients of `V` in the dual basis of `basis` in the form of
%       multidimensional matrix of size `basis.sz`-by-`basis.sz`.
%
% Description
%    If `V` is a matrix of size `basis.ct`-by-`basis.ct`, `B` is the
%    change-of-basis matrix of `basis`, and `x` is a multidimensional matrix
%    of size `basis.sz`-by-`basis.sz`, the function calculates
%
%       X = (B * B')^(-1) * B * V * B' * (B * B')^(-1)
%
%    where the rows of `B`, rows of 'X', and columns of `X` are read as
%    vectorized arrays.
%
% See also
%    basis_expand_t, basis_mat_evaluate_t

% Author
%    Joakim Anden <janden@flatironinstitute.org>

function X = basis_mat_expand_t(basis, V)
    fun = @(v)(basis_expand_t(basis, v));

    X = mdim_mat_fun_conj(V, 1, numel(basis.sz), fun);
end
