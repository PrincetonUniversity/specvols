% BASIS_MAT_EXPAND Expand array matrix in basis
%
% Usage
%    V = basis_mat_expand(basis, X);
%
% Input
%    basis: A basis object (for example, obtained from `dirac_basis`) in which
%       the array is to be expanded.
%    X: An multidimensional matrix of size `basis.sz`-by-`basis.sz` to be
%       expanded in.
%
% Output
%    V: The coefficients of `X` in `basis` in the form of a matrix of size
%       `basis.ct`-by-`basis.ct`.
%
% Description
%    If `V` is a matrix of size `basis.ct`-by-`basis.ct`, `B` is the
%    change-of-basis matrix of `basis`, and `x` is a multidimensional matrix
%    of size `basis.sz`-by-`basis.sz`, the function calculates
%
%       V = (B' * B)^(-1) * B' * X * B * (B' * B)^(-1)
%
%    where the rows of `B`, rows of 'X', and columns of `X` are read as
%    vectorized arrays.
%
% See also
%    basis_expand, basis_mat_evaluate

% Author
%    Joakim Anden <janden@flatironinstitute.org>

function V = basis_mat_expand(basis, X)
    fun = @(x)(basis_expand(basis, x));

    V = mdim_mat_fun_conj(X, numel(basis.sz), 1, fun);
end
