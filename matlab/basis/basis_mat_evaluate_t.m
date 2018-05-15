% BASIS_MAT_EVALUATE_T Evaluate coefficient matrix in dual basis
%
% Usage
%    V = basis_mat_evaluate_t(basis, X);
%
% Input
%    basis: A basis object (for example, obtained from `dirac_basis`) with
%       respect to whose dual the coefficients are defined.
%    X: The coefficient array of size `basis.sz`-by-`basis.sz` to be
%       evaluated.
%
% Output
%    V: The evaluation of `X` in the dual basis of `basis`. This is a matrix
%       of size `basis.count`-by-`basis.count`.
%
% Description
%    If `V` is a matrix of size `basis.ct`-by-`basis.ct`, `B` is the
%    change-of-basis matrix of `basis`, and `x` is a multidimensional matrix
%    of size `basis.sz`-by-`basis.sz`, the function calculates
%
%       V = B' * X * B
%
%    where the rows of `B`, rows of 'X', and columns of `X` are read as
%    vectorized arrays.
%
% See also
%    basis_evaluate_t, basis_expand_t

% Author
%    Joakim Anden <janden@flatironinstitute.org>

function V = basis_mat_evaluate_t(basis, X)
    fun = @(x)(basis_evaluate_t(basis, x));

    V = mdim_mat_fun_conj(X, numel(basis.sz), 1, fun);
end
