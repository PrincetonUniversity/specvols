% BASIS_MAT_EVALUATE Evaluate coefficient matrix in basis
%
% Usage
%    X = basis_mat_evaluate(basis, V);
%
% Input
%    basis: A basis object (for example, obtained from `dirac_basis`) with
%       respect to which the coefficients are defined.
%    V: A coefficient matrix of size `basis.count`-by-`basis.count` to be
%       evaluated.
%
% Output
%    X: A multidimensional matrix of size `basis.sz`-by-`basis.sz`
%       corresponding to the evaluation of `V` in the basis.
%
% Description
%    If `V` is a matrix of size `basis.ct`-by-`basis.ct`, `B` is the
%    change-of-basis matrix of `basis`, and `x` is a multidimensional matrix
%    of size `basis.sz`-by-`basis.sz`, the function calculates
%
%       X = B * V * B'
%
%    where the rows of `B`, rows of 'X', and columns of `X` are read as
%    vectorized arrays.
%
% See also
%    basis_evaluate, basis_mat_expand

% Author
%    Joakim Anden <janden@flatironinstitute.org>

function X = basis_mat_evaluate(basis, V)
    fun = @(v)(basis_evaluate(basis, v));

    X = mdim_mat_fun_conj(V, 1, numel(basis.sz), fun);
end
