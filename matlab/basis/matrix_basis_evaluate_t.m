% MATRIX_BASIS_EVALUATE_T Evaluate coefficient in dual matrix basis
%
% Usage
%    v = matrix_basis_evaluate_t(basis, x);
%
% Input/Output
%    See documentation for `basis_evaluate_t`.
%
% See also
%    basis_evaluate_t

% Author
%    Joakim Anden <janden@flatironinstitute.org>

function v = matrix_basis_evaluate_t(basis, x)
    [x, sz_roll] = unroll_dim(x, numel(basis.sz)+1);

    x = reshape(x, [prod(basis.sz) size(x, numel(basis.sz)+1)]);

    v = basis.B'*x;

    v = roll_dim(v, sz_roll);
end
