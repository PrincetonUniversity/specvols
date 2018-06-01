% MATRIX_BASIS_EXPAND_T Expand array in dual matrix basis
%
% Usage
%    x = matrix_basis_expand_t(basis, v);
%
% Input/Output
%    See documentation for `basis_expand_t`.
%
% See also
%    basis_expand_t

% Author
%    Joakim Anden <janden@flatironinstitute.org>

function x = matrix_basis_expand_t(basis, v)
    [v, sz_roll] = unroll_dim(v, 2);

    x = basis.B'\v;

    x = reshape(x, [basis.sz size(x, 2)]);

    x = roll_dim(x, sz_roll);
end
