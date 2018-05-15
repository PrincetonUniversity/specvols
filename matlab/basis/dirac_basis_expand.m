% DIRAC_BASIS_EXPAND Expand array in Dirac basis
%
% Usage
%    v = dirac_basis_expand(basis, x);
%
% Input/Output
%    See documentation for `basis_expand`.
%
% See also
%    basis_expand

% Author
%    Joakim Anden <janden@flatironinstitute.org>

function v = dirac_basis_expand(basis, x)
    [x, sz_roll] = unroll_dim(x, numel(basis.sz)+1);

    x = reshape(x, [prod(basis.sz) size(x, numel(basis.sz)+1)]);

    v = x(basis.mask,:);

    v = roll_dim(v, sz_roll);
end
