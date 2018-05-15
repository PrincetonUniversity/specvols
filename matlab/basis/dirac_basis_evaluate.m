% DIRAC_BASIS_EVALUATE Evaluate coefficient vector in Dirac basis
%
% Usage
%    x = dirac_basis_evaluate(basis, v);
%
% Input/Output
%    See documentation for `basis_evaluate`.
%
% See also
%    basis_evaluate

% Author
%    Joakim Anden <janden@flatironinstitute.org>

function x = dirac_basis_evaluate(basis, v)
    [v, sz_roll] = unroll_dim(v, 2);

    x = zeros([prod(basis.sz) size(v, 2)], class(v));

    x(basis.mask,:) = v;

    x = reshape(x, [basis.sz size(x, 2)]);

    x = roll_dim(x, sz_roll);
end
