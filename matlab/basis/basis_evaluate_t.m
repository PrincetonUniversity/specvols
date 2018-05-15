% BASIS_EVALUATE_T Evaluate coefficient in dual basis
%
% Usage
%    v = basis_evaluate_t(basis, x);
%
% Input
%    basis: A basis object (for example, obtained from `dirac_basis`) with
%       respect to whose dual the coefficients are defined.
%    x: The coefficient array to be evaluated. The first dimensions must
%       equal `basis.sz`.
%
% Output
%    v: The evaluation of the coefficient array `x` in the dual basis of
%       `basis`. This is an array of vectors whose first dimension equals
%       `basis.count` and whose remaining dimensions correspond to higher
%       dimensions of `x`.
%
% Description
%    If `v` is a matrix of size `basis.ct`-by-..., `B` is the change-of-basis
%    matrix of `basis`, and `x` is a matrix of size `basis.sz`-by-..., the
%    function calculates
%
%       v = B' * x
%
%    where the rows of `B` and columns of `x` are read as vectorized arrays.
%
% See also
%    basis_expand_t

% Author
%    Joakim Anden <janden@flatironinstitute.org>

function v = basis_evaluate_t(basis, x)
    sz_x = size(x);

    if numel(sz_x) < numel(basis.sz) || ...
        ~all(sz_x(1:numel(basis.sz)) == basis.sz)
        error('First dimensions of `x` must match `basis.sz`.');
    end

    if basis.type == basis_type_dirac()
        v = dirac_basis_expand(basis, x);
    elseif basis.type == basis_type_fb()
        v = fb_basis_evaluate_t(basis, x);
    else
        error('Invalid basis type.');
    end
end
