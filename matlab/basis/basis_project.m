% BASIS_PROJECT Project array into basis
%
% Usage
%    y = basis_project(basis, x);
%
% Input
%    basis: A basis object onto whose span `x` is to be projected.
%    x: An array of size `basis.sz`.
%
% Output
%    y: The array `x` expanded in `basis` (using `basis_expand`), and then
%       evaluated (using `basis_evaluate`).
%
% See also
%    basis_expand, basis_evaluate

% Author
%    Joakim Anden <janden@flatironinstitute.org>

function y = basis_project(basis, x)
    y = basis_evaluate(basis, basis_expand(basis, x));
end
