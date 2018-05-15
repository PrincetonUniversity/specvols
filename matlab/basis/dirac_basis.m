% DIRAC_BASIS Construct a Dirac basis object
%
% Usage
%    basis = dirac_basis(sz, mask);
%
% Input
%    sz: The size of the vectors for which to define the basis.
%    mask: A boolean mask of size sz indicating which coordinates to include
%       in the basis (default true(sz)).
%
% Output
%    basis: A Dirac basis object corresponding to the parameters.
%
% Description
%    The returned basis object can be used to project a vector x in the
%    standard basis onto the basis or to map a coefficient vector v into
%    the standard basis. This is achieved using the basis.expand and
%    basis.evaluate functions, respectively. The fields basis.sz denotes
%    the size of the vectors in the standard basis while basis.count denotes
%    the number of vectors in the basis.
%
%    For example, we can generate a random vector in the standard basis and
%    project it onto the Dirac basis through
%
%       basis = dirac_basis(sz, mask);
%       x = randn([sz 1]);
%       v = basis_expand(basis, x);
%
%    Likewise, we can map a random coefficient vector into the standard basis
%    using
%
%       v = randn(basis.count, 1);
%       x = basis_evaluate(basis, v);

% Author
%    Joakim Anden <janden@flatironinstitute.org>

function basis = dirac_basis(sz, mask)
    if nargin < 2 || isempty(mask)
        mask = true(sz);
    end

    basis = struct();

    basis.type = basis_type_dirac();

    basis.sz = sz;
    basis.count = sum(mask(:));

    basis.mask = mask;
end
