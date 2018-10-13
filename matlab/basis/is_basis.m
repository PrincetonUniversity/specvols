% IS_BASIS Checks whether object is a valid basis
%
% Usage
%    b = is_basis(basis);
%
% Input
%    basis: The object to be tested.
%
% Output
%    b: True is `basis` is a basis object, such as those returned by
%       `dirac_basis` or `fb_basis`.

% Author
%    Joakim Anden <janden@flatironinstitute.org>

function b = is_basis(basis)
    b = true;

    b = b && isstruct(basis);

    required_fields = {'type', 'sz', 'count'};

    b = b && all(isfield(basis, required_fields));

    b = b && ndims(basis.sz) == 2 && size(basis.sz, 1) == 1 && ...
        all(floor(basis.sz) == basis.sz) && all(basis.sz > 0);
    b = b && numel(basis.count) == 1 && ...
        floor(basis.count) == basis.count && basis.count > 0;

    allowed_types = [basis_type_dirac() ...
        basis_type_fb() ...
        basis_type_matrix() ...
        basis_type_ffb()];

    b = b && any(ismember(basis.type, allowed_types));
end
