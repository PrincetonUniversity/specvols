% SCALAR_FILTER Create scalar filter
%
% Usage
%    filter = scalar_filter(value, dim);
%
% Input
%    value: The multiplier value of the filter (default 1).
%    dim: The dimension of the filter (default 2).
%
% Output
%    filter: A filter structure corresponding to the scalar filter which
%       multiplies a signal by a certain value.
%
% See also
%    identity_filter

% Author
%    Joakim Anden <janden@flatironinstitute.org>

function filter = scalar_filter(value, dim)
    if nargin < 2 || isempty(dim)
        dim = 2;
    end

    if nargin < 1 || isempty(value)
        value = 1;
    end

    filter.type = filter_type_scalar();
    filter.dim = dim;
    filter.radial = true;

    filter.value = value;
end
