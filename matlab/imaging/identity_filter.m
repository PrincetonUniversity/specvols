% IDENTITY_FILTER Create identity filter
%
% Usage
%    filter = identity_filter(dim);
%
% Input
%    dim: The dimension of the filter (default 2).
%
% Output
%    filter: A filter structure corresponding to the identity filter which
%       leaves a signal unchanged.
%
% See also
%    scalar_filter

% Author
%    Joakim Anden <janden@flatironinstitute.org>

function filter = identity_filter(dim)
    if nargin < 1 || isempty(dim)
        dim = 2;
    end

    filter = scalar_filter(1, dim);
end
