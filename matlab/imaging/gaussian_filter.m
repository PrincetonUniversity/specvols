% GAUSSIAN_FILTER Create Gaussian filter
%
% Usage
%    filter = gaussian_filter(scale, dim);
%
% Input
%    scale: The scale of the filter in the frequency domain (default pi).
%    dim: The dimension of the filter (default 2).
%
% Output
%    filter: A filter structure corresponding to the Gaussian filter whose
%       transfer function equals
%
%          exp(-k^2/(2*scale^2)),
%
%       where `k` is the radial frequency normalized so that the Nyquist
%       frequency is `pi`.
%
% See also
%    scalar_filter

% Author
%    Joakim Anden <janden@flatironinstitute.org>

function filter = gaussian_filter(scale, dim)
    if nargin < 2 || isempty(dim)
        dim = 2;
    end

    if nargin < 1 || isempty(scale)
        scale = pi;
    end

    filter.type = filter_type_gaussian();
    filter.dim = dim;
    filter.radial = true;

    filter.scale = scale;
end
