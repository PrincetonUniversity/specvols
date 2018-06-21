% ARRAY_FILTER Create filter from an array
%
% Usage
%    filter = array_filter(filter_f);
%
% Input
%    filter_f: The transfer function of the filter in the form of an array of
%       one or two dimensions.
%
% Output
%    filter: A filter structure corresponding to the filter with the specified
%       transfer function.

% Author
%    Joakim Anden <janden@flatironinstitute.org>

function filter = array_filter(filter_f)
    dim = ndims(filter_f);

    if dim == 2 && size(filter_f, 2) == 1
        dim = 1;
    end

    if dim ~= 1 && dim ~= 2
        error('Only dimensions 1 and 2 supported.');
    end

    filter = struct();
    filter.type = filter_type_array();
    filter.dim = dim;
    filter.radial = false;

    sz = size(filter_f);

    if dim >= 1 && mod(size(filter_f, 1), 2) == 0
        filter_f = cat(1, filter_f, filter_f(1,end:-1:1));
    end

    if dim >= 2 && mod(size(filter_f, 2), 2) == 0
        filter_f = cat(2, filter_f, filter_f(end:-1:1,1));
    end

    filter.filter_f = filter_f;
    filter.sz = sz;
    filter.scale = 1;
end
