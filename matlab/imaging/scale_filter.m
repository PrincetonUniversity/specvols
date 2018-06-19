% SCALE_FILTER Scale filter by a constant factor
%
% Usage
%    filter = scale_filter(filter, c);
%
% Input
%    filter: The filter object (or array of filter objects) to be scaled.
%    c: The scaling factor. For c < 1, it dilates the filter(s) in frequency,
%       while for c > 1, it compresses (default 1).
%
% Output
%    filter: The same filter object (or array of filter objects), scaled by
%       `c`.

function filter = scale_filter(filter, c)
    if nargin < 2 || isempty(c)
        c = 1;
    end

    if numel(filter) > 1
        filter = arrayfun(@(f)(scale_filter(f, c)), filter);
        return;
    end

    if filter.type == filter_type_scalar()
        return;
    elseif filter.type == filter_type_ctf()
        filter.ctf_params.pixel_size = filter.ctf_params.pixel_size*c;
        return;
    elseif filter.type == filter_type_power()
        filter.original_filter = scale_filter(filter.original_filter, c);
    elseif filter.type == filter_type_gaussian()
        filter.scale = c*filter.scale;
    end
end
