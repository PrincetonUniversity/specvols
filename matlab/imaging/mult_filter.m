% MULT_FILTER Multiplies two filters
%
% Usage
%    filter = mult_filter(original_filter1, original_filter2);
%
% Input
%    original_filter1, original_filter2: The filter objects whose product is
%       to be computed. The first filter, `original_filter1`, may be an array
%       of filters, win which case the output will be an array as well, with
%       each filter being multiplied by `original_filter2`.
%
% Output
%    filter: The filter obtained by taking `original_filter1` multiplied by
%       `original_filter2`. If the value of `original_filter1` at `k` is
%       `h1(k)` and that of `original_filter2` is `h2(k)`, `filter` will have
%       value `h1(k) * h2(k)`.

% Author
%    Joakim Anden <janden@flatironinstitute.org>

function filter = mult_filter(original_filter1, original_filter2)
    filter = struct();

    filter(1:numel(original_filter1)) = struct();

    for k = 1:numel(original_filter1)
        filter(k).type = filter_type_mult();
        filter(k).dim = original_filter1(k).dim;
        filter(k).radial = ...
            original_filter1(k).radial && original_filter2.radial;
        filter(k).original_filter1 = original_filter1(k);
        filter(k).original_filter2 = original_filter2;
    end
end
