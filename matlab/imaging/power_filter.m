% POWER_FILTER Takes the power of a filter
%
% Usage
%    filter = power_filter(original_filter, p);
%
% Input
%    original_filter: The filter object (or array of filter objects) whose
%       power is to be computed.
%    p: The exponent.
%
% Output
%    filter: The filter obtained by taking `original_filter` to the pth power.
%       If the value of `original_filter` at `k` is `h(k)`, `filter` will have
%       the value `h(k)^p`.

% Author
%    Joakim Anden <janden@flatironinstitute.org>

function filter = power_filter(original_filter, p)
    if nargin < 2 || isempty(p)
        p = 1;
    end

    filter = struct();

    filter(1:numel(original_filter)) = struct();

    for k = 1:numel(original_filter)
        filter(k).type = filter_type_power();
        filter(k).dim = original_filter(k).dim;
        filter(k).radial = original_filter(k).radial;
        filter(k).original_filter = original_filter(k);
        filter(k).p = p;
    end
end
