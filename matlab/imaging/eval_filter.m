% EVAL_FILTER Evaluates a filter
%
% Usage
%    h = eval_filter(filter, omega);
%
% Input
%    filter: A filter structure obtained from `identity_filter` or
%       `ctf_filter`. This can also be a struct array, in which case all
%       filters in the array are evaluated.
%    omega: An array of size filter.dim-by-n representing the frequencies at
%       which the filter is to be evaluated. These are normalized so that
%       pi is equal to the Nyquist frequency.
%
% Output
%    h: The value of the filter at the n specified frequencies. If `filter` is
%       an array of length K, `h` is of size n-by-K, where each column corre-
%       sponds to a different filter.

% Author
%    Joakim Anden <janden@flatironinstitute.org>

function h = eval_filter(filter, omega)
    if ~all(isfield(filter, {'type', 'dim', 'radial'}))
        error('`filter` must be a filter structure');
    end

    if numel(filter) > 1
        h = arrayfun(@(f)(eval_filter(f, omega)), filter, ...
            'uniformoutput', false);

        h = cat(2, h{:});

        return;
    end

    if size(omega, 1) ~= filter.dim
        error(sprintf('`omega` must of size %d-by-n.', filter.dim));
    end

    if filter.radial
        omega = sqrt(sum(omega.^2, 1));

        [omega, ~, idx] = unique(omega);
        omega = cat(1, omega, zeros(size(omega), class(omega)));
    end

    if filter.type == filter_type_scalar
        h = filter.value*ones(size(omega));
    elseif filter.type == filter_type_ctf
        h = eval_filter_ctf(filter, omega);
    elseif filter.type == filter_type_power
        h = eval_filter(filter.original_filter, omega).^filter.p;
    elseif filter.type == filter_type_gaussian
        h = exp(-omega(1,:).^2/(2*filter.scale^2));
    elseif filter.type == filter_type_array
        sz = filter.sz;

        omega = omega/filter.scale;
        omega = bsxfun(@times, omega/(2*pi), sz');
        omega = bsxfun(@plus, omega, floor(sz'/2)+1);

        omega = num2cell(omega, 2)';

        h = interpn(filter.filter_f, omega{:}, 'linear', 0);
    else
        error(['`filter.type` must be a valid filter type']);
    end

    if filter.radial
        h = h(idx);
    end

    h = h(:);
end
