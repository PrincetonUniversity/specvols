% WHITEN_SRC Apply whitenin filter to source
%
% Usage
%    src = whiten_src(original_src, noise_psd);
%
% Input
%    original_src: A source whose images are to be whitened.
%    noise_psd: The noise PSD of the images, typically estimated using
%       `estimate_noise`.
%
% Output
%    src: A source containing the images as the source, but filtered so that
%       their noise is approximately white.

% Author
%    Joakim Anden <janden@flatironinstitute.org>

function src = whiten_src(original_src, noise_psd)
    src = struct();

    src.type = src_type_filtered();

    src.original_src = original_src;

    src.L = original_src.L;
    src.n = original_src.n;
    src.precision = original_src.precision;

    whiten_filter = power_filter(noise_psd, -1/2);

    src.filters = whiten_filter;
    src.filter_idx = ones(1, src.n);

    src.params = original_src.params;

    src.params.filters = mult_filter(src.params.filters, whiten_filter);
end
