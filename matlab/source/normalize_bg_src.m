% NORMALIZE_BG_SRC Normalize background of source
%
% Usage
%    src = normalize_bg_src(original_src, bg_radius, do_ramp);
%
% Input
%    original_src: A source whose images are to have their background
%       normalized.
%    bg_radius: The radius of the disk whose complement specifies the
%       background noise domain (default 1).
%    do_ramp: Fit a ramping background to the data and subtract. If false, a
%       constant background level is used instead (default true).
%
% Output
%    src: A source containing the same images as `original_src`, but with their
%       backgrounds normalized using `im_normalize_bg`.
%
% See also
%    im_normalize_bg

% Author
%    Joakim Anden <janden@flatironinstitute.org>

function src = normalize_bg_src(original_src, bg_radius, do_ramp)
    if nargin < 2 || isempty(bg_radius)
        bg_radius = 1;
    end

    if nargin < 3 || isempty(do_ramp)
        do_ramp = true;
    end

    src = struct();

    src.type = src_type_normalized_bg();

    src.original_src = original_src;

    src.L = original_src.L;
    src.n = original_src.n;
    src.precision = original_src.precision;

    src.params = original_src.params;

    src.bg_radius = bg_radius;
    src.do_ramp = do_ramp;
end
