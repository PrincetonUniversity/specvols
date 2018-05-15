% IM_TRANSLATE Translate image by shifts
%
% Usage
%    im_translated = im_translate(im, shifts);
%
% Input
%    im: An array of size L-by-L-by-n containing images to be translated.
%    shifts: An array of size 2-by-n specifying the shifts in pixels.
%       Alternatively, it can be a column vector of length 2, in which case
%       the same shifts is applied to each image.
%
% Output
%    im_translated: The images translated by the shifts, with periodic
%       boundaries.

% Author
%    Joakim Anden <janden@flatironinstitute.org>

function im_translated = im_translate(im, shifts)
    n_im = size(im, 3);

    n_shifts = size(shifts, 2);

    if size(shifts, 1) ~= 2
        error('Input `shifts` must be of size 2-by-n.');
    end

    if n_shifts ~= 1 && n_im ~= n_shifts
        error(['The number of shifts must be 1 or match the number of ' ...
            'images.']);
    end

    if size(im, 1) ~= size(im, 2)
        error('Images must be square.');
    end

    L = size(im, 1);

    grid1d = ifftshift(ceil(-L/2:L/2-1))*2*pi/L;

    [om_x, om_y] = ndgrid(grid1d, grid1d);

    phase_shifts = ...
        bsxfun(@times, om_x, permute(-shifts(1,:), [1 3 2])) + ...
        bsxfun(@times, om_y, permute(-shifts(2,:), [1 3 2]));

    mult_f = exp(-i*phase_shifts);

    im_f = fft2(im);

    im_translated_f = bsxfun(@times, im_f, mult_f);

    im_translated = ifft2(im_translated_f);
    im_translated = real(im_translated);
end
