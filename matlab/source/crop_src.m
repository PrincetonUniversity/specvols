% CROP_SRC Crop all images in a source
%
% Usage
%    src = crop_src(original_src, upper_left, lower_right);
%
% Input
%    original_src: The source containing images to be cropped.
%    upper_left: Upper left corner of cropping box.
%    lower_right: Lower right corner of cropping box.
%
% Output
%    src: The source with each image cropped according to `upper_left` and
%       `lower_right`.

% Author
%    Joakim Anden <janden@flatironinstitute.org>

function src = crop_src(original_src, upper_left, lower_right)
    if lower_right(1)-upper_left(1) ~= lower_right(2)-upper_left(2)
        error('crop_src: Cropping must be to a square');
    end

    src = struct();

    src.type = src_type_cropped();

    src.original_src = original_src;

    src.L = lower_right(1)-upper_left(1)+1;
    src.n = original_src.n;
    src.precision = original_src.precision;

    src.coords = [upper_left lower_right];

    src.params = original_src.params;
end
