% LOAD_IMAGE Open, extract slices from, and close an image file
%
% Usage
%    x = load_image(filename, start, num);
%
% Input
%    filename: The filename of the MRC/MRCS or SPIDER file.
%    start: The starting number of the slices to be extracted (default 1).
%    num: The number of slices to extract. If set to Inf, all slices after
%       start will be returned (default Inf).
%
% Output
%    x: The slices, arranged in an array of size N(1)-by-N(2)-by-n, where
%       N are the slice dimensions and n is the number of slices extracted.
%
% Description
%    This function provides a wrapper for the `load_mrc` (called for
%    extensions '.mrcs' and '.mrc') and `load_spider` (called for all other
%    extensions) functions.

% Author
%    Joakim Anden <janden@flatironinstitute.org>

function x = load_image(filename, start, num)
    [~, ~, ext] = fileparts(filename);

    if any(strcmp(ext, {'.mrc', '.mrcs'}))
        x = load_mrc(filename, start, num);
    else
        x = load_spider(filename, start, num);
    end
end
