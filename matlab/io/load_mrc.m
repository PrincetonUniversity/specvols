% LOAD_MRC Open, extract slices from, and close an MRC file
%
% Usage
%    [x, header] = load_mrc(filename, start, num);
%
% Input
%    filename: The filename of the MRC file.
%    start: The starting number of the slices to be extracted (default 1).
%    num: The number of slices to extract. If set to Inf, all slices after
%       start will be returned (default Inf).
%
% Output
%    x: The slices, arranged in an array of size N(1)-by-N(2)-by-n, where
%       N are the slice dimensions from header.N and n is the number of slices
%       extracted.
%    header: The header structure obtained from the MRC file.
%
% Description
%    This function provides a wrapper for the mrc_open, mrc_skip, mrc_read
%    and mrc_close functions.

% Author
%    Joakim Anden <janden@flatironinstitute.org>

function [x, header] = load_mrc(filename, start, num)
    if nargin < 2 || isempty(start)
        start = 1;
    end

    if nargin < 3 || isempty(num)
        num = Inf;
    end
    
    mrc = mrc_open(filename);

    mrc_skip(mrc, start-1);

    x = mrc_read(mrc, num);

    mrc_close(mrc);

    header = mrc.header;
end
