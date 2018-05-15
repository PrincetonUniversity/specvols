% MRC_SKIP Skips slices in an MRC file
%
% Usage
%    mrc_skip(mrc, s);
%
% Input
%    mrc: The MRC structure obtained from mrc_open.
%    s: The number of slices to skip.

% Author
%    Joakim Anden <janden@flatironinstitute.org>

function mrc_skip(mrc, s)
    N = mrc.header.N(1:2);

    fseek(mrc.fd, prod(N)*s*mrc.data_width, 0);
end
