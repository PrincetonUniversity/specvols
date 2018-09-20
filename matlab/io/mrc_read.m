% MRC_READ Read from MRC file
%
% Usage
%    x = mrc_read(mrc, n);
%
% Input
%    mrc: The MRC structure obtained from mrc_open.
%    n: The number of slices to read in. If set to Inf, all remaining slices
%       in the file are read (default Inf).
%
% Output
%    x: An array of size N(1)-by-N(2)-by-n, where N are the slice dimensions
%       from the MRC header found in mrc.header.N and n is the number of
%       slices acutally read.

% Author
%    Joakim Anden <janden@flatironinstitute.org>

function x = mrc_read(mrc, n)
    if nargin < 2 || isempty(n)
        n = Inf;
    end
    N = mrc.header.N(1:2);

    if isinf(n)
        n = mrc.header.N(3);
    end

    [x, count] = fread(mrc.fd, prod(N)*n, ['*' mrc.data_type]);

    x = reshape(x, [N' n]);
end
