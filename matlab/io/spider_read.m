% SPIDER_READ Read from a SPIDER file
%
% Usage
%    x = spider_read(spider, n);
%
% Input
%    spider: The SPIDER file structure obtained from spider_open.
%    n: The number of slices to read in. If set to Inf, all remaining slices
%       in the file are read (default Inf).
%
% Output
%    x: An array of size nx-by-ny-by-n, where nx and nx are the slice
%       dimensions from the SPIDER header and n is the number of slices
%       actually read.

function x = spider_read(spider, n)
    if nargin < 2 || isempty(n)
        n = Inf;
    end

    N = [spider.header.nx spider.header.ny];

    [x, count] = fread(spider.fd, prod(N)*n, 'float32=>float32');

    x = reshape(x, [N count/prod(N)]);
end
