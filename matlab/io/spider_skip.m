% SPIDER_SKIP Skips slices in a SPIDER file
%
% Usage
%    spider_skip(spider, s);
%
% Input
%    spider: The SPIDER file structure obtained from spider_open.
%    s: The number of slices to skip.

function spider_skip(spider, s)
    N = [spider.header.nx spider.header.ny];

    fseek(spider.fd, prod(N)*s*4, 0);
end
