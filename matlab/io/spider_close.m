% SPIDER_CLOSE Closes a SPIDER file
%
% Usage
%    spider_close(spider);
%
% Input
%    spider: The SPIDER file structure obtained from spider_open.

function spider_close(spider)
    fclose(spider.fd);
end
