% IMAGE_INFO Extract information from image file
%
% Usage
%    info = image_info(filename);
%
% Input
%    filename: The filename of the MRC/MRCS or SPIDER file.
%
% Output
%    info: A struct containing the fields:
%           - 'sz': The size of the image stack contained in the file.
%           - 'precision': The datatype of the images in the file.

% Author
%    Joakim Anden <janden@flatironinstitute.org>

function info = image_info(filename)
    [~, ~, ext] = fileparts(filename);

    info = struct();

    if any(strcmp(lower(ext), {'.mrc', '.mrcs'}))
        mrc = mrc_open(filename);

        info.sz = mrc.header.N;
        info.precision = mrc.data_type;

        mrc_close(mrc);
    else
        spider = spider_open(filename);

        info.sz = [spider.header.nx spider.header.ny spider.header.nz];
        info.precision = 'single';

        spider_close(spider);
    end
end
