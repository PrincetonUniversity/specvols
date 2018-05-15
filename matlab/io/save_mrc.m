% SAVE_MRC Save a data to an MRC file
%
% Usage
%    save_mrc(filename, x, pixel_size, header);
%
% Input
%    filename: The filename of the MRC file.
%    x: The data (image stack, volume, etc.) to be written.
%    pixel_size: The pixel size in angstroms (default 0).
%    header: A header structure for the MRC file (default generated from x and
%       pixel_size).

% Author
%    Joakim Anden <janden@flatironinstitute.org>

function save_mrc(filename, x, pixel_size, header)
    if nargin < 3 || isempty(pixel_size)
        pixel_size = 0;
    end

    if nargin < 4 || isempty(header)
        header = default_header(x, pixel_size);
    end

    fd = fopen(filename, 'w', 'ieee-le');

    data_types = {'int8', 'int16', 'single'};
    data_widths = [1 2 4];

    data_type = data_types{header.mode+1};
    data_width = data_widths(header.mode+1);

    write_header(fd, header);

    fwrite(fd, x(:), data_type);

    fclose(fd);
end

function header = default_header(x, pixel_size)
    header.N = size(x);

    if strcmp(class(x), 'int8')
        header.mode = 0;
    elseif strcmp(class(x), 'int16')
        header.mode = 1;
    elseif strcmp(class(x), 'single')
        header.mode = 2;
    else
        error(['Invalid precision. Map must be ''int8'', ''int16'', or ' ...
            '''single''.']);
    end

    header.nstart = -floor(size(x)/2);

    header.m = size(x);

    header.cella = pixel_size*size(x);

    header.cellb = 90*ones(1, 3);

    header.mapcrs = 1:3;

    header.dmin = min(x(:));
    header.dmax = max(x(:));
    header.dmean = mean(x(:));

    % TODO: Check for volume/stack/stack of volumes.
    header.ispg = 0;

    header.nsymbt = 0;

    header.extra = zeros(1, 25);

    header.origin = zeros(1, 3);

    % TODO: Check for endianess and change accordingly.
    header.map = 'MAP ';
    header.machst = typecast(uint8([68 65 0 0]), 'uint32');

    header.rms = std(x(:));

    header.nlabl = 0;

    header.label = char(zeros(10, 80));
end

function write_header(fd, header)
    fwrite(fd, header.N, 'uint32');

    fwrite(fd, header.mode, 'uint32');

    fwrite(fd, header.nstart, 'int32');

    fwrite(fd, header.m, 'uint32');

    fwrite(fd, header.cella, 'float32');

    fwrite(fd, header.cellb, 'float32');

    fwrite(fd, header.mapcrs, 'uint32');

    fwrite(fd, header.dmin, 'float32');
    fwrite(fd, header.dmax, 'float32');
    fwrite(fd, header.dmean, 'float32');

    fwrite(fd, header.ispg, 'uint32');

    fwrite(fd, header.nsymbt, 'uint32');

    fwrite(fd, header.extra, 'uint32');

    fwrite(fd, header.origin, 'float32');

    fwrite(fd, header.map, 'uchar');

    fwrite(fd, header.machst, 'uint32');

    fwrite(fd, header.rms, 'float32');

    fwrite(fd, header.nlabl, 'uint32');

    fwrite(fd, header.label', 'uchar');
end
