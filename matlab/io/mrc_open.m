% MRC_OPEN Opens an MRC file for reading
%
% Usage
%    mrc = mrc_open(filename);
%
% Input
%    filename: MRC filename.
%
% Output
%    mrc: A structure representing the MRC file containing the file descriptor
%       and header.

% Author
%    Joakim Anden <janden@flatironinstitute.org>

function mrc = mrc_open(filename)
    fd = fopen(filename, 'r', 'ieee-le');

    h = mrc_header(fd);

    check_header(h);

    mrc.fd = fd;
    mrc.header = h;

    datatypes = {'int8', 'int16', 'single'};

    mrc.data_type = datatypes{mrc.header.mode+1};

    data_widths = [1 2 4];

    mrc.data_width = data_widths(mrc.header.mode+1);
end

function h = mrc_header(f)
    h.N = fread(f, 3, 'uint32');

    h.mode = fread(f, 1, 'uint32');

    h.nstart = fread(f, 3, 'int32');

    h.m = fread(f, 3, 'uint32');

    h.cella = fread(f, 3, 'float32');

    h.cellb = fread(f, 3, 'float32');

    h.mapcrs = fread(f, 3, 'uint32');

    h.dmin = fread(f, 1, 'float32');
    h.dmax = fread(f, 1, 'float32');
    h.dmean = fread(f, 1, 'float32');

    h.ispg = fread(f, 1, 'uint32');

    h.nsymbt = fread(f, 1, 'uint32');

    h.extra = fread(f, 25, 'uint32');

    h.origin = fread(f, 3, 'float32');

    h.map = char(fread(f, 4, 'uchar'))';

    h.machst = fread(f, 1, 'uint32');

    h.rms = fread(f, 1, 'float32');

    h.nlabl = fread(f, 1, 'uint32');

    h.label = char(fread(f, 10*80, 'uchar'));
    h.label = reshape(h.label, [80 10])';
end

function check_header(h)
    % Check against specification.
    if ~ismember(h.mode, 0:4)
        error('MODE must be in the range 0 through 4');
    end

    if (~all(ismember(h.mapcrs, 1:3)) || numel(unique(h.mapcrs)) ~= 3) && ...
        ~all(h.mapcrs == 0)

        error(['MAPC, MAPR and MAPS must form a permutation of {1, 2, 3} ' ...
            'or all equal 0.']);
    end

    % Check against current capabilities.
    if ismember(h.mode, 3:4)
        error('Transform data (MODE 3 and 4) currently not supported');
    end

    if h.nsymbt > 0
        error('Files with symmetry data currently not supported');
    end
end
