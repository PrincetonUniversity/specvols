% SPIDER_OPEN Opens a SPIDER file for reading
%
% Usage
%    spider = spider_open(filename);
%
% Input
%    filename: A SPIDER filename.
%
% Output
%    spider: A structure representing the SPIDER file containing the file
%        descriptor and header.

function spider = spider_open(filename)
    fd = fopen(filename, 'r', 'ieee-le');

    h = spider_header(fd);

    if ~ismember(h.iform, [1 3 -11 -12 -21 -22])
        fclose(fd);

        fd = fopen(filename, 'r', 'ieee-be');

        h = spider_header(fd);
    end

    fread(fd, h.labbyt-1024, 'char');

    spider = struct();
    spider.fd = fd;
    spider.header = h;
end

function h = spider_header(fd)
    h = struct();

    h.nz = fread(fd, 1, 'float32');
    h.ny = fread(fd, 1, 'float32');
    h.irec = fread(fd, 1, 'float32');
    fread(fd, 1, 'float32');
    h.iform = fread(fd, 1, 'float32');
    h.imami = fread(fd, 1, 'float32');
    h.fmax = fread(fd, 1, 'float32');
    h.fmin = fread(fd, 1, 'float32');
    h.av = fread(fd, 1, 'float32');
    h.sig = fread(fd, 1, 'float32');
    fread(fd, 1, 'float32');
    h.nx = fread(fd, 1, 'float32');
    h.labrec = fread(fd, 1, 'float32');
    h.iangle = fread(fd, 1, 'float32');
    h.tilt = fread(fd, 3, 'float32');
    h.offset = fread(fd, 3, 'float32');
    h.scale = fread(fd, 1, 'float32');
    h.labbyt = fread(fd, 1, 'float32');
    h.lenbyt = fread(fd, 1, 'float32');
    h.istack = fread(fd, 1, 'float32');
    fread(fd, 1, 'float32');
    h.maxim = fread(fd, 1, 'float32');
    h.imgnum = fread(fd, 1, 'float32');
    h.lastindx = fread(fd, 1, 'float32');
    fread(fd, 2, 'float32');
    h.kangle = fread(fd, 1, 'float32');
    h.rot1 = fread(fd, 3, 'float32');
    h.rot2 = fread(fd, 3, 'float32');
    h.pixsiz = fread(fd, 1, 'float32');
    h.ev = fread(fd, 1, 'float32');
    h.proj = fread(fd, 1, 'float32');
    h.mic = fread(fd, 1, 'float32');
    h.num = fread(fd, 1, 'float32');
    h.glonum = fread(fd, 1, 'float32');
    fread(fd, 57, 'float32');
    h.rot3 = fread(fd, 3, 'float32');
    h.langle = fread(fd, 1, 'float32');
    fread(fd, 107, 'float32');
    h.cdat = char(fread(fd, 12, 'char'))';
    h.cdat = h.cdat(1:11);
    h.ctim = char(fread(fd, 8, 'char'))';
    h.ctit = char(fread(fd, 160, 'char'))';
end
