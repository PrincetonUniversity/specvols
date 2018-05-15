% FPEEK Peek ahead in the input stream
%
% Usage
%    c = fpeek(fid);
%
% Input
%    fid: The file identifier to read from.
%
% Output
%    c: The next character to be read.
%

% Author
%    Joakim Anden <janden@flatironinstitute.org>

function c = fpeek(fid)
    c = fread(fid, 1, 'char');
    fseek(fid, -1, 0);

    if ~isempty(c)
        c = char(c);
    end
end
