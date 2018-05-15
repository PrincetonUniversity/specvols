% IS_ABSOLUTE_PATHNAME Check whether a path is absolute
%
% Usage
%    ret = is_absolute_pathname(pathname);
%
% Input
%    pathname: A path.
%
% Output
%    ret: `true` if the path is absolute, `false` if not.
%
% Note
%    The implementation of this function is taken from that of
%    `is_absolute_filename` in GNU Octave, which in turn is adapted from that
%    of GNU Bash.

% Author
%    Joakim Anden <janden@flatironinstitute.org>

function ret = is_absolute_pathname(pathname)
    ret = false;

    if length(pathname) > 0 && pathname(1) == filesep
        ret = true;
    end

    if ispc()
        if length(pathname) == 2 && ...
            is_alpha(pathname(1)) && pathname(2) == ':'

            ret = true;
        end

        if length(pathname) > 2 && ...
            is_alpha(pathname(1)) && pathname(2) == ':' && ...
            pathname(3) == filesep

            ret = true;
        end
    end
end
