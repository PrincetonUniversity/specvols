% IS_DOUBLE Determine whether string is double or not
%
% Usage
%    b = is_double(str);
%
% Input
%    str: The string to check.
%
% Output
%    b: True if the string represents a double, false if not.

% Author
%    Joakim Anden <janden@flatironinstitute.org>

function b = is_double(str)
    b = all((str>='0' & str<='9') | str=='.' | ...
             str=='e' | str=='E' | ...
             str=='+' | str=='-');
end
