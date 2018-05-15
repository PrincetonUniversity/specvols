% FAST_SPLIT Split string fast
%
% Usage
%    parts = fast_split(str, ch);
%
% Input
%    str: The string to be split.
%    ch: A cell array of characters which should be used to split the string.
%
% Output
%    parts: A cell array of strings containing the parts of the string sep-
%       arated by the delimiters in ch. Adjacent instances of delimiters are
%       treated as a single delimiter.

% Author
%    Joakim Anden <janden@flatironinstitute.org>

function parts = fast_split(str, ch)
    if nargin < 2
        ch = {' '};
    end

    inds = [];
    for k = 1:length(ch)
        inds = [inds find(str==ch{k})];
    end

    inds1 = [0 inds];
    inds2 = [inds length(str)+1];

    mask = (inds2-inds1)>1;
    inds1 = inds1(mask);
    inds2 = inds2(mask);

    parts = cell(1,length(inds1));

    for k = 1:length(parts)
        parts{k} = str(inds1(k)+1:inds2(k)-1);
    end
end
