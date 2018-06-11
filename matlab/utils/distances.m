% DISTANCES Calculate distances between vectors
%
% Usage
%    d = distances(vecs);
%
% Input
%    vecs: A p-by-n array of vectors, arranged as columns.
%
% Output
%    d: An n-by-n distance matrix whose ijth entry corresponds to the
%       distance between the ith and the jth columns of `vecs`.

% Author
%    Joakim Anden <janden@flatironinstitute.org>

function d = distances(vecs)
    n2 = sum(abs(vecs).^2, 1);

    d2 = bsxfun(@plus, n2, n2') - 2*vecs'*vecs;

    d = sqrt(max(0, d2));
end
