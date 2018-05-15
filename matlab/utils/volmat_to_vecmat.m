% VOLMAT_TO_VECMAT Unroll volume matrices to vector matrices
%
% Usage
%    vecmat = volmat_to_vecmat(volmat);
%
% Input
%    volmat: A volume "matrix" of size L1-by-L1-by-L2-by-L2-by-L2-by-
%       L2-by-... .
%
% Output
%    vecmat: A vector matrix of size L1^3-by-L2^3-by-... .
%
% See also
%    vecmat_to_volmat

% Author
%    Joakim Anden <janden@flatironinstitute.org>

function vecmat = volmat_to_vecmat(volmat)
    sz = size(volmat);

    if length(sz) < 6
        error('Volume matrix must be at least six-dimensional.');
    end

    L1 = sz(1);
    L2 = sz(4);

    if all(sz(1:6) ~= [L1*ones(1, 3) L2*ones(1, 3)])
        error('Volumes must be cubic.');
    end

    vecmat = reshape(volmat, [L1^3 L2^3 sz(7:end)]);
end
