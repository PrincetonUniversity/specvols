% IMMAT_TO_VECMAT Unroll image matrices to vector matrice
%
% Usage
%    vecmat = immat_to_vecmat(immat);
%
% Input
%    immat: An image "matrix" of size N1-by-N1-by-N2-by-N2-by-... .
%
% Output
%    vecmat: A vector matrix of size N1^2-by-N2^2-by-... .

% Author
%    Joakim Anden <janden@flatironinstitute.org>

function vecmat = immat_to_vecmat(immat)
    sz = size(immat);

    if length(sz) < 4
        error('image matrix must be at least four-dimensional');
    end

    N1 = sz(1);
    N2 = sz(3);

    if all(sz(1:4) ~= [N1*ones(1, 2) N2*ones(1, 2)])
        error('images must be square');
    end

    vecmat = reshape(immat, [N1^2 N2^2 sz(5:end)]);
end
