% VECMAT_TO_IMMAT Converts a vectorized image matrix into a image matrix
%
% Usage
%    immat = vec_to_immat(vec);
%
% Input
%    vec: A vector of size N^4-by-... .
%
% Output
%    immat: An image "matrix" of size N-by-N-by-N-by-N-by-... .

% Author
%    Joakim Anden <janden@flatironinstitute.org>

function immat = vec_to_immat(vec)
    immat = vecmat_to_immat(vec_to_mat(vec));
end
