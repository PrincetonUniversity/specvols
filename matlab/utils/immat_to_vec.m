% IMMAT_TO_VEC Convert a image matrix into vectorized form
%
% Usage
%    vec = immat_to_vec(immat);
%
% Input
%    volmat: An image "matrix" of size N-by-N-by-N-by-N-by-... .
%
% Output
%   vec: A vector of size N^4-by-... .

% Author
%    Joakim Anden <janden@flatironinstitute.org>

function vec = immat_to_vec(immat)
    vec = mat_to_vec(immat_to_vecmat(immat));
end
