% VECMAT_TO_VOLMAT Converts a vectorized volume matrix into a volume matrix
%
% Usage
%    volmat = vec_to_volmat(vec);
%
% Input
%    vec: A vector of size N^6-by-... .
%
% Output
%    volmat: A volume "matrix" of size N-by-N-by-N-by-N-by-N-by-N-by-... .

% Author
%    Joakim Anden <janden@flatironinstitute.org>

function volmat = vec_to_volmat(vec)
    volmat = vecmat_to_volmat(vec_to_mat(vec));
end
