% VECMAT_TO_VOLMAT Roll up vector matrices into volume matrices
%
% Usage
%    volmat = vecmat_to_volmat(vecmat);
%
% Input
%    vecmat: A vector matrix of size L1^3-by-L2^3-by-... .
%
% Output
%    volmat: A volume "matrix" of size L1-by-L1-by-L1-by-L2-by-L2-by-
%       L2-by-... .
%
% See also
%    volmat_to_vecmat

% Author
%    Joakim Anden <janden@flatironinstitute.org>

function volmat = vecmat_to_volmat(vecmat)
    sz = size(vecmat);

    L1 = round(sz(1)^(1/3));
    L2 = round(sz(2)^(1/3));

    if L1^3 ~= sz(1) || L2^3 ~= sz(2)
        error('Volumes must be cubic.');
    end

    volmat = reshape(vecmat, [L1*ones(1,3) L2*ones(1,3) sz(3:end)]);
end
