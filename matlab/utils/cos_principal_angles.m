% COS_PRINCIPAL_ANGLES Cosines of principal angles between subspaces
%
% Usage
%    cos_theta = cos_principal_angles(A, B);
%
% Input
%    A, B: Two p-by-k matrices containing bases for two k-dimensional
%       subspaces of R^p.
%
% Output
%    cos_theta: The cosines of the principal angles between of the subspaces
%       defined by A and B.

% Author
%    Joakim Anden <janden@flatironinstitute.org>

function cos_theta = cos_principal_angles(A, B)
    [A, ~] = qr(A, 0);
    [B, ~] = qr(B, 0);

    cos_theta = svd(A'*B);
end
