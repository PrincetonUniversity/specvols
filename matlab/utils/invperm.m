% INVPERM Invert permutation
%
% Usage
%    iperm = invperm(perm);
%
% Input
%    perm: A permutation in the form of an array of length n containing all
%       integers between 1 and n exactly once.
%
% Output
%    iperm: The inverse of `perm`, in the sense that
%
%       iperm(perm) == 1:n
%       perm(iperm) == 1:n .

% Author
%    Joakim Anden <janden@flatironinstitute.org>

function iperm = invperm(perm)
    n = numel(perm);

    [~, iperm] = ismember(1:n, perm);
end
