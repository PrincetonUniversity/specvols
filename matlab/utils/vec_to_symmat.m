% VEC_TO_SYMMAT Convert packed lower triangular vector to symmetric matrix
%
% Usage
%    mat = vec_to_symmat(vec);
%
% Input
%    vec: A vector of size N*(N+1)/2-by-... describing a symmetric (or
%       Hermitian) matrix.
%
% Output
%    mat: An array of size N-by-N-by-... which indexes symmetric/Hermitian
%       matrices that occupy the first two dimensions. The lower triangular
%       parts of these matrices consists of the corresponding vectors in vec.

% Author
%    Joakim Anden <janden@flatironinstitute.org>

function mat = vec_to_symmat(vec)
    % Calculate size of matrix.
    N = round(sqrt(2*size(vec, 1)+1/4)-1/2);

    if size(vec, 1) ~= N*(N+1)/2 || N == 0
        error('vector must be of size N*(N+1)/2 for some N > 0');
    end

    [vec, sz_roll] = unroll_dim(vec, 2);

    % Create index matrix.
    I = tril(ones(N*ones(1, 2)));
    I = cumsum(I(:));
    I = tril(reshape(I, N*ones(1, 2)));
    I = I+I'-diag(diag(I));

    if isreal(vec)
        mat = vec(I(:),:);
    else
        % For complex, need to conjugate upper triangle for Hermitian
        % symmetry.
        Jl = tril(true(N*ones(1, 2)));
        Ju = triu(true(N*ones(1, 2)));
        Il = tril(I);
        Iu = triu(I);

        Jl = find(Jl(:));
        Ju = find(Ju(:));
        Il = Il(Il~=0);
        Iu = Iu(Iu~=0);

        mat = zeros([N^2 size(vec, 2)], class(vec));

        mat(Jl,:) = vec(Il(:),:);
        mat(Ju,:) = conj(vec(Iu(:),:));

        % Ensure real diagonal.
        Jd = 1:N+1:N^2;

        mat(Jd,:) = real(mat(Jd,:));
    end
    mat = reshape(mat, [N*ones(1, 2) size(mat, 2)]);

    mat = roll_dim(mat, sz_roll);
end
