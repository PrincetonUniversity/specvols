% FUN_SVDS Find top singular values and vectors of a linear mapping
%
% Usage
%    [U, S, V] = fun_svds(A_fun, At_fun, N, M, K, opts);
%
% Input
%    A_fun: A function handle applying a matrix A to a column vector of length
%       N.
%    At_fun: A function handle applying the transpose of A to a column vector
%       of length M.
%    N, M: The dimenions of the matrix A (N columns and M rows).
%    K: The number of desired singular values and vectors (default 6).
%    opts: An options structure to be passed to eigs (default empty).
%
% Output
%    U: The top K left singular vectors of A arranged as columns.
%    S: The top K singular values of A arranged in a diagonal matrix.
%    V: The top K right singular vectors of A arranged as columns.

function [U, S, V] = fun_svds(A_fun, At_fun, N, M, K, opts)
    if nargin < 5 || isempty(K)
        K = 6;
    end

    if nargin < 6 || isempty(opts)
        opts = struct();
    end

    fun = @(x)(cat(1, A_fun(x(M+1:M+N,:)), At_fun(x(1:M,:))));

    [V0, D0] = eigs(fun, N+M, K, 'la', opts);

    U = V0(1:M,1:K);
    U = bsxfun(@times, U, 1./sqrt(sum(abs(U).^2, 1)));

    V = V0(M+1:M+N,1:K);
    V = bsxfun(@times, V, 1./sqrt(sum(abs(V).^2, 1)));

    sigma = diag(abs(D0));

    [sigma, I] = sort(sigma, 'descend');

    S = diag(sigma);
    U = U(:,I);
    V = V(:,I);
end