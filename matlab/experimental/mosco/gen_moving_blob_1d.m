function outmat = gen_moving_blob_1d(L,sigma)
    blob = fspecial('gaussian', [1 L], sigma);
    outmat = zeros([L L]);
    for i=1:L
        outmat(i,:) = circshift(blob, i-1);
    end
    outmat = outmat + normrnd(0,0.01,[L L]);
end

