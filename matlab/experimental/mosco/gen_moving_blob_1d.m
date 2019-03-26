function outmat = gen_moving_blob_1d(L,blobwidth, noise_sigma)
    blob = fspecial('gaussian', [1 L], blobwidth);
    outmat = zeros([L L]);
    for i=1:L
        outmat(i,:) = circshift(blob, i-1);
    end
    outmat = outmat + normrnd(0,noise_sigma,[L L]);
    outmat = outmat - sum(outmat,1)/L;
end

