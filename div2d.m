function d = div2d(G)

[numRow, ~] = size(G);
f = zeros(numRow,1);

f(1) = -1;
f(2) = 1;

F = diag(fft(f));
gradX = ifft(F*fft(squeeze(G(:,:,1))),[], 1);

gradY = ifft(F*fft(squeeze(G(:,:,2))',[],1))';

d = gradX + gradY;
end