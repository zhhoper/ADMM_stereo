function G = grad2d(X)
% compute gradient of a 2-dim matrix using fft

[numRow, ~] = size(X);
f = zeros(numRow,1);

f(1) = -1;
f(end) = 1;

F = diag(fft(f));
gradX = ifft(F*fft(X),[], 1);

gradY = ifft(F*fft(X'),[],1)';

G = cat(3, gradX, gradY);
end


