% test to see whether we can get back the results of x given y
x = rand(5);
%y = grad2d(x);
[dx, dy] = gradient_xy(x);

y = recover_img(dx, dy);
ccc = 0;
% % try to get back x, y
% base = eye(5);
% Fbase = fft2(base);
% 
% grad_filter = zeros(5,5);
% grad_filter(1,1) = -1;
% grad_filter(end,1) = 1;
% F1 = fft2(grad_filter);
% 
% tmpx = ifft2(F1.*fft2(x));
% 
% F2 = fft2(conj(grad_filter));
% 
% 
% grad_filter1 = zeros(5,5);
% grad_filter1(1,1) = -1;
% grad_filter1(1, end) = 1;
% F12 = fft2(grad_filter1);
% F22 = fft2(conj(grad_filter1));
% 
% tmpy = ifft2(F12.*fft2(x));
% 
% tmp = (F2.*F1 + F22.*F12)
% D = 1./(tmp);
% D(tmp == 0) = 0;
% tmpy = Fbase'*D.*Fbase*x;
% 
% 
% D = F2*F1;
% D = pinv(D);
% Z = Fbase'*D*Fbase*x

