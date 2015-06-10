% written by debingzhang
% if you have any questions, please fell free to contact
% debingzhangchina@gmail.com

% version 1.0, 2012.11.16








clear;clc;
%% for experiment use
pic_list = {'re1.jpg','re2.jpg','re3.jpg','re4.jpg','re5.jpg','re6.jpg','re7.jpg','re8.jpg','re9.jpg','re10.jpg','re11.jpg'};

imagenum = 10;
pic_name = pic_list{imagenum};

cd pic
Xfull = double(imread(pic_name));
cd ..

sizem = size(Xfull,1);
sizen = size(Xfull,2);

ind = randint(sizem,sizen,10,2011);

oldind = ind;

iter = zeros(10,3);
err = zeros(10,3);
tim = zeros(10,3);

kk = 5; %       
disp('50% entries are missing');
ind = (oldind < kk);   % 50% entries are missing, controlled by kk.

mask(:,:,1)=ind;mask(:,:,2)=ind;mask(:,:,3)=ind;

to_create_eps_pic = 1;
lower_R = 9; upper_R = 9;  % as we test r=9 maybe the best for re10.jpg, for other images, sometimes, it is needed to test all small rs.
fprintf('now is running admm optimization method to recovery an image with missing pixels\n');

cd TNNR-admm
tic;
[admmret psnradmm]= admm_pic(pic_name,Xfull,mask,0,lower_R,upper_R);
cd ..
%toc;
admm_num_iteration = max(admmret.iterations);
admm_psnr_error = admmret.psnr;
admm_time_cost = toc;

fprintf('\n TNNR-admm: psnr(%.3f), time(%.1fs), iterations(%d)\n',admm_psnr_error,admm_time_cost,admm_num_iteration);

fprintf('\n\nnow is running apgl optimization method to recovery an image with missing pixels\n');
cd TNNR-apgl
tic;
[apglret psnrapgl]= apgl_pic(pic_name,Xfull,mask,0,lower_R,upper_R);
cd ..
%toc;
apgl_num_iteration = max(apglret.iterations);
apgl_psnr_error = apglret.psnr;
apgl_time_cost = toc;

fprintf('\n TNNR-apgl: psnr(%.3f), time(%.1fs), iterations(%d)\n',apgl_psnr_error,apgl_time_cost,apgl_num_iteration);


