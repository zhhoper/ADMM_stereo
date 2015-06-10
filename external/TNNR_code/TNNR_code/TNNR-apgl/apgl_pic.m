function [ret fullpsnr] = apgl_pic(image_name,matrix_pic,matrix_mask,to_create_eps_pic,lower_R,upper_R)
%% nothing needs to be changed below
%Xfull = double(imread(pic_name));
%Xfull = imresize(Xfull,[resizedm,resizedn]);
Xfull = matrix_pic;
mask = matrix_mask;
Xmiss = Xfull.*mask;
[m, n, dim] = size(Xfull);
known = Xmiss(:,:,1) > 0;
%missing = ones(m,n) - known;
Psnr = zeros(50,1);
time_cost = zeros(50,1);
iterations_cost = zeros(50,1);

best_psnr = 0;
best_rank = 0;
%X_rec = zeros(m,n,
Xrecover = zeros(m,n,3,50);
number_of_out_iter = 10;
%eps = 100;
eps = 0.1;
lambda = 0.06;
fullpsnr = [];

figure(2);
subplot(2,1,1);
imshow(Xmiss/255); % show the missing image.
xlabel('image with missing pixels');


fprintf('now is running under rank(r)=                   ');
for R = lower_R:upper_R
    tic;
    for i = 1:3
        fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b%2d, channel(rgb) %1d',R,i);
        X = Xmiss(:,:,i);
        M = X;
        X_rec = zeros(size(Xmiss));
        for out_iter = 1:number_of_out_iter
            %out_iter
            [u, sigma, v] = svd(X);
            A = u(:,1:R)'; B = v(:,1:R)';
            
            [X_rec(:,:,out_iter),iter_count] = my_apgl(A,B,X,M,known,eps,lambda);
            %fullpsnr = [fullpsnr localpsnr];
            iterations_cost(R) = iterations_cost(R) + iter_count;
            
            if(out_iter>=2 && norm(X_rec(:,:,out_iter)-X_rec(:,:,out_iter-1),'fro')/norm(M,'fro')<0.01)
                X = X_rec(:,:,out_iter);
                break;
            end
            X = X_rec(:,:,out_iter);
        end
        Xrecover(:,:,i,R) = X;
    end
    time_cost(R) = toc;
    tem_recover=max(Xrecover(:,:,:,R),0);
    tem_recover = min(tem_recover,255);
    Psnr(R) = PSNR(Xfull,tem_recover,ones(size(mask))-mask);
    if(best_psnr<Psnr(R))
        best_psnr = Psnr(R);
        best_rank = R;
    end
    
end


subplot(2,1,2); % show the recovery image
Xrecover=max(Xrecover(:,:,:,best_rank),0);
Xrecover = min(Xrecover,255);
imshow(Xrecover/255);
xlabel('recovered image by TNNR-apgl');





if(to_create_eps_pic)
        figure;
        Xrecover=max(Xrecover(:,:,:,best_rank),0);
        Xrecover = min(Xrecover,255);
        imshow(Xrecover/255);
        str_R = sprintf('%d',best_rank);
        str_error = sprintf('%.2f',best_psnr);
        cd result
        saveas(gcf,strcat(image_name,'Rankis',str_R,'PSNRis',num2str(str_error),'APGL','.eps'),'psc2');
        cd ..
        close all;
end
ret.time = time_cost;
ret.iterations = iterations_cost;
ret.psnr = best_psnr;
ret.Psnr = Psnr;
ret.rank = best_rank;
end