clear all; clc
list=dir('result_new/'); list=list(3:end);
cd('result_new/');



base=[]; admm=[]; miss=0:10:80;

%sort_id=[6 2 3 4 5 1];

for i=1:length(list)
    load(list(i).name);
    disp('....')
    disp(list(i).name)
    
    err_base=main_result.base_error.depth;
    
    
    err_admm=[];
    for j=1:length(main_result.admm)        
        err_admm=[err_admm main_result.admm(j).error.depth];
    end
    
    [Err,I]=min(err_admm);
    
    main_result.admm_error=main_result.admm(I).error;
    main_result.admm_final=main_result.admm(I).final;
    
    index=sum(err_admm<err_base);
    
    fprintf('Baseline Error : %.4g ; ADMM error : %.4g ; Index : %d\n',err_base,Err,index);
    
    base=[base err_base];
    admm=[admm Err];
    
    fprintf('Baseline : depth %.4g ; normal %.4g ; albedo %.4g ; light %.4g ; matrix completion %.4g \n',...
    main_result.base_error.depth,main_result.base_error.normal,main_result.base_error.albedo,...
    main_result.base_error.light,main_result.base_error.M);

    fprintf('ADMM : depth %.4g ; normal %.4g ; albedo %.4g ; light %.4g ; matrix completion %.4g \n',...
    main_result.admm_error.depth,main_result.admm_error.normal,main_result.admm_error.albedo,...
    main_result.admm_error.light,main_result.admm_error.M(end));


    disp('...');

    %save([ list(i).name],'main_result','-v7.3');

end
%base=base(sort_id); admm=admm(sort_id);

plot(miss,base,'b*-',miss,admm,'r*-','LineWidth',4);
xlabel('Missing Data in percentage','FontSize', 20);
ylabel('Depth Reconstruction Error','FontSize', 20);
title('No Noise','FontSize', 20);
legend('Baseline','ADMM','FontSize', 20);

