clear;  close all; clc;
addpath('../solvers');
addpath('../functions');
addpath('../external');
addpath('../simulate')

    
noise=0;
know=0.5:-0.1:0.2;

parameter=combvec(noise,know);

for i=1:size(parameter,2)
    
    Sim = generate_data(parameter(1,i),parameter(2,i));
    
    % initialzie data
    [X, Y, lam, l3, Z, lambda, ini_SGRB,e,r] = initialize_ADMM(Sim);
    
    main_result.base_error=e; main_result.base_final=r;
    main_result.sim=Sim;

    var = struct;
    var.M = Sim.MNoise;
    var.W = Sim.INC;
    var.X = X;
    var.Y = Y;
    var.lam = lam;
    var.l3 = l3;
    var.Z = Z;
    var.lambda = lambda;
    [var.img_row, var.img_col] = size(Sim.Z);

% what should be the correct parameter?
    
    C2=[1 10]; TAU=[100];

    para_exp=combvec(C2,TAU);
    
    
    tim=0;
    for t=1:size(para_exp,2)
    
        tic
        var.c2 = para_exp(1,t);
        var.tau = para_exp(2,t);

        result = ADMM_opt(var);

        numPixel = var.img_row*var.img_col;

        depth = result.var.Z;
        light = -[result.var.X(3:end, 1:2), result.var.l3];
        %normal = [-result.var.X(1:2, 3:end)', ones(numPixel,1)];
        normal = [result.var.X(1:2, 3:end)', -ones(numPixel,1)];
        albedo = -(sqrt(sum(normal.^2, 2)))./result.var.lam;
        normal = my_normr(normal);

        %depth=reconstructDepthMap_adapted(-normal',[var.img_row,var.img_col]);

        [LGBR,ZGBR,SGBR,RGBR,G] = applyGBR(Sim.Z, Sim.R, depth, light, normal, albedo);

        
        
        
        
        result.final.light=LGBR;
        result.final.depth=ZGBR;
        result.final.normal=SGBR;
        result.final.albedo=RGBR;
        
        
        Mask=ones(var.img_row, var.img_col); %Mask(2:var.img_row-1, 2:var.img_col-1)=1; %Mask that deletes boundary
        Mask1=Mask(:);

        recZ=reshape(ZGBR,[var.img_row, var.img_col]);
        result.error.depth = norm((Sim.Z - recZ),'fro')/numPixel;
        result.error.light = norm(Sim.L - LGBR,'fro');
        result.error.normal = norm((Sim.S - SGBR).*repmat(Mask1,[1,3]),'fro')/numPixel;
        result.error.albedo = norm((Sim.R - RGBR).*Mask1)/numPixel;
        result.error.M=result.out.obj1;
    
        time=toc; tim=time+tim;
        disp(['Done Run ' num2str(i) ' parameter : ' num2str(t) ' Time : ' num2str(time)]);
        
        main_result.admm(t)=result;
        
        
        save(['result_new_3/Result_know_' num2str(parameter(2,i)) '_noise_' num2str(parameter(1,i)) '.mat'],'main_result','-v7.3');
        
    end
    
    %save main result %
    
    
    save(['result_new_3/Result_know_' num2str(parameter(2,i)) '_noise_' num2str(parameter(1,i)) '.mat'],'main_result','-v7.3');
    
    disp(['Total Time : ' num2str(tim)]);
    
    
    main_result=[];
end
    