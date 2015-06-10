% simulate

close all; clear all; clear;
addpath('../solvers'); addpath('../functions'); addpath('../external'); addpath('../simulate');

know=1;



%seed = round(rand()*1e6);
%% Test Case 1 : no Noise, no missing data
seed = 58597;
rng(seed);

size = [40,40];
noiseStrength = 0.0;
knowledgeRatio = know;  % change this accordingly
fprintf('Processing : noise %d, missing data %d \n',noiseStrength,1-knowledgeRatio);
rows = size(1);
cols = size(2);
lights = 20;
pixels = prod(size);

scenario.object = 'new1';          % random || peaks     ||  ellipse
scenario.lighting = 'random';       % random    ||  fitted
scenario.albedo = 'random';         % number || 'random'


Sim = Simulator(lights,size,scenario,noiseStrength, knowledgeRatio);
Sim.simulate();

showNormals(Sim.S, Sim.Z, 'Ground Truth');
f = showImages(Sim.MNoise, size);
set(f,'Name','Input Images');

%load a file 

load(['/fs/jacobsvideo/Photo-stereo/ADMM_stereo_0507/ADMM/result_new/1Result_know_' num2str(know) '_noise_0.mat']);

    err_base=main_result.base_error.depth;
    
    
    err_admm=[];
    for j=1:length(main_result.admm)        
        err_admm=[err_admm main_result.admm(j).error.depth];
    end
    
    [Err,I]=min(err_admm);
    
    % processing depth %
    
    Depth_base=100*err_base*rows*cols/norm(Sim.Z);
    Depth_admm=100*Err*rows*cols/norm(Sim.Z);
    
    figure;
    subplot(1,2,1)
    imagesc(abs(Sim.Z-main_result.base_final.depth));
    title('Depth Reconstruction Error for Baseline','FontSize',15);
    h = colorbar; set(h, 'ylim')
    
    subplot(1,2,2)
    imagesc(abs(Sim.Z-main_result.admm(I).final.depth));
    title('Depth Reconstruction Error for ADMM','FontSize',15);
    h = colorbar; set(h, 'ylim')
    
    
    
    fprintf(' Depth Reconstruction Error : Baseline = %d ; ADMM = %d \n',Depth_base,Depth_admm);
    
    % Processing Albedo %
    
    Alb_base=100*main_result.base_error.albedo*rows*cols/norm(Sim.R);
    Alb_admm=100*main_result.admm(I).error.albedo*rows*cols/norm(Sim.R);
    
    fprintf(' Albedo Reconstruction Error : Baseline = %d ; ADMM = %d \n',Alb_base,Alb_admm);
    
    figure;
    subplot(1,2,1)
    
    imagesc(reshape(abs(Sim.R-main_result.base_final.albedo),[Sim.resolutionX,Sim.resolutionY]));
    title('Albedo Reconstruction Error for Baseline','FontSize',15);
    colorbar
    
    subplot(1,2,2)    
    imagesc(reshape(abs(Sim.R-main_result.admm(I).final.albedo),[Sim.resolutionX,Sim.resolutionY]));
    title('Albedo Reconstruction Error for ADMM','FontSize',15);
    colorbar
    
    
    % Process Surface normal %
    
    Surf_base=main_result.base_final.normal;
    Surf_admm=main_result.admm(I).final.normal;
    
    S_base_err=mean(abs(acos(sum(Surf_base.*Sim.S,2))));
    S_admm_err=mean(abs(acos(sum(Surf_admm.*Sim.S,2))));
    
    fprintf(' Surface Normal Reconstruction Error (angle) : Baseline = %.3d ; ADMM = %.3d \n',rad2deg(S_base_err),rad2deg(S_admm_err));
    
    figure;
subplot(3,3,1),imagesc((reshape(Sim.S(:,1), [rows, cols])));
subplot(3,3,2),imagesc((reshape(Surf_base(:,1), [rows, cols])));
subplot(3,3,3),imagesc((reshape(Surf_admm(:,1), [rows, cols])));   

subplot(3,3,4),imagesc((reshape(Sim.S(:,2), [rows, cols])));
subplot(3,3,5),imagesc((reshape(Surf_base(:,2), [rows, cols])));
subplot(3,3,6),imagesc((reshape(Surf_admm(:,2), [rows, cols])));

subplot(3,3,7),imagesc((reshape(Sim.S(:,3), [rows, cols])));
subplot(3,3,8),imagesc((reshape(Surf_base(:,3), [rows, cols])));
subplot(3,3,9),imagesc((reshape(Surf_admm(:,3), [rows, cols])));

    
    % Process Light %
 
    
    Light_base=main_result.base_final.light;
    Light_admm=main_result.admm(I).final.light;
    
    n_Light_base=Light_base./repmat(sqrt(sum(Light_base.^2,2)),[1,3]);
    n_Light_admm=Light_admm./repmat(sqrt(sum(Light_admm.^2,2)),[1,3]);
    n_Light_L=Sim.L./repmat(sqrt(sum(Sim.L.^2,2)),[1,3]);
    
    L_base_err=mean(abs(acos(sum(n_Light_base.*n_Light_L,2))));
    L_admm_err=mean(abs(acos(sum(n_Light_admm.*n_Light_L,2))));
    
    fprintf(' Light Reconstruction Error (angle) : Baseline = %.3d ; ADMM = %.3d \n',rad2deg(L_base_err),rad2deg(L_admm_err));
    
    h = showLights(Sim.L, Light_admm,Light_base);
    
    
    % Convergence Plots %
    
    figure;
subplot(2,2,1)
semilogy(main_result.admm(I).out.obj,'LineWidth',2.5); title('Convergence plot of overall objective function'); xlabel('Iterations'); ylabel('Objective Value');
subplot(2,2,2)
semilogy(main_result.admm(I).out.obj1,'LineWidth',2.5); title('Convergence plot of Data consistency term'); xlabel('Iterations'); ylabel('Objective Value');
subplot(2,2,3)
semilogy(main_result.admm(I).out.obj2,'LineWidth',2.5); title('Convergence plot of Integrability Constraint term'); xlabel('Iterations'); ylabel('Objective Value');
subplot(2,2,4)
semilogy(main_result.admm(I).out.obj3,'LineWidth',2.5); title('Covergence plot of Constraint violation'); xlabel('Iterations'); ylabel('Objective Value');
    
