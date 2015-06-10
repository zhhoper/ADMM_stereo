%% Default params for Sigmax
function params = SigmaxDefaultParams(imgsize,mode)
    params = struct();          
    
    
    
    params.mustImproveAfterIters = 20;
    params.errDecayThreshold = 0.002;
    params.spectralErrThreshold = 1e-5;
    params.dataErrThreshold = 1e-5;
    params.initWith = 'M';
    
    switch mode 
        case 'direct'      
            params.maxIter = 200;              
            params.maxErrIn = 10* prod(imgsize);
            params.sigmaTailIn = -1;  
            params.boundSpecErr = 1e-2;
            params.mode = 'direct';
        case 'direct-fast'
            params.maxIter = 50;        
            params.mustImproveAfterIters = 5;
            params.maxErrIn = 10* prod(imgsize);
            params.sigmaTailIn = -1;  
            params.boundSpecErr = 1e-2;
            params.mode = 'direct';
        case 'augmented1-10'
            params.maxIter = 50;
            params.mustImproveAfterIters = 10;
            params.muInit = 1;
            params.muInc = 10;
            params.muMaxIter = 6;
            params.lambdaMaxIter = 6;
            params.lambdaMaxValue = 1e6;
            params.mode = 'augmented';
        case 'augmented1-10-random'
            params.maxIter = 50;
            params.mustImproveAfterIters = 10;
            params.muInit = 1;
            params.muInc = 10;
            params.muMaxIter = 6;
            params.lambdaMaxIter = 6;
            params.lambdaMaxValue = 1e6;
            params.initWith='random';
            params.mode = 'augmented';
        case 'augmented10-2'
            params.spectralErrThreshold = 1e-5;
            params.dataErrThreshold = 1e-4;
            params.maxIter = 50;
            params.mustImproveAfterIters = 10;
            params.muInit = 10;
            params.muInc = 2;
            params.muMaxIter = 6;
            params.lambdaMaxIter = 6;
            params.lambdaMaxValue = 1e6;
            params.initWith = 'random';
            params.mode = 'augmented';
        case 'augmented10-2-random'
            params.maxIter = 50;
            params.mustImproveAfterIters = 5;
            params.muInit = 10;
            params.muInc = 2;
            params.muMaxIter = 6;
            params.lambdaMaxIter = 3;
            params.lambdaMaxValue = 1e6;
            params.initWith = 'random';
            params.mode = 'augmented';
        otherwise
            error('Unknown mode...');
    end
    params.imgsize = imgsize;
%     params.mode = mode;
end
