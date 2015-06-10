%% Completes the matrix using the start values stored in the Solver
%%
%% params: maxErrIn,    spectralErrThreshold,   boundSpecErr, GT
%%
%% result:  iterations, ANew, 
%%          specErr,    distErr
%%          breakReason,
%%          LNew,       RNew,       ZNew,       ANew
function [result] = SigmaxSolver(M, params)
    global errors;  
    global Zs;
    Zs = [];
    errors = [0,0];
    myprintf('Starting matrix completion\n');        

    [params.Init.ANew, params.Init.MNew,params.Init.LNew,params.Init.RNew,params.Init.ZNew] = initialize(M,params);
    AInit = params.Init.ANew;
    MInit = params.Init.MNew;
    
    
    params.M = M;
    params.knownInds = find(M>0);
    
    [distErrNew, specErrNew] = getError(AInit, M, params.Init.LNew, params.Init.RNew, params.knownInds);
    myprintf('#0: Distance error = %.3g, spectral error %.3g\n', distErrNew, specErrNew);
    
    

    switch params.mode
        case 'direct'
            params.maxErrIn = 1e5*norm(getIndContents(M-MInit,params.knownInds));
            result = completeMatrixDirect(AInit, params);
        case 'augmented'
            %params.muInit = max(1e3* norm(getIndContents(M-MInit,params.knownInds)), params.muInit);
            result = completeMatrixAugmented(AInit, params);    
        otherwise
            error(strcat('SigmaxSolver doesn''t know this mode ', params.mode));
    end
%     [result.ANew, result.MNew,result.LNew,result.RNew,result.ZNew] = initialize(M,params);
    Vars.AEnc = result.ANew;
    Vars.LEnc = result.LNew;
    Vars.LambdaEnc = result.RNew;
    Vars.ZEnc = result.ZNew;
    [result.M, result.L, result.S, result.R, result.Z] = decode(Vars);
    
    result.errors = errors;
    result.Zs = Zs;
    
end
        

function [result] = completeMatrixDirect(AInit, params)
    
    
    [res] = completeMatrix(AInit, params);
    
    iter1 = res.iterations;
    % secone call, optimize data error
    if res.specErr >= params.spectralErrThreshold
        params.sigmaTailIn = res.specErr * 10000;
    else
        params.sigmaTailIn = params.boundSpecErr;
    end

    if ~strcmp(res.breakReason,'infeasable')
        params.maxErrIn = -1;
        res = completeMatrix(res.ANew,params);
        iter2 = res.iterations;
    end
    result = res;
    result.iterations = [iter1, iter2];
end

   
%% params:  muInit,     muInc
%%          muMaxIter,  lambdaMaxIter,  lambdaMaxValue
function [result] = completeMatrixAugmented(ACurr, params) 
    mu = params.muInit;
    result = struct();
    result.errors = zeros(1,2);
    result.ANew = ACurr;
    for i = 1:params.muMaxIter
        lambda = 0;
        for j = 1:params.lambdaMaxIter
            myprintf('Updating multipliers: lambda = %.2g, mu = %.2g\n\n',lambda, mu);
            params.muIn = mu;
            params.lambdaIn = lambda;
            resultOld = result;
            result = completeMatrix(result.ANew, params);
            [distErrNew, ~] = getError(result.ANew, params.M, result.LNew, result.RNew, params.knownInds);
            lambda = lambda + distErrNew * mu;
            if strcmp(result.breakReason,'threshold')
                return;
            end
            if lambda > params.lambdaMaxValue 
                break;
            end
        end
        mu = mu * params.muInc;
    end
end
        
        
        
%% Completes the matrix minimizing either sigmaTail or the dataError
%%
%% params:  M,          knownInds,  imgsize,
%%          mode = 'direct' || 'augmented',
%%          maxErrIn,   sigmaTailIn,
%%          muIn,       lambdaIn,
%%          mustImproveAfterIters,  maxIter
%% result:  iterations, ANew, 
%%          specErr,    distErr
%%          breakReason,
%%          LNew,       RNew,       ZNew,       ANew,       status
function [result] = completeMatrix(ACurr, params)
            
    % Initialization
    assignVars(params); % Make solver params directly accessible in this function

    prevErrs = inf(1, 2*mustImproveAfterIters);
    iter = 0;

    result.ANew = ACurr;
    result.iterations = 0;
    
    while true
        iter = iter +1;
        resultOld = result;
        [result] = optimize(result.ANew, params);
       
        
        result.iterations = iter;
        
        [distErrNew, specErrNew] = getError(result.ANew, params.M, result.LNew, result.RNew, params.knownInds);
        result.specErr = specErrNew;
        result.distErr = distErrNew;
        
        switch result.status
            case 0
                myprintf('#%d: Distance error = %.3g, spectral error %.3g\n', iter, distErrNew, specErrNew);
            case 1
                result.breakReason = 'infeasable';
                warning('Infeasable Problem. Stopping');
                break;
            case 9
                if iter > 1 
                    result.ANew = resultOld.ANew;
                    result.ZNew = resultOld.ZNew;
                    result.LNew = resultOld.LNew;
                    result.RNew = resultOld.RNew;
                end
                result.breakReason = 'stalling';
                break;
                if strcmp(params.mode,'direct')
                    if (sigmaTailIn == -1) 
                        myprintf('Notice: Increasing maxErrIn by *10.\n');
                        maxErrIn = maxErrIn*10;
                    else
                        myprintf('Notice: Increasing sigmaTailIn by *10.\n');
                        sigmaTailIn = sigmaTailIn*10;
                    end
                    if iter < maxIter 
                        continue;
                    else
                        break;
                    end
                elseif strcmp(params.mode,'augmented')
                    break;
                end
            otherwise
                result.breakReason = 'unknown';
                warning('Unknown problem encountered');
        end

        
        switch params.mode
            case 'direct'
                % If the spectral error reached the desired threshold, break.
                if maxErrIn >= 0 && specErrNew <= spectralErrThreshold 
                    result.breakReason = 'spectralErr';
                    myprintf('\nReached spectral error threshold. Breaking.\n\n');
                    break;
                end

                % If the data error reached the desired threshold, break.
                if sigmaTailIn >= 0 && distErrNew <= dataErrThreshold 
                    result.breakReason = 'distErr';
                    myprintf('\nReached distance error threshold. Breaking.\n\n');
                    break;
                end
            case 'augmented'
                if specErrNew <= spectralErrThreshold && distErrNew <= dataErrThreshold
                    result.breakReason = 'threshold';
                    myprintf('\nReached error threshold. Breaking.\n\n');
                    break;
                end
        end


        % If the last mustImproveAfterIters iterations didn't improve over the
        % previous ones, break.
        switch params.mode
            case 'direct'                
                if sigmaTailIn < 0
                    errNew = specErrNew;
                else
                    errNew = distErrNew;
                end
            case 'augmented'
                errNew = specErrNew + distErrNew;
        end
        prevErrs = [prevErrs(2:end) errNew];

        if min(prevErrs(mustImproveAfterIters+1:end)) >= (1-errDecayThreshold) * min(prevErrs(1:mustImproveAfterIters))
            result.breakReason = 'noImprovement';
            myprintf('\n Last %d iterations made no improvement. Breaking.\n\n', mustImproveAfterIters);
            break;
        end

        % If we've just finished running the last iteration, break.
        if iter >= maxIter
            result.breakReason = 'maxIter';
            myprintf('\nReached iteration limit. Breaking.\n\n');
            break;
        end

    end
end % function complete

        
        
%% Solves a convex optimization problem by looking only at 
%% the halfspace obtained by projection of the current A onto the 
%% 
%% params:  M,          knownInds,  imgsize,    
%%          mode = 'direct' || 'augmented',
%%          maxErrIn,   sigmaTailIn,
%%          muIn,       lambdaIn
%% result:  LNew,       RNew,       ZNew,       ANew,       status
function [result] = optimize(ACurr, params)

    assignVars(params);
    
    switch params.mode
        case 'direct'
            if sigmaTailIn < 0
                % objective sigma tail
                sigmaTail = sdpvar(1);
                Objective = sigmaTail;
                maxErr = maxErrIn;
                %myprintf('Optimize sigmaTails. maxErr= %d\n', maxErr);
            else
                % objective max error
                maxErr = sdpvar(1);
                Objective = maxErr;
                sigmaTail = sigmaTailIn;
                %myprintf('Optimize maxErr. sigmaTailBound = %d\n', sigmaTail);
            end
        case 'augmented'
            sigmaTail = sdpvar(1);
            maxErr = sdpvar(1);
            Objective = sigmaTail + muIn*0.5 * maxErr^2 + lambdaIn *maxErr;
    end

    


    r = 2;
    [lights, pixels] = size(M);
    m = lights+r; 
    n = pixels+r;
    rows = imgsize(1);
    cols = imgsize(2);


    %% Add all the optimization variables
    
    
    AOpt = sdpvar(m, n, 'full');            % Rank 2 matrix
    ROpt = sdpvar(1,pixels,'full');         % Albedo + Normalization vector
    L3Opt = sdpvar(lights,1,'full');        % Third row of the light vectors
    MOpt = sdpvar(m + n);                   % Symmetric, positive-semidefinite matrix to access the nuclear norm of A
    ZOpt = sdpvar(rows, cols, 'full');      % Optimization variable holding all the depth values

    %% Add all the constraints 
    Constr = [];
    Constr = [Constr, ROpt(:) <= -1];                % Avoid the trivial solution            
    Constr = [Constr, AOpt(1:r, 1:r) == eye(r)];    % Ensure identity matrix in AOpt

    %% Half space constraint
    % substituting nuclear norm constraint with a trace (Fazel's thesis)
    
    Constr = [Constr, MOpt >= 0];                   % positive semi definite
    %Constr = [Constr, MOpt(n+1:end, 1:n) == AOpt];  % coherence with A
    Constr = [Constr, MOpt(1:m,m+1:end) == AOpt];  % coherence with A

    [Ua, ~, Va] = svd(ACurr); 
    Wa = zeros(m, n); 
    Wa(1:r, 1:r) = eye(r);
    Constr = [Constr, 0.5 * trace(MOpt) - trace((Ua*Wa*Va')'*AOpt) <= sigmaTail];    % half space constraint

    %% Data Error constraint
    Constr = [Constr, cone(getIndContents( AOpt(r+1:end,r+1:end) - L3Opt * ones(1,pixels) -M * diag(ROpt), knownInds) ,maxErr)]; % proximity to known data constraint

    %% Integrability constraint

    
    % x derivatives
    for j=1:cols
        for i=1:rows-1
            ind = sub2ind([rows, cols],i,j);
            Constr = [Constr, ZOpt(i+1,j) - ZOpt(i,j) == AOpt(1,r+ind)];
        end                
    end
    % y derivatives
    for i=1:rows
        for j=1:cols-1
            ind = sub2ind([rows, cols],i,j);
            Constr = [Constr, ZOpt(i,j+1) - ZOpt(i,j) == AOpt(2,r+ind)];
        end                
    end

    
    %% Calling Yalmip to solve the optimization problem
    if isfield(params,'GT') 
        GT = params.GT;
        assign(AOpt,params.Init.ANew);
        assign(L3Opt,params.Init.LNew);
        assign(ROpt,params.Init.RNew');
        assign(ZOpt,params.Init.ZNew);    
        [U,S,V] = svd(params.Init.ANew,'econ');
        M = [U*S*U',params.Init.ANew; params.Init.ANew', V*S*V'];
        assign(MOpt,M);
        if strcmp(params.mode,'augmented')
            assign(sigmaTail,1e-10);
            assign(maxErr,1e-10);
        elseif sigmaTailIn == -1
            assign(sigmaTail,1e-10);
        else 
            assign(maxErr,1e-10);
        end
        %checkset(Constr);
        
    end
    
    %myprintf('Calling yalmip\n');
    options = sdpsettings('cachesolvers',1, 'verbose', 0, 'debug', 0, 'solver', 'mosek');
    sol = solvesdp(Constr,Objective,options);

    %checkset(Constr);
    
    result = struct();
    result.status = sol.problem;
    
    if (sol.problem ~= 0)
        sol.info
        yalmiperror(sol.problem);
        %warning('YALMIP error #%d: Something went wrong in the optimization process!\n Returning previous results', sol.problem);               
        result.ANew = ACurr;
        result.ZNew = double(ZOpt);
        result.LNew = double(L3Opt);
        result.RNew = double(ROpt);
    else
        result.ANew = double(AOpt);
        result.ZNew = double(ZOpt);
        result.LNew = double(L3Opt);
        result.RNew = double(ROpt);
    end
    
    % Evaluate depth (after GBR) with respect to GT
    global Zs;
    if isfield(params,'ZGT')
        [~, ZRes] = getGBR(params.ZGT,result.ZNew);
        err = norm(ZRes-params.Init.ZNew,'fro');
        Zs = [Zs, err];
    end

end % function optimize
        
 
   
        

    
%% Initialization of the algorithm
function [AInit,MInit, LInit, LambdaInit,ZInit] = initialize(M,params)
    r = 2;
    imgsize = params.imgsize;
    pixels = prod(imgsize);
    lights = size(M,1);
    
    if isfield(params,'GT')
        %% Initialize with ground truth
        fprintf('Using ground truth\n');
        Vars = params.GT;
        sanitytests(Vars);
    else
        %% Initialize without ground truth
        
        switch params.initWith
            case 'random'
                L = rand(lights,3);
                S = rand(pixels,3);
            case 'M'
                knownInds = M > 0;
                fprintf('Completing matrix\n');
                try 
                    [L, S] = damped_wiberg_new(M, knownInds, 3);
                catch 
                    [~,L,S] = l2_rpca_mask_alm_fast(M,knownInds, 3);
                end
        end
              
        Vars = struct();    
        Vars.M = L*S';
        Vars.L = L;
        Vars.S = normr(S); 
        Vars.R = sqrt(sum(S'.^2))';
        Vars.Z = zeros(imgsize);
        
%         Sol=SimpleSolver(M,imgsize,3, struct());
%         Sol.solve();
% 
%         [Vars.M, Vars.L, Vars.S, Vars.R, Z] = Sol.extract();
%         Vars.Z = reshape(Z,imgsize);


    end
    
    [AInit, MInit,LInit,LambdaInit,ZInit] = encode(Vars);
end
        
function [AEnc, MEnc,LEnc,LambdaEnc,ZEnc] = encode(Vars)
    lights = size(Vars.M,1);    
    pixels = size(Vars.M,2);
    
    if size(Vars.S,1) ==3 
        Vars.S = Vars.S';
    end
        
    r = 2;
    LEnc = Vars.L(:,3);
    
    NEnc = [-diag(1./Vars.S(:,3)) * Vars.S(:,1:r), -ones(pixels,1)];
    LambdaEnc = -sqrt(sum(NEnc'.^2))' ./ Vars.R;
    AEnc = zeros(size(Vars.M)+2);
    AEnc(1:r, 1:r) = eye(r);
    AEnc(r+1:end,r+1:end) = Vars.M * diag(LambdaEnc) + Vars.L(:,3)*ones(1,pixels);
    AEnc(1:r,r+1:end) = NEnc(:,1:r)';
    AEnc(r+1:end,1:r) = Vars.L(:,1:r);
    ZEnc = Vars.Z;    
    MEnc = Vars.M;
    
    
        
end

function [MDec,LDec,SDec,RDec,ZDec] = decode(Vars)
    r = 2;
    lights = size(Vars.AEnc,1)-r;
    pixels = size(Vars.AEnc,2)-r;
    LDec = [Vars.AEnc(r+1:end,1:r),Vars.LEnc];
    MDec = (Vars.AEnc(r+1:end,r+1:end) - Vars.LEnc*ones(1,pixels)) / diag(Vars.LambdaEnc) ;
    SDec = LDec \ MDec;
    RDec = sqrt(sum(SDec.^2))';
    SDec = normc(SDec)';
    %SDec = SDec';
    ZDec = Vars.ZEnc;
    
end

function [] = sanitytests(Vars)
fprintf('Running sanitytests\n');
r = 2;
eeppss = 1e-10;

[AEnc, MEnc,LEnc,LambdaEnc,ZEnc] = encode(Vars);

%% test rank of A

if (rank(AEnc) > r)
    warning('Rank of AInit is greater than 2!');
    svd(AEnc);
end

%% test factorization of A
S = AEnc(1:r,r+1:end);
L = AEnc(r+1:end,1:r);
M = AEnc(r+1:end,r+1:end);
if norm(L*S - M,'fro') > eeppss
    warning('Not a consistent factorization in A');
    L*S-M
end


%% test derivatives of A
Zx = diff(ZEnc,1);
Zy = diff(ZEnc,1,2);

Sx = reshape(AEnc(1,r+1:end),size(ZEnc));
Sx = Sx(1:end-1,:);
Sy = reshape(AEnc(2,r+1:end), size(ZEnc));
Sy = Sy(:,1:end-1);

if norm(Zx - Sx,'fro') > eeppss || norm(Zy - Sy,'fro') > eeppss
    warning('Something with the derivatives in A is not correct');
    Zx
    Sx
    Zx-Sx
end



end

        
      
%% Calculate the current data and specular error
function [distErr, specErr] = getError(A, M, L, R, knownInds)
    r = 2;
    [lights, pixels] = size(M);
    distErr = norm(getIndContents(A(r+1:end,r+1:end) - L * ones(1,pixels) - M *diag(R), knownInds), 'fro');
    s = svd(A);
    specErr = norm(s(r+1:end))/norm(s);
    
    global errors;
    errors = [errors; distErr, specErr];
end



