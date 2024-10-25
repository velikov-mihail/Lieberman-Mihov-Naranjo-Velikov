function [res] = runGMM(Factors, Rf, TestAssets, HAC_Lags)  

    % Store some dimensions
    [nMonths, nTestAssets] = size(TestAssets);
    nFactors = size(Factors, 2);

    % Subtract risk-free rate from test asset returns
    rptdRF = repmat(Rf, 1, nTestAssets);
    Re = TestAssets - rptdRF;

    % Store the long/short return
    lsRet = mean(Re(:,nTestAssets)-Re(:,1), 1, 'omitnan');
    
    % Start with identity weighting matrix
    W0 = eye(nTestAssets+nFactors);

    % Initial guesses for the parameters
    initialRiskPrices = zeros(1, nFactors); 
    initialFactorMeans = zeros(1, nFactors);
    p0 = [initialRiskPrices initialFactorMeans]';

    % Optimization options for fminunc
    options = optimoptions(@fminunc, 'Algorithm', 'quasi-newton', ...
                                     'Display', 'off', ...
                                     'MaxIterations', 20000);
    
    % Minimize the GMM objective function
    [p1, ~] = fminunc(@(params) gmm_objective(params, Re, Factors, W0), p0, options);
    
    W_opt1 = gmm_weight(Re, Factors, p1, HAC_Lags);

    [theta, ~] = fminunc(@(params) gmm_objective(params, Re, Factors, W_opt1), p1, options);
    
    W = gmm_weight(Re, Factors, theta, HAC_Lags);
    
    res.theta = theta;
    riskPrices = theta(1:nFactors);
    factorMeans = theta(nFactors+1:end);

    % Get the gradient and parameter standard errors
    [Moment_V, se] = getParamsSE(Re, Factors, factorMeans, riskPrices, W);

    % Store the standard errors and t-stats
    res.se = se;
    res.t = res.theta./res.se;

    % Get the pricing errors
    [alphas, talphas] = getAlphas(Re, Factors, riskPrices, factorMeans, Moment_V);

    % Get the mean absolute error
    MAE = mean(abs(alphas), 'omitnan');

    % Correct the L/S standard error
    R = zeros(1, cols(Moment_V));
    R(1,1) = sign(-lsRet);
    R(1,nTestAssets) = sign(lsRet);
    temp_VAR = R*Moment_V*R';
    alphas(nTestAssets+1) = (alphas(nTestAssets)-alphas(1)) * sign(lsRet);
    talphas(nTestAssets+1) = alphas(nTestAssets+1)./sqrt(temp_VAR/nMonths);
    
    betas = NaN(nTestAssets+1, 1);
    tbetas = nan(nTestAssets+1, 1); 
    for c=1:nTestAssets+1
        if c<=nTestAssets
            tempRet = Re(:,c);
        else
            tempRet = sign(lsRet)*(Re(:,nTestAssets) - Re(:,1));
        end
        [mu_COV_Rf2, SE_COV_Rf2] = GMM_Covariance(tempRet, Factors(:,end), theta(end));
        betas(c) = mu_COV_Rf2(end);
        tbetas(c) = mu_COV_Rf2(end)/SE_COV_Rf2(end);
    end
        
    res.MAE     = 100*MAE;
    res.alphas = alphas;
    res.talphas = talphas;
    res.betas = betas;
    res.tbetas = tbetas;
end

function [Moment_V, se] = getParamsSE(Re, Factors, factorMeans, riskPrices, W)
    % Initialize the sizes of the inputs
    [nMonths, nTestAssets] = size(Re); % Get the number of months and test assets
    nFactors = size(Factors, 2); % Number of factors
    nParams = 2 * nFactors; % Total number of parameters (2 for each factor)
    
    % Preallocate the gradient matrix G
    G = zeros(nTestAssets+nFactors, nParams, nMonths);
    
    % Constructing the gradient for first moment condition
    for i = 1:nTestAssets
        % Calculate the partial derivative w.r.t. the risk price   
        % d/db E_T[Re - Re(f'-Ef')b] = E_T[-Re(f'-Ef')]
        for lambda = 1:nFactors
            for t = 1:nMonths        
                G(i, lambda, t) = -Re(t, i) * (Factors(t, lambda) - factorMeans(lambda));
            end 
        end 
        
        % Calculate the partial derivative w.r.t. the factor mean
        % d/dEf E_T[Re-Re(f'-Ef')b] = E_T[Re*b]
        for mu = 1:nFactors
            for t = 1:nMonths
                G(i, nFactors+mu, t) = Re(t, i) * riskPrices(mu);
            end 
        end 
    end 
    
    % Constructing the gradient for the second moment condition
    for mu = 1:size(Factors, 2)
        for t = 1:nMonths
            % Calculate the partial derivative w.r.t. the factor mean
            % d/dEf E_T[f'-Ef'] = -1
            G(nTestAssets+mu, nFactors+mu, t) = -1;
        end
    end 
    
    % Average the gradients over all months
    G = mean(G, 3, 'omitnan');
    
    % Calculate the variance-covariance matrix V 
    V = inv(G' * W * G); 

    % Calculate the standard errors 
    se = sqrt(diag(V) / nMonths);

    % Calculate the moment variance covariance matrix
    Moment_V = inv(W)-G/(G'*W*G)*G';
    
end


function [alphas, talphas] = getAlphas(Re, Factors, riskPrices, factorMeans, Moment_V)

    [nMonths, nTestAssets] = size(Re);
        
    PriceError = NaN(nTestAssets, nMonths);
    for t=1:nMonths
        for i=1:nTestAssets
            PriceError(i,t) = Re(t, i) - Re(t, i) * ... 
                                 (Factors(t, :) - factorMeans') * riskPrices;
        end
    end  

    % Get the alphas their t-stats
    alphas = [mean(PriceError, 2, 'omitnan')];
    talphas = alphas./sqrt(diag(Moment_V(1:nTestAssets,1:nTestAssets))/nMonths);
end



function W = gmm_weight(Re, Factors, theta,lags)
% This function calculates the optimal GMM weighting matrix.
% Moment_Cond: The moment conditions for each observation.
% theta: Parameter estimates from the first step.
% lags: Number of lags for Newey-West standard errors.

    
    [nMonths, nTestAssets] = size(Re); 
    nFactors = size(Factors,2); 
    
    riskPrices = theta(1:nFactors);
    factorMeans = theta(nFactors+1:end);


    % Preallocate GAMMA for speed.
    GAMMA = zeros(nTestAssets+nFactors, nTestAssets+nFactors, lags+1);
        
    fTilde = Factors - repmat(factorMeans', nMonths, 1);
    u = [Re' - Re' .* repmat((fTilde*riskPrices)', nTestAssets, 1);
         fTilde'];

    lagU = lag(u', 1, 0)';

    % Compute the mean of f across all time periods.    
    f = [mean(u, 2) mean(lagU, 2)];

    % Calculate the autocovariance matrices GAMMA.
    for k=0:lags
        V = zeros(nTestAssets+nFactors, nTestAssets+nFactors, nMonths);
        for i=1:nTestAssets
            for j=1:nTestAssets
                for t=1+k:nMonths
                    V(i,j,t) = ...
                        ((Re(t, i) -...
                        Re(t, i)*(Factors(t, :)-factorMeans')*riskPrices)-f(i,1+k))*...
                        ((Re(t-k, j)-...
                        Re(t-k, j)*(Factors(t-k, :)-factorMeans')*riskPrices)-f(j,1+k));     
                end
            end 
            for fact = 1:nFactors
                for t=1+k:nMonths
                    V(i,nTestAssets+fact,t) = ...
                        ((Re(t, i)-...
                        Re(t, i)*(Factors(t, :) - factorMeans')*riskPrices) - f(i,1+k))*...
                        ((Factors(t-k, fact)-factorMeans(fact))-f(nFactors+fact));
                    V(nTestAssets+fact,i,t) = ...
                        ((Re(t, i)-...
                        Re(t, i)*(Factors(t, :) - factorMeans')*riskPrices) - f(i,1+k))*...
                        ((Factors(t-k, fact) - factorMeans(fact)) - f(nFactors+fact));
                end
            end
        end 
        for fact1 = 1:nFactors
            for fact2 = 1:nFactors
                for t=1+k:nMonths
                    V(nTestAssets+fact1,nTestAssets+fact2,t) = ...
                        ((Factors(t, fact1) - factorMeans(fact1)) - f(nFactors+fact1,1+k))*...
                        ((Factors(t, fact2) - factorMeans(fact2)) - f(nFactors+fact2,1+k));
                end
            end
        end 
        
        % Average V over T to get GAMMA for each lag.
        GAMMA(:,:,k+1) = mean(V, 3, 'omitnan');
    end

    % Calculate weights 
    w = zeros(lags+1,1);
    for j=0:lags
        w(j+1,1)=1-(j/(lags+1));
    end

    % Compute the spectral density matrix S.
    S = GAMMA(:,:,1);
    for lag_num = 1:lags
        % Add weighted autocovariances.
        S = S + w(lag_num+1)*(GAMMA(:,:,lag_num+1)+GAMMA(:,:,lag_num+1)');
    end

    % The optimal weighting matrix W is the inverse of S.    
    W = inv(S);
end

function [mu, SE, V] = GMM_Covariance(Y, X, muX)
    % Calculate covariance and standard errors for GMM estimation
    % Inputs:
    % X - Vector of returns for the factor
    % Y - Vector of returns for the test asset
    % muX - Mean return of the factor X

    % Number of observations
    n = size(Y, 1);
    
    % Calculate mean of Y
    muY = mean(Y);

    % Calculate covariance between X and Y
    % sigma_XY = (1/(n-1))*sum(X.*Y - muX*muY);
    sigma_XY = (1/(n-1))*sum((X - muX) .* (Y - muY));
    
    % Form matrix f of demeaned X and product of demeaned X and Y
    f = [(Y - muY) (X.*Y - muX*muY)];
    
    % Define matrix d for linear transformation
    d = [-1 0; -muX -1];
    
    % Calculate sample covariance matrix S of f
    S = (1/n)*(f' * f);
    
    % Invert S through pre-multiplication and post-multiplication by d and its transpose
    % to get covariance matrix V of the parameter estimates
    V = inv(d' * inv(S) * d);
    
    % Standard errors for mean of X and covariance between X and Y
    SE_muY = sqrt(V(1,1)/n);
    SE_sigmaXY = sqrt(V(2,2)/n);
    
    % Package the means and standard errors into output variables
    mu = [muY; sigma_XY];
    SE = [SE_muY; SE_sigmaXY];
end



function moments = moment_conditions(params, Re, f)               
    %  This function calculates moment conditions following Cochrane
    %  (2005), Chapter 13.2, pg. 257: 
    %  E_T[Re - Re(f'-Ef')b]
    %  E_T(f-Ef)
    
    [T, ~] = size(Re);
    nFactors = size(f, 2);

    riskPrices = params(1:nFactors)';
    factorMeans = params(nFactors+1:end)';

    % Calculate the expected excess returns
    rptdFactorMeans = repmat(factorMeans, T, 1);
    expected_excess_return = Re .* ((f - rptdFactorMeans) * riskPrices');
    
    % Calculate deviations for each asset at each time
    deviations = Re - expected_excess_return;
    
    % Time-series average of these deviations for each asset 
    moments_assets = mean(deviations);
    
    % Additional moment condition for the average of the factors
    moments_factors = mean(f) - factorMeans;
    
    % Combine the moments from the assets and the factors
    moments = [moments_assets, moments_factors];
end

function J = gmm_objective(params, Re, f, W)
    % Define the GMM objective function
    moments = moment_conditions(params, Re, f);
    
    % Quadratic form of the moments
    J = moments * W * moments'; 
end



