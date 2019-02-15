function cumDistrNGauss = calcCumDistrNGauss(param,abscissa,variableMean,...
    variableStd,logData,gaussParamIn,ratioTol)
%CALCCUMDISTRNGAUSS calculates the cumulative distribution of N Gaussians
%
%SYNOPSIS cumDistrNGauss = calcCumDistrNGauss(param,abscissa,variableMean,...
%     variableStd,logData,gaussParamIn)
%
%INPUT  param       : Vector of parameters indicating the means,
%                     stds and amplitudes of the N Gaussians.
%       abscissa    : Abscissa values at which the cumulative
%                     distribution is calculated.
%       variableMean: Flag with potential values:
%                     - 0 if assuming the fixed relationship
%                     (mean of nth Gaussian) = n * (mean of 1st Gaussian).
%                     - 1 if there is no relationship between the means of
%                     different Gaussians.
%                     - 2, 3, etc. if assuming the same fixed relationship
%                     as 0 but that the first detected Gaussian is actually
%                     the 2nd, 3rd, etc. Gaussian in the relationship.
%                     Optional. Default: 1.
%       variableStd : Flag with potential values:
%                     - 0 if assuming that all Gaussians have the same
%                     standard deviation. 
%                     - 1 if there is no relationship
%                     between the standard deviations of different
%                     Gaussians.
%                     - 2 if assuming the relationship
%                     (std of nth Gaussian) = sqrt(n) * (std of 1st Gaussian).
%                     This relationship is generalized if variableMean > 1.
%                     *** variableStd can equal 2 only if variableMean is
%                     not 1. ***
%                     - 3 if assuming the relationship
%                     (std of nth Gaussian) = n * (std of 1st Gaussian).
%                     This relationship is generalized if variableMean > 1.
%                     *** variableStd can equal 3 only if variableMean is
%                     not 1. ***
%                     Optional. Default: 1.
%       logData     : 1 for log normal data, where the log(data) is being
%                     fitted to a normal distribution, 0 otherwise. Note
%                     that data are passed to this function already after
%                     taking the log.
%                     Optional. Default: 0.
%       gaussParamIn: Matrix with number of rows equal to number of
%                     modes and two columns indicating the mean and
%                     standard deviation of each mode. If input, the
%                     specified mode parameters are used, and only the
%                     mode amplitudes are determined by data fitting. In
%                     this case, the input variableMean and variableStd
%                     are not used.
%                     Optional. Default: [].
%       ratioTol    : Tolerance for ratio between mean/std of 1st Gaussian
%                     and mean/std of subsequent Gaussians.
%                     If 0, ratio is taken strictly.
%                     If > 0, ratio is allowed to wiggle by ratioTol about
%                     the theoretial ratio. 
%                     Optional. Default: 0.
%                     Option currently implemented only for 3 cases:
%                     variableMean ~= 1 and variableStd = 1, 2 or 3.
%
%OUTPUT cumDistrNGauss: Values of the resulting cumulative distribution
%                       given the input abscissa values.
%
%Khuloud Jaqaman, August 2006; major updates in 2014 and 2015

%% Output

cumDistrNGauss = [];

%% Input

%check whether correct number of input arguments was used
if nargin < 2
    disp('--calcCumDistrNGauss: Incorrect number of input arguments!');
    return
end

if nargin < 3 || isempty(variableMean)
    variableMean = 1;
end

if nargin < 4 || isempty(variableStd)
    variableStd = 1;
end

if nargin < 5 || isempty(logData)
    logData = 0;
end
if logData && (variableMean==1&&variableStd~=1)
    disp('--fitHistWithGaussians: For log-normal fit, no current implementation for constrained std but variable mean. Exiting.')
    return
end

if nargin < 6 || isempty(gaussParamIn)
    gaussParamIn = [];
else
    variableMean = -1;
    variableStd = -1;
end

if nargin < 7 || isempty(ratioTol)
    ratioTol = 0;
end

%% Calculating the cumulative distribution

%get the means, standard deviations and amplitudes of the Gaussians from the input
%parameter vector
switch variableMean
    
    case -1 %if mean is given
        
        switch variableStd
            
            case -1 %if std is given
                
                %get number of Gaussians
                numGauss = length(param);
                
                %get their means, stds and amplitudes
                gaussMean = gaussParamIn(:,1);
                gaussStd  = gaussParamIn(:,2);
                gaussAmp  = param;
                
        end
        
    case 1 %if mean is variable

        switch variableStd

            case 0 %if std is constrained to all stds are equal

                %get number of Gaussians
                numGauss = floor(length(param)/2);

                %get their means, stds and amplitudes
                gaussMean = param(1:numGauss);
                gaussStd  = repmat(param(numGauss+1),numGauss,1);
                gaussAmp  = param(numGauss+2:end);

            case 1 %if std is variable

                %get number of Gaussians
                numGauss = length(param)/3;

                %get their means, stds and amplitudes
                gaussMean = param(1:numGauss);
                gaussStd  = param(numGauss+1:2*numGauss);
                gaussAmp  = param(2*numGauss+1:end);

        end %(switch variableStd)

    otherwise %if mean is constrained to mean_n = n * mean_1
        
        %get relationship of first fitted Gaussian to real first
        %Gaussian in series
        firstGauss = max(variableMean,1);

        switch variableStd
            
            case 0 %if std is constrained to all stds are equal
                
                %get number of Gaussians
                numGauss = length(param)-2;
                
                %get their means, stds and amplitudes
                if ~logData
                    tmpMean = param(1)/firstGauss;
                    gaussMean = (firstGauss:numGauss+firstGauss-1)' * tmpMean;
                    gaussStd  = repmat(param(2),numGauss,1);
                else
                    dataMean1 = exp(param(1)+param(2)^2/2);
                    dataMean1 = dataMean1 / firstGauss;
                    dataVar1 = exp(param(2)^2+2*param(1))*(exp(param(2)^2)-1);
                    dataMeanN = (firstGauss:numGauss+firstGauss-1)'*dataMean1;
                    dataVarN = repmat(dataVar1,numGauss,1);
                    gaussMean = log(dataMeanN.^2./sqrt(dataVarN+dataMeanN.^2));
                    gaussStd = sqrt(log(dataVarN./dataMeanN.^2+1));
                end
                gaussAmp  = param(3:end);
                
            case 1 %if std is variable
                
                if sum(ratioTol) == 0 %if strict ratio
                    
                    %get number of Gaussians
                    numGauss = floor(length(param)/2);
                    
                    %get their means, stds and amplitudes
                    if ~logData
                        tmpMean = param(1)/firstGauss;
                        gaussMean = (firstGauss:numGauss+firstGauss-1)' * tmpMean;
                    else
                        dataMean1= exp(param(1)+param(2)^2/2);
                        dataMean1 = dataMean1 / firstGauss;
                        dataMeanN = (firstGauss:numGaussT+firstGauss-1)' * dataMean1;
                        gaussMean = log(dataMeanN) - param(2:numGaussT+1).^2/2;
                    end
                    gaussStd  = param(2:numGauss+1);
                    gaussAmp  = param(numGauss+2:end);
                    
                else %if there is wiggle room
                
                    %get number of Gaussians
                    numGauss = length(param)/3;
                    
                    %get their means, stds and amplitudes
                    param = reshape(param,numGauss,3);
                    if ~logData
                        tmpMean = param(1,1) / firstMode;
                        gaussMean = [firstGauss; param(2:end,1)] * tmpMean;
                    else
                        dataMean1 = exp(param(1,1)+param(1,2)^2/2);
                        dataMean1 = dataMean1 / firstGauss;
                        dataMeanN = [firstGauss; param(2:end,1)] * dataMean1;
                        gaussMean = log(dataMeanN) - param(:,2).^2/2;
                    end
                    gaussStd  = param(:,2);
                    gaussAmp  = param(:,3);
                    
                end
                
            case 2 %if std is constrained to std_n = sqrt(n)*std_1
                
                if sum(ratioTol) == 0 %if strict ratio
                    
                    %get number of Gaussians
                    numGauss = length(param)-2;
                    
                    %get their means, stds and amplitudes
                    if ~logData
                        tmpMean = param(1)/firstGauss;
                        gaussMean = (firstGauss:numGauss+firstGauss-1)' * tmpMean;
                        tmpStd = param(2) / sqrt(firstGauss);
                        gaussStd  = sqrt(firstGauss:numGauss+firstGauss-1)' * tmpStd;
                    else
                        dataMean1 = exp(param(1)+param(2)^2/2);
                        dataMean1 = dataMean1 / firstGauss;
                        dataVar1 = exp(param(2)^2+2*param(1))*(exp(param(2)^2)-1);
                        dataVar1 = dataVar1 / firstGauss;
                        dataMeanN = (firstGauss:numGauss+firstGauss-1)' * dataMean1;
                        dataVarN = (firstGauss:numGauss+firstGauss-1)' * dataVar1;
                        gaussMean = log(dataMeanN.^2./sqrt(dataVarN+dataMeanN.^2));
                        gaussStd = sqrt(log(dataVarN./dataMeanN.^2+1));
                    end
                    gaussAmp  = param(3:end);
                    
                else %if there is wiggle room
                    
                    %get number of Gaussians
                    numGauss = length(param)/3;
                    
                    %get their means, stds and amplitudes
                    param = reshape(param,numGauss,3);
                    if ~logData
                        tmpMean = param(1,1) / firstGauss;
                        gaussMean = [firstGauss; param(2:end,1)] * tmpMean;
                        tmpStd = param(1,2) / sqrt(firstGauss);
                        gaussStd = sqrt([firstGauss; param(2:end,2)]) * tmpStd;
                    else
                        dataMean1 = exp(param(1,1)+param(1,2)^2/2);
                        dataMean1 = dataMean1 / firstGauss;
                        dataVar1 = exp(param(1,2)^2+2*param(1,1))*(exp(param(1,2)^2)-1);
                        dataVar1 = dataVar1 / firstGauss;
                        dataMeanN = [firstGauss; param(2:end,1)] * dataMean1;
                        dataVarN = [firstGauss; param(2:end,2)] * dataVar1;
                        gaussMean = log(dataMeanN.^2./sqrt(dataVarN+dataMeanN.^2));
                        gaussStd = sqrt(log(dataVarN./dataMeanN.^2+1));
                    end
                    gaussAmp = param(:,3);
                    
                end

            case 3 %if std is constrained to std_n = n*std_1
                
                if sum(ratioTol) == 0 %if strict ratio
                    
                    %get number of Gaussians
                    numGauss = length(param)-2;
                    
                    %get their means, stds and amplitudes
                    if ~logData
                        tmpMean = param(1)/firstGauss;
                        gaussMean = (firstGauss:numGauss+firstGauss-1)' * tmpMean;
                        tmpStd = param(2) / firstGauss;
                        gaussStd  = (firstGauss:numGauss+firstGauss-1)' * tmpStd;
                    else
                        dataMean1 = exp(param(1)+param(2)^2/2);
                        dataMean1 = dataMean1 / firstGauss;
                        dataVar1 = exp(param(2)^2+2*param(1))*(exp(param(2)^2)-1);
                        dataVar1 = dataVar1 / (firstGauss^2);
                        dataMeanN = (firstGauss:numGauss+firstGauss-1)' * dataMean1;
                        dataVarN = (firstGauss:numGauss+firstGauss-1)'.^2 * dataVar1;
                        gaussMean = log(dataMeanN.^2./sqrt(dataVarN+dataMeanN.^2));
                        gaussStd = sqrt(log(dataVarN./dataMeanN.^2+1));
                    end
                    gaussAmp  = param(3:end);
                    
                else %if there is wiggle room
                    
                    %get number of Gaussians
                    numGauss = length(param)/3;
                    
                    %get their means, stds and amplitudes
                    param = reshape(param,numGauss,3);
                    if ~logData
                        tmpMean = param(1,1) / firstGauss;
                        gaussMean = [firstGauss; param(2:end,1)] * tmpMean;
                        tmpStd = param(1,2) / firstGauss;
                        gaussStd = [firstGauss; param(2:end,2)] * tmpStd;
                    else
                        dataMean1 = exp(param(1,1)+param(1,2)^2/2);
                        dataMean1 = dataMean1 / firstGauss;
                        dataVar1 = exp(param(1,2)^2+2*param(1,1))*(exp(param(1,2)^2)-1);
                        dataVar1 = dataVar1 / (firstGauss^2);
                        dataMeanN = [firstGauss; param(2:end,1)] * dataMean1;
                        dataVarN = [firstGauss; param(2:end,2)].^2 * dataVar1;
                        gaussMean = log(dataMeanN.^2./sqrt(dataVarN+dataMeanN.^2));
                        gaussStd = sqrt(log(dataVarN./dataMeanN.^2+1));
                    end
                    gaussAmp = param(:,3);
                    
                end
                
        end %(switch variableStd)
        
end %(switch variableMean)

%calculate the cumulative distribution
cumDistrNGauss = zeros(size(abscissa));
for i=1:numGauss
    cumDistrNGauss = cumDistrNGauss + gaussAmp(i)*normcdf(abscissa,...
        gaussMean(i),gaussStd(i));
end

%% %%% ~~ the end ~~ %%%%%

