function retval = bramila_gPPI(response, source, task, hrf)
% bml_gppi(response, source, task) performs general linear model based
% gPPI regression analysis
%
% Inputs: 
%   response    = response time series (Nx1 matrix)
%   source      = source time series (Nx1 matrix)
%   task        = task time series (Nx1 matrix for M number of tasks)
%   (hrf)       = (OPTIONAL) hrf model to convolve with task
%
%   If the input parameter hrf is not provided, the task time series are 
%   expected to be in their final form, i.e. convolved with the approrpiate 
%   HRF if necessary, not in their raw block format
% 
%
% The return variable is a structure with the following fields:
%
% beta      = Coefficient estimates b, beta(4) is the gPPI interaction
% dfe       = Degrees of freedom for error
% sfit      = Estimated dispersion parameter
% s         = Theoretical or estimated dispersion parameter
% estdisp   = 0 when the 'estdisp' name-value pair argument value is 'off' 
%           and 1 when the 'estdisp' name-value pair argument value is 'on'
% covb      = Estimated covariance matrix for b
% se        = Vector of standard errors of the coefficient estimates b
% coeffcorr = Correlation matrix for b
% t         = t statistics for b
% p         = p-values for b
% resid     = Vector of residuals
% residp    = Vector of Pearson residuals
% residd    = Vector of deviance residuals
% resida    = Vector of Anscombe residuals
%
% 
% The coefficient estimates (retval.beta) are given in the following order:
%      
%       Estimate        Description
%       1               Baseline regressor (all 1s)
%       2               Context independent connectivity estimator (beta_0)
%       3 to M+2        Task activity estimators (beta_1 to beta_M)
%       M+3 to 2M+2     Context dependent estimators (beta_M+1 to beta_2M)
%
%       Where M = number of tasks

% Get the length of the time series (n) and the number of tasks (m)
[n,m] = size(task); %#ok<ASGLU>


% Check if hrf was provided as an input parameter
if nargin > 3
    % If so, then convolve task input with the hrf model
    for i = 1:m
        convtask(:,m) = conv(hrf,task(:,m)); %#ok<AGROW>
    end
    % Remove the excess data points from the end of the time series due to
    % convolution
    task = convtask(1:end-(length(hrf)-1),:);
end

% Get the binarized form of the task time series (all values above 0 get a
% value of 1, otherwise value = 0)
bintask = task>0;

% Generate the baseline regressor
%baseline = ones(n,1); Not required as glmfit does this automatically

% Generate context dependent regressors by convolving the source time
% series with the binarized task time series
contextdep = repmat(source,1,m).*bintask;

% Combine the regressors into one matrix
glminput = horzcat(source,task,contextdep);

% See result description at the beginning of this function
[~,~,retval] = glmfit(glminput,response);                         
                                       