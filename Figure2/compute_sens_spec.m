function s = compute_sens_spec(target_sd, p, n)
    % Computes required sensitivity when specificity = sensitivity
    % Inputs:
    %   target_sd - Target standard deviation in p (scalar)
    %   p         - True gene drive carrier frequency (scalar)
    %   n         - Sample size (scalar)
    %
    % Output:
    %   s - Required sensitivity (scalar)

    % Define the function to solve
    std_dev_function = @(s_guess) sqrt(((1 - s_guess) + (2 * s_guess - 1) * p) * ...
                          (s_guess - (2 * s_guess - 1) * p) / (n * (2 * s_guess - 1)^2)) - target_sd;
    
    % Sensitivity must be between 0.5 and 1 to keep (2s - 1) positive
    s_lower = 0.50001;  % Slightly above 0.5 to avoid division by zero
    s_upper = 1;

    % Ensure the function changes sign in the interval before calling fzero
    if std_dev_function(s_lower) * std_dev_function(s_upper) > 0
        error('No solution found for the given parameters.');
    end

    % Solve for sensitivity
    s = fzero(std_dev_function, [s_lower, s_upper]);
end

% function s = compute_sens_spec(SE_target, y_true, n)
% % compute_sens_spec: Computes sensitivity and specificity when they are equal,
% % given the target standard deviation, true prevalence, and sample size.
% %
% % Inputs:
% %   SE_target - Desired standard deviation (scalar)
% %   y_true    - True prevalence (scalar between 0 and 1)
% %   n         - Sample size (positive integer)
% %
% % Output:
% %   s         - Required Sensitivity = Specificity (scalar). 
% %                Returns NaN if no feasible solution.
% 
%     % Compute intermediate term
%     term = sqrt(y_true * (1 - y_true) / (n * SE_target^2));
% 
%     % Compute sensitivity/specificity
%     s = (1 + term) / 2;
% 
%     % Check feasibility
%     if s > 1
%         warning('No feasible solution: Required sensitivity/specificity exceeds 1.');
%         s = NaN;
%     elseif s < 0
%         warning('Computed sensitivity/specificity is negative, setting to NaN.');
%         s = NaN;
%     end
% end
