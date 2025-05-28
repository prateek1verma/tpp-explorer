function s = compute_sensitivity(target_sd, p, n, specificity)
    % Computes the required sensitivity 's' to achieve a target standard deviation.
    % 
    % Inputs:
    %   target_sd   - Target standard deviation (scalar)
    %   p           - True gene drive carrier frequency (scalar, between 0 and 1)
    %   n           - Sample size (scalar, positive integer)
    %   specificity - Test specificity per mosquito (scalar, between 0 and 1)
    %
    % Output:
    %   s - Required sensitivity (scalar between specificity and 1)

    % Calculate false positive rate
    f = 1 - specificity;

    % Define the standard error function in terms of sensitivity
    std_dev_function = @(s_guess) sqrt((f + (s_guess - f) * p) * ...
                          (1 - f - (s_guess - f) * p) / (n * (s_guess - f)^2)) - target_sd;

    % Define the search interval for sensitivity (must be greater than specificity)
    s_lower = f + 1e-6;  % Slightly above specificity to avoid division by zero
    s_upper = 1;

    % Check if a valid solution exists in the interval
    if std_dev_function(s_lower) * std_dev_function(s_upper) > 0
        error('No solution found for the given parameters. Try adjusting inputs.');
    end

    % Numerically solve for sensitivity using fzero
    s = fzero(std_dev_function, [s_lower, s_upper]);
end
