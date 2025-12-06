function y = smooth_custom(x, span)
    % smooth_custom - Simple moving average smoother
    % x    : Input data vector
    % span : Number of points to average over (must be odd)
    
    if nargin < 2
        span = 5; % default window size
    end

    if mod(span, 2) == 0
        error('Span must be an odd number.');
    end

    n = length(x);
    y = zeros(size(x));

    halfSpan = floor(span / 2);

    for i = 1:n
        % Determine window range
        startIdx = max(1, i - halfSpan);
        endIdx   = min(n, i + halfSpan);

        % Compute average
        y(i) = mean(x(startIdx:endIdx));
    end
end

