function xq = quantize_signal(x, levels)
    % Quantize each sample in x using the provided levels
    xq = zeros(size(x));
    for i = 1:numel(x)
        [~, idx] = min(abs(x(i) - levels)); % Find the closest quantization level
        xq(i) = levels(idx); % Assign the quantized value
    end
end
