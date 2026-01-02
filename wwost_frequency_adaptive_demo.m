function wwost_frequency_adaptive_demo()
%WWOST_FREQUENCY_ADAPTIVE_DEMO Frequency-wise window-width optimization.
%   Implements the per-frequency WWOST compactness search described in the
%   accompanying paper. For each candidate exponent p in p_list, the script
%   computes the S-transform S_p(t, f) (using sigma(f) = 1 / |f|^p, matching
%   ST.m), evaluates a compactness measure CM(f, p) = 1 / L_q(f, p), and
%   selects p_opt(f) = argmax_p CM(f, p). The optimized map S^{p_opt(f)} is
%   then assembled and plotted alongside the time-domain waveform and the
%   p_opt(f) profile.
%
%   Key steps:
%     1) Define a three-component signal x1(t) sampled at Fs = 256 Hz.
%     2) Optionally zero-pad (and taper edges) to compare padding effects.
%     3) Sweep p_list, accumulate S_p, and compute CM(f, p) via an L_q norm.
%     4) Build S_popt using the best p at each frequency.
%     5) Visualize x1(t), |S^{p_opt(f)}(t, f)| with axis xy, and p_opt(f).

    % Sampling and base signal definition.
    Fs = 256;
    dt = 1 / Fs;
    t = 0:dt:(1 - dt); % 1 second, uniformly sampled

    % Three-component test signal x1(t).
    tone = 0.8 * sin(2 * pi * 18 * t);
    chirp_up = 0.6 * chirp(t, 12, t(end), 96, 'quadratic');
    burst = 1.1 * exp(-((t - 0.6) / 0.07).^2) .* cos(2 * pi * 72 * t);
    x1 = tone + chirp_up + burst;

    % Optional padding/tapering to mirror paper-style preprocessing.
    zero_pad_factor = 1.0; % set > 1 (e.g., 2) to append zeros
    taper_edges = true;    % apply gentle cosine ramps at boundaries
    taper_fraction = 0.04; % fraction of samples tapered on each edge

    [t_tf, x_tf] = pad_and_taper(t, x1, zero_pad_factor, taper_edges, taper_fraction);

    % Compactness search parameters.
    p_list = 0.01:0.01:1.0;
    q = 0.2; % L_q exponent for compactness
    num_p = numel(p_list);

    % Compute S_p for all p values and evaluate CM(f, p).
    [S_first, f] = ST(x_tf, t_tf, p_list(1));
    num_f = numel(f);
    num_t = numel(t_tf);
    dt_tf = t_tf(2) - t_tf(1);

    S_stack = complex(zeros(num_f, num_t, num_p));
    CM = zeros(num_f, num_p);

    S_stack(:, :, 1) = S_first;
    CM(:, 1) = 1 ./ max(sum(abs(S_first) .^ q, 2) * dt_tf, eps);

    for idx = 2:num_p
        p = p_list(idx);
        [S_p, f_cur] = ST(x_tf, t_tf, p);

        if numel(f_cur) ~= num_f || max(abs(f_cur - f)) > 1e-12
            error('Frequency grids do not match across p values.');
        end

        S_stack(:, :, idx) = S_p;
        CM(:, idx) = 1 ./ max(sum(abs(S_p) .^ q, 2) * dt_tf, eps);
    end

    % Frequency-wise optimal p and assembled map.
    [~, opt_idx] = max(CM, [], 2);
    p_opt = p_list(opt_idx);

    S_popt = zeros(num_f, num_t);
    for ff = 1:num_f
        S_popt(ff, :) = S_stack(ff, :, opt_idx(ff));
    end

    % Visualization.
    figure('Name', 'Frequency-Adaptive WWOST Demo');
    tiledlayout(3, 1, 'Padding', 'compact', 'TileSpacing', 'compact');

    % Time waveform (non-padded signal).
    nexttile;
    plot(t, x1, 'LineWidth', 1.25);
    grid on;
    xlabel('Time (s)');
    ylabel('x_1(t)');
    title('Three-Component Signal x_1(t)');

    % Optimized time-frequency magnitude.
    nexttile;
    imagesc(t_tf, f, abs(S_popt));
    axis xy;
    xlabel('Time (s)');
    ylabel('Frequency (Hz)');
    title('|S^{p_{opt}(f)}(t, f)|');
    colorbar;

    % Frequency-dependent optimal p.
    nexttile;
    plot(f, p_opt, 'LineWidth', 1.25);
    grid on;
    xlabel('Frequency (Hz)');
    ylabel('p_{opt}(f)');
    title(sprintf('p_{opt}(f) using CM(f, p) with q = %.2f', q));

    fprintf('Computed frequency-wise p_{opt} over %d frequency bins. Example range: [%.3g, %.3g]\n', ...
        num_f, min(p_opt), max(p_opt));
end

function [t_out, x_out] = pad_and_taper(t, x, pad_factor, taper_edges, taper_fraction)
%PAD_AND_TAPER Optional zero padding and edge tapering.
    if pad_factor < 1
        error('pad_factor must be >= 1.');
    end

    dt = t(2) - t(1);
    N = numel(t);
    N_pad = ceil(N * pad_factor);

    t_out = (0:N_pad-1) * dt;
    x_out = zeros(1, N_pad);
    x_out(1:N) = x;

    if taper_edges && taper_fraction > 0
        taper_len = max(1, round(taper_fraction * N_pad));
        ramp = 0.5 * (1 - cos(pi * (0:taper_len-1) / taper_len));

        x_out(1:taper_len) = x_out(1:taper_len) .* ramp;
        x_out(end-taper_len+1:end) = x_out(end-taper_len+1:end) .* fliplr(ramp);
    end
end
