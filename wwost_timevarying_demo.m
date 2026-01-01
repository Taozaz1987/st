function wwost_timevarying_demo()
%WWOST_TIMEVARYING_DEMO Time-varying window-width optimization.
%   Demonstrates selecting the optimal exponent p per time index using
%   compactness measure CM_t(p) = 1 / L1_t(p), where L1_t integrates the
%   S-transform magnitude over frequency at each time sample.
%
%   Steps:
%     1) Compute the baseline S-transform at p = 1 and its energy.
%     2) For each p in p_list, compute S_p, normalize using the baseline
%        energy, and evaluate CM_t(p).
%     3) Derive p_opt(t) = argmax_p CM_t(p).
%     4) Build the WWOST time-frequency map by selecting columns from the
%        best-performing S_p at each time index.
%     5) Visualize p_opt(t) and the TF magnitude maps for p = 1 and WWOST.

    [t, s] = demo_signal();
    dt = t(2) - t(1);

    % Candidate exponents (keep modest length to limit memory footprint).
    p_list = [0.75, 1, 1.25, 1.5, 1.75];
    num_p = numel(p_list);

    % Baseline transform at p = 1.
    [S_baseline, f] = ST(s, t, 1);
    df = f(2) - f(1);
    baseline_energy = sqrt(sum(abs(S_baseline).^2, 'all') * dt * df);
    S_baseline = S_baseline / baseline_energy;

    num_f = numel(f);
    num_t = numel(t);

    % Storage for normalized transforms (kept small for demo).
    S_stack = complex(zeros(num_f, num_t, num_p));
    CM_t = zeros(num_p, num_t);

    for idx = 1:num_p
        p = p_list(idx);
        [S_p, f_cur] = ST(s, t, p);
        if numel(f_cur) ~= num_f || max(abs(f_cur - f)) > 1e-9
            error('Frequency grids do not match across p values.');
        end

        % Normalize using baseline energy to enable fair comparison.
        S_p = S_p / baseline_energy;

        S_stack(:, :, idx) = S_p;
        L1_t = df * sum(abs(S_p), 1);
        CM_t(idx, :) = 1 ./ max(L1_t, eps);
    end

    % Time-varying optimal p and assembled WWOST map.
    [~, opt_idx_t] = max(CM_t, [], 1);
    p_opt_t = p_list(opt_idx_t);

    S_wwost = zeros(num_f, num_t);
    for tt = 1:num_t
        S_wwost(:, tt) = S_stack(:, tt, opt_idx_t(tt));
    end

    clim_max = max([max(abs(S_baseline), [], 'all'), max(abs(S_wwost), [], 'all')]);

    figure('Name', 'Time-Varying WWOST Demo');
    tiledlayout(3, 1, 'Padding', 'compact', 'TileSpacing', 'compact');

    % p_opt(t) trajectory.
    nexttile;
    plot(t, p_opt_t, 'LineWidth', 1.5);
    grid on;
    xlabel('Time (s)');
    ylabel('p_{opt}(t)');
    title('Time-Varying Optimal Exponent');

    % Baseline TF map.
    nexttile;
    imagesc(t, f, abs(S_baseline));
    axis xy;
    xlabel('Time (s)');
    ylabel('Frequency (Hz)');
    title('Baseline p = 1 (normalized by baseline energy)');
    colorbar;
    caxis([0, clim_max]);

    % WWOST TF map.
    nexttile;
    imagesc(t, f, abs(S_wwost));
    axis xy;
    xlabel('Time (s)');
    ylabel('Frequency (Hz)');
    title('WWOST with p_{opt}(t) (baseline-normalized)');
    colorbar;
    caxis([0, clim_max]);

    fprintf('Computed p_opt(t) over %d samples. Example values: min=%.3g, max=%.3g\n', ...
        num_t, min(p_opt_t), max(p_opt_t));
end

function [t, s] = demo_signal()
%DEMO_SIGNAL Create a broadband test signal for WWOST examples.
    T = 2; % seconds
    fs = 500; % Hz
    t = (0:1/fs:T-1/fs);

    tone1 = 1.2 * sin(2 * pi * 30 * t);
    chirp1 = 0.8 * chirp(t, 10, T, 120, 'quadratic');
    am = (1 + 0.2 * sin(2 * pi * 0.5 * t));
    fm = cos(2 * pi * 5 * t + 0.4 * sin(2 * pi * 1.5 * t));
    burst = am .* cos(2 * pi * (60 + 5 * fm) .* t) .* exp(-((t - 1).^2) / 0.08);

    s = tone1 + chirp1 + burst;
end

