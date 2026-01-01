function wwost_global_demo()
%WWOST_GLOBAL_DEMO Global window-width optimization for the S-transform.
%   Demonstrates selection of a single optimal window-width exponent p that
%   maximizes a global compactness measure CM(p) = 1 / L1(p) where L1 is
%   the L1 norm of a unit-energy S-transform magnitude.
%
%   The routine:
%     1) Generates a test signal.
%     2) Defines a candidate list of p values.
%     3) Computes S_p for each p, normalizes it to unit L2 energy
%        (integrated over time and frequency), and evaluates CM(p) using
%        the positive-frequency (half-spectrum) with non-DC bins doubled to
%        approximate full-spectrum energy.
%     4) Picks p_opt = argmax_p CM(p).
%     5) Recomputes and visualizes the p = 1 and p = p_opt transforms with
%        shared color limits, and annotates p_opt.

    % Test signal and sampling.
    [t, s] = demo_signal();
    dt = t(2) - t(1);

    % Candidate window-width exponents (limit list length for memory safety).
    p_list = [0.5, 0.75, 1, 1.25, 1.5, 2];

    num_p = numel(p_list);
    S_norm = cell(1, num_p);
    CM = zeros(1, num_p);

    f_ref = [];
    df = [];

    % Compute transforms, normalize to unit L2 energy, and evaluate CM.
    for idx = 1:num_p
        p = p_list(idx);
        [S_p_full, f_full] = ST(s, t, p);
        [f, S_p, weights] = positive_spectrum(f_full, S_p_full);

        if isempty(f_ref)
            f_ref = f;
            df = f(2) - f(1);
        else
            if numel(f) ~= numel(f_ref) || max(abs(f - f_ref)) > 1e-9
                error('Frequency grids do not match across p values.');
            end
        end

        % Normalize using the positive-frequency (half-spectrum) energy,
        % doubling non-DC bins to approximate full-spectrum energy.
        energy = sqrt(sum(weights(:) .* sum(abs(S_p).^2, 2)) * dt * df);
        if energy == 0
            error('Zero transform energy encountered.');
        end
        S_norm{idx} = S_p / energy;

        L1 = sum(weights(:) .* sum(abs(S_norm{idx}), 2)) * dt * df;
        CM(idx) = 1 ./ max(L1, eps);
    end

    % Select optimal p.
    [~, opt_idx] = max(CM);
    p_opt = p_list(opt_idx);

    % Recompute canonical transforms for plotting.
    [S_p1_full, f_full] = ST(s, t, 1);
    [f, S_p1, weights] = positive_spectrum(f_full, S_p1_full);

    [S_popt_full, ~] = ST(s, t, p_opt);
    [~, S_popt] = positive_spectrum(f_full, S_popt_full);

    energy_p1 = sqrt(sum(weights(:) .* sum(abs(S_p1).^2, 2)) * dt * df);
    energy_popt = sqrt(sum(weights(:) .* sum(abs(S_popt).^2, 2)) * dt * df);
    S_p1 = S_p1 / energy_p1;
    S_popt = S_popt / energy_popt;

    clim_max = max([max(abs(S_p1), [], 'all'), max(abs(S_popt), [], 'all')]);

    % Visualization.
    figure('Name', 'Global WWOST Demo');
    tiledlayout(2, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

    % (a) CM vs p.
    nexttile([1 2]);
    plot(p_list, CM, 'LineWidth', 1.5);
    hold on;
    plot(p_opt, CM(opt_idx), 'ro', 'MarkerFaceColor', 'r');
    grid on;
    xlabel('p (window-width exponent)');
    ylabel('CM(p) = 1 / L1(p)');
    title(sprintf('Compactness Measure, p_{opt} = %.3g', p_opt));

    % (b) TF magnitude maps for p = 1 and p = p_opt.
    nexttile;
    imagesc(t, f, abs(S_p1));
    axis xy;
    xlabel('Time (s)');
    ylabel('Frequency (Hz)');
    title('p = 1 (unit L2 energy)');
    colorbar;
    caxis([0, clim_max]);

    nexttile;
    imagesc(t, f, abs(S_popt));
    axis xy;
    xlabel('Time (s)');
    ylabel('Frequency (Hz)');
    title(sprintf('p = %.3g (unit L2 energy)', p_opt));
    colorbar;
    caxis([0, clim_max]);

    fprintf('Optimal global p: %.3g\n', p_opt);
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

function [f_pos, S_pos, weights] = positive_spectrum(f, S)
%POSITIVE_SPECTRUM Extract non-negative frequencies with dimension safety.
    pos_idx = f >= 0;
    f_pos = f(pos_idx);

    subs = repmat({':'}, 1, ndims(S));
    subs{1} = pos_idx;
    S_pos = S(subs{:});

    % Double non-DC bins to approximate full-spectrum energy from half-spectrum.
    weights = ones(size(f_pos));
    weights(f_pos > 0) = 2;
end
