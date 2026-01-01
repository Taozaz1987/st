function [stTFR, f] = ST(s, t, p)
%ST Compute the forward S-transform of a 1-D real signal.
%   [stTFR, f] = ST(s, t) returns the Stockwell transform (S-transform)
%   time-frequency representation stTFR with size
%   [num_frequencies x num_time_samples] and the corresponding frequency
%   vector f. The frequency resolution is consistent with the sampling rate
%   implied by the time vector t.
%
%   Inputs:
%       s - Real-valued input signal vector.
%       t - Time vector (must match the length of s and have uniform
%           spacing).
%       p - Window-width exponent (scalar, positive). The Gaussian standard
%           deviation follows sigma(f; p) = 1 / |f|^p, producing the window
%           factor |f|.^p .* exp(-0.5 * (|f|.^(2*p)) .* t^2). Default: 1.
%
%   Outputs:
%       stTFR - Complex S-transform matrix (frequency x time).
%       f     - Frequency vector corresponding to rows of stTFR.

% Ensure inputs are row vectors for broadcasting convenience.
s = s(:).';
t = t(:).';
N = numel(s);

if numel(t) ~= N
    error('Signal and time vectors must have the same length.');
end

if nargin < 3 || isempty(p)
    p = 1;
end
validateattributes(p, {'numeric'}, {'scalar', 'real', 'positive'}, mfilename, 'p');

% Sampling characteristics.
dt = mean(diff(t));
if any(abs(diff(t) - dt) > eps(max(abs(t)))*10)
    warning('Time vector is not uniformly sampled; results may be inaccurate.');
end
fs = 1 / dt;

% Two-sided frequency vector centered at zero; df matches sampling rate.
f = (-floor(N/2):ceil(N/2)-1).' * (fs / N);

% Centered time vector used for zero-phase Gaussian construction.
t_centered = ((0:N-1) - floor(N/2)) * dt;

% Modulate the signal for each frequency.
modulation = exp(-1j * 2 * pi * (f * t));
s_mod = modulation .* s;

% Frequency-dependent Gaussian windows (zero-centered) with exponent p.
abs_f = abs(f);
gaussian = (abs_f .^ p) .* exp(-0.5 * ((abs_f .^ (2 * p)) .* (t_centered .^ 2)));

% Handle the zero-frequency row explicitly (limit of the Gaussian window).
gaussian(abs_f == 0, :) = 1;

% Align zero time to the first element for circular convolution via FFT.
gaussian_shifted = ifftshift(gaussian, 2);

% FFT-based convolution for all frequencies simultaneously.
S_mod_fft = fft(s_mod, [], 2);
gaussian_fft = fft(gaussian_shifted, [], 2);
stTFR = dt * ifft(S_mod_fft .* gaussian_fft, [], 2);

% Zero-frequency row equals the local average of the signal.
stTFR(f == 0, :) = mean(s) * ones(1, N);
