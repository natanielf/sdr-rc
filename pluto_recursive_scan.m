function detected = pluto_recursive_scan(device, fStart, fStop, maxSpanHz, sampleRate, samplesPerCapture, threshold_dB, minSpanHz)
% Returns detected ranges (Nx2 array) or single frequencies when minSpanHz reached.
% Example:
% detected = pluto_recursive_scan('Pluto',2.4e9,2.4835e9,20e6,2e6,2^14,-40,100e3)

if nargin<9, minSpanHz = 100e3; end
if nargin<8, threshold_dB = -40; end

% Create receiver template
rx = sdrrx(device, ...
    'BasebandSampleRate', sampleRate, ...
    'SamplesPerFrame', samplesPerCapture, ...
    'OutputDataType','double');

% Recursive worker
detected = scan_range(fStart, fStop);

release(rx);

    function out = scan_range(flo, fhi)
        span = fhi - flo;
        if span <= maxSpanHz
            % Single-center capture covering this span (use sampleRate <= span)
            center = (flo+fhi)/2;
            rx.CenterFrequency = center;
            pause(0.08); % settle after retune
            s = rx();
            if isempty(s), out = []; return; end
            % PSD
            nfft = max(2^nextpow2(4096), 2^nextpow2(round(span/sampleRate*4096)));
            [pxx,f] = pwelch(s, hann(1024), 512, nfft, sampleRate, 'centered');
            fAbs = f + center;
            % Limit PSD to the flo..fhi window, compute max level
            idx = fAbs>=flo & fAbs<=fhi;
            if any(idx)
                maxdB = max(10*log10(pxx(idx)));
            else
                maxdB = -Inf;
            end
            if maxdB >= threshold_dB
                if span <= minSpanHz
                    out = [flo, fhi]; % resolved segment (or treat as frequency)
                else
                    % subdivide
                    mid = (flo+fhi)/2;
                    left = scan_range(flo, mid);
                    right = scan_range(mid, fhi);
                    out = vertcat(left, right);
                end
            else
                out = [];
            end
        else
            % span larger than maxSpanHz: split into chunks and scan each
            n = ceil(span / maxSpanHz);
            edges = linspace(flo, fhi, n+1);
            out = [];
            for k=1:n
                outk = scan_range(edges(k), edges(k+1));
                if ~isempty(outk), out = vertcat(out, outk); end
            end
        end
    end
end

