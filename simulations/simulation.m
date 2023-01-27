
function simulation()

    function dc_fft()
    % interference signal with DC-FFT s7
        for i = 1:90
            Fs = 1024;                % sampling freq
            T = 1/Fs;                 % sampling period
            L = 1024;                 % length of signal
            k = ((0:L-1)*T)';         % wavenumber sampling vector
            A = 10;                   % amplitude
            dz = 50;                  % pathlength difference (um)
            DC = A + i;               % DC

            I_k = A*cos(2*pi*k*dz);        % interference signal
            I_k = I_k + ones([L 1]) .* DC; % interference signal with DC
            I_k = I_k .* gausswin(L);      % interference signal with Gauss
            I_k = I_k ./ max(I_k(:));      % normed interference signal

            I_z = abs(fft(I_k));      % fft signal
            I_z = I_z ./ max(I_z(:)); % normed fft signal

            plot(I_k);
            xlim([1 L]);

            plot(I_z);
            xlim([1 L/2]);
            end
    end

    function gaus_window_fft()
    % Gaussian windowing s11
        for i =1:100
            Fs = 1024;                % sampling freq
            T = 1/Fs;                 % sampling period
            L = 1024;                 % length of signal
            k = ((0:L-1)*T)';         % wavenumber sampling vector
            A = 10;                   % amplitude
            dz = 50;                  % pathlength difference (um)
            I_k = A*cos(2*pi*k*dz);   % interference signal
            DC = A;                   % DC
            alpha = 2.5 + (0.1*i);    % width of gauss window

            I_k = I_k + ones([L 1]) .* DC;   % interference signal with DC
            I_k = I_k .* gausswin(L, alpha); % interference signal with Gauss
            I_k = I_k ./ max(I_k(:));        % normed interference signal

            I_z = abs(fft(I_k));      % fft signal
            I_z = I_z ./ max(I_z(:)); % normed fft signal

            plot(I_k);
            xlim([1 L]);

            plot(I_z);
            xlim([1 L/2]);
        end
    end

    function extractor()
    % s13
        Fs = 1024;                % sampling freq
        T = 1/Fs;                 % sampling period
        L = 1024;                 % length of signal
        k = ((0:L-1)*T)';         % wavenumber sampling vector
        A = 10;                   % amplitude
        dz = 50;                  % pathlength difference (um)
        DC = ones([L 1])*10;      % DC
        I_k = A*cos(2*pi*k*dz);   % inferference signal
        I_k = I_k + DC;           % inferference signal with DC
        I_k = I_k .* gausswin(L); % interference signal with Gaussian
        I_z = fft(I_k);           % fft signal

        % binary window
        BinWin = zeros([L,1]);
        BinWin(44:57) = 1;

        % I_z x binary window
        I_z_BinWin = I_z .* BinWin;

        % spectrum amplitude & phase
        I_k_BinWin = ifft(I_z_BinWin);        % extracted inferference fringe
        I_k_Spectrum = abs(I_k_BinWin);       % spectrum amplitude
        I_k_Phase = angle(I_k_BinWin);        % spectrum phase
        I_k_Phase_unwrap = unwrap(I_k_Phase); % spectrum phase - unwrap
    end

    function opd()
    % optical path-length difference s15
        for i = 1:50
            Fs = 1024;                % sampling freq
            T = 1/Fs;                 % sampling period
            L = 1024;                 % length of signal
            k = ((0:L-1)*T)';         % wavenumber sampling vector
            A = 10;                   % amplitude
            DC = ones([L 1])*10;      % DC

            dz = 10 + (i-1)*10;       % pathlength difference (um)

            I_k = A*cos(2*pi*k*dz);   % inferference signal
            I_k = I_k + DC;           % inferference signal with DC
            I_k = I_k .* gausswin(L); % interference signal with Gaussian

            I_z = abs(fft(I_k));

            plot(I_k, 'color', 'b');
            xlim([1 L]);

            plot(I_z, 'color', 'r');
            xlim([1 L/2]);
        end
    end

    function dispersion_2nd()
    % s19
        for i = 1:50
            Fs = 1024;                % sampling freq
            T = 1/Fs;                 % sampling period
            L = 1024;                 % length of signal
            k = ((0:L-1)*T)';         % wavenumber sampling vector
            A = 10;                   % amplitude
            DC = ones([L 1])*10;      % DC

            disp2nd = i .* (linspace(-1,1,L).^2); % 2nd order dispersion

            I_k = A*cos(2*pi*k*dz);   % inferference signal
            I_k = I_k + DC;           % inferference signal with DC
            I_k = I_k .* gausswin(L); % interference signal with Gaussian

            I_z = abs(fft(I_k));

            plot(disp2nd);
            xlim([1 L]);
            ylim([0 50]);

            plot(abs(I_z));
            xlim([1 L/2]);
        end
    end

    function dispersion_3rd()
    % s23
        for i = 1:50
            Fs = 1024;                % sampling freq
            T = 1/Fs;                 % sampling period
            L = 1024;                 % length of signal
            k = ((0:L-1)*T)';         % wavenumber sampling vector
            A = 10;                   % amplitude
            DC = ones([L 1])*10;      % DC

            disp3rd = i .* (linspace(-1,1,L).^3); % 2nd order dispersion

            I_k = A*cos(2*pi*k*dz);   % inferference signal
            I_k = I_k + DC;           % inferference signal with DC
            I_k = I_k .* gausswin(L); % interference signal with Gaussian

            I_z = fft(I_k);

            plot(disp3rd);
            xlim([1 L]);
            ylim([-50 50]);

            plot(abs(I_z));
            xlim([1 L/2]);
        end
    end

    function non_linr_sample()
    % 29
        Fs = 1024;                % sampling freq
        T = 1/Fs;                 % sampling period
        L = 1024;                 % length of signal
        A = 10;                   % amplitude
        dz = 50;                  % pathlength difference (um)
        DC = ones([L 1])*10;      % DC
        k_linear = ((0:L-1)*T)';  % linear k sampling vector
        k_Nlinear = flip((1 - (k_linear' .^1.2)))'; % non-linear k sampling vec

        I_k_lin = A*cos(2*pi*k_linear*dz);
        I_k_lin = I_k_lin + DC;
        I_k_lin = I_k_lin .* gausswin(L);
        I_k_Nlin = interp1(k_linear, I_k_lin, k_Nlinear, 'spline');

        I_z_lin = fft(I_k_lin);
        I_z_Nlin = fft(I_k_Nlin);
    end

end