function reference_processing()
    
    fn = '/media/fred/working_drive/2022W2/bmeg591c/simulations/data/project_a/RawOCT.mat';
    %fn = 'RawOCT.mat';
    raw_data = load(fn);
    rawOCT = getfield(raw_data,'rawOCT');

    ref_rawData = rawOCT(:,:,250);
    ref_fftData = fft(hilbert(ref_rawData));

    % binary filter
    calSigMin = [29,40];
    winFunc = zeros(size(ref_fftData));
    winFunc(calSigMin(1):calSigMin(2),:) = 1;
    ref_cal_fftData = ref_fftData .* winFunc;

    % get calibrated raw data
    ref_cal_rawData = ifft(ref_cal_fftData);

    % resampling
    ref_rawData_rescaled = reSampling_CalSig(ref_fftData, ref_cal_rawData);
    ref_fftData_rescaled = fft(ref_rawData_rescaled);

    % spectra shift / phase compensation
    cplxConjX = ref_fftData_rescaled .* ...
        repmat( conj(ref_fftData_rescaled(:,1)), [1 size(ref_fftData_rescaled,2)] );

    depthIdx = 34;

    ref_fftData_1D = ref_fftData_rescaled(:,1);
    ref_rawData_PhaseComp = compPhaseShift(ref_fftData_1D, ref_fftData_rescaled, depthIdx);
    ref_fftData_PhaseComp = fft(ref_rawData_PhaseComp);

    % Fixed pattern noise (FPN) removal
    ref_rawData_FPNSub = ref_rawData_PhaseComp...
        - (repmat(median(real(ref_rawData_PhaseComp), 2), [1,size(ref_rawData_PhaseComp,2)]) ...
           +1j .* repmat(median(imag(ref_rawData_PhaseComp), 2), [1, size(ref_rawData_PhaseComp, 2)]));
    ref_fftData_FPNSub = fft(ref_rawData_FPNSub);

    % hanning windowing
    ref_rawData_hanWin = ref_rawData_FPNSub ...
        .* repmat(hann(size(ref_rawData_FPNSub,1)), [1 size(ref_rawData_FPNSub,2)]);
    ref_fftData_hanWin = fft(ref_rawData_hanWin);

    % dispersion estimation
    maxDispOrders = 5;
    coeffRange = 10;
    depthROI = [40,300]; % guessed, removed dark region (use any image before)

    dispCoeffs = setDispCoeffs(ref_rawData_hanWin, depthROI, maxDispOrders, coeffRange);
    disp("Dispersion compensation coefficients are:");
    disp(dispCoeffs);
    ref_rawData_dispComp = compDisPhase(ref_rawData_hanWin, maxDispOrders, dispCoeffs);
end
