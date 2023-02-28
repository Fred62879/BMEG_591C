function reference_processing()
    
    fn = '/media/fred/working_drive/2022W2/bmeg591c/simulations/data/project_a/RawOCT.mat';
    %fn = 'RawOCT.mat';
    raw_data = load(fn);
    rawOCT = getfield(raw_data,'rawOCT');

    ref_rawData = rawOCT(:,:,250);
    ref_fftData = fft(ref_rawData);
    ref_hilData = hilbert(ref_rawData);
    ref_fhilData = fft(ref_hilData);

%     subplot(1,6,1); plot(ref_rawData(:,:));
%     subplot(1,6,2); plot(ref_rawData(:,1));
%     subplot(1,6,3); plot(abs(ref_fftData(:,:))); ylim([0 4e6]);
%     subplot(1,6,4); plot(abs(ref_hilData(:,1)));
%     subplot(1,6,5); plot(imag(ref_hilData(:,1)));
%     subplot(1,6,6); plot(abs(ref_fhilData(:,1))); ylim([0 4e6]);

    ref_fftData = fft(hilbert(ref_rawData));

%     subplot(1,2,1); plot(abs(ref_fftData(1:end/2,:))); ylim([0 4e6]);
%     subplot(1,2,2); plot(abs(ref_fftData(1:end/2,1))); ylim([0 4e6]);
%     subplot(1,2,1); imagesc( imadjust(mat2gray(abs(ref_fftData(1:end/2,:))))); colormap(gray);
%     subplot(1,2,2); imagesc( imadjust(mat2gray(20 .* log10(...
%         abs(ref_fftData(1:end/2,:)))))); colormap(gray);

    % binary filter
    calSigMin = [29,40];
    winFunc = zeros(size(ref_fftData));
    winFunc(calSigMin(1):calSigMin(2),:) = 1;
    ref_cal_fftData = ref_fftData .* winFunc;

%     figure(); imagesc( imadjust(mat2gray(20 .* log10(...
%         abs(ref_cal_fftData(1:end/2,:)))))); colormap(gray);hold on;

    % get calibrated raw data
    ref_cal_rawData = ifft(ref_cal_fftData);
%     figure(); plot(real(ref_cal_rawData));
%     figure(); plot(real(ref_cal_rawData(1:end/2,1)));

    % resampling
    ref_rawData_rescaled = reSampling_CalSig(ref_fftData, ref_cal_rawData);
    ref_fftData_rescaled = fft(ref_rawData_rescaled);

%     plot(real(ref_rawData_rescaled));
%     subplot(1,2,1); imagesc( imadjust(mat2gray(20 .* log10(...
%         abs(ref_fftData(1:end/2,:)))))); colormap(gray);
%     subplot(1,2,2); imagesc( imadjust(mat2gray(20 .* log10(...
%         abs(ref_fftData_rescaled(1:end/2,:)))))); colormap(gray);

    % spectra shift / phase compensation
    cplxConjX = ref_fftData_rescaled .* ...
        repmat( conj(ref_fftData_rescaled(:,1)), [1 size(ref_fftData_rescaled,2)] );

    % plot(abs(cplxConjX(:,1:10))); xlim([ 10 100]) % find optimal depthIdx
    depthIdx = 34;

    ref_fftData_1D = ref_fftData_rescaled(:,1);
    ref_rawData_PhaseComp = compPhaseShift(ref_fftData_1D, ref_fftData_rescaled, depthIdx);
    ref_fftData_PhaseComp = fft(ref_rawData_PhaseComp);

%     subplot(1,2,1); imagesc( imadjust(mat2gray(20 .* log10(...
%          abs(ref_fftData(1:end/2,:)))))); colormap(gray);hold on;
%     subplot(1,2,2); imagesc( imadjust(mat2gray(20 .* log10(...
%          abs(ref_fftData_PhaseComp(1:end/2,:)))))); colormap(gray);

    % Fixed pattern noise (FPN) removal
    ref_rawData_FPNSub = ref_rawData_PhaseComp...
        - (repmat(median(real(ref_rawData_PhaseComp), 2), [1,size(ref_rawData_PhaseComp,2)]) ...
           +1j .* repmat(median(imag(ref_rawData_PhaseComp), 2), [1, size(ref_rawData_PhaseComp, 2)]));
    ref_fftData_FPNSub = fft(ref_rawData_FPNSub);

%     subplot(1,2,1); imagesc( imadjust(mat2gray(20 .* log10(...
%          abs(ref_fftData(1:end/2,:)))))); colormap(gray);hold on;
%     subplot(1,2,2); imagesc( imadjust(mat2gray(20 .* log10(...
%          abs(ref_fftData_FPNSub(1:end/2,:)))))); colormap(gray);

    % hanning windowing
    ref_rawData_hanWin = ref_rawData_FPNSub ...
        .* repmat(hann(size(ref_rawData_FPNSub,1)), [1 size(ref_rawData_FPNSub,2)]);
    ref_fftData_hanWin = fft(ref_rawData_hanWin);

%     subplot(1,2,1); plot(real(ref_rawData_FPNSub));
%     subplot(1,2,2); plot(real(ref_rawData_hanWin));

%     subplot(1,2,1); imagesc( imadjust(mat2gray(20 .* log10(...
%          abs(ref_fftData_FPNSub(1:end/2,:)))))); colormap(gray);
%     subplot(1,2,2); imagesc( imadjust(mat2gray(20 .* log10(...
%          abs(ref_fftData_hanWin(1:end/2,:)))))); colormap(gray);

    % dispersion estimation
    maxDispOrders = 5;
    coeffRange = 10;
    depthROI = [40,300]; % guessed, removed dark region (use any image before)

    dispCoeffs = setDispCoeffs(ref_rawData_hanWin, depthROI, maxDispOrders, coeffRange);
    disp("Dispersion compensation coefficients are:");
    disp(dispCoeffs);
    ref_rawData_dispComp = compDisPhase(ref_rawData_hanWin, maxDispOrders, dispCoeffs);
    ref_fftData_dispComp = fft(ref_rawData_dispComp);

%     subplot(1,2,1); imagesc( imadjust(mat2gray(20 .* log10(...
%          abs(ref_fftData(1:end/2,:)))))); colormap(gray);hold on;
%     subplot(1,2,2); imagesc( imadjust(mat2gray(20 .* log10(...
%          abs(ref_fftData_dispComp(1:end/2,:)))))); colormap(gray);
end
