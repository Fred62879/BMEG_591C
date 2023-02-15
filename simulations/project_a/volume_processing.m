function volume_processing()
    fn = '/media/fred/working_drive/2022W2/bmeg591c/simulations/data/RawOCT.mat';
    raw_data = load(fn);
    rawOCT = getfield(raw_data,'rawOCT');

    dispCoeffs = -1;

    for frame_num = 1:2 %size(rawOCT, 3)
        raw_data = rawOCT(:,:,frame_num);

        fft_Data = fft(hilbert(raw_data));

        % binary filter
        calSigMin = [29,40];
        winFunc = zeros(size(fft_Data));
        winFunc(calSigMin(1):calSigMin(2),:) = 1;
        cal_fftData = fft_Data .* winFunc;

        % get calibrated raw data
        cal_rawData = ifft(cal_fftData);

        % resampling
        rawData_rescaled = reSampling_CalSig(fft_Data, cal_rawData);
        fftData_rescaled = fft(rawData_rescaled);

        % phase compensation
        cplxConjX = fftData_rescaled .* ...
            repmat( conj(fftData_rescaled(:,1)), [1 size(fftData_rescaled,2)] );

        depthIdx = 34;

        fftData_1D = fftData_rescaled(:,1);
        rawData_PhaseComp = compPhaseShift(fftData_1D, fftData_rescaled, depthIdx);
        fftData_PhaseComp = fft(rawData_PhaseComp);

        % Fixed pattern noise (FPN) removal
        rawData_FPNSub = rawData_PhaseComp...
            - (repmat(median(real(rawData_PhaseComp), 2), [1,size(rawData_PhaseComp,2)]) ...
               +1j .* repmat(median(imag(rawData_PhaseComp), 2), [1, size(rawData_PhaseComp, 2)]));
        fftData_FPNSub = fft(rawData_FPNSub);

        % hanning windowing
        rawData_hanWin = rawData_FPNSub ...
            .* repmat(hann(size(rawData_FPNSub,1)), [1 size(rawData_FPNSub,2)]);
        fftData_hanWin = fft(rawData_hanWin);

        % dispersion compensation
        if dispCoeffs == -1
            disp("calculating dispersion coefficient")
            maxDispOrders = 3;
            coeffRange = 10;
            depthROI = [45,275]; % guessed, removed dark region (use any image before)
            dispCoeffs = setDispCoeffs(rawData_FPNSub, depthROI, maxDispOrders, coeffRange);
        end

        rawData_dispComp = compDisPhase(rawData_FPNSub, maxDispOrders, dispCoeffs);
        fftData_dispComp = fft(rawData_dispComp);

        % collect result
        procd_data(:,:,frame_num) = fftData_dispComp(1:size(raw_data, 1)/2,:);

    end

    subplot(1,2,1); imagesc( imadjust(mat2gray(20 .* log10(...
         abs(procd_data(:,:,1)))))); colormap(gray);
    subplot(1,2,2); imagesc( imadjust(mat2gray(20 .* log10(...
         abs(procd_data(:,:,2)))))); colormap(gray);
    imwrite(procd_data(:,:,1),"../res.tiff")
    imwrite(procd_data(:,:,2),"../res.tiff","WriteMode","append")

    a=imread('../res.tiff');
    subplot(1,2,1); imagesc( imadjust(mat2gray(20 .* log10(...
         abs(a(:,:,1)))))); colormap(gray);
    subplot(1,2,2); imagesc( imadjust(mat2gray(20 .* log10(...
         abs(a(:,:,2)))))); colormap(gray);

end