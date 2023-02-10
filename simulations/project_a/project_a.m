function project_a
    fn = '/media/fred/working_drive/2022W2/bmeg591c/simulations/data/RawOCT.mat';
    raw_data = load(fn);
    rawOCT = getfield(raw_data,'rawOCT');

    ref_rawData = rawOCT(:,:,250);
    ref_fftData = fft(hilbert(ref_rawData));

    subplot(1,2,1); plot(abs(ref_fftData(1:end/2,:))); ylim([0 4e6]);
    subplot(1,2,2); plot(abs(ref_fftData(1:end/2,1))); ylim([0 4e6]);
    subplot(1,2,1); imagesc( imadjust(mat2gray(abs(ref_fftData(1:end/2,:))))); colormap(gray);
    subplot(1,2,2); imagesc( imadjust(mat2gray(20 .* log10(...
        abs(ref_fftData(1:end/2,:)))))); colormap(gray);

    % binary filter
    calSigMin = [29,40];
    winFunc = zeros(size(ref_fftData));
    winFunc(calSigMin(1):calSigMin(2),:) = 1;
    ref_cal_fftData = ref_fftData .* winFunc;

    figure(); imagesc( imadjust(mat2gray(20 .* log10(...
        abs(ref_cal_fftData(1:end/2,:)))))); colormap(gray);hold on;

    % get calibrated raw data
    ref_cal_rawData = ifft(ref_cal_fftData);
    figure(); plot(real(ref_cal_rawData));
    figure(); plot(real(ref_cal_rawData(1:end/2,1)));

    % resampling
    ref_rawData_rescaled = reSampling_CalSig(ref_fftData, ref_cal_rawData);
    ref_fftData_rescaled = fft(ref_rawData_rescaled);

    plot(real(ref_rawData_rescaled));
    subplot(1,2,1); imagesc( imadjust(mat2gray(20 .* log10(...
        abs(ref_fftData(1:end/2,:)))))); colormap(gray);
    subplot(1,2,2); imagesc( imadjust(mat2gray(20 .* log10(...
        abs(ref_fftData_rescaled(1:end/2,:)))))); colormap(gray);

    % phase compensation
    % todo: figure out best depthIdx
    cplxConjX = ref_fftData_rescaled .* ...
        repmat( conj(ref_fftData_rescaled(:,1)), [1 size(ref_fftData_rescaled,2)] );
    for i=1:10
        subplot(1,10,i); plot(abs(cplxConjX(i,:)));
    end
    for i = 1:10
        subplot(2,5,i); plot(angle(cplxConjX(i,:)));
    end
    subplot(1,2,2); plot(real(cal_RawData()));

    depthIdx = 47;
    ref_fftData_1D = ref_fftData_rescaled(:,1);
    ref_rawData_PhaseComp = compPhaseShift(ref_fftData_1D, ref_fftData_rescaled, depthIdx);
    ref_fftData_PhaseComp = fft(ref_rawData_PhaseComp);

    subplot(1,2,1); imagesc( imadjust(mat2gray(20 .* log10(...
         abs(ref_fftData(1:end/2,:)))))); colormap(gray);hold on;
    subplot(1,2,2); imagesc( imadjust(mat2gray(20 .* log10(...
         abs(ref_fftData_PhaseComp(1:end/2,:)))))); colormap(gray);

    % Fixed pattern noise (FPN) removal
    ref_rawData_FPNSub = ref_rawData_rescaled...
        - (repmat(median(real(ref_rawData_rescaled), 2), [1,size(ref_rawData_rescaled,2)]) ...
           +1j .* repmat(median(imag(ref_rawData_rescaled), 2), [1, size(ref_rawData_rescaled, 2)]));
    ref_fftData_FPNSub = fft(ref_rawData_FPNSub);

    subplot(1,2,1); imagesc( imadjust(mat2gray(20 .* log10(...
         abs(ref_fftData(1:end/2,:)))))); colormap(gray);hold on;
    subplot(1,2,2); imagesc( imadjust(mat2gray(20 .* log10(...
         abs(ref_fftData_FPNSub(1:end/2,:)))))); colormap(gray);

    % dispersion estimation
    % TODO: visualize image to find depthROI
    maxDispOrders = 3;
    coeffRange = 10;
    depthROI = [45,275];

    dispCoeffs = setDispCoeffs(ref_rawData_FPNSub, depthROI, maxDispOrders, coeffRange);
    ref_rawData_dispComp = compDisPhase(ref_rawData_FPNSub, maxDispOrders, dispCoeffs);
    ref_fftData_dispComp = fft(ref_rawData_dispComp);

    subplot(1,2,1); imagesc( imadjust(mat2gray(20 .* log10(...
         abs(ref_fftData(1:end/2,:)))))); colormap(gray);hold on;
    subplot(1,2,2); imagesc( imadjust(mat2gray(20 .* log10(...
         abs(ref_fftData_dispComp(1:end/2,:)))))); colormap(gray);

    % function volume_processing(rawOCT)
    %     for frame_num = 1:size(rawOCT, 3)
    %         raw_data = rawOCT(:,:,frame_num);

    %         procd_data(:,:,frame_num) = FFTData(1:size(raw_data, 1)/2,:);
    %     end
    % end
end