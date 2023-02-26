function volume_processing()
    % process data in batches
    procd_data = zeros(260,500,500);
    lo = 1;
    hi = 500;

    for frame_num = lo:hi %size(rawOCT, 3)
        if mod(frame_num, 2) == 0
            disp(frame_num);
        end
        
        raw_data = rawOCT(:,:,frame_num);
        fft_Data = fft(hilbert(raw_data));

        % binary windowing
        cal_fftData = fft_Data .* winFunc;

        % get calibrated raw data
        cal_rawData = ifft(cal_fftData);

        % resampling
        rawData_rescaled = reSampling_CalSig(fft_Data, cal_rawData);
        fftData_rescaled = fft(rawData_rescaled);

        % phase compensation (can used the ref)
        rawData_PhaseComp = compPhaseShift(ref_fftData_1D, fftData_rescaled, depthIdx);

        % Fixed pattern noise (FPN) removal
        rawData_FPNSub = rawData_PhaseComp...
            - (repmat(median(real(rawData_PhaseComp), 2), [1,size(rawData_PhaseComp,2)]) ...
               +1j .* repmat(median(imag(rawData_PhaseComp), 2), [1, size(rawData_PhaseComp, 2)]));

        % hanning windowing
        rawData_hanWin = rawData_FPNSub ...
            .* repmat(hann(size(rawData_FPNSub,1)), [1 size(rawData_FPNSub,2)]);

        % dispersion compensation
        rawData_dispComp = compDisPhase(rawData_hanWin, maxDispOrders, dispCoeffs);
        fftData_dispComp = fft(rawData_dispComp);

        % collect result
        procd_data(:,:,frame_num) = fftData_dispComp(40:299,:);
        % procd_data(:,:,frame_num) = fftData_dispComp(1:size(raw_data, 1)/2,:);
    end

    subplot(1,2,1); imagesc( imadjust(mat2gray(20 .* log10(...
         abs(procd_data(:,:,1)))))); colormap(gray);
    subplot(1,2,2); imagesc( imadjust(mat2gray(20 .* log10(...
         abs(procd_data(:,:,2)))))); colormap(gray);

    % save data
    batches = 5;
    bsz = round(500 / batches);
    
    for i=0:batches-1
        lo = i*bsz+1
        hi = min(lo + bsz - 1,500)
        batch_data = procd_data(:,:,lo:hi);
        fname = append('../procd_', int2str(i), '.mat');
        save(fname,'batch_data');
    end

    % load data
    loaded_procd_data = zeros(260,500,500);
    for i=0:batches-1
        lo = i*bsz + 1;
        hi = min(lo + bsz - 1, 500);
        fname = append('../procd_', int2str(i), '.mat');
        batch_data = load(fname);
        loaded_procd_data(:,:,lo:hi) = getfield(batch_data,'batch_data');
    end
end