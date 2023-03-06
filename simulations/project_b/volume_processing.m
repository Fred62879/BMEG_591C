function procd_data = volume_processing(rawOCT, num_splits, depthIdx, depthROI, maxDispOrders, load_cache)
    % process data in batches
    num_frames = size(rawOCT, 3);
    lo = 1;
    hi = num_frames;
    bsz = 400;
    procd_data = zeros(...
        depthROI(2) - depthROI(1) + 1,...
        size(rawOCT,2),...
        num_frames*(num_splits + 1));

    if load_cache
        batches = round(num_frames*(num_splits+1) / bsz);
        for i=0:batches-1
            lo = i*bsz + 1;
            hi = min(lo + bsz - 1, num_frames*(num_splits+1));
            fname = append('../data/project_b/procd_', int2str(i+1), '.mat');
            batch_data = load(fname);
            procd_data(:,:,lo:hi) = getfield(batch_data,'batch_data');
        end
        return
    end

    for frame_num = lo:hi
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
        maxDispOrders = size(dispCoeffs,2) + 1;
        rawData_dispComp = compDisPhase(rawData_hanWin, maxDispOrders, dispCoeffs);
        fftData_dispComp = fft(rawData_dispComp);

        % collect full data
        id_full = (frame_num - 1)*(num_splits+1) + 1;
        procd_data(:,:,id_full) = fftData_dispComp(depthROI(1):depthROI(2),:);

        % split raw data and disp comp and collect
        rawData_hanWin_split = split(rawData_hanWin);
        for i=1:num_splits
            cur = fft(compDisPhase(rawData_hanWin_split(:,:,i), maxDispOrders, dispCoeffs));
            procd_data(:,:,id_full+i) = cur(depthROI(1):depthROI(2),:);
        end

        if mod(frame_num, 50) == 0
            disp(frame_num);
        end

        if mod(frame_num*(num_splits+1), bsz) == 0
            disp("saving data");
            i = frame_num*(num_splits+1) / bsz;
            batch_data = procd_data(:,:,(i-1)*bsz+1 : i*bsz);
            fname = append('../data/project_b/procd_', int2str(i), '.mat');
            save(fname,'batch_data');
        end
    end

    % visualize sequence
%     for i = 1:4:num_frames*(num_splits+1) %size(procd_data,3)
%         imagesc( imadjust(mat2gray(20 .* log10( ...
%             abs(procd_data(:,:,i)))))); colormap(gray);
%         pause(0.001);
%     end

    % load validate
%     for i = 1:4:4000 %size(procd_data,3)
%         imagesc( imadjust(mat2gray(20 .* log10( ...
%             abs(loaded_procd_data(:,:,i)))))); colormap(gray);
%         pause(0.001);
%     end
% 
%     for i=1:num_splits+1
%         subplot(1,num_splits+1,i); imagesc( imadjust(mat2gray(20 .* log10(...
%             abs(loaded_procd_data(:,:,i+36)))))); colormap(gray);
%     end
%     subplot(1,4,1); imagesc( imadjust(mat2gray(20 .* log10(...
%          abs(loaded_procd_data(:,:,1)))))); colormap(gray);
%     subplot(1,4,2); imagesc( imadjust(mat2gray(20 .* log10(...
%          abs(loaded_procd_data(:,:,49)))))); colormap(gray);
%     subplot(1,4,3); imagesc( imadjust(mat2gray(20 .* log10(...
%          abs(loaded_procd_data(:,:,97)))))); colormap(gray);
%     subplot(1,4,4); imagesc( imadjust(mat2gray(20 .* log10(...
%          abs(loaded_procd_data(:,:,117)))))); colormap(gray);
end