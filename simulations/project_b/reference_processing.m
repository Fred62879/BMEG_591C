function [depthIdx, depthROI, winFunc, ref_fftData_1D, dispCoeffs, ref_fftData_dispComp, ref_fftData_dispComp_splits] = ...
    reference_processing(ref_rawData, dispCoeffs, num_splits)

    ref_fftData = fft(hilbert(ref_rawData));
    
    %subplot(1,2,1); plot(abs(ref_fftData(1:end/2,:))); ylim([0 4e6]);
    %subplot(1,2,2); plot(abs(ref_fftData(1:end/2,1))); ylim([0 4e6]);
    %imagesc( imadjust(mat2gray(20 .* log10(abs(ref_fftData(1:end/2,:)))))); colormap(gray);

    % binary filter
    calSigMin = [29,40];
    winFunc = zeros(size(ref_fftData));
    winFunc(calSigMin(1):calSigMin(2),:) = 1;
    ref_cal_fftData = ref_fftData .* winFunc;
    %imagesc( imadjust(mat2gray(20 .* log10(abs(ref_cal_fftData(1:end/2,:)))))); colormap(gray);

    % get calibrated raw data
    ref_cal_rawData = ifft(ref_cal_fftData);
    %plot(real(ref_cal_rawData(1:end/2,:)));

    % resampling
    ref_rawData_rescaled = reSampling_CalSig(ref_fftData, ref_cal_rawData);
    ref_fftData_rescaled = fft(ref_rawData_rescaled);
    %imagesc( imadjust(mat2gray(20 .* log10(abs(ref_fftData_rescaled(1:end/2,:)))))); colormap(gray);    

    % spectra shift / phase compensation
    cplxConjX = ref_fftData_rescaled .* ...
        repmat( conj(ref_fftData_rescaled(:,1)), [1 size(ref_fftData_rescaled,2)] );
    %plot(abs(cplxConjX(:,1:10))); xlim([ 10 100]) % find optimal depthIdx

    depthIdx = 34;
    ref_fftData_1D = ref_fftData_rescaled(:,1);
    ref_rawData_PhaseComp = compPhaseShift(ref_fftData_1D, ref_fftData_rescaled, depthIdx);
    %ref_fftData_PhaseComp = fft(ref_rawData_PhaseComp);
    %imagesc( imadjust(mat2gray(20 .* log10(abs(ref_fftData_PhaseComp(1:end/2,:)))))); colormap(gray);

    % Fixed pattern noise (FPN) removal
    ref_rawData_FPNSub = ref_rawData_PhaseComp...
        - (repmat(median(real(ref_rawData_PhaseComp), 2), [1,size(ref_rawData_PhaseComp,2)]) ...
           +1j .* repmat(median(imag(ref_rawData_PhaseComp), 2), [1, size(ref_rawData_PhaseComp, 2)]));
    ref_fftData_FPNSub = fft(ref_rawData_FPNSub);
    %imagesc( imadjust(mat2gray(20 .* log10(abs(ref_fftData_FPNSub(1:end/2,:)))))); colormap(gray);

    % hanning windowing
    ref_rawData_hanWin = ref_rawData_FPNSub ...
        .* repmat(hann(size(ref_rawData_FPNSub,1)), [1 size(ref_rawData_FPNSub,2)]);
    ref_fftData_hanWin = fft(ref_rawData_hanWin);

    % split raw data
    ref_rawData_hanWin_split = split(ref_rawData_hanWin);

    % visualize fringe data and image of splitted signal
    subplot(1,num_splits+1,1); plot(real(ref_rawData_hanWin));title('full spectrum OCT fringe data')
    for i=1:num_splits
       subplot(1,num_splits+1,i+1); plot(real(ref_rawData_hanWin_split(:,:,i))); ylim([-5000 5000]); ...
           title(strcat('split spectrum OCT fringe data - ', int2str(i)))
    end

    %subplot(1,num_splits+1,1); imagesc( imadjust(mat2gray(20 .* log10(...
    %     abs(ref_fftData_hanWin(1:end/2,:)))))); colormap(gray);
    %for i=1:num_splits
    %    cur_ref_fftData_split = fft(ref_rawData_hanWin_split(:,:,i));
    %    subplot(1,num_splits+1,i+1); imagesc( imadjust(mat2gray(20 .* log10(...
    %        abs(cur_ref_fftData_split(1:end/2,:)))))); colormap(gray);
    %end

    % dispersion estimation
    maxDispOrders = 5;
    coeffRange = 10;
    depthROI = [40,250]; % guessed, removed dark region (use any image before)

    if dispCoeffs == -1
        dispCoeffs = setDispCoeffs(ref_rawData_hanWin, depthROI, maxDispOrders, coeffRange);
    end
    ref_rawData_dispComp = compDisPhase(ref_rawData_hanWin, maxDispOrders, dispCoeffs);
    ref_fftData_dispComp = fft(ref_rawData_dispComp);

    % disp comp splitted signal
    ref_fftData_dispComp_splits = zeros([depthROI(2)-depthROI(1)+1 size(ref_fftData_FPNSub,2) num_splits]);
    for i=1:num_splits
        % calculate disp coeff for each split individually
        %if dispCoeffs_split(i,1) == -1
        %    dispCoeffs_split(i,:) = setDispCoeffs(ref_rawData_hanWin_split(:,:,i), depthROI, maxDispOrders, coeffRange);
        %end
        %cur_ref_rawData_dispComp = compDisPhase(ref_rawData_hanWin_split(:,:,i), maxDispOrders, dispCoeffs_split(i,:));

        % use disp coeff from full data for compensation
        cur = fft(compDisPhase(ref_rawData_hanWin_split(:,:,i), maxDispOrders, dispCoeffs));
        ref_fftData_dispComp_splits(:,:,i) = cur(depthROI(1):depthROI(2),:);
    end

    ref_fftData_dispComp = ref_fftData_dispComp(depthROI(1):depthROI(2),:);
end
