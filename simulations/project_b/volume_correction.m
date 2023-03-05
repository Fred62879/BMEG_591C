function [axialShift, OCT_corr] = volume_correction(OCT_data, axialShift, estimate_shift)

    axialShift = zeros(numFrames, 1);

    OCT_mcorr = OCT_data;

    for I = 1:numFrames
        if estimate_shift
            [output, ~] = dftregistration(...
                fft2( 20 .* log10(abs(OCT_data(:,:,round(numFrames)./2)))),...
                fft2(20 .* log10(abs(OCT_data(:,:,I)))),...
                usfac);
            axialShift(I) = round(output(3));
        end

        OCT_mcorr(:,:,I) = circshift(OCT_data(:,:,I), [axialShift(I), 0]);
    end

    %[axialShift procd_OCT_ROI_mcorr] = motion_correction(procd_OCT_ROI);
    %numLines = size(procd_OCT_ROI_mcorr, 2);

    % apply the fitted tilt estimation
    axialshift = tilt_estimation(procd_OCT_ROI_mcorr);
    plot(axialshift);

    % fit a curve to the axialshift curve
    x = [1:numLines]';
    p = polyfit(x,axialshift,2);

    x1 = linspace(0,500,500);
    y1 = polyval(p,x1);
    subplot(1,2,1); plot(axialshift); xlim([0 500]); ylim([-30 10]); title('Tilting Slope Estimation');
    subplot(1,2,2); plot(x1,y1); title('Interpolated Tilt Correction Parameter'); % y1 is the smoothed axialshift

    % tilt correction
    procd_OCT_ROI_tcorr = tilt_correction(procd_OCT_ROI_mcorr, y1);

    subplot(1,2,1); imagesc( imadjust(mat2gray(20 .* log10( ...
            abs(squeeze(procd_OCT_ROI_mcorr(:,:,5))))))) ; colormap(gray);
    subplot(1,2,2); imagesc( imadjust(mat2gray(20 .* log10( ...
            abs(squeeze(procd_OCT_ROI_tcorr(:,:,5))))))) ; colormap(gray);

    % fast and slow scan direction
    subplot(1,2,1); imagesc( imadjust(mat2gray(20 .* log10( ...
            abs(squeeze(procd_OCT_ROI_tcorr(:,:,1))))))) ; colormap(gray); title('Fast Scan Direction'); axis xy;
    subplot(1,2,2); imagesc( imadjust(mat2gray(20 .* log10( ...
            abs(squeeze(procd_OCT_ROI_tcorr(:,1,:))))))) ; colormap(gray); title('Slow Scan Direction'); axis xy;

    % en face image
    frames = [50 60 70 80 90 100 110 120 160 200 250];
    frames = 110:1:120;
    n = size(frames,2);
    for i=1:n
        frame = frames(1,i);
        subplot(1,n,i); imagesc( imadjust(mat2gray(20 .* log10( ...
            abs(squeeze(procd_OCT_ROI_tcorr(frame,:,:))))))) ; colormap(gray); title('en face'); axis xy;
    end
end