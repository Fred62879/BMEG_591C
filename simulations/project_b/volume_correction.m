function [axialShift, OCT_tcorr] = volume_correction(usfac, OCT_data, axialShift, estimate_shift)

    num_frames = size(OCT_data,3);
    axialShift = zeros(num_frames, 1);

    OCT_mcorr = OCT_data;

    % fast b scan
%     for i = 1:num_frames
%         imagesc( imadjust(mat2gray(20 .* log10( ...
%             abs(OCT_data(:,:,i)))))); colormap(gray);
%         pause(0.001);
%     end

    for I = 1:num_frames
        if estimate_shift
            [output, ~] = dftregistration(...
                fft2( 20 .* log10(abs(OCT_data(:,:,round(num_frames)./2)))),...
                fft2(20 .* log10(abs(OCT_data(:,:,I)))),...
                usfac);
            axialShift(I) = round(output(3));
        end

        OCT_mcorr(:,:,I) = circshift(OCT_data(:,:,I), [axialShift(I), 0]);
    end

    % slow b scan (very bumpy vs smoothed)
%     subplot(1,2,1); imagesc( imadjust(mat2gray(20 .* log10( ...
%             abs(squeeze(OCT_data(:,250,:))))))) ; colormap(gray);
%     subplot(1,2,2); imagesc( imadjust(mat2gray(20 .* log10( ...
%             abs(squeeze(OCT_mcorr(:,250,:))))))) ; colormap(gray);

    % estimate tile axial shift
    numLines = size(OCT_mcorr, 2);
    tilt_axialshift = zeros(1, num_frames);
    for I = 1:numLines
        [output,~] = dftregistration( fft2(squeeze(OCT_mcorr(:,round(numLines./2),:) )), ...
                                      fft2(squeeze(OCT_mcorr(:,I,:))), usfac );
        tilt_axialshift(I) = round(output(3));
    end
%     plot(tilt_axialshift);

    % fit a curve to the axialshift curve
    numLines = size(OCT_mcorr, 2);
    x = [1:numLines]';
    p = polyfit(x,tilt_axialshift,2);

    x1 = linspace(0,500,500);
    y1 = polyval(p,x1);
%     subplot(1,2,1); plot(tilt_axialshift); title('Tilting Slope Estimation');
%     subplot(1,2,2); plot(x1,y1); title('Interpolated Tilt Correction Parameter'); % y1 is the smoothed axialshift

    % tilt correction
    OCT_tocorr = zeros(size(OCT_mcorr));
    for I = 1:numLines
        OCT_tcorr(:,I,:) = circshift(OCT_mcorr(:,I,:), [round(tilt_axialshift(1,I)), 0]);
    end

%     subplot(1,2,1); imagesc( imadjust(mat2gray(20 .* log10( ...
%             abs(squeeze(OCT_mcorr(:,:,5))))))) ; colormap(gray);
%     subplot(1,2,2); imagesc( imadjust(mat2gray(20 .* log10( ...
%             abs(squeeze(OCT_tcorr(:,:,5))))))) ; colormap(gray);
% 
%     % fast and slow scan direction
%     subplot(1,2,1); imagesc( imadjust(mat2gray(20 .* log10( ...
%             abs(squeeze(procd_OCT_ROI_tcorr(:,:,1))))))) ; colormap(gray); title('Fast Scan Direction'); axis xy;
%     subplot(1,2,2); imagesc( imadjust(mat2gray(20 .* log10( ...
%             abs(squeeze(procd_OCT_ROI_tcorr(:,1,:))))))) ; colormap(gray); title('Slow Scan Direction'); axis xy;

    % en face image
%     frames = [50 60 70 80 90 100 110 120 160 200 250];
%     frames = 110:1:120;
%     n = size(frames,2);
%     for i=1:n
%         frame = frames(1,i);
%         subplot(1,n,i); imagesc( imadjust(mat2gray(20 .* log10( ...
%             abs(squeeze(procd_OCT_ROI_tcorr(frame,:,:))))))) ; colormap(gray); title('en face'); axis xy;
%     end
end