function [motion_axialshift, tilt_axialshift, OCT_tcorr] = volume_correction(...
    usfac, OCT_data, motion_axialshift, tilt_axialshift, estimate_shift, plot_figure)

    num_frames = size(OCT_data,3);
    OCT_mcorr = OCT_data;

    % fast b scan
%     for i = 1:num_frames
%         imagesc( imadjust(mat2gray(20 .* log10( ...
%             abs(OCT_data(:,:,i)))))); colormap(gray);
%         pause(0.001);
%     end

    if estimate_shift
        motion_axialshift = zeros(num_frames, 1);
    end
    for I = 1:num_frames
        if estimate_shift
            [output, ~] = dftregistration(...
                fft2(20 .* log10(abs(OCT_data(:,:,round(num_frames)./2)))),...
                fft2(20 .* log10(abs(OCT_data(:,:,I)))),...
                usfac);
            motion_axialshift(I) = round(output(3));
        end

        OCT_mcorr(:,:,I) = circshift(OCT_data(:,:,I), [motion_axialshift(I), 0]);
    end

    % slow b scan (very bumpy vs smoothed)
    if plot_figure
        subplot(1,2,1); imagesc( imadjust(mat2gray(20 .* log10( ...
            abs(squeeze(OCT_data(:,250,:))))))) ; colormap(gray);
        subplot(1,2,2); imagesc( imadjust(mat2gray(20 .* log10( ...
            abs(squeeze(OCT_mcorr(:,250,:))))))) ; colormap(gray);
    end

    % estimate tile axial shift
    numLines = size(OCT_mcorr, 2);
    if estimate_shift
        tilt_axialshift = zeros(1, num_frames);
        for I = 1:numLines
            [output,~] = dftregistration( fft2(squeeze(OCT_mcorr(:,round(numLines./2),:) )), ...
                                          fft2(squeeze(OCT_mcorr(:,I,:))), usfac );
            tilt_axialshift(I) = round(output(3));
        end
        % plot(tilt_axialshift);

        % fit a curve to the axialshift curve
        numLines = size(OCT_mcorr, 2);
        x = [1:numLines]';
        p = polyfit(x,tilt_axialshift,2);
    
        x1 = linspace(0,500,500);
        y1 = polyval(p,x1);
        tilt_axialshift = y1;
        % subplot(1,2,1); plot(tilt_axialshift); title('Tilting Slope Estimation');
        % subplot(1,2,2); plot(x1,y1); title('Interpolated Tilt Correction Parameter'); % y1 is the smoothed axialshift
    end

    % tilt correction
    OCT_tcorr = zeros(size(OCT_mcorr));
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