function volume_correction()
    % load
%     fn = '/media/fred/working_drive/2022W2/bmeg591c/simulations/data/procd_data.mat';
%     procd_data = load(fn);
%     procd_OCT_ROI = getfield(procd_data,'procd_data');

    % or run on the fly
    procd_OCT_ROI = procd_data(:,:,:);

    % fast b scan
    for i = 1:2 %size(procd_OCT_ROI, 3)
        imagesc( imadjust(mat2gray(20 .* log10( ...
            abs(procd_data(:,:,1)))))); colormap(gray);
        pause(0.1);
    end

    % slow b scan (very bumpy)
    imagesc( imadjust(mat2gray(20 .* log10( ...
            abs(squeeze(procd_OCT_ROI(:,250,:))))))) ; colormap(gray);

    % motion correction
    [axialShift procd_OCT_ROI_mcorr] = motion_correction(procd_OCT_ROI);
    numLines = size(procd_OCT_ROI_mcorr, 2);
    subplot(1,2,1); imagesc( imadjust(mat2gray(20 .* log10( ...
            abs(squeeze(procd_OCT_ROI(:,250,:))))))) ; colormap(gray);
    subplot(1,2,2); imagesc( imadjust(mat2gray(20 .* log10( ...
            abs(squeeze(procd_OCT_ROI_mcorr(:,250,:))))))) ; colormap(gray);
    plot(axialShift); xlim([0 500]); title('Axial Motion Estimation');

    % apply the fitted tilt estimation
    axialshift = tilt_estimation(procd_OCT_ROI_mcorr);

    % fit a curve to the axialshift curve
    x = [1:numLines]';
    p = polyfit(x,axialshift,2);

    x1 = linspace(0,500,500);
    y1 = polyval(p,x1);
    subplot(1,2,1); plot(axialshift); xlim([0 500]); ylim([-30 10]); title('Tilting Slope Estimation');
    subplot(1,2,2); plot(x1,y1); title('Interpolated Tilt Correction Parameter'); % y1 is the smoothed axialshift

    % tilt correction
    procd_OCT_ROI_tcorr = tilt_correction(procd_OCT_ROI_mcorr, y1);

    % fast and slow scan direction
    subplot(1,2,1); imagesc( imadjust(mat2gray(20 .* log10( ...
            abs(squeeze(procd_OCT_ROI_tcorr(:,:,1))))))) ; colormap(gray); title('Fast Scan Direction'); axis xy;
    subplot(1,2,2); imagesc( imadjust(mat2gray(20 .* log10( ...
            abs(squeeze(procd_OCT_ROI_tcorr(:,1,:))))))) ; colormap(gray); title('Slow Scan Direction'); axis xy;

    % save locally as tiff
    scaled = imadjust(mat2gray(20 .* log10(abs(squeeze(procd_OCT_ROI_tcorr(1,:,:))))));
    imwrite(scaled,"../data/procd_OCT_ROI_tcorr.tiff");

    for i = 2:size(procd_OCT_ROI_tcorr,1)
        scaled = imadjust(mat2gray(20 .* log10(abs(squeeze(procd_OCT_ROI_tcorr(i,:,:))))));
        imwrite(scaled,"../data/procd_OCT_ROI_tcorr.tiff","WriteMode","append")
    end

end