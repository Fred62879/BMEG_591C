function volume_correction()
    procd_OCT_ROI = procd_data(:,:,:);

    % fast b scan
    for i = 1:size(procd_OCT_ROI, 3)
        imagesc( imadjust(mat2gray(20 .* log10( ...
            abs(procd_data(:,:,1)))))); colormap(gray);
        pause(0.1);
    end

    % slow b scan (very bumpy)
    imagesc( imadjust(mat2gray(20 .* log10( ...
            abs(squeeze(procd_OCT_ROI(:,250,:))))))) ; colormap(gray);

    % motion correction
    procd_OCT_ROI_mcorr = motion_correction(procd_OCT_ROI);
    numLines = size(procd_OCT_ROI_mcorr, 2);
    subplot(1,2,1); imagesc( imadjust(mat2gray(20 .* log10( ...
            abs(squeeze(procd_OCT_ROI(:,250,:))))))) ; colormap(gray);
    subplot(1,2,2); imagesc( imadjust(mat2gray(20 .* log10( ...
            abs(squeeze(procd_OCT_ROI_mcorr(:,250,:))))))) ; colormap(gray);

    % apply the fitted tilt estimation
    axialshift = tilt_estimation(procd_OCT_ROI_mcorr);
    plot(axialshift);

    % fit a curve to the axialshift curve
    x = [1:numLines]';
    p = polyfit(x,axialshift,2);

    x1 = linspace(0,500,500);
    y1 = polyval(p,x1);
    subplot(1,2,1); plot(axialshift);
    subplot(1,2,2); plot(x1,y1); % y1 is the smoothed axialshift

    % tilt correction
    procd_OCT_ROI_tcorr = tilt_correction(procd_OCT_ROI_mcorr, y1);

    subplot(1,2,1); imagesc( imadjust(mat2gray(20 .* log10( ...
            abs(squeeze(procd_OCT_ROI_mcorr(:,:,5))))))) ; colormap(gray);
    subplot(1,2,2); imagesc( imadjust(mat2gray(20 .* log10( ...
            abs(squeeze(procd_OCT_ROI_tcorr(:,:,5))))))) ; colormap(gray);

    % save locally
    

end