
function angiography

    procd_OCT_BM_ROI = load(); % no motion or tilt correction

    usfac = 1;
    numFrames = size(procd_OCT_BM_ROI, 3);
    repBScans = 2;

    for i = 1:size(procd_OCT_BM_ROI)
        imagesc( :,:,i);
        pause(0.01);
    end

    xShift_local = zeros([numFrames 1]);
    yShift_local = zeros([numFrames 1]);

    % local motion (axial and lateral) correction
    cplxOCT_mcorr_local = local_motion_correction(procd_OCT_BM_ROI);

    plot(xShift_local) % lateral

    % oct average & octa
    [Var, Sub, Dec] = process_oct_a(cplxOCT_mcorr_local);

    % thresholding
    avgOCT_log = 20.*log10(abs(flipud(avgOCT)));

    % OCTA with intensity thresholding
    OCTA_Var = flipud(Var); OCTA_Var(avgOCT_log<90)=min(Var(:));
    OCTA_Sub = flipud(Var); OCTA_Sub(avgOCT_log<90)=min(Sub(:));
    OCTA_Dec = flipud(Var); OCTA_Dec(avgOCT_log<90)=min(Dec(:));

    % smoothing & global motion correction
    avgOCTA_var = mov2DAvg(OCTA_Var, [2,2]);
    avgOCTA_sub = mov2DAvg(OCTA_Sub, [2,2]);
    avgOCTA_dec = mov2DAvg(OCTA_Dec, [2,2]);

    % global motion correction
    for I = 1:numFrames
        [output,~] = dftregistration(fft2(avgOCT_log(:,:,numFrames/2)),...
                                     fft2(avgOCT_log(:,:,I)),usfac);

        yShift_global(I) = round(output(3));
        avgOCT_mcorr(:,:,I) = circshift(avgOCT_log(:,:,I), [output(3),0]);
        varOCT_mcorr(:,:,I) = circshift(avgOCT_var(:,:,I), [output(3),0]);
        subOCT_mcorr(:,:,I) = circshift(avgOCT_sub(:,:,I), [output(3),0]);
        decOCT_mcorr(:,:,I) = circshift(avgOCT_dec(:,:,I), [output(3),0]);
    end

    % tilt correction

    % ssada processing
end