
function ref_angiography

    % load procd data
    bsz = 100;
    num_frames = 1000;
    batches = num_frames/bsz;
    loaded_procd_data = zeros(260,500,num_frames);
    for i=0:batches-1
        lo = i*bsz + 1;
        hi = min(lo + bsz - 1, num_frames);
        fname = append('../data/project_b/procd_', int2str(i+1), '.mat');
        batch_data = load(fname);
        loaded_procd_data(:,:,lo:hi) = getfield(batch_data,'batch_data');
    end

    % visualize
    usfac = 1;
    repBScans = 2;
    procd_OCT_BM_ROI = loaded_procd_data; % no motion or tilt correction

    for i = 1:30 %size(procd_OCT_BM_ROI,3)
        imagesc( imadjust(mat2gray(20 .* log10( ...
            abs(loaded_procd_data(:,:,i)))))); colormap(gray);
        pause(0.01);
    end

    % load raw data
    fn = '/media/fred/working_drive/2022W2/bmeg591c/simulations/data/project_b/RawOCT_BM.mat';
    raw_data = load(fn);
    rawOCT = getfield(raw_data,'rawOCT_BM');
    ref_RawData_Full = rawOCT(:,:,250);

    % process reference data (and split)dispCoeffs
    dispCoeffs = -1;
    num_splits = 3;
    [depthIdx, depthROI, winFunc, ref_fftData_1D, dispCoeffs, ref_procd_Full, ref_procd_splits] = ...
        reference_processing(ref_RawData_Full, dispCoeffs, num_splits);

    % visualize
    subplot(1,num_splits+1,1); imagesc( imadjust(mat2gray(20 .* log10(...
         abs(ref_procd_Full(:,:)))))); colormap(gray);
    for i=1:num_splits
        subplot(1,num_splits+1,i+1); imagesc( imadjust(mat2gray(20 .* log10(...
            abs(ref_procd_splits(:,:,i) ))))); colormap(gray);
    end

    % local motion (axial and lateral) correction
    % use procd_data_full (higher resolution) to estimate shift amount
    xShift_local = zeros([num_frames 1]);
    yShift_local = zeros([num_frames 1]);
    cplxOCT_mcorr_local = local_motion_correction(procd_OCT_BM_ROI);

    plot(xShift_local) % lateral

    % octa (mean, var, sub, decorr)
    % oct average & octa
    [Var, Sub, Dec] = process_oct_a(cplxOCT_mcorr_local);

    % bulk phase offset (do it for full but do we do it for the 3 split?)


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

    % ssada processing (3.2.3)
end