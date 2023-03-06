
function angiography

    usfac = 1;
    num_frames = 500;
    repBScans = 2;
    num_splits = 3;

    % load raw data
    % fn = '/media/fred/working_drive/2022W2/bmeg591c/simulations/data/project_b/RawOCT_BM.mat';
    fn = '/scratch/academic/2022w2/bmeg591c/simulations/data/project_b/RawOCT_BM.mat';
    raw_data = load(fn);
    rawOCT = getfield(raw_data,'rawOCT_BM');
    ref_RawData_Full = rawOCT(:,:,500);

    % process reference data (and split)dispCoeffs
    dispCoeffs = -1;
    [depthIdx, depthROI, winFunc, ref_fftData_1D, dispCoeffs, ref_procd_Full, ref_procd_splits] = ...
        reference_processing(ref_RawData_Full, dispCoeffs, num_splits);

    % visualize
    subplot(1,num_splits+1,1); imagesc( imadjust(mat2gray(20 .* log10(...
         abs(ref_procd_Full(:,:)))))); colormap(gray);
    for i=1:num_splits
        subplot(1,num_splits+1,i+1); imagesc( imadjust(mat2gray(20 .* log10(...
            abs(ref_procd_splits(:,:,i) ))))); colormap(gray);
    end

    % volume processing
    load_cache = true;
    maxDispOrders = size(dispCoeffs)+1;
    procd_data = volume_processing(rawOCT, num_splits, depthIdx, depthROI, maxDispOrders, load_cache);
    procd_OCT_BM_ROI = procd_data(:,:,1:4:end);

    % local motion (axial and lateral) correction
    % use procd_data_full (higher resolution) to estimate shift amount
    [xShift_local, yShift_local, cplxOCT_mcorr_local] = local_motion_correction(...
        procd_OCT_BM_ROI, -1, -1, true);
%     for i = 1:4:num_frames*repBScans
%         imagesc( imadjust(mat2gray(20 .* log10( ...
%             abs(cplxOCT_mcorr_local(:,:,i)))))); colormap(gray);
%         pause(0.001);
%     end
    %subplot(1,2,1); plot(xShift_local); ylim([-0.2 0.2]); % lateral
    %subplot(1,2,2); plot(yShift_local); ylim([-0.2 0.2]); % lateral

    % local motion correction for the splitted data
    n = num_frames * repBScans;
    cplx_OCT_mcorr_local_split = zeros([size(procd_data,1) size(procd_data,2) num_splits*n]);
    for i=1:num_splits
        [~, ~, cur_mcorr_local] = local_motion_correction(...
            procd_data(:,:,i+1:4:end), xShift_local, yShift_local, false);
        cplx_OCT_mcorr_local_split(:,:,(i-1)*n+1:i*n) = cur_mcorr_local;
    end
%     for i = 1:4:num_frames*repBScans
%         imagesc( imadjust(mat2gray(20 .* log10( ...
%             abs(cplx_OCT_mcorr_local_split(:,:,i+1)))))); colormap(gray);
%         pause(0.001);
%     end

    % octa (mean, var, sub, decorr)
    % bulk phase offset (do it for full but no need for ssada)
    [avgOCT, Var, Sub, Dec] = process_oct_a(repBScans, cplxOCT_mcorr_local);
    Dec_ssada = decorrelate_ssada(repBScans, num_frames, num_splits, cplx_OCT_mcorr_local_split);
    for i = 1:num_frames
        imagesc( imadjust(mat2gray(20 .* log10( ...
            abs(Dec_ssada(:,:,i)))))); colormap(gray);
        pause(0.001);
    end

    % thresholding
    avgOCT_log = 20.*log10(abs(flipud(avgOCT)));

    % OCTA with intensity thresholding
    OCTA_Var = flipud(Var); OCTA_Var(avgOCT_log<90)=min(Var(:));
    OCTA_Sub = flipud(Var); OCTA_Sub(avgOCT_log<90)=min(Sub(:));
    OCTA_Dec = flipud(Var); OCTA_Dec(avgOCT_log<90)=min(Dec(:));
    OCTA_Dec_ssada = flipud(Var); OCTA_Dec_ssada(avgOCT_log<90)=min(Dec_ssada(:));

    % smoothing
    avgOCTA_var = mov2DAvg(OCTA_Var, [2,2]);
    avgOCTA_sub = mov2DAvg(OCTA_Sub, [2,2]);
    avgOCTA_dec = mov2DAvg(OCTA_Dec, [2,2]);

    % global motion correction
    usfac = 1;
    avgOCT_corr = volume_correction(usfac, avgOCT_log, -1, true);
    varOCT_corr = volume_correction(usfac, avgOCTA_var, -1, true);
    subOCT_corr = volume_correction(usfac, avgOCTA_sub, -1, false);
    decOCT_corr = volume_correction(usfac, avgOCTA_dec, -1, false);
    decOCT_corr_ssada(:,:,I) = volume_correction(OCTA_Dec_ssada(:,:,I), -1, false);

    % save locally as tiff
    scaled = imadjust(mat2gray(20 .* log10(abs(squeeze(varOCT_corr(1,:,:))))));
    imwrite(scaled,"../data/project_b/procd_OCT_corr_var.tiff");

    for i = 2:size(varOCT_corr,1)
        scaled = imadjust(mat2gray(20 .* log10(abs(squeeze(varOCT_corr(i,:,:))))));
        imwrite(scaled,"../data/project_b/procd_OCT_corr_var.tiff","WriteMode","append")
    end

    % en face image
    frames = [50 60 70 80 90 100 110 120 160 200 250];
    frames = 70:2:80;
    n = size(frames,2);
    for i=1:n
        frame = frames(1,i);
        subplot(1,n,i); imagesc( imadjust(mat2gray(20 .* log10( ...
            abs(squeeze(avgOCTA_sub(frame,:,:))))))) ; colormap(gray); title('en face'); axis xy;
    end
end