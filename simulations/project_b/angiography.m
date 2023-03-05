
function angiography

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
    num_frames = 500;
    repBScans = 2;
    procd_OCT_BM_ROI = loaded_procd_data; % no motion or tilt correction

    for i = 1:30 %size(procd_OCT_BM_ROI,3)
        imagesc( imadjust(mat2gray(20 .* log10( ...
            abs(loaded_procd_data(:,:,i)))))); colormap(gray);
        pause(0.01);
    end

    % load raw data
    % fn = '/media/fred/working_drive/2022W2/bmeg591c/simulations/data/project_b/RawOCT_BM.mat';
    fn = '/scratch/academic/2022w2/bmeg591c/simulations/data/project_b/RawOCT_BM.mat';
    raw_data = load(fn);
    rawOCT = getfield(raw_data,'rawOCT_BM');
    ref_RawData_Full = rawOCT(:,:,500);

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

    % volume processing
    load = true;
    procd_data = volume_processing(rawOCT, num_splits, depthIdx, depthROI, maxDispOrders, load);
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
    avgOCT_mcorr = volume_correction(avgOCT_log, -1, true);
    varOCT_mcorr = volume_correction(avgOCTA_var, -1, false);
    subOCT_mcorr = volume_correction(avgOCTA_sub, -1, false);
    decOCT_mcorr = volume_correction(avgOCTA_dec, -1, false);
    decOCT_mcorr_ssada(:,:,I) = volume_correction(avgOCTA_Dec_ssada(:,:,I), -1, false);

    % save locally as tiff
    scaled = imadjust(mat2gray(20 .* log10(abs(squeeze(procd_OCT_ROI_tcorr(1,:,:))))));
    imwrite(scaled,"../data/procd_OCT_ROI_tcorr.tiff");

    for i = 2:size(procd_OCT_ROI_tcorr,1)
        scaled = imadjust(mat2gray(20 .* log10(abs(squeeze(procd_OCT_ROI_tcorr(i,:,:))))));
        imwrite(scaled,"../data/procd_OCT_ROI_tcorr.tiff","WriteMode","append")
    end

end