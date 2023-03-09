
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
    disp("loaded")

    % process reference data (and split)dispCoeffs
    dispCoeffs = -1;
    [depthIdx, depthROI, winFunc, ref_fftData_1D, dispCoeffs, ref_procd_Full, ref_procd_splits] = ...
        reference_processing(ref_RawData_Full, dispCoeffs, num_splits);
    disp("process ref")

    % visualize
%     subplot(1,num_splits+1,1); imagesc( imadjust(mat2gray(20 .* log10(...
%          abs(ref_procd_Full(:,:)))))); colormap(gray); title('full spectrum')
%     for i=1:num_splits
%         subplot(1,num_splits+1,i+1); imagesc( imadjust(mat2gray(20 .* log10(...
%             abs(ref_procd_splits(:,:,i) ))))); colormap(gray); title(strcat('spectrum', int2str(i)));
%     end

    % volume processing
    load_cache = true;
    maxDispOrders = size(dispCoeffs)+1;
    procd_data = volume_processing(rawOCT, num_splits, depthIdx, depthROI, maxDispOrders, load_cache);
    procd_OCT_BM_ROI = procd_data(:,:,1:4:end);
    %for i = 1:4:num_frames*repBScans*(num_splits+1)
%     for i = 1:num_frames*repBScans
%         imagesc( imadjust(mat2gray(20 .* log10( ...
%             abs(procd_OCT_BM_ROI(:,:,i)))))); colormap(gray);
%             %abs(procd_data(:,:,i)))))); colormap(gray);
%         pause(0.001);
%     end
    disp("loaded volume")

    % local motion (axial and lateral) correction
    % use procd_data_full (higher resolution) to estimate shift amount
    [xShift_local, yShift_local, cplxOCT_mcorr_local] = local_motion_correction(...
        procd_OCT_BM_ROI, -1, -1, true);

%     for i = 1:4:num_frames*repBScans
%         imagesc( imadjust(mat2gray(20 .* log10( ...
%             abs(cplxOCT_mcorr_local(:,:,i)))))); colormap(gray);
%         pause(0.001);
%     end
%     subplot(1,2,1); plot(xShift_local); title('axial'); % lateral
%     subplot(1,2,2); plot(yShift_local); title('lateral'); % lateral

    % local motion correction for the splitted data
    n = num_frames * repBScans;
    cplx_OCT_mcorr_local_split = zeros([size(procd_data,1) size(procd_data,2) num_splits*n]);
    for i=1:num_splits
        [~, ~, cur_mcorr_local] = local_motion_correction(...
            procd_data(:,:,i+1:4:end), xShift_local, yShift_local, false);
        cplx_OCT_mcorr_local_split(:,:,(i-1)*n+1:i*n) = cur_mcorr_local;
    end

    % for i = 1:n
    %     imagesc( imadjust(mat2gray(20 .* log10( ...
    %         abs(cplx_OCT_mcorr_local_split(:,:,i)))))); colormap(gray);
    %     pause(0.001);
    % end

    % octa (mean, var, sub, decorr)
    % bulk phase offset (do it for full but no need for ssada)
    [avgOCT, Var, Sub, Dec] = process_oct_a(repBScans, cplxOCT_mcorr_local);
    disp("got oct a")
    Dec_ssada = decorrelate_ssada(repBScans, num_frames, num_splits, cplx_OCT_mcorr_local_split);

    for i = 1:num_frames
        %imagesc( imadjust(mat2gray(20 .* log10(abs(Sub(:,:,i)))))); colormap(gray);
        imagesc( imadjust(mat2gray(abs(Dec(:,:,i))))); colormap(gray);
        pause(0.001);
    end

    % thresholding
    avgOCT_log = 20.*log10(abs(flipud(avgOCT)));

    % OCTA with intensity thresholding
    OCTA_Var = flipud(Var); OCTA_Var(avgOCT_log<90)=min(Var(:));
    OCTA_Sub = flipud(Sub); OCTA_Sub(avgOCT_log<90)=min(Sub(:));
    OCTA_Dec = flipud(Dec); OCTA_Dec(avgOCT_log<90)=min(Dec(:));
    OCTA_Dec_ssada = flipud(Dec_ssada); OCTA_Dec_ssada(avgOCT_log<90)=min(Dec_ssada(:));

    % smoothing
    avgOCTA_var = mov2DAvg(OCTA_Var, [2,2]);
    avgOCTA_sub = mov2DAvg(OCTA_Sub, [2,2]);
    avgOCTA_dec = mov2DAvg(OCTA_Dec, [2,2]);

    % figs = subplot(1,3,1);imagesc( imadjust(mat2gray(abs(avgOCTA_sub(:,:,1))))); colormap(gray);
    % figs = subplot(1,3,2);imagesc( imadjust(mat2gray(abs(avgOCTA_var(:,:,1))))); colormap(gray);
    % figs = subplot(1,3,3);imagesc( imadjust(mat2gray(abs(avgOCTA_dec(:,:,1))))); colormap(gray); saveas(figs, "b_scan.png")
    % return

    % for i = 1:num_frames
    %     imagesc( imadjust(mat2gray(abs(avgOCTA_var(:,:,i))))); colormap(gray);
    %     pause(0.001);
    % end

    % global motion correction
    usfac = 1;
    [global_mshift, global_tshift, avgOCT_corr] = volume_correction(usfac, avgOCT_log, -1, -1, true, false);
    [~, ~, varOCT_corr] = volume_correction(usfac, avgOCTA_var, global_mshift, global_tshift, false, false);
    [~, ~, subOCT_corr] = volume_correction(usfac, avgOCTA_sub, global_mshift, global_tshift, false, false);
    [~, ~, decOCT_corr] = volume_correction(usfac, avgOCTA_dec, global_mshift, global_tshift, false, false);
    [~, ~, decOCT_corr_ssada] = volume_correction(usfac, OCTA_Dec_ssada, global_mshift, global_tshift, false, false);

    % save locally as tiff
    save_tiff(varOCT_corr, "../data/project_b/procd_OCT_corr_var.tiff");
    save_tiff(subOCT_corr, "../data/project_b/procd_OCT_corr_sub.tiff");
    save_tiff(decOCT_corr, "../data/project_b/procd_OCT_corr_dec.tiff");
    save_tiff(decOCT_corr_ssada, "../data/project_b/procd_OCT_corr_dec_ssada.tiff");

    % en face image
    % frames = [50 60 70 80 90 100 110 120 160 200];
    % frames = 70:2:80;
    % n = size(frames,2);
    % for i=1:n
    %     frame = frames(1,i);
    %     subplot(1,n,i); imagesc( imadjust(mat2gray(20 .* log10( ...
    %         abs(squeeze(avgOCTA_sub(frame,:,:))))))) ; colormap(gray); title('en face'); axis xy;
    %     subplot(1,n,i); imagesc( imadjust(mat2gray(...
    %         abs(squeeze(avgOCTA_sub(frame,:,:)))))); colormap(gray); title('en face'); axis xy;
    % end
end
