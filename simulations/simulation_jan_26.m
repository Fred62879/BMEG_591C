
function simulation_jan_26()

    % read data
    fn = 'RawData.mat';
    % fn = 'RawDataArray.mat'
    raw_data = load(fn);
    raw_data = getfield(raw_data,'rawData'); %_array');

    % get complex data
    cplx_data = hilbert(raw_data);

%     for i = 1:12
%         subplot(3,4,i);plot(real(cplx_data(:,i)));xlim([0 1024]);ylim([1e4,5.5e4]);
%     end

    % get fft data
    fft_data = fft(cplx_data);
%     for i = 1:12
%         subplot(3,4,i);plot(abs(fft_data(:,i)));xlim([0 512]);ylim([0,2e7]);
%     end

    fft_data_ref = fft_data(:,3);
    plot(abs(fft_data_ref));xlim([100 200]);hold on;
    rectangle('position',[100 0 40 2e7]);

    % binary filter
    depth_roi = [20 60];
    win_func = zeros([size(fft_data_ref,1) 1]);
    win_func(depth_roi(1) : depth_roi(2), 1) = 1;
    fft_data_cal = fft_data_ref .* win_func;

    % calculate raw data from fft data
    raw_data_cal = ifft(fft_data_cal);
    plot(real(raw_data_cal));xlim([0 1024]);ylim([-2e4,2e4]);% s7 right figure

    % extract phase info
    phase = unwrap(angle(raw_data_cal));
    phase_norm = phase - min(phase(:));
    phase_norm = phase_norm ./ max(phase_norm(:));

    phase_linear = linspace(0,1,length(raw_data_cal))';
    plot(phase_linear,color='blue');hold on;plot(phase_norm,color='red');xlim([0 1024]); % s8 left

    % non-linear sampling
    idx_non_linr_k = (1:length(raw_data_cal))'; % linear sampling gives non-linear phase
    poly_fit = fit(phase_norm, idx_non_linr_k, 'cubicinterp');
    idx_linr_k = poly_fit(phase_linear); % non-linear sampling gives linear phase

    plot(idx_linr_k,color='blue');hold on;plot(idx_non_linr_k,color='red');xlim([0 1024]);ylim([0 1024]); % s8 right

    % resampling to make sure phase of sampled data is linear
    [n, d] = size(cplx_data);
    for i = 1:size(cplx_data,2)
        cplx_data_rescaled(:,i) = ...
            interp1( [1:n]', real(cplx_data(:,i)), idx_linr_k, 'spline') ...
            + 1j .* interp1( [1:n]', imag(cplx_data(:,i)), idx_linr_k, 'spline' );
    end

    % *** Resampled data ***
    % fft of linear data
    fft_data_rescaled = fft(cplx_data_rescaled);
    fft_data_rescaled_ref = fft_data_rescaled(:,1);
    plot(abs(fft_data_rescaled_ref));

    % binary filter
    depth_roi = [20 60];
    win_func = zeros([size(fft_data_rescaled_ref,1) 1]);
    win_func(depth_roi(1) : depth_roi(2), 1) = 1;
    fft_data_rescaled_cal = fft_data_rescaled_ref .* win_func;

    % calculate raw data from fft data
    raw_data_rescaled_cal = ifft(fft_data_rescaled_cal);

    % before and after resampling- fringe
    % (s9) left figure
    plot(real(raw_data_rescaled_cal),color='blue');hold on;plot(real(raw_data_cal),color='red');xlim([0 1024]);ylim([-2e4,2e4]);

    % before and after resampling - amplitude
    % (s9) middle figure
    plot(abs(raw_data_rescaled_cal),color='blue');hold on;plot(abs(raw_data_cal),color='red');xlim([0 1024]);

    % extract phase info
    phase_rescaled = unwrap(angle(raw_data_rescaled_cal));
    phase_rescaled_norm = phase_rescaled - min(phase_rescaled(:));
    phase_rescaled_norm = phase_rescaled_norm ./ max(phase_rescaled_norm(:));

    % before and after resampling - phase
    % (s9) right figure
    plot(phase_rescaled_norm,color='blue');hold on;plot(phase_norm,color='red');xlim([0 1024]); % s8 left

    phase_rescaled = unwrap(angle(cplx_data_rescaled));
    phase_rescaled_norm = phase_rescaled - min(phase_rescaled(:));
    plot(phase_rescaled_norm,color='blue');hold on;plot(phase_norm,color='red');xlim([0 1024]); % s9 right


    % *** plot all data ***
    for k = 1:12
        depth_roi = [20 + (k-1)*40 20 + k*40];
        win_func = zeros([size(fft_data(:,k),1) 1]);
        win_func(depth_roi(1) : depth_roi(2), 1) = 1;
        fft_data_cal = fft_data(:,k) .* win_func;
    
        % calculate raw data from fft data
        raw_data_cal = ifft(fft_data_cal);
    
        % extract phase info
        phase = unwrap(angle(raw_data_cal));
        phase_norm = phase - min(phase(:));
        phase_norm = phase_norm ./ max(phase_norm(:));
    
        phase_linear = linspace(0,1,length(raw_data_cal))';
       
        % non-linear sampling
        idx_non_linr_k = (1:length(raw_data_cal))'; % linear sampling gives non-linear phase
        poly_fit = fit(phase_norm, idx_non_linr_k, 'cubicinterp');
        idx_linr_k = poly_fit(phase_linear); % non-linear sampling gives linear phase

        % resampling to make sure phase of sampled data is linear
        n = size(cplx_data,1);
        cplx_data_rescaled(:,k) = ...
            interp1( [1:n]', real(cplx_data(:,k)), idx_linr_k, 'spline') ...
            + 1j .* interp1( [1:n]', imag(cplx_data(:,k)), idx_linr_k, 'spline' );
    
        % fft of linear data
        fft_data_rescaled(:,k) = fft(cplx_data_rescaled(:,k));
    end

    % original fft signals, s10 left figure
    for i = 1:12
        subplot(1,2,1);plot(abs(fft_data(:,i)),color='red');hold on;xlim([0 512]);ylim([0 2e7]);
    end

    % resampled fft signal, s10 right figure
    for i = 1:12
        subplot(1,2,2);plot(abs(fft_data_rescaled(:,i)),color='blue');hold on;xlim([0 512]);ylim([0 2e7]);
    end

end
