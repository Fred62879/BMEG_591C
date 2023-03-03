
function rawData_rescaled = reSampling_CalSig(fftData, cal_RawData)
    ifftData = ifft(fftData); % equivalent to apply only hilbert without fft
    rawData_real = real(ifftData);
    rawData_imag = imag(ifftData);
    rawData_rescaled = zeros(size(fftData));

    idx = (1:length(cal_RawData))';

    for i = 1:size(fftData,2)
        phase = unwrap(angle(cal_RawData(:,i)));
        phase_norm = phase - min(phase(:));
        phase_norm = phase_norm ./ max(phase_norm(:));

        p = fit(phase_norm, idx, 'cubicinterp');

        phase_lin = linspace(0, 1, length(cal_RawData))';

        idx_linK = p(phase_lin);

        rawData_rescaled(:,i) = (interp1(...
            [1:size(rawData_real,1)]', rawData_real(:,i), idx_linK,'spline'))...
            + 1j.*(interp1([1:size(rawData_imag,1)]',rawData_imag(:,i),idx_linK,'spline'));
    end
end
