function [avgOCT, Var, Sub, Dec] = process_oct_a(BM, cplxOCT_mcorr_local)
% 1536x500x1000->1536x500x500

    numFrames = size(cplxOCT_mcorr_local,3);
    for I = 1:BM:numFrames
        K = ((I-1)/BM)+1;

        Xconj = cplxOCT_mcorr_local(:,:,I+1).*conj(cplxOCT_mcorr_local(:,:,I));
        BulkOff = repmat(angle(sum(Xconj,1)), [size(Xconj,1) 1]);

        Bscan_1 = cplxOCT_mcorr_local(:,:,I);
        Bscan_2 = cplxOCT_mcorr_local(:,:,I+1) .* exp(-1j*BulkOff);

        % average oct
        avgOCT(:,:,K) = (Bscan_1+Bscan_2)./2;

        % variance
        Var(:,:,K) = abs(var(cat(3,Bscan_1,Bscan_2),0,3));

        % subtraction
        Sub(:,:,K) = abs(Bscan_1 - Bscan_2);

        % decorrelation
        Dec(:,:,K) = 1 - ((abs(Bscan_1).*abs(Bscan_2))...
                          ./((abs(Bscan_1).^2 + abs(Bscan_2).^2)./2));
    end
end
