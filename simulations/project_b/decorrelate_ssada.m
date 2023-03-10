function Dec_ssada = decorrelate_ssada(repBScans, num_frames, num_splits, cplx_OCT_mcorr_local_split)
% cplx_OCT_mcorr_local_split [.,.,3000]
% 1-500 split spectrum 1 for first frame, 501-1000 split spectrum 1 for
% second frame
% 1001-2000 split spectrum 2 for first frame, ...

    [a, b, ~] = size(cplx_OCT_mcorr_local_split);
    Dec_ssada = zeros([a b num_frames]);
    for i=1:num_frames
       decorr_sum = zeros([a b]);
       for j=1:repBScans-1
           for k=1:num_splits
               id1 = (k-1)*repBScans*num_frames+1 + (j-1)*num_frames + i-1;
               id2 = id1 + num_frames;
               %disp(id1);
               %disp(id2);
               Bscan_1 = cplx_OCT_mcorr_local_split(:,:,id1);
               Bscan_2 = cplx_OCT_mcorr_local_split(:,:,id2);
               cur_decorr = ((abs(Bscan_1).*abs(Bscan_2))...
                          ./((abs(Bscan_1).^2 + abs(Bscan_2).^2)./2));
               decorr_sum = decorr_sum + cur_decorr;
           end
       end
       Dec_ssada(:,:,i) = 1 - decorr_sum / ((repBScans-1)*num_splits);
    end
end
