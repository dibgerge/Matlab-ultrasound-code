function [ output_args ] = getNoiseLevel( sig, threshold )
% getNoiseLevel sig is a three dimensional array representing all A-Scans
% threshold < 1, ratio of considering a signal a noise.

N1 = size(sig, 2);
N2 = size(sig, 3);

amp = max(max(max(abs(sig))));
noise_level = 0;
tot = [];

mx = squeeze(max(abs(sig)));

ind = find(mx(:) < threshold*amp);
s = sig(:,:);
s(:,ind) = [];


% 
% for i=1:N1
%     for j=1:N2
%         mx = max(abs(sig(:,i,j)));
%         if mx < threshold*amp
%             noise_level = noise_level + mean(sig(:,i,j));
%             tot = [tot; sig(:,i,j)];
%         %    plot(sig(:,i,j)); hold on
%         end
%     end
% end


 hist(s(:),100);