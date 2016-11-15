% function Mat = CCf(Mat,dr,NthOrd)
% %
% %    Mat = df(Mat,dr,NthOrd)
% %
% % Calculates Nth order (NthOrd) cross correlation like operation: 
% %   S(n)=S(n)-[S(n+dr)+S(n-dr)]/2
% % Applied N times to vector/matrix of size [r,c] by subtracting
% % data located dr points on opposite sides of each point.
% % If I/P is a matrix, the operation is performed col.-wise.
% % Notes: 
% % Initial and end dr number of data points will be corrupted
% %
% % also see df.m
% % [S.B.; last modified 200901]
% 
% if nargin == 1; dr=1; NthOrd=1; end;
% if nargin == 2; NthOrd=1; end;
% [r,c] = size(Mat); 
% tr=0; if isequal(r,1); tr=1; Mat=transpose(Mat); end; [r,c] = size(Mat);
% 
% %enlarge Mat
% % Mat = [Mat(end-4:end, :); Mat; Mat(1:4,:)];
% 
% for n=1:NthOrd
%    Mat = [Mat; ones(dr,1)*Mat(r,:)]-[ones(dr,1)*Mat(1,:); Mat]/2-[Mat(dr+1:r,:); ones(2*dr,1)*Mat(r,:)]/2;
%    Mat = Mat(1:r,:);
%    Mat(1:dr,:) = ones(dr,1)*Mat(dr+1,:);
%    Mat(r-dr+1:r,:) = ones(dr,1)*Mat(r-dr,:);
% end;
% % Mat=Mat(5:end-5,:);
% if isequal(tr,1); Mat=transpose(Mat); end;

function Mat = CCf(Mat,dr,NthOrd)
%
%    Mat = df(Mat,dr,NthOrd)
%
% Calculates Nth order (NthOrd) cross correlation like operation: 
%   S(n)=S(n)-[S(n+dr)+S(n-dr)]/2
% Applied N times to vector/matrix of size [r,c] by subtracting
% data located dr points on opposite sides of each point.
% If I/P is a matrix, the operation is performed col.-wise.
% Notes: 
% Initial and end dr number of data points will be corrupted
%
% also see df.m
% [S.B.; last modified 200901]

if nargin == 1; dr=1; NthOrd=1; end;
if nargin == 2; NthOrd=1; end;
[r,c] = size(Mat); 
tr=0; if isequal(r,1); tr=1; Mat=transpose(Mat); end; [r,c] = size(Mat);


for n=1:NthOrd
    Mat = [Mat; Mat(1:dr,:)]-[Mat((r-dr+1):r,:); Mat]/2-[Mat(dr+1:r,:); Mat(1:2*dr,:)]/2;
    Mat = Mat(1:r,:);
end

if isequal(tr,1); Mat=transpose(Mat); end;