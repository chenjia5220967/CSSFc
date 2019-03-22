function patches=im2col_forSpatiotemporalFusion_onlyCenterPoints(X, winsize, resRate)

if nargin<3
	resRate = 500/30;	% default for modis/landsat
end

[height,width] = size(X);
% winsize=[17 17];	% an odd number
blockHeight = winsize(1); halfBlockHeight = floor(blockHeight/2); blockHeight = halfBlockHeight*2+1;
blockWidth = winsize(2); halfBlockWidth = floor(blockWidth/2); blockWidth = halfBlockWidth*2+1;
if any([height width] < [blockHeight blockWidth]) % if neighborhood is larger than image
	patches = zeros(blockHeight*blockWidth,0);
	return
end

% centers of exactly mapped blocks
centers_x = round(resRate/2+1:resRate:width);	
centers_y = round(resRate/2+1:resRate:height);
centers_x(centers_x-halfBlockWidth<1 | centers_x+halfBlockWidth>width)=[];
centers_y(centers_y-halfBlockHeight<1 | centers_y+halfBlockHeight>height)=[];

patches=zeros(blockHeight*blockWidth, length(centers_y)*length(centers_x));
for iy=1:length(centers_y)
	for ix=1:length(centers_x)
		block = X(centers_y(iy)-halfBlockHeight:centers_y(iy)+halfBlockHeight, centers_x(ix)-halfBlockWidth:centers_x(ix)+halfBlockWidth);
		patches(:,ix+(iy-1)*length(centers_x)) = block(:);
	end
end
