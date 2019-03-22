function [patches]=im2col_forSpatiotemporalFusion(X, winsize, step, resRate)

% step must be no larger than window size to avoid missing blocks
if step(1)>winsize(1)
	step(1) = winsize(1);
end
if step(2)>winsize(2)
	step(2) = winsize(2);
end

if nargin<4
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

% centers_y = [halfBlockHeight+1 centers_y height-halfBlockHeight];
points_y = [];	% centers of all blocks
for i=1:length(centers_y)-1
% 	c = (centers_y(i)+centers_y(i+1))/2;
% 	fr = c-step(1)/2:-step(1):centers_y(i)+2;
% 	r = floor([fr(end:-1:1) c+step(1)/2:step(1):centers_y(i+1)-2]);
	r = centers_y(i)+step(1):step(1):centers_y(i+1)-1;
	if (~isempty(r))
		offset = step(1)-1+centers_y(i+1)-r(end)-1;
		offset = ceil(offset/2);
		r = r - step(1) + 1 + offset;
	end
	points_y = [points_y centers_y(i) r];
end
fr = centers_y(1):-step(1):halfBlockHeight+1;
points_y = [fr(end:-1:2) points_y centers_y(end):step(1):height-halfBlockHeight];	% add first and last
if points_y(1)>halfBlockHeight+1
	points_y = [halfBlockHeight+1 points_y];	% make sure that all pixels are divided into blocks
end
if points_y(end)<height-halfBlockHeight
	points_y = [points_y height-halfBlockHeight];
end


% centers_x = [halfBlockWidth+1 centers_x width-halfBlockWidth];	% centers of exactly mapped blocks
points_x = [];	% centers of all blocks
for i=1:length(centers_x)-1
% 	c = (centers_x(i)+centers_x(i+1))/2;
% 	fr = c-step(2)/2:-step(2):centers_x(i)+2;
% 	r = floor([fr(end:-1:1) c+step(2)/2:step(2):centers_x(i+1)-2]);
	r = centers_x(i)+step(2):step(2):centers_x(i+1)-1;
	if (~isempty(r))
		offset = step(2)-1+centers_x(i+1)-r(end)-1;
		offset = ceil(offset/2);
		r = r - step(2) + 1 + offset;
	end
	points_x = [points_x centers_x(i) r];
end
fr = centers_x(1):-step(2):halfBlockWidth+1;
points_x = [fr(end:-1:2) points_x centers_x(end):step(2):width-halfBlockWidth];	% add first and last
if points_x(1)>halfBlockWidth+1
	points_x = [halfBlockWidth+1 points_x];	% make sure that all pixels are divided into blocks
end
if points_x(end)<width-halfBlockWidth
	points_x = [points_x width-halfBlockWidth];
end

patches=zeros(1+2+blockHeight*blockWidth, length(points_y)*length(points_x));
for iy=1:length(points_y)
	for ix=1:length(points_x)
		block = X(points_y(iy)-halfBlockHeight:points_y(iy)+halfBlockHeight, points_x(ix)-halfBlockWidth:points_x(ix)+halfBlockWidth);
		patches(:,ix+(iy-1)*length(points_x)) = [0; points_y(iy); points_x(ix); block(:)];
	end
end

%mark the central points of mapped blocks
pos = 1;
for iy=1:length(centers_y)
	for ix=1:length(centers_x)
		while (pos<length(points_y)*length(points_x))
			if (patches(2:3,pos) == [centers_y(iy); centers_x(ix)])
				patches(1,pos) = 1;
				break;
			end
			pos = pos + 1;
		end
	end
end

