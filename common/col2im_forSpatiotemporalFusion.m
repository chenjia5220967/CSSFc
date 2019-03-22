function im = col2im_forSpatiotemporalFusion(patches, winsize)

blockHeight = winsize(1); halfBlockHeight = floor(blockHeight/2); blockHeight = halfBlockHeight*2+1;
blockWidth = winsize(2); halfBlockWidth = floor(blockWidth/2); blockWidth = halfBlockWidth*2+1;

height = patches(2, end) + halfBlockHeight;
width = patches(3, end) + halfBlockWidth;

im = zeros(height, width);
count = zeros(height, width);

for i=1:size(patches, 2)
	center_y = patches(2,i);
	center_x = patches(3,i);
	block = reshape(patches(4:end,i), winsize);
	region_y = [center_y-halfBlockHeight:center_y+halfBlockHeight];
	region_x = [center_x-halfBlockWidth:center_x+halfBlockWidth];
	im(region_y, region_x) = im(region_y, region_x) + block;
	count(region_y, region_x) = count(region_y, region_x) + 1;
end

if (any(count<1))
	fprintf('Warning: not all pixels are divided into patches. Check im2col function for reason.');
	count(count<1) = 1;
end
im = im ./ count;