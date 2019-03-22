addpath('../common');

band=1;	% adjustable
bandName={'red' 'green' 'nir'};

resRate=10/3;	% adjustable

% read data
Y01 = double(imread('../data/MOD09GHK.05-24-01.r-g-nir.tif'));
X01 = double(imread('../data/L7SR.05-24-01.r-g-nir.tif'));
Y02 = double(imread('../data/MOD09GHK.07-11-01.r-g-nir.tif'));
X02 = double(imread('../data/L7SR.07-11-01.r-g-nir.tif'));
Y03 = double(imread('../data/MOD09GHK.08-12-01.r-g-nir.tif'));
X03 = double(imread('../data/L7SR.08-12-01.r-g-nir.tif'));

Y01 = Y01(:,:,band);
Y02 = Y02(:,:,band);
Y03 = Y03(:,:,band);
X01 = X01(:,:,band);
X02 = X02(:,:,band);
X03 = X03(:,:,band);

[height, width] = size(X01);
regionHeight=400;
regionWidth=400;
nVertRegions=floor(height/regionHeight);
nHoriRegions=floor(width/regionWidth);
region(nVertRegions*nHoriRegions, 1)=struct('stRow', [], 'endRow', [], 'stCol', [], 'endCol', []);
for y=1:floor(height/regionHeight)
	for x=1:nHoriRegions
		iRegion = (y-1)*nHoriRegions+x;
		region(iRegion).stRow = (y-1)*regionHeight+1;
		region(iRegion).endRow = y*regionHeight;
		region(iRegion).stCol = (x-1)*regionWidth+1;
		region(iRegion).endCol = x*regionWidth;
	end
end

iRegion=1;	% adjustable
XB1=X01(region(iRegion).stRow:region(iRegion).endRow, region(iRegion).stCol:region(iRegion).endCol,:);
XB2=X02(region(iRegion).stRow:region(iRegion).endRow, region(iRegion).stCol:region(iRegion).endCol,:);
XB3=X03(region(iRegion).stRow:region(iRegion).endRow, region(iRegion).stCol:region(iRegion).endCol,:);
YB1=Y01(region(iRegion).stRow:region(iRegion).endRow, region(iRegion).stCol:region(iRegion).endCol,:);
YB2=Y02(region(iRegion).stRow:region(iRegion).endRow, region(iRegion).stCol:region(iRegion).endCol,:);
YB3=Y03(region(iRegion).stRow:region(iRegion).endRow, region(iRegion).stCol:region(iRegion).endCol,:);

%%%%%%%%%% begin of generating medium image %%%%%%%%%%%%%%%%%%
% filepath=['../model/' 'KMeans_5x5_32_Factor3_2'];
param.win=5;
param.K=128;
nClass=32;
filepath=['../model/Params_' bandName{band} '_' num2str(param.win) 'x' num2str(param.K) 'x' num2str(nClass)];
fprintf('param file: %s\n', filepath);
load(filepath, 'dictVec', 'param');

param.step=2;
param.dictVec = dictVec;

% filepath=['../model/' 'Dict_SR_Factor3.mat'];
filepath=['../model/DictSet_' bandName{band} '_' num2str(param.win) 'x' num2str(param.K) 'x' num2str(nClass)];
fprintf('dictionary file: %s\n', filepath);
load(filepath, 'Dict');
param.Dict=Dict;

param.lambda_3 = 1;
param.lambda_4 = 10000;

param.resRate = 500/30/resRate;

% downsample to middle resolution
Y1=imresize(YB1, 1.0/resRate);
Y2=imresize(YB2, 1.0/resRate);
Y3=imresize(YB3, 1.0/resRate);
X1=imresize(XB1, 1.0/resRate);
X2=imresize(XB2, 1.0/resRate);	% for psnr evaluation only
X3=imresize(XB3, 1.0/resRate);

% Measurement matrix
B=im2col_forSpatiotemporalFusion(imresize(X01-X03, 1.0/resRate), [param.win param.win], [param.step param.step], resRate);
B=B(4:end,:);
B2=im2col_forSpatiotemporalFusion(imresize(Y01-Y03, 1.0/resRate), [param.win param.win], [param.step param.step], resRate);
B2=B2(end-floor(param.win*param.win/2),:);
[U]=pca(B');
[~, pos] = min(sum((B'*U-repmat(B2',1,size(U,2))).^2, 1));
param.M0 = reshape(U(:, pos), 1, []);

param.originalImage=X2-X1;	% only for psnr comparison
XB21m = ebscdl_interp(Y2-Y1, param);%, refRec-X1);

param.originalImage=X2-X3;	% only for psnr comparison
XB23m = ebscdl_interp(Y2-Y3, param);%, refRec-X3);

X2m = ((X1+XB21m)+(X3+XB23m))/2;

[psnr, mse, maxerr]=psnr_mse_maxerr(X2m(:), X2(:));
fprintf('\nResult of middle image after regularization: PSNR:%f rmse:%f\n', psnr, sqrt(mse));

%%%%%%%%%% end of generating medium image %%%%%%%%%%%%%%%%%%
%%
X2m = imresize(X2m, size(XB1));

param.win=5;
param.K=128;
nClass=32;
% nBands=length(bands);
filepath=['../model/Params_' bandName{band} '_' num2str(param.win) 'x' num2str(param.K) 'x' num2str(nClass)];
fprintf('param file: %s\n', filepath);
load(filepath, 'dictVec', 'param');
param.dictVec = dictVec;

param.step = param.win-1;
param.lambda1=0.02;
param.lambda2=0;
param.lambda_3 = 1;
param.lambda_4 = 20;

param.resRate = 500/30;%resRate;

filepath=['../model/DictSet_' bandName{band} '_' num2str(param.win) 'x' num2str(param.K) 'x' num2str(nClass)];
fprintf('dictionary file: %s\n', filepath);
load(filepath, 'Dict');
param.Dict=Dict;

% Measurement matrix
B=im2col_forSpatiotemporalFusion(X01-X03, [param.win param.win], [param.step param.step], param.resRate);
B=B(4:end,:);
B2=im2col_forSpatiotemporalFusion(Y01-Y03, [param.win param.win], [param.step param.step], param.resRate);
B2=B2(end-floor(param.win*param.win/2),:);
[U]=pca(B');
[mindif, pos] = min(sum((B'*U-repmat(B2',1,size(U,2))).^2, 1));
param.M0 = reshape(U(:, pos), 1, []);	

param.originalImage=XB2-XB1;	% only for psnr comparison
XB21 = ebscdl_interp(YB2-YB1, param, X2m-XB1);

param.originalImage=XB2-XB3;	% only for psnr comparison
X2B3 = ebscdl_interp(YB2-YB3, param, X2m-XB3);
%%
XB2p = ((XB1+XB21)+(XB3+X2B3))/2;

[psnr, mse, maxerr]=psnr_mse_maxerr(XB2p(:), XB2(:));
fprintf('\nResult of CSSCDL:  PSNR:%f rmse:%f\n', psnr, sqrt(mse));

imwrite(uint16(XB2p),['../result/x2p-CSSF-region' num2str(iRegion) '-rate' num2str(resRate) '-' bandName{band}  '-' num2str(psnr) '.tif']);

rmpath('../common');