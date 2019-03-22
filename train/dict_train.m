
% addpath('../../spams-matlab/build');
addpath('../common');

% Parameters Setting
param = SetSCDLParams();
nClass = 32;
param.win = 5;
param.K = 128;

Y1 = double(imread('../data/MOD09GHK.05-24-01.r-g-nir.tif'));
X1 = double(imread('../data/L7SR.05-24-01.r-g-nir.tif'));
Y3 = double(imread('../data/MOD09GHK.08-12-01.r-g-nir.tif'));
X3 = double(imread('../data/L7SR.08-12-01.r-g-nir.tif'));

band=1;
bandName={'red' 'green' 'nir'};
Y1 = Y1(:,:,band);
X1 = X1(:,:,band);
Y3 = Y3(:,:,band);
X3 = X3(:,:,band);
X31 = X3-X1;
Y31 = Y3-Y1;

XH = im2col(conv2(X31, param.psf, 'same') - X31, [param.win param.win], 'sliding');

[cls_idx0,dictVec] = kmeans_1(XH, nClass);	
dictVec = dictVec';

filepath=['../model/Params_' bandName{band} '_' num2str(param.win) 'x' num2str(param.K) 'x' num2str(nClass)];
save(filepath, 'dictVec','param');

XH = im2col(X31, [param.win param.win], 'sliding');
XL = im2col(Y31, [param.win param.win], 'sliding');

% SCDL
Dict=[];
for iClass = 1 : nClass	
	XH_t = XH(:,cls_idx0==iClass);
	XL_t = XL(:,cls_idx0==iClass);
	XH_t = XH_t - repmat(mean(XH_t), [param.win^2 1]);
	XL_t = XL_t - repmat(mean(XL_t), [param.win^2 1]);
	fprintf('dictionary learning: Cluster: %d (%d)\n', iClass, size(XH_t,2));
	D = mexTrainDL([XH_t;XL_t], param.lassoParam);
	Dh = D(1:param.win^2,:);
	Dl = D(param.win^2+1:end,:);
	W = eye(size(Dl, 2));
	alphaH = mexLasso([XH_t;XL_t], D, param.lassoParam);
	alphaL = alphaH;
	fprintf('Semi-Coupled dictionary learning: Cluster: %d (%d)\n', iClass, size(XH_t,2));
	[Dh, Dl, W] = coupled_DL(alphaH, XH_t, XL_t, Dh, Dl, W, param);
	Dict.DH{iClass} = Dh;
	Dict.DL{iClass} = Dl;
	Dict.W{iClass} = W;
end

filepath=['../model/DictSet_' bandName{band} '_' num2str(param.win) 'x' num2str(param.K) 'x' num2str(nClass)];
save(filepath, 'Dict');

rmpath('../common');