function hdImage = ebscdl_interp(ldImage, param, initOutImage)
% ldImage: up-sampled low-resolution image

if (nargin<3)
    initOutImage = ldImage;
end
hdImage = initOutImage;

lassoParam=param.lassoParam;
Dict=param.Dict;

nOuterLoop = 5;
nInnerLoop = 1;

cls_num = size(param.dictVec,1);
[height, width, nBands] = size(ldImage);
nBands = 1;

M = param.M0;
param.M0=0;	%||M-M0||->||M||

% Used for calculate CSNR
param.height = height;
param.width = width;

if exist('param.psf', 'var')
	psf = param.psf;
else
	psf = fspecial('gaussian', param.win+2, 2.2);
end

ldPatches = im2col_forSpatiotemporalFusion(ldImage, [param.win param.win], [param.step param.step], param.resRate); 
nBlocks = size(ldPatches, 2);

nAtoms = size(Dict.DH{1},2);
param.K=nAtoms;
AL = zeros(nAtoms*nBands, nBlocks);
AH = zeros(nAtoms*nBands, nBlocks);

if (size(param.originalImage,1)>1)
	[psnr, mse, maxerr]=psnr_mse_maxerr(hdImage(:), param.originalImage(:));
	fprintf('initial psnr: %f; rmse: %f\n', psnr, sqrt(mse));
end

for outIter = 1 : nOuterLoop
for iter = 1 : nInnerLoop
	fprintf('Iter: %d x %d\n   ', outIter,iter);
	
	% clustering
	if iter == 1
		XH0 = hdImage;
% 		XH = data2patch(conv2(XH0, psf, 'same') - XH0, param);
		XH = im2col_forSpatiotemporalFusion(conv2(XH0, psf, 'same') - XH0, [param.win param.win], [param.step param.step], param.resRate); 
		cls_idx = setPatchIdx(XH(4:end,:), param.dictVec');
		clear XH XH0;
	end
	
% 	%normalize	
	[m, stdev] = GetMeanStdVar(hdImage(:));
	hdImage = (hdImage-m)/stdev;
	ldImage = (ldImage-m)/stdev;
	
	% image to blocks
	XH = im2col_forSpatiotemporalFusion(hdImage, [param.win param.win], [param.step param.step], param.resRate); % X: the high-resolution image
	XL = im2col_forSpatiotemporalFusion(ldImage, [param.win param.win], [param.step param.step], param.resRate); % Y: the low-resolution image
	
	meanX = repmat(mean(XH(4:end,:),1), [param.win^2 1]);
	XL(4:end,:) = XL(4:end,:) - meanX;
	XH(4:end,:) = XH(4:end,:) - meanX; 
	
	for iClass = 1 : cls_num
		idx_cluster = find(cls_idx == iClass);
		length_idx = length(idx_cluster);
		fprintf(' class %d(%d):', iClass, length_idx);
		if (length_idx==0)
			continue;
		end
		start_idx = [1:10000:length_idx, length_idx];
		for j  = 1 : length(start_idx) - 1
			idx_temp = idx_cluster(start_idx(j):start_idx(j+1));
			Xl	= double(XL(4:end, idx_temp));
			Xh	= double(XH(4:end, idx_temp));
			Dl	= Dict.DL{iClass};
			Dh	= Dict.DH{iClass};
			W	= Dict.W{iClass};
			
			if (iter == 1)
				alphaL = mexLasso(Xl, Dl, lassoParam);
				alphaH = W * alphaL;
				Xh = Dh * alphaH;
			else
				alphaH = AH(:, idx_temp);
			end
			
			alphaL = mexLasso([Xl;param.sqrtmu * full(alphaH)], [Dl; param.sqrtmu * W], lassoParam);

			alphaH = mexLasso([Xh;param.sqrtmu * W * full(alphaL)], [Dh; param.sqrtmu * eye(size(alphaH, 1))], lassoParam);
			
			AL(:, idx_temp) = alphaL;
			AH(:, idx_temp) = alphaH;
			flag_center = logical(XH(1, idx_temp)==1);
			
			if (any(flag_center))
				%centers:
				y = Xl(ceil(end/2),flag_center);
				B = Dh*alphaH(:,flag_center)+param.lambda_3*M'*y;
				A = eye(size(M,2))+param.lambda_3*(M'*M);
				Xh = A\B;
				XH(4:end, idx_temp(flag_center)) = Xh;
			end
			if (any(~flag_center))
				%non-centers:
				Xh = Dh * alphaH(:,~flag_center);
				XH(4:end, idx_temp(~flag_center)) = Xh;
			end
		end
	end
	fprintf('\n');
	
% % 	update Measurement matrix
% 	flag_center = logical(XH(1, :)==1);
% 	y = XL(end-floor(param.win*param.win/2), flag_center);
% 	X = XH(4:end, flag_center);
% 	M = (param.lambda_4*param.M0+y*X')/(X*X'+param.lambda_4*eye(size(X,1)));
	
	XH(4:end,:) = XH(4:end,:) + meanX;
	for iBand=1:nBands
		hdImage(:,:,iBand) = col2im_forSpatiotemporalFusion(XH, [param.win param.win]);
	end
	
	hdImage = hdImage*stdev+m;	
	
%%%%%%%%%%%%%%%%%% regularization %%%%%%%%%%%%%%%%%%%%%
% 	fprintf('  start of regularization\n');

	%nonlocal regulaziation
% 	[N, ~] = Compute_NLM_Matrix( hdImage, 5, param);
% 	NTN	   = N'*N*0.05;
% 	im_f = sparse(double(reshape(hdImage, [], 1)));
% 	for i = 1 : fix(60 / iter.^2)
% 		im_f = im_f  - NTN*im_f;
% 	end
% 	hdImage = reshape(full(im_f), height, width);	

% 	fprintf('  end of regularization\n');	
%%%%%%%%%%%%%%%%%% regularization %%%%%%%%%%%%%%%%%%%%%		
	
	if (size(param.originalImage,1)>1)
		[psnr, mse, maxerr]=psnr_mse_maxerr(hdImage(:), param.originalImage(:));
		fprintf('  mapping psnr: %f; rmse: %f\n', psnr, sqrt(mse));
	end
end
end

