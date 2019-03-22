function [DH, DL, W] = coupled_DL(alphaH, XH, XL, DH, DL, W, par)
% Semi-Coupled Dictionary Learning
% Shenlong Wang
% Reference: S Wang, L Zhang, Y Liang and Q. Pan, "Semi-coupled Dictionary Learning with Applications in Super-resolution and Photo-sketch Synthesis", CVPR 2012

[dimX, numX]        =       size(XH);
dimY                =       size(alphaH, 1);
numD                =       size(DH, 2);
rho                 =       par.rho;
lambda1             =       par.lambda1;
lambda2             =       par.lambda2;
mu                  =       par.mu;
sqrtmu              =       sqrt(mu);
nu                  =       par.nu;
nIter               =       par.nIter;
t0                  =       par.t0;
epsilon             =       par.epsilon;
dictSize = par.K;

%lasso param
param.lambda       = lambda1; % not more than 20 non-zeros coefficients
param.lambda2       = lambda2;
param.mode          = 2;       % penalized formulation
param.approx=0;
param.K = par.K;
param.L = par.L;

f = 0;	% cost
for t = 1 : nIter
	
	% Alphat = mexLasso(Xt,D,param);
	f_prev = f;
	
	alphaL = mexLasso([XL;sqrtmu * full(alphaH)], [DL; sqrtmu * W], param);
	alphaH = mexLasso([XH;sqrtmu * W * full(alphaL)], [DH; sqrtmu * eye(size(alphaH, 1))], param);
	
	% Update D with K-SVD
	for i=1:dictSize
		ai        =    alphaL(i,:);
		Y         =    XL-DL*alphaL+DL(:,i)*ai;
		di        =    Y*ai';
		di        =    di./(norm(di,2) + eps);
		DL(:,i)    =    di;
	end
	for i=1:dictSize
		ai        =    alphaH(i,:);
		Y         =    XH-DH*alphaH+DH(:,i)*ai;
		di        =    Y*ai';
		di        =    di./(norm(di,2) + eps);
		DH(:,i)    =    di;
	end
	
	% Update W
	%     Ws = Alphap * Alphas' * inv(Alphas * Alphas' + par.nu * eye(size(Alphas, 1))) ;
	%     Wp = Alphas * Alphap' * inv(Alphap * Alphap' + par.nu * eye(size(Alphap, 1))) ;
	W = (1 - rho) * W  + rho * alphaH * alphaL' * inv(alphaL * alphaL' + par.nu * eye(size(alphaL, 1))) ;
% 	Wp = (1 - rho) * Wp  + rho * alphaL * alphaH' * inv(alphaH * alphaH' + par.nu * eye(size(alphaH, 1))) ;

	% Alpha = pinv(D' * D + lambda2 * eye(numD)) * D' * X;
	P1 = XH - DH * alphaH;
	P1 = P1(:)'*P1(:) / 2;
	P2 = lambda1 *  norm(alphaH, 1);
% 	P3 = alphaL - Wp * alphaH;
% 	P3 = P3(:)'*P3(:) / 2;
% 	P4 = nu * norm(Wp, 'fro');
	fp = 1 / 2 * P1 + P2;% + mu * (P3 + P4);
	
	P1 = XL - DL * alphaL;
	P1 = P1(:)'*P1(:) / 2;
	P2 = lambda1 *  norm(alphaL, 1);
	P3 = alphaH - W * alphaL;
	P3 = P3(:)'*P3(:) / 2;
	P4 = nu * norm(W, 'fro');
	fs = 1 / 2 * P1 + P2 + mu * (P3 + P4);
	f = fp + fs;
	if (abs(f_prev - f) / f < epsilon)
		break;
	end
	%     fprintf('Energy: %d\n',f);
	%     save tempDict_SR_NL Ds Dp Ws Wp par param i;
	% fprintf('Iter: %d, E1 : %d, E2 : %d, E : %d\n', t, mu * (P1 + P2), (1 - mu) * (P3 + P4), E);
end
