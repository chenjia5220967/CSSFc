function param = SetSCDLParams()
param.win				=	5;
param.step			   =	floor(param.win/2);
param.rho = 5e-2;
param.lambda1		 =	   0.01;
param.lambda2		 =	   0.1;
param.mu			  =	   0.01;
param.sqrtmu		  =	   sqrt(param.mu);
param.nu			  =	   0.1;
param.nIter		   =	  10;
param.epsilon		 =	   5e-3;
param.t0			  =	   5;
param.K			   =	   256;%512;
param.L			   =	   param.win * param.win;
param.psf = fspecial('gaussian', param.win+2, 2.2);
param.nClass = 40;
param.nKmeansIters = 1000;

lassoParam.K = param.K;
lassoParam.lambda = param.lambda1;
lassoParam.iter=100;
lassoParam.L = param.win * param.win;
lassoParam.verbose = false;
lassoParam.numThreads = -1; % mt

param.lassoParam = lassoParam;
