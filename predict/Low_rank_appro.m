function   [dim, wei, U_arr]  =  Low_rank_appro(nim, par, blk_arr, U_arr, it, flag)
b            =   par.win;
[h  w ch]    =   size(nim);
N            =   h-b+1;
M            =   w-b+1;
r            =   [1:N];
c            =   [1:M];
X       =     Im2Patch( nim, par );%%%%%X是图像块 每个图像块拉成列  (b^2, M*N)
Ys           =   zeros( size(X) );  %%%%%% (b^2, M*N)
W            =   zeros( size(X) ); %%%%%% (b^2, M*N)
L            =   size(blk_arr,2);%%%%%%
T            =   8;
for  i  =  1 : L
	B          =   X(:, blk_arr(:, i));  %%%%%  %%%%输入的B是相似的 图像块 的列
	if it==1 || mod(it, T)==0
		[tmp_y, tmp_w, U_arr(:,i)]   =    Weighted_SVT( double(B), par.c1, par.nSig^2, flag, par.c0 ); %mB,
	else
		[tmp_y, tmp_w]    =   Weighted_SVT_fast( double(B), par.c1, par.nSig^2, U_arr(:,i), flag, par.c0 ); %mB,
	end
	
	Ys(:, blk_arr(:,i))   =    Ys(:, blk_arr(:,i)) + tmp_y;%%%%辅助变量
	W(:, blk_arr(:,i))    =    W(:, blk_arr(:,i)) + tmp_w;%%%%辅助变量
end

dim     =  zeros(h,w);
wei     =  zeros(h,w);
k       =   0;
for i  = 1:b
	for j  = 1:b
		k    =  k+1;
		dim(r-1+i,c-1+j)  =  dim(r-1+i,c-1+j) + reshape( Ys(k,:)', [N M]);  %%%%%
		wei(r-1+i,c-1+j)  =  wei(r-1+i,c-1+j) + reshape( W(k,:)', [N M]);   %%%%%%
	end
end
return;


function  [X W U]   =   Weighted_SVT( Y, c1, nsig2,  flag, c0 )   %m, Y补丁
c1                =   c1*sqrt(2);
[U0,Sigma0,V0]    =   svd(full(Y),'econ');%%%%%% 图像块 U=9*9  再拉成列
Sigma0            =   diag(Sigma0); %%%%%对角化
if flag==1
	S                 =   max( Sigma0.^2/size(Y, 2) - nsig2, 0 );
	
	% 当目标函数为log时，
	thr               =   c1*nsig2./ ( sqrt(S) + eps );
	% 当目标函数为lp时，
	% ppppp = 0.05; thr = c1*nsig2.*ppppp.*( sqrt(S) + eps ).^(ppppp-1);
	
	% 当目标函数为geman时，
	%geman=0.01; thr  =   c1*nsig2*geman./ (( sqrt(S) + geman ).^2);
	
	% 当目标函数为lap时，
	%lap=0.01; thr  =   c1*nsig2/lap.*exp(-sqrt(S)./lap);
	
	% 当目标函数为核范数时，
	%thr=1;
	
	S                 =   soft(Sigma0, thr);
else
	S                 =   soft(Sigma0, c0*nsig2);
end
r                 =   sum( S>0 );%%%%求秩的大小
U                 =   U0(:,1:r);   %%%%
V                 =   V0(:,1:r);
X                 =   U*diag(S(1:r))*V'; %%%%%相似列低秩
if r==size(Y,1)
	wei           =   1/size(Y,1);
else
	wei           =   (size(Y,1)-r)/size(Y,1);
end
W                 =   wei*ones( size(X) );
X                 =   (X)*wei;% + m
U                 =   U0(:);%%%%U  方阵  U=9*9  再拉成列
return;


function  [X W]   =   Weighted_SVT_fast( Y, c1, nsig2,  U0, flag, c0 )%m,
c1                =   c1*sqrt(2);
n                 =   sqrt(length(U0));
U0                =   reshape(U0, n, n);
A                 =   U0'*Y;
Sigma0            =   zeros(size(A,1),1);
for i = 1:size(A,1)
	aa = A(i,:);
	Sigma0(i,1) = sqrt(aa*aa');
end
V0                =   (diag(1./Sigma0)*A)';

if flag==1
	S                 =   max( Sigma0.^2/size(Y, 2) - nsig2, 0 );
	
	% 当目标函数为log时，
	thr               =   c1*nsig2./ ( sqrt(S) + eps );
	
	% 当目标函数为lp时，
	% ppppp = 0.05; thr = c1*nsig2.*ppppp.*( sqrt(S) + eps ).^(ppppp-1);
	
	% 当目标函数为geman时，
	%geman=0.01; thr  =   c1*nsig2*geman./ (( sqrt(S) + geman ).^2);
	
	% 当目标函数为lap时，
	%lap=0.01; thr  =   c1*nsig2/lap.*exp(-sqrt(S)./lap);
	
	% 当目标函数为核范数时，
	%thr=1;
	
	S                 =   soft(Sigma0, thr);
else
	S                 =   soft(Sigma0, c0*nsig2);
end
r                 =   sum( S>0 ); %%%%%求和 秩的大小
P                 =   find(S);%%%%%%找出
X                 =   U0(:,P)*diag(S(P))*V0(:,P)';%%式子（17）
if r==size(Y,1)
	wei           =   1/size(Y,1);
else
	wei           =   (size(Y,1)-r)/size(Y,1);
end
W                 =   wei*ones( size(X) );
X                 =   (X )*wei;%+ m
return;


