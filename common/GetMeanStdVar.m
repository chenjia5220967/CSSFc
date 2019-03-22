function [m,stdVar] = GetMeanStdVar(A)
depths=size(A,3);
A=double(A);
for color=1:depths
	S=A(:,:,color);
	m(color)=mean(S(:));
	stdVar(color) = sqrt(sum(sum((S-m(color)).^2))/numel(S));
% 	fprintf('      band %d: mean:%f std_var:%f\n', color, m(color), stdVar(color));
end