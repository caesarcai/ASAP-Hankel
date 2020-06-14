function [U_out,V_out] = trim(U,Sig,V,threshold_U,threshold_V)

row_norm_square_U = sum(U.^2,2);
big_rows_U = row_norm_square_U > threshold_U;
U(big_rows_U,:) = bsxfun(@times,U(big_rows_U,:),sqrt(threshold_U./(row_norm_square_U(big_rows_U))));   %This is A in paper

row_norm_square_V = sum(V.^2,2);
big_rows_V = row_norm_square_V > threshold_V;
V(big_rows_V,:) = bsxfun(@times,V(big_rows_V,:),sqrt(threshold_V./(row_norm_square_V(big_rows_V))));   %This is B in paper

% The following 2 QR-decomp can be computed parallelly
[Q1,R1] = qr(U,0);
[Q2,R2] = qr(V,0);
[U_temp,~,V_temp] = svd(R1*Sig*R2');

% The following 2 matrices multiplications can be computed parallelly
U_out = Q1*U_temp;
V_out = Q2*V_temp;

end