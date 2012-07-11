%tests the factorization generated by matrix-factor
mat_name = 'test_matrices\stokes64';
base = 'C:\Users\Paul\Documents\My Dropbox\Online Resources\My Homework\UBC 2011W\NSERC\matrix-factor-recode\';
file = strcat(base, mat_name, '.mtx');
A = mmread(file);
l = mmread(strcat(base, 'output_matrices\outL.mtx'));
d = mmread(strcat(base, 'output_matrices\outD.mtx'));
p = mmread(strcat(base, 'output_matrices\outPerm.mtx'));
p = diag(p)+1;

avg_col_norm = norm(A,1)/size(A,1);
%the residual only provides a heuristic. e.g. lfil = 5, droptol = 0.1 for
%testmat 5 gives 4.5*10^-2 residual but lfil = 5, droptol = 0.01 gives
%1.8*10^-1 as a residual. the latter is a better preconditioner, but
%measures further from the original matrix.
fprintf('The residual is %d.\n', norm(A(p,p) - l*d*l', 1)/norm(A, 1));spy(A(p,p)); figure; spy(abs(l*d*l') > 0.01*avg_col_norm);

e = ones(size(A,1),1);
gmres(A(p,p),e,5,1e-6,50,l*d+0.01*avg_col_norm*speye(size(A)), l');