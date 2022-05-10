//-*- mode: c++ -*-
#ifndef _SOLVER_FGMRES_H_
#define _SOLVER_FGMRES_H_

#include <string>
#include <algorithm>
#include <cmath>

#define  epsmac  1.0e-16

template<class el_type>
el_type DNRM2(int n, vector<el_type>& x, int incx) {
	el_type res=0.;
	for (int i=0; i<n; i+=incx)
		res+=x[i]*x[i];
	return std::sqrt(res);
}

template<class el_type>
el_type DDOT(int n, vector<el_type>& x, int incx, vector<el_type>& y, int incy) {
	el_type res=0.;
	for (int i=0; i<n; i+=incx)
		res+=x[i]*y[i];
	return res;
}

template<class el_type>
void DAXPY(int n, el_type t, vector<el_type>& x, int incx, vector<el_type>& y, int incy) {
	for (int i=0; i<n; i+=incx) y[i] += t*x[i];
	return;
}

template<class el_type>
void DSCAL(int n, el_type t, vector<el_type>& x, int incx) {
	for (int i=0; i<n; i+=incx) x[i]*=t;
	return;
}

template<class el_type, class mat_type >
void solver<el_type, mat_type> :: fgmres(int max_iter, int kdim, double tol) {
	//Zero out solution vector
	int n = A.n_rows();
	sol_vec.resize(n, 0);

    int i, i1, ii, j, k, k1, its, kdim_p1, pti, pti1, ptih=0, retval, one = 1;
	double *hh, *c, *s, *rs, t;
	double negt, beta, eps1=0, gam;
	kdim_p1 = kdim+1;

	el_type norm_rhs = norm(rhs, 2.0);

    vector<el_type> tmp1(n);
    vector<el_type> tmp2(n);

    vector< vector<el_type> > vv(kdim_p1, vector<el_type>(n));
    vector< vector<el_type> > z(kdim, vector<el_type>(n));
	//vv = (double *)malloc(kdim_p1*n*sizeof(double));
	//z  = (double *)malloc(kdim*n*sizeof(double));
	kdim_p1 = kdim+1;
	hh = (double *)malloc((kdim_p1*(kdim+3))*sizeof(double));
	c  = hh+kdim_p1*kdim ; s  = c+kdim_p1;  rs = s+kdim_p1;
	/*-------------------- outer loop starts here */
	retval = 0;
	its = 0;
	/*-------------------- Outer loop */
	while (its < max_iter) {
		/*-------------------- compute initial residual vector */
		A.multiply(sol_vec, vv[0]);
		//Amat->matvec(Amat, sol, vv);
		for (j=0; j<n; j++)
			vv[0][j] = rhs[j] - vv[0][j];    /*  vv[0]= initial residual */

		beta = DNRM2(n, vv[0], one);
		if (beta == 0.0)
			break;
		t = 1.0 / beta;
		/*--------------------   normalize:  vv    =  vv   / beta */
		DSCAL(n, t, vv[0], one);
		if (its == 0)
			eps1 = tol*beta;
		/*--------------------initialize 1-st term  of rhs of hessenberg mtx */
		rs[0] = beta;
		i = 0;
		/*-------------------- Krylov loop*/
		i = -1;
		pti=pti1=0;
		while((i < kdim-1) && (beta > eps1) && (its++ < max_iter))  {
			i++;
			i1   = i+1;
			pti  = i*n;
			pti1 = i1*n;
			/*------------------------------------------------------------
			  |  (Right) Preconditioning Operation   z_{j} = M^{-1} v_{j}
			  +-----------------------------------------------------------*/
			L.backsolve(vv[i], tmp1);
            D.solve(tmp1, tmp2);
            L.forwardsolve(tmp2, z[i]);

			/*-------------------- matvec operation w = A z_{j} = A M^{-1} v_{j} */
			A.multiply(z[i], vv[i1]);
			//Amat->matvec(Amat, &z[pti], &vv[pti1]);
			/*-------------------- modified gram - schmidt...
			  |     h_{i,j} = (w,v_{i});
			  |     w  = w - h_{i,j} v_{i}
			+------------------------------------------------------------*/
			ptih=i*kdim_p1;
			for (j=0; j<=i; j++) {
				t = DDOT(n, vv[j], one, vv[i1], one);
				hh[ptih+j] = t;
				negt = -t;
				DAXPY(n, negt, vv[j], one, vv[i1], one);
			}
			/*-------------------- h_{j+1,j} = ||w||_{2}    */
			t = DNRM2(n, vv[i1], one);
			hh[ptih+i1] = t;
			if (t == 0.0)
				return;
			t = 1.0/t;
			/*-------------------- v_{j+1} = w / h_{j+1,j}  */
			DSCAL(n, t, vv[i1], one);
			/*-------- done with modified gram schimdt/arnoldi step
			  | now  update factorization of hh.
			  | perform previous transformations  on i-th column of h
			  +-------------------------------------------------------*/
			for (k=1; k<=i; k++) {
				k1 = k-1;
				t = hh[ptih+k1];
				hh[ptih+k1] = c[k1]*t + s[k1]*hh[ptih+k];
				hh[ptih+k] = -s[k1]*t + c[k1]*hh[ptih+k];
			}
			gam = sqrt(pow(hh[ptih+i],2) + pow(hh[ptih+i1],2) );
			/*-------------------- check if gamma is zero */
			if (gam == 0.0) gam = epsmac;
			/*-------------------- get  next plane rotation    */
			c[i] = hh[ptih+i]/gam;
			s[i] = hh[ptih+i1]/gam;
			rs[i1] = -s[i]*rs[i];
			rs[i] =  c[i]*rs[i];
			/*-------------------- get residual norm + test convergence*/
			hh[ptih+i] = c[i]*hh[ptih+i] + s[i]*hh[ptih+i1];
			beta = fabs(rs[i1]);
			/*-------------------- end [inner] while loop [Arnoldi] */
			printf("iter %d, relative residual %e.\n", its, beta/norm_rhs);
		}
		/*-------------------- now compute solution. 1st, solve upper
		  triangular system*/
		rs[i] = rs[i]/hh[ptih+i];
		for (ii=i-1; ii>=0; ii--) {
			t=rs[ii];
			for (j=ii+1; j<=i; j++)
				t -= hh[j*kdim_p1+ii]*rs[j];
			rs[ii] = t/hh[ii*kdim_p1+ii];
		}
		/*---------- linear combination of z_j's to get sol. */
		for (j=0; j<= i; j++)
			DAXPY(n, rs[j], z[j], one, sol_vec, one);
		/*--------------------  restart outer loop if needed */
		if (beta < eps1)
			break;
		else
			if (its >= max_iter)
				retval = 1;
		/*---------- end main [outer] while loop */
	}
	/*-------------------- prepare to return */
	//free(vv);
	//free(z);
	free(hh);

	A.multiply(sol_vec, tmp1);
	for (j=0; j<n; j++)
		tmp2[j] = rhs[j] - tmp1[j];
	el_type norm_num = DNRM2(n,tmp2,one);
	el_type norm_x = DNRM2(n,sol_vec,one);
	el_type norm_A = A.fronorm();
	el_type norm_b = norm_rhs;
	el_type norm_den = norm_b + norm_A*norm_x;
	el_type nrbe = norm_num/norm_den;

    std::string iter_str = "iterations";
    if (k-1 == 1) iter_str = "iteration";

	if (msg_lvl) printf("FGMRES took %i %s and got down to relative residual %e nrbe %e\n", its, iter_str.c_str(), beta/norm_rhs, nrbe);
	return;
}

#endif // _SOLVER_FGMRES_H_
