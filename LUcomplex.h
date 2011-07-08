/* 
 * File:   LUcomplex.h
 * Author: veal
 *
 * Created on March 11, 2011, 1:29 PM
 */

#ifndef LUCOMPLEX_H
#define	LUCOMPLEX_H

typedef complex<Doub> Comp;
#ifdef	__cplusplus
extern "C" {
#endif

struct LUcomplex
{
	Int n;
	NRmatrix<Comp> lu;
	VecInt indx;
	Comp d;
	LUcomplex(const NRmatrix<Comp> &a);
	void solve(const NRvector<Comp> &b, NRvector<Comp> &x);
	void solve(const NRmatrix<Comp> &b, NRmatrix<Comp> &x);
	void inverse(NRmatrix<Comp> &ainv);
	Comp det();
	void mprove(const NRvector<Comp> &b, NRvector<Comp> &x);
	const NRmatrix<Comp> &aref;
};
LUcomplex::LUcomplex(const NRmatrix<Comp> &a) : n(a.nrows()), lu(a), aref(a), indx(n) {
	const Doub TINY=1.0e-40;
	Int i,imax,j,k;
	Doub big,temp;
        Comp myTemp;
	VecDoub vv(n);
	d=Comp(1.0, 0.0);
	for (i=0;i<n;i++) {
		big=0.0;
		for (j=0;j<n;j++)
			if ((temp=abs(lu[i][j])) > big) big=temp;
		if (big == 0.0) throw("Singular matrix in LUdcmp");
		vv[i]=1.0/big;
	}
	for (k=0;k<n;k++) {
		big=0.0;
		for (i=k;i<n;i++) {
			temp=vv[i]*abs(lu[i][k]);
			if (temp > big) {
				big=temp;
				imax=i;
			}
		}
		if (k != imax) {
			for (j=0;j<n;j++) {
				myTemp=lu[imax][j];
				lu[imax][j]=lu[k][j];
				lu[k][j]=myTemp;
			}
			d = -d;
			vv[imax]=vv[k];
		}
		indx[k]=imax;
		if (abs(lu[k][k]) == 0.0) lu[k][k]=Comp(TINY, 0.0);
		for (i=k+1;i<n;i++) {
			myTemp=lu[i][k] /= lu[k][k];
			for (j=k+1;j<n;j++)
				lu[i][j] -= myTemp*lu[k][j];
		}
	}
}
void LUcomplex::solve(const NRvector<Comp> &b, NRvector<Comp> &x)
{
	Int i,ii=0,ip,j;
	Comp sum;
	if (b.size() != n || x.size() != n)
		throw("LUdcmp::solve bad sizes");
	for (i=0;i<n;i++) x[i] = b[i];
	for (i=0;i<n;i++) {
		ip=indx[i];
		sum=x[ip];
		x[ip]=x[i];
		if (ii != 0)
			for (j=ii-1;j<i;j++) sum -= lu[i][j]*x[j];
		else if (abs(sum) != 0.0)
			ii=i+1;
		x[i]=sum;
	}
	for (i=n-1;i>=0;i--) {
		sum=x[i];
		for (j=i+1;j<n;j++) sum -= lu[i][j]*x[j];
		x[i]=sum/lu[i][i];
	}
}

void LUcomplex::solve(const NRmatrix<Comp> &b, NRmatrix<Comp> &x)
{
	int i,j,m=b.ncols();
	if (b.nrows() != n || x.nrows() != n || b.ncols() != x.ncols())
		throw("LUdcmp::solve bad sizes");
	NRvector<Comp> xx(n);
	for (j=0;j<m;j++) {
		for (i=0;i<n;i++) xx[i] = b[i][j];
		solve(xx,xx);
		for (i=0;i<n;i++) x[i][j] = xx[i];
	}
}
void LUcomplex::inverse(NRmatrix<Comp> &ainv)
{
	Int i,j;
	ainv.resize(n,n);
	for (i=0;i<n;i++) {
		for (j=0;j<n;j++) ainv[i][j] = 0.;
		ainv[i][i] = 1.;
	}
	solve(ainv,ainv);
}
Comp LUcomplex::det()
{
	Comp dd = d;
	for (Int i=0;i<n;i++) dd *= lu[i][i];
	return dd;
}
void LUcomplex::mprove(const NRvector<Comp> &b, NRvector<Comp> &x)
{
	Int i,j;
	NRvector<Comp> r(n);
	for (i=0;i<n;i++) {
		Comp sdp = -b[i];
		for (j=0;j<n;j++)
			sdp += aref[i][j] * x[j];
		r[i]=sdp;
	}
	solve(r,r);
	for (i=0;i<n;i++) x[i] -= r[i];
}



#ifdef	__cplusplus
}
#endif

#endif	/* LUCOMPLEX_H */

