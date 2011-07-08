#include "sys_param.h"
#include "eigen_sym.h"
#include "LUcomplex.h"
void cd_Invert(const Comp a[N_LAT*N_Band][N_LAT*N_Band], Comp aa[N_LAT*N_Band][N_LAT*N_Band])
{
    NRmatrix<Comp> temp(N_LAT*N_Band, N_LAT*N_Band);
    for (int i = 0; i < N_LAT*N_Band; i++) {
        for (int j = 0; j < N_LAT*N_Band; j++) {
            temp[i][j] = a[i][j];
        }
    }

    LUcomplex ytrewq(temp);
    NRmatrix<Comp> result(N_LAT*N_Band, N_LAT*N_Band);
    ytrewq.inverse(result);

    for (int i = 0; i < N_LAT*N_Band; i++) {
        for (int j = 0; j < N_LAT*N_Band; j++) {
            aa[i][j] = result[i][j];
        }
    }

//      Comp a_ll, c_ll[N_LAT*N_Band-1],pr[N_LAT*N_Band];
//      Comp (*cl)[N_LAT*N_Band] = new Comp[N_LAT*N_Band][N_LAT*N_Band];
//      Comp (*cl1)[N_LAT*N_Band] = new Comp[N_LAT*N_Band][N_LAT*N_Band];
//
//      if ((real(a[0][0]) == 0.0)&&(imag(a[0][0]) == 0.0))
//	  {
//		  cout << "Error in Cd_Invert: a(0,0)=0";
//		  return;
//	  }
//	for (int i = 0; i < N_LAT*N_Band; i++)
//	{
//		aa[i][i] = a[i][i];
//		for (int j = i+1; j < N_LAT*N_Band; j++)
//		{
//			aa[i][j] = a[i][j];
//			aa[j][i] = a[j][i];
//			cl1[j][i] = 0;
//		}
//	}
//	for(int l = 0; l < N_LAT*N_Band; l++)
//	{
//		if ((real(aa[l][l]) == 0.0) && (imag(aa[l][l]) == 0.0))
//		{
//			aa[l][l] = Comp(1.0e-31, 1.0e-31);
//		}
//		a_ll = 1.0/aa[l][l];
//		for(int i = l+1; i < N_LAT*N_Band; i++)
//		{
//			c_ll[i-1]=aa[i][l]*a_ll*(-1.0);
//		}
//		for(int j = 0; j < N_LAT*N_Band; j++)
//		{
//			pr[j] = aa[l][j];
//			aa[l][j] = a_ll*pr[j];
//		}
//		for(int i = l+1; i < N_LAT*N_Band; i++)
//		{
//			for(int j = 0; j < N_LAT*N_Band; j++)
//			{
//				aa[i][j] = pr[j]*c_ll[i-1] + aa[i][j];
//			}
//		}
//		cl[l][l]=a_ll;
//		for(int j = 0; j < l-1; j++)
//		{
//			cl[l][j]=a_ll*cl1[l][j];
//		}
//		for(int i = l+1; i < N_LAT*N_Band; i++)
//		{
//			cl[i][l]=c_ll[i];
//			for(int j = 0; j < l-1; j++)
//			{
//				cl[i][j]=cl1[i][j]+cl1[l][j]*c_ll[i-1];
//			}
//		}
//		for(int i = l; i < N_LAT*N_Band; i++)
//		{
//			for(int j = 0; j < i; j++)
//			{
//				cl1[i][j]=cl[i][j];
//			}
//		}
//	}
//	cl[N_LAT*N_Band-1][N_LAT*N_Band-1] = 1.0;
//	for(int i = N_LAT*N_Band-2; i >= 1; i--)
//	{
//		cl[i][i] = 1.0;
//		for(int j = i+1; j < N_LAT*N_Band; j++)
//		{
//			cl[j][i] = 0.0;
//			cl[i][j] = 0.0;
//			for(int k = i+1; k <= j; k++)
//			{
//				cl[i][j] = cl[i][j] + cl[k][j]*(-1.0)*aa[i][k];
//			}
//		}
//	}
//	for (int i = 0; i < N_LAT*N_Band; i++)
//	{
//		for(int j = 0; j < N_LAT*N_Band; j++)
//		{
//			aa[i][j] = 0.0;
//			for(int k = j; k < N_LAT*N_Band; k++)
//			{
//				aa[i][j]=aa[i][j]+cl[i][k]*cl1[k][j];
//			}
//		}
//	}
//        delete[] cl;
//        delete[] cl1;
}

void cd_Invert(const Comp** a, Comp** aa, int dimension) {
    Comp* ta = (Comp*)a;
    Comp* taa = (Comp*)aa;
    NRmatrix<Comp> temp(dimension, dimension);
    for (int i = 0; i < dimension; i++) {
        for (int j = 0; j < dimension; j++) {
            temp[i][j] = ta[i*dimension+j];
        }
    }

    LUcomplex ytrewq(temp);
    NRmatrix<Comp> result(dimension, dimension);
    ytrewq.inverse(result);

    for (int i = 0; i < dimension; i++) {
        for (int j = 0; j < dimension; j++) {
            taa[i*dimension+j] = result[i][j];
        }
    }
}

void matmul(Comp A[N_LAT*N_Band][N_LAT*N_Band], Comp B[N_LAT*N_Band][N_LAT*N_Band],
			Comp res[N_LAT*N_Band][N_LAT*N_Band])
{
	for(int i = 0; i < N_LAT*N_Band; i++)
	{
		for(int j = 0; j < N_LAT*N_Band; j++)
		{
			res[i][j] = 0;
			for (int it = 0; it < N_LAT*N_Band; it++)
			{
				res[i][j] += A[i][it]*B[it][j];
			}
		}
	}
}
void matmul3(Comp A[N_Alpha][N_Alpha], Comp B[N_Alpha][N_Alpha], Comp res[N_Alpha][N_Alpha])
{
	for(int i = 0; i < N_LAT*N_Band; i++)
	{
		for(int j = 0; j < N_LAT*N_Band; j++)
		{
			res[i][j] = 0;
			for (int it = 0; it < N_LAT*N_Band; it++)
			{
				res[i][j] += A[i][it]*B[it][j];
			}
		}
	}
}
void matmul36(Comp A[N_LAT*N_Alpha][N_LAT*N_Alpha], Comp B[N_LAT*N_Alpha][N_LAT*N_Alpha],
			  Comp res[N_LAT*N_Alpha][N_LAT*N_Alpha])
{
	for(int i = 0; i < N_LAT*N_Alpha; i++)
	{
		for(int j = 0; j < N_LAT*N_Alpha; j++)
		{
			res[i][j] = 0;
			for (int it = 0; it < N_LAT*N_Alpha; it++)
			{
				res[i][j] += A[i][it]*B[it][j];
			}
		}
	}
}
void matmul4(Comp A[N_Band][N_Band], Comp B[N_Band][N_Band], Comp res[N_Band][N_Band])
{
	for(int i = 0; i < N_Band; i++)
	{
		for(int j = 0; j < N_Band; j++)
		{
			res[i][j] = 0;
			for (int it = 0; it < N_Band; it++)
			{
				res[i][j] += A[i][it]*B[it][j];
			}
		}
	}
}
void matinv(Comp B[N_LAT*N_Band][N_LAT*N_Band],Comp res[N_LAT*N_Band][N_LAT*N_Band])
{
	Comp (*A)[N_LAT*N_Band] = new Comp[N_LAT*N_Band][N_LAT*N_Band];
	for(int i = 0; i < N_LAT*N_Band; i++)
	{
		for(int j = 0; j < N_LAT*N_Band; j++)
		{
			res[i][j] = (i == j)?1:0;
			A[i][j] = B[i][j];
		}
	}
	for(int it = 0; it < N_LAT*N_Band; it++)
	{
		if (abs(A[it][it]) == 0)
		{
			cout << it << '\n';
			A[it][it] = 1e-20;
//			return;
		}
		Comp del = A[it][it];
		for(int it1 = 0; it1 < N_LAT*N_Band; it1++)
		{
			res[it][it1] /= del;
			A[it][it1] /= del;
		}
		for(int it1 = 1; it1 < N_LAT*N_Band-it; it1++)
		{
			Comp dell = A[it+it1][it];
			if ((it+it1 < 0) || (it+it1 > 119))
			{
				cout << it << ' ' << it1 << '\n';
			}
			for(int it2 = 0; it2 < N_LAT*N_Band; it2++)
			{
				res[it+it1][it2] -= res[it][it2]*dell;
				A[it+it1][it2] -= A[it][it2]*dell;
			}
		}
	}
	for(int it = 0; it < N_LAT*N_Band; it++)
	{
		for(int it1 = 1; it1 < N_LAT*N_Band-it; it1++)
		{
			Comp dell = A[N_LAT*N_Band-1-it-it1][N_LAT*N_Band-1-it];
			for(int it2 = 0; it2 < N_LAT*N_Band; it2++)
			{
				if ((N_LAT*N_Band-1-it-it1 < 0) || (N_LAT*N_Band-1-it-it1 > 119))
				{
					cout << it << ' ' << it1 << '\n';
				}
				res[N_LAT*N_Band-1-it-it1][it2] -= res[N_LAT*N_Band-1-it][it2]*dell;
				A[N_LAT*N_Band-1-it-it1][it2] -= A[N_LAT*N_Band-1-it][it2]*dell;
			}
		}
	}
	delete [] A;
}
void matinv36(Comp B[N_LAT*N_Alpha][N_LAT*N_Alpha],Comp res[N_LAT*N_Alpha][N_LAT*N_Alpha])
{
	Comp (*A)[N_LAT*N_Alpha] = new Comp[N_LAT*N_Alpha][N_LAT*N_Alpha];
	for(int i = 0; i < N_LAT*N_Alpha; i++)
	{
		for(int j = 0; j < N_LAT*N_Alpha; j++)
		{
			res[i][j] = (i == j)?1:0;
			A[i][j] = B[i][j];
		}
	}
	for(int it = 0; it < N_LAT*N_Alpha; it++)
	{
		if (abs(A[it][it]) == 0)
		{
			cout << it << '\n';
			A[it][it] = 1e-20;
//			return;
		}
		Comp del = A[it][it];
		for(int it1 = 0; it1 < N_LAT*N_Alpha; it1++)
		{
			res[it][it1] /= del;
			A[it][it1] /= del;
		}
		for(int it1 = 1; it1 < N_LAT*N_Alpha-it; it1++)
		{
			Comp dell = A[it+it1][it];
			if ((it+it1 < 0) || (it+it1 > 119))
			{
				cout << it << ' ' << it1 << '\n';
			}
			for(int it2 = 0; it2 < N_LAT*N_Alpha; it2++)
			{
				res[it+it1][it2] -= res[it][it2]*dell;
				A[it+it1][it2] -= A[it][it2]*dell;
			}
		}
	}
	for(int it = 0; it < N_LAT*N_Alpha; it++)
	{
		for(int it1 = 1; it1 < N_LAT*N_Alpha-it; it1++)
		{
			Comp dell = A[N_LAT*N_Alpha-1-it-it1][N_LAT*N_Alpha-1-it];
			for(int it2 = 0; it2 < N_LAT*N_Alpha; it2++)
			{
				if ((N_LAT*N_Alpha-1-it-it1 < 0) || (N_LAT*N_Alpha-1-it-it1 > 119))
				{
					cout << it << ' ' << it1 << '\n';
				}
				res[N_LAT*N_Alpha-1-it-it1][it2] -= res[N_LAT*N_Alpha-1-it][it2]*dell;
				A[N_LAT*N_Alpha-1-it-it1][it2] -= A[N_LAT*N_Alpha-1-it][it2]*dell;
			}
		}
	}
	delete [] A;
}
void matinv4(Comp B[N_Band][N_Band],Comp res[N_Band][N_Band])
{
	Comp (*A)[N_Band] = new Comp[N_Band][N_Band];
	for(int i = 0; i < N_Band; i++)
	{
		for(int j = 0; j < N_Band; j++)
		{
			res[i][j] = (i == j)?1:0;
			A[i][j] = B[i][j];
		}
	}
	for(int it = 0; it < N_Band; it++)
	{
		if (abs(A[it][it]) == 0)
		{
			cout << it << '\n';
			A[it][it] = 1e-20;
//			return;
		}
		Comp del = A[it][it];
		for(int it1 = 0; it1 < N_Band; it1++)
		{
			res[it][it1] /= del;
			A[it][it1] /= del;
		}
		for(int it1 = 1; it1 < N_Band-it; it1++)
		{
			Comp dell = A[it+it1][it];
			if ((it+it1 < 0) || (it+it1 > 119))
			{
				cout << it << ' ' << it1 << '\n';
			}
			for(int it2 = 0; it2 < N_Band; it2++)
			{
				res[it+it1][it2] -= res[it][it2]*dell;
				A[it+it1][it2] -= A[it][it2]*dell;
			}
		}
	}
	for(int it = 0; it < N_Band; it++)
	{
		for(int it1 = 1; it1 < N_Band-it; it1++)
		{
			Comp dell = A[N_Band-1-it-it1][N_Band-1-it];
			for(int it2 = 0; it2 < N_Band; it2++)
			{
				if ((N_Band-1-it-it1 < 0) || (N_Band-1-it-it1 > 119))
				{
					cout << it << ' ' << it1 << '\n';
				}
				res[N_Band-1-it-it1][it2] -= res[N_Band-1-it][it2]*dell;
				A[N_Band-1-it-it1][it2] -= A[N_Band-1-it][it2]*dell;
			}
		}
	}
	delete [] A;
}
void matinv3(Comp B[N_Alpha][N_Alpha],Comp res[N_Alpha][N_Alpha])
{
	Comp (*A)[N_Alpha] = new Comp[N_Alpha][N_Alpha];
	for(int i = 0; i < N_Alpha; i++)
	{
		for(int j = 0; j < N_Alpha; j++)
		{
			res[i][j] = (i == j)?1:0;
			A[i][j] = B[i][j];
		}
	}
	for(int it = 0; it < N_Alpha; it++)
	{
		if (abs(A[it][it]) == 0)
		{
			cout << it << '\n';
			A[it][it] = 1e-20;
//			return;
		}
		Comp del = A[it][it];
		for(int it1 = 0; it1 < N_Alpha; it1++)
		{
			res[it][it1] /= del;
			A[it][it1] /= del;
		}
		for(int it1 = 1; it1 < N_Alpha-it; it1++)
		{
			Comp dell = A[it+it1][it];
			if ((it+it1 < 0) || (it+it1 > 119))
			{
				cout << it << ' ' << it1 << '\n';
			}
			for(int it2 = 0; it2 < N_Alpha; it2++)
			{
				res[it+it1][it2] -= res[it][it2]*dell;
				A[it+it1][it2] -= A[it][it2]*dell;
			}
		}
	}
	for(int it = 0; it < N_Alpha; it++)
	{
		for(int it1 = 1; it1 < N_Alpha-it; it1++)
		{
			Comp dell = A[N_Alpha-1-it-it1][N_Alpha-1-it];
			for(int it2 = 0; it2 < N_Alpha; it2++)
			{
				if ((N_Alpha-1-it-it1 < 0) || (N_Alpha-1-it-it1 > 119))
				{
					cout << it << ' ' << it1 << '\n';
				}
				res[N_Alpha-1-it-it1][it2] -= res[N_Alpha-1-it][it2]*dell;
				A[N_Alpha-1-it-it1][it2] -= A[N_Alpha-1-it][it2]*dell;
			}
		}
	}
	delete [] A;
}
void Trans_Hk(Comp S[N_LAT][N_LAT][N_Band][N_Band], Comp H[N_LAT][N_LAT][N_Band][N_Band],
			  Comp H1[N_LAT][N_LAT][N_Band][N_Band], bool isOutputNeeded) {

    MatDoub temp(2 * N_LAT * N_Band, 2 * N_LAT * N_Band);
    Comp (*hamilt)[N_LAT*N_Band] = new Comp[N_LAT*N_Band][N_LAT*N_Band];

    for (int i = 0; i < N_LAT; i++) {
        for (int j = 0; j < N_LAT; j++) {
            for (int i1 = 0; i1 < N_Band; i1++) {
                for (int j1 = 0; j1 < N_Band; j1++) {
                    temp[i * N_Band + i1][j * N_Band + j1] = real(S[i][j][i1][j1]);
                    hamilt[i * N_Band + i1][j * N_Band + j1] = H[i][j][i1][j1];
//                    if (conj(H[i][j][i1][j1]) != H[j][i][j1][i1])
//                        cout << i * N_Band + i1 << "  " << j * N_Band + j1 << '\n';
                }
            }
        }
    }

    for (int i = 0; i < N_LAT; i++) {
        for (int j = 0; j < N_LAT; j++) {
            for (int i1 = 0; i1 < N_Band; i1++) {
                for (int j1 = 0; j1 < N_Band; j1++) {
                    temp[i * N_Band + i1][N_LAT*N_Band + j * N_Band + j1] = -imag(S[i][j][i1][j1]);
                }
            }
        }
    }

    for (int i = 0; i < N_LAT; i++) {
        for (int j = 0; j < N_LAT; j++) {
            for (int i1 = 0; i1 < N_Band; i1++) {
                for (int j1 = 0; j1 < N_Band; j1++) {
                    temp[N_LAT*N_Band + i * N_Band + i1][j * N_Band + j1] = imag(S[i][j][i1][j1]);
                }
            }
        }
    }
    for (int i = 0; i < N_LAT; i++) {
        for (int j = 0; j < N_LAT; j++) {
            for (int i1 = 0; i1 < N_Band; i1++) {
                for (int j1 = 0; j1 < N_Band; j1++) {
                    temp[N_LAT*N_Band + i * N_Band + i1][N_LAT*N_Band + j * N_Band + j1] = real(S[i][j][i1][j1]);
                }
            }
        }
    }

    Symmeig qwerty(temp);

    Comp (*middle)[N_LAT * N_Band] = new Comp[N_LAT * N_Band][N_LAT * N_Band];
    Comp (*middle2)[N_LAT * N_Band] = new Comp[N_LAT * N_Band][N_LAT * N_Band];
    Comp (*first)[N_LAT * N_Band] = new Comp[N_LAT * N_Band][N_LAT * N_Band];
    Comp (*last)[N_LAT * N_Band] = new Comp[N_LAT * N_Band][N_LAT * N_Band];
    Comp (*temp_matr)[N_LAT * N_Band] = new Comp[N_LAT * N_Band][N_LAT * N_Band];

    ofstream file("Rezults/S_k.dat", ios::app);
    file.precision(5);

    for (int i = 0; i < N_LAT * N_Band; i++) {
        for (int i1 = 0; i1 < N_LAT * N_Band; i1++) {
            first[i][i1] = Comp(qwerty.z[i][2*i1], qwerty.z[N_LAT*N_Band+i][2*i1]);
            middle[i][i1] = (i == i1) ? Comp(1.0, 0.0)/sqrt(Comp(qwerty.d[2*i], 0.0)) : 0.0;
            middle2[i][i1] = (i == i1) ? sqrt(Comp(qwerty.d[2*i], 0.0)) : 0.0;
            last[i1][i] = conj(first[i][i1]);
        }
        if (isOutputNeeded)
            file << scientific << qwerty.d[2*i] << '\t';
    }
    if (isOutputNeeded)
        file << '\n';
    file.close();

    for (int i = 0; i < N_LAT * N_Band; i++) {
        for (int i1 = 0; i1 < N_LAT * N_Band; i1++) {
            Comp result = 0.0;
            for (int j = 0; j < N_LAT * N_Band; j++) {
                result += middle2[i][j] * last[j][i1];
            }
            temp_matr[i][i1] = result;
        }
    }

    for (int i = 0; i < N_LAT * N_Band; i++) {
        for (int i1 = 0; i1 < N_LAT * N_Band; i1++) {
            Comp result = 0.0;
            for (int j = 0; j < N_LAT * N_Band; j++) {
                result += first[i][j] * temp_matr[j][i1];
            }
            middle[i][i1] = result;
        }
    }

    for (int i = 0; i < N_LAT; i++) {
        for (int k = 0; k < N_LAT; k++) {
            for (int i1 = 0; i1 < N_Band; i1++) {
                for (int k1 = 0; k1 < N_Band; k1++) {
                    S[i][k][i1][k1] = middle[i*N_Band + i1][k*N_Band + k1];
//                    Comp result = 0.0;
//                    for (int j = 0; j < N_LAT; j++) {
//                        for (int j1 = 0; j1 < N_Band; j1++) {
//                            result += middle[i*N_Band + i1][j*N_Band + j1] * middle[j*N_Band + j1][k*N_Band + k1];
//                        }
//                    }
//                    if(abs(S[i][k][i1][k1] - result) > 0.01*abs(S[i][k][i1][k1])/*)//*/ && abs(S[i][k][i1][k1]) != 0.0)
//                        cout << i << " " << k << " " << i1 << " " << k1 << " " <<
//                                abs(S[i][k][i1][k1] - result)/abs(S[i][k][i1][k1]) << '\n';
////                                real(result) << "  " << imag(result) << '\n';
                }
            }
        }
    }

//    checkIfMatrixIsHermitian(S, 0.01);

//    (S-1/2) matrix complete

    for (int i = 0; i < N_LAT * N_Band; i++) {
        for (int i1 = 0; i1 < N_LAT * N_Band; i1++) {
            Comp result = 0.0;
            for (int j = 0; j < N_LAT * N_Band; j++) {
                result += hamilt[i][j] * middle[j][i1];
            }
            last[i][i1] = result;
        }
    }
//    H*(S-1/2) matrix complete

    for (int i = 0; i < N_LAT * N_Band; i++) {
        for (int i1 = 0; i1 < N_LAT * N_Band; i1++) {
            Comp result = 0.0;
            for (int j = 0; j < N_LAT * N_Band; j++) {
                result += middle[i][j] * last[j][i1];
            }
            hamilt[i][i1] = result;
        }
    }
//    (S-1/2)*H*(S-1/2) matrix complete

    for (int i = 0; i < N_LAT; i++) {
        for (int j = 0; j < N_LAT; j++) {
            for (int i1 = 0; i1 < N_Band; i1++) {
                for (int j1 = 0; j1 < N_Band; j1++) {
                    H1[i][j][i1][j1] = hamilt[i * N_Band + i1][j * N_Band + j1];
                }
            }
        }
    }

    delete first;
    delete middle;
    delete last;
    delete temp_matr;
    delete hamilt;

//    checkIfMatrixIsHermitian(H1, 0.01);

//	Comp (*X_tem)[N_LAT*N_Band] = new Comp[N_LAT*N_Band][N_LAT*N_Band];
//	Comp (*H_tem)[N_LAT*N_Band] = new Comp[N_LAT*N_Band][N_LAT*N_Band];
//	Comp (*T_tem)[N_LAT*N_Band] = new Comp[N_LAT*N_Band][N_LAT*N_Band];
//	Comp (*Y_tem)[N_LAT*N_Band] = new Comp[N_LAT*N_Band][N_LAT*N_Band];
//	Comp (*X)[N_LAT][N_Band][N_Band] = new Comp[N_LAT][N_LAT][N_Band][N_Band];
//	Comp (*Y)[N_LAT][N_Band][N_Band] = new Comp[N_LAT][N_LAT][N_Band][N_Band];
//	for(int n = 0; n < N_LAT; n++)
//	{
//		for(int n1 = 0; n1 < N_LAT; n1++)
//		{
//			for(int gm = 0; gm < N_Band; gm++)
//			{
//				for(int gm1 = 0; gm1 < N_Band; gm1++)
//				{
//					X[n][n1][gm][gm1] = S[n][n1][gm][gm1];
////					X_tem[n*N_Band+gm][n1*N_Band+gm1]=X[n][n1][gm][gm1];
//					T_tem[n*N_Band+gm][n1*N_Band+gm1]=X[n][n1][gm][gm1];
//					H_tem[n*N_Band+gm][n1*N_Band+gm1]=H[n][n1][gm][gm1];
//				}
//			}
//		}
//	}
////	matinv(X_tem, T_tem);
//	matmul(T_tem, H_tem, X_tem);
//	matmul(H_tem, T_tem, Y_tem);
//	for(int n = 0; n < N_LAT; n++)
//	{
//		for(int n1 = 0; n1 < N_LAT; n1++)
//		{
//			for(int gm = 0; gm < N_Band; gm++)
//			{
//				for(int gm1 = 0; gm1 < N_Band; gm1++)
//				{
//					H1[n][n1][gm][gm1] -= 0.5*(X_tem[n*N_Band+gm][n1*N_Band+gm1]+Y_tem[n*N_Band+gm][n1*N_Band+gm1]);
//				}
//			}
//		}
//	}
//	delete [] X_tem;
//	delete [] H_tem;
//	delete [] T_tem;
//	delete [] Y_tem;
//	delete [] X;
//	delete [] Y;
	
}

void Eigen_values(Comp S[N_LAT][N_LAT][N_Band][N_Band], double* Eigen_val, bool isOutputNeeded) {

    MatDoub temp(2 * N_LAT * N_Band, 2 * N_LAT * N_Band);

    for (int i = 0; i < N_LAT; i++) {
        for (int j = 0; j < N_LAT; j++) {
            for (int i1 = 0; i1 < N_Band; i1++) {
                for (int j1 = 0; j1 < N_Band; j1++) {
                    temp[i * N_Band + i1][j * N_Band + j1] = real(S[i][j][i1][j1]);
                }
            }
        }
    }

    for (int i = 0; i < N_LAT; i++) {
        for (int j = 0; j < N_LAT; j++) {
            for (int i1 = 0; i1 < N_Band; i1++) {
                for (int j1 = 0; j1 < N_Band; j1++) {
                    temp[i * N_Band + i1][N_LAT*N_Band + j * N_Band + j1] = -imag(S[i][j][i1][j1]);
                }
            }
        }
    }

    for (int i = 0; i < N_LAT; i++) {
        for (int j = 0; j < N_LAT; j++) {
            for (int i1 = 0; i1 < N_Band; i1++) {
                for (int j1 = 0; j1 < N_Band; j1++) {
                    temp[N_LAT*N_Band + i * N_Band + i1][j * N_Band + j1] = imag(S[i][j][i1][j1]);
                }
            }
        }
    }

    for (int i = 0; i < N_LAT; i++) {
        for (int j = 0; j < N_LAT; j++) {
            for (int i1 = 0; i1 < N_Band; i1++) {
                for (int j1 = 0; j1 < N_Band; j1++) {
                    temp[N_LAT*N_Band + i * N_Band + i1][N_LAT*N_Band + j * N_Band + j1] = real(S[i][j][i1][j1]);
                }
            }
        }
    }

    Symmeig qwerty(temp);

    ofstream file("Rezults/Matrix_values.dat", ios::app);
    file.precision(5);

    for (int i = 0; i < N_LAT * N_Band; i++) {
        if (Eigen_val != NULL)
            Eigen_val[i] = qwerty.d[2*i];
        if (isOutputNeeded)
            file << scientific << qwerty.d[2*i] << '\t';
//        cout << scientific << qwerty.d[2*i] << '\n';
    }
    if (isOutputNeeded)
        file << '\n';
    file.close();

}

Doub Fun_Fermi(Doub E, Doub Ef, Doub T)
{
	Doub res;
	if (T == 0.0)
	{
		if (E-Ef < 0.0) 
		{
			res = 1.0;
		}
		else
		{
			res = 0.0;
		}
	}
	else
	{
		Doub x = 1.0/(0.634552e-5*T);
		res = 1.0/(1.0 + exp((E-Ef)*x));
	}
	return res;
}
Doub DFun_Fermi(Doub E, Doub Ef, Doub kT)
{
	Doub x, temp;
	if (kT == 0.0)
	{
		cout << "Error: T=0 & dF/dE -> infinity" << '\n';
	}

	x = 1/(0.634552e-5*kT); // 1/(Kb*T) v Ry

	if ((E-Ef)*x > 150.0 || (E-Ef)*x < -150.0)
	{
		return 0.0;
	}
	else
	{
		temp = (1 + exp((E-Ef)*x))*(1 + exp((E-Ef)*x));
		return -x*exp((E-Ef)*x)/(temp);
	}
}
void Int_Simpson(Doub &Int, int i, Doub y_2, Doub y_1, Doub y, Doub h)
{
//**************
	if((i % 2) != 0.0)
	{
		Int += dI(y_2,y_1,y,h);
	}
}
Doub dI(Doub y_2, Doub y_1, Doub y, Doub h)
{
	return (y_2 + 4.0*y_1 + y)*h/3.0;
}

void checkIfMatrixIsHermitian(Comp A[N_LAT][N_LAT][N_Band][N_Band], Doub ACCURACY) {
    for (int i = 0; i < N_LAT; i++) {
        for (int j = 0; j < N_LAT; j++) {
            for (int i1 = 0; i1 < N_Band; i1++) {
                for (int j1 = 0; j1 < N_Band; j1++) {
                    if (abs(real(A[i][j][i1][j1]) - real(A[j][i][j1][i1])) > ACCURACY*abs(real(A[i][j][i1][j1]))) {
                        if(abs(A[i][j][i1][j1]) > 10e-10) {
                            cout << "Matrix is not Hermitian" << '\n';
                            exit(-1);
                        } else {
                            A[i][j][i1][j1] = 0.0;
                        }
                    }
                    if (abs(imag(A[i][j][i1][j1]) + imag(A[j][i][j1][i1])) > ACCURACY*abs(imag(A[i][j][i1][j1]))) {
                        if(abs(A[i][j][i1][j1]) > 10e-10) {
                            cout << "Matrix is not Hermitian" << '\n';
                            exit(-1);
                        } else {
                            A[i][j][i1][j1] = 0.0;
                        }
                    } else {
                        A[i][j][i1][j1] = conj(A[j][i][j1][i1]);
                    }

                }
            }
        }
    }
}