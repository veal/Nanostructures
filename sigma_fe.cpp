#include "sys_param.h"
void Sigma_fe(double Ef, double T, complex<Doub> Gs[N_En][3][3][N_LAT][N_LAT][N_Band][N_Band],
			  complex<Doub> Sig_fe[3][3][N_LAT][N_LAT][N_Alpha][N_Alpha], double ee[N_En], double Cl[N_Com][N_LAT])
{

	int i = 0;
	int niE_S = 0, niE_E = 0;
	double Kb=0.6346e-5;
	double (*V_f)[N_Com][N_Alpha][N_Band][N_Band] = new double[3][N_Com][N_Alpha][N_Band][N_Band];
	complex<Doub> (*Sigma_fe)[3][3][N_LAT][N_LAT][N_Alpha][N_Alpha] = new complex<Doub>[3][3][3][N_LAT][N_LAT][N_Alpha][N_Alpha];
	double (*CL)[N_Com][N_LAT] = new double[3][N_Com][N_LAT];
	double (*ee_th)[N_En] = new double[3][N_En];
	Sig_fe_param* param = new Sig_fe_param[3];
	string line;
	ifstream file_in("H:\\MyBComp\\Study\\Clean Nanotube\\Progy\\mine\\Matrix_f\\C.dat");
	getline(file_in, line);
	for (int iter2 = 0; iter2 < N_Alpha; iter2++)
	{
		for (int iter1 = 0; iter1 < N_Band; iter1++)
		{
			file_in >> V_f[0][0][iter2][iter1][0] >> V_f[0][0][iter2][iter1][1] >> V_f[0][0][iter2][iter1][2]
			>> V_f[0][0][iter2][iter1][3] >> V_f[0][0][iter2][iter1][4] >> V_f[0][0][iter2][iter1][5]
			>> V_f[0][0][iter2][iter1][6] >> V_f[0][0][iter2][iter1][7] >> V_f[0][0][iter2][iter1][8]
			>> V_f[0][0][iter2][iter1][9];
		}
	}
	file_in.close();
	file_in.open("H:\\MyBComp\\Study\\Clean Nanotube\\Progy\\mine\\Matrix_f\\Cr.dat");
	getline(file_in, line);
	for (int iter2 = 0; iter2 < N_Alpha; iter2++)
	{
		for (int iter1 = 0; iter1 < N_Band; iter1++)
		{
			file_in >> V_f[0][1][iter2][iter1][0] >> V_f[0][1][iter2][iter1][1] >> V_f[0][1][iter2][iter1][2]
			>> V_f[0][1][iter2][iter1][3] >> V_f[0][1][iter2][iter1][4] >> V_f[0][1][iter2][iter1][5]
			>> V_f[0][1][iter2][iter1][6] >> V_f[0][1][iter2][iter1][7] >> V_f[0][1][iter2][iter1][8] 
			>> V_f[0][1][iter2][iter1][9];
		}
	}
	file_in.close();
	for (int it  = 0; it < 2; it++)
	{
		for (int iter2 = 0; iter2 < N_Alpha; iter2++)
		{
			for (int iter1 = 0; iter1 < N_Band; iter1++)
			{
				for (int ib = 0; ib < N_Band; ib++)
				{
					V_f[it+1][0][iter2][iter1][ib] = V_f[0][0][iter2][iter1][ib];
					V_f[it+1][1][iter2][iter1][ib] = V_f[0][1][iter2][iter1][ib];
				}
			}
		}
	}
	while (1.0/(exp((ee[i]-Ef)/(Kb*T))+1.0) > 1.0e-10)
	{
		double d  = 1.0/(exp((ee[i]-Ef)/(Kb*T))+1.0);
		niE_E = i;
		i++;
	}
	for (int it = 0; it < 3; it++)
	{
		for (int i = 0; i < N_Com; i++)
		{
			for (int j = 0; j < N_LAT; j++)
			{
				CL[it][i][j] = Cl[i][j];
			}
		}
		for (int i = 0; i < N_En; i++)
		{
			ee_th[it][i] = ee[i];
		}
		param[it].num = it;
		param[it].niE_S = niE_S + (niE_E - niE_S)*it/3;
		param[it].niE_E = niE_S + (niE_E - niE_S)*(it+1)/3;
		param[it].V_f = V_f[it];
		param[it].Sigma = Sigma_fe[it];
		param[it].Cl = CL[it];
		param[it].G_s = Gs;
		param[it].ee = ee_th[it];
		param[it].T = T;
		param[it].Ef = Ef;
	}
	pthread_t aThread[3];

	for (int iter = 0; iter < 3; iter++)
	{
		pthread_create(&aThread[iter], NULL, &Sigma_fe_thread, NULL);
	}

        for (int iter = 0; iter < 3; iter++)
	{
		pthread_join(aThread[iter], NULL);
	}
	
	for (int n = 0; n < 3; n++)
	{
		for (int ns = 0; ns < 3; ns++)
		{
			for (int il = 0; il < N_LAT; il++)
			{
				for (int jl = 0; jl < N_LAT; jl++)
				{
					for (int ia = 0; ia < N_Alpha; ia++)
					{
						for (int ja = 0; ja < N_Alpha; ja++)
						{
							Sig_fe[n][ns][il][jl][ia][ja] = Sigma_fe[0][n][ns][il][jl][ia][ja] + 
								Sigma_fe[1][n][ns][il][jl][ia][ja] + Sigma_fe[2][n][ns][il][jl][ia][ja];

						}
					}
				}
			}
		}
	}
	delete [] V_f;
	delete [] Sigma_fe;
	delete [] CL;
	delete [] ee_th;
	delete [] param;
}
void* Sigma_fe_thread(void* param)
{
	complex<Doub> Matr1(0.0, 0.0);
	complex<Doub> Matr2(0.0, 0.0);
	complex<Doub> Matr3(0.0, 0.0);
	complex<Doub> Matr4(0.0, 0.0);

	Sig_fe_param* pr = (Sig_fe_param*)param;

	double Kb=0.6346e-5;
	double Koef1;
	double T = pr->T;
	double* ee = pr->ee;
	int niE_S = pr->niE_S;
	int niE_E = pr->niE_E;
	double Ef = pr->Ef;
	double (*V_f)[N_Alpha][N_Band][N_Band] = pr->V_f;
	complex<Doub> (*Gs)[3][3][N_LAT][N_LAT][N_Band][N_Band] = pr->G_s;
	complex<Doub>(*Sig_fe)[3][N_LAT][N_LAT][N_Alpha][N_Alpha] = pr->Sigma;
	double (*Cl)[N_LAT] = pr->Cl;
	for (int n = 0; n < 3; n++)
	{
		for (int ns = 0; ns < 3; ns++)
		{
			for (int il = 0; il < N_LAT; il++)
			{
				for (int jl = 0; jl < N_LAT; jl++)
				{
					for (int ia = 0; ia < N_Alpha; ia++)
					{
						for (int ja = 0; ja < N_Alpha; ja++)
						{
							for (int i = 0; i < N_Com; i++)
							{
								for (int j = 0; j < N_Com; j++)
								{
									complex<Doub> Sig_fe_old(0.0, 0.0);
									complex<Doub> Sig_fe_new(0.0, 0.0);
									complex<Doub> Sig_feij(0.0, 0.0);
									double E = 0.0;
									int N_Sig_fe = niE_E - niE_S;
									double dE = (ee[niE_E] - ee[niE_S])/N_Sig_fe;
									if (T == 0)
									{
										if (ee[niE_S] > Ef) {Koef1 = 0.0;}
										if (ee[niE_S] < Ef) {Koef1 = 1.0;}
									}
									else
									{
										Koef1 = 1.0/(exp((ee[niE_S]-Ef)/(Kb*T))+1.0);
									}

									for (int ib1 = 0; ib1 < N_Band; ib1++)
									{
										for (int ib2 = 0; ib2 < N_Band; ib2++)
										{
											for (int ib3 = 0; ib3 < N_Band; ib3++)
											{
												for (int ib4 = 0; ib4 < N_Band; ib4++)
												{
													Matr1 = 0.0;
													Matr2 = conj(Gs[0][ns][n][jl][il][ib3][ib1])*conj(Gs[0][n][ns][il][jl][ib2][ib4]);
													Matr3 = Gs[0][n][ns][il][jl][ib1][ib3]*Gs[0][ns][n][jl][il][ib4][ib2];
													Matr4 = 0.0;
													Sig_fe_old -= (CI/(2.0*Pi))*Koef1*(V_f[i][ia][ib2][ib1]*
													(Matr1-Matr2+Matr3-Matr4)*V_f[j][ja][ib3][ib4]);
													double b = 0;
												}
											}
										}
									}
									for (int nE2 = niE_S; nE2 < niE_E; nE2++)
									{
										E = ee[niE_S] + nE2*dE;
										if (T == 0)
										{
											 if (E > Ef) {Koef1 = 0.0;}
											 if (E < Ef) {Koef1 = 1.0;}
										}
										else
										{
											Koef1 = 1.0/(exp((E-Ef)/(Kb*T))+1.0);
										}
										for (int ib1 = 0; ib1 < N_Band; ib1++)
										{
											for (int ib2 = 0; ib2 < N_Band; ib2++)
											{
												for (int ib3 = 0; ib3 < N_Band; ib3++)
												{
													for (int ib4 = 0; ib4 < N_Band; ib4++)
													{
														Matr1 = 0.0;
														Matr2 = conj(Gs[nE2][ns][n][jl][il][ib3][ib1])*conj(Gs[nE2][n][ns][il][jl][ib2][ib4]);
														Matr3 = Gs[nE2][n][ns][il][jl][ib1][ib3]*Gs[nE2][ns][n][jl][il][ib4][ib2];
														Matr4 = 0.0;
														Sig_fe_new -= (CI/(2.0*Pi))*Koef1*(V_f[i][ia][ib2][ib1]*
														(Matr1-Matr2+Matr3-Matr4)*V_f[j][ja][ib3][ib4]);
													}
												}
											}
										}
										Sig_feij += 0.5*dE*(Sig_fe_old + Sig_fe_new);

										Sig_fe_old = Sig_fe_new;
										Sig_fe_new = 0.0;
									}
									Sig_fe[n][ns][il][jl][ia][ja] += Cl[i][il]*Cl[j][jl]*Sig_feij;
								}
							}
						}
					}
				}
			}
		}
	}
}