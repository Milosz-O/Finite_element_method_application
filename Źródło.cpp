#include <fstream>
#include <iostream>
#include <vector>
#include <string>
#include <iomanip>

using namespace std;

double dN_dksi_1(double eta)
{
	return -1.0 / 4.0 * (1 - eta);
};
double dN_dksi_2(double eta)
{
	return 1.0 / 4.0 * (1 - eta);
};
double dN_dksi_3(double eta)
{
	return 1.0 / 4.0 * (1 + eta);
};
double dN_dksi_4(double eta)
{
	return -1.0 / 4.0 * (1 + eta);
};

double dN_deta_1(double ksi)
{
	return -1.0 / 4.0 * (1 - ksi);
};
double dN_deta_2(double ksi)
{
	return -1.0 / 4.0 * (1 + ksi);
};
double dN_deta_3(double ksi)
{
	return 1.0 / 4.0 * (1 + ksi);
};
double dN_deta_4(double ksi)
{
	return 1.0 / 4.0 * (1 - ksi);
};

struct Global_data
{
	int simulationTime;
	int simulationStepTime;
	int conductivity;
	int alfa;
	int tot;
	int initialTemp;
	int density;
	int specificHeat;
	int nodesNumber;
	int elementsNumber;
};
struct Element
{
	double** tab_dN_po_dksi;
	double** tab_dN_po_deta;

	Element(int n)
	{
		int h = n * n;		//iloœæ punktów ca³kowania

		tab_dN_po_dksi = new double* [h];
		for (int i = 0; i < h; i++)
		{
			tab_dN_po_dksi[i] = new double[4];
		}

		tab_dN_po_deta = new double* [h];
		for (int i = 0; i < h; i++)
		{
			tab_dN_po_deta[i] = new double[4];
		}

		if (n == 2)
		{
			for (int i = 0; i < h; i++)
			{
				double eta[4] = { -1.0 / sqrt(3.0) , -1.0 / sqrt(3.0) , 1.0 / sqrt(3.0) , 1.0 / sqrt(3.0) };

				tab_dN_po_dksi[i][0] = dN_dksi_1(eta[i]);
				tab_dN_po_dksi[i][1] = dN_dksi_2(eta[i]);
				tab_dN_po_dksi[i][2] = dN_dksi_3(eta[i]);
				tab_dN_po_dksi[i][3] = dN_dksi_4(eta[i]);
			}

			for (int i = 0; i < h; i++)
			{
				double ksi[4] = { -1.0 / sqrt(3.0) , 1.0 / sqrt(3.0) , -1.0 / sqrt(3.0) , 1.0 / sqrt(3.0) };

				tab_dN_po_deta[i][0] = dN_deta_1(ksi[i]);
				tab_dN_po_deta[i][1] = dN_deta_2(ksi[i]);
				tab_dN_po_deta[i][2] = dN_deta_3(ksi[i]);
				tab_dN_po_deta[i][3] = dN_deta_4(ksi[i]);
			}
		}
		if (n == 3)
		{
			for (int i = 0; i < h; i++)
			{
				double eta[9] = { -sqrt(3.0 / 5.0) , -sqrt(3.0 / 5.0) , -sqrt(3.0 / 5.0) , 0 , 0 , 0 , sqrt(3.0 / 5.0) , sqrt(3.0 / 5.0) , sqrt(3.0 / 5.0) };

				tab_dN_po_dksi[i][0] = dN_dksi_1(eta[i]);
				tab_dN_po_dksi[i][1] = dN_dksi_2(eta[i]);
				tab_dN_po_dksi[i][2] = dN_dksi_3(eta[i]);
				tab_dN_po_dksi[i][3] = dN_dksi_4(eta[i]);
			}

			for (int i = 0; i < h; i++)
			{
				double ksi[9] = { -sqrt(3.0 / 5.0) , 0 , sqrt(3.0 / 5.0) , -sqrt(3.0 / 5.0) , 0 , sqrt(3.0 / 5.0) , -sqrt(3.0 / 5.0) , 0 , sqrt(3.0 / 5.0) };

				tab_dN_po_deta[i][0] = dN_deta_1(ksi[i]);
				tab_dN_po_deta[i][1] = dN_deta_2(ksi[i]);
				tab_dN_po_deta[i][2] = dN_deta_3(ksi[i]);
				tab_dN_po_deta[i][3] = dN_deta_4(ksi[i]);
			}
		}
		if (n == 4)
		{
			for (int i = 0; i < h; i++)
			{
				double eta[16] = { -0.861136 , -0.861136 , -0.861136 , -0.861136 , -0.339981 , -0.339981 , -0.339981 , -0.339981 , 0.339981 , 0.339981 , 0.339981 , 0.339981 , 0.861136 , 0.861136 , 0.861136 , 0.861136 };

				tab_dN_po_dksi[i][0] = dN_dksi_1(eta[i]);
				tab_dN_po_dksi[i][1] = dN_dksi_2(eta[i]);
				tab_dN_po_dksi[i][2] = dN_dksi_3(eta[i]);
				tab_dN_po_dksi[i][3] = dN_dksi_4(eta[i]);
			}

			for (int i = 0; i < h; i++)
			{
				double ksi[16] = { -0.861136 , -0.339981 , 0.339981 , 0.861136 , -0.861136 , -0.339981 , 0.339981 , 0.861136 , -0.861136 , -0.339981 , 0.339981 , 0.861136 , -0.861136 , -0.339981 , 0.339981 , 0.861136 };

				tab_dN_po_deta[i][0] = dN_deta_1(ksi[i]);
				tab_dN_po_deta[i][1] = dN_deta_2(ksi[i]);
				tab_dN_po_deta[i][2] = dN_deta_3(ksi[i]);
				tab_dN_po_deta[i][3] = dN_deta_4(ksi[i]);
			}
		}
	};
};
struct Node
{
	double x, y = 0;
	int BC = 0;
};
struct ElementID
{
	int ID[4] = { 0, 0, 0, 0 };
};
struct Grid
{
	vector<Node> ND;
	vector<ElementID> EL;
};
struct Uklad
{
	double** macierz_duza;
	double** macierz_Hbc;
	double** macierz_suma;
	double** macierz_C;
	double** macierz_H_plus_C_przez_deltaT;
	double* wektorP_duzy;
	double* wektor_t0;
	double* wektor_do_ostatniego_rownania;

	Uklad(double*** macierze_h_dla_elementow, double*** macierze_Hbc_dla_elementow, Grid grid, double** wektoryP_dla_elementow, Global_data gd, double*** macierze_C_dla_elementow)
	{
		macierz_duza = new double* [gd.nodesNumber];
		macierz_Hbc = new double* [gd.nodesNumber];
		macierz_suma = new double* [gd.nodesNumber];
		macierz_C = new double* [gd.nodesNumber];
		macierz_H_plus_C_przez_deltaT = new double* [gd.nodesNumber];
		wektorP_duzy = new double[gd.nodesNumber];
		wektor_t0 = new double[gd.nodesNumber];
		wektor_do_ostatniego_rownania = new double[gd.nodesNumber];

		for (int i = 0; i < gd.nodesNumber; i++)
		{
			macierz_duza[i] = new double[gd.nodesNumber];
			macierz_Hbc[i] = new double[gd.nodesNumber];
			macierz_suma[i] = new double[gd.nodesNumber];
			macierz_C[i] = new double[gd.nodesNumber];
			macierz_H_plus_C_przez_deltaT[i] = new double[gd.nodesNumber];
		}
		for (int i = 0; i < gd.nodesNumber; i++)
		{
			for (int j = 0; j < gd.nodesNumber; j++)
			{
				macierz_duza[i][j] = 0;
				macierz_Hbc[i][j] = 0;
				macierz_suma[i][j] = 0;
				macierz_C[i][j] = 0;
				macierz_H_plus_C_przez_deltaT[i][j] = 0;
			}
			cout << endl;
		}

		int id_temp = 0;

		for (int u = 0; u < gd.elementsNumber; u++)
		{
			for (int i = 0; i < 4; i++)
			{
				id_temp = grid.EL[u].ID[i] - 1;

				for (int j = 0; j < 4; j++)
				{
					macierz_duza[id_temp][grid.EL[u].ID[j] - 1] += macierze_h_dla_elementow[u][i][j];
					macierz_Hbc[id_temp][grid.EL[u].ID[j] - 1] += macierze_Hbc_dla_elementow[u][i][j];
					macierz_C[id_temp][grid.EL[u].ID[j] - 1] += macierze_C_dla_elementow[u][i][j];
				}

				id_temp = 0;
			}
		}

		for (int i = 0; i < gd.nodesNumber; i++)
		{
			wektorP_duzy[i] = 0;
			wektor_t0[i] = gd.initialTemp;
			wektor_do_ostatniego_rownania[i] = 0;
		}

		for (int i = 0; i < gd.elementsNumber; i++)
		{
			for (int j = 0; j < 4; j++)
			{
				id_temp = grid.EL[i].ID[j] - 1;
				wektorP_duzy[id_temp] += wektoryP_dla_elementow[i][j];
			}
		}
		/*
		cout << "DUZY WEKTOR P" << endl << endl;
		for (int i = 0; i < gd.nodesNumber; i++)
		{
			cout << wektorP_duzy[i] << "	";
		}
		cout << endl << endl;

		cout << "DUZA MACIERZ H" << endl << endl;
		for (int i = 0; i < gd.nodesNumber; i++)
		{
			for (int j = 0; j < gd.nodesNumber; j++)
			{
				cout << macierz_duza[i][j] << " ";
			}
			cout << endl;
		}
		cout << endl << endl;

		cout << "DUZA MACIERZ Hbc" << endl << endl;
		for (int i = 0; i < gd.nodesNumber; i++)
		{
			for (int j = 0; j < gd.nodesNumber; j++)
			{
				cout << macierz_Hbc[i][j] << " ";
			}
			cout << endl;
		}
		cout << endl << endl;
		*/

		for (int i = 0; i < gd.nodesNumber; i++)
		{
			for (int j = 0; j < gd.nodesNumber; j++)
			{
				macierz_suma[i][j] = macierz_Hbc[i][j] + macierz_duza[i][j];
			}
		}
		/*
		cout << "ZSUMOWANE MACIERZE" << endl << endl;
		for (int i = 0; i < gd.nodesNumber; i++)
		{
			for (int j = 0; j < gd.nodesNumber; j++)
			{
				cout << macierz_suma[i][j] << " ";
			}
			cout << endl;
		}
		*/
		/*
		cout << "DUZA MACIERZ C" << endl << endl;
		for (int i = 0; i < gd.nodesNumber; i++)
		{
			for (int j = 0; j < gd.nodesNumber; j++)
			{
				cout << macierz_C[i][j] << " ";
			}
			cout << endl;
		}
		*/
	}
};

double* obliczanie_rownania(double** macierz, double* wektor, int n)
{
	int i, j, k;
	cout.precision(4);        
	cout.setf(ios::fixed);

	double* x = new double[n];

	double** a = new double* [n];
	for (int i = 0; i < n; i++)
	{
		a[i] = new double[n + 1];
	}

	for (i = 0; i < n; i++)
	{
		for (j = 0; j < n; j++)
		{
			a[i][j] = macierz[i][j];
		}
	}

	for (i = 0; i < n; i++)
	{
		a[i][n] = wektor[i];
	}

	for (i = 0; i < n; i++)                   
		for (k = i + 1; k < n; k++)
			if (abs(a[i][i]) < abs(a[k][i]))
				for (j = 0; j <= n; j++)
				{
					double temp = a[i][j];
					a[i][j] = a[k][j];
					a[k][j] = temp;
				}

	for (i = 0; i < n - 1; i++)           
		for (k = i + 1; k < n; k++)
		{
			double t = a[k][i] / a[i][i];
			for (j = 0; j <= n; j++)
				a[k][j] = a[k][j] - t * a[i][j];   
		}


	for (i = n - 1; i >= 0; i--)                
	{                       
		x[i] = a[i][n];                
		for (j = i + 1; j < n; j++)
			if (j != i)                                             
				x[i] = x[i] - a[i][j] * x[j];
		x[i] = x[i] / a[i][i];           
	}
	
	double temp_min = x[0];
	double temp_max = x[0];
	
	

	for (int i = 0; i < (n - 1); i++)
	{
		if (temp_min > x[i + 1])
		{
			temp_min = x[i + 1];
		}
		if (temp_max < x[i + 1])
		{
			temp_max = x[i + 1];
		}
	}

	/*
	cout << "ROZWIAZANIE" << endl;
	for (int i = 0; i < n; i++)
	{
		cout << x[i] << endl;
	}
	*/

	cout << "TEMP MIN: " << temp_min << " " << "TEMP MAX: " << temp_max << endl;

	return x;
}
double** tworzenie_Hbc_P_C(Element element, int element_number, Global_data gd, Grid grid, int n, double** wektoryP_dla_elementow, double*** macierze_C_dla_elementow)
{
	//biorê element
	//przechodzê po jego œcianach, je¿eli dla danej œciany jest warunek brzegowy, to liczê jej macierz
	//lcizê jakobian dla tej œciany, obliczam ca³¹ tê macierz razem a alf¹ i jakobianem
	//sumujê te macierze dla scian, bo mog¹ byæ dwie, które w danym elemencie maj¹ warunek graniczny
	//otrzymujê macierz Hbc dla tego konkretnego elementu, teraz muszê j¹ jakoœ wypluæ z funkcji, zebraæ je w mainie wszystkie i wrzuciæ funkcji, która je zagreguje
	//zagregowan¹ macierz zsumowaæ z du¿¹, zagregowan¹ uprzednio macierz¹ H 

	double** tab_pomocnicze_dolna = new double* [n];
	double** tab_pomocnicze_prawa = new double* [n];
	double** tab_pomocnicze_gorna = new double* [n];
	double** tab_pomocnicze_lewa = new double* [n];

	for (int i = 0; i < n; i++)
	{
		tab_pomocnicze_dolna[i] = new double[2];
		tab_pomocnicze_prawa[i] = new double[2];
		tab_pomocnicze_gorna[i] = new double[2];
		tab_pomocnicze_lewa[i] = new double[2];
	}

	double* tab_dN_po_dksi_eta;
	tab_dN_po_dksi_eta = new double[n];

	if (n == 2)
	{
		tab_dN_po_dksi_eta[0] = -1.0 / sqrt(3);
		tab_dN_po_dksi_eta[1] = 1.0 / sqrt(3);
	}
	if (n == 3)
	{
		tab_dN_po_dksi_eta[0] = -sqrt(3.0 / 5.0);
		tab_dN_po_dksi_eta[1] = 0;
		tab_dN_po_dksi_eta[2] = sqrt(3.0 / 5.0);
	}
	if (n == 4)
	{
		tab_dN_po_dksi_eta[0] = -0.861136;
		tab_dN_po_dksi_eta[1] = -0.339981;
		tab_dN_po_dksi_eta[2] = 0.339981;
		tab_dN_po_dksi_eta[3] = 0.861136;
	}

	for (int i = 0; i < n; i++)
	{
		tab_pomocnicze_dolna[i][0] = tab_dN_po_dksi_eta[i];
		tab_pomocnicze_dolna[i][1] = -1;

		tab_pomocnicze_prawa[i][0] = 1;
		tab_pomocnicze_prawa[i][1] = tab_dN_po_dksi_eta[i];

		tab_pomocnicze_gorna[i][0] = tab_dN_po_dksi_eta[i];
		tab_pomocnicze_gorna[i][1] = 1;

		tab_pomocnicze_lewa[i][0] = -1;
		tab_pomocnicze_lewa[i][1] = tab_dN_po_dksi_eta[i];
	}

	double** tab_macierze_sciana_dolna = new double* [n];
	double** tab_macierze_sciana_prawa = new double* [n];
	double** tab_macierze_sciana_gorna = new double* [n];
	double** tab_macierze_sciana_lewa = new double* [n];

	for (int i = 0; i < n; i++)
	{
		tab_macierze_sciana_dolna[i] = new double[4];
		tab_macierze_sciana_prawa[i] = new double[4];
		tab_macierze_sciana_gorna[i] = new double[4];
		tab_macierze_sciana_lewa[i] = new double[4];
	}

	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			tab_macierze_sciana_dolna[i][j] = 0;
			tab_macierze_sciana_prawa[i][j] = 0;
			tab_macierze_sciana_gorna[i][j] = 0;
			tab_macierze_sciana_lewa[i][j] = 0;
		}
	}

	int tab_for_flags[4] = { 0,0,0,0 };
	for (int i = 0; i < 4; i++)
	{
		if (grid.ND[grid.EL[element_number].ID[i] - 1].BC == 1)
			tab_for_flags[i] = 1;
	}

	double l = 0.0;
	double a = 0.0;
	double l_ost = 0.0;
	double det_J = 0.0;


	//wagi
	double w[4];

	if (n == 2)
	{
		w[0] = 1.0;
		w[1] = 1.0;
	}
	if (n == 3)
	{
		w[0] = 5.0 / 9.0;
		w[1] = 8.0 / 9.0;
		w[2] = 5.0 / 9.0;
	}
	if (n == 4)
	{
		w[0] = 0.347855;
		w[1] = 0.652145;
		w[2] = 0.652145;
		w[3] = 0.347855;
	}

	double H_bc_g[4][4];
	double H_bc_l[4][4];
	double H_bc_d[4][4];
	double H_bc_p[4][4];

	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			H_bc_g[i][j] = 0;
			H_bc_l[i][j] = 0;
			H_bc_d[i][j] = 0;
			H_bc_p[i][j] = 0;
		}
	}

	double wektorP_g[4];
	double wektorP_l[4];
	double wektorP_d[4];
	double wektorP_p[4];

	for (int i = 0; i < 4; i++)
	{
		wektorP_g[i] = 0.0;
		wektorP_l[i] = 0.0;
		wektorP_d[i] = 0.0;
		wektorP_p[i] = 0.0;
	}

	//N1 = 1/4 * (1 + ksi) * (1 + eta)
	//N2 = 1/4 * (1 - ksi) * (1 + eta)
	//N3 = 1/4 * (1 - ksi) * (1 - eta)
	//N4 = 1/4 * (1 + ksi) * (1 - eta)

	if (tab_for_flags[0] == 1 && tab_for_flags[1] == 1)
	{
		for (int i = 0; i < n; i++)
		{
			tab_macierze_sciana_gorna[i][0] = 1.0 / 4.0 * (1.0 + tab_pomocnicze_gorna[i][0]) * (1.0 + tab_pomocnicze_gorna[i][1]);
			tab_macierze_sciana_gorna[i][1] = 1.0 / 4.0 * (1.0 - tab_pomocnicze_gorna[i][0]) * (1.0 + tab_pomocnicze_gorna[i][1]);
			tab_macierze_sciana_gorna[i][2] = 1.0 / 4.0 * (1.0 - tab_pomocnicze_gorna[i][0]) * (1.0 - tab_pomocnicze_gorna[i][1]);
			tab_macierze_sciana_gorna[i][3] = 1.0 / 4.0 * (1.0 + tab_pomocnicze_gorna[i][0]) * (1.0 - tab_pomocnicze_gorna[i][1]);
		}

		l = grid.ND[grid.EL[element_number].ID[0] - 1].x - grid.ND[grid.EL[element_number].ID[1] - 1].x;
		a = grid.ND[grid.EL[element_number].ID[0] - 1].y - grid.ND[grid.EL[element_number].ID[1] - 1].y;
		l_ost = sqrt(l * l + a * a);
		det_J = l_ost / 2;

		for (int i = 0; i < n; i++)
		{
			for (int j = 0; j < 4; j++)
			{
				for (int o = 0; o < 4; o++)
				{
					H_bc_g[j][o] += tab_macierze_sciana_gorna[i][j] * tab_macierze_sciana_gorna[i][o] * w[i];
				}
			}
		}

		for (int i = 0; i < n; i++)
		{
			for (int j = 0; j < 4; j++)
			{
				wektorP_g[j] += tab_macierze_sciana_gorna[i][j] * gd.tot * w[i];
			}
		}

		for (int i = 0; i < 4; i++)
		{
			for (int j = 0; j < 4; j++)
			{
				H_bc_g[i][j] = H_bc_g[i][j] * gd.alfa * det_J;
			}
		}

		for (int i = 0; i < 4; i++)
		{
			wektorP_g[i] = wektorP_g[i] * gd.alfa * det_J;
		}
	}
	if (tab_for_flags[1] == 1 && tab_for_flags[2] == 1)
	{
		for (int i = 0; i < n; i++)
		{
			tab_macierze_sciana_lewa[i][0] = 1.0 / 4.0 * (1.0 + tab_pomocnicze_lewa[i][0]) * (1.0 + tab_pomocnicze_lewa[i][1]);
			tab_macierze_sciana_lewa[i][1] = 1.0 / 4.0 * (1.0 - tab_pomocnicze_lewa[i][0]) * (1.0 + tab_pomocnicze_lewa[i][1]);
			tab_macierze_sciana_lewa[i][2] = 1.0 / 4.0 * (1.0 - tab_pomocnicze_lewa[i][0]) * (1.0 - tab_pomocnicze_lewa[i][1]);
			tab_macierze_sciana_lewa[i][3] = 1.0 / 4.0 * (1.0 + tab_pomocnicze_lewa[i][0]) * (1.0 - tab_pomocnicze_lewa[i][1]);
		}

		l = grid.ND[grid.EL[element_number].ID[2] - 1].y - grid.ND[grid.EL[element_number].ID[1] - 1].y;
		a = grid.ND[grid.EL[element_number].ID[2] - 1].x - grid.ND[grid.EL[element_number].ID[1] - 1].x;
		l_ost = sqrt(l * l + a * a);
		det_J = l_ost / 2;

		cout << "JAKOBIAN DLA SCIANY LEWEJ: " << det_J << endl;

		for (int i = 0; i < n; i++)
		{
			for (int j = 0; j < 4; j++)
			{
				for (int o = 0; o < 4; o++)
				{
					H_bc_l[j][o] += tab_macierze_sciana_lewa[i][j] * tab_macierze_sciana_lewa[i][o] * w[i];
				}
			}
		}

		for (int i = 0; i < n; i++)
		{
			for (int j = 0; j < 4; j++)
			{
				wektorP_l[j] += tab_macierze_sciana_lewa[i][j] * gd.tot * w[i];
			}
		}

		for (int i = 0; i < 4; i++)
		{
			wektorP_l[i] = wektorP_l[i] * gd.alfa * det_J;
		}

		for (int i = 0; i < 4; i++)
		{
			for (int j = 0; j < 4; j++)
			{
				H_bc_l[i][j] = H_bc_l[i][j] * gd.alfa * det_J;
			}
		}
	}
	if (tab_for_flags[2] == 1 && tab_for_flags[3] == 1)
	{
		for (int i = 0; i < n; i++)
		{
			tab_macierze_sciana_dolna[i][0] = 1.0 / 4.0 * (1.0 + tab_pomocnicze_dolna[i][0]) * (1.0 + tab_pomocnicze_dolna[i][1]);
			tab_macierze_sciana_dolna[i][1] = 1.0 / 4.0 * (1.0 - tab_pomocnicze_dolna[i][0]) * (1.0 + tab_pomocnicze_dolna[i][1]);
			tab_macierze_sciana_dolna[i][2] = 1.0 / 4.0 * (1.0 - tab_pomocnicze_dolna[i][0]) * (1.0 - tab_pomocnicze_dolna[i][1]);
			tab_macierze_sciana_dolna[i][3] = 1.0 / 4.0 * (1.0 + tab_pomocnicze_dolna[i][0]) * (1.0 - tab_pomocnicze_dolna[i][1]);
		}

		l = grid.ND[grid.EL[element_number].ID[3] - 1].x - grid.ND[grid.EL[element_number].ID[2] - 1].x;
		a = grid.ND[grid.EL[element_number].ID[3] - 1].y - grid.ND[grid.EL[element_number].ID[2] - 1].y;
		l_ost = sqrt(l * l + a * a);
		det_J = l_ost / 2;

		cout << "JAKOBIAN DLA SCIANY DOLNEJ: " << det_J << endl;

		for (int i = 0; i < n; i++)
		{
			for (int j = 0; j < 4; j++)
			{
				for (int o = 0; o < 4; o++)
				{
					H_bc_d[j][o] += tab_macierze_sciana_dolna[i][j] * tab_macierze_sciana_dolna[i][o] * w[i];
				}
			}
		}

		for (int i = 0; i < n; i++)
		{
			for (int j = 0; j < 4; j++)
			{
				wektorP_d[j] += tab_macierze_sciana_dolna[i][j] * gd.tot * w[i];
			}
		}

		for (int i = 0; i < 4; i++)
		{
			for (int j = 0; j < 4; j++)
			{
				H_bc_d[i][j] = H_bc_d[i][j] * gd.alfa * det_J;
			}
		}

		for (int i = 0; i < 4; i++)
		{
			wektorP_d[i] = wektorP_d[i] * gd.alfa * det_J;
		}
	}
	if (tab_for_flags[3] == 1 && tab_for_flags[0] == 1)
	{
		for (int i = 0; i < n; i++)
		{
			tab_macierze_sciana_prawa[i][0] = 1.0 / 4.0 * (1.0 + tab_pomocnicze_prawa[i][0]) * (1.0 + tab_pomocnicze_prawa[i][1]);
			tab_macierze_sciana_prawa[i][1] = 1.0 / 4.0 * (1.0 - tab_pomocnicze_prawa[i][0]) * (1.0 + tab_pomocnicze_prawa[i][1]);
			tab_macierze_sciana_prawa[i][2] = 1.0 / 4.0 * (1.0 - tab_pomocnicze_prawa[i][0]) * (1.0 - tab_pomocnicze_prawa[i][1]);
			tab_macierze_sciana_prawa[i][3] = 1.0 / 4.0 * (1.0 + tab_pomocnicze_prawa[i][0]) * (1.0 - tab_pomocnicze_prawa[i][1]);
		}

		l = grid.ND[grid.EL[element_number].ID[0] - 1].y - grid.ND[grid.EL[element_number].ID[3] - 1].y;
		a = grid.ND[grid.EL[element_number].ID[0] - 1].x - grid.ND[grid.EL[element_number].ID[3] - 1].x;
		l_ost = sqrt(l * l + a * a);
		det_J = l_ost / 2;

		cout << "JAKOBIAN DLA SCIANY PRAWEJ: " << det_J << endl;

		for (int i = 0; i < n; i++)
		{
			for (int j = 0; j < 4; j++)
			{
				for (int o = 0; o < 4; o++)
				{
					H_bc_p[j][o] += tab_macierze_sciana_prawa[i][j] * tab_macierze_sciana_prawa[i][o] * w[i];
				}
			}
		}

		for (int i = 0; i < n; i++)
		{
			for (int j = 0; j < 4; j++)
			{
				wektorP_p[j] += tab_macierze_sciana_prawa[i][j] * gd.tot * w[i];
			}
		}

		for (int i = 0; i < 4; i++)
		{
			for (int j = 0; j < 4; j++)
			{
				H_bc_p[i][j] = H_bc_p[i][j] * gd.alfa * det_J;
			}
		}

		for (int i = 0; i < 4; i++)
		{
			wektorP_p[i] = wektorP_p[i] * gd.alfa * det_J;
		}

		//wypisywanie œcian
		/*
		cout << "SCIANA GORNA" << endl;
		for (int i = 0; i < 4; i++)
		{
			for (int j = 0; j < 4; j++)
			{
				cout << H_bc_g[i][j] << "	";
			}
			cout << endl;
		}
		cout << endl;
		cout << endl;
		cout << "SCIANA DOLNA" << endl;
		for (int i = 0; i < 4; i++)
		{
			for (int j = 0; j < 4; j++)
			{
				cout << H_bc_d[i][j] << "	";
			}
			cout << endl;
		}
		cout << endl;
		cout << endl;
		cout << "SCIANA LEWA" << endl;
		for (int i = 0; i < 4; i++)
		{
			for (int j = 0; j < 4; j++)
			{
				cout << H_bc_l[i][j] << "	";
			}
			cout << endl;
		}
		cout << endl;
		cout << endl;
		cout << "SCIANA PRAWA" << endl;
		for (int i = 0; i < 4; i++)
		{
			for (int j = 0; j < 4; j++)
			{
				cout << H_bc_p[i][j] << "	";
			}
			cout << endl;
		}
		cout << endl;
		cout << endl;
		*/
	}

	double** H_bc;
	H_bc = new double* [4];

	H_bc[0] = new double[4];
	H_bc[1] = new double[4];
	H_bc[2] = new double[4];
	H_bc[3] = new double[4];

	double* wektorP_finalny = new double[4];

	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			H_bc[i][j] = 0;
		}
	}

	for (int i = 0; i < 4; i++)
	{
		wektorP_finalny[i] = 0;
	}

	for (int i = 0; i < 4; i++)
	{
		wektorP_finalny[i] = wektorP_g[i] + wektorP_l[i] + wektorP_d[i] + wektorP_p[i];
	}

	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			H_bc[i][j] = H_bc_g[i][j] + H_bc_l[i][j] + H_bc_d[i][j] + H_bc_p[i][j];
		}
	}

	cout << "WEKTOR P DLA ELEMENTU " << element_number << endl;
	for (int i = 0; i < 4; i++)
	{
		cout << wektorP_finalny[i] << "	";
	}
	cout << endl << endl;

	cout << "MACIERZ H_bc DLA ELEMENTU " << element_number << endl;
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			cout << H_bc[i][j] << "	";
		}
		cout << endl;
	}
	cout << endl << endl;

	wektoryP_dla_elementow[element_number] = wektorP_finalny;

	double* tab_dN_po_dksi_C;
	tab_dN_po_dksi_C = new double[n * n];

	double* tab_dN_po_deta_C;
	tab_dN_po_deta_C = new double[n * n];

	double* wagi_ksi;
	wagi_ksi = new double[n * n];

	double* wagi_eta;
	wagi_eta = new double[n * n];

	if (n == 2)
	{
		tab_dN_po_dksi_C[0] = -1.0 / sqrt(3);
		tab_dN_po_dksi_C[1] = 1.0 / sqrt(3);
		tab_dN_po_dksi_C[2] = -1.0 / sqrt(3);
		tab_dN_po_dksi_C[3] = 1.0 / sqrt(3);

		tab_dN_po_deta_C[0] = -1.0 / sqrt(3);
		tab_dN_po_deta_C[1] = -1.0 / sqrt(3);
		tab_dN_po_deta_C[2] = 1.0 / sqrt(3);
		tab_dN_po_deta_C[3] = 1.0 / sqrt(3);

		wagi_ksi[0] = 1.0;
		wagi_ksi[1] = 1.0;
		wagi_ksi[2] = 1.0;
		wagi_ksi[3] = 1.0;

		wagi_eta[0] = 1.0;
		wagi_eta[1] = 1.0;
		wagi_eta[2] = 1.0;
		wagi_eta[3] = 1.0;
	}
	if (n == 3)
	{
		tab_dN_po_dksi_C[0] = -sqrt(3.0 / 5.0);
		tab_dN_po_dksi_C[1] = 0;
		tab_dN_po_dksi_C[2] = sqrt(3.0 / 5.0);
		tab_dN_po_dksi_C[3] = -sqrt(3.0 / 5.0);
		tab_dN_po_dksi_C[4] = 0;
		tab_dN_po_dksi_C[5] = sqrt(3.0 / 5.0);
		tab_dN_po_dksi_C[6] = -sqrt(3.0 / 5.0);
		tab_dN_po_dksi_C[7] = 0;
		tab_dN_po_dksi_C[8] = sqrt(3.0 / 5.0);

		tab_dN_po_deta_C[0] = -sqrt(3.0 / 5.0);
		tab_dN_po_deta_C[1] = -sqrt(3.0 / 5.0);
		tab_dN_po_deta_C[2] = -sqrt(3.0 / 5.0);
		tab_dN_po_deta_C[3] = 0;
		tab_dN_po_deta_C[4] = 0;
		tab_dN_po_deta_C[5] = 0;
		tab_dN_po_deta_C[6] = sqrt(3.0 / 5.0);
		tab_dN_po_deta_C[7] = sqrt(3.0 / 5.0);
		tab_dN_po_deta_C[8] = sqrt(3.0 / 5.0);

		wagi_ksi[0] = 5.0 / 9.0;
		wagi_ksi[1] = 8.0 / 9.0;
		wagi_ksi[2] = 5.0 / 9.0;
		wagi_ksi[3] = 5.0 / 9.0;
		wagi_ksi[4] = 8.0 / 9.0;
		wagi_ksi[5] = 5.0 / 9.0;
		wagi_ksi[6] = 5.0 / 9.0;
		wagi_ksi[7] = 8.0 / 9.0;
		wagi_ksi[8] = 5.0 / 9.0;

		wagi_eta[0] = 5.0 / 9.0;
		wagi_eta[1] = 5.0 / 9.0;
		wagi_eta[2] = 5.0 / 9.0;
		wagi_eta[3] = 8.0 / 9.0;
		wagi_eta[4] = 8.0 / 9.0;
		wagi_eta[5] = 8.0 / 9.0;
		wagi_eta[6] = 5.0 / 9.0;
		wagi_eta[7] = 5.0 / 9.0;
		wagi_eta[8] = 5.0 / 9.0;
	}
	if (n == 4)
	{
		tab_dN_po_dksi_C[0] = -0.861136;
		tab_dN_po_dksi_C[1] = -0.339981;
		tab_dN_po_dksi_C[2] = 0.339981;
		tab_dN_po_dksi_C[3] = 0.861136;
		tab_dN_po_dksi_C[4] = -0.861136;
		tab_dN_po_dksi_C[5] = -0.339981;
		tab_dN_po_dksi_C[6] = 0.339981;
		tab_dN_po_dksi_C[7] = 0.861136;
		tab_dN_po_dksi_C[8] = -0.861136;
		tab_dN_po_dksi_C[9] = -0.339981;
		tab_dN_po_dksi_C[10] = 0.339981;
		tab_dN_po_dksi_C[11] = 0.861136;
		tab_dN_po_dksi_C[12] = -0.861136;
		tab_dN_po_dksi_C[13] = -0.339981;
		tab_dN_po_dksi_C[14] = 0.339981;
		tab_dN_po_dksi_C[15] = 0.861136;

		tab_dN_po_deta_C[0] = -0.861136;
		tab_dN_po_deta_C[1] = -0.861136;
		tab_dN_po_deta_C[2] = -0.861136;
		tab_dN_po_deta_C[3] = -0.861136;
		tab_dN_po_deta_C[4] = -0.339981;
		tab_dN_po_deta_C[5] = -0.339981;
		tab_dN_po_deta_C[6] = -0.339981;
		tab_dN_po_deta_C[7] = -0.339981;
		tab_dN_po_deta_C[8] = 0.339981;
		tab_dN_po_deta_C[9] = 0.339981;
		tab_dN_po_deta_C[10] = 0.339981;
		tab_dN_po_deta_C[11] = 0.339981;
		tab_dN_po_deta_C[12] = 0.861136;
		tab_dN_po_deta_C[13] = 0.861136;
		tab_dN_po_deta_C[14] = 0.861136;
		tab_dN_po_deta_C[15] = 0.861136;

		wagi_ksi[0] = 0.347855;
		wagi_ksi[1] = 0.652145;
		wagi_ksi[2] = 0.652145;
		wagi_ksi[3] = 0.347855;
		wagi_ksi[4] = 0.347855;
		wagi_ksi[5] = 0.652145;
		wagi_ksi[6] = 0.652145;
		wagi_ksi[7] = 0.347855;
		wagi_ksi[8] = 0.347855;
		wagi_ksi[9] = 0.652145;;
		wagi_ksi[10] = 0.652145;
		wagi_ksi[11] = 0.347855;
		wagi_ksi[12] = 0.347855;
		wagi_ksi[13] = 0.652145;
		wagi_ksi[14] = 0.652145;
		wagi_ksi[15] = 0.347855;

		wagi_eta[0] = 0.347855;
		wagi_eta[1] = 0.347855;
		wagi_eta[2] = 0.347855;
		wagi_eta[3] = 0.347855;
		wagi_eta[4] = 0.652145;
		wagi_eta[5] = 0.652145;
		wagi_eta[6] = 0.652145;
		wagi_eta[7] = 0.652145;
		wagi_eta[8] = 0.652145;
		wagi_eta[9] = 0.652145;
		wagi_eta[10] = 0.652145;
		wagi_eta[11] = 0.652145;
		wagi_eta[12] = 0.347855;
		wagi_eta[13] = 0.347855;
		wagi_eta[14] = 0.347855;
		wagi_eta[15] = 0.347855;
	}

	double** tab_N_dla_C;
	tab_N_dla_C = new double* [n * n];

	for (int i = 0; i < n * n; i++)
	{
		tab_N_dla_C[i] = new double[4];
	}

	for (int i = 0; i < n * n; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			tab_N_dla_C[i][j] = 0;
		}
	}

	for (int i = 0; i < n * n; i++)
	{
		tab_N_dla_C[i][0] = 1.0 / 4.0 * (1.0 - tab_dN_po_dksi_C[i]) * (1.0 - tab_dN_po_deta_C[i]);
		tab_N_dla_C[i][1] = 1.0 / 4.0 * (1.0 + tab_dN_po_dksi_C[i]) * (1.0 - tab_dN_po_deta_C[i]);
		tab_N_dla_C[i][2] = 1.0 / 4.0 * (1.0 + tab_dN_po_dksi_C[i]) * (1.0 + tab_dN_po_deta_C[i]);
		tab_N_dla_C[i][3] = 1.0 / 4.0 * (1.0 - tab_dN_po_dksi_C[i]) * (1.0 + tab_dN_po_deta_C[i]);
	}

	double x00;
	double y01;
	double x10;
	double y11;

	double x1 = 0.0;
	double y1 = 0.0;

	double x2 = 0.0;
	double y2 = 0.0;

	double x3 = 0.0;
	double y3 = 0.0;

	double x4 = 0.0;
	double y4 = 0.0;

	x1 = grid.ND[grid.EL[element_number].ID[0] - 1].x;
	y1 = grid.ND[grid.EL[element_number].ID[0] - 1].y;

	x2 = grid.ND[grid.EL[element_number].ID[1] - 1].x;
	y2 = grid.ND[grid.EL[element_number].ID[1] - 1].y;

	x3 = grid.ND[grid.EL[element_number].ID[2] - 1].x;
	y3 = grid.ND[grid.EL[element_number].ID[2] - 1].y;

	x4 = grid.ND[grid.EL[element_number].ID[3] - 1].x;
	y4 = grid.ND[grid.EL[element_number].ID[3] - 1].y;

	//cout << x00 << "	" << y01 << endl << x10 << "	" << y11 << endl << endl;

	double jakobian = 0;

	double** macierz_C_dla_elementu;

	macierz_C_dla_elementu = new double* [4];

	macierz_C_dla_elementu[0] = new double[4];
	macierz_C_dla_elementu[1] = new double[4];
	macierz_C_dla_elementu[2] = new double[4];
	macierz_C_dla_elementu[3] = new double[4];

	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			macierz_C_dla_elementu[i][j] = 0;
		}
	}


	for (int u = 0; u < n * n; u++)
	{
		x00 = element.tab_dN_po_dksi[u][0] * x1 + element.tab_dN_po_dksi[u][1] * x2 + element.tab_dN_po_dksi[u][2] * x3 + element.tab_dN_po_dksi[u][3] * x4;
		y01 = element.tab_dN_po_dksi[u][0] * y1 + element.tab_dN_po_dksi[u][1] * y2 + element.tab_dN_po_dksi[u][2] * y3 + element.tab_dN_po_dksi[u][3] * y4;
		x10 = element.tab_dN_po_deta[u][0] * x1 + element.tab_dN_po_deta[u][1] * x2 + element.tab_dN_po_deta[u][2] * x3 + element.tab_dN_po_deta[u][3] * x4;
		y11 = element.tab_dN_po_deta[u][0] * y1 + element.tab_dN_po_deta[u][1] * y2 + element.tab_dN_po_deta[u][2] * y3 + element.tab_dN_po_deta[u][3] * y4;

		jakobian = x00 * y11 - (y01 * x10);

		for (int i = 0; i < 4; i++)
		{
			for (int j = 0; j < 4; j++)
			{
				macierz_C_dla_elementu[i][j] += tab_N_dla_C[u][i] * tab_N_dla_C[u][j] * jakobian * wagi_ksi[u] * wagi_eta[u] * gd.specificHeat * gd.density;
			}
		}
	}

	cout << "MACIERZ C DLA ELEMENTU " << element_number << endl;
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			cout << macierz_C_dla_elementu[i][j] << "	";
		}
		cout << endl;
	}
	cout << endl << endl;

	macierze_C_dla_elementow[element_number] = macierz_C_dla_elementu;

	return H_bc;
}
double** tworzenie_macierzy_H(double x1, double y1, double x2, double y2, double x3, double y3, double x4, double y4, Element element, int point_number, Global_data gd, Grid grid, int n)
{
	double x00;
	double y01;
	double x10;
	double y11;

	x00 = element.tab_dN_po_dksi[point_number][0] * x1 + element.tab_dN_po_dksi[point_number][1] * x2 + element.tab_dN_po_dksi[point_number][2] * x3 + element.tab_dN_po_dksi[point_number][3] * x4;
	y01 = element.tab_dN_po_dksi[point_number][0] * y1 + element.tab_dN_po_dksi[point_number][1] * y2 + element.tab_dN_po_dksi[point_number][2] * y3 + element.tab_dN_po_dksi[point_number][3] * y4;
	x10 = element.tab_dN_po_deta[point_number][0] * x1 + element.tab_dN_po_deta[point_number][1] * x2 + element.tab_dN_po_deta[point_number][2] * x3 + element.tab_dN_po_deta[point_number][3] * x4;
	y11 = element.tab_dN_po_deta[point_number][0] * y1 + element.tab_dN_po_deta[point_number][1] * y2 + element.tab_dN_po_deta[point_number][2] * y3 + element.tab_dN_po_deta[point_number][3] * y4;

	//cout << x00 << "	" << y01 << endl << x10 << "	" << y11 << endl << endl;

	double jakobian = x00 * y11 - (y01 * x10);

	//potem do obliczeñ trzeba zmieniæ macierz x00, y01, na y11, -y01
	//										   x10, y11     -x10, x00

	double tab_dN_po_dx[4];
	double tab_dN_po_dy[4];

	for (int o = 0; o < 4; o++)
	{
		tab_dN_po_dx[o] = (1 / jakobian * y11 * element.tab_dN_po_dksi[point_number][o]) + (1 / jakobian * (-y01) * element.tab_dN_po_deta[point_number][o]);
		tab_dN_po_dy[o] = (1 / jakobian * (-x10) * element.tab_dN_po_dksi[point_number][o]) + (1 / jakobian * x00 * element.tab_dN_po_deta[point_number][o]);
	}

	//cout << "[" << tab_dN_po_dx[0] << ", " << tab_dN_po_dx[1] << ", " << tab_dN_po_dx[2] << ", " << tab_dN_po_dx[3] << "]" << "		" << "[" << tab_dN_po_dy[0] << ", " << tab_dN_po_dy[1] << ", " << tab_dN_po_dy[2] << ", " << tab_dN_po_dy[3] << "]" << endl << endl;

	double** tab_nowa_macierz_x;
	tab_nowa_macierz_x = new double* [4];

	tab_nowa_macierz_x[0] = new double[4];
	tab_nowa_macierz_x[1] = new double[4];
	tab_nowa_macierz_x[2] = new double[4];
	tab_nowa_macierz_x[3] = new double[4];

	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			tab_nowa_macierz_x[i][j] = tab_dN_po_dx[i] * tab_dN_po_dx[j];
		}
	}

	double** tab_nowa_macierz_y;
	tab_nowa_macierz_y = new double* [4];

	tab_nowa_macierz_y[0] = new double[4];
	tab_nowa_macierz_y[1] = new double[4];
	tab_nowa_macierz_y[2] = new double[4];
	tab_nowa_macierz_y[3] = new double[4];

	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			tab_nowa_macierz_y[i][j] = tab_dN_po_dy[i] * tab_dN_po_dy[j];
		}
	}

	double** tab_macierz_h_dla_punktu;
	tab_macierz_h_dla_punktu = new double* [4];

	tab_macierz_h_dla_punktu[0] = new double[4];
	tab_macierz_h_dla_punktu[1] = new double[4];
	tab_macierz_h_dla_punktu[2] = new double[4];
	tab_macierz_h_dla_punktu[3] = new double[4];

	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			tab_macierz_h_dla_punktu[i][j] = (tab_nowa_macierz_x[i][j] + tab_nowa_macierz_y[i][j]) * gd.conductivity * jakobian;
		}
	}
	/*
	cout << "MACIERZ H:" << endl;

	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			cout << tab_macierz_h_dla_punktu[i][j] << "		";
		}
		cout << endl;
	}

	cout << "--------------------------------------------------------------------------" << endl;
	*/

	return tab_macierz_h_dla_punktu;
}
bool isNumber(const string& str)
{
	for (char const& c : str)
	{
		if (isdigit(c) == 0) return false;
	}
	return true;
}
void printAll(Global_data g, Grid gr)
{
	cout << "simulationTime = " << g.simulationTime << endl;
	cout << "simulationStepTime = " << g.simulationStepTime << endl;
	cout << "conductivity = " << g.conductivity << endl;
	cout << "alfa = " << g.alfa << endl;
	cout << "tot = " << g.tot << endl;
	cout << "initialTemp = " << g.initialTemp << endl;
	cout << "density = " << g.density << endl;
	cout << "specificHeat = " << g.specificHeat << endl;
	cout << "nodesNumber = " << g.nodesNumber << endl;
	cout << "elementsNumber = " << g.elementsNumber << endl << endl;

	for (int i = 0; i < g.nodesNumber; i++)
	{
		cout << "Node " << i + 1 << " | x=" << gr.ND[i].x << " | y=" << gr.ND[i].y << " | BC=" << gr.ND[i].BC << endl;
	}

	cout << endl;

	for (int i = 0; i < g.elementsNumber; i++)
	{
		cout << "Element " << i + 1 << " | Nodes: " << gr.EL[i].ID[0] << ", " << gr.EL[i].ID[1] << ", " << gr.EL[i].ID[2] << ", " << gr.EL[i].ID[3] << endl;
	}
}
int readToStruct(Global_data* gd, Grid* gr)
{
	ifstream myfile("siatka.txt");

	vector<int> values;

	string input;

	while (myfile >> input)
	{
		if (input == "*Node")
			break;

		if (isNumber(input))
			values.push_back(stoi(input));
	}

	Global_data Gdata;

	Gdata.simulationTime = values[0];
	Gdata.simulationStepTime = values[1];
	Gdata.conductivity = values[2];
	Gdata.alfa = values[3];
	Gdata.tot = values[4];
	Gdata.initialTemp = values[5];
	Gdata.density = values[6];
	Gdata.specificHeat = values[7];
	Gdata.nodesNumber = values[8];
	Gdata.elementsNumber = values[9];

	Grid grid1;

	int counter = 0;
	double tempX, tempY = 0;
	int actual = 0;

	for (int i = 0; i < Gdata.nodesNumber * 3; i++)
	{
		counter++;
		myfile >> input;

		if (counter % 3 == 2)
		{
			tempX = stod(input);
		}

		if (counter % 3 == 0)
		{
			tempY = stod(input);

			grid1.ND.push_back(Node());

			grid1.ND[actual].x = tempX;
			grid1.ND[actual].y = tempY;

			actual++;
		}
	}

	myfile >> input;
	myfile >> input;

	counter = 0;
	actual = 0;

	int tempID[4] = { 0 };

	for (int i = 0; i < Gdata.elementsNumber * 5; i++)
	{
		counter++;
		myfile >> input;

		if (counter % 5 == 2)
		{
			tempID[0] = stoi(input, 0, 10);
		}

		if (counter % 5 == 3)
		{
			tempID[1] = stoi(input, 0, 10);
		}

		if (counter % 5 == 4)
		{
			tempID[2] = stoi(input, 0, 10);
		}

		if (counter % 5 == 0)
		{
			tempID[3] = stoi(input, 0, 10);

			grid1.EL.push_back(ElementID());

			grid1.EL[actual].ID[0] = tempID[0];
			grid1.EL[actual].ID[1] = tempID[1];
			grid1.EL[actual].ID[2] = tempID[2];
			grid1.EL[actual].ID[3] = tempID[3];

			actual++;
		}
	}

	myfile >> input;

	while (myfile >> input)
	{
		grid1.ND[stoi(input, 0, 10) - 1].BC = 1;
	}

	myfile.close();

	printAll(Gdata, grid1);

	*gd = Gdata;
	*gr = grid1;

	return 0;
}

int main()
{
	Grid grid;
	Global_data globalData;

	readToStruct(&globalData, &grid);

	//zmienna odpowiedzialna za to ilu punktowy ma byæ schemat ca³kowania
	int n = 3;

	double x1 = 0.0;
	double y1 = 0.0;

	double x2 = 0.0;
	double y2 = 0.0;

	double x3 = 0.0;
	double y3 = 0.0;

	double x4 = 0.0;
	double y4 = 0.0;

	double w[4];

	Element element1(n);

	double** tab_macierz_H_dla_punktow[16];
	double*** tab_macierze_H_dla_elementow = new double** [globalData.elementsNumber];

	for (int u = 0; u < globalData.elementsNumber; u++)
	{
		x1 = grid.ND[grid.EL[u].ID[0] - 1].x;
		y1 = grid.ND[grid.EL[u].ID[0] - 1].y;

		x2 = grid.ND[grid.EL[u].ID[1] - 1].x;
		y2 = grid.ND[grid.EL[u].ID[1] - 1].y;

		x3 = grid.ND[grid.EL[u].ID[2] - 1].x;
		y3 = grid.ND[grid.EL[u].ID[2] - 1].y;

		x4 = grid.ND[grid.EL[u].ID[3] - 1].x;
		y4 = grid.ND[grid.EL[u].ID[3] - 1].y;

		for (int i = 0; i < n * n; i++)
		{
			tab_macierz_H_dla_punktow[i] = tworzenie_macierzy_H(x1, y1, x2, y2, x3, y3, x4, y4, element1, i, globalData, grid, n);
		}

		double** tab_macierz_H_finalna;
		tab_macierz_H_finalna = new double* [4];

		tab_macierz_H_finalna[0] = new double[4];
		tab_macierz_H_finalna[1] = new double[4];
		tab_macierz_H_finalna[2] = new double[4];
		tab_macierz_H_finalna[3] = new double[4];

		//WAGI I LICZENIE MACIERZY H
		if (n == 2)
		{
			w[0] = 1.0;
			w[1] = 1.0;

			for (int i = 0; i < 4; i++)
			{
				for (int j = 0; j < 4; j++)
				{
					tab_macierz_H_finalna[i][j] = tab_macierz_H_dla_punktow[0][i][j] * w[0] * w[0] + tab_macierz_H_dla_punktow[1][i][j] * w[1] * w[0] + tab_macierz_H_dla_punktow[2][i][j] * w[0] * w[1] + tab_macierz_H_dla_punktow[3][i][j] * w[1] * w[1];
				}
			}
		}
		if (n == 3)
		{
			w[0] = 5.0 / 9.0;
			w[1] = 8.0 / 9.0;
			w[2] = 5.0 / 9.0;

			for (int i = 0; i < 4; i++)
			{
				for (int j = 0; j < 4; j++)
				{
					tab_macierz_H_finalna[i][j] = tab_macierz_H_dla_punktow[0][i][j] * w[0] * w[0] + tab_macierz_H_dla_punktow[1][i][j] * w[1] * w[0] + tab_macierz_H_dla_punktow[2][i][j] * w[2] * w[0] + tab_macierz_H_dla_punktow[3][i][j] * w[0] * w[1] + tab_macierz_H_dla_punktow[4][i][j] * w[1] * w[1] + tab_macierz_H_dla_punktow[5][i][j] * w[2] * w[1] + tab_macierz_H_dla_punktow[6][i][j] * w[0] * w[2] + tab_macierz_H_dla_punktow[7][i][j] * w[1] * w[2] + tab_macierz_H_dla_punktow[8][i][j] * w[2] * w[2];
				}
			}
		}
		if (n == 4)
		{
			w[0] = 0.347855;
			w[1] = 0.652145;
			w[2] = 0.652145;
			w[3] = 0.347855;

			for (int i = 0; i < 4; i++)
			{
				for (int j = 0; j < 4; j++)
				{
					tab_macierz_H_finalna[i][j] = tab_macierz_H_dla_punktow[0][i][j] * w[0] * w[0] + tab_macierz_H_dla_punktow[1][i][j] * w[1] * w[0] + tab_macierz_H_dla_punktow[2][i][j] * w[2] * w[0] + tab_macierz_H_dla_punktow[3][i][j] * w[3] * w[0] + tab_macierz_H_dla_punktow[4][i][j] * w[0] * w[1] + tab_macierz_H_dla_punktow[5][i][j] * w[1] * w[1] + tab_macierz_H_dla_punktow[6][i][j] * w[2] * w[1] + tab_macierz_H_dla_punktow[7][i][j] * w[3] * w[1] + tab_macierz_H_dla_punktow[8][i][j] * w[0] * w[2] + tab_macierz_H_dla_punktow[9][i][j] * w[1] * w[2] + tab_macierz_H_dla_punktow[10][i][j] * w[2] * w[2] + tab_macierz_H_dla_punktow[11][i][j] * w[3] * w[2] + tab_macierz_H_dla_punktow[12][i][j] * w[0] * w[3] + tab_macierz_H_dla_punktow[13][i][j] * w[1] * w[3] + tab_macierz_H_dla_punktow[14][i][j] * w[2] * w[3] + tab_macierz_H_dla_punktow[15][i][j] * w[3] * w[3];
				}
			}
		}

		tab_macierze_H_dla_elementow[u] = tab_macierz_H_finalna;
	}

	double*** tab_macierze_Hbc_dla_elementow = new double** [globalData.elementsNumber];
	double*** tab_macierze_C_dla_elementow = new double** [globalData.elementsNumber];

	double** wektory_P_dla_elementow = new double* [globalData.elementsNumber];

	for (int i = 0; i < globalData.elementsNumber; i++)
	{
		tab_macierze_Hbc_dla_elementow[i] = tworzenie_Hbc_P_C(element1, i, globalData, grid, n, wektory_P_dla_elementow, tab_macierze_C_dla_elementow);
	}

	cout << "WEKTORY" << endl;
	for (int i = 0; i < globalData.elementsNumber; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			cout << wektory_P_dla_elementow[i][j] << "	";
		}
		cout << endl;
	}
	cout << endl << endl;


	Uklad uklad1(tab_macierze_H_dla_elementow, tab_macierze_Hbc_dla_elementow, grid, wektory_P_dla_elementow, globalData, tab_macierze_C_dla_elementow);

	obliczanie_rownania(uklad1.macierz_suma, uklad1.wektorP_duzy, globalData.nodesNumber);





	//tworzê macierz [H] + [C]/deltaT

	for (int i = 0; i < globalData.nodesNumber; i++)
	{
		for (int j = 0; j < globalData.nodesNumber; j++)
		{
			uklad1.macierz_H_plus_C_przez_deltaT[i][j] = uklad1.macierz_suma[i][j] + (uklad1.macierz_C[i][j] / globalData.simulationStepTime);
		}
	}

	for (int u = globalData.simulationStepTime; u <= globalData.simulationTime; )
	{
		//[C]/deltaT * {t0} + {P}

		for (int i = 0; i < globalData.nodesNumber; i++)
		{
			for (int j = 0; j < globalData.nodesNumber; j++)
			{
				uklad1.wektor_do_ostatniego_rownania[i] += (uklad1.macierz_C[i][j] / globalData.simulationStepTime) * uklad1.wektor_t0[j];
			}
		}

		for (int i = 0; i < globalData.nodesNumber; i++)
		{
			uklad1.wektor_do_ostatniego_rownania[i] += uklad1.wektorP_duzy[i];
		}

		cout << "TIME: " << u << " ";
		uklad1.wektor_t0 = obliczanie_rownania(uklad1.macierz_H_plus_C_przez_deltaT, uklad1.wektor_do_ostatniego_rownania, globalData.nodesNumber);

		for (int i = 0; i < globalData.nodesNumber; i++)
		{
			uklad1.wektor_do_ostatniego_rownania[i] = 0;
		}
		u += globalData.simulationStepTime;
	}
	
	return 0;
}

