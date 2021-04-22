
#include<iostream>
#include <omp.h> 

using namespace std;

void make_matrix(int** s, int n, int m) 
{
	for (int i = 0; i < n; i++) 
		s[i] = new int[m];
	for (int i = 0; i < n; i++)
		for (int j = 0; j < m; j++)
			s[i][j] = 1;
}

void print_matrix(int** x, int z, int c) {
	for (int i = 0; i < z; i++) 
	{
		for (int j = 0; j < c; j++)
			cout << x[i][j] << " ";
		cout << endl;
	}
}

void pr_matrix_posled(int** a1, int** b1,int**c1, int n1, int m1, int m2) {
	for (int i = 0; i < n1; i++) {
		for (int j = 0; j < m2; j++) {
			c1[i][j] = 0;
			for (int l = 0; l < m1; l++)
				c1[i][j] = c1[i][j] + (a1[i][l] * b1[l][j]);
		}
			
	}
}

int sum_elem_matrix(int** a2, int n, int m) {
	long int sum = 0;
	for (int i = 0; i < n; i++) 
		for (int j = 0; j < m; j++) 
			sum += a2[i][j];	
	return sum;
}


void main()
{
	setlocale(LC_ALL, "rus");
	int  n1, m1;
	int  n2, m2;
	

	cout << "Введите размерность первой матрицы: ";
	cin >> n1 >> m1;
	cout << "Введите размерность второй матрицы: ";
	cin >> n2 >> m2;
	if (m1 != n2)
		cout << " Умножение невозможно";
	else {

		int** a = new int* [n1];
		make_matrix(a, n1, m1);
	

		int** b = new int* [n2];
		make_matrix(b, n2, m2);
		

		int** c = new int* [n1];
		make_matrix(c, n1, m2);
		for (int i = 0; i < n1; i++)
			for (int j = 0; j < m2; j++)
				c[i][j] = 0;

		double timein, timeout;
		timein = omp_get_wtime();
		pr_matrix_posled(a, b, c, n1, m1, m2);
		timeout = omp_get_wtime();
		
		cout << "Последовательное: " << "\n";
		cout << "Сумма элементов матрицы, полученая после умножения: " << sum_elem_matrix(c, n1, m2) << "\n";
		cout << "Замер времени: " << fixed << timeout - timein << "\n" << "\n";
		for (int i = 0; i < n1; i++)
			for (int j = 0; j < m2; j++)
				c[i][j] = 0;



		// Ленточное разбиение
		double timein2, timeout2;
		timein2 = omp_get_wtime();
#pragma omp parallel num_threads(4) shared(a, b, c)
		{
			int numt = omp_get_num_threads();
			int q = n1 / numt;

			int tid = omp_get_thread_num();

			for (int i = 0; i < q; ++i)
				for (int j = 0; j < n1; ++j)
					for (int k = 0; k < n1; ++k)
					{
#pragma omp atomic
						c[i * numt + tid][j] += a[i * numt + tid][k] * b[k][j];
					}
		}
		timeout2 = omp_get_wtime();
		cout << "Ленточное разбиение: " << "\n";
		cout << "Сумма элементов матрицы, полученая после умножения: " << sum_elem_matrix(c, n1, n1) << "\n";
		cout << "Замер времени: " << fixed << timeout2 - timein2 << "\n" << "\n";
		for (int i = 0; i < n1; ++i)
			delete[] c[i];
		delete[] c;
		c = new int* [n1];
		for (int i = 0; i < n1; ++i)
			c[i] = new int[n1];

		int i, j, k;
		for (int i = 0; i < n1; i++)
			for (int j = 0; j < m2; j++)
				c[i][j] = 0;
		
		// Блочное разбиение
		double timein3, timeout3;
		timein3 = omp_get_wtime();
#pragma omp parallel num_threads(4) private(i, j, k) shared(a, b,c)
		{
			int numt = omp_get_num_threads();
			int q = n1 / 2;
			int p = n1 / 2;

			int tid1 = omp_get_thread_num() % 2;
			int tid2 = omp_get_thread_num() % 2;

			for (i = 0; i < q; ++i)
				for (j = 0; j < n1; ++j)
					for (k = 0; k < p; ++k)
					{
						#pragma omp atomic
						c[i * 2 + tid1][j] += a[i * 2 + tid1][k * 2 + tid2] * b[k * 2 + tid2][j];
					}
		}
		timeout3 = omp_get_wtime();
		cout << "Блочное разбиение: " << "\n";
		cout << "Сумма элементов матрицы, полученая после умножения: " << sum_elem_matrix(c, n1, n1) << "\n";
		cout << "Замер времени: " << fixed << timeout3 - timein3 << "\n" << "\n";
	}
}



