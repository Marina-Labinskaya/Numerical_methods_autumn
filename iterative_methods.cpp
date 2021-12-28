#include <iostream>
#include <locale.h>
#include <conio.h>
#include <iomanip>
#include <vector>
#include <math.h>

const double a = -1.0;
const double b = 1.0;
const double c = -1.0;
const double d = 1.0;
const double delta_u = -4.0;

struct solution_characteristics {
	std::vector<std::vector<double>> V;
	std::vector<std::vector<double>> U;
	std::vector<double> R_N;
	double eps_n;
	int N_steps;
	double R_N_2;
};

double f1(double x) {
	return (-x * x);
}

double f2(double y) {
	return (-y * y);
}

double u(double x, double y) {
	return (1 - x * x - y * y);
}

void print_grid(const std::vector<std::vector<double>>& V, int i) {
	int par = 12;
	double r = 1e+7;
	std::cout << i << "-ое приближение" << std::endl;
	for (int i = V.size() - 1; i >= 0; i--) {
		for (size_t j = 0; j < V.size(); j++) {
			std::cout << std::setw(par) << round(V[i][j] * r) / r;
		}
		std::cout << std::endl;
	}

}

solution_characteristics seidel_method(const std::vector<std::vector<double>>& V, int Nmax, double eps) {
	solution_characteristics res;
	res.V = V;
	res.eps_n = 1.0;

    int N = (res.V[0].size() - 2) * (res.V[0].size() - 2);
	double h = (b - a) / (res.V[0].size() - 1);
	double k = (d - c) / (res.V[0].size() - 1);
	double A = -2 * (1 / (h * h) + 1 / (k * k));
	std::vector<double> F;
	F.assign(N, delta_u);

	for (size_t i = 0; i < res.V[0].size() - 2; i++) {
		F[i] -= f1(a + (i + 1) * h) / (k * k);
		F[N - 1 - i] -= f1(b - (i + 1) * h) / (k * k);
		F[i * (res.V[0].size() - 2)] -= f2(c + (i + 1) * k) / ( h * h);
		F[i * (res.V[0].size() - 2) + (res.V[0].size() - 3)] -= f2(c + (i + 1) * k) / (h * h);
	}

	/*for (size_t i = 0; i < F.size(); i++)
		std::cout << std::setw(10) << F[i];
	std::cout << std::endl;
	std::cout << std::endl;*/

	std::vector<std::vector<double>> V_prev;
	V_prev = V;

	std::vector<double> v(N, 0.0);
	std::vector<std::vector<double>> A_kudr;
	A_kudr.assign(N, v);

	std::vector<double> R_N(N);

	for (size_t i = 0; i < A_kudr[0].size(); i++) {
		A_kudr[i][i] = A;

	    if (i % (res.V[0].size() - 2) == 0)
			A_kudr[i + 1][i] = 1 / (h * h);
		if (i % (res.V[0].size() - 2) == res.V[0].size() - 3)
			A_kudr[i - 1][i] = 1 / (h * h);
		if ((i % (res.V[0].size() - 2) != 0) && (i % (res.V[0].size() - 2) != res.V[0].size() - 3)) {
			A_kudr[i + 1][i] = 1 / (h * h);
			A_kudr[i - 1][i] = 1 / (h * h);
		}
		if (int(i / (res.V[0].size() - 2)) == 0)
			A_kudr[i][i + (res.V[0].size() - 2)] = 1 / (k * k);
		if (int(i / (res.V[0].size() - 2)) == res.V[0].size() - 3)
			A_kudr[i][i - (res.V[0].size() - 2)] = 1 / (k * k);
		if ((int(i / (res.V[0].size() - 2)) != 0) && (int(i / (res.V[0].size() - 2)) != res.V[0].size() - 3)) {
			A_kudr[i][i + (res.V[0].size() - 2)] = 1 / (k * k);
			A_kudr[i][i - (res.V[0].size() - 2)] = 1 / (k * k);
		}
	}

	std::vector<double> temp_v;

	/*for (size_t i = 0; i < A_kudr[0].size(); i++) {
		for (size_t j = 0; j < A_kudr[0].size(); j++)
			std::cout << std::setw(10) << A_kudr[i][j];
		std::cout << std::endl;
	}*/

	int iter = 0;

	print_grid(res.V, iter);
	iter++;

	while (iter <= Nmax && res.eps_n >= eps) {
		V_prev = res.V;
		for (size_t j = 1; j < res.V[0].size() - 1; j++) {
			for (size_t i = 1; i < res.V[0].size() - 1; i++) {
				if (i != 1 && i != res.V[0].size() - 2) {
					if (j != 1 && j != res.V[0].size() - 2)
					    res.V[j][i] = (F[(j - 1) * (res.V[0].size() - 2) + (i - 1)] - res.V[j][i + 1] / (h * h) - res.V[j][i - 1] / (h * h) - res.V[j + 1][i] / (k * k) - res.V[j - 1][i] / (k * k)) / A;
					if (j == 1) {
						res.V[j][i] = (F[(j - 1) * (res.V[0].size() - 2) + (i - 1)] - res.V[j][i + 1] / (h * h) - res.V[j][i - 1] / (h * h) - res.V[j + 1][i] / (k * k)) / A;
					}
					if (j == res.V[0].size() - 2) {
						res.V[j][i] = (F[(j - 1) * (res.V[0].size() - 2) + (i - 1)] - res.V[j][i + 1] / (h * h) - res.V[j][i - 1] / (h * h) - res.V[j - 1][i] / (k * k)) / A;
					}
				}
				if (i == 1) {
					if (j == 1)
						res.V[j][i] = (F[(j - 1) * (res.V[0].size() - 2) + (i - 1)] - res.V[j][i + 1] / (h * h) - res.V[j + 1][i] / (k * k)) / A;
					if (j == res.V[0].size() - 2)
						res.V[j][i] = (F[(j - 1) * (res.V[0].size() - 2) + (i - 1)] - res.V[j][i + 1] / (h * h) - res.V[j - 1][i] / (k * k)) / A;
					if (j != 1 && j != res.V[0].size() - 2)
						res.V[j][i] = (F[(j - 1) * (res.V[0].size() - 2) + (i - 1)] - res.V[j][i + 1] / (h * h) - res.V[j - 1][i] / (k * k) - res.V[j + 1][i] / (k * k)) / A;
				}
				if (i == res.V[0].size() - 2) {
					if (j == 1)
						res.V[j][i] = (F[(j - 1) * (res.V[0].size() - 2) + (i - 1)] - res.V[j][i - 1] / (h * h) - res.V[j + 1][i] / (k * k)) / A;
					if (j == res.V[0].size() - 2)
						res.V[j][i] = (F[(j - 1) * (res.V[0].size() - 2) + (i - 1)] - res.V[j][i - 1] / (h * h) - res.V[j - 1][i] / (k * k)) / A;
					if (j != 1 && j != res.V[0].size() - 2)
						res.V[j][i] = (F[(j - 1) * (res.V[0].size() - 2) + (i - 1)] - res.V[j][i - 1] / (h * h) - res.V[j + 1][i] / (k * k) - res.V[j - 1][i] / (k * k)) / A;
				}
			}
		}
		for (size_t i = 1; i < res.V[0].size() - 1; i++)
			for (size_t j = 1; j < res.V[0].size() - 1; j++) {
				res.eps_n = 0;
				if (abs(V_prev[i][j] - res.V[i][j]) >= res.eps_n)
					res.eps_n = abs(V_prev[i][j] - res.V[i][j]);
			}
		print_grid(res.V, iter);
		iter++;
	}
	res.N_steps = iter - 1;
	
	std::cout << std::endl;
	for (size_t j = 1; j < res.V[0].size() - 1; j++)
		for (size_t i = 1; i < res.V[0].size() - 1; i++) {
			temp_v.push_back(res.V[j][i]);
		}

	/*for (size_t i = 0; i < temp_v.size(); i++)
		std::cout << std::setw(12) << temp_v[i];
	std::cout << std::endl;*/

	for (size_t i = 0; i < A_kudr[0].size(); i++) {
		for (size_t j = 0; j < A_kudr[0].size(); j++) {
			R_N[i] += A_kudr[i][j] * temp_v[j];
		}
	}

	/*std::cout << std::endl;
	for (size_t i = 0; i < R_N.size(); i++)
		std::cout << std::setw(12) << R_N[i];*/

	for (size_t i = 0; i < R_N.size(); i++)
		R_N[i] -= F[i];
	
	/*std::cout << std::endl;
	for (size_t i = 0; i < R_N.size(); i++)
		std::cout << std::setw(12) << R_N[i];*/

	double sum = 0.0;

	for (size_t i = 0; i < R_N.size(); i++)
		sum += R_N[i] * R_N[i];
	res.R_N_2 = std::sqrt(sum);

	return res;
}

solution_characteristics seidel_method(const std::vector<std::vector<double>>& V, int Nmax, double eps, const std::vector<std::vector<double>>& U) {
	solution_characteristics res;
	res.V = V;
	res.eps_n = 1.0;

	int N = (res.V[0].size() - 2) * (res.V[0].size() - 2);
	double h = (b - a) / (res.V[0].size() - 1);
	double k = (d - c) / (res.V[0].size() - 1);
	double A = -2 * (1 / (h * h) + 1 / (k * k));
	std::vector<double> F;
	F.assign(N, delta_u);

	for (size_t i = 0; i < res.V[0].size() - 2; i++) {
		F[i] -= f1(a + (i + 1) * h) / (k * k);
		F[N - 1 - i] -= f1(b - (i + 1) * h) / (k * k);
		F[i * (res.V[0].size() - 2)] -= f2(c + (i + 1) * k) / (h * h);
		F[i * (res.V[0].size() - 2) + (res.V[0].size() - 3)] -= f2(c + (i + 1) * k) / (h * h);
	}

	std::vector<std::vector<double>> V_prev;
	V_prev = V;

	/*std::vector<double> v(N, 0.0);
	std::vector<std::vector<double>> A_kudr;
	A_kudr.assign(N, v);*/

	std::vector<double> R_N(N);

	/*for (size_t i = 0; i < A_kudr[0].size(); i++) {
		A_kudr[i][i] = A;

		if (i % (res.V[0].size() - 2) == 0)
			A_kudr[i + 1][i] = 1 / (h * h);
		if (i % (res.V[0].size() - 2) == res.V[0].size() - 3)
			A_kudr[i - 1][i] = 1 / (h * h);
		if ((i % (res.V[0].size() - 2) != 0) && (i % (res.V[0].size() - 2) != res.V[0].size() - 3)) {
			A_kudr[i + 1][i] = 1 / (h * h);
			A_kudr[i - 1][i] = 1 / (h * h);
		}
		if (int(i / (res.V[0].size() - 2)) == 0)
			A_kudr[i][i + (res.V[0].size() - 2)] = 1 / (k * k);
		if (int(i / (res.V[0].size() - 2)) == res.V[0].size() - 3)
			A_kudr[i][i - (res.V[0].size() - 2)] = 1 / (k * k);
		if ((int(i / (res.V[0].size() - 2)) != 0) && (int(i / (res.V[0].size() - 2)) != res.V[0].size() - 3)) {
			A_kudr[i][i + (res.V[0].size() - 2)] = 1 / (k * k);
			A_kudr[i][i - (res.V[0].size() - 2)] = 1 / (k * k);
		}
	}

	std::vector<double> temp_v;*/

	/*for (size_t i = 0; i < A_kudr[0].size(); i++) {
		for (size_t j = 0; j < A_kudr[0].size(); j++)
			std::cout << std::setw(10) << A_kudr[i][j];
		std::cout << std::endl;
	}*/

	int iter = 0;

	print_grid(res.V, iter);
	iter++;
	double eps1 = abs(U[1][1] - res.V[1][1]);
	for (size_t i = 1; i < res.V[0].size() - 1; i++)
		for (size_t j = 1; j < res.V[0].size() - 1; j++) {
			if (abs(U[i][j] - res.V[i][j]) >= res.eps_n)
				eps1 = abs(U[i][j] - res.V[i][j]);
		}

	while (iter <= Nmax && eps1 >= eps) {
		V_prev = res.V;
		for (size_t j = 1; j < res.V[0].size() - 1; j++) {
			for (size_t i = 1; i < res.V[0].size() - 1; i++) {
				if (i != 1 && i != res.V[0].size() - 2) {
					if (j != 1 && j != res.V[0].size() - 2)
						res.V[j][i] = (F[(j - 1) * (res.V[0].size() - 2) + (i - 1)] - res.V[j][i + 1] / (h * h) - res.V[j][i - 1] / (h * h) - res.V[j + 1][i] / (k * k) - res.V[j - 1][i] / (k * k)) / A;
					if (j == 1) {
						res.V[j][i] = (F[(j - 1) * (res.V[0].size() - 2) + (i - 1)] - res.V[j][i + 1] / (h * h) - res.V[j][i - 1] / (h * h) - res.V[j + 1][i] / (k * k)) / A;
					}
					if (j == res.V[0].size() - 2) {
						res.V[j][i] = (F[(j - 1) * (res.V[0].size() - 2) + (i - 1)] - res.V[j][i + 1] / (h * h) - res.V[j][i - 1] / (h * h) - res.V[j - 1][i] / (k * k)) / A;
					}
				}
				if (i == 1) {
					if (j == 1)
						res.V[j][i] = (F[(j - 1) * (res.V[0].size() - 2) + (i - 1)] - res.V[j][i + 1] / (h * h) - res.V[j + 1][i] / (k * k)) / A;
					if (j == res.V[0].size() - 2)
						res.V[j][i] = (F[(j - 1) * (res.V[0].size() - 2) + (i - 1)] - res.V[j][i + 1] / (h * h) - res.V[j - 1][i] / (k * k)) / A;
					if (j != 1 && j != res.V[0].size() - 2)
						res.V[j][i] = (F[(j - 1) * (res.V[0].size() - 2) + (i - 1)] - res.V[j][i + 1] / (h * h) - res.V[j - 1][i] / (k * k) - res.V[j + 1][i] / (k * k)) / A;
				}
				if (i == res.V[0].size() - 2) {
					if (j == 1)
						res.V[j][i] = (F[(j - 1) * (res.V[0].size() - 2) + (i - 1)] - res.V[j][i - 1] / (h * h) - res.V[j + 1][i] / (k * k)) / A;
					if (j == res.V[0].size() - 2)
						res.V[j][i] = (F[(j - 1) * (res.V[0].size() - 2) + (i - 1)] - res.V[j][i - 1] / (h * h) - res.V[j - 1][i] / (k * k)) / A;
					if (j != 1 && j != res.V[0].size() - 2)
						res.V[j][i] = (F[(j - 1) * (res.V[0].size() - 2) + (i - 1)] - res.V[j][i - 1] / (h * h) - res.V[j + 1][i] / (k * k) - res.V[j - 1][i] / (k * k)) / A;
				}
			}
		}
		res.eps_n = abs(V_prev[1][1] - res.V[1][1]);
		for (size_t i = 1; i < res.V[0].size() - 1; i++)
			for (size_t j = 1; j < res.V[0].size() - 1; j++) {
				if (abs(V_prev[i][j] - res.V[i][j]) >= res.eps_n)
					res.eps_n = abs(V_prev[i][j] - res.V[i][j]);
			}

		eps1 = abs(U[1][1] - res.V[1][1]);
		for (size_t i = 1; i < res.V[0].size() - 1; i++)
			for (size_t j = 1; j < res.V[0].size() - 1; j++) {
				if (abs(U[i][j] - res.V[i][j]) >= res.eps_n)
					eps1 = abs(U[i][j] - res.V[i][j]);
			}

		print_grid(res.V, iter);
		iter++;
	}
	res.N_steps = iter - 1;

	std::cout << std::endl;
	/*for (size_t j = 1; j < res.V[0].size() - 1; j++)
		for (size_t i = 1; i < res.V[0].size() - 1; i++) {
			temp_v.push_back(res.V[j][i]);
		}*/

	/*for (size_t i = 0; i < temp_v.size(); i++)
		std::cout << std::setw(12) << temp_v[i];
	std::cout << std::endl;*/

	/*for (size_t i = 0; i < A_kudr[0].size(); i++) {
		for (size_t j = 0; j < A_kudr[0].size(); j++) {
			R_N[i] += A_kudr[i][j] * temp_v[j];
		}
	}*/


	for (size_t j = 1; j < res.V[0].size() - 1; j++) {
		for (size_t i = 1; i < res.V[0].size() - 1; i++) {
			if (i != 1 && i != res.V[0].size() - 2) {
				if (j != 1 && j != res.V[0].size() - 2)
					R_N[(j - 1) * (res.V[0].size() - 2) + (i - 1)] = res.V[j][i + 1] / (h * h) + res.V[j][i - 1] / (h * h) + res.V[j + 1][i] / (k * k) + res.V[j - 1][i] / (k * k) + res.V[j][i] * A;
				if (j == 1) {
					R_N[(j - 1) * (res.V[0].size() - 2) + (i - 1)] = res.V[j][i + 1] / (h * h) + res.V[j][i - 1] / (h * h) + res.V[j + 1][i] / (k * k) + res.V[j][i] * A;
				}
				if (j == res.V[0].size() - 2) {
					R_N[(j - 1) * (res.V[0].size() - 2) + (i - 1)] = res.V[j][i + 1] / (h * h) + res.V[j][i - 1] / (h * h) + res.V[j - 1][i] / (k * k) + res.V[j][i] * A;
				}
			}
			if (i == 1) {
				if (j == 1)
					R_N[(j - 1) * (res.V[0].size() - 2) + (i - 1)] = res.V[j][i + 1] / (h * h) + res.V[j + 1][i] / (k * k) + res.V[j][i] * A;
				if (j == res.V[0].size() - 2)
					R_N[(j - 1) * (res.V[0].size() - 2) + (i - 1)] = res.V[j][i + 1] / (h * h) + res.V[j - 1][i] / (k * k) + res.V[j][i] * A;
				if (j != 1 && j != res.V[0].size() - 2)
					R_N[(j - 1) * (res.V[0].size() - 2) + (i - 1)] = res.V[j][i + 1] / (h * h) + res.V[j - 1][i] / (k * k) - res.V[j + 1][i] / (k * k) + res.V[j][i] * A;
			}
			if (i == res.V[0].size() - 2) {
				if (j == 1)
					R_N[(j - 1) * (res.V[0].size() - 2) + (i - 1)] = res.V[j][i - 1] / (h * h) + res.V[j + 1][i] / (k * k) + res.V[j][i] * A;
				if (j == res.V[0].size() - 2)
					R_N[(j - 1) * (res.V[0].size() - 2) + (i - 1)] = res.V[j][i - 1] / (h * h) + res.V[j - 1][i] / (k * k) + res.V[j][i] * A;
				if (j != 1 && j != res.V[0].size() - 2)
					R_N[(j - 1) * (res.V[0].size() - 2) + (i - 1)] = res.V[j][i - 1] / (h * h) + res.V[j + 1][i] / (k * k) - res.V[j - 1][i] / (k * k) + res.V[j][i] * A;
			}
		}
	}


	/*std::cout << std::endl;
	for (size_t i = 0; i < R_N.size(); i++)
		std::cout << std::setw(12) << R_N[i];*/

	for (size_t i = 0; i < R_N.size(); i++)
		R_N[i] -= F[i];

	/*std::cout << std::endl;
	for (size_t i = 0; i < R_N.size(); i++)
		std::cout << std::setw(12) << R_N[i];*/

	double sum = 0.0;

	for (size_t i = 0; i < R_N.size(); i++)
		sum += R_N[i] * R_N[i];
	res.R_N_2 = std::sqrt(sum);

	return res;
}

int main()
{

	setlocale(LC_ALL, "Rus");
	int N = 3;
	//int N;

	int Nmax;
	double eps;
	double Z_N = 0.0;

	std::cout << "Применение итерационного метода Зейделя. u(x,y) = 1 - x * x - y * y, -1 < x < 1, -1 < y < 1" << std::endl;
	std::cout << "delta_u = -4" << std::endl;
	std::cout << "u(-1, y) = - y * y" << std::endl;
	std::cout << "u(1, y) = - y * y" << std::endl;
	std::cout << "u(x, -1) = - x * x" << std::endl;
	std::cout << "u(x, 1) = - x * x" << std::endl;
	std::cout << "Сетка 3 X 3, начальное приближение - нулевое" << std::endl;

	//std::cout << "Введите размерность сетки N: ";
	//std::cin >> N;
	std::cout << "Введите Nmax: ";
	std::cin >> Nmax;
	std::cout << "Введите eps: ";
	std::cin >> eps;
	std::cout << std::endl;

	const double h = (b - a) / N;
	const double k = (d - c) / N;
	std::vector<double> v(N + 1);

	std::vector<std::vector<double>> U;
	U.assign(N + 1, v);

	for (int i = 0; i < N + 1; i++)
		for (int j = 0; j < N + 1; j++)
			U[i][j] = u(a + i * h, c + j * k);

	std::cout << "Точное решение задачи в узлах сетки:" << std::endl;
	for (int i = 0; i < N + 1; i++) {
		for (int j = 0; j < N + 1; j++)
			std::cout << std::setw(10) << U[i][j];
		std::cout << std::endl;
	}

	std::cout << std::endl;

	std::vector<std::vector<double>> V;
	V.assign(N + 1, v);

	V[0][0] = u(a, c);
	V[0][N] = u(a, d);
	V[N][0] = u(b, c);
	V[N][N] = u(b, d);

	for (int i = 1; i < N; i++) {
		V[i][0] = f1(a + i * h);
		V[i][N] = f1(a + i * h);
	}

	for (int j = 1; j < N; j++) {
		V[0][j] = f2(c + j * k);
		V[N][j] = f2(c + j * k);
	}

	for (int i = 1; i < N; i++)
		for (int j = 1; j < N; j++)
			V[i][j] = 0.0;

	solution_characteristics res = seidel_method(V, Nmax, eps, U);

	for (size_t i = 0; i < U[0].size(); i++)
	    for (size_t j = 0; j < U[0].size(); j++) {
			if (abs(U[i][j] - res.V[i][j]) > Z_N)
				Z_N = abs(U[i][j] - res.V[i][j]);
	}

	std::cout << std::endl;
	std::cout << "Выполнено N = " << res.N_steps << " шагов" << std::endl;
	std::cout << "Точность на выходе eps_N: " << res.eps_n << std::endl;
	std::cout << "Невязка на выходе имеет норму || R_N || _2 = " << res.R_N_2 << std::endl;
	std::cout << "Погрешность решения СЛАУ составит по норме || Z_N || _oo = " << Z_N;
	std::cout << std::endl;

	_getch();
}
