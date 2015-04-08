#include "Mie.h"
using namespace std;

void SphereMie(double x, Complex m, pComplex s1, pComplex s2, double *q_sca, double *q_ext)
{
	int n_max;
	if (x < 8)
		n_max = (int)(x + 4 * pow(x, 1 / 3) + 1);
	else if (x < 4200)
		n_max = (int)(x + 4.05*pow(x, 1 / 3) + 2);
	else
		n_max = (int)(x + 4 * pow(x, 1 / 3) + 2);
	n_max++;

	//Calculate D(n, x) downward
	//using recurrence formula D(n-1, x)=n/x-1/(D(n,x)+n/x)
	Complex mx = x*m;
	pComplex dn = (pComplex)malloc(sizeof(Complex)*(n_max + 1));
	dn[n_max] = Complex(0, 0);
	for (int i = n_max; i > 0; i--)
		dn[i - 1] = Complex(i) / mx - Complex(1) / (dn[i] + (Complex(i) / mx));

	//Calculate Spherical Bessel functions 
	double *jn = (double*)malloc(sizeof(double)*n_max);
	double *yn = (double*)malloc(sizeof(double)*n_max);
	jn[0] = sin(x) / x;
	jn[1] = sin(x) / (x*x) - cos(x) / x;
	yn[0] = 0 - cos(x) / x;
	yn[1] = 0 - cos(x) / (x*x) - sin(x) / x;
	for (int i = 2; i < n_max; i++)
	{
		jn[i] = ((2 * i - 1) / x)*jn[i - 1] - jn[i - 2];
		yn[i] = ((2 * i - 1) / x)*yn[i - 1] - yn[i - 2];
	}

	double qst = 0;

	//Calculate scattering coefficients a(n) and b(n)
	pComplex an = (pComplex)malloc(sizeof(Complex)*n_max);
	pComplex bn = (pComplex)malloc(sizeof(Complex)*n_max);
	for (int i = 1; i < n_max; i++)
	{
		double psi1 = x*jn[i], psi0 = x*jn[i - 1];
		Complex xi1 = Complex(psi1, x*yn[i]), xi0 = Complex(psi0, x*yn[i - 1]);
		an[i] = (psi1*(dn[i] / m + Complex(i / x)) - psi0) / (xi1*(dn[i] / m + Complex(i / x)) - xi0);
		bn[i] = (psi1*(dn[i] * m + Complex(i / x)) - psi0) / (xi1*(dn[i] * m + Complex(i / x)) - xi0);

		double a_s = an[i].Abs(), b_s = bn[i].Abs();
		qst += (2 * i + 1)*(a_s*a_s + b_s*b_s);
	}

	free(dn); free(jn); free(yn);

	//Calculate Scattering matrix elements for each angle
	double pi = 4 * atan(1);
	double interval = pi / NUM_ANGLE;
	int count = 0;
	for (double angle = 0; angle < pi; angle += interval)
	{
		//Calculate pi and tao
		double *pn = (double*)malloc(sizeof(double)*n_max);
		double *tn = (double*)malloc(sizeof(double)*n_max);
		pn[0] = 0; pn[1] = 1;
		for (int i = 1; i < n_max; i++)
		{
			if (i>1)
				pn[i] = ((2 * i - 1) / (i - 1))*cos(angle)*pn[i - 1] - (i / (i - 1))*pn[i - 2];
			tn[i] = i*cos(angle)*pn[i] - (i + 1)*pn[i - 1];
		}

		Complex sum1 = Complex(0, 0), sum2 = Complex(0, 0);
		for (int i = 1; i < n_max; i++)
		{
			double coe = (double)(2 * i + 1) / (double)(i*i + i);
			sum1 = sum1 + coe*(pn[i] * an[i] + tn[i] * bn[i]);
			sum2 = sum2 + coe*(tn[i] * an[i] + pn[i] * bn[i]);
		}
		s1[count] = sum1;
		s2[count] = sum2;
		count++;
		free(pn); free(tn);
	}
	free(an); free(bn);

	*q_ext = (4.0 / (x*x))*(s1[0].real);
	*q_sca = (2.0 / (x*x))*qst;
}