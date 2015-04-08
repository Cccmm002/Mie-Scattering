#ifndef _COMPLEX_H_
#define _COMPLEX_H_

#define _COMPLEX_DEFINED
#include <math.h>
#include <string>

typedef struct Complex
{
	double real;
	double imagine;

	Complex(double r = 0, double i = 0)
	{
		real = r;
		imagine = i;
	}

	std::string Complex::ToString()
	{
		if (imagine == 0) return std::to_string(real);
		else if (imagine > 0)
			return std::to_string(real) + "+" + std::to_string(imagine) + "i";
		else
			return std::to_string(real) + std::to_string(imagine) + "i";
	}

	inline double Complex::Abs()
	{
		return sqrt(real*real + imagine*imagine);
	}
} *pComplex;

inline Complex operator + (Complex &a, Complex &b)
{
	return Complex(a.real + b.real, a.imagine + b.imagine);
}

inline Complex operator - (Complex &a, Complex &b)
{
	return Complex(a.real - b.real, a.imagine - b.imagine);
}

inline Complex operator * (Complex &a, Complex &b)
{
	return Complex(a.real*b.real - a.imagine*b.imagine, a.real*b.imagine + a.imagine*b.real);
}

inline Complex operator * (double &a, Complex &b)
{
	return Complex(a*b.real, a*b.imagine);
}

inline Complex operator + (Complex &a, double &b)
{
	return Complex(a.real + b, a.imagine);
}

inline Complex operator - (Complex &a, double &b)
{
	return Complex(a.real - b, a.imagine);
}

inline Complex operator - (double &a, Complex &b)
{
	return Complex(a - b.real, 0 - b.imagine);
}

inline Complex operator / (Complex &a, Complex &b)
{
	double denominator = b.real*b.real + b.imagine*b.imagine;
	double r = a.real*b.real + a.imagine*b.imagine;
	double i = a.imagine*b.real - a.real*b.imagine;
	return Complex(r / denominator, i / denominator);
}

#endif