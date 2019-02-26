#pragma once
#include "stdafx.h"
class LeastSquare
{
public:
	LeastSquare(int length_array,double *x, double *y);
	~LeastSquare(void);
	double a;
	double b;
	double get_Y(double x) const;
	void print(void) const;
	void error(double *x, double *y,double y_);
	double error_reg;
	double Rxy;
	double Rxy2;
	double Rxy3;
	int array_length;
	double navi_angle;
private:
	void cal_hxj(double *x, double *y);
};

