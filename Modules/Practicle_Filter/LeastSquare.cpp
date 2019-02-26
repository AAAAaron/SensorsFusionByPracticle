
// #include "LeastSquare.h"
#include "Practicle_Filter.h"

LeastSquare::LeastSquare(int length_array,double *x, double *y)
	: a(0)
	, b(0)
	, error_reg(0)
	, Rxy(0)
	, array_length(length_array)
	, navi_angle(0)
{
	    double t1=0, t2=0, t3=0, t4=0,t5=0;  
        for(int i=0; i<array_length; ++i)  
        {  
            t1 += x[i]*x[i];  
            t2 += x[i];  
            t3 += x[i]*y[i];  
            t4 += y[i];
			t5 += y[i]*y[i];
			//cout<<x[i]<<",";
			//cout<<y[i]<<",";
        }  
		//cout<<endl;
		if (t1*array_length - t2*t2 <1e-2)
		{
			Rxy = 0.0;
			Rxy3 =1.0;
			a = 0.0;
			b = 0.0;
			if (y[array_length-1]-y[0]>0){
				navi_angle = 0.0;
			}
			else
			{
				navi_angle = -M_PI;
			}

		}
		else if(t5*array_length - t4*t4<1e-2){
			Rxy = 0.0;
			Rxy3 =1.0;
			a = 0.0;
			b = 0.0;
			if (x[array_length-1]-x[0]>0){
				navi_angle = M_PI*0.5;
			}
			else
			{
				navi_angle = -M_PI*0.5;
			}

		}
		else{	
		Rxy = (t3*array_length - t2*t4)/sqrt(abs((t1*array_length - t2*t2)*(t5*array_length - t4*t4)));
		Rxy3 =Rxy;
		//double x_mean = t2/x.size();
		//double y_mean = t4/x.size();
        a = (t3*array_length - t2*t4) / (t1*array_length - t2*t2);  
        b = (t4 - a*t2) /array_length;  
        //a = (t1*t4 - t2*t3) / (t1*x.size() - t2*t2); 
		error(x,y,t4/array_length);
		cal_hxj(x,y);
		}

		return;
}


LeastSquare::~LeastSquare(void)
{
}


double LeastSquare::get_Y(double x) const
{
	 return a*x + b;  
}


void LeastSquare::print(void) const
{
	cout<<"y = "<<a<<"x + "<<b<<",error ="<<this->error_reg<<",rxy="<<Rxy<<"angle="<<navi_angle<<"\n"; 
}
        //vector<double> x;  
        //ifstream in(argv[1]);  
        //for(double d; in>>d; )  
        //    x.push_back(d);  
        //int sz = x.size();  
        //vector<double> y(x.begin()+sz/2, x.end());  
        //x.resize(sz/2);  
        //LeastSquare ls(x, y);  
//ls.print(); 

void LeastSquare::error(double *x, double *y,double y_)
{	
	double sst=0.0;
	for (int i = 0; i < array_length; ++i)
	{
		double g_y = get_Y(x[i]);
		double err = g_y-y[i];
		error_reg += err *err;
		sst += (y[i] - y_)*(y[i] - y_);
	}
	this->Rxy2 =1- error_reg/sst;
	return;
}



void LeastSquare::cal_hxj(double *x, double *y)
{

	navi_angle=atan(a);

	double cosresult= x[array_length-1]-x[0]+(y[array_length-1]-y[0])*a;
	if (cosresult<0){
		navi_angle=-M_PI/2-navi_angle;
	}   
    else{
        navi_angle=M_PI/2-navi_angle; 
	}

	return;
}
