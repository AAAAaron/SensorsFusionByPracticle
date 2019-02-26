#include "Practicle_Filter.h"


#define RANDI(a,b) rand()%(b-a)+a//[a,b)的随机整数
#define Eur_distance(x1,y1,z1,x2,y2,z2) sqrt((x1-x2) *(x1-x2) + (y1-y2)*(y1-y2) +(z1-z2)*(z1-z2))
double round(double r)
{
    return (r > 0.0) ? floor(r + 0.5) : ceil(r - 0.5);
}

Practicle::Practicle()
	: practicle_num()
	, practicle_dimension()
	, weighted_position()
	, geometric_center_position()
	, practicle_diversity(-1000)
	, practicle_diversity_thershold(0.5)
	,wifi_weight_unit(1.0)
	,mag_weight_unit(0.5)
	,mag_match_num(7)
	,startx(-10)
	,starty(-10)
	,endx(80)
	,endy(33)
	, cannot_find_value(-1)
	, se_error(0)
{
	practicle_diversity_thershold=0.1;

}
//类里面写死了粒子滤波的状态只能是5个，索引了1,2,3

Practicle::~Practicle(void)
{
	//mkl_free(this->practicle_matrix);
  
}

void Practicle::Practicle_progration(double step_length, double deta_heading, Noise_practicle step_noise, Noise_practicle deta_heading_noise)
{
	for(int practicle_index =0 ;practicle_index<practicle_num ;practicle_index++){
		double practicle_step_length =step_length + step_noise.miu + step_length*(rand() / double(RAND_MAX)-0.5)  * step_noise.sigma*step_noise.sigma;
// 		double practicle_step_length =step_length *(practicle_index/practicle_num*0.5*2-0.5);
		double heading = -deta_heading + practicle_matrix(3,practicle_index) +deta_heading_noise.miu+ (rand() / double(RAND_MAX)-0.5)  *deta_heading_noise.sigma*deta_heading_noise.sigma;
		//float heading = - deta_heading + practicle_matrix(3,practicle_index) +deta_heading_noise.miu+ (gaussrand()-0.5)  *deta_heading_noise.sigma*deta_heading_noise.sigma;
		practicle_matrix(3,practicle_index) = heading;
		practicle_matrix(0,practicle_index) += practicle_step_length * sin(heading);
		practicle_matrix(1,practicle_index) += practicle_step_length * cos(heading);
	}

	calculate_weighted_position();
	return;
}

void Practicle::Init_practicle(int start_x,int start_y,int end_x,int end_y, Point init_point, double init_heading, int practicle_num, int practicle_dimension)
{
		
	practicle_diversity_thershold *= practicle_num;
		 vector<double> vec_ini(4,0.0);
	
		for (int i = 0; i < mag_match_num; i++)
		{
			vector<double> vectemp(vec_ini);
			mag_ob_buf.push_back(vectemp);
		}
		for (int i = 0; i < practicle_num; i++)
		{
			vector<vector<double > > vectemp(mag_ob_buf);
			mag_pro_buf.push_back(vectemp);

		}

		this->startx = start_x;
		this->starty = start_y;
		this->endx = end_x;
		this->endy = end_y;

		srand((unsigned)time(NULL));
	  //5*n
		this->practicle_dimension =practicle_dimension;
		this->practicle_num =practicle_num;
		this->practicle_matrix.setZero(practicle_dimension,practicle_num);
		for (int i = 0; i < practicle_num; i++){
			this->practicle_matrix(0,i) = init_point.x+(rand()/double(RAND_MAX)-0.5)*2;
			this->practicle_matrix(1,i) = init_point.y+(rand()/double(RAND_MAX)-0.5)*2;
			this->practicle_matrix(2,i) = init_point.z+(rand()/double(RAND_MAX)-0.5)*2;
			this->practicle_matrix(3,i) = init_heading;
			this->practicle_matrix(4,i) = 1.0/practicle_num;

		}  
		
		this->weighted_position = init_point;
		this->geometric_center_position = init_point;
		return;
}

void Practicle::Init_practicle(int start_x,int start_y,int end_x,int end_y, vector<wifi_match_result> wifidata_list, double init_heading, int practicle_num, int practicle_dimension)
{
	practicle_diversity_thershold *= practicle_num;
		this->startx = start_x;
		this->starty = start_y;
		this->endx = end_x;
		this->endy = end_y;

		//srand((unsigned)time(NULL));
	  //5*n
		this->practicle_dimension =practicle_dimension;
		this->practicle_num =practicle_num;
		this->practicle_matrix.setZero(practicle_dimension,practicle_num);
		for (int i = 0; i < practicle_num; i++){
			int wifi_index =practicle_num%wifidata_list.size();
			this->practicle_matrix(0,i) = wifidata_list[wifi_index].x + (rand()/double(RAND_MAX)-0.5)*wifidata_list[wifi_index].promise_length;
			this->practicle_matrix(1,i) = wifidata_list[wifi_index].y+(rand()/double(RAND_MAX)-0.5)*wifidata_list[wifi_index].promise_length;
			this->practicle_matrix(2,i) = wifidata_list[wifi_index].z+(rand()/double(RAND_MAX)-0.5)*wifidata_list[wifi_index].promise_length;
			this->practicle_matrix(3,i) = init_heading;
			this->practicle_matrix(4,i) = 1.0/practicle_num;

		}  
		
		calculate_weighted_position();
		return;
	
}

void Practicle::Init_practicle( vector<wifi_match_result> wifidata_list)
{
		//srand((unsigned)time(NULL));
		for (int i = 0; i < practicle_num; i++){
			int wifi_index =practicle_num%wifidata_list.size();
			this->practicle_matrix(0,i) = wifidata_list[wifi_index].x + (rand()/double(RAND_MAX)-0.5)*wifidata_list[wifi_index].promise_length;
			this->practicle_matrix(1,i) = wifidata_list[wifi_index].y+(rand()/double(RAND_MAX)-0.5)*wifidata_list[wifi_index].promise_length;
			this->practicle_matrix(2,i) = wifidata_list[wifi_index].z+(rand()/double(RAND_MAX)-0.5)*wifidata_list[wifi_index].promise_length;
			this->practicle_matrix(4,i) = 1.0/practicle_num;

		}  
		calculate_weighted_position();
		return;
	
}

void Practicle::Init_practicle(double heading)
{
		for (int i = 0; i < practicle_num; i++){
			this->practicle_matrix(3,i) = heading;
		}
		return;

}

void Practicle::Init_practicle(Point point)
{
		for (int i = 0; i < practicle_num; i++){
			this->practicle_matrix(0,i) = point.x;
			this->practicle_matrix(1,i) = point.y;
			this->practicle_matrix(2,i) = point.z;
		}  
		this->weighted_position = point;
		this->geometric_center_position = point;
		return;
}

void Practicle::calculate_weighted_position(int result_mode)
{
	//mode=1:默认权值结果，mode=2:几何平均，mode=3最大值ot(practicle_num,practicle_matrix,1,&practicle_matrix[practicle_num*4],1) ;
	//this->weighted_position.y = cblas_ddot(practicle_num,&practicle_matrix[practicle_num*1],1,&practicle_matrix[practicle_num*4],1) ;
	//this->weighted_position.z = cblas_ddot(practicle_num,&practicle_matrix[practicle_num*2],1,&practicle_matrix[practicle_num*4],1) ;
	//this->geometric_center_position.x = cblas_dasum(practicle_num,practicle_matrix,1)/practicle_num;
	//this->geometric_center_position.y = cblas_dasum(practicle_num,&practicle_matrix[practicle_num*1],1)/practicle_num;
	//this->geometric_center_position.z = cblas_dasum(practicle_num,&practicle_matrix[practicle_num*2],1)/practicle_num;
	//int maxweight_index = int(cblas_idamax (practicle_num,&(practicle_matrix[practicle_num*4]),1));
	//this->max_weight_position.x = practicle_matrix[practicle_num*0+maxweight_index];
	//this->max_weight_position.y = practicle_matrix[practicle_num*1+maxweight_index];
	//this->max_weight_position.z = practicle_matrix[practicle_num*2+maxweight_index];
	//this->se_error = SE_error(0.65,practicle_matrix,&(practicle_matrix[practicle_num]),practicle_num);
	if (result_mode==1)
	{
		this->weighted_position.x = practicle_matrix.row(0).dot(practicle_matrix.row(4)) ;	
		this->weighted_position.y = practicle_matrix.row(1).dot(practicle_matrix.row(4)) ;
		this->weighted_position.z = practicle_matrix.row(2).dot(practicle_matrix.row(4)) ;
	}
	if (result_mode==2)
	{
		this->geometric_center_position.x = practicle_matrix.row(0).mean();
		this->geometric_center_position.y = practicle_matrix.row(1).mean();
		this->geometric_center_position.z = practicle_matrix.row(2).mean();
	}


	int maxweight_index;
	double max_weight=practicle_matrix.row(4).maxCoeff(&maxweight_index);
// 	cout<<"max"<<practicle_matrix<<endl;
	this->max_weight_position.x = practicle_matrix(0,maxweight_index);
	this->max_weight_position.y = practicle_matrix(1,maxweight_index);
	this->max_weight_position.z = practicle_matrix(2,maxweight_index);
	calculate_se(0.65f);
	return;
}

void Practicle::calculate_se(float se_percent)
{
	vector<double > error_array;
	for (int i = 0; i < practicle_num; i++)
	{
		double x = this->practicle_matrix(0,i)-this->weighted_position.x;
		double y = this->practicle_matrix(1,i)-this->weighted_position.y;
		error_array.push_back(sqrt(x*x+y*y));
	}
	this->se_error = calculate_se_error(se_percent,error_array);
}

void Practicle::Rample_max2(void)
{
	practicle_diversity =0.0;
	practicle_diversity = practicle_matrix.row(4).dot(practicle_matrix.row(4));
	//for (int item = 0; item < practicle_num; item++)
	//{
	//	
	//	practicle_diversity += practicle_matrix(4,item) * practicle_matrix(4,item);
	//}
	//cout<<practicle_matrix.row(4)<<endl;
	if (practicle_diversity>=1000)
	{
		for (int i = 0; i < practicle_num; i++)
		{
			cout<<practicle_matrix(4,i)<<endl;
		}
	}
	if (practicle_diversity==0)
	{
		cout<<"all weight are 0,lost control"<<endl;
		weight_normalization(true);
	}
	else
	{
		practicle_diversity = 1.0/practicle_diversity;
		double w_thre=0.00001;
		double sum_weight = 0.0;
		//for (int item = 0; item < practicle_num; item++)
		//{
		//	sum_weight += practicle_matrix(4,item);
		//}
		sum_weight = practicle_matrix.row(4).sum();
		if (practicle_diversity < practicle_diversity_thershold){
			//cout<<"多样性不满足配置"<<endl;
			double weight_max = practicle_matrix.row(4).maxCoeff();
			weight_max*=2*rand()/double(RAND_MAX);
			for (int i = 0; i < practicle_num; i++)
			{
				if (practicle_matrix(4,i)<w_thre)
				{
					int index = RANDI(0,practicle_num);
					while (weight_max > practicle_matrix(4,index))
					{
						weight_max -= practicle_matrix(4,index);
						index +=1;
						if (index>=practicle_num)
						{
							index =0;
						}
					}
					practicle_matrix(0,i) = practicle_matrix(0,index);
					practicle_matrix(1,i) = practicle_matrix(1,index);
					practicle_matrix(2,i) = practicle_matrix(2,index);
					practicle_matrix(3,i) = practicle_matrix(3,index);
					practicle_matrix(4,i) = (1-sum_weight)/practicle_num;
				}
				else
				{
					practicle_matrix(4,i) *= sum_weight;
				}
			}
			
		}	
		
		weight_normalization();
		
	}
	calculate_weighted_position();
	//cout<<practicle_diversity<<"resample"<<endl;
	return;

}

void Practicle::Rample_max(void)
{
	practicle_diversity =0.0;
	for (int item = 0; item < practicle_num; item++)
	{
		
		practicle_diversity += practicle_matrix(4,item) * practicle_matrix(4,item);
	}

	if (practicle_diversity==0)
	{
		cout<<"全是0则重置，失锁"<<endl;
		weight_normalization(true);
	}
	else
	{
		practicle_diversity = 1.0/practicle_diversity;

		if (practicle_diversity < practicle_diversity_thershold)
		{
			for (int k = 0; k < practicle_num; k++)
			{
				double u = rand() /double(RAND_MAX);
				double qtempsum =0.0;
				for (int index = 0; index < practicle_num; index++)
				{
					qtempsum += practicle_matrix(4,k);
					if (qtempsum>=u)
					{
						practicle_matrix(0,k) = practicle_matrix(0,index);
						practicle_matrix(1,k) = practicle_matrix(1,index);
						practicle_matrix(2,k) = practicle_matrix(2,index);
						practicle_matrix(3,k) = practicle_matrix(3,index);
					}
				}
					


			}

			
		}	
		
		weight_normalization();
		
	}
	calculate_weighted_position();
	//cout<<practicle_diversity<<"resample"<<endl;
	return;

}

void Practicle::weight_normalization(bool reset)
{//reset代表全部重置为1/n
	
	if (reset == true){
		double weight =1.0 /practicle_num;
		for (int i = 0; i < practicle_num; i++){
			this->practicle_matrix(4,i) = weight;
		}
	}
	else
	{
		double sum_weight = practicle_matrix.row(4).sum();
		for (int item = 0; item < practicle_num; item++)
		{
			practicle_matrix(4,item) =  practicle_matrix(4,item)/sum_weight;
		}


	}
	return;
}


void Practicle::load_fp(string filename ,string fp_name,string floor_id)
{
// 		int fp_row = 0;
// 		int fp_col =0;
// 		fp_row = endy - starty ;
// 		fp_col = endx - startx ;
// 		    // 读文件  
// 		ifstream inFile(filename, ios::in);  
// 		string lineStr;  
// 
// 		vector<vector<double>> strArray;  
// 		while (getline(inFile, lineStr))  
// 		{  
// 			// 打印整行字符串  
// 			//cout << lineStr << endl;  
// 			// 存成二维表结构  
// 			stringstream ss(lineStr);  
// 			string str;  
// 			vector<double> lineArray;  
// 			// 按照逗号分隔  
// 			while (getline(ss, str, ','))  
// 				lineArray.push_back(atof(str.c_str()));  
// 			strArray.push_back(lineArray);  
// 		}  
// 		
// 		if (fp_name == "space"){
// 			//this->var_id_floorid.insert(make_pair(this->space_fp.size(),floor_id));
// 			this->var_id_floorid[floor_id+fp_name] =this->space_fp.size();
// 			this->space_fp.push_back(strArray);
// 			return;
// 		}
// 		if (fp_name == "mag1")
// 		{
// 			this->var_id_floorid[floor_id+fp_name] =this->mag1_fp.size();
// 			this->mag1_fp.push_back(strArray);
// 			return;
// 		}
// 		if (fp_name == "mag2")
// 		{
// 			this->var_id_floorid[floor_id+fp_name] =this->mag2_fp.size();
// 			this->mag2_fp.push_back(strArray);
// 			return;
// 		}
// 		if (fp_name == "mag3")
// 		{
// 			this->var_id_floorid[floor_id+fp_name] =this->mag3_fp.size();
// 			this->mag3_fp.push_back(strArray);
// 			return;
// 		}
// 		if (fp_name == "mag4")
// 		{
// 			this->var_id_floorid[floor_id+fp_name] =this->mag4_fp.size();
// 			this->mag4_fp.push_back(strArray);
// 			return;
// 		}
// 		return;

}

void Practicle::test_loadfp()
{
	//for (unsigned int i =0;i<fpdata.size();i++){
	//	for (unsigned int j = 0; j < fpdata[i].size(); j++)
	//	{
	//		cout<<fpdata[i][j];
	//	}
	//	cout<<endl;
	//}
	map<string,int>::iterator iter;  
	for(iter = var_id_floorid.begin(); iter != var_id_floorid.end(); iter++) { 
  
       cout<<iter->first<<' '<<iter->second<<endl; 
	}
	//cout<<"space has"<<this->space_fp.size()<<endl;
	//cout<<"mag has"<<this->mag1_fp.size()<<endl;
	return;
}

//类似于传入时序的蓝牙信息，那么就使用这种方式
double Practicle::weight_from_single_signal(double x,double y,vector<wifi_match_result> signal_list)
{
	double weight_unit =  1.0 / practicle_num ;
	double weight =0;
	for (unsigned int sig_index = 0; sig_index < signal_list.size(); sig_index++)
	{
			double p_p_distance = Eur_distance(signal_list[sig_index].x,signal_list[sig_index].y,signal_list[sig_index].z,x,y,0);
		
			weight += signal_list[sig_index].probability * weight_unit * 1 / pow(1.1, p_p_distance);
	}
	
	return weight;
	//Rample_max2();

}

//最好先用空间判断权值，后面的权值都可以根据这个调整完之后，如果是0就不判断了
void Practicle::adjust_weight_from_space(string floor_id)
{
	this->startx = int(all_fp[floor_id+"xxyy"][0][0]);
	this->starty = int(all_fp[floor_id+"xxyy"][0][2]);
	this->endx = int(all_fp[floor_id+"xxyy"][0][1]);
	this->endy  = int(all_fp[floor_id+"xxyy"][0][3]);
	//乘法表示在空间上的一票否决制
	for (int index = 0; index < practicle_num; index++)
	{
		double x =practicle_matrix(0,index);
		double y = practicle_matrix(1,index);
		if (x>=startx&&x<=endx&&y>=starty&&y<=endy){
		practicle_matrix(4,index) *= weight_from_space(x,y,floor_id);
		}
		else
		{
			practicle_matrix(4,index) =0;
		}
	}
	return;
}


vector<double> Practicle::magnetic_match(double x,double y,string floor_id,vector<double> mag_info_obs)
{
	int indi;
    int indj;
	double mag_p1;
    double mag_p2;
	double mag_p3;
	double mag_p4;
	double probability;
	vector<double> result;

	


	indi=int(round(y-starty));
	indj=int(round(x-startx));
	mag_p1=this->all_fp[floor_id+"mag1"][indi][indj];
	mag_p2=this->all_fp[floor_id+"mag2"][indi][indj];
	mag_p3=this->all_fp[floor_id+"mag3"][indi][indj];
	mag_p4=this->all_fp[floor_id+"mag4"][indi][indj];

	double indexexp = max(sqrt((mag_info_obs[1]-mag_p2)*(mag_info_obs[1]-mag_p2)+(mag_info_obs[2]-mag_p3)*(mag_info_obs[2]-mag_p3)),0.0);
	probability = exp(-indexexp/4);
	if (mag_p1<1e-5&&mag_p2<1e-5&&mag_p3<1e-5&&mag_p4<1e-5){
		probability=0.0;
	}
	result.push_back(probability);
	result.push_back(mag_p1);
	result.push_back(mag_p2);
	result.push_back(mag_p3);
	result.push_back(mag_p4);
	
	return result;
}


double Practicle::magnetic_match_s(int practicle_index,double p0)
{
		//double mag_ob1_av = average(mag_ob_buf[0],mag_match_num);
		//double mag_ob2_av = average(mag_ob_buf[1],mag_match_num);
		//double mag_ob3_av = average(mag_ob_buf[2],mag_match_num);
		//double mag_ob4_av = average(mag_ob_buf[3],mag_match_num);
		//double mag_pro1_av=average(mag_pro_buf[practicle_index][0],mag_match_num);
		//double mag_pro2_av=average(mag_pro_buf[practicle_index][1],mag_match_num);
		//double mag_pro3_av=average(mag_pro_buf[practicle_index][2],mag_match_num);
		//double mag_pro4_av=average(mag_pro_buf[practicle_index][3],mag_match_num);
		double *mag_obs;
		mag_obs= (double *)malloc(mag_match_num * 4 * sizeof(double));
		double *mag_pro;
		mag_pro= (double *)malloc(mag_match_num * 4 * sizeof(double));
		
		for (int j = 0; j < 4; j++)
		{
			double av_temp=0.0,av2_temp=0.0;
			for (int i = 0; i < mag_match_num; i++)
			{
				av_temp+=mag_ob_buf[i][j];
				av2_temp+=mag_pro_buf[practicle_index][i][j];
			}
			av_temp/=mag_match_num;
			av2_temp/=mag_match_num;

			for (int i = 0; i < mag_match_num; i++)
			{
				mag_obs[i+mag_match_num*j] = mag_ob_buf[i][j] - av_temp;
				mag_pro[i+mag_match_num*j] = mag_pro_buf[practicle_index][i][j] - av2_temp;
			}
		}
		double w_mfactor;

		LeastSquare ls1 (mag_match_num,&(mag_obs[mag_match_num]), &(mag_pro[mag_match_num]));
		LeastSquare ls2 (mag_match_num,&(mag_obs[mag_match_num*2]), &(mag_pro[mag_match_num*2]));
		LeastSquare ls3 (mag_match_num,&(mag_obs[mag_match_num*3]),&(mag_pro[mag_match_num*3]));
		//p1=corr(mag_obv_series2,mag_fin_series2,'type','pearson');  
		//p2=corr(mag_obv_series3,mag_fin_series3,'type','pearson');  
		//p3=corr(mag_obv_series4,mag_fin_series4,'type','pearson'); 
		
		w_mfactor=p0+max(ls3.Rxy,0.0)*max(ls1.Rxy,0.0)*max(ls2.Rxy,0.0)*1.5;
		if (w_mfactor<0.2)
		{
			w_mfactor = 0.0;
		}
								if (std::isnan(w_mfactor))
								{
								  cout<<" 系数不存在"<<endl;
									for (int i = 0; i < practicle_num; i++)
									{
										cout<<practicle_matrix(4,i)<<endl;
									}
								}
		free(mag_obs);
		free(mag_pro);
	return w_mfactor;
}


void Practicle::adjust_weight_mag_space(int step_index, string floor_id,vector<double> mag_info_obs)
{
	this->startx = int(all_fp[floor_id+"xxyy"][0][0]);
	this->starty = int(all_fp[floor_id+"xxyy"][0][2]);
	this->endx = int(all_fp[floor_id+"xxyy"][0][1]);
	this->endy  = int(all_fp[floor_id+"xxyy"][0][3]);

	mag_ob_buf.push_back(mag_info_obs);
	mag_ob_buf.erase(mag_ob_buf.begin());
	

	for (int prac_index = 0; prac_index < practicle_num; prac_index++)
	{
		double sapce_weight ,p0;
		double x =practicle_matrix(0,prac_index);
		double y = practicle_matrix(1,prac_index);
		vector<double> temp_result ;

			if (x>=startx&&x<=endx&&y>=starty&&y<=endy)
			{

				sapce_weight=weight_from_space(x,y,floor_id); 
				practicle_matrix(4,prac_index) *= sapce_weight;



				if (practicle_matrix(4,prac_index)<1e-8){
					practicle_matrix(4,prac_index)=0.0;
					//continue;
				}
				else
				{//TODO 确定顺序是否正确

				

				temp_result	= magnetic_match(x,y,floor_id,mag_info_obs);
				p0 = temp_result[0];
				temp_result.erase(temp_result.begin());
				mag_pro_buf[prac_index].push_back(temp_result);
				mag_pro_buf[prac_index].erase(mag_pro_buf[prac_index].begin());
				
						
				if (step_index>=mag_match_num&&p0!=0.0)
				{
					p0 = magnetic_match_s(prac_index,p0);
				}

				temp_result.clear();
				practicle_matrix(4,prac_index) +=  p0*1.0/practicle_num;
				}
				//cout<<practicle_matrix(4,prac_index)<<endl;
			}
			else
			{

				practicle_matrix(4,prac_index)=0.0;
			}


	}
	return;

}

void Practicle::adjust_weight_mag_space_wifi(int step_index, string floor_id,vector<double> mag_info_obs,vector<wifi_match_result> wifidata_list)
{
	this->startx = int(all_fp[floor_id+"xxyy"][0][0]);
	this->starty = int(all_fp[floor_id+"xxyy"][0][2]);
	this->endx = int(all_fp[floor_id+"xxyy"][0][1]);
	this->endy  = int(all_fp[floor_id+"xxyy"][0][3]);
	mag_ob_buf.push_back(mag_info_obs);
	mag_ob_buf.erase(mag_ob_buf.begin());
	for (int prac_index = 0; prac_index < practicle_num; prac_index++)
	{
		double sapce_weight ,wifi_weight,p0;
		double x =practicle_matrix(0,prac_index);
		double y = practicle_matrix(1,prac_index);
		vector<double> temp_result ;

			if (x>=startx&&x<=endx&&y>=starty&&y<=endy)
			{

				sapce_weight=weight_from_space(x,y,floor_id); 

				practicle_matrix(4,prac_index) *= sapce_weight;



				if (practicle_matrix(4,prac_index)<1e-8){
					practicle_matrix(4,prac_index)=0.0;
					//continue;
				}
				else
				{//TODO 确定顺序是否正确
				wifi_weight =1/practicle_num * weight_from_wifi(x,y,wifidata_list);		
				practicle_matrix(4,prac_index) += wifi_weight;
				

				temp_result	= magnetic_match(x,y,floor_id,mag_info_obs);
				p0 = temp_result[0];
				temp_result.erase(temp_result.begin());
				mag_pro_buf[prac_index].push_back(temp_result);
				mag_pro_buf[prac_index].erase(mag_pro_buf[prac_index].begin());
				
						
				if (step_index>=mag_match_num&&p0!=0.0)
				{
					p0 = magnetic_match_s(prac_index,p0);

				}

				temp_result.clear();
				practicle_matrix(4,prac_index) +=  p0 *this->mag_weight_unit/practicle_num;
				}
				//cout<<practicle_matrix(4,prac_index)<<endl;
			}
			else
			{

				practicle_matrix(4,prac_index)=0.0;
			}


	}

	return;

}


//前面导致了多次循环，所以重写成单独的函数

int Practicle::Var_from_floor_id(string floor_id)
{
	map<string,int>::iterator fp_now;
	fp_now = this->var_id_floorid.find(floor_id);
	if (fp_now==var_id_floorid.end())
	{
		cout<<floor_id<<"not found"<<endl;
		return cannot_find_value;
	}
	else
	{
		return fp_now->second;
	}
}

double Practicle:: weight_from_space(double x,double y,string floor_id)
{
	double weight = 0.0;
	this->startx = int(all_fp[floor_id+"xxyy"][0][0]);
	this->starty = int(all_fp[floor_id+"xxyy"][0][2]);
	this->endx = int(all_fp[floor_id+"xxyy"][0][1]);
	this->endy  = int(all_fp[floor_id+"xxyy"][0][3]);
	if (x>=startx&&x<=endx&&y>=starty&&y<=endy)
	{
		int indi;
		int indj;
		indi=int(round(y-starty));
		indj=int(round(x-startx));
		weight = this->all_fp[floor_id+"space"][indi][indj];

	}
	if (weight>1)
	{
		weight =0;
	}
	return weight;
}

double Practicle::weight_from_wifi(double x,double y,vector<wifi_match_result> wifidata_list)
{
	double mblen=0.4;
	for (unsigned int index = 0; index < wifidata_list.size(); index++)
	{
		double distance = Eur_distance(x,y,0,wifidata_list[index].x,wifidata_list[index].y,0);
		if (distance<wifidata_list[index].promise_length)
		{
			if (wifidata_list[index].probability>mblen)
			{
				mblen =wifidata_list[index].probability;
			} 
		}
		
	}
	if (mblen==0.4)//小于0.4直接截断
	{
		return 0.0*this->wifi_weight_unit;
	}
	else
	{
		return mblen*this->wifi_weight_unit;
	}
	
}


void Practicle::adjust_weight_space_wifi(int step_index, string floor_id, vector<wifi_match_result> wifidata_list)
{
  
	for (int prac_index = 0; prac_index < practicle_num; prac_index++)
	{
		practicle_matrix(4,prac_index) *= weight_from_space(practicle_matrix(0,prac_index),practicle_matrix(1,prac_index),floor_id);
		if (practicle_matrix(4,prac_index) >1e-5)
		{
		  practicle_matrix(4,prac_index) += 1/practicle_num* weight_from_wifi(practicle_matrix(0,prac_index),practicle_matrix(1,prac_index),wifidata_list);	
		}

	}
	return;

}

void Practicle::adjust_weight_space_single_signal(string floor_id, vector<wifi_match_result> signal_list)
{
	for (int prac_index = 0; prac_index < practicle_num; prac_index++)
	{
				practicle_matrix(4,prac_index) *= weight_from_space(practicle_matrix(0,prac_index),practicle_matrix(1,prac_index),floor_id);
				if (practicle_matrix(4,prac_index) >1e-5)
				{
				practicle_matrix(4,prac_index) += 0.3/practicle_num* weight_from_single_signal(practicle_matrix(0,prac_index),practicle_matrix(1,prac_index),signal_list);	
				}

	}
	return;

}

double Practicle::calculate_se_error(double SE_para, vector<double > error_array)
{

	sort(error_array.begin(),error_array.end());	

	return error_array[int(SE_para*100)];
}


void Practicle::reset(void)
{
	//mkl_free(this->practicle_matrix);
}


double Practicle::gaussrand(void)
{
	static double V1, V2, S;
    static int phase = 0;
    double X;
     
    if ( phase == 0 ) {
        do {
            double U1 = (double)rand() / RAND_MAX;
            double U2 = (double)rand() / RAND_MAX;
             
            V1 = 2 * U1 - 1;
            V2 = 2 * U2 - 1;
            S = V1 * V1 + V2 * V2;
        } while(S >= 1 || S == 0);
         
        X = V1 * sqrt(-2 * log(S) / S);
    } else
        X = V2 * sqrt(-2 * log(S) / S);
         
    phase = 1 - phase;
 
    return X;
}
