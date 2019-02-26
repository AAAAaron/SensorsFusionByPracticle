
#include "stdafx.h"
#include "cnew.h"
// 根据地磁调整
//根据wifi调整
//暴露在外的函数用大写
class Practicle
{
public:
	//公共参数
	//Practicle();
	Practicle();
	~Practicle(void);
	//double *practicle_matrix;
	MatrixXd practicle_matrix;
	Point weighted_position;
	Point geometric_center_position;
	Point max_weight_position;
	map<string,int> var_id_floorid;
	double practicle_diversity;
	double practicle_diversity_thershold;
	int mag_match_num;
	double se_error;
	map<string,vector<vector<double> > > all_fp;//维护的指纹,map 赋值是引用吗，因为多线程调用的时候不确定会不会造成错误
	//共有函数
	void Practicle_progration(double step_length, double deta_heading, Noise_practicle step_noise, Noise_practicle deta_heading_noise);
	void Init_practicle(int start_x,int start_y,int end_x,int end_y,Point init_point, double init_heading, int praciticle_num, int practicle_dimension);
	void Init_practicle(double heading);
	void Init_practicle(Point init_point);
	void Init_practicle(int start_x,int start_y,int end_x,int end_y, vector<wifi_match_result> wifidata_list, double init_heading, int practicle_num, int practicle_dimension);
	void Init_practicle( vector<wifi_match_result> wifidata_list);

	void Rample_max2(void);
	void Rample_max(void);
	//

	void weight_normalization(bool reset=false);
	double weight_from_single_signal(double x,double y,vector<wifi_match_result> signal_list);
	void calculate_weighted_position(int result_mode=1);
	void load_fp(string filename,string fp_name,string floor_id);
	void test_loadfp();
	void adjust_weight_from_space(string floor_id);
	double magnetic_match_s(int practicle_index,double p0);
	void adjust_weight_mag_space(int step_index,string floor_id,vector<double> mag_info_obs);
	void adjust_weight_mag_space_wifi(int step_index, string floor_id,vector<double> mag_info_obs,vector<wifi_match_result> wifidata_list);
	double weight_from_space(double x,double y,string floor_id);
	void adjust_weight_space_wifi(int step_index, string floor_id, vector<wifi_match_result> wifidata_list);
	void adjust_weight_space_single_signal(string floor_id, vector<wifi_match_result> signal_list);
	double calculate_se_error(double SE_para, vector<double > error_array);
	void calculate_se(float se_percent);
	void reset(void);

private:
	//参数私有
	int practicle_num;
	int practicle_dimension;
	vector<vector<vector<double > > > mag1_fp;//0,1,2,3 地磁的四个floor 0,1,2,3;floor,0,1,2,3....
	vector<vector<vector<double> > > mag2_fp;
	vector<vector<vector<double> > > mag3_fp;
	vector<vector<vector<double> > > mag4_fp;
	vector<vector<vector<double> > > space_fp;//空间指纹
	
	vector<vector<double> > mag_ob_buf;//观测序列7-4
	vector<vector<vector<double> > > mag_pro_buf;//历史粒子的概率[100][7][4]
	int startx;
	int starty;
	int endx;
	int endy;
	int cannot_find_value;

	//私有函数
	vector<double> magnetic_match(double x,double y,string floor_id,vector<double> mag_info_obs);
	int Var_from_floor_id(string floor_id);
	double weight_from_wifi(double x,double y,vector<wifi_match_result> wifidata_list);
	double gaussrand(void);
};

