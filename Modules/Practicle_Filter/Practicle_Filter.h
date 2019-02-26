#pragma once
#include "stdafx.h"
#define Point_List_length 10

#include <algorithm>
#include <vector> 
#include <map>



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

class Practicle
{
public:
	//��������
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
	double wifi_weight_unit;//һ����ԪȨֵ����
	double mag_weight_unit;
	int mag_match_num;
	double se_error;
	map<string,vector<vector<double> > > all_fp;//ά����ָ��,map ��ֵ����������Ϊ���̵߳��õ�ʱ��ȷ���᲻����ɴ���
	//���к���
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
	//����˽��
	int practicle_num;
	int practicle_dimension;
	vector<vector<vector<double > > > mag1_fp;//0,1,2,3 �شŵ��ĸ�floor 0,1,2,3;floor,0,1,2,3....
	vector<vector<vector<double> > > mag2_fp;
	vector<vector<vector<double> > > mag3_fp;
	vector<vector<vector<double> > > mag4_fp;
	vector<vector<vector<double> > > space_fp;//�ռ�ָ��
	
	vector<vector<double> > mag_ob_buf;//�۲�����7-4
	vector<vector<vector<double> > > mag_pro_buf;//��ʷ���ӵĸ���[100][7][4]
	int startx;
	int starty;
	int endx;
	int endy;
	int cannot_find_value;

	//˽�к���
	vector<double> magnetic_match(double x,double y,string floor_id,vector<double> mag_info_obs);
	int Var_from_floor_id(string floor_id);
	double weight_from_wifi(double x,double y,vector<wifi_match_result> wifidata_list);
	double gaussrand(void);
};

class WifiModule
{
public:
	unsigned int array_length;
	int offset_rssi;
	map<string,vector<vector<string> > > wifi_list_fp;//wifiָ���ж�ά�б���str
	map<string,vector<vector<double> > >  wifi_rssi_fp;//���¹���wifi����
	vector<string> region_list;//�ڽ��������б�
	wifi_match_result wifi_aftermatch_result;
	vector<wifi_match_result> wifi_cal_list;
	WifiModule(void);
	~WifiModule(void);
	void WifiMach(vector<wifi> wifi_receive,string region_name);
	void WifiMach(vector<wifi> wifi_receive);

private:
	void vec2matrix(void);
};

class Practicle_Filter
{
public:
	vector<double> x_list;
	vector<double> y_list;
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
	WifiModule wifi_obj;
	map<string,vector<vector<double> > > all_fp;//ά����ָ��,map ��ֵ����������Ϊ���̵߳��õ�ʱ��ȷ���᲻����ɴ���
	map<string,vector<vector<string> > > wifi_list_fp;//wifiָ���ж�ά�б���str
	map<string,vector<vector<double> > > wifi_rssi_fp;//���¹���wifi����
	vector<string> region_list;//�ڽ��������б�
	map<string ,Practicle> practicle_map;//ά�������Ӷ��󣬶��,�ⲿ�������Ƿ��ظ����ظ�ֱ�Ӹ���
	Practicle_Filter(void);
	~Practicle_Filter(void);
	void load_fp_all(string filedir);
	void InitPracClass(string prac_name);
	void main_free(string prac_name);
	void main_practicle_init(string prac_name,string floor_name,double heading,int practicle_num,int practicle_dis,vector<wifi_match_result> wifidata_list);
	void main_practicle_init(string prac_name,string floor_name,double x,double y,double z,double heading,int practicle_num,int practicle_dis);
	void main_practicle_progration(string prac_name,double step_length, double step_deta_angle, double step_miu,double step_sigma, double heading_miu,double heading_sigma);
	void adjust_practicle_angle(string prac_name);
	Point main_progration(string prac_name,string floor_id,double step_length,double step_deta_angle,double Heading);
	Point main_progration(string prac_name,int step_index,string floor_id, double step_length,double step_deta_angle,double Heading, double mag_norm,double mag_z,double mag_y,double inclination);//�شſռ�
	Point main_progration(string prac_name,int step_index,string floor_id,double step_length,double step_deta_angle,double Heading, vector<wifi_match_result> wifidata_list);//����wifi�źŵ�����
	Point main_progration(string prac_name,int step_index,string floor_id,double step_length,double step_deta_angle,double Heading, double mag_norm,double mag_z,double mag_y,double inclination,vector<wifi_match_result> wifidata_list);

private:
	vector<string> getFiles( string cate_dir) ;
	
	void test_loadfp();

	
};