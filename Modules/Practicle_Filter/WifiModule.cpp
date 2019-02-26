
// #include "WifiModule.h"
#include "Practicle_Filter.h"

WifiModule::WifiModule()
{
	array_length=5;
	offset_rssi=110;
}


WifiModule::~WifiModule(void)
{
}

//����ָ�Ƶ�ָ�������ƥ��
void WifiModule::WifiMach(vector<wifi> wifi_receive,string region_name)
{
	vector<vector<string> > wifilist=this->wifi_list_fp[region_name+"wifi_list"];
	RowVectorXd slice_match =RowVectorXd::Zero(wifilist.size());

	//��1������ά����
	for (unsigned int i = 0; i < wifi_receive.size(); i++)
	{
		for (unsigned int j = 0; j < wifilist.size(); j++)
		{
			if (wifi_receive[i].mac==wifilist[j][1])
			{
				slice_match[j]+=offset_rssi+wifi_receive[i].rssi;//��������������
				break;
			}
		}
	}
	slice_match/=slice_match.dot(slice_match);//��һ��
	//��2��ά���������ֱ�ӱȽ�
	//if(wholefp)
	vector<vector<double> > tmatrix=this->wifi_rssi_fp[region_name+"wifi"];

	vector<wifi_match_result> max_cos_x_y_vec;
	max_cos_x_y_vec.reserve(this->array_length);
	int min_index=-1;
	double vec_min=1;
	for (unsigned int rowindex = 0; rowindex < tmatrix.size(); rowindex++)
	{
		RowVectorXd thiswifi=RowVectorXd::Zero(tmatrix[rowindex].size()-2);
		for (int col_vec_index = 0; col_vec_index < thiswifi.cols(); col_vec_index++)
		{
			thiswifi(col_vec_index) = tmatrix[rowindex][col_vec_index];
		}
		double real_cos=slice_match.dot(thiswifi)/thiswifi.dot(thiswifi);
		if (real_cos<=0.4)
		{
			continue;
		}
		if( max_cos_x_y_vec.size()<this->array_length)//С�ڸ���
		{
			wifi_match_result twmr;
			twmr.probability=real_cos;
			twmr.x=tmatrix[rowindex][0];
			twmr.y=tmatrix[rowindex][1];
			twmr.region_id=region_name;
			max_cos_x_y_vec.push_back(twmr);
			if (real_cos<vec_min)//ˢ����С�ĺ�����
			{
				vec_min=real_cos;
				min_index=max_cos_x_y_vec.size();
			}
		}
		else
		{
			if (real_cos>vec_min)
			{
				wifi_match_result twmr;
				twmr.probability=real_cos;
				twmr.x=tmatrix[rowindex][0];
				twmr.y=tmatrix[rowindex][1];
				twmr.region_id=region_name;
				max_cos_x_y_vec.push_back(twmr);
				max_cos_x_y_vec[min_index]=twmr;
				//����ѭ���õ���Сֵ������
				min_index=0;
				vec_min=max_cos_x_y_vec[0].probability;
				for (unsigned int loop_min_index	= 1; loop_min_index < max_cos_x_y_vec.size(); loop_min_index++)
				{
					if (max_cos_x_y_vec[loop_min_index].probability<vec_min)
					{
						vec_min=max_cos_x_y_vec[loop_min_index].probability;
						min_index=loop_min_index;
					}
				}

			}
		}
		
	}
	//��3.�����ֵ���ڵ�һ��,��Сֵ�������һ�����Ա�ʹ�ú�ֵ�˳�
	int max_vec_index=0;
	for (unsigned int loop_max_index = 1; loop_max_index < max_cos_x_y_vec.size(); loop_max_index++)
	{
		if (max_cos_x_y_vec[loop_max_index].probability>max_cos_x_y_vec[max_vec_index].probability)
		{
			max_vec_index=loop_max_index;
		}
	}
	swap(max_cos_x_y_vec[max_vec_index],max_cos_x_y_vec[0]);
	swap(max_cos_x_y_vec[min_index],max_cos_x_y_vec[max_cos_x_y_vec.size()-1]);
	this->wifi_cal_list = max_cos_x_y_vec;
	
}
//���������ƥ�䣬�������µĲ��
void WifiModule::WifiMach(vector<wifi> wifi_receive )
{
	vector <wifi_match_result> wifiresult;
	vector<wifi_match_result> max_cos_x_y_vec;
	max_cos_x_y_vec.reserve(this->array_length);
	int min_index=-1;
	double vec_min=1;

	for (unsigned int i = 0; i < this->region_list.size(); i++)
	{
		string region_name=region_list[i];
		vector<vector<string> > wifilist=this->wifi_list_fp[region_name+"wifi_list"];
		RowVectorXd slice_match =RowVectorXd::Zero(wifilist.size());

		//��1������ά����
		for (unsigned int i = 0; i < wifi_receive.size(); i++)
		{
			for (unsigned int j = 0; j < wifilist.size(); j++)
			{
				if (wifi_receive[i].mac==wifilist[j][1])
				{
					slice_match[j]+=offset_rssi+wifi_receive[i].rssi;//��������������
					break;
				}
			}
		}
		slice_match/=slice_match.dot(slice_match);//��һ��
		//��2��ά���������ֱ�ӱȽ�
		//if(wholefp)
		vector<vector<double> > tmatrix=this->wifi_rssi_fp[region_name+"wifi"];



		for (unsigned int rowindex = 0; rowindex < tmatrix.size(); rowindex++)
		{
			RowVectorXd thiswifi=RowVectorXd::Zero(tmatrix[rowindex].size()-2);
			for (int col_vec_index = 0; col_vec_index < thiswifi.cols(); col_vec_index++)
			{
				thiswifi[col_vec_index] = tmatrix[rowindex][col_vec_index];
			}
			double real_cos=slice_match.dot(thiswifi)/thiswifi.dot(thiswifi);
			if (real_cos<=0.4)
			{
				continue;
			}
			if( max_cos_x_y_vec.size()<this->array_length)//С�ڸ���
			{
				wifi_match_result twmr;
				twmr.probability=real_cos;
				twmr.x=tmatrix[rowindex][0];
				twmr.y=tmatrix[rowindex][1];
				twmr.region_id=region_name;
				max_cos_x_y_vec.push_back(twmr);
				if (real_cos<vec_min)//ˢ����С�ĺ�����
				{
					vec_min=real_cos;
					min_index=max_cos_x_y_vec.size();
				}
			}
			else
			{
				if (real_cos>vec_min)
				{
					wifi_match_result twmr;
					twmr.probability=real_cos;
					twmr.x=tmatrix[rowindex][0];
					twmr.y=tmatrix[rowindex][1];
					twmr.region_id=region_name;
					max_cos_x_y_vec.push_back(twmr);
					max_cos_x_y_vec[min_index]=twmr;
					//����ѭ���õ���Сֵ������
					min_index=0;
					vec_min=max_cos_x_y_vec[0].probability;
					for (unsigned int loop_min_index	= 1; loop_min_index < max_cos_x_y_vec.size(); loop_min_index++)
					{
						if (max_cos_x_y_vec[loop_min_index].probability<vec_min)
						{
							vec_min=max_cos_x_y_vec[loop_min_index].probability;
							min_index=loop_min_index;
						}
					}

				}
			}
		
		}

	}
	//��3.�����ֵ���ڵ�һ��,��Сֵ�������һ�����Ա�ʹ�ú�ֵ�˳�
	int max_vec_index=0;
	for (unsigned int loop_max_index = 1; loop_max_index < max_cos_x_y_vec.size(); loop_max_index++)
	{
		if (max_cos_x_y_vec[loop_max_index].probability>max_cos_x_y_vec[max_vec_index].probability)
		{
			max_vec_index=loop_max_index;
		}
	}
	swap(max_cos_x_y_vec[max_vec_index],max_cos_x_y_vec[0]);
	swap(max_cos_x_y_vec[min_index],max_cos_x_y_vec[max_cos_x_y_vec.size()-1]);
	this->wifi_cal_list = max_cos_x_y_vec;
}

void WifiModule::vec2matrix(void)
{
}
