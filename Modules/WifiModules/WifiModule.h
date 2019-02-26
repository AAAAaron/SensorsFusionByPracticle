#pragma once
#include "stdafx.h"

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

