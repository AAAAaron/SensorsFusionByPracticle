// stdafx.h : ��׼ϵͳ�����ļ��İ����ļ���
// ���Ǿ���ʹ�õ��������ĵ�
// �ض�����Ŀ�İ����ļ�
//
// TODO: �ڴ˴����ó�����Ҫ������ͷ�ļ�
extern "C" 
#pragma once

#define M_PI 3.14159265358979323846
 
#include <iostream>  

using namespace std;
#include <Eigen/Dense> 
using namespace Eigen;

// #include <iosfwd>
#include <fstream>
#include <sstream>

#include <time.h>
#include <string>
#include <math.h> 


struct Point{      //����һ���ṹ������Point 
 double x;         
 double y;   
 double z;
 string region_id;
};

struct Noise_practicle{
	string mode;
	double miu;
	double sigma;
};
struct signal_match_result:Point
{

	double probability;

};

struct wifi_match_result:signal_match_result
{
	double promise_length;

};

struct RFsignal
{
	double time;
	string ssid;
	string mac;
	int rssi;
};

struct ble:RFsignal
{
	string Major;
	string Minor;

};

struct wifi:RFsignal
{

};