// stdafx.h : 标准系统包含文件的包含文件，
// 或是经常使用但不常更改的
// 特定于项目的包含文件
//
// TODO: 在此处引用程序需要的其他头文件
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


struct Point{      //声明一个结构体类型Point 
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