
#include "Practicle_Filter.h"
#ifdef linux  
#include <unistd.h>  
#include <dirent.h>  
#endif  
#ifdef WIN32  
#include <direct.h>  
#include <io.h>  
#endif 

Practicle_Filter::Practicle_Filter(void)
{
	//使用流程
	//·1，载入指纹
	//load_fp_all("D:/指纹文件/30188888");
	//并将指纹分别赋给对应的对象
	// ·2新建一个粒子对象
	//InitPracClass(string prac_name)
	//·3初始粒子矩阵
	//main_practicle_init(prac_name,startx,starty,endx,endy,init_x,init_y,init_z, init_heading,practicle_num,practicle_dis);
	//·4选择不同的更新方式
	//main_progration()
	//·5必要时释放掉该粒子对象
	
}

Practicle_Filter::~Practicle_Filter(void)
{
}

//读取固定文件夹下所有的文件
vector<string> Practicle_Filter::getFiles(string cate_dir)  
{  
    vector<string> files;
#ifdef WIN32  
    _finddata_t file;  
    long lf;  
    //输入文件夹路径  
    if ((lf=_findfirst(cate_dir.c_str(), &file)) == -1) {  
        cout<<cate_dir<<" not found!!!"<<endl;  
    } else {  
        while(_findnext(lf, &file) == 0) {  
            //输出文件名  
            //cout<<file.name<<endl;  
            if (strcmp(file.name, ".") == 0 || strcmp(file.name, "..") == 0)  
                continue;  
            files.push_back(file.name);  
        }  
    }  
    _findclose(lf);  
#endif  

#ifdef linux  
    DIR *dir;  
    struct dirent *ptr;  


    if ((dir=opendir(cate_dir.c_str())) == NULL)  
        {  
	  cout<<"opendir is wrong"<<endl;
        }  

    while ((ptr=readdir(dir)) != NULL)  
    {  
        if(strcmp(ptr->d_name,".")==0 || strcmp(ptr->d_name,"..")==0)    ///current dir OR parrent dir  
                continue;  
        else{ 
	  switch(ptr->d_type){
	    case 8:///file  
	    {//printf("d_name:%s/%s\n",basePath,ptr->d_name); 
	      files.push_back(cate_dir+"/"+ptr->d_name); 
	      break;
	    }
	    case 10:    ///link file  
	    {  //printf("d_name:%s/%s\n",basePath,ptr->d_name);  
	      continue;
	      break;	  
	    }
	    case 4:    ///dir    
	    {  
	      vector<string> tempfiles;
	      tempfiles=getFiles(cate_dir+"/"+ptr->d_name);
	      for(unsigned int i=0;i<tempfiles.size();i++)
	      {
		files.push_back(tempfiles[i]); 
	      }

	      break;
	    }  
	    default:
	      continue;
	      
	  }
	  
	}
	     
         
    }  
    closedir(dir);  
#endif   
    return files;
} 

//载入给定文件夹中的所有指纹
void Practicle_Filter::load_fp_all(string filedir)
{
	vector<string> files;     
    //获取该路径下的所有文件  ，载入指纹
	files=getFiles(filedir); 
// 	for(unsigned int i=0;i<files.size();i++)
// 	{
// 	  cout<<files[i]<<endl;
// 	}
	for (unsigned int i = 0; i <files.size(); i++)
	{
		int startfilename =files[i].find_last_of("/");
		int re_second_sep = files[i].substr(0,startfilename).find_last_of("/");
		//int startfilename =files[i].find_last_of("/");//据说这样可以很好的在两个系统中运行
		int strindex = files[i].find_last_of(".");
		string FName =files[i].substr(re_second_sep+1,startfilename-re_second_sep-1);//目录中倒数第二个，F1等
		string sName =files[i].substr(startfilename+1,strindex-startfilename-1);//文件名

		if (region_list.size()>0)
		{
			vector<string>::iterator it;
			it=find(region_list.begin(),region_list.end(),FName);
			if (it!=region_list.end())
			{
				//找到
			}
			else
			{
				//未找到
				region_list.push_back(FName);
			}
		}
		else
		{
			region_list.push_back(FName);
		}
		if (sName=="wifi_list")
		{	
			ifstream inFile(files[i].c_str());  
			string lineStr;  
			vector<vector<string> > strArray;  
			while (getline(inFile, lineStr))  
			{  
				// 打印整行字符串  
				//cout << lineStr << endl;  
				// 存成二维表结构  
				stringstream ss(lineStr);  
				string str;  
				vector<string> lineArray;  
				// 按照逗号分隔  
				while (getline(ss, str, ','))  
					lineArray.push_back(str.c_str());  
				strArray.push_back(lineArray);  
			} 
			
			this->wifi_list_fp[FName+sName].swap(strArray);
			inFile.close();
		}
		else
		{
			
			ifstream inFile(files[i].c_str());  
			string lineStr;  
			vector<vector<double> > strArray;  
			while (getline(inFile, lineStr))  
			{  
				// 打印整行字符串  
				//cout << lineStr << endl;  
				// 存成二维表结构  
				stringstream ss(lineStr);  
				string str;  
				vector<double> lineArray;  
				// 按照逗号分隔  
				while (getline(ss, str, ','))  
					lineArray.push_back(atof(str.c_str()));  
				strArray.push_back(lineArray);  
			}
			if (sName=="wifi")
			{
				this->wifi_rssi_fp[FName+sName].swap(strArray);
			}
			else
			{
				this->all_fp[FName+sName].swap(strArray);
			}
			inFile.close();
			
		}
  

		
		//this->practicle_map[prac_name].load_fp(files[i],sName,FName);
	}

	this->wifi_obj.wifi_list_fp=this->wifi_list_fp;
	this->wifi_obj.wifi_rssi_fp=this->wifi_rssi_fp;
	test_loadfp();
}

void Practicle_Filter::test_loadfp()
{


	cout<<"test all fp data"<<endl;
	map<string,vector<vector<double> > >::iterator iter;  
	for(iter = this->all_fp.begin(); iter != this->all_fp.end(); iter++) { 
		cout<<iter->first<<' '<<iter->second.size()<<endl; 
	}
// 	cout<<this->all_fp.size()<<endl;
	cout<<"wifi fp------"<<endl;
	cout<<this->wifi_rssi_fp.size()<<endl;
	cout<<this->wifi_list_fp.size()<<endl;
 
	
	//cout<<"space has"<<this->space_fp.size()<<endl;
	//cout<<"mag has"<<this->mag1_fp.size()<<endl;
	return;
}

//希望可以完成多线程的调用
void Practicle_Filter::InitPracClass(string prac_name)
{
	Practicle Practicle_Filter;
	this->practicle_map[prac_name]=Practicle_Filter;
	this->practicle_map[prac_name].all_fp=this->all_fp;
	cout<<"sucess creat practicle class :"<<prac_name<<endl;
}

//直接使用wifi结果做初定位，获取该路径下的所有文件  
void Practicle_Filter::main_practicle_init(string prac_name,string floor_name,double heading,int practicle_num,int practicle_dis,vector<wifi_match_result> wifidata_list)
{
    map<string,vector<vector<double> > >::iterator iter; 
	iter = this->all_fp.find(floor_name+"xxyy");  
    if(iter != this->all_fp.end())  
	{
       //cout<<"Find "<<endl; 
		vector<vector<double> > xxyy=iter->second;
		int start_x = int(xxyy[0][0]);
		int start_y = int(xxyy[0][2]);
		int end_x = int(xxyy[0][1]);
		int end_y  = int(xxyy[0][3]);
		this->practicle_map[prac_name].Init_practicle(start_x,start_y,end_x,end_y,wifidata_list,heading,practicle_num,practicle_dis);
	}
    else  
	{
       cout<<"Do not Find"<<endl; 
	}
}
//使用已知位置初始化
void Practicle_Filter::main_practicle_init(string prac_name,string floor_name,double x,double y,double z,double heading,int practicle_num,int practicle_dis =5)
{
    map<string,vector<vector<double> > >::iterator iter; 
	iter = this->all_fp.find(floor_name+"xxyy");  
    if(iter != this->all_fp.end())  
	{
       //cout<<"Find "<<endl; 
		vector<vector<double> > xxyy=iter->second;
		int start_x = int(xxyy[0][0]);
		int start_y = int(xxyy[0][2]);
		int end_x = int(xxyy[0][1]);
		int end_y  = int(xxyy[0][3]);
		Point initpoint = {x,y,z};
		this->practicle_map[prac_name].Init_practicle(start_x,start_y,end_x,end_y,initpoint,heading,practicle_num,practicle_dis);
	}
    else  
	{
       cout<<"Do not Find"<<endl; 
	}
}
//粒子更新过程
void Practicle_Filter::main_practicle_progration(string prac_name,double step_length, double step_deta_angle, double step_miu=0 ,double step_sigma=0.5, double heading_miu =0,double heading_sigma=0.2)
{
	Noise_practicle step_noise ={"uniform",step_miu,step_sigma};
	Noise_practicle heading_noise ={"uniform",heading_miu,heading_sigma};
	this->practicle_map[prac_name].Practicle_progration(step_length, step_deta_angle,step_noise,heading_noise);
}

void Practicle_Filter::adjust_practicle_angle(string prac_name)
{
		double nxl[Point_List_length];
		double nyl[Point_List_length];
		for (int i = 0; i < Point_List_length; i++)
		{
			nxl[i]=x_list[i];
			nyl[i]=y_list[i];
		}
	    LeastSquare ls (Point_List_length,nxl, nyl);         
		if (ls.Rxy3>0.93 || ls.error_reg<0.01){
				this->practicle_map[prac_name].Init_practicle(ls.navi_angle);
				ls.print();
			}
		x_list.clear();
		y_list.clear();
}

//单独更新，均不扩散
Point Practicle_Filter::main_progration(string prac_name,string floor_id,double step_length,double step_deta_angle,double Heading)
{
	main_practicle_progration(prac_name,step_length, step_deta_angle,0,0,0,0);
// 	this->practicle_map[prac_name].adjust_weight_from_space(floor_id);
	this->practicle_map[prac_name].Rample_max2();
// 			// 处理修正角度
// 		if (x_list.size()%(Point_List_length)==0 && x_list.size()!=0){
// 			adjust_practicle_angle(prac_name);
// 		}
// 		if (step_length>0.2)
// 		{
// 			x_list.push_back( this->practicle_map[prac_name].weighted_position.x);
// 			y_list.push_back( this->practicle_map[prac_name].weighted_position.y);
// 		}
// 
// 		if (x_list.size()>Point_List_length)
// 		{
// 			x_list.erase(x_list.begin());
// 			y_list.erase(y_list.begin());
// 		}
	return this->practicle_map[prac_name].weighted_position;
	
	
}
//地磁，空间更新权值
Point Practicle_Filter::main_progration(string prac_name,int step_index,string floor_id, double step_length,double step_deta_angle,double Heading, double mag_norm,double mag_z,double mag_y,double inclination)
{
	
		main_practicle_progration(prac_name,step_length, step_deta_angle,0.1,sqrt(3/3.0),0,sqrt(2.0*12.0/180.0*M_PI));
		//main_practicle_progration(step_length, step_deta_angle,0.1,sqrt(1/3.0),0,sqrt(1));
		//this->practicle_map[prac_name].adjust_weight_from_space(floor_id);
		vector<double> mag_info_obs;
		mag_info_obs.push_back(mag_norm);
		mag_info_obs.push_back(mag_z);
		mag_info_obs.push_back(mag_y);
		mag_info_obs.push_back(inclination);
		this->practicle_map[prac_name].adjust_weight_mag_space(step_index,floor_id,mag_info_obs);

		this->practicle_map[prac_name].Rample_max2();
		//cout<<count<<","<<this->practicle_map[prac_name].practicle_diversity<<endl;
		if (this->practicle_map[prac_name].practicle_diversity ==0.0)
		{
			this->practicle_map[prac_name].Init_practicle(Heading);
		}
		//cout << this->practicle_map[prac_name].weighted_position.x<<","<<this->practicle_map[prac_name].weighted_position.y<<","<< endl;  


		 //处理修正角度
		if (x_list.size()%(Point_List_length)==0 && x_list.size()!=0){
			adjust_practicle_angle(prac_name);
		}
		if (step_length>0.2)
		{
			x_list.push_back( this->practicle_map[prac_name].weighted_position.x);
			y_list.push_back( this->practicle_map[prac_name].weighted_position.y);
		}
		if (x_list.size()>=Point_List_length)
		{
			x_list.erase(x_list.begin());
			y_list.erase(y_list.begin());
		}

		return this->practicle_map[prac_name].weighted_position;
}
// //这个用于蓝牙定位
// Point Practicle_Filter::main_progration(string prac_name,int step_index,string floor_id,double step_length,double step_deta_angle,double Heading, vector<wifi_match_result> bt_list)
// {
// 	main_practicle_progration(prac_name,step_length, step_deta_angle,0,sqrt(1/3.0),0,sqrt(2.0*12.0/180.0*M_PI));
// 	this->practicle_map[prac_name].adjust_weight_space_single_signal(floor_id, bt_list);
// 	this->practicle_map[prac_name].Rample_max2();
// 
// 		if (x_list.size()%(Point_List_length)==0 && x_list.size()!=0){
// 			adjust_practicle_angle(prac_name);
// 		}
// 		if (step_length>0.2)
// 		{
// 			x_list.push_back( this->practicle_map[prac_name].weighted_position.x);
// 			y_list.push_back( this->practicle_map[prac_name].weighted_position.y);
// 		}
// 		if (x_list.size()>=Point_List_length)
// 		{
// 			x_list.erase(x_list.begin());
// 			y_list.erase(y_list.begin());
// 		}
// 
// 	if (this->practicle_map[prac_name].practicle_diversity==0.0)
// 	{
// 		this->practicle_map[prac_name].Init_practicle(bt_list);
// 		return this->practicle_map[prac_name].max_weight_position;
// 	}
// 	else
// 	{
// 		return this->practicle_map[prac_name].weighted_position;
// 	}
// }
//wifi,空间更新权值和重采样
Point Practicle_Filter::main_progration(string prac_name,int step_index,string floor_id,double step_length,double step_deta_angle,double Heading, vector<wifi_match_result> wifidata_list)
{
		main_practicle_progration(prac_name,step_length, step_deta_angle,0,1,0,sqrt(2.0*12.0/180.0*M_PI));
		//main_practicle_progration(step_length, step_deta_angle,0,0,0,0);

		//this->practicle_map[prac_name].adjust_weight_mag_space_wifi(step_index,floor_id,mag_info_obs,wifidata_list);
		this->practicle_map[prac_name].adjust_weight_space_wifi(step_index,floor_id,wifidata_list);
		//this->practicle_map[prac_name].adjust_weight_from_space(floor_id);
		//this->practicle_map[prac_name].adjust_weight_mag_space(step_index,floor_id,mag_info_obs);
		//this->practicle_map[prac_name].weight_normalization(true);

		this->practicle_map[prac_name].Rample_max2();

		// 处理修正角度
		if (x_list.size()%(Point_List_length)==0 && x_list.size()!=0){
			adjust_practicle_angle(prac_name);
		}
		if (step_length>0.2)
		{
			x_list.push_back( this->practicle_map[prac_name].weighted_position.x);
			y_list.push_back( this->practicle_map[prac_name].weighted_position.y);
		}

		if (x_list.size()>Point_List_length)
		{
			x_list.erase(x_list.begin());
			y_list.erase(y_list.begin());
		}



		if (this->practicle_map[prac_name].practicle_diversity==0.0)
		{
			//this->practicle_map[prac_name].Init_practicle(startx,starty,endx,endy,wifidata_list,Heading,practicle_num,practicle_dis);
			this->practicle_map[prac_name].Init_practicle(wifidata_list);
			cout<<"toggle wifi init"<<endl;
// 			for(unsigned int i=1;i<wifidata_list.size();i++)
// 			{
// 			  cout<<wifidata_list[i].x<<"--"<<wifidata_list[i].y<<endl;
// 			}
			return this->practicle_map[prac_name].max_weight_position;
		}
		else
		{
			return this->practicle_map[prac_name].weighted_position;
		}
		

}
//wifi，地磁，空间更新权值和重采样
Point Practicle_Filter::main_progration(string prac_name,int step_index,string floor_id,double step_length,double step_deta_angle,double Heading, double mag_norm,double mag_z,double mag_y,double inclination,vector<wifi_match_result> wifidata_list)
{
		main_practicle_progration(prac_name,step_length, step_deta_angle,0,sqrt(1/3.0),0,sqrt(2.0*12.0/180.0*M_PI));
		//main_practicle_progration(step_length, step_deta_angle,0,0,0,0);
				vector<double> mag_info_obs;
		mag_info_obs.push_back(mag_norm);
		mag_info_obs.push_back(mag_z);
		mag_info_obs.push_back(mag_y);
		mag_info_obs.push_back(inclination);
		this->practicle_map[prac_name].adjust_weight_mag_space_wifi(step_index,floor_id,mag_info_obs,wifidata_list);
// 		this->practicle_map[prac_name].adjust_weight_space_wifi(step_index,floor_id,wifidata_list);
		//this->practicle_map[prac_name].adjust_weight_from_space(floor_id);
		//this->practicle_map[prac_name].adjust_weight_mag_space(step_index,floor_id,mag_info_obs);
		//this->practicle_map[prac_name].weight_normalization(true);

		this->practicle_map[prac_name].Rample_max2();
		 

		// 处理修正角度
		if (x_list.size()%(Point_List_length)==0 && x_list.size()!=0){
			adjust_practicle_angle(prac_name);
		}
		if (step_length>0.2)
		{
			x_list.push_back( this->practicle_map[prac_name].weighted_position.x);
			y_list.push_back( this->practicle_map[prac_name].weighted_position.y);
		}

		if (x_list.size()>Point_List_length)
		{
			x_list.erase(x_list.begin());
			y_list.erase(y_list.begin());
		}



		if (this->practicle_map[prac_name].practicle_diversity==0.0)
		{
			//this->practicle_map[prac_name].Init_practicle(startx,starty,endx,endy,wifidata_list,Heading,practicle_num,practicle_dis);
			this->practicle_map[prac_name].Init_practicle(wifidata_list);
			
			cout<<"toggle wifi init"<<wifidata_list.size()<<endl;
			return this->practicle_map[prac_name].max_weight_position;
		}
		else
		{
			return this->practicle_map[prac_name].weighted_position;
		}
		

}

//释放某个对象
void Practicle_Filter::main_free(string prac_name)
{
	int n = this->practicle_map.erase(prac_name);//如果删除了会返回1，否则返回0  
	if (n==1)
	{
		cout<<"delete object sucess!"<<endl;
	}
	else
	{
		cout<<"delete object fail!"<<endl;
	}
}
