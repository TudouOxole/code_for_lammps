#include "Data.h"

//编辑区
const int folder = 9;
const int Nfile = 1601;
vector<int> deg = {500};//链长
const int N_chain = 1;
const double deltaT = 0.01;//定义脚本文件步长
const int Nstep = 50;//定义Nstep
const string data_in = "D:\\data\\non_uniformity_obstacles\\test\\calculate2";
const string data_out = data_in + "/data3";
string str0 = data_in + "/N";
string str1 = "/";
string str2 = "/CAL_PO_";
string str3 = ".DAT";

void Write(vector<vector<Point>> dat, int deg, string Name)//输出每一个构象的数据
{	
	ofstream ofs;
    string file = data_out + "/N" + to_string(deg) + "-" + Name + "-"+ to_string(Nstep) + ".txt";//输出文件名
    ofs.open(file , ios::out);
    ofs <<"t"<< "\t";
	for (int i = 0; i < dat.size(); i++)
	{
		ofs <<Name<<to_string(i+1)<< "\t";
	}
	ofs << endl;
    double val,val1;
	for (int i = 1; i < dat[0].size(); i++) {
		val1 = Nstep * deltaT * i;
		ofs << val1  << "\t";
		for (int j = 0; j < dat.size(); j++)
		{
			val = (dat[j][i].addPoint());
			ofs << val  << "\t";
		}
        ofs << endl;
    }
    ofs.close();
} 
void Write(vector<double> dat, int deg, string Name)//输出所有构象平均的数据
{	
	ofstream ofs;
    string file = data_out + "/N" + to_string(deg) + "-" + Name + "-"+ to_string(Nstep) + ".txt";//输出文件名
    ofs.open(file , ios::out);
    ofs <<"t"<< "\t"<<Name<< "\t"<< endl;
    double val1,val2;
	for (int i = 1; i < dat.size(); i++) {
		val1 = Nstep * deltaT * i;
        val2 = dat[i];
		ofs << val1  << "\t"<< val2  << "\t"<< endl;
    }
    ofs.close();
} 

void test01(){
    clock_t start_time = clock();
    for(int degree = 0; degree < deg.size(); degree++){
        printf("The N%d system was calculating...\n", deg[degree]);
        vector<vector<Point>> data;
        vector<double> aveDat;
        aveDat.resize(Nfile);
        for (int N_folder = 0; N_folder < folder; N_folder++){
            Data *d = new Data(Nfile,deg[degree],N_chain,deltaT,Nstep,N_folder,data_in,data_out);
            d->calculateCMMSD();
            data.push_back(d->g3);
            delete d;
        }
        for (int i = 0; i < data.size(); i++) {
            for (int j = 0; j < data[i].size(); j++){
                aveDat[j] += data[i][j].addPoint();
            }
        }
        printf("Writing data...\n");
        Write(aveDat,deg[degree],"CMMSD2");

    }
    clock_t end_time = clock();
    printf("Finished!!! Using %d s!\n",(double)((end_time - start_time) / CLOCKS_PER_SEC));
}
int main(){
    test01();
    return 0;
}