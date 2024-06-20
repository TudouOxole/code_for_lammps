//构象文件处理
#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <vector>
#include <iomanip>
#include <stdexcept>
#include <cstdio>
#include "PolyChain.h"
using namespace std;
using std::vector;
const int folder = 20;//平行样本数
vector<int> deg = {300};//链长 
const double del_mol = 0.97 ;//粒子间距
const int del_Tube = 4 ;//柱子的间距
const double L_box = 120 ;//盒子边长
const string path_in = "D:\\data\\non_uniformity_obstacles\\buquanlizi\\origin_form\\N";//输入文件位置
const string file_in = "poly.chain";//输入文件名
const string path_out = "D:\\data\\non_uniformity_obstacles\\buquanlizi\\form\\N";//输出文件位置(需提前创建好文件夹)
const string file_out = "poly.chain";//输出文件名
void test01(){
    PolyChain polyin;
    for(int degree = 0; degree < deg.size(); degree++){
        printf("The N%d system was setting...\n", deg[degree]);
        for (int N_folder = 0; N_folder < folder; N_folder++){
            //读取文件(针对lammps格式的构象文件)
            polyin.deg = deg[degree];
            polyin.path = path_in + to_string(polyin.deg) + "/" + to_string(N_folder + 1) + "/";
            polyin.name = file_in;//输入文件名
            polyin.dTube = del_Tube;
            polyin.Initialize_Polychain();
            printf("The %d folder's data was read.\n", (N_folder+1));
            PolyChain polyout(polyin);

            //修改区
            polyout.moveChain(Point(1,0,0));
            printf("Setting success!\n");

            //输出文件(针对lammps格式的构象文件)
            polyout.path = path_out + to_string(deg[degree]) + "/" + to_string(N_folder + 1) + "/";
            polyout.name = file_out;//输出文件名
            polyout.WritePolychain();
            printf("Writing complete!\n");
        }
    }
}
int main(){
    test01();
}