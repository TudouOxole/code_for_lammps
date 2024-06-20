#pragma once//防止头文件重复包含
#include <iostream>//输入输出流头文件
#include <string>
#include <vector>
#include <fstream>
#include <random>
#include <cmath>
#include <algorithm>
#include "Point.h"
#include "Atom.h"
using namespace std;
class Box{//盒子的两个点
    public:
    Point p1,p2;
};
class Masses{//粒子类型、质量
    public:
    int type,mass;
};
class Bonds{//键ID、键类型、键链粒子1、2
    public:
    int ID,type,AID1,AID2;
};
class PolyChain{//构象文件的信息
    public:
    void Initialize_Polychain();
    void ReadPolychain();
    void WritePolychain();
    void Set0();//把跑的太远的链拉回初始盒子
    void Unwrap();//周期坐标转换成真实坐标
    void SortID();//按粒子ID重新排序
    void SortBonds();//键排序
    void Setunifor(double mid , double standard_deviation);
    void Setunifor2(double min , double max);
    void Disassemble(int deg2);//拆链
    void SetMess();
    void Set2DObs(double L_box, double del_mol);
    void mySet();
    void moveChain(Point m);

    string path,name;//文件位置,名
    int N_atoms,N_atomtypes,N_bonds,N_bondtypes;//粒子数、粒子类型、键数、键类型
    int N_obs_chain;//柱子数
    int deg;//体系链长
    double dTube;
    Box box;//盒子信息
    vector<Atoms> a;//粒子信息
    vector<Atoms> TestAtoms;//测试链粒子信息
    vector<Atoms> ObstaclesAtoms;//柱子粒子信息
    vector<Masses> m;//质量信息
    vector<Bonds> b;//键信息
};