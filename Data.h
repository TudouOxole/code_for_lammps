#include <string>
#include <vector>
#include <fstream>
#include <algorithm>
#include "Atom.h"
class Data{
    public:

    Data(int Nfile,int deg,int N_chain,double deltaT,int Nstep,int N_folder,string data_in,string data_out);

    void Read();
    void Write(const vector<Point>dat,string Name);
    void Write(const vector<double>dat,string Name);
    void calculateRcm();
    void calculateDeltaRcm();
    void calculateRg();
    void calculateRe();
    void calculateCMMSD();
    void calculateMONOMSD();
    void calculateCt();
    void debugWrite(vector<vector<Point>> dat,string Name);

    int Nfile;
    int N_folder,deg,N_atoms,N_chain,Nstep;
    double deltaT; 
    string pathin,pathout,fileName,fileWrite;
    vector<vector<Atoms>> R;//粒子信息
    vector<vector<Point>> Rcm,Rg_chain,Re_chain;
    vector<Point> deltaRcm;
    vector<Point> Rg,Re;
    vector<Point> g1,g3;
    vector<Point> ptp0;
    vector<double> ptp0xy,ptp0xyz;
    
};