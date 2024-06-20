#include "Data.h"
Data::Data(int Nfile,int deg,int N_chain,double deltaT,int Nstep,int N_folder,string data_in,string data_out){
	this->Nfile = Nfile;
	this->N_atoms = deg * N_chain;
	this->deg = deg;
	this->N_chain = N_chain;
	this->Nstep = Nstep;
	this->N_folder = N_folder;
	this->pathin = data_in;
	this->pathout = data_out;
	this->deltaT = deltaT;
	this->Read();
	printf("The %d folder's data was read.\n", (N_folder+1));
	this->calculateRcm();
	this->calculateDeltaRcm();
	this->calculateRg();
	this->calculateRe();
}

void Data::Read(){
	this->R.resize(this->Nfile);
	for (int N_time = 0; N_time < this->Nfile; N_time++) {
		R[N_time].resize(this->N_atoms);
		long long T = this->Nstep * N_time;
		int ID;
		this->fileName = this->pathin + "/N" + to_string(this->N_atoms/this->N_chain) + "/" + to_string(this->N_folder+1) + "/CAL_PO_" + to_string((long long)this->Nstep * N_time) + ".DAT";
		ifstream file(this->fileName);
		if (!file.is_open()) {
			throw runtime_error("file open failed: " + this->fileName + "\n");
		}
		string line;
		while (getline(file, line)){
			if (line == "ITEM: ATOMS id xu yu zu ") {
				for (int i = 0; i < this->N_atoms; i++) {
					file >> ID;
					file >> this->R[N_time][i].p.x;
					file >> this->R[N_time][i].p.y;
					file >> this->R[N_time][i].p.z;
				}
				break;
			}
		}
	}
}
void Data::Write(const vector<Point>dat,string Name){
	ofstream ofs;
    this->fileWrite = this->pathout + "/N" + to_string(this->N_atoms/this->N_chain) + "-" + Name + "-"+ to_string(this->Nstep) + ".txt";//输出文件名
    ofs.open(fileWrite , ios::out);
        ofs <<         "t"        << "    \t" <<                Name                << "    \t" <<         Name<<"xy"        << "    \t" <<  Name<<"z"  << endl;
    for (int i = 1; i < this->Nfile; i++) {
        ofs << this->Nstep * this->deltaT * i << "    \t" << (dat[i].x+dat[i].y+dat[i].z) << "    \t" << (dat[i].x+dat[i].y)/2 << "    \t" <<  dat[i].z << endl;
    }
    ofs.close();
}
void Data::Write(const vector<double>dat,string Name){
	ofstream ofs;
    this->fileWrite = this->pathout + "/N" + to_string(this->N_atoms/this->N_chain) + "-" + Name + "-"+ to_string(this->Nstep) + ".txt";//输出文件名
    ofs.open(fileWrite , ios::out);
        ofs <<         "t"        << "    \t" <<                Name                << "    \t" <<         Name<<"xy"        << "    \t" <<  Name<<"z"  << endl;
    for (int i = 0; i < Nfile; i++) {
        ofs << Nstep * deltaT * i << "    \t" << (this->ptp0xyz[i]) << "    \t" << this->ptp0xy[i] << "    \t" <<  this->ptp0[i].z << endl;
    }
    ofs.close();
}

void Data::calculateRcm(){
	this->Rcm.resize(this->Nfile);
	for (int N_time = 0; N_time < this->Nfile; N_time++) {
		Rcm[N_time].resize(this->N_chain);
		for (int i = 0; i < this->N_chain; i++) {
			for (int j = 0; j < this->deg; j++) {
				int k = i * this->deg + j;
				this->Rcm[N_time][i] += this->R[N_time][k].p / (this->deg*N_chain);
			}	
		}
	}
}

void Data::calculateRg(){
	this->Rg.resize(this->Nfile);
	this->Rg_chain.resize(this->Nfile);
	for (int N_time = 0; N_time < this->Nfile; N_time++) {
		Rg_chain[N_time].resize(this->N_chain);
		for (int i = 0; i < this->N_chain; i++) {
			for (int j = 0; j < this->deg; j++) {
				int k = i * this->deg + j;
					this->Rg_chain[N_time][i] += (this->R[N_time][k].p - this->Rcm[N_time][i]).scaleByPower(2) / this->deg;
			}
			this->Rg[N_time] += this->Rg_chain[N_time][i] / (this->N_chain);
		}
	}
}

void Data::calculateRe(){
	this->Re.resize(this->Nfile);
	this->Re_chain.resize(this->Nfile);
	for (int N_time = 0; N_time < this->Nfile; N_time++) {
		Re_chain[N_time].resize(this->N_chain);
		for (int i = 0; i < this->N_chain; i++) {
			int k = i * this->deg ;
			this->Re_chain[N_time][i] = (this->R[N_time][k].p - this->R[N_time][k + this->deg - 1].p);
			this->Re[N_time] += this->Re_chain[N_time][i].scaleByPower(2) / this->N_chain;
		}
	}
}

void Data::calculateDeltaRcm(){
	this->deltaRcm.resize(this->Nfile);
	for (int N_time = 0; N_time < this->Nfile; N_time++) {
		for (int deltat = 0; deltat < (this->Nfile - N_time); deltat++){
			if (this->N_chain != 1){ 
				for (int i = 0; i < this->N_chain; i++) {
					this->deltaRcm[N_time] += (this->Rcm[N_time+deltat][i] - this->Rcm[deltat][i]) / (this->N_chain*(this->Nfile - deltat));
				}
			}else{
				for (int i = 0; i < this->N_chain; i++) {
					this->deltaRcm[N_time] = 0;
				}
			}
		}
	}
}

void Data::calculateCMMSD(){
	this->g3.resize(this->Nfile);
	for (int N_time = 0; N_time < this->Nfile; N_time++) {
		for (int deltat = 0; deltat < (this->Nfile - N_time); deltat++){
			for (int i = 0; i < this->N_chain; i++) {
				this->g3[N_time] += (this->Rcm[N_time+deltat][i] - this->Rcm[deltat][i] - this->deltaRcm[N_time]).scaleByPower(2) /((this->Nfile - N_time) * this->N_chain);
			}
		}
	}
	printf("CMMSD completed...\n");
}

void Data::calculateMONOMSD(){
	this->calculateDeltaRcm();
	this->g1.resize(this->Nfile);
	int mono = this->deg / 5 ;
	for (int N_time = 0; N_time < this->Nfile; N_time++) {
		for (int deltat = 0; deltat < (this->Nfile - N_time); deltat++){
			for (int i = 0; i < this->N_chain; i++) {
				for (int j = 0; j < mono ; j++) {
					int k = i * this->deg + j + (this->deg - mono)/2;
					this->g1[N_time] += (this->R[N_time+deltat][k].p - this->R[deltat][k].p - this->deltaRcm[N_time]).scaleByPower(2) / ((this->Nfile - N_time)* this->N_chain * mono);
				}
			}
		}
	}
	printf("MONOMSD completed...\n");
}

void Data::calculateCt(){
	this->ptp0.resize(this->Nfile);
	this->ptp0xy.resize(this->Nfile);
	this->ptp0xyz.resize(this->Nfile);
	for (int N_time = 0; N_time < this->Nfile; N_time++) {
		for (int deltat = 0; deltat < (this->Nfile - N_time); deltat++){
			for (int i = 0; i < this->N_chain; i++) {
				this->ptp0xyz[N_time] += (this->Re_chain[N_time+deltat][i].x * this->Re_chain[deltat][i].x
											+ this->Re_chain[N_time+deltat][i].y * this->Re_chain[deltat][i].y
											+ this->Re_chain[N_time+deltat][i].z * this->Re_chain[deltat][i].z)
											/(this->N_chain * (this->Nfile - N_time)
											* (this->Re[deltat].x+this->Re[deltat].y+this->Re[deltat].z));

				this->ptp0xy[N_time] += (this->Re_chain[N_time+deltat][i].x * this->Re_chain[deltat][i].x 
											+ this->Re_chain[N_time+deltat][i].y * this->Re_chain[deltat][i].y)
											/(this->N_chain * (this->Nfile - N_time) * (this->Re[deltat].x + this->Re[deltat].y));

				this->ptp0[N_time] += (this->Re_chain[N_time+deltat][i] * this->Re_chain[deltat][i])/(this->Re[deltat] * this->N_chain * (this->Nfile - N_time));
			}
		}
	}
	printf("C(t) completed...\n");
}


void Data::debugWrite(vector<vector<Point>> dat,string Name)
{	
	ofstream ofs;
    this->fileWrite = this->pathout + "/N" + to_string(this->N_atoms/this->N_chain) + "-" + Name + "-"+ to_string(this->Nstep) + ".txt";//输出文件名
    ofs.open(fileWrite , ios::out);
    ofs <<"t"<< "\t";
	for (int i = 0; i < dat.size(); i++)
	{
		ofs <<Name<<to_string(i+1)<< "\t";
	}
	ofs << endl;
    double val,val1;
	for (int i = 1; i < this->Nfile; i++) {
		val1 = this->Nstep * this->deltaT * i;
		ofs << val1  << "\t";
		for (int j = 0; j < dat.size(); j++)
		{
			val = dat[j][i].addPoint();
			ofs << val  << "\t";
		}
        ofs << endl;
    }
    ofs.close();
}








