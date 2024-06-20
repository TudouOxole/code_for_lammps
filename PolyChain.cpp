#include "PolyChain.h"
void PolyChain::Initialize_Polychain(){//初始化读取构象文件
    this->ReadPolychain();
    this->N_obs_chain = (this->box.p2.x * 2/this->dTube) * (this->box.p2.y * 2/this->dTube);
    //cout << "N_obs_chain = "<<N_obs_chain <<endl;
    this->Set0();
    this->Unwrap();
    this->SortID();
    this->SortBonds();
    //提取测试链和柱子的信息
    for(int i = 0; i < this->N_atoms; i++){
        if(this->a[i].ID <= this->deg){
            this->TestAtoms.push_back(this->a[i]);
        }else{this->ObstaclesAtoms.push_back(this->a[i]);}
    }
    printf("Polymer data initializing success!\n");
}
void PolyChain::ReadPolychain(){//读取文件(针对lammps格式的构象文件)
    string file0 = this->path + this->name;
    string line;
    ifstream file;
    file.open(file0);
    if (!file.is_open()) {//如果文件打开失败直接退出
        std::cerr << "file open failed: " << file0 << "\n";
        exit(-1);
    }
    getline(file, line);// 第一行是注释，直接读取并忽略
    getline(file, line);//2
    file >> this->N_atoms;//3
    getline(file, line);//3
    file >> this->N_atomtypes;//4
    getline(file, line);//4
    file >> this->N_bonds;//5
    getline(file, line);//5
    file >> this->N_bondtypes;//6
    getline(file, line);//6
    getline(file, line);//7
    file >> this->box.p1.x >> this->box.p2.x;//8
    getline(file, line);//8
    file >> this->box.p1.y >> this->box.p2.y;//9
    getline(file, line);//9
    file >> this->box.p1.z >> this->box.p2.z;//10
    getline(file, line);//10
    getline(file, line);//11
    getline(file, line);//12
    getline(file, line);//13
    this->m.resize(this->N_atomtypes);//vector给容器重新定义大小，不然会栈溢出
    for (int i = 0; i < this->N_atomtypes; i++) {//14-17
        file >> this->m[i].type;
        file >> this->m[i].mass;
    }
    getline(file, line);//18
    getline(file, line);//19
    getline(file, line);//20
    this->a.resize(this->N_atoms);
    for (int i = 0; i < this->N_atoms; i++) {//21+
        file >> this->a[i].ID;
        file >> this->a[i].ChainID;
        file >> this->a[i].type;
        file >> this->a[i].p.x;
        file >> this->a[i].p.y;
        file >> this->a[i].p.z;
        file >> this->a[i].nx;
        file >> this->a[i].ny;
        file >> this->a[i].nz;
    }
    getline(file, line);
    getline(file, line);//Velocity
    getline(file, line);
    for (int i = 0; i < this->N_atoms; i++) {
        file >> this->a[i].ID;
        file >> this->a[i].v.x;
        file >> this->a[i].v.y;
        file >> this->a[i].v.z;
    }
    getline(file, line);
    getline(file, line);//Bonds
    getline(file, line);
    this->b.resize(this->N_bonds);
    for (int i = 0; i < this->N_bonds; i++) {
        file >> this->b[i].ID;
        file >> this->b[i].type;
        file >> this->b[i].AID1;
        file >> this->b[i].AID2;
    }
    printf("Polymer data reading success!\n");
}
void PolyChain::WritePolychain(){//输出文件(针对lammps格式的构象文件)
    ofstream ofs;
    string outfile = this->path + this->name;//输出文件名
    ofs.open(outfile, ios::out);
    ofs << "LAMMPS data file via obs.cpp ,by xcm" << endl;//1
    ofs << endl;//2
    ofs << this->N_atoms << " atoms" << endl;//3
    ofs << this->N_atomtypes << " atom types" << endl;//4
    ofs << this->N_bonds << " bonds" << endl;//5
    ofs << this->N_bondtypes << " bond types" << endl;//6
    ofs << endl;//7
    ofs << this->box.p1.x<<" "<<this->box.p2.x<<" xlo xhi" << endl;//8
    ofs << this->box.p1.y<<" "<<this->box.p2.y<<" ylo yhi" << endl;//9
    ofs << this->box.p1.z<<" "<<this->box.p2.z<<" zlo zhi" << endl;//10
    ofs << endl;//11
    ofs << "Masses" << endl;//12
    ofs << endl;//13
    for(int i = 0; i < this->N_atomtypes; i++){//14-17
        ofs << this->m[i].type<<" "<<this->m[i].mass << endl;
    }
    ofs << endl;//18
    ofs << "Atoms # bond" << endl;//19输出粒子ID、链ID、类型、坐标
    ofs << endl; //20
    for (int i = 0; i < this->N_atoms; i++) {//21+输出测试链，粒子ID 链ID 粒子类型 x y z nx ny nz
        ofs<< this->a[i].ID<<" "<<this->a[i].ChainID<<" "<<this->a[i].type<<" "<<this->a[i].p.x<<" "<<this->a[i].p.y<<" "<<this->a[i].p.z<<" "<<this->a[i].nx<<" "<<this->a[i].ny<<" "<<this->a[i].nz<<endl;
    }
    ofs << endl;
    ofs << "Velocities" << endl;
    ofs << endl;
    for (int i = 0; i < this->N_atoms; i++) {
        ofs<< this->a[i].ID<<" "<<this->a[i].v.x<<" "<<this->a[i].v.y<<" "<<this->a[i].v.z<<endl;
    }
    ofs << endl;
    ofs << "Bonds" << endl;
    ofs << endl;
    for (int i = 0; i < this->N_bonds; i++) {
        ofs<< i+1<<" "<<this->b[i].type<<" "<<this->b[i].AID1<<" "<<this->b[i].AID2<<endl;
    }
    printf("Polymer data Writing success!\n");
}
void PolyChain::Set0(){//把跑的太远的链拉回初始盒子
    int nx,ny,nz;
    for(int i = 0; i < this->N_atoms; i++){
        if (this->a[i].ChainID == 1){//链ID=1指测试链  
            nx = this->a[i].nx;
            ny = this->a[i].ny;
            nz = this->a[i].nz;
            break;
        }
    }
    for(int i = 0; i < this->N_atoms; i++){
        if (this->a[i].ChainID == 1){//链ID=1指测试链  
            this->a[i].nx = this->a[i].nx - nx;
            this->a[i].ny = this->a[i].ny - ny;
            this->a[i].nz = this->a[i].nz - nz;
        }
    }
    printf("Set0 Setting success!\n");
}
void PolyChain::Unwrap(){//周期坐标转换成真实坐标
    for(int i = 0; i < this->N_atoms; i++){
        if (this->a[i].ChainID == 1){//链ID=1指测试链  
            this->a[i].p.x += this->a[i].nx * 2 * this->box.p2.x;
            this->a[i].nx = 0;
            this->a[i].p.y += this->a[i].ny * 2 * this->box.p2.y;
            this->a[i].ny = 0;
            this->a[i].p.z += this->a[i].nz * 2 * this->box.p2.z;
            this->a[i].nz = 0;
        }
    }
    printf("Unwrap Setting success!\n");
}
void PolyChain::SortID(){//按粒子ID重新排序
    vector<Atoms> temp;
    temp.resize(this->N_atoms);
    for(int i = 0; i < this->N_atoms; i++){
        temp[this->a[i].ID-1] = this->a[i];
    }
    for(int i = 0; i < this->N_atoms; i++){
        this->a[i] = temp[i];
    }
    printf("SortID Setting success!\n");
}
void PolyChain::SortBonds(){//键排序
    for(int i = 0; i < (this->deg-1); i++){
        this->b[i].type = 1;
        this->b[i].AID1 = i+1;
        this->b[i].AID2 = i+2;
    }
    for(int i = 0; i < this->N_obs_chain; i++){
        for(int j = 0; j < (this->dTube-1); j++){
            int k = (this->deg-1) + (this->dTube-1) * i + j;
            this->b[k].type = 2;
            this->b[k].AID1 = k + 2 + i;
            this->b[k].AID2 = k + 3 + i; 
        }
    }
    printf("SortBonds Setting success!\n");
}
void PolyChain::Setunifor(double mid, double standard_deviation){//柱子非均匀化non-uniformity-obstacles
    //生成高斯分布随机数
    default_random_engine gen;
    normal_distribution<double> n(mid,standard_deviation);
    vector<double> p;
    int num = N_obs_chain * 2;
    double sum=0;
    double variance=0;
    for (int i = 0; i < num; i++){
        p.push_back(n(gen));
        if(p[i]>=(this->dTube/2)){
            p.pop_back();
            i--;
        }else if(p[i]<=(-this->dTube/2)){
            p.pop_back();
            i--;
        }
    }
    for (int i = 0; i < p.size(); i++){
        sum += p[i]/p.size();
    }
    cout << "sum = " << sum << endl;
    for (int i = 0; i < p.size();i++){
        variance += pow((p[i]-sum),2)/p.size();
    }
    cout << "variance = " << variance << endl;
    //加到柱子坐标上
    for(int i = 0; i < this->N_atoms; i++){
        if (this->a[i].ID > this->deg){//链ID=0指obs 
            int k = (this->a[i].ID - this->deg - 1)/this->dTube; 
            this->a[i].p.x += p[2*k]; 
            this->a[i].p.y += p[2*k+1]; 
        }
    }
    printf("Setunifor Setting success!\n");
}
void PolyChain::Setunifor2(double min, double max){
    //生成均匀分布随机数
    
    default_random_engine gen;
    uniform_real_distribution<double> n(min,max);
    vector<double> p;
    int num = this->N_obs_chain * 2;
    double sum=0;
    double variance=0;
    for (int i = 0; i < num; i++){
        p.push_back(n(gen));
    }
    for (int i = 0; i < p.size(); i++){
        sum += p[i]/p.size();
    }
    cout << "sum = " << sum << endl;
    for (int i = 0; i < p.size();i++){
        variance += pow((p[i]-sum),2)/p.size();
    }
    cout << "variance = " << variance << endl;
    //加到柱子坐标上
    for(int i = 0; i < this->N_atoms; i++){
        if (this->a[i].ID > this->deg){//链ID=0指obs 
            int k = (this->a[i].ID - this->deg - 1)/5; 
            this->a[i].p.x += p[2*k]; 
            this->a[i].p.y += p[2*k+1]; 
        }
    }
    printf("Setunifor2 Setting success!\n");
}
void PolyChain::Disassemble(int deg2){ // 拆链
    this->N_atoms = this->N_atoms + deg2 - this->deg;
    this->N_bonds = this->N_bonds + deg2 - this->deg;
    this->deg = deg2;
    this->a.clear();
    for (int i = 0; i < this->N_atoms; i++){
        if(i < deg2){
            this->a.push_back(this->TestAtoms[i]);
            if(i == (deg2-1)){
                this->a[i].type = 3;
            }
        }
        else{this->a.push_back(this->ObstaclesAtoms[i-deg2]);}
        this->a[i].ID = i+1;
    }
    this->SortBonds();
}
void PolyChain::SetMess(){   
    this->m.resize(this->N_atomtypes);//vector给容器重新定义大小，不然会栈溢出
    for (int i = 0; i < this->N_atomtypes; i++) {
        this->m[i].type = i + 1;
        this->m[i].mass = 1 ;
    }
}
void PolyChain::Set2DObs(double L_box,double del_mol){
    int obs_deg = this->dTube/del_mol + 1;             //柱子链长(3d:L_box;2d:this->dTube)
    double L_box_z = this->dTube;                      //z方向柱子长度(3d:this->obs_deg * del_mol;2d:this->dTube)
    int N_obs_x = (this->box.p2.x * 2/this->dTube);                 //x一行柱子数
    int N_obs_y = (this->box.p2.y * 2/this->dTube);                 //y一行柱子数
    double obs_del = del_mol;                       //柱子内粒子间距
    int N_atoms_obs = this->N_obs_chain * obs_deg;  //总粒子数（未计入高分子链的粒子，后面写的时候把高分子链的粒子加上去了）
    int N_atom_types_obs = 5;                       //总粒子类型（4,5为柱子,5为柱子末端; 1、2、3为高分子链; 1，3为高分子链的两末端）
    int N_bonds_obs = N_obs_x * N_obs_y*(obs_deg-1);//总键数（未计入高分子链的键，后面写的时候把高分子链的键加上去了）
    int N_bond_types_obs = 2;                       //总键类型(2为柱子; 1为高分子链)

    this->N_atoms = N_atoms_obs + this->deg;
    this->N_atomtypes = N_atom_types_obs;
    this->SetMess();
    this->N_bonds = N_bonds_obs + this->deg-1;
    this->N_bondtypes = N_bond_types_obs;

    this->a.resize(this->N_atoms);
    this->b.resize(this->N_bonds);
    for(int i = 0; i < this->N_obs_chain; i++){//柱子建模
        int x = i / N_obs_x;
        int y = i % N_obs_y;
        for(int j = 0; j < obs_deg; j++){
            int k = this->deg + j + i * obs_deg;
            this->a[k].ID = this->deg + i * obs_deg + j + 1;
            this->a[k].ChainID = 0;
            if (j == 0 || j ==(obs_deg-1)){
                this->a[k].type = 5;
            }else{
                this->a[k].type = 4;
            }
            this->a[k].p.x = ((1 + (2 * x)) * this->dTube / 2) - (L_box / 2) + 0.5;
            this->a[k].p.y = ((1 + (2 * y)) * this->dTube / 2) - (L_box / 2) + 0.5;
            this->a[k].p.z = j * obs_del - this->dTube / 2 + 0.06;
            this->a[k].nx = 0;
            this->a[k].ny = 0;
            this->a[k].nz = 0;
        }
    }
    for(int i = 0; i < (this->deg-1); i++){
        this->b[i].type = 1;
        this->b[i].AID1 = i+1;
        this->b[i].AID2 = i+2;
    }
    for(int i = 0; i < this->N_obs_chain; i++){
        for(int j = 0; j < (this->dTube); j++){
            int k = (this->deg-1) + (this->dTube) * i + j;
            this->b[k].type = 2;
            this->b[k].AID1 = k + 2 + i;
            this->b[k].AID2 = k + 3 + i; 
        }
    }
    printf("SortBonds Setting success!\n");
}

void PolyChain::mySet()
{
    int obs_deg = this->dTube + 1;             //柱子链长(3d:L_box;2d:this->dTube)
    int N_obs_x = (this->box.p2.x * 2/this->dTube);                 //x一行柱子数
    int N_obs_y = (this->box.p2.y * 2/this->dTube);                 //y一行柱子数
    int N_atoms_obs = this->N_obs_chain * obs_deg;  //总粒子数（未计入高分子链的粒子，后面写的时候把高分子链的粒子加上去了）
    int N_atom_types_obs = 5;                       //总粒子类型（4,5为柱子,5为柱子末端; 1、2、3为高分子链; 1，3为高分子链的两末端）
    int N_bonds_obs = N_obs_x * N_obs_y*(obs_deg-1);//总键数（未计入高分子链的键，后面写的时候把高分子链的键加上去了）
    int N_bond_types_obs = 2;                       //总键类型(2为柱子; 1为高分子链)

    this->N_atoms = N_atoms_obs + this->deg;
    this->N_atomtypes = N_atom_types_obs;
    this->SetMess();
    this->N_bonds = N_bonds_obs + this->deg-1;
    this->N_bondtypes = N_bond_types_obs;

    this->a.resize(this->N_atoms);
    this->b.resize(this->N_bonds);
    for(int i = 0; i < this->N_obs_chain; i++){//柱子建模
        for(int j = 0; j < obs_deg; j++){
            int k = this->deg + j + i * obs_deg;
            this->a[k].ID = this->deg + i * obs_deg + j + 1;
            this->a[k].ChainID = 0;
            if (j == 0 || j ==(obs_deg-1)){
                this->a[k].type = 5;
            }else{
                this->a[k].type = 4;
            }
            this->a[k].p.x = ObstaclesAtoms[(i+1)*4-1].p.x;
            this->a[k].p.y = ObstaclesAtoms[(i+1)*4-1].p.y;
            this->a[k].p.z = j * 0.97 - this->dTube / 2 + 0.06;
            this->a[k].nx = 0;
            this->a[k].ny = 0;
            this->a[k].nz = 0;
        }
    }
    for(int i = 0; i < (this->deg-1); i++){
        this->b[i].type = 1;
        this->b[i].AID1 = i+1;
        this->b[i].AID2 = i+2;
    }
    for(int i = 0; i < this->N_obs_chain; i++){
        for(int j = 0; j < (this->dTube); j++){
            int k = (this->deg-1) + (this->dTube) * i + j;
            this->b[k].type = 2;
            this->b[k].AID1 = k + 2 + i;
            this->b[k].AID2 = k + 3 + i; 
        }
    }
}

void PolyChain::moveChain(Point m){
    for (int i = 0; i < TestAtoms.size(); i++){
        this->TestAtoms[i].p += m;
        this->a[i].p = this->TestAtoms[i].p;
    }
}