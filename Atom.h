#include "Point.h"
class Atoms{
    public:
    int ID,ChainID,type;//粒子ID、链ID、粒子类型
    Point p,v;//粒子坐标、速度
    int nx,ny,nz;//粒子所在的盒子周期
};