#pragma once//防止头文件重复包含
#include <iostream>  
#include <cmath>  
#include <stdexcept>  
using namespace std;
class Point {  
public:  
    // 构造函数  
    Point(double x = 0, double y = 0, double z = 0) : x(x), y(y), z(z) {}  
    Point operator+(const Point& other) const;  //加法运算符
    Point& operator+=(const Point& other);      //加法赋值运算符
    Point operator-(const Point& other) const;  //减法运算符
    Point& operator-=(const Point& other);      //减法赋值运算符
    Point operator*(double scalar) const;       //标量乘法运算符
    Point operator*(int scalar) const;
    Point operator*(const Point& scalar) const;
    Point& operator*=(double scalar);           //标量乘法赋值运算符
    Point operator/(double scalar) const;       //标量除法运算符
    Point operator/(const Point& scalar) const;
    Point& operator/=(double scalar);           //标量除法赋值运算符
    Point scaleByPower(double exponent) const;
    double addPoint();
    // 成员变量 
    double x, y, z;  
    // 输出Point对象  
    friend ostream& operator<<(ostream& os, const Point& p) {  
        os << "(" << p.x << ", " << p.y << ", " << p.z << ")";  
        return os;  
    }  
};  
class Vector : public Point {  
public:  
    using Point::Point; // 继承Point的构造函数  
    // 点乘运算符（需要另一个Vector对象）  
    double operator*(const Vector& other) const {  
        return this->x * other.x + this->y * other.y + this->z * other.z;  
    }  
    // 叉乘（cross product）

    // 长度（magnitude）或模（norm）

};  