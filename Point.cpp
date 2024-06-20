#include "Point.h"
Point Point::operator+(const Point& other) const {// 加法运算符
    return Point(this->x + other.x, this->y + other.y, this->z + other.z);  
}
Point& Point::operator+=(const Point& other) {// 加法赋值运算符
    this->x += other.x;  
    this->y += other.y;  
    this->z += other.z;  
    return *this;  
}
Point Point::operator-(const Point& other) const {// 减法运算符
    return Point(this->x - other.x, this->y - other.y, this->z - other.z);  
}
Point& Point::operator-=(const Point& other) {// 减法赋值运算符
    this->x -= other.x;  
    this->y -= other.y;  
    this->z -= other.z;  
    return *this;  
}
Point Point::operator*(double scalar) const {// 标量乘法运算符
    return Point(this->x * scalar, this->y * scalar, this->z * scalar);  
}
Point Point::operator*(int scalar) const {// 标量乘法运算符
    return Point(this->x * scalar, this->y * scalar, this->z * scalar);  
}
Point Point::operator*(const Point& scalar) const {// 标量乘法运算符
    return Point(this->x * scalar.x, this->y * scalar.y, this->z * scalar.z);  
}
Point& Point::operator*=(double scalar) {// 标量乘法赋值运算符
    this->x *= scalar;  
    this->y *= scalar;  
    this->z *= scalar;  
    return *this;  
}
Point Point::operator/(double scalar) const {// 标量除法运算符
    if (scalar == 0) {  
        throw invalid_argument("Division by zero is not allowed.");  
    }  
    return Point(this->x / scalar, this->y / scalar, this->z / scalar);  
}
Point Point::operator/(const Point& scalar) const {// 标量除法运算符
    if (this->x==0||this->y==0||this->z==0) {  
        throw invalid_argument("Division by zero is not allowed.");  
    }  
    return Point(this->x / scalar.x, this->y / scalar.y, this->z / scalar.z);  
}
Point& Point::operator/=(double scalar) {// 标量除法赋值运算符
    if (scalar == 0) {  
        throw invalid_argument("Division by zero is not allowed.");  
    }  
    this->x /= scalar;  
    this->y /= scalar;  
    this->z /= scalar;  
    return *this;  
}
Point Point::scaleByPower(double exponent) const {//分别幂运算
    return Point(pow(this->x, exponent), pow(this->y, exponent), pow(this->z, exponent));  
}
double Point::addPoint(){
    return this->x+this->y+this->z;
}
