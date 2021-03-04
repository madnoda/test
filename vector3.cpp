/*

  Copyright 2021、野田篤司、Twitter: @madnoda 、All right reserved.

  このプログラムは、MITライセンスとする

  プログラムに対する改変を加えた場合でも、加えない場合でも、ソース・
  コードとバイナリーのいずれの形式における再配布や使用は、下記の条件を
  満たした場合にのみ許可される。
  (条件 1)
  ソース・コードの再配布には、必ず上記の著作権表示と、この条件リスト、
  および以下の免責事項を改変せずにファイルの最初の行に入れなければなら
  ない。
  (条件 2)
  バイナリ−形式での再配布には、必ず上記の著作権表示と、この条件リスト
  および以下の免責事項をドキュメント中、もしくは配布と共にアーカイブ
  などのパッケージに入れて、提供しなければならない。

  (免責事項)
  このソフトウエアは、著作者により『現状どおり』で提供するものである。
  ソフトウエアに市場性があり、『売れる製品』であるような含みがあったり、
  『一定の目的に適うような製品』であるような、いかなる明示的な保証は
  無い。また、何らかの保証があるような含みのある一切の保証の類からの責任
  も全て拒否する。
  著作者は、直接的・間接的・偶発的・突発的あるいは必然的等、いかなる
  場合においてもいかなる損害に一切の責任は負わない。その損害には、この
  ソフトウエアの代替となりうるソフトウエアやサービスへの調達、効果、
  情報や利益の低下、商取引自体の損害を含む。また、いかなる倫理観に
  おいても、契約に関するものであろうと無かろうと、いかなる責務において
  も、不正行為があろうが無かろうと、このソフトウエアを使用することによっ
  て生じるいかなる損害の、一切の責任から著作者は免れる。たとえ、それが
  事前に当該の損害への可能性が示唆されていた場合においても同様である。
*/

#include	<string.h>
#include	<stdio.h>
#include	<stdlib.h>
#include	<math.h>
#include	"vector3.h"

double Vector3::abs() const
{
  double w = x * x + y * y + z * z;
  if (w > 0.0 ) {
    return sqrt(w);
  } else {
    return 0.0;
  }
}

Vector3 Vector3::operator-() const
{
  Vector3 r;
  r.x = -x;
  r.y = -y;
  r.z = -z;
  return r;
}

Vector3 Vector3::operator+=(const Vector3& a)
{
  x += a.x;
  y += a.y;
  z += a.z;
  return *this;
}
Vector3 Vector3::operator-=(const Vector3& a)
{
  x -= a.x;
  y -= a.y;
  z -= a.z;
  return *this;
}
Vector3 Vector3::operator*=(double w)
{
  x *= w;
  y *= w;
  z *= w;
  return *this;
}
Vector3 Vector3::operator/=(double w)
{
  x /= w;
  y /= w;
  z /= w;
  return *this;
}

Vector3 Vector3::unit()
{
  double w;
  w = abs();
  x /= w;
  y /= w;
  z /= w;
  return *this;
}

Vector3 operator^(const Vector3&a,const Vector3&b)
{
  Vector3 r;
  r.x = a.y * b.z - a.z * b.y;
  r.y = a.z * b.x - a.x * b.z;
  r.z = a.x * b.y - a.y * b.x;
  return r;
}

double operator*(const Vector3& a,const Vector3&b)
{
  return a.x * b.x + a.y * b.y + a.z * b.z;
}

Vector3 operator+(const Vector3&a,const Vector3&b)
{
  Vector3 r;
  r.x = a.x + b.x;
  r.y = a.y + b.y;
  r.z = a.z + b.z;
  return r;
}

Vector3 operator-(const Vector3&a,const Vector3&b)
{
  Vector3 r;
  r.x = a.x - b.x;
  r.y = a.y - b.y;
  r.z = a.z - b.z;
  return r;
}

Vector3 operator*(double a,const Vector3& b)
{
  Vector3 r;
  r.x = a * b.x;
  r.y = a * b.y;
  r.z = a * b.z;
  return r;
}

Vector3 operator/(const Vector3& b,double a)
{
  Vector3 r;
  r.x = b.x / a;
  r.y = b.y / a;
  r.z = b.z / a;
  return r;
}

double & Matrix33::operator()(int i,int j)
{
  if ((i == 1) && (j == 1)) return e11;
  if ((i == 1) && (j == 2)) return e12;
  if ((i == 1) && (j == 3)) return e13;
  if ((i == 2) && (j == 1)) return e21;
  if ((i == 2) && (j == 2)) return e22;
  if ((i == 2) && (j == 3)) return e23;
  if ((i == 3) && (j == 1)) return e31;
  if ((i == 3) && (j == 2)) return e32;
  if ((i == 3) && (j == 3)) return e33;
  return e11;
}

Matrix33 Matrix33::operator*=(double w)
{
  e11 *= w;
  e12 *= w;
  e13 *= w;
  e21 *= w;
  e22 *= w;
  e23 *= w;
  e31 *= w;
  e32 *= w;
  e33 *= w;
  return *this;
}

Vector3 Matrix33::GetX() const
{
  Vector3 r;
  r.x = e11;
  r.y = e21;
  r.z = e31;
  return r;
}

Vector3 Matrix33::GetY() const
{
  Vector3 r;
  r.x = e12;
  r.y = e22;
  r.z = e32;
  return r;
}

Vector3 Matrix33::GetZ() const
{
  Vector3 r;
  r.x = e13;
  r.y = e23;
  r.z = e33;
  return r;
}

void Matrix33::SetXYZ(const Vector3&ix,const Vector3&iy,const Vector3&iz)
{
  e11 = ix.x;
  e21 = ix.y;
  e31 = ix.z;
  e12 = iy.x;
  e22 = iy.y;
  e32 = iy.z;
  e13 = iz.x;
  e23 = iz.y;
  e33 = iz.z;
}
Matrix33 Matrix33::t() const
{
  Matrix33 r;
  r.e11 = e11;
  r.e12 = e21;
  r.e13 = e31;
  r.e21 = e12;
  r.e22 = e22;
  r.e23 = e32;
  r.e31 = e13;
  r.e32 = e23;
  r.e33 = e33;
  return r;
}

Matrix33 operator*(const Matrix33 & a,const Matrix33& b)
{
  Matrix33 r;
  r.e11 = a.e11 * b.e11 + a.e12 * b.e21 + a.e13 * b.e31;
  r.e12 = a.e11 * b.e12 + a.e12 * b.e22 + a.e13 * b.e32;
  r.e13 = a.e11 * b.e13 + a.e12 * b.e23 + a.e13 * b.e33;
  r.e21 = a.e21 * b.e11 + a.e22 * b.e21 + a.e23 * b.e31;
  r.e22 = a.e21 * b.e12 + a.e22 * b.e22 + a.e23 * b.e32;
  r.e23 = a.e21 * b.e13 + a.e22 * b.e23 + a.e23 * b.e33;
  r.e31 = a.e31 * b.e11 + a.e32 * b.e21 + a.e33 * b.e31;
  r.e32 = a.e31 * b.e12 + a.e32 * b.e22 + a.e33 * b.e32;
  r.e33 = a.e31 * b.e13 + a.e32 * b.e23 + a.e33 * b.e33;
  return r;
}

Vector3 operator*(const Matrix33 & a,const Vector3& b)
{
  Vector3 r;
  r.x = a.e11 * b.x + a.e12 * b.y + a.e13 * b.z;
  r.y = a.e21 * b.x + a.e22 * b.y + a.e23 * b.z;
  r.z = a.e31 * b.x + a.e32 * b.y + a.e33 * b.z;
  return r;
}

Quaternion::Quaternion(const Vector3& v)
{
  double a = sqrt(v.x * v.x + v.y * v.y + v.z * v.z);
  q0 = cos(a / 2.0);
  if (a != 0.0) {
    q1 = v.x * sin(a / 2.0) / a;
    q2 = v.y * sin(a / 2.0) / a;
    q3 = v.z * sin(a / 2.0) / a;
  } else {
    q1 = q2 = q3 = 0.0;
  }
}

Quaternion Quaternion::operator =(const Quaternion & s)
{
  q0 = s.q0;
  q1 = s.q1;
  q2 = s.q2;
  q3 = s.q3;
  kikakuka();
  return *this;
}

Quaternion Quaternion::operator =(const Matrix33 & s)
{
  const double e = 0.001;
  q0 = sqrt(  s.e11 + s.e22 + s.e33 + 1.0) / 2.0;
  if (q0 > e) {
    q1 = (s.e32 - s.e23)/(4.0 * q0);
    q2 = (s.e13 - s.e31)/(4.0 * q0);
    q3 = (s.e21 - s.e12)/(4.0 * q0);
  } else {
    q0 = 0.0;
    q1 = sqrt(  s.e11 - s.e22 - s.e33 + 1.0) / 2.0;
    if (q1 > e) {
      q2 = (s.e12 + s.e21)/(4.0 * q1);
      q3 = (s.e13 + s.e31)/(4.0 * q1);
    } else {
      q1 = 0.0;
      q2 = sqrt(- s.e11 + s.e22 - s.e33 + 1.0) / 2.0;
      if (q2 > e) {
	q3 = (s.e23 + s.e32)/(4.0 * q2);
      } else {
	q2 = 0.0;
	q3 = sqrt(- s.e11 - s.e22 + s.e33 + 1.0) / 2.0;
      }
    }
  }
  kikakuka();
  return *this;
}

void Quaternion::kikakuka()
{
  if (q0 < 0.0) {
    q0 = - q0;
    q1 = - q1;
    q2 = - q2;
    q3 = - q3;
  }
  double a = sqrt(q0 *q0 + q1 * q1 + q2 * q2 + q3 * q3);
  q0 /= a;
  q1 /= a;
  q2 /= a;
  q3 /= a;
}

Quaternion Quaternion::operator -()
{
  Quaternion q;
  q.q0 = q0;
  q.q1 = - q1;
  q.q2 = - q2;
  q.q3 = - q3;
  return q;
}
void Quaternion::GetV(double &a,double &x,double &y,double &z)
{
  double si = sqrt(q1 * q1 + q2 * q2 + q3 * q3);
  if (si == 0.0) {
    a = y = z = 0.0;
    x = 1.0;

  } else {
    a = 2.0 * atan2(si,q0);
    x = q1 / si;
    y = q2 / si;
    z = q3 / si;
  }
}

Quaternion operator*(const Quaternion& n,const Quaternion& q)
{
  Quaternion r;
  r.q0 = n.q0 * q.q0 - n.q1 * q.q1 - n.q2 * q.q2 - n.q3 * q.q3;
  r.q1 = n.q0 * q.q1 + n.q1 * q.q0 + n.q2 * q.q3 - n.q3 * q.q2;
  r.q2 = n.q0 * q.q2 - n.q1 * q.q3 + n.q2 * q.q0 + n.q3 * q.q1;
  r.q3 = n.q0 * q.q3 + n.q1 * q.q2 - n.q2 * q.q1 + n.q3 * q.q0;
  return r;
}

//	四元数から座標変換行列Aを求める。
Matrix33 A(Quaternion q)
{
  Matrix33 a;
  a.e11 = q.q0 * q.q0 + q.q1 * q.q1 - q.q2 * q.q2 - q.q3 * q.q3;
  a.e12 = 2.0 * (q.q1 * q.q2 - q.q0 * q.q3);
  a.e13 = 2.0 * (q.q1 * q.q3 + q.q0 * q.q2);
  a.e21 = 2.0 * (q.q1 * q.q2 + q.q0 * q.q3);
  a.e22 = q.q0 * q.q0 - q.q1 * q.q1 + q.q2 * q.q2 - q.q3 * q.q3;
  a.e23 = 2.0 * (q.q2 * q.q3 - q.q0 * q.q1);
  a.e31 = 2.0 * (q.q1 * q.q3 - q.q0 * q.q2);
  a.e32 = 2.0 * (q.q2 * q.q3 + q.q0 * q.q1);
  a.e33 = q.q0 * q.q0 - q.q1 * q.q1 - q.q2 * q.q2 + q.q3 * q.q3;
  return a;
}

//	四元数から座標変換行列ATを求める。
Matrix33 AT(Quaternion q)
{
  Matrix33 a;
  a.e11 = q.q0 * q.q0 + q.q1 * q.q1 - q.q2 * q.q2 - q.q3 * q.q3;
  a.e12 = 2.0 * (q.q1 * q.q2 + q.q0 * q.q3);
  a.e13 = 2.0 * (q.q1 * q.q3 - q.q0 * q.q2);
  a.e21 = 2.0 * (q.q1 * q.q2 - q.q0 * q.q3);
  a.e22 = q.q0 * q.q0 - q.q1 * q.q1 + q.q2 * q.q2 - q.q3 * q.q3;
  a.e23 = 2.0 * (q.q2 * q.q3 + q.q0 * q.q1);
  a.e31 = 2.0 * (q.q1 * q.q3 + q.q0 * q.q2);
  a.e32 = 2.0 * (q.q2 * q.q3 - q.q0 * q.q1);
  a.e33 = q.q0 * q.q0 - q.q1 * q.q1 - q.q2 * q.q2 + q.q3 * q.q3;
  return a;
}

Quaternion operator+(const Quaternion& a,const Vector3& v)
{
  return Quaternion(v) * a;
}

Quaternion Quaternion::operator+=(const Vector3& v)
{
  static double a,aq0,aq1,aq2,aq3,rq0,rq1,rq2,rq3;
  a = sqrt(v.x * v.x + v.y * v.y + v.z * v.z);
  aq0 = cos(a / 2.0);
  if (a != 0.0) {
    aq1 = v.x * sin(a / 2.0) / a;
    aq2 = v.y * sin(a / 2.0) / a;
    aq3 = v.z * sin(a / 2.0) / a;
  } else {
    aq1 = aq2 = aq3 = 0.0;
  }
  rq0 = aq0 * q0 - aq1 * q1 - aq2 * q2 - aq3 * q3;
  rq1 = aq0 * q1 + aq1 * q0 + aq2 * q3 - aq3 * q2;
  rq2 = aq0 * q2 - aq1 * q3 + aq2 * q0 + aq3 * q1;
  rq3 = aq0 * q3 + aq1 * q2 - aq2 * q1 + aq3 * q0;
  q0 = rq0;
  q1 = rq1;
  q2 = rq2;
  q3 = rq3;
  return *this;
}

