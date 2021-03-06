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

#define STDCALLbyMAD __stdcall


struct Vector3
{
  double x,y,z;
  Vector3(){ x = y = z = 0.0;};
  Vector3(const Vector3& s) {
    x = s.x;
    y = s.y;
    z = s.z;
  };
  Vector3(double xin,double yin,double zin){
    x = xin;
    y = yin;
    z = zin;
  };
  double abs() const;
	       Vector3 unit();
  Vector3 operator=(const Vector3&s) {
    x = s.x;
    y = s.y;
    z = s.z;
    return *this;
  };
  Vector3 operator-() const;
		      Vector3 operator+=(const Vector3&);
  Vector3 operator-=(const Vector3&);
  Vector3 operator*=(double);
  Vector3 operator/=(double);
};

Vector3 operator^(const Vector3&,const Vector3&);
double operator*(const Vector3& ,const Vector3&);
Vector3 operator+(const Vector3& ,const Vector3&);
Vector3 operator-(const Vector3& ,const Vector3&);
Vector3 operator*(double ,const Vector3&);
Vector3 operator/(const Vector3&,double);

struct Matrix33
{
  double e11,e12,e13;
  double e21,e22,e23;
  double e31,e32,e33;
  Matrix33() {
    e11 = e12 = e13 = 0.0;
    e21 = e22 = e23 = 0.0;
    e31 = e32 = e33 = 0.0;
  };
  Matrix33(const Matrix33&s) {
    e11 = s.e11;
    e12 = s.e12;
    e13 = s.e13;
    e21 = s.e21;
    e22 = s.e22;
    e23 = s.e23;
    e31 = s.e31;
    e32 = s.e32;
    e33 = s.e33;
  };
  Matrix33(double i11,double i12,double i13,double i21,double i22,double i23,double i31,double i32,double i33) {
    e11 = i11;
    e12 = i12;
    e13 = i13;
    e21 = i21;
    e22 = i22;
    e23 = i23;
    e31 = i31;
    e32 = i32;
    e33 = i33;
  };
  Matrix33 operator =(const Matrix33&s) {
    e11 = s.e11;
    e12 = s.e12;
    e13 = s.e13;
    e21 = s.e21;
    e22 = s.e22;
    e23 = s.e23;
    e31 = s.e31;
    e32 = s.e32;
    e33 = s.e33;
    return *this;
  };
  double & operator()(int i,int j);
  Matrix33 operator*=(double);
  Vector3 GetX() const;
  Vector3 GetY() const;
  Vector3 GetZ() const;
  void SetXYZ(const Vector3&,const Vector3&,const Vector3&);
  Matrix33 t() const;
};

Matrix33 operator*(const Matrix33 &,const Matrix33&);
Vector3 operator*(const Matrix33 &,const Vector3&);

struct Quaternion
{
  double q0,q1,q2,q3;
  Quaternion(){
    q0 = 1.0;
    q1 = q2 = q3 = 0.0;
  };
  Quaternion(double i0,double i1,double i2,double i3){
    q0 = i0;
    q1 = i1;
    q2 = i2;
    q3 = i3;
  };
  Quaternion(const Quaternion&s){
    q0 = s.q0;
    q1 = s.q1;
    q2 = s.q2;
    q3 = s.q3;
  };
  Quaternion(const Vector3&);
  Quaternion operator =(const Quaternion&s);
  Quaternion operator =(const Matrix33 & s);
  void kikakuka();
  Quaternion operator -();
  Quaternion operator +=(const Vector3&);
  void GetV(double &,double &,double &,double &);
};

Quaternion operator*(const Quaternion& ,const Quaternion&);

Matrix33 A(Quaternion q);	//	四元数から座標変換行列Aを求める。
Matrix33 AT(Quaternion q);	//	四元数から座標変換逆行列ATを求める。

Quaternion operator+(const Quaternion& ,const Vector3&);

