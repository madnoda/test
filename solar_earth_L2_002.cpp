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

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdarg.h>
#include "vector3.h"

const double GMs = 1.32712438e20;
const double GMe = 3.9860044e14;
const double AU = 1.4959787e11;
double gcs,gce,a_e,omega_e,t_e,rl2;
const int n_satellite = 7;

struct dXV;
struct XV {
    Vector3 x;
    Vector3 v;
    XV();
    XV(const XV& s);
    XV operator= (const XV& s);
    void operator+= (const dXV& s);
};

struct dXV {
    Vector3 x;
    Vector3 v;
    dXV();
    dXV(const dXV& s);
    dXV operator= (const dXV& s);
    void operator+= (const dXV& s);
    void operator*= (double s);
    void operator/= (double s);
};

XV::XV()
{
    x.x = x.y = x.z = 0.0;
    v.x = v.y = v.z = 0.0;
}

XV::XV(const XV& s)
{
    x = s.x;
    v = s.v;
}

XV XV::operator= (const XV& s)
{
    x = s.x;
    v = s.v;
    return *this;
}

void XV::operator+= (const dXV& s)
{
    x += s.x;
    v += s.v;
}

dXV::dXV()
{
    x.x = x.y = x.z = 0.0;
    v.x = v.y = v.z = 0.0;
}

dXV::dXV(const dXV& s)
{
    x = s.x;
    v = s.v;
}

dXV dXV::operator= (const dXV& s)
{
    x = s.x;
    v = s.v;
    return *this;
}

void dXV::operator+= (const dXV& s)
{
    x += s.x;
    v += s.v;
}

void dXV::operator*= (double s)
{
    x.x *= s;
    x.y *= s;
    x.z *= s;
    v.x *= s;
    v.y *= s;
    v.z *= s;
}

void dXV::operator/= (double s)
{
    x.x /= s;
    x.y /= s;
    x.z /= s;
    v.x /= s;
    v.y /= s;
    v.z /= s;
}

class Satellite
{
    void f(dXV&, const XV&,double);
public:
    XV xv;
    XV runge(double, double);
    void errPuts();
};


void Satellite::f(dXV& dx, const XV& x, double t)
{
    static Vector3 solor;
    static Vector3 earth;
    solor.x = -gcs * cos(t * omega_e);
    solor.y = -gcs * sin(t * omega_e);
    solor.z = 0.0;
    earth.x = gce * cos(t * omega_e);
    earth.y = gce * sin(t * omega_e);
    earth.z = 0.0;
    double r2s = (x.x - solor) * (x.x - solor);
    double r3s = r2s * sqrt(r2s);
    double r2e = (x.x - earth) * (x.x - earth);
    double r3e = r2e * sqrt(r2e);
    dx.v = (-GMs / r3s) * (x.x - solor);
    dx.v += (-GMe / r3e) * (x.x - earth);
    dx.x = x.v;

}

XV Satellite::runge(double t, double h)
{
    static dXV k1, k2, k3, k4, kwork;
    static XV work;
    f(k1, xv, t);
    k1 *= h;
    work = xv;
    kwork = k1;
    kwork /= 2.0;
    work += kwork;
    f(k2, work, t + h / 2.0);
    k2 *= h;
    work = xv;
    kwork = k2;
    kwork /= 2.0;
    work += kwork;
    f(k3, work, t + h / 2.0);
    k3 *= h;
    work = xv;
    work += k3;
    f(k4, work, t + h);
    k4 *= h;
    k2 *= 2.0;
    k3 *= 2.0;
    k1 += k2;
    k1 += k3;
    k1 += k4;
    k1 /= 6.0;
    work = xv;
    work += k1;
    return work;
//    xv += k1;
}

double get_time_in_limit_range(double te,double x,double vy, double z,double limit_range)
{
  Satellite satellite;
  satellite.xv.x.x = x;
  satellite.xv.x.y = 0.0;
  satellite.xv.x.z = z;
  satellite.xv.v.x = 0.0;
  satellite.xv.v.y = vy;
  satellite.xv.v.z = 0.0;
  double step = 60.0;
  double t = 0.0;
  double min_r = -1.0;
  while (t < te) {
	satellite.xv = satellite.runge(t,step);
	t += step;
    // 回転座標系にする
    double rx = satellite.xv.x.x * cos(t * omega_e) + satellite.xv.x.y * sin(t * omega_e) - gce - rl2;
    double ry = - satellite.xv.x.x * sin(t * omega_e) + satellite.xv.x.y * cos(t * omega_e);
    double rz = satellite.xv.x.z;
    if (limit_range > 0.0) {
      if (limit_range * limit_range < rx * rx + ry * ry + rz) {
        return t;
      }
    } else {
      if (min_r < rx * rx + ry * ry + rz) {
        min_r = rx * rx + ry * ry + rz;
      }
    }
  }
  if (limit_range > 0.0) {
    return t;
  }
  if (min_r > 0.0) return sqrt(min_r);
  return 0.0;
}

double get_minmin(double te,double x,double vy0,double vyw,int ny, double z,double limit_range)
{
  double vy_h = vyw / double(ny - 1);
  double vy = - vyw / 2.0;
  double minmin = -1.0;
  double min_vy = -1.0;
  for (int i = 0;i < ny;i++) {
    double r = get_time_in_limit_range(te,x,vy0 + vy,z,limit_range);
    if (limit_range > 0.0) {
      if (minmin < r) {
        minmin = r;
        min_vy = vy;
      }
      printf("%d/%d vy:%f r:%f day\n",i,ny,vy0 + vy,r / (24.0 * 3500.0));
    } else {
      if ((minmin < 0.0) || (minmin > r)) {
        minmin = r;
        min_vy = vy;
      }
      printf("%d/%d vy:%f r:%f\n",i,ny,vy0 + vy,r / 1.0e7);
    }
    vy += vy_h;
  }
  return vy0 + min_vy;
}

double get_from_init_position(double te,double x,double vy, double z)
{
  Satellite satellite;
  satellite.xv.x.x = x;
  satellite.xv.x.y = 0.0;
  satellite.xv.x.z = z;
  satellite.xv.v.x = 0.0;
  satellite.xv.v.y = vy;
  satellite.xv.v.z = 0.0;
  double step = 60.0;
  double t = 0.0;
  double min_r = -1.0;
  while (t < te * 0.9) {
	satellite.xv = satellite.runge(t,step);
	t += step;
  }
  // 回転座標系にする
  double ry = - satellite.xv.x.x * sin(t * omega_e) + satellite.xv.x.y * cos(t * omega_e);
  XV xv_next;
  for (int i = 0;i < 30;i++) {
    while (ry > 0.0) {
      xv_next = satellite.runge(t,step);
      ry = - xv_next.x.x * sin(t * omega_e) + xv_next.x.y * cos(t * omega_e);
      if (ry > 0.0) {
        t += step;
        satellite.xv = xv_next;
      }
    }
//    printf("%d %fday %f \n",i,t/(3600.0 * 24.0),ry);
    ry = 100.0;
    step /= 2.0;
    t -= step;
  }
  double dx = satellite.xv.x.x * cos(t * omega_e) + satellite.xv.x.y * sin(t * omega_e) - x;
  double dy = - satellite.xv.x.x * sin(t * omega_e) + satellite.xv.x.y * cos(t * omega_e);
  double dz = satellite.xv.x.z - z;
//  printf("%fday %f %f %f \n",t/(3600.0 * 24.0),dx/1.0e7,dy,dz/1.0e7);
  return dx * dx + dz * dz;
}

void get_minminmin(double te,double x,double & vy0,double vyw,int ny, double & z0, double zw, int nz)
{
  double vy_h = 0.0;
  if (ny > 1) {
    vy_h = vyw / double(ny - 1);
  } else {
    vyw = 0.0;
  }
  double vy = - vyw / 2.0;
  double z_h = 0.0;
  if (nz > 1) {
    z_h = zw / double(nz - 1);
  } else {
    zw = 0.0;
  }
  double minmin = -1.0;
  double min_vy = -1.0;
  double min_z = -1.0;
  for (int i = 0;i < ny;i++) {
    double z = - zw / 2.0;
    for (int j = 0;j < nz;j++) {
      double r = get_from_init_position(te,x,vy0 + vy,z0 + z);
      if ((minmin < 0.0) || (minmin > r)) {
        minmin = r;
        min_vy = vy;
        min_z = z;
      }
      printf("%d/%d %d/%d vy:%f z:%f r:%f\n",i,ny,j,nz,vy0 + vy,z0 + z,sqrt(r) / 1.0e7);
      z += z_h;
    }
    vy += vy_h;
  }
  vy0 += min_vy;
  z0 += min_z;
}

double get_a(double r) {
  return (gce + r) * omega_e * omega_e - GMs/(AU + r)/(AU + r) - GMe/r/r;
}

double get_da(double r) {
  return omega_e * omega_e + 2.0 * GMs/(AU + r)/(AU + r)/(AU + r) + 2.0 * GMe/r/r/r;
}


main (int argc, char* argv[]) {
  char distance[128];
  char filename[128];
  distance[0] = 0;
  double distance_f;
  filename[0] = 0;;
  int i = 1;
  while (i < argc) {
    int j = 0;
    int k = 0;
    if (i == 1) {
      j = 0;
      int k = 0;
      while (distance[j++] = argv[i][k]) {
        k++;
      }
      if (k) {
        sscanf(distance,"%lf",&distance_f);
      }
    } else if (i == 2) {
      j = 0;
      int k = 0;
      while (filename[j++] = argv[i][k]) {
        k++;
       }
    }
    i++;
  }
  printf("%s %s %f\n",distance,filename,distance_f);
  gcs = AU * GMe / (GMs + GMe);
  gce = AU * GMs / (GMs + GMe);
  a_e = GMs / gce / gce;
  omega_e = sqrt(a_e/gce);
  t_e = 2.0 * M_PI / omega_e;
  printf("%f %f\n",t_e,t_e / 24.0 / 3600.0);
  rl2 = 120e7;
  for (int i = 0;i < 50;i++) {
    rl2 -= get_a(rl2) / get_da(rl2);
  }

  double vyw = distance_f * omega_e * 200.0;
  double vy = (gce + rl2 + distance_f) * omega_e;
  double z = distance_f ;
  double limit_range = distance_f * 7.0;
  vy = get_minmin(t_e,gce + rl2 + distance_f,vy,vyw,400, z,limit_range);
  printf(" vy:%f z:%f\n",vy,z);
  vyw /= 10.0;
  for (int i = 0;i < 5;i++) {
    vy = get_minmin(t_e,gce + rl2 + distance_f,vy,vyw,20, z,limit_range);
    printf(" vy:%f z:%f\n",vy,z);
    if (z < 0.0) {
      z = -z;
    }
    vyw /= 10.0;
  }
  double zw = z;
  get_minminmin(t_e / 2.0,gce + rl2 + distance_f,vy,vyw,10, z,zw,10);
  printf(" vy:%f z:%f\n",vy,z);
  vyw /= 10.0;
  zw /= 10.0;
  for (int i = 0;i < 5;i++) {
    get_minminmin(t_e / 2.0,gce + rl2 + distance_f,vy,vyw,10, z,zw,10);
    printf(" vy:%f z:%f\n",vy,z);
    if (z < 0.0) {
      z = -z;
    }
    vyw /= 10.0;
    zw /= 10.0;
  }

  Satellite satellites[n_satellite];
  double x000 = gce + rl2 + distance_f;
  double z000 = z;
  double vy000 = vy;
  
  for (int j = 0;j < n_satellite;j++) {
    double x = x000 + 100.0 * double(j);
    if (j > n_satellite / 2) {
      x = x000 - 100.0 * double(j - int(n_satellite / 2));
    }
    z = z000;
    zw = z000 / 100.0;
    vy = vy000;
    vyw = 10.0;
    get_minminmin(t_e / 2.0,x,vy,vyw,10, z,zw,10);
    printf(" vy:%f z:%f\n",vy,z);
    vyw /= 10.0;
    zw /= 10.0;
    for (int i = 0;i < 5;i++) {
      get_minminmin(t_e / 2.0,x,vy,vyw,10, z,zw,10);
      printf(" vy:%f z:%f\n",vy,z);
      if (z < 0.0) {
        z = -z;
      }
      vyw /= 10.0;
      zw /= 10.0;
    }
    satellites[j].xv.x.x = x;
    satellites[j].xv.x.y = 0.0;
    satellites[j].xv.x.z = z;
    satellites[j].xv.v.x = 0.0;
    satellites[j].xv.v.y = vy;
    satellites[j].xv.v.z = 0.0;
  }
  for (int j = 0;j < n_satellite;j++) {
    printf("%d %e %e %e\n",j,satellites[j].xv.x.x,satellites[j].xv.x.z,satellites[j].xv.v.y);
    printf("%d %f %f %f\n",j,satellites[j].xv.x.x - gce - rl2,satellites[j].xv.x.z,satellites[j].xv.v.y);
  }
  FILE* fp;
  if ((fp = fopen(filename, "w")) == NULL) {
     fprintf(stderr, "Can not open file :%s\n",filename);
  }
  
  double step = 60.0;
  double t = 0.0;
  i = 0;
  while (t < t_e / 2.0) {
	if ((i % 1000 == 0) && (fp != NULL)) {
      // 回転座標系にする
      double rx0 = satellites[0].xv.x.x * cos(t * omega_e) + satellites[0].xv.x.y * sin(t * omega_e) - gce - rl2;
      double ry0 = - satellites[0].xv.x.x * sin(t * omega_e) + satellites[0].xv.x.y * cos(t * omega_e);
      double rz0 = satellites[0].xv.x.z;
      int fc = fprintf(fp,"%f %f %f %f",t,rx0 / 1.0e7,ry0 / 1.0e7,rz0 / 1.0e7);
      if (fc < 0) {
        fprintf(stderr, "Can not write file %s.\n",filename);
        fclose(fp);
        fp = NULL;
      }
      for (int j = 1;j < n_satellite;j++) {
        double rx = satellites[j].xv.x.x * cos(t * omega_e) + satellites[j].xv.x.y * sin(t * omega_e) - gce - rl2;
        double ry = - satellites[j].xv.x.x * sin(t * omega_e) + satellites[j].xv.x.y * cos(t * omega_e);
        double rz = satellites[j].xv.x.z;
        fc = fprintf(fp," %f %f %f",rx / 1.0e7,ry / 1.0e7,rz / 1.0e7);
        if (fc < 0) {
          fprintf(stderr, "Can not write file %s.\n",filename);
          fclose(fp);
          fp = NULL;
        }
        fc = fprintf(fp," %f %f %f",rx - rx0,ry - ry0,rz - rz0);
        if (fc < 0) {
          fprintf(stderr, "Can not write file %s.\n",filename);
          fclose(fp);
          fp = NULL;
        }
      }
      fc = fprintf(fp,"\n");
      if (fc < 0) {
        fprintf(stderr, "Can not write file %s.\n",filename);
        fclose(fp);
        fp = NULL;
      }

    }
    for (int j = 0;j < n_satellite;j++) {
      satellites[j].xv = satellites[j].runge(t,step);
    }
	t += step;
    i++;
  }
  if (fp != NULL) {
    fclose(fp);
  }
}

