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

#include <quadmath.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdarg.h>
#include "vector3q.h"

const __float128 step_base = 600.0q;

const __float128 GMs = 1.32712438e20q;
const __float128 GMe = 3.9860044e14q;
const __float128 AU = 1.4959787e11q;
__float128 gcs,gce,a_e,omega_e,t_e,rl2;
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
    void operator*= (__float128 s);
    void operator/= (__float128 s);
};

XV::XV()
{
    x.x = x.y = x.z = 0.0q;
    v.x = v.y = v.z = 0.0q;
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
    x.x = x.y = x.z = 0.0q;
    v.x = v.y = v.z = 0.0q;
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

void dXV::operator*= (__float128 s)
{
    x.x *= s;
    x.y *= s;
    x.z *= s;
    v.x *= s;
    v.y *= s;
    v.z *= s;
}

void dXV::operator/= (__float128 s)
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
    void f(dXV&, const XV&,__float128);
public:
    XV xv;
    XV runge(__float128, __float128);
    void errPuts();
};


void Satellite::f(dXV& dx, const XV& x, __float128 t)
{
    static Vector3 solor;
    static Vector3 earth;
    solor.x = -gcs * cosq(t * omega_e);
    solor.y = -gcs * sinq(t * omega_e);
    solor.z = 0.0q;
    earth.x = gce * cosq(t * omega_e);
    earth.y = gce * sinq(t * omega_e);
    earth.z = 0.0q;
    __float128 r2s = (x.x - solor) * (x.x - solor);
    __float128 r3s = r2s * sqrtq(r2s);
    __float128 r2e = (x.x - earth) * (x.x - earth);
    __float128 r3e = r2e * sqrtq(r2e);
    dx.v = (-GMs / r3s) * (x.x - solor);
    dx.v += (-GMe / r3e) * (x.x - earth);
    dx.x = x.v;

}

XV Satellite::runge(__float128 t, __float128 h)
{
    static dXV k1, k2, k3, k4, kwork;
    static XV work;
    f(k1, xv, t);
    k1 *= h;
    work = xv;
    kwork = k1;
    kwork /= 2.0q;
    work += kwork;
    f(k2, work, t + h / 2.0q);
    k2 *= h;
    work = xv;
    kwork = k2;
    kwork /= 2.0q;
    work += kwork;
    f(k3, work, t + h / 2.0q);
    k3 *= h;
    work = xv;
    work += k3;
    f(k4, work, t + h);
    k4 *= h;
    k2 *= 2.0q;
    k3 *= 2.0q;
    k1 += k2;
    k1 += k3;
    k1 += k4;
    k1 /= 6.0q;
    work = xv;
    work += k1;
    return work;
//    xv += k1;
}

__float128 get_period(__float128 te,__float128 x,__float128 vy, __float128 z)
{
  Satellite satellite;
  satellite.xv.x.x = x;
  satellite.xv.x.y = 0.0q;
  satellite.xv.x.z = z;
  satellite.xv.v.x = 0.0q;
  satellite.xv.v.y = vy;
  satellite.xv.v.z = 0.0q;
  __float128 step = step_base;
  __float128 t = 0.0q;

  while (t < te * 0.7q) {
    satellite.xv = satellite.runge(t,step);
    t += step;
  }

  __float128 ry = - 100.0q;
  XV xv_next;
  while (ry < 0.0q) {
    xv_next = satellite.runge(t,step);
    // 回転座標系にする
    ry = - xv_next.x.x * sinq(t * omega_e) + xv_next.x.y * cosq(t * omega_e);
    if (ry < 0.0q) {
      t += step;
      if (t > te * 1.5q) {
        return -1.0q;
      }
      satellite.xv = xv_next;
    }
  }
  ry = -100.0q;
  step /= 2.0q;
  bool flg = true;
  for (int i = 0;flg &&(i < 100);i++) {
    int j = 0;
    while (flg && (j < 1000) && (ry < 0.0q)) {
      xv_next = satellite.runge(t,step);
      // 回転座標系にする
      ry = - xv_next.x.x * sinq(t * omega_e) + xv_next.x.y * cosq(t * omega_e);
      if (ry < 0.0q) {
        t += step;
        if (t > te * 1.5q) {
          return -1.0q;
        }
        satellite.xv = xv_next;
      }
      j++;
    }
    ry = -100.0q;
    step /= 2.0q;
  }
//  printf("period:%fday\n",(double)(t / (24.0q * 3600.0q)));
  return t;
}

__float128 get_time_in_limit_range(__float128 te,__float128 x,__float128 vy, __float128 z,__float128 limit_range)
{
  Satellite satellite;
  satellite.xv.x.x = x;
  satellite.xv.x.y = 0.0q;
  satellite.xv.x.z = z;
  satellite.xv.v.x = 0.0q;
  satellite.xv.v.y = vy;
  satellite.xv.v.z = 0.0q;
  __float128 step = step_base;
  __float128 t = 0.0q;
  __float128 min_r = -1.0q;
  while (t < te) {
    satellite.xv = satellite.runge(t,step);
    t += step;
    // 回転座標系にする
    __float128 rx = satellite.xv.x.x * cosq(t * omega_e) + satellite.xv.x.y * sinq(t * omega_e) - gce - rl2;
    __float128 ry = - satellite.xv.x.x * sinq(t * omega_e) + satellite.xv.x.y * cosq(t * omega_e);
    __float128 rz = satellite.xv.x.z;
    if (limit_range > 0.0q) {
      if (limit_range * limit_range < rx * rx + ry * ry + rz) {
        return t;
      }
    } else {
      if (min_r < rx * rx + ry * ry + rz) {
        min_r = rx * rx + ry * ry + rz;
      }
    }
  }
  if (limit_range > 0.0q) {
    return t;
  }
  if (min_r > 0.0q) return sqrtq(min_r);
  return 0.0q;
}

__float128 get_min(__float128 te,__float128 x,__float128 vy0,__float128 vyw,int ny, __float128 z,__float128 limit_range)
{
  __float128 vy_h = vyw / __float128(ny - 1);
  __float128 vy = - vyw / 2.0q;
  __float128 minmin = -1.0q;
  __float128 min_vy = -1.0q;
  for (int i = 0;i < ny;i++) {
    __float128 r = get_time_in_limit_range(te,x,vy0 + vy,z,limit_range);
    if (limit_range > 0.0q) {
      if (minmin < r) {
        minmin = r;
        min_vy = vy;
      }
      printf("%d/%d vy:%f z:%f r:%f day\n",i,ny,(double)(vy0 + vy),(double)z,(double)(r / (24.0q * 3500.0q)));
    } else {
      if ((minmin < 0.0q) || (minmin > r)) {
        minmin = r;
        min_vy = vy;
      }
      printf("%d/%d vy:%f z:%f r:%f\n",i,ny,(double)(vy0 + vy),(double)z,(double)(r / 1.0e7));
    }
    vy += vy_h;
  }
  return vy0 + min_vy;
}

bool get_from_init_position(__float128 te,__float128 x,__float128 vy, __float128 z,__float128 &dx,__float128 &dz)
{
  Satellite satellite;
  satellite.xv.x.x = x;
  satellite.xv.x.y = 0.0q;
  satellite.xv.x.z = z;
  satellite.xv.v.x = 0.0q;
  satellite.xv.v.y = vy;
  satellite.xv.v.z = 0.0q;
  __float128 step = step_base;
  __float128 t = 0.0q;
  __float128 min_r = -1.0q;

  while (t < te * 0.7q) {
    satellite.xv = satellite.runge(t,step);
    t += step;
  }

  __float128 ry = - 100.0q;
  XV xv_next;
  while (ry < 0.0q) {
    xv_next = satellite.runge(t,step);
    // 回転座標系にする
    ry = - xv_next.x.x * sinq(t * omega_e) + xv_next.x.y * cosq(t * omega_e);
    if (ry < 0.0q) {
      t += step;
      if (t > te * 1.5q) {
        dx = -1.0q;
        return false;
      }
      satellite.xv = xv_next;
    }
  }
  ry = -100.0q;
  step /= 2.0q;
  bool flg = true;
  for (int i = 0;flg &&(i < 100);i++) {
    int j = 0;
    while (flg && (j < 1000) && (ry < 0.0q)) {
      xv_next = satellite.runge(t,step);
      // 回転座標系にする
      ry = - xv_next.x.x * sinq(t * omega_e) + xv_next.x.y * cosq(t * omega_e);
      if (ry < 0.0q) {
        t += step;
        if (t > te * 1.5q) {
          dx = -1.0q;
          return false;
        }
        satellite.xv = xv_next;
//        if (ry > -1.0e-6) {
//          flg = false;
//        }
      }
      j++;
    }
    ry = -100.0q;
    step /= 2.0q;
  }
  dx = fabsq(satellite.xv.x.x * cosq(t * omega_e) + satellite.xv.x.y * sinq(t * omega_e) - x);
  dz = satellite.xv.x.z - z;
  return true;
}

void get_minmin(bool flgOptX,__float128 te,__float128 x,__float128 & vy0,__float128 vyw,int ny, __float128 & z0, __float128 zw, int nz,__float128 &dx,__float128 &dz)
{
  int i = 0;
  int l = 0;
  __float128 minmin_old = 0.0q;
  __float128 min_dx_old = 0.0q;
  __float128 min_dz_old = 0.0q;

  while (i < 30){
    __float128 vy_h = 0.0q;
    if (ny > 1) {
      vy_h = vyw / __float128(ny - 1);
    } else {
      vyw = 0.0q;
    }
    __float128 vy = - vyw / 2.0q;
    __float128 z_h = 0.0q;
    if (nz > 1) {
      z_h = zw / __float128(nz - 1);
    } else {
      zw = 0.0q;
    }
    __float128 minmin = -1.0q;
    __float128 min_vy = -1.0q;
    __float128 min_z = -1.0q;
    int min_j = -1;
    int min_k = -1;
    __float128 min_dx = -1.0q;
    __float128 min_dz = -1.0q;
    for (int j = 0;j < ny;j++) {
      __float128 z = - zw / 2.0q;
      for (int k = 0;k < nz;k++) {
        __float128 dx,dz;
        bool flg = get_from_init_position(te,x,vy0 + vy,z0 + z,dx,dz);
        __float128 r = fabsq(dx);
        if ((flg) && ((minmin < 0.0q) || (minmin > r))) {
          minmin = r;
          min_vy = vy;
          min_z = z;
          min_j = j;
          min_k = k;
          min_dx = dx;
          min_dz = dz;
        }
        if (flg)
          printf("%d %d/%d %d/%d vy:%f z:%f r:%f dx:%f dz:%f\n",i,j,ny,k,nz,(double)(vy0 + vy),(double)(z0 + z),(double)(r / 1.0e7),(double)min_dx,(double)(min_dz / 1.0e7));
        z += z_h;
      }
      vy += vy_h;
    }
    vy0 += min_vy;
    z0 += min_z;
    if (z0 > 0.0q) {
      z0 = -z0;
    }
    if (fabsq(minmin - minmin_old) / minmin < 0.01q) {
      l++;
      if (l > 3) i = 100;
    } else {
      l = 0;
    }
    printf("l:%d / %f\n",l,(double)(fabsq(minmin - minmin_old) / minmin));
    minmin_old = minmin;
    min_dx_old = min_dx;
    min_dz_old = min_dz;
    if ((ny > 1) && (nz > 1)) {
      if ((min_j != 0) && (min_j != ny - 1) && (min_k != 0) && (min_k != nz - 1)) {
        vyw /= 3.0q;
        zw /= 3.0q;
        i++;
      } else {
        if ((min_j == 0) || (min_j == ny - 1)) {
          vyw *= 3.0q;
        }
        if ((min_k == 0) || (min_k == nz - 1)) {
          zw *= 3.0q;
        }
      }
    } else if (ny < 2) {
      if ((min_k != 0) && (min_k != nz - 1)) {
        zw /= 3.0q;
        i++;
      } else {
        zw *= 3.0q;
      }
    } else {
      if ((min_j != 0) && (min_j != ny - 1)) {
        vyw /= 3.0q;
        i++;
      } else {
        vyw *= 3.0q;
      }
    }
    printf("min_j:%d min_k:%d vyw:%f zw:%f min_r:%f dx:%f dz:%f\n",min_j,min_k,(double)vyw,(double)zw,(double)(minmin / 1.0e7),(double)min_dx,(double)(min_dz / 1.0e7));
    printf("vy:%f z:%f min:%f\n",(double)vy0,(double)z0,(double)sqrtq(minmin));
  }
  dx = min_dx_old;
  dz = min_dz_old;
  if (flgOptX) {
    __float128 dx0,dz0;
    __float128 d_vy = 0.01q / 100.0q;
    __float128 old_r = -1.0q;
    bool flg = get_from_init_position(te,x,vy0,z0,dx0,dz0);
    for (int i = 0;i < 10;i++) {
      old_r = dx0;
      __float128 dx_vy,dz_vy;
      flg = get_from_init_position(te,x,vy0 + d_vy,z0,dx_vy,dz_vy);
      __float128 a = (dx_vy - dx0) / d_vy;
      const __float128 TINY_A = 0.0000000001q; // 1.0e-10
      if (fabsq(a) > TINY_A) {
        __float128 yyy = (-dx0 / a) * 2.0q;
        int j = 0;
        do {
          yyy /= 2.0q;
          flg = get_from_init_position(te,x,vy0 + yyy,z0,dx_vy,dz_vy);
          printf("%f %f\n",(double)old_r,(double)dx_vy);
          j++;
        } while ((j < 20) && (fabsq(old_r) < fabsq(dx_vy)));
        vy0 += yyy;
      }
      flg = get_from_init_position(te,x,vy0,z0,dx0,dz0);
    }
    dx = dx0;
    dz = dz0;
  }
}

void get_minminmin(__float128 te,__float128 x,__float128 & vy0, __float128 & z0, __float128 &dx0, __float128 &dz0)
{
//  __float128 dx0,dz0;
  __float128 d_vy = 0.01q / 100.0q;
  __float128 d_z = 10000000.0q / 100.0q;
  bool flg;
  __float128 old_r = -1.0q;
  int k = 0;
  flg = get_from_init_position(te,x,vy0,z0,dx0,dz0);
  for (int i = 0;i < 10000;i++) {
    if (sqrtq(dx0 * dx0 + dz0 * dz0) < 10000.0) {
      printf("i:%d vy:%f z:%f dx:%f dz:%f r:%f(m) \n",i,(double)vy0,(double)z0,(double)dx0,(double)dz0,(double)(sqrtq(dx0 * dx0 + dz0 * dz0)));
    } else {
      printf("i:%d vy:%f z:%f dx:%f dz:%f r:%f \n",i,(double)vy0,(double)z0,(double)dx0,(double)dz0,(double)(sqrtq(dx0 * dx0 + dz0 * dz0) / 1.0e7));
    }
    old_r = dx0 * dx0 + dz0 * dz0;
    __float128 dx_vy,dz_vy;
    flg = get_from_init_position(te,x,vy0 + d_vy,z0,dx_vy,dz_vy);
    __float128 dx_z,dz_z;
    flg = get_from_init_position(te,x,vy0,z0 + d_z, dx_z,dz_z);
  
    __float128 a = (dx_vy - dx0) / d_vy;
    __float128 c = (dz_vy - dz0) / d_vy;
    __float128 b = (dx_z - dx0) / d_z;
    __float128 d = (dz_z - dz0) / d_z;
    __float128 yyy = -(d * dx0 - b * dz0) / (a * d - b * c);
    __float128 zzz = -(-c * dx0 + a * dz0) / (a * d - b * c);
    yyy /= 2.0q;
    zzz /= 2.0q;
    printf("yyy: %f zzz:%f:\n",(double)yyy, (double)zzz);
    const __float128 TINY = 0.000000000001q; // 1.0e-12
    if (d_vy > fabsq(yyy) / 100.0q) {
      d_vy = fabsq(yyy) / 100.0q;
      if (d_vy < TINY) {
        d_vy = TINY;
      }
    }
    if (d_z > fabsq(zzz) / 100.0q) {
      d_z = fabsq(zzz) / 100.0q;
      if (d_z < TINY) {
        d_z = TINY;
      }
    }
    flg = get_from_init_position(te,x,vy0 + yyy,z0 + zzz,dx0,dz0);
    __float128 new_r  = dx0 * dx0 + dz0 * dz0;
    int j = 0;
    while ((j < 1000) && (new_r > old_r)) {
      yyy /= 2.0q;
      zzz /= 2.0q;
      flg = get_from_init_position(te,x,vy0 + yyy,z0 + zzz,dx0,dz0);
      new_r  = dx0 * dx0 + dz0 * dz0;
      j++;
    }
    if (d_vy > fabsq(yyy) / 100.0q) {
      d_vy = fabsq(yyy) / 100.0q;
      if (d_vy < TINY) {
        d_vy = TINY;
      }
    }
    if (d_z > fabsq(zzz) / 100.0q) {
      d_z = fabsq(zzz) / 100.0q;
      if (d_z < TINY) {
        d_z = TINY;
      }
    }
    old_r = new_r;
    printf("d_vy: %e d_z:%e j:%d\n",(double)d_vy, (double)d_z,j);
    printf("yyy: %e zzz:%e:\n",(double)yyy, (double)zzz);
    vy0 += yyy;
    z0 += zzz;
    const __float128 treshold = 0.001q;
    const __float128 treshold2 = 0.01q;
    if (sqrtq(dx0 * dx0 + dz0 * dz0) < treshold) {
      printf("i:%d vy:%f z:%f dx:%f dz:%f r:%f(m) \n",i,(double)vy0,(double)z0,(double)dx0,(double)dz0,(double)(sqrtq(dx0 * dx0 + dz0 * dz0)));
      i = 10001;
    } else {
      if (sqrtq(dx0 * dx0 + dz0 * dz0) < treshold2) {
        k++;
        if (k > 20) {
          if (sqrtq(dx0 * dx0 + dz0 * dz0) < 1000.0q) {
            printf("i:%d vy:%f z:%f dx:%f dz:%f r:%f (m)\n",i,(double)vy0,(double)z0,(double)dx0,(double)dz0,(double)(sqrtq(dx0 * dx0 + dz0 * dz0)));
          } else {
            printf("i:%d vy:%f z:%f dx:%f dz:%f r:%f (km)\n",i,(double)vy0,(double)z0,(double)dx0,(double)dz0,(double)(sqrtq(dx0 * dx0 + dz0 * dz0) / 1000.0));
          }
          i = 10001;
        }
      } else {
        k = 0;
      }
    }
  }
}

bool get_minminminmin(__float128 d_q,__float128 &vy,__float128 &z,__float128 &dx,__float128 &dz)
{
  if (fabsq(z) < 0.1q) {
    __float128 vyw = d_q * omega_e * 200.0q;
    vy = (gce + rl2 + d_q) * omega_e;
    z = - d_q;
    __float128 limit_range = d_q * 4.0q;
    vy = get_min(t_e,gce + rl2 - d_q,vy,vyw,400, z,limit_range);
    printf(" vy:%f z:%f\n",(double)vy,(double)z);
    vyw /= 10.0q;
    for (int i = 0;i < 5;i++) {
      vy = get_min(t_e,gce + rl2 - d_q,vy,vyw,20, z,limit_range);
      printf(" vy:%f z:%f\n",(double)vy,(double)z);
      if (z > 0.0q) {
        z = -z;
      }
      vyw /= 10.0q;
    }
    __float128 zw = -z * 3.0q;
    vyw *= 1000.0q;
    get_minmin(false,t_e / 2.0q,gce + rl2 - d_q,vy,vyw,21, z,zw,1,dx,dz);
    printf(" vy:%f z:%f r%f:\n",(double)vy,(double)z,(double)(sqrtq(dx * dx + dz * dz) / 1.0e7));
  }
  get_minminmin(t_e / 2.0q,gce + rl2 - d_q,vy,z,dx,dz);
  return true;
}

__float128 get_a(__float128 r) {
  return (gce + r) * omega_e * omega_e - GMs/(AU + r)/(AU + r) - GMe/r/r;
}

__float128 get_da(__float128 r) {
  return omega_e * omega_e + 2.0q * GMs/(AU + r)/(AU + r)/(AU + r) + 2.0q * GMe/r/r/r;
}

main (int argc, char* argv[]) {
  char distance[128];
  char filename[128];
  distance[0] = 0;
  double distance_d;
  __float128 distance_q;
  filename[0] = 0;;
  int i = 1;
  printf("argc:%d\n",argc);
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
        sscanf(distance,"%lf",&distance_d);
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
  distance_q = distance_d;
  printf("%s %s %f\n",distance,filename,(double)distance_q);
  gcs = AU * GMe / (GMs + GMe);
  gce = AU * GMs / (GMs + GMe);
  a_e = GMs / gce / gce;
  omega_e = sqrtq(a_e/gce);
  t_e = 2.0q * M_PI / omega_e;
  printf("%f %f\n",(double)t_e,(double)(t_e / 24.0q / 3600.0q));
  rl2 = 120e7;
  for (int i = 0;i < 50;i++) {
    rl2 -= get_a(rl2) / get_da(rl2);
  }
  __float128 vy,z,dx,dz;
  z = 0.0q;
  get_minminminmin(distance_q,vy,z,dx,dz);
//  printf(" vy:%f z:%f\n",(double)vy,(double)z);
  char buf[128];
  printf(" vy:%f z:%f r:%f\n",(double)vy,(double)z,(double)(sqrtq(dx * dx + dz * dz) / 1.0e7));
  __float128 t_period = get_period(t_e / 2.0,gce + rl2 - distance_q,vy,z);
  int n = quadmath_snprintf(buf, sizeof buf, "%+-#*.20Qe", 46, vy);
  if ((size_t)n < sizeof buf) printf ("vy:%s\n", buf);
  n = quadmath_snprintf(buf, sizeof buf, "%+-#*.20Qe", 46, z);
  if ((size_t)n < sizeof buf) printf ("z:%s\n", buf);

  Satellite satellites[n_satellite];

  satellites[0].xv.x.x = gce + rl2 - distance_q;
  satellites[0].xv.x.y = 0.0;
  satellites[0].xv.x.z = z;
  satellites[0].xv.v.x = 0.0;
  satellites[0].xv.v.y = vy;
  satellites[0].xv.v.z = 0.0;

  for (int j = 1;j < n_satellite;j++) {
    __float128 d_q = distance_q + 100.0 * __float128(j);
    if (j > n_satellite / 2) {
      d_q = distance_q - 100.0 * __float128(j - int(n_satellite / 2));
    }
    get_minminminmin(d_q,vy,z,dx,dz);
    satellites[j].xv.x.x = gce + rl2 - d_q;
    satellites[j].xv.x.y = 0.0;
    satellites[j].xv.x.z = z;
    satellites[j].xv.v.x = 0.0;
    satellites[j].xv.v.y = vy;
    satellites[j].xv.v.z = 0.0;
  }
  for (int j = 0;j < n_satellite;j++) {
    printf("%d %e %e %e\n",j,(double)satellites[j].xv.x.x,(double)satellites[j].xv.x.z,(double)satellites[j].xv.v.y);
    printf("%d %f %f %f\n",j,(double)satellites[j].xv.x.x - gce - rl2,(double)satellites[j].xv.x.z,(double)satellites[j].xv.v.y);
  }
  FILE* fp;
  if ((fp = fopen(filename, "w")) == NULL) {
     fprintf(stderr, "Can not open file :%s\n",filename);
  }
  
  __float128 step = step_base;
  __float128 t = 0.0;
  i = 0;
  while (t < t_e / 2.0) {
	if ((i % 1000 == 0) && (fp != NULL)) {
      // 回転座標系にする
      __float128 rx0 = satellites[0].xv.x.x * cosq(t * omega_e) + satellites[0].xv.x.y * sinq(t * omega_e) - gce - rl2;
      __float128 ry0 = - satellites[0].xv.x.x * sinq(t * omega_e) + satellites[0].xv.x.y * cosq(t * omega_e);
      __float128 rz0 = satellites[0].xv.x.z;
      int fc = fprintf(fp,"%f %f %f %f",(double)t,(double)(rx0 / 1.0e7),(double)(ry0 / 1.0e7),(double)(rz0 / 1.0e7));
      if (fc < 0) {
        fprintf(stderr, "Can not write file %s.\n",filename);
        fclose(fp);
        fp = NULL;
      }
      for (int j = 1;j < n_satellite;j++) {
        __float128 rx = satellites[j].xv.x.x * cosq(t * omega_e) + satellites[j].xv.x.y * sinq(t * omega_e) - gce - rl2;
        __float128 ry = - satellites[j].xv.x.x * sinq(t * omega_e) + satellites[j].xv.x.y * cosq(t * omega_e);
        __float128 rz = satellites[j].xv.x.z;
        fc = fprintf(fp," %f %f %f",(double)(rx / 1.0e7),(double)(ry / 1.0e7),(double)(rz / 1.0e7));
        if (fc < 0) {
          fprintf(stderr, "Can not write file %s.\n",filename);
          fclose(fp);
          fp = NULL;
        }
        fc = fprintf(fp," %f %f %f",(double)(rx - rx0),(double)(ry - ry0),(double)(rz - rz0));
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

