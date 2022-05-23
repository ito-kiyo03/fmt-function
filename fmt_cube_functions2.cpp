#include <cmath>
#include <iostream>
#include<random>
#include <fstream>
#include <string>
#include <utility>
#include "fmt_cube2.h"
#include "intde2.h"
using namespace::std;
const int lenaw_global = 8000;
double r1_global, r2_global, aw_global[lenaw_global];
_parameter prm_global;

// 与えられた関数 f(prm, r[]) について、r[1] = given_r1 のときのf(prm, r[])= 0 の解 r[0] を range[0] から range[1] の範囲で求める関数。
// 評価座標での関数値の絶対値が eps 未満なら、その評価座標を解と見なす。
// 領域を二等分して関数値の符号が反転する領域を見つけ、その領域をさらに二等分して探すので、最初の3点（range[0]、range[1]とその中点）での関数値が
// 同符号なら解が見つからない。このときは nan を戻す
bool solve_f_eq_0_binary(const _parameter prm, const double given_r1, double (*f)(const _parameter prm2, const double r2[]), const double range[], const int max_i, const double eps, double& solution)
{
    bool found = false;
    double tmp_r[3][2], fval[3];
    for (int i1 = 0; i1 < 3; i1++)
    {
        tmp_r[i1][1] = given_r1;
    }
    tmp_r[0][0] = range[0];//領域左端
    tmp_r[2][0] = range[1];//領域右端
    tmp_r[1][0] = 0.5 * (tmp_r[0][0] + tmp_r[2][0]);//領域中央
    for (int i1 = 0; i1 < 3; i1++)
    {
        fval[i1] = f(prm, tmp_r[i1]);
    }
    for (int i1 = 0; i1 < max_i; i1++)
    {
        if (fval[0] * fval[1] <= 0)
        {
            found = true;
            tmp_r[2][0] = tmp_r[1][0];
            fval[2] = fval[1];
        }
        else if (fval[1] * fval[2] <= 0)
        {
            found = true;
            tmp_r[0][0] = tmp_r[1][0];
            fval[0] = fval[1];
        }
        else
        {
            break;
        }
        tmp_r[1][0] = 0.5 * (tmp_r[0][0] + tmp_r[2][0]);
        fval[1] = f(prm, tmp_r[1]);
        if (abs(fval[1]) < eps)
        {
            break;
        }
    }
    solution = tmp_r[1][0];
    return found;
}

// 与えられた関数 f(prm, r[]) について、r[1] = given_r1 のときのf(prm, r[])= 0 の解 r[0] を range[0] から range[1] の範囲で求める関数。
// 領域を Ndiv 等分して関数値の符号が反転する領域を見つけ、その領域を result_range[] に代入する。
// その領域が複数ある場合は、range[0] に最も近い領域だけが result_range[] に代入される。
bool solve_f_eq_0_naive(const _parameter prm, const double given_r1, double (*f)(const _parameter prm2, const double r2[]), const double range[], const int Ndiv, double result_range[])
{
    bool found = false;
    double tmp_r[2][2], fval[2], delta = (range[1] - range[0]) / Ndiv;
    result_range[0] = result_range[1] = nan("");
    for (int i = 0; i < 2; i++)
    {
        tmp_r[i][1] = given_r1;
    }
    tmp_r[0][0] = range[0];
    fval[0] = f(prm, tmp_r[0]);
    for (int i1 = 1; i1 < Ndiv; i1++)
    {
        tmp_r[1][0] = tmp_r[0][0] + delta;
        fval[1] = f(prm, tmp_r[1]);
        if (fval[0] * fval[1] <= 0)
        {
            found = true;
            result_range[0] = tmp_r[0][0];
            result_range[1] = tmp_r[1][0];
            break;
        }
        else
        {
            fval[0] = fval[1];
            tmp_r[0][0] = tmp_r[1][0];
        }
    }
    return found;
}

// 3次元でカットオフ距離（関数 cutoff_r_for_rho を参照）をニュートンラフソン法で求めるための補助関数
void get_delta_and_slope_3D(const struct _parameter prm, const double x, double& delta, double& slope)
{
    double y = prm.sqrt_gamma * x;
    delta = prm.eta * (2 * y * exp(-y * y) / sqrt(M_PI) + erfc(y));
    slope = -prm.eta * pow(prm.gamma / M_PI, 1.5) * 4 * M_PI * exp(-pow(y, 2)) * x * x;
}

// カットオフ距離、つまりこれより遠い全ての格子点を無視したとき、誤差が eps になる距離を戻す関数。
// 2次元では解析的な結果を戻すが、3次元ではニュートンラフソン法で求めた計算結果を戻す。
double cutoff_r_for_rho(const int dim, const double eps, const struct _parameter prm)
{
    double val = nan("");
    if (dim == 2)
    {
        val = sqrt(log(prm.eta / eps)) / prm.sqrt_gamma;
    }
    else if (dim == 3)
    {
        double delta, slope, x, newx = 1;
        for (int i = 0; i < 100; i++)
        {
            x = newx;
            get_delta_and_slope_3D(prm, x, delta, slope);
            if (delta < eps)
            {
                break;
            }
            newx = x + (eps - delta) / slope;
            if (abs(newx - x) < 1.e-8)
            {
                break;
            }
        }
        val = newx;
    }
    return val;
}

// 計算に必要なパラメータをまとめた構造体を設定する関数
// vacancy < 0 （格子点より粒子数が多い場合）は許容するが > 1 は許容しない（戻り値 err = true、つまりエラーが生じた扱いにする）
bool set_parameter(const int dim, const double eta, const double gamma, const double vacancy, const double eps, _parameter& prm)
{
    bool err = false;
    if (vacancy <= 1.0)
    {
        prm.eta = eta;
        prm.gamma = gamma;
        prm.vacancy = vacancy;
        prm.lambda = pow((1 - vacancy) / eta, 1.0 / dim);
        if (prm.lambda < 1.0)//粒子一辺より格子定数が小さいとき
        {
            err = true;
        }
        prm.coeff_rho = eta * pow(prm.lambda, dim) * pow(gamma / M_PI, dim / 2.0);// rho の係数
        prm.coeff_rhoSm = eta * prm.lambda * sqrt(gamma / M_PI);;// rhoSm の係数ηλ√(γ/π)
        prm.coeff_t = 0.5 * sqrt(M_PI / gamma);// t() の係数
        prm.sqrt_gamma = sqrt(gamma);// sqrt(gamma)
        // cutoff_r_for_rho より遠方は、 rho を計算するときに無視しても誤差 eps 以下
        prm.cutoff_r_for_rho = cutoff_r_for_rho(dim, eps, prm);
        prm.cutoff_m_for_rho = 1 + int(ceil(prm.cutoff_r_for_rho / prm.lambda));
    }
    else
    {
        err = true;
    }
    return err;
}

// 自由エネルギー最小 Fmin とその位置を記録するための配列 Fmin_at[] を初期化する関数
void refresh_F_minimum(const int i_vacancy, const int i_gamma, const double newF, int Fmin_at[2], double& Fmin)
{
    Fmin_at[0] = i_vacancy;
    Fmin_at[1] = i_gamma;
    Fmin = newF;
}

// これまでの自由エネルギー最小 Fmin と新たな自由エネルギー newF を比較し newF < Fmin なら Fmin 更新。
// 自由エネルギー最小の位置 Fmin_at[] も {i_vacancy, i_gamma} に更新。
void renew_F_minimum(const int i_vacancy, const int i_gamma, const double newF, int Fmin_at[2], double& Fmin)
{
    if (newF < Fmin)
    {
        Fmin_at[0] = i_vacancy;
        Fmin_at[1] = i_gamma;
        Fmin = newF;
    }
}

double t(const struct _parameter prm, const double x)
{
  return prm.coeff_t * (erf(prm.sqrt_gamma * (x + 0.5)) - erf(prm.sqrt_gamma * (x - 0.5)));
}

double T(const struct _parameter prm, const double x)
{
    int max_i = prm.cutoff_m_for_rho, ix = int(round(x / prm.lambda));
    double T = 0;
    for (int i = ix - max_i; i <= ix + max_i; i++)
    {
        T += t(prm, x - prm.lambda * i);
    }
    return T;
}

double zexp(const struct _parameter prm, const double x)
{
    return (exp(-prm.gamma * pow(x + 0.5, 2.0)) + exp(-prm.gamma * pow(x - 0.5, 2.0))) / 2.0;
}

double Zexp(const struct _parameter prm, const double x)
{
    int max_i = prm.cutoff_m_for_rho, ix = int(round(x / prm.lambda));
    double Z = 0;
    for (int i = ix - max_i; i <= ix + max_i; i++)
    {
        Z += zexp(prm, x - prm.lambda * i);
    }
    return Z;
}

// 単位胞 [- lambda/2, lambda/2] を n - 1 分割した分点 r[] を求める関数
void divide_unit_cell(const struct _parameter prm, const int n, double rdiv[])
{
    for (int i = 0; i < n; i++)
    {
        rdiv[i] = (-0.5 + double(i) / (n - 1)) * prm.lambda;
    }
}

// 点 r[] = (x, y) における二次元荷重関数
double n0_2D(const _parameter prm, const double r[])
{
    return prm.coeff_rho * Zexp(prm, r[0]) * Zexp(prm, r[1]);
}

// 点 r[] = (x, y) における二次元荷重関数
//ベクトル n1 のx成分のみ
double n1x_2D(const _parameter prm, const double r[])
{
    return prm.coeff_rho * Zexp(prm, r[0]) * T(prm, r[1]);
}

// 点 r[] = (x, y) における二次元荷重関数
double n2_2D(const _parameter prm, const double r[])
{
    return prm.coeff_rho * T(prm, r[0]) * T(prm, r[1]);
}

//二次元で、余剰自由エネルギー密度の特異点に関係する関数。この関数の零点で自由エネルギー密度は発散
double one_minus_n2_2D(const _parameter prm, const double r[])
{
    return 1 - n2_2D(prm, r);
}

// 点 r[] = (x, y, z) における3次元荷重関数
double n0_3D(const _parameter prm, const double r[])
{
    return prm.coeff_rho * Zexp(prm, r[0]) * Zexp(prm, r[1]) * Zexp(prm, r[2]);
}

// 点 r[] = (x, y, z) における3次元荷重関数
//ベクトル n1 のx成分のみ
double n1x_3D(const _parameter prm, const double r[])
{
    return prm.coeff_rho * T(prm, r[0]) * Zexp(prm, r[1]) * Zexp(prm, r[2]);
}

//ベクトル n1 のy成分のみ
double n1y_3D(const _parameter prm, const double r[])
{
    return prm.coeff_rho * Zexp(prm, r[0]) * T(prm, r[1]) * Zexp(prm, r[2]);
}

//ベクトル n1 のz成分のみ
double n1z_3D(const _parameter prm, const double r[])
{
    return prm.coeff_rho * Zexp(prm, r[0]) * Zexp(prm, r[1]) * T(prm, r[2]);
}

//ベクトルn2のx成分
double n2x_3D(const _parameter prm, const double r[])
{
    return prm.coeff_rho * Zexp(prm, r[0]) * T(prm, r[1]) * T(prm, r[2]);
}

//ベクトルn2のy成分
double n2y_3D(const _parameter prm, const double r[])
{
    return prm.coeff_rho * T(prm, r[0]) * Zexp(prm, r[1]) * T(prm, r[2]);
}

//ベクトルn2のz成分
double n2z_3D(const _parameter prm, const double r[])
{
    return prm.coeff_rho * T(prm, r[0]) * T(prm, r[1]) * Zexp(prm, r[2]);
}

//
double n3_3D(const _parameter prm, const double r[])
{
    return prm.coeff_rho * T(prm, r[0]) * T(prm, r[1]) * T(prm, r[2]);
}

//3次元で、余剰自由エネルギー密度の特異点に関係する関数。この関数の零点で自由エネルギー密度は発散
double one_minus_n3_3D(const _parameter prm, const double r[])
{
    return 1 - n3_3D(prm, r);
}

// 2次元。
// y 座標 given_r1 が与えられたとき、1 - n2 = 0 となる x 座標を戻す関数
// 解が見つからないときは nan を戻す。
// 2次元正方結晶（一辺 lambda）の対称性から、解となる x 座標を探す領域は [given_r1, lambda / 2] で十分である。
// まずこの領域を n 個の小領域に等分して関数 1 - n2(x, y) の符号が反転する小領域を探し、
// さらにその小領域を 2 等分しながら解が存在する領域を絞り込む作業を m 回だけ繰り返す。
// |1 - n2| < eps になったら終了。
double solution_n2_eq_1_2D(const _parameter prm, const double given_r1, const double eps, const int n, const int m)
{
    bool found = false;
    double solution, range[2] = { given_r1, 0.5 * prm.lambda }, result_range[2];
    found = solve_f_eq_0_naive(prm, given_r1, one_minus_n2_2D, range, n, result_range);
    if (found)
    {
        found = solve_f_eq_0_binary(prm, given_r1, one_minus_n2_2D, result_range, m, eps, solution);
    }
    else
    {
        solution = nan("");
    }
    return solution;
}

double solution_n3_eq_1_3D(const _parameter prm, const double given_r1, const double eps, const int n, const int m)
{
    bool found = false;
    double solution, range[2] = { given_r1, 0.5 * prm.lambda }, result_range[2];
    found = solve_f_eq_0_naive(prm, given_r1, one_minus_n3_3D, range, n, result_range);
    if (found)
    {
        found = solve_f_eq_0_binary(prm, given_r1, one_minus_n3_3D, result_range, m, eps, solution);
    }
    else
    {
        solution = nan("");
    }
    return solution;
}

// 2次元。
// 1 - n2(r[0], r[1]) = 0 が解を持つかどうか判定する関数。
// 解あり、および無しの場合、それぞれ true、および false を戻す
// 2次元正方結晶（一辺 lambda）の対称性から、調べる領域は 
// r[1] の範囲： [0, lambda / 2]
// r[0] の範囲： [r[1], lambda / 2]
// で十分である。
// まずこの r[1] の範囲を n - 1 個の小領域に等分する（小領域の幅 lambda / (2 * (n - 1))）。
// r[1] はその値に固定し、r[0] の範囲を同じ幅 lambda / (2 * (n - 1)) で分割し、 1 - n2 の符号反転を探すことで解の有無を調べる
bool solution_n2_eq_1_2D_exists(const _parameter prm, const int n)
{
    bool found = false;
    double r[2], result_range[2], delta = 0.5 * prm.lambda / (n - 1), range[2] = { 0.0, 0.5 * prm.lambda };
    for (int i1 = 0; i1 < n; i1++)
    {
        r[1] = i1 * delta;
        range[0] = r[1];
        found = solve_f_eq_0_naive(prm, r[1], one_minus_n2_2D, range, n - i1, result_range);
        if (found)
        {
            break;
        }
    }
    return found;
}

bool solution_n3_eq_1_3D_exists(const _parameter prm, const int n)
{
    bool found = false;
    double r[2], result_range[2], delta = 0.5 * prm.lambda / (n - 1), range[2] = { 0.0, 0.5 * prm.lambda };
    for (int i1 = 0; i1 < n; i1++)
    {
        r[1] = i1 * delta;
        range[0] = r[1];
        found = solve_f_eq_0_naive(prm, r[1], one_minus_n3_3D, range, n - i1, result_range);
        if (found)
        {
            break;
        }
    }
    return found;
}

// 点 r[] = (x, y) における2次元の余剰自由エネルギー密度
// 計算量を減らすため、定義した荷重関数 n0_2D, n1_2D, n2_2D を使っていない
double Phi_2D(const struct _parameter prm, const double r[])
{
    double Tval[2], Zval[2], n0, n1xn1y, n2;
    for (int i1 = 0; i1 < 2; i1++)
    {
        Tval[i1] = T(prm, r[i1]);
        Zval[i1] = Zexp(prm, r[i1]);
    }
    n0 = prm.coeff_rho * Zval[0] * Zval[1];
    n2 = prm.coeff_rho * Tval[0] * Tval[1];
    n1xn1y = n0 * n2;
    return -n0 * log(abs(1 - n2)) + n1xn1y / abs(1 - n2);
}

// 点 r[] = (x, y) における2次元の余剰自由エネルギー密度
// 1次元積分のルーチン intde2 を用いる都合上、パラメータおよび y 座標にそれぞれグローバル変数 prm_global および r1_global を利用して
// 1 変数関数として定義する。
double Phi_2D(double r0)
{
    extern double r1_global;
    extern _parameter prm_global;
    double Tval[2], Zval[2], n0, n1xn1y, n2, r[2] = { r0, r1_global };
    for (int i1 = 0; i1 < 2; i1++)
    {
        Tval[i1] = T(prm_global, r[i1]);
        Zval[i1] = Zexp(prm_global, r[i1]);
    }
    n0 = prm_global.coeff_rho * Zval[0] * Zval[1];
    n2 = prm_global.coeff_rho * Tval[0] * Tval[1];
    n1xn1y = n0 * n2;
    return -n0 * log(abs(1 - n2)) + n1xn1y / abs(1 - n2);
}

// 点 r[] = (x, y, z) における3次元の余剰自由エネルギー密度
// 計算量を減らすため、定義した荷重関数 n0_3D, n1_3D, n2_3D を使っていない
double Phi_3D(const struct _parameter prm, const double r[])
{
    double Tval[3], Zval[3], n0, n1[3], n2[3], n1n2 = 0.0, n2_all = 1.0, n3;
    for (int i1 = 0; i1 < 3; i1++)
    {
        Tval[i1] = T(prm, r[i1]);
        Zval[i1] = Zexp(prm, r[i1]);
    }
    n0 = prm.coeff_rho * Zval[0] * Zval[1] * Zval[2];
    n1[0] = prm.coeff_rho * Tval[0] * Zval[1] * Zval[2];
    n1[1] = prm.coeff_rho * Zval[0] * Tval[1] * Zval[2];
    n1[2] = prm.coeff_rho * Zval[0] * Zval[1] * Tval[2];
    n2[0] = prm.coeff_rho * Zval[0] * Tval[1] * Tval[2];
    n2[1] = prm.coeff_rho * Tval[0] * Zval[1] * Tval[2];
    n2[2] = prm.coeff_rho * Tval[0] * Tval[1] * Zval[2];
    n3 = prm.coeff_rho * Tval[0] * Tval[1] * Tval[2];
    for (int i = 0; i < 3; i++)
    {
        n1n2 += n1[i] * n2[i];
        n2_all *= n2[i];
    }
    return -n0 * log(abs(1.0 - n3)) + n1n2 / abs(1.0 - n3) + n2_all / pow(1.0 - n3, 2);
}

// 点 r[] = (x, y, z) における3次元の余剰自由エネルギー密度
// 1次元積分のルーチン intde2 を用いる都合上、パラメータおよび y 座標にそれぞれグローバル変数 prm_global および r2_global を利用して
// 1 変数関数として定義する。
double Phi_3D(double r0)
{
    extern double r1_global, r2_global;
    extern _parameter prm_global;
    double Tval[3], Zval[3], n0, n1[3], n2[3], n1n2 = 0.0, n2_all = 1.0, n3, r[3] = { r0, r1_global, r2_global };
    for (int i1 = 0; i1 < 3; i1++)
    {
        Tval[i1] = T(prm_global, r[i1]);
        Zval[i1] = Zexp(prm_global, r[i1]);
    }
    n0 = prm_global.coeff_rho * Zval[0] * Zval[1] * Zval[2];
    n1[0] = prm_global.coeff_rho * Tval[0] * Zval[1] * Zval[2];
    n1[1] = prm_global.coeff_rho * Zval[0] * Tval[1] * Zval[2];
    n1[2] = prm_global.coeff_rho * Zval[0] * Zval[1] * Tval[2];
    n2[0] = prm_global.coeff_rho * Zval[0] * Tval[1] * Tval[2];
    n2[1] = prm_global.coeff_rho * Tval[0] * Zval[1] * Tval[2];
    n2[2] = prm_global.coeff_rho * Tval[0] * Tval[1] * Zval[2];
    n3 = prm_global.coeff_rho * Tval[0] * Tval[1] * Tval[2];
    for (int i = 0; i < 3; i++)
    {
        n1n2 += n1[i] * n2[i];
        n2_all *= n2[i];
    }
    return -n0 * log(abs(1.0 - n3)) + n1n2 / abs(1.0 - n3) + n2_all / pow(1.0 - n3, 2);
}

// 点 r[] = (x, y) における密度 rho を求める関数
double rho_2D(const struct _parameter prm, const double r[])
{
    int max_i = prm.cutoff_m_for_rho;
    double rho[2] = {0, 0};

    for (int i = 0; i < 2; i++)
    {
        int ir = round(r[i] / prm.lambda);
        for (int m = -max_i + ir; m <= max_i + ir; m++)
        {
            rho[i] += exp(-prm.gamma * pow(r[i] - prm.lambda * m, 2));
        }
    }
    //return prm.eta;//d=0
    return prm.coeff_rho * rho[0] * rho[1];//d=2
}

// 点 r[] = (x, y, z) における密度 rho を求める関数
double rho_3D(const struct _parameter prm, const double r[])
{
    int max_i = prm.cutoff_m_for_rho;
    double rho[3] = { 0, 0 ,0 };

    for (int i = 0; i < 3; i++)
    {
        int ir = round(r[i] / prm.lambda);
        for (int m = -max_i + ir; m <= max_i + ir; m++)
        {
            rho[i] += exp(-prm.gamma * pow(r[i] - prm.lambda * m, 2));
        }
    }
    //return prm.eta;//d=0
    return prm.coeff_rho * rho[0] * rho[1] * rho[2];//d=2
}

// 点 r[] = (x, y) における理想自由エネルギー密度
double Fid_density_2D(const struct _parameter prm, const double r[])
{
    double rho = rho_2D(prm, r);
    return rho * (log(rho) - 1);
}

// 点 r[] = (x, y, z) における理想自由エネルギー密度
double Fid_density_3D(const struct _parameter prm, const double r[])
{
    double rho = rho_3D(prm, r);
    return rho * (log(rho) - 1);
}

// 1次元台形公式で配列の積分を求める
double trapezoidal_integration(const int n, const double integrand[], const double delta)
{
    double val = 0;
    for (int i = 0; i < n; i++)
    {
        val += integrand[i];
    }
    return delta * (val - 0.5 * (integrand[0] + integrand[n - 1]));
}

// 2次元の単位体積あたりの理想自由エネルギー。
// 台形公式で求める。
double Fid_density_2D(const struct _parameter prm, const int n)
{
    double rdiv[n], r[2], integrand_x[n], integrand_y[n], delta, val;
    divide_unit_cell(prm, n, rdiv);
    delta = rdiv[1] - rdiv[0];
    for (int i0 = 0; i0 < n; i0++)
    {
        r[0] = rdiv[i0];
        for (int i1 = 0; i1 < n; i1++)
        {
            r[1] = rdiv[i1];
            integrand_y[i1] = Fid_density_2D(prm, r);
        }
        integrand_x[i0] = trapezoidal_integration(n, integrand_y, delta);
    }
    val = trapezoidal_integration(n, integrand_x, delta);//ここまでで単位胞あたりの理想自由エネルギー
    return val / pow(prm.lambda, 2);// lambda^2 で割って単位体積当たりの値に変換
}

// 3次元の単位体積あたりの理想自由エネルギー。
// 台形公式で求める。
double Fid_density_3D(const struct _parameter prm, const int n)
{
    double rdiv[n], r[3], integrand_x[n], integrand_y[n], integrand_z[n], delta, val;
    divide_unit_cell(prm, n, rdiv);
    delta = rdiv[1] - rdiv[0];
    for (int i0 = 0; i0 < n; i0++)
    {
        r[0] = rdiv[i0];
        for (int i1 = 0; i1 < n; i1++)
        {
            r[1] = rdiv[i1];
            for (int i2 = 0; i2 < n; i2++)
            {
                r[2] = rdiv[i2];
                integrand_z[i2] = Fid_density_3D(prm, r);
            }
            integrand_y[i1] = trapezoidal_integration(n, integrand_z, delta);
        }
        integrand_x[i0] = trapezoidal_integration(n, integrand_y, delta);
    }
    val = trapezoidal_integration(n, integrand_x, delta);//ここまでで単位胞あたりの理想自由エネルギー
    return val / pow(prm.lambda, 3);// lambda^3 で割って単位体積当たりの値に変換
}

//一粒子当たりの理想自由エネルギー。
//gamma が大きい場合の近似
double Fid_density_approx(const int dim, const struct _parameter prm)
{
    double Fid_density = nan("");
    if (dim == 2)
    {
        Fid_density = prm.eta * pow(prm.lambda, 2) * (log(prm.eta * pow(prm.lambda, 2) * prm.gamma / M_PI) - 2);
    }
    else if (dim == 3)
    {
        Fid_density = 0.5 * prm.eta * pow(prm.lambda, 3) * (2 * log(prm.eta * pow(prm.lambda, 3) * pow(prm.gamma / M_PI, 1.5)) - 5);
    }//ここまでで単位胞あたりの理想自由エネルギー
    return Fid_density / pow(prm.lambda, dim);//単位体積あたりの値に変換
}

bool isnan_array1(const int n, const double arry[])
{
    bool rslt = false;
    for (int i = 0; i < n; i++)
    {
        if (isnan(arry[i]))
        {
            rslt = true;
            break;
        }
    }
    return rslt;
}

// この関数は、次の関数 r0range_n2_eq_1_2D 内で使うことを想定している。
// 2次元の余剰自由エネルギー密度 Phi_2D の発散に寄与する項（1 - n2 を含む項）の絶対値が、微小量 epsilon より小さいとき true、大きいとき false を戻す関数。
// ただし、|1 - n2| > eps_n2 である場合に限る。
// このとき |1 - n2| > eps_n2 なので  |(ln|1 - n2|)| < |ln(eps_n2)|、および 1/|1 - n2| < 1/eps_n2。
// したがって n0 < |ln(eps_n2)| * epsilon、かつ n1x * n1y < epsilon /eps_n2 なら Phi_2D の各項は epsilon 未満
bool Phi_2D_terms_are_small(const _parameter prm, const double r[], const double eps_n2, const double epsilon)
{
    bool val;
    double n0 = n0_2D(prm, r), n1x = n1x_2D(prm, r);
    double rflip[2] = { r[1], r[0] }, n1y = n1x_2D(prm, rflip);
    if ((n0 < epsilon * abs(log(eps_n2))) && (n1x * n1y < epsilon / eps_n2))
    {
        val = true;
    }
    else
    {
        val = false;
    }
    return val;
}

bool Phi_3D_terms_are_small(const _parameter prm, const double r[], const double eps_n3, const double epsilon)
{
    bool val;
    double n0 = n0_3D(prm, r), n1[3], n2[3], n1n2 = 0.0, n2_all = 1.0;

    n1[0] = n1x_3D(prm, r);
    n1[1] = n1y_3D(prm, r);
    n1[2] = n1z_3D(prm, r);
    n2[0] = n2x_3D(prm, r);
    n2[1] = n2y_3D(prm, r);
    n2[2] = n2z_3D(prm, r);

    for (int i = 0; i < 3; i++)
    {
        n1n2 += n1[i] * n2[i];
        n2_all *= n2[i];
    }

    if ((n0 < epsilon * abs(log(eps_n3))) && (n1n2 < epsilon / eps_n3) && (n2_all < epsilon / eps_n3))
    {
        val = true;
    }
    else
    {
        val = false;
    }
    return val;
}

// パラメータが prm であり、y座標が r1 の場合の x 領域 [y, 0.5 * prm.lambda] のうち、 
// abs(1 - n2) > eps_n2 となる x の領域 r0range[2] を求める関数。この領域で Phi_2D が数値計算可能と見なせるように eps_n2 を選ぶ。
// Phi_2D が積分可能であるためには、r0range[2] 外で Phi_2D が零であることが必要だが、
// その確認は abs(n2 - 1) ~= 0 となる側の r0range[2] の端点で、Phi_2D が十分に小さいことで代用する。
// abs(n2 - 1) > eps_n2 となる x の領域が [r0min, r0max] のように連続した1つの範囲で表せて、かつ端点で十分に小さいときは戻り値 true、
// [r0min1, r0max1] および [r0min2, r0max2] のように複数になる場合または端点で小さくない場合は戻り値 false
bool r0range_n2_eq_1_2D(const _parameter prm, const double r1, double r0range[])
{
    bool val = true, found_lower = false, found_higher = false;
    int Ndiv = 1000;
    double delta = (0.5 * prm.lambda - r1) / (Ndiv - 1), eps_n2 = 1.0e-13, epsilon = eps_n2;
    double r[2];
    r[1] = r1;
    for (int ix = 0; ix < Ndiv; ix++)
    {
        r[0] = r1 + ix * delta;
        if (abs(one_minus_n2_2D(prm, r)) > eps_n2)
        {
            if (found_higher) // r0range 上端が確定したのにふたたび |1 - n2| > eps_n2 となったら、戻り値 false
            {
                val = false;
                break;
            }
            else {
                if (!found_lower) // r0range 上端が確定しておらず r0range下端が確定していない場合、この r[0] が r0range 下端
                {
                    r0range[0] = r[0];
                    r0range[1] = r[1]; // r0range 下端と上端が一致する場合に備えて上端を設定
                    found_lower = true; // r0range 下端が確定した
                }
                else //r0range下端が見つかっている場合、この r[0] が r0range 上端の候補。さらに大きい r[0] で |1 - n2| > eps_n2 となる場合は更新する
                {
                    r0range[1] = r[0];
                }
            }
        }
        else
        {
            if (found_lower) // r0range下端が見つかっていて |1 - n2| < eps_n2 なら、上端が確定したことになる
            {
                found_higher = true;
            }
        }
    }
    if ((found_lower) && (!found_higher)) // r0range 下端が確定していて上端が確定していない場合は lambda / 2 が r0range 上端
    {
        r0range[1] = 0.5 * prm.lambda;
        found_higher = true;
    }
    if (val) // この時点で val = false なら、ここは飛ばして val = false を戻す
    {
        r[0] = r1;
        if (abs(one_minus_n2_2D(prm, r)) < eps_n2)//この場合、r0range[0] で Phi_2D が十分に小さいかどうか確認が必要
        {
            r[0] = r0range[0];
            val = Phi_2D_terms_are_small(prm, r, eps_n2, epsilon);
        }
        if (val)
        {
            r[0] = 0.5 * prm.lambda;
            if (abs(one_minus_n2_2D(prm, r)) < eps_n2)//この場合、r0range[1] で Phi_2D が十分に小さいかどうか確認が必要
            {
                r[0] = r0range[1];
                val = Phi_2D_terms_are_small(prm, r, eps_n2, epsilon);
            }
        }
    }
    return val;
}

bool r0range_n3_eq_1_3D(const _parameter prm, const double r1, double r0range[])
{
    bool val = true, found_lower = false, found_higher = false;
    int Ndiv = 1000;
    double delta = (0.5 * prm.lambda - r1) / (Ndiv - 1), eps_n2 = 1.0e-13, epsilon = eps_n2;
    double r[2];
    r[1] = r1;
    for (int ix = 0; ix < Ndiv; ix++)
    {
        r[0] = r1 + ix * delta;
        if (abs(one_minus_n3_3D(prm, r)) > eps_n2)
        {
            if (found_higher) // r0range 上端が確定したのにふたたび |1 - n3| > eps_n2 となったら、戻り値 false
            {
                val = false;
                break;
            }
            else {
                if (!found_lower) // r0range 上端が確定しておらず r0range下端が確定していない場合、この r[0] が r0range 下端
                {
                    r0range[0] = r[0];
                    r0range[1] = r[1]; // r0range 下端と上端が一致する場合に備えて上端を設定
                    found_lower = true; // r0range 下端が確定した
                }
                else //r0range下端が見つかっている場合、この r[0] が r0range 上端の候補。さらに大きい r[0] で |1 - n2| > eps_n2 となる場合は更新する
                {
                    r0range[1] = r[0];
                }
            }
        }
        else
        {
            if (found_lower) // r0range下端が見つかっていて |1 - n3| < eps_n2 なら、上端が確定したことになる
            {
                found_higher = true;
            }
        }
    }
    if ((found_lower) && (!found_higher)) // r0range 下端が確定していて上端が確定していない場合は lambda / 2 が r0range 上端
    {
        r0range[1] = 0.5 * prm.lambda;
        found_higher = true;
    }
    if (val) // この時点で val = false なら、ここは飛ばして val = false を戻す
    {
        r[0] = r1;
        if (abs(one_minus_n3_3D(prm, r)) < eps_n2)//この場合、r0range[0] で Phi_3D が十分に小さいかどうか確認が必要
        {
            r[0] = r0range[0];
            val = Phi_3D_terms_are_small(prm, r, eps_n2, epsilon);
        }
        if (val)
        {
            r[0] = 0.5 * prm.lambda;
            if (abs(one_minus_n3_3D(prm, r)) < eps_n2)//この場合、r0range[1] で Phi_3D が十分に小さいかどうか確認が必要
            {
                r[0] = r0range[1];
                val = Phi_3D_terms_are_small(prm, r, eps_n2, epsilon);
            }
        }
    }
    return val;
}

// パラメータが prm_global であり、y座標が r1 の場合の、Phi_2D の x 積分を求める関数。
// 積分に1変数積分ルーチン intde2 を使うので、被積分関数を 1 変数積分にするため prm_global, r1_global はグローバル変数にする。
double integrate_Phi_2D_wrt_x(const double r1)
{
    extern const int lenaw_global;
    extern double r1_global, aw_global[lenaw_global];
    extern _parameter prm_global;
    const int Ndiv = 100, Mdiv = 50;
    double integ_result[2] = { 0.0, 0.0 }, integ_err[2] = { 0.0, 0.0 }, r0_range[2] = { r1, 0.5 * prm_global.lambda };
    double solution, val = nan("");
    bool integrable;
    r1_global = r1;
    // n2 = 1 の解があるかどうか探す
    solution = solution_n2_eq_1_2D(prm_global, r1, 1.e-13, Ndiv, Mdiv);
    if (isnan(solution)) // n2 = 1 の解がない場合、通常の積分で問題ない
    {
        // r0 積分の範囲は r0_range[0] から r0_range[1]。(prm_global.lambda)
        intde(Phi_2D, r0_range[0], r0_range[1], aw_global, &integ_result[0], &integ_err[0]);
        val = integ_result[0];
    }
    else // n2 = 1 の解 solution がある場合、r0 積分領域を solution で分割して和を取る
    {
        intde(Phi_2D, r0_range[0], solution - 1.0e-8, aw_global, &integ_result[0], &integ_err[0]);
        intde(Phi_2D, solution + 1.0e-8, r0_range[1], aw_global, &integ_result[1], &integ_err[1]);
        if (integ_err[0] < 0 || integ_err[1] < 0 || integ_err[0] > 1.0 || integ_err[1] > 1.0) // 分割した積分が異常値の場合
        {
            integrable = r0range_n2_eq_1_2D(prm_global, r1, r0_range);
            // n2 ~= 1 となる領域のすぐ外側で Phi の各項が十分に 0 に近いとき、その領域からの積分の寄与は 0 とみなし、n2 != 1 の領域だけで積分
            if (integrable)
            {
                intde(Phi_2D, r0_range[0], r0_range[1], aw_global, &integ_result[0], &integ_err[0]);
                val = integ_result[0];
            }
            else
            {
                val = nan("");
            }
        }
        else // 分割した積分が正常値の場合
        {
            val = integ_result[0] + integ_result[1];
        }
    }
    return val;
}

double integrate_Phi_3D_wrt_yz(const double r2)
{
    extern const int lenaw_global;
    extern double r2_global, aw_global[lenaw_global];
    extern _parameter prm_global;
    const int Ndiv = 100, Mdiv = 50;
    double integ_result[2] = { 0.0, 0.0 }, integ_err[2] = { 0.0, 0.0 }, r2_range[2] = { r2, 0.5 * prm_global.lambda };
    double solution, val = nan("");
    bool integrable;
    r2_global = r2;
    // n3 = 1 の解があるかどうか探す
    solution = solution_n3_eq_1_3D(prm_global, r2, 1.e-13, Ndiv, Mdiv);
    if (isnan(solution)) // n3 = 1 の解がない場合、通常の積分で問題ない
    {
        // r0 積分の範囲は r0_range[0] から r0_range[1]。(prm_global.lambda)
        intde(Phi_3D, r2_range[0], r2_range[1], aw_global, &integ_result[0], &integ_err[0]);
        val = integ_result[0];
    }
    else // n3 = 1 の解 solution がある場合、r0 積分領域を solution で分割して和を取る
    {
        intde(Phi_3D, r2_range[0], solution - 1.0e-8, aw_global, &integ_result[0], &integ_err[0]);
        intde(Phi_3D, solution + 1.0e-8, r2_range[1], aw_global, &integ_result[1], &integ_err[1]);
        if (integ_err[0] < 0 || integ_err[1] < 0 || integ_err[0] > 1.0 || integ_err[1] > 1.0) // 分割した積分が異常値の場合
        {
            integrable = r0range_n3_eq_1_3D(prm_global, r2, r2_range);
            // n3 ~= 1 となる領域のすぐ外側で Phi の各項が十分に 0 に近いとき、その領域からの積分の寄与は 0 とみなし、n3 != 1 の領域だけで積分
            if (integrable)
            {
                intde(Phi_3D, r2_range[0], r2_range[1], aw_global, &integ_result[0], &integ_err[0]);
                val = integ_result[0];
            }
            else
            {
                val = nan("");
            }
        }
        else // 分割した積分が正常値の場合
        {
            val = integ_result[0] + integ_result[1];
        }
    }
    return val;
}

double integrate_Phi_3D_wrt_xy(const double r1)
{
    extern const int lenaw_global;
    extern double r1_global, aw_global[lenaw_global];
    extern _parameter prm_global;
    const int Ndiv = 100, Mdiv = 50;
    double integ_result[2] = { 0.0, 0.0 }, integ_err[2] = { 0.0, 0.0 }, r1_range[2] = { r1, 0.5 * prm_global.lambda };
    double solution, val = nan("");
    bool integrable;
    r1_global = r1;
    intde(integrate_Phi_3D_wrt_yz, r1_range[0], r1_range[1], aw_global, &integ_result[0], &integ_err[0]);
    // r1 積分の範囲は r1_range[0] から r1_range[1]。(prm_global.lambda)
    //intde(Phi_3D, r1_range[0], r1_range[1], aw_global, &integ_result[0], &integ_err[0]);
    val = integ_result[0];

    return val;
}

/*double integrate_Phi_3D_wrt_xy(const double r1)
{
    extern const int lenaw_global;
    extern double r1_global, aw_global[lenaw_global];
    extern _parameter prm_global;
    const int Ndiv = 100, Mdiv = 50;
    double integ_result[2] = { 0.0, 0.0 }, integ_err[2] = { 0.0, 0.0 }, r1_range[2] = { r1, 0.5 * prm_global.lambda };
    double solution, val = nan("");
    bool integrable;
    r1_global = r1;
    intde(integrate_Phi_3D_wrt_yz, r1_range[0], r1_range[1], aw_global, &integ_result[0], &integ_err[0]);
    // n3 = 1 の解があるかどうか探す
    solution = solution_n3_eq_1_3D(prm_global, r1, 1.e-13, Ndiv, Mdiv);
    if (isnan(solution)) // n3 = 1 の解がない場合、通常の積分で問題ない
    {
        // r0 積分の範囲は r0_range[0] から r0_range[1]。(prm_global.lambda)
        intde(Phi_3D, r1_range[0], r1_range[1], aw_global, &integ_result[0], &integ_err[0]);
        val = integ_result[0];
    }
    else // n3 = 1 の解 solution がある場合、r0 積分領域を solution で分割して和を取る
    {
        intde(Phi_3D, r1_range[0], solution - 1.0e-8, aw_global, &integ_result[0], &integ_err[0]);
        intde(Phi_3D, solution + 1.0e-8, r1_range[1], aw_global, &integ_result[1], &integ_err[1]);
        if (integ_err[0] < 0 || integ_err[1] < 0 || integ_err[0] > 1.0 || integ_err[1] > 1.0) // 分割した積分が異常値の場合
        {
            integrable = r0range_n3_eq_1_3D(prm_global, r1, r1_range);
            // n3 ~= 1 となる領域のすぐ外側で Phi の各項が十分に 0 に近いとき、その領域からの積分の寄与は 0 とみなし、n3 != 1 の領域だけで積分
            if (integrable)
            {
                intde(Phi_3D, r1_range[0], r1_range[1], aw_global, &integ_result[0], &integ_err[0]);
                val = integ_result[0];
            }
            else
            {
                val = nan("");
            }
        }
        else // 分割した積分が正常値の場合
        {
            val = integ_result[0] + integ_result[1];
        }
    }
    return val;
}*/

double Fex_density_2D(const struct _parameter prm, const double rel_error_req)
{
  extern const int lenaw_global;
  extern double aw_global[lenaw_global];
  extern _parameter prm_global;
  const double tiny = 1.e-307;
  double integ_result, integ_err, r0_range[2] = {0.0, 0.5 * prm.lambda};
  prm_global = prm;
  intdeini(lenaw_global, tiny, rel_error_req, aw_global);
  intde(integrate_Phi_2D_wrt_x, r0_range[0], r0_range[1], aw_global, &integ_result, &integ_err);
  // ここまでで、integ_result は単位胞の 1/8 あたりの余剰自由エネルギー
  return 8 * integ_result / pow(prm.lambda, 2);//単位体積あたりの余剰自由エネルギーに変換
}

double Fex_density_3D(const struct _parameter prm, const double rel_error_req)
{
    extern const int lenaw_global;
    extern double aw_global[lenaw_global];
    extern _parameter prm_global;
    const double tiny = 1.e-307;
    double integ_result, integ_err, r1_range[2] = { 0.0, 0.5 * prm.lambda };
    prm_global = prm;
    intdeini(lenaw_global, tiny, rel_error_req, aw_global);
    intde(integrate_Phi_3D_wrt_xy, r1_range[0], r1_range[1], aw_global, &integ_result, &integ_err);
    // ここまでで、integ_result は単位胞の 1/8 あたりの余剰自由エネルギー
    return 32 * integ_result / pow(prm.lambda, 3);//単位体積あたりの余剰自由エネルギーに変換
}

void F_density_dD(const int dim, const struct _parameter prm, const double rel_error_req, double F_density_set[])
{
    //理想自由エネルギー
    if (dim == 2)
    {
        if (prm.cutoff_r_for_rho < 0.5)
        {
            F_density_set[0] = Fid_density_approx(dim, prm);//近似
        }
        else
        {
            F_density_set[0] = Fid_density_2D(prm, 100);
        }
        //余剰自由エネルギー
        F_density_set[1] = Fex_density_2D(prm, rel_error_req);
    }
    else if (dim == 3)
    {
        if (prm.cutoff_r_for_rho < 0.5)
        {
            F_density_set[0] = Fid_density_approx(dim, prm);//近似
        }
        else
        {
            F_density_set[0] = Fid_density_3D(prm, 100);
        }
        //余剰自由エネルギー
        F_density_set[1] = Fex_density_3D(prm, rel_error_req);
    }
    else 
    {
        cout << "false";
    }
    //id+ex
    F_density_set[2] = F_density_set[0] + F_density_set[1];
}

// 2次元の単位体積あたりの余剰自由エネルギー。
// Phi_2D の特異点がないものと仮定して台形公式で求める。
double Fex_density_2D_naive(const struct _parameter prm, const int n)
{
    double rdiv[n], r[2], integrand_x[n], integrand_y[n], delta, val;
    divide_unit_cell(prm, n, rdiv);
    delta = rdiv[1] - rdiv[0];
    for (int i0 = 0; i0 < n; i0++)
    {
        r[0] = rdiv[i0];
        for (int i1 = 0; i1 < n; i1++)
        {
            r[1] = rdiv[i1];
            integrand_y[i1] = Phi_2D(prm, r);
        }
        integrand_x[i0] = trapezoidal_integration(n, integrand_y, delta);
    }
    val = trapezoidal_integration(n, integrand_x, delta);//ここまでで単位胞あたりの理想自由エネルギー
    return val / pow(prm.lambda, 2);// lambda^2 で割って単位体積当たりの値に変換
}

double Fex_density_3D_naive(const struct _parameter prm, const int n)
{
    double rdiv[n], r[3], integrand_x[n], integrand_y[n], integrand_z[n], delta, val;
    divide_unit_cell(prm, n, rdiv);
    delta = rdiv[1] - rdiv[0];
    for (int i0 = 0; i0 < n; i0++)
    {
        r[0] = rdiv[i0];
        for (int i1 = 0; i1 < n; i1++)
        {
            r[1] = rdiv[i1];
            for (int i2 = 0; i2 < n; i2++)
            {
                r[2] = rdiv[i2];
                integrand_z[i2] = Phi_3D(prm, r);
            }
            integrand_y[i1] = trapezoidal_integration(n, integrand_z, delta);
        }
        integrand_x[i0] = trapezoidal_integration(n, integrand_y, delta);
    }
    val = trapezoidal_integration(n, integrand_x, delta);//ここまでで単位胞あたりの理想自由エネルギー
    return val / pow(prm.lambda, 3);// lambda^3 で割って単位体積当たりの値に変換
}

void write_setting(ofstream& ofs, const int dim, const int Neta, const int Nvacancy, const int Ngamma, const double eps, const double eta_set[], const double vacancy_set[], const double gamma_set[])
{
    ofs << "dim: " << dim << ", Neta: " << Neta << ", Nvacancy: " << Nvacancy << ", Ngamma: " << Ngamma << endl;
    ofs << "eps: " << eps << endl;
    ofs << "eta_set: " << eta_set[0] << ", " << eta_set[1] << endl;
    ofs << "vacancy_set: " << vacancy_set[0] << ", " << vacancy_set[1] << endl;
    ofs << "gamma_set: " << gamma_set[0] << ", " << gamma_set[1] << endl;
}

// ofs で指定されるファイルに、いくつかのパラメータを出力する。
// write_title が true なら見出しも出力。
// pre_string に与えられた文字を行頭に書く。
// たとえば gnuplot で処理するファイルに出力する場合は pre_string = "#" とすると gnuplot にとってのコメントアウトになる。
void write_parameter(ofstream& ofs, const _parameter prm, const bool write_title, const string pre_string)
{
    if (write_title)
    {
        ofs << pre_string << "eta, gamma, vacancy, lambda, coeff_rho, cutoff_r_for_rho" << endl;
    }
    ofs << pre_string << prm.eta << ", " << prm.gamma << ", " << prm.vacancy << ", " << prm.lambda << ", " << prm.coeff_rho << " , " << prm.cutoff_r_for_rho << endl;
}

// 密度プロファイル、荷重密度、Phi を、xy 平面を各軸 Ndiv - 1 個に離散化して描く
// 第一象限だけ描く場合は quarter = true にする
void write_profile_2D(ofstream& ofs, const _parameter prm, const int Ndiv, const bool write_title, const bool quarter, const string pre_string)
{
    double r[2], rdiv[Ndiv];
    // 描画範囲
    double range[2] = { -0.5 * prm.lambda, 0.5 * prm.lambda };
    if (quarter)
    {
        range[0] = 0.0;
    }
    // xy 平面離散化
    for (int i = 0; i < Ndiv; i++)
    {
        rdiv[i] = range[0] + (range[1] - range[0]) * i / (Ndiv - 1);
    }
    write_parameter(ofs, prm, write_title, pre_string);
    for (int i = 0; i < Ndiv; i++)
    {
        r[0] = rdiv[i];
        for (int j = 0; j < Ndiv; j++)
        {
            r[1] = rdiv[j];
            ofs << r[0] << " " << r[1] << " " << rho_2D(prm, r) << endl;// << " " << n0_2D(prm, r) << " " << n1x_2D(prm, r) << " " << n2_2D(prm, r) << " " << Phi_2D(prm, r) << endl;
        }
        ofs << endl;
    }
    ofs << endl;
}

double n0_2DSm(const _parameter prm, const double r[])
{
    return prm.coeff_rhoSm * Zexp(prm, r[0]);
}

double n1x_2DSm(const _parameter prm, const double r[])
{
    return prm.coeff_rhoSm * Zexp(prm, r[0]);
}

double n1y_2DSm(const _parameter prm, const double r[])
{
    return prm.coeff_rhoSm * T(prm, r[0]);
}

double n2_2DSm(const _parameter prm, const double r[])
{
    return prm.coeff_rhoSm * T(prm, r[0]);
}


double one_minus_n2_2DSm(const _parameter prm, const double r[])
{
    return 1 - n2_2DSm(prm, r);
}

double Phi_2DSm(const struct _parameter prm, const double r[])
{
    double Tval[2], Zval[2], n0, n1xn1y, n2;
    for (int i1 = 0; i1 < 2; i1++)
    {
        Tval[i1] = T(prm, r[i1]);
        Zval[i1] = Zexp(prm, r[i1]);
    }
    n0 = prm.coeff_rhoSm * Zval[0];
    n2 = prm.coeff_rhoSm * Tval[0];
    n1xn1y = n0 * n2;
}

double Phi_2DSm(double r0)
{
    extern double r1_global;
    extern _parameter prm_global;
    double Tval[2], Zval[2], n0, n1xn1y, n2, r[2] = { r0, r1_global };
    for (int i1 = 0; i1 < 2; i1++)
    {
        Tval[i1] = T(prm_global, r[i1]);
        Zval[i1] = Zexp(prm_global, r[i1]);
    }
    n0 = prm_global.coeff_rhoSm * Zval[0];
    n2 = prm_global.coeff_rhoSm * Tval[0];
    n1xn1y = n0 * n2;
    return -n0 * log(abs(1 - n2)) + n1xn1y / abs(1 - n2);
}

double rho_2DSm(const struct _parameter prm, const double r[])
{
    int max_i = prm.cutoff_m_for_rho;
    double rho = 0;

    int ir = round(r[0] / prm.lambda);
    for (int m = -max_i + ir; m <= max_i + ir; m++)
    {
        rho += exp(-prm.gamma * pow(r[0] - prm.lambda * m, 2));
    }

    return prm.coeff_rhoSm * rho * prm.lambda;//
}

void write_profile_3D(ofstream& ofs, const _parameter prm, const int Ndiv, const bool write_title, const bool quarter, const string pre_string)
{
    double r[3], rdiv[Ndiv];
    // 描画範囲
    double range[2] = { -0.5 * prm.lambda, 0.5 * prm.lambda };
    if (quarter)
    {
        range[0] = 0.0;
    }
    // xy 平面離散化
    for (int i = 0; i < Ndiv; i++)
    {
        rdiv[i] = range[0] + (range[1] - range[0]) * i / (Ndiv - 1);
    }
    write_parameter(ofs, prm, write_title, pre_string);
    for (int i = 0; i < Ndiv; i++)
    {
        r[0] = 0.0;
        r[1] = rdiv[i];
        for (int j = 0; j < Ndiv; j++)
        {
            r[2] = rdiv[j];
            for (int k = 0; k < Ndiv; k++)
            {
                r[2] = rdiv[k];
                ofs << r[0] << " " << r[1] << " " << r[2] << " " << Phi_3D(prm, r) << endl;
            }
            ofs << endl;
        }
        ofs << endl;
    }
    ofs << endl;
}

void write_profile_3Dna(ofstream& ofs, const _parameter prm, const int Ndiv, const bool write_title, const bool quarter, const string pre_string)
{
    double r[3], rdiv[Ndiv];
    // 描画範囲
    double range[2] = { -0.5 * prm.lambda, 0.5 * prm.lambda };
    if (quarter)
    {
        range[0] = 0.0;
    }
    // xy 平面離散化
    for (int i = 0; i < Ndiv; i++)
    {
        rdiv[i] = range[0] + (range[1] - range[0]) * i / (Ndiv - 1);
    }
    write_parameter(ofs, prm, write_title, pre_string);
    for (int i = 0; i < Ndiv; i++)
    {
        r[2] = 0.0;
        r[0] = rdiv[i];
        for (int j = 0; j < Ndiv; j++)
        {
            r[1] = rdiv[j];
            double n0 = n0_3D(prm, r);
            double n1[3], n2[3], n1n2 = 0.0, n2_all = 1.0;
            double n3 = n3_3D(prm, r);

            n1[0] = n1x_3D(prm, r);
            n1[1] = n1y_3D(prm, r);
            n1[2] = n1z_3D(prm, r);
            n2[0] = n2x_3D(prm, r);
            n2[1] = n2y_3D(prm, r);
            n2[2] = n2z_3D(prm, r);
            for (int a = 0; a < 3; a++)
            {
                n1n2 += n1[a] * n2[a];
                n2_all *= n2[a];
            }
            ofs << r[0] << " " << r[1] << " " << r[2] << " " << n0 << " " << n1n2 << " " << n2_all << " " << n3;
            /*for (int k = 0; k < Ndiv; k++)
            {
                r[2] = rdiv[k];
                double n0 = n0_3D(prm, r);
                double n1[3], n2[3], n1n2 = 0.0, n2_all = 1.0;
                double n3 = n3_3D(prm, r);

                n1[0] = n1x_3D(prm, r);
                n1[1] = n1y_3D(prm, r);
                n1[2] = n1z_3D(prm, r);
                n2[0] = n2x_3D(prm, r);
                n2[1] = n2y_3D(prm, r);
                n2[2] = n2z_3D(prm, r);
                for (int a = 0; a < 3; a++)
                {
                    n1n2 += n1[a] * n2[a];
                    n2_all *= n2[a];
                }
                //ofs << r[0] << " " << r[1] << " " << r[2] << " " << n0 << " " << n1n2 << " " << n2_all << " " << n3 << endl;
            }*/
            ofs << endl;
        }
        ofs << endl;
    }
    ofs << endl;
}

void write_profile_3DPhi(ofstream& ofs, const _parameter prm, const int Ndiv, const bool write_title, const bool quarter, const string pre_string)
{
    double r[3], rdiv[Ndiv];
    // 描画範囲
    double range[2] = { -0.5 * prm.lambda, 0.5 * prm.lambda };
    if (quarter)
    {
        range[0] = 0.0;
    }
    // xy 平面離散化
    for (int i = 0; i < Ndiv; i++)
    {
        rdiv[i] = range[0] + (range[1] - range[0]) * i / (Ndiv - 1);
    }
    write_parameter(ofs, prm, write_title, pre_string);
    for (int i = 0; i < Ndiv; i++)
    {
        r[2] = 0.0;
        r[0] = rdiv[i];
        for (int j = 0; j < Ndiv; j++)
        {
            r[1] = rdiv[j];
            double n0 = n0_3D(prm, r);
            double n1[3], n2[3], n1n2 = 0.0, n2_all = 1.0;
            double n3 = n3_3D(prm, r);

            n1[0] = n1x_3D(prm, r);
            n1[1] = n1y_3D(prm, r);
            n1[2] = n1z_3D(prm, r);
            n2[0] = n2x_3D(prm, r);
            n2[1] = n2y_3D(prm, r);
            n2[2] = n2z_3D(prm, r);
            for (int a = 0; a < 3; a++)
            {
                n1n2 += n1[a] * n2[a];
                n2_all *= n2[a];
            }
            ofs << r[0] << " " << r[1] << " " << r[2] << " " << n0 * log(abs(1.0 - n3)) << " " << n1n2 / abs(1.0 - n3) << " " << n2_all / pow(1 - n3, 3);
            /*for (int k = 0; k < Ndiv; k++)
            {
                r[2] = rdiv[k];
                double n0 = n0_3D(prm, r);
                double n1[3], n2[3], n1n2 = 0.0, n2_all = 1.0;
                double n3 = n3_3D(prm, r);

                n1[0] = n1x_3D(prm, r);
                n1[1] = n1y_3D(prm, r);
                n1[2] = n1z_3D(prm, r);
                n2[0] = n2x_3D(prm, r);
                n2[1] = n2y_3D(prm, r);
                n2[2] = n2z_3D(prm, r);
                for (int a = 0; a < 3; a++)
                {
                    n1n2 += n1[a] * n2[a];
                    n2_all *= n2[a];
                }
                //ofs << r[0] << " " << r[1] << " " << n0 * log(abs(1.0 - n3)) << " " << n1n2 / abs(1.0 - n3) << " " << n2_all / pow(1 - n3, 3) << endl;
            }*/
            ofs << endl;
        }
        ofs << endl;
    }
    ofs << endl;
}

void write_profile_3Drho(ofstream& ofs, const _parameter prm, const int Ndiv, const bool write_title, const bool quarter, const string pre_string)
{
    double r[3], rdiv[Ndiv];
    // 描画範囲
    double range[2] = { -0.5 * prm.lambda, 0.5 * prm.lambda };
    if (quarter)
    {
        range[0] = 0.0;
    }
    // xy 平面離散化
    for (int i = 0; i < Ndiv; i++)
    {
        rdiv[i] = range[0] + (range[1] - range[0]) * i / (Ndiv - 1);
    }
    write_parameter(ofs, prm, write_title, pre_string);
    for (int i = 0; i < Ndiv; i++)
    {
        r[2] = 0.0;
        r[0] = rdiv[i];
        for (int j = 0; j < Ndiv; j++)
        {
            r[1] = rdiv[j];
            ofs << r[0] << " " << r[1] << " " << r[2] << " " << rho_3D(prm, r);
            /*for (int k = 0; k < Ndiv; k++)
            {
                r[2] = rdiv[k];
                ofs << r[0] << " " << r[1] << " " << r[2] << " " << rho_3D(prm, r) << endl;
            }*/
            ofs << endl;
        }
        ofs << endl;
    }
    ofs << endl;
}

void write_profile_2DSm(ofstream& ofs, const _parameter prm, const int Ndiv, const bool write_title, const bool quarter, const string pre_string)
{
    double r[2], rdiv[Ndiv];
    // 描画範囲
    double range[2] = { -0.5 * prm.lambda, 0.5 * prm.lambda };
    if (quarter)
    {
        range[0] = 0.0;
    }
    // xy 平面離散化
    for (int i = 0; i < Ndiv; i++)
    {
        rdiv[i] = range[0] + (range[1] - range[0]) * i / (Ndiv - 1);
    }
    write_parameter(ofs, prm, write_title, pre_string);
    for (int i = 0; i < Ndiv; i++)
    {
        r[0] = rdiv[i];
        for (int j = 0; j < Ndiv; j++)
        {
            r[1] = rdiv[j];
            double n0 = n0_2DSm(prm, r);
            double n1x = n1x_2DSm(prm, r), n1y = n2_2DSm(prm, r);
            double n2 = n2_2DSm(prm, r);
            ofs << r[0] << " " << r[1] << " " << n0 << " " << n1x * n1y << " " << n2 << endl;
        }
        ofs << endl;
    }
    ofs << endl;
}

void write_profile_2DSmn2(ofstream& ofs, const _parameter prm, const int Ndiv, const bool write_title, const bool quarter, const string pre_string)
{
    double r[2], rdiv[Ndiv];
    // 描画範囲
    double range[2] = { -0.5 * prm.lambda, 0.5 * prm.lambda };
    if (quarter)
    {
        range[0] = 0.0;
    }
    // xy 平面離散化
    for (int i = 0; i < Ndiv; i++)
    {
        rdiv[i] = range[0] + (range[1] - range[0]) * i / (Ndiv - 1);
    }
    write_parameter(ofs, prm, write_title, pre_string);
    for (int i = 0; i < Ndiv; i++)
    {
        r[0] = rdiv[i];
        for (int j = 0; j < Ndiv; j++)
        {
            r[1] = rdiv[j];
            double n0 = n0_2DSm(prm, r);
            double n1x = n1x_2DSm(prm, r), n1y = n1y_2DSm(prm, r);
            double n2 = n2_2DSm(prm, r);
            ofs << r[0] << " " << r[1] << " " << n0 * log(abs(1.0 - n2)) << " " << n1x * n1y / (1.0 - n2) << endl;
        }
        ofs << endl;
    }
    ofs << endl;
}

void write_profile_2DrhoSm(ofstream& ofs, const _parameter prm, const int Ndiv, const bool write_title, const bool quarter, const string pre_string)
{
    double r[2], rdiv[Ndiv];
    // 描画範囲
    double range[2] = { -0.5 * prm.lambda, 0.5 * prm.lambda };
    if (quarter)
    {
        range[0] = 0.0;
    }
    // xy 平面離散化
    for (int i = 0; i < Ndiv; i++)
    {
        rdiv[i] = range[0] + (range[1] - range[0]) * i / (Ndiv - 1);
    }
    write_parameter(ofs, prm, write_title, pre_string);
    for (int i = 0; i < Ndiv; i++)
    {
        r[0] = rdiv[i];
        for (int j = 0; j < Ndiv; j++)
        {
            r[1] = rdiv[j];
            ofs << r[0] << " " << r[1] << " " << rho_2DSm(prm, r) << endl;
        }
        ofs << endl;
    }
    ofs << endl;
}

double Fid_density_2DSm(const struct _parameter prm, const double r[])
{
    double rho = rho_2DSm(prm, r);
    return rho * (log(rho) - 1);
}

double Fid_density_2DSm(const struct _parameter prm, const int n)
{
    double rdiv[n], r[2], integrand_x[n], integrand_y[n], delta, val;
    divide_unit_cell(prm, n, rdiv);
    delta = rdiv[1] - rdiv[0];
    for (int i0 = 0; i0 < n; i0++)
    {
        r[0] = rdiv[i0];
        for (int i1 = 0; i1 < n; i1++)
        {
            r[1] = rdiv[i1];
            integrand_y[i1] = Fid_density_2DSm(prm, r);
        }
        integrand_x[i0] = trapezoidal_integration(n, integrand_y, delta);
    }
    val = trapezoidal_integration(n, integrand_x, delta);//ここまでで単位胞あたりの理想自由エネルギー
    return val / pow(prm.lambda, 2);// lambda^2 で割って単位体積当たりの値に変換
}

bool Phi_2DSm_terms_are_small(const _parameter prm, const double r[], const double eps_n2, const double epsilon)
{
    bool val;
    double n0 = n0_2DSm(prm, r), n1x = n1x_2DSm(prm, r), n1y = n2_2DSm(prm, r);
    if ((n0 < epsilon * abs(log(eps_n2))) && (n1x * n1y < epsilon / eps_n2))
    {
        val = true;
    }
    else
    {
        val = false;
    }
    return val;
}

double solution_n2_eq_1_2DSm(const _parameter prm, const double given_r1, const double eps, const int n, const int m)
{
    bool found = false;
    double solution, range[2] = { given_r1, 0.5 * prm.lambda }, result_range[2];
    found = solve_f_eq_0_naive(prm, given_r1, one_minus_n2_2DSm, range, n, result_range);
    if (found)
    {
        found = solve_f_eq_0_binary(prm, given_r1, one_minus_n2_2DSm, result_range, m, eps, solution);
    }
    else
    {
        solution = nan("");
    }
    return solution;
}

bool solution_n2_eq_1_2DSm_exists(const _parameter prm, const int n)
{
    bool found = false;
    double r[2], result_range[2], delta = 0.5 * prm.lambda / (n - 1), range[2] = { 0.0, 0.5 * prm.lambda };
    for (int i1 = 0; i1 < n; i1++)
    {
        r[1] = i1 * delta;
        range[0] = r[1];
        found = solve_f_eq_0_naive(prm, r[1], one_minus_n2_2DSm, range, n - i1, result_range);
        if (found)
        {
            break;
        }
    }
    return found;
}

bool r0range_n2_eq_1_2DSm(const _parameter prm, const double r1, double r0range[])
{
    bool val = true, found_lower = false, found_higher = false;
    int Ndiv = 1000;
    double delta = (0.5 * prm.lambda - r1) / (Ndiv - 1), eps_n2 = 1.0e-13, epsilon = eps_n2;
    double r[2];
    r[1] = r1;
    for (int ix = 0; ix < Ndiv; ix++)
    {
        r[0] = r1 + ix * delta;
        if (abs(one_minus_n2_2DSm(prm, r)) > eps_n2)
        {
            if (found_higher) // r0range 上端が確定したのにふたたび |1 - n2| > eps_n2 となったら、戻り値 false
            {
                val = false;
                break;
            }
            else {
                if (!found_lower) // r0range 上端が確定しておらず r0range下端が確定していない場合、この r[0] が r0range 下端
                {
                    r0range[0] = r[0];
                    r0range[1] = r[1]; // r0range 下端と上端が一致する場合に備えて上端を設定
                    found_lower = true; // r0range 下端が確定した
                }
                else //r0range下端が見つかっている場合、この r[0] が r0range 上端の候補。さらに大きい r[0] で |1 - n2| > eps_n2 となる場合は更新する
                {
                    r0range[1] = r[0];
                }
            }
        }
        else
        {
            if (found_lower) // r0range下端が見つかっていて |1 - n2| < eps_n2 なら、上端が確定したことになる
            {
                found_higher = true;
            }
        }
    }
    if ((found_lower) && (!found_higher)) // r0range 下端が確定していて上端が確定していない場合は lambda / 2 が r0range 上端
    {
        r0range[1] = 0.5 * prm.lambda;
        found_higher = true;
    }
    if (val) // この時点で val = false なら、ここは飛ばして val = false を戻す
    {
        r[0] = r1;
        if (abs(one_minus_n2_2DSm(prm, r)) < eps_n2)//この場合、r0range[0] で Phi_2DSm が十分に小さいかどうか確認が必要
        {
            r[0] = r0range[0];
            val = Phi_2DSm_terms_are_small(prm, r, eps_n2, epsilon);
        }
        if (val)
        {
            r[0] = 0.5 * prm.lambda;
            if (abs(one_minus_n2_2DSm(prm, r)) < eps_n2)//この場合、r0range[1] で Phi_2DSm が十分に小さいかどうか確認が必要
            {
                r[0] = r0range[1];
                val = Phi_2DSm_terms_are_small(prm, r, eps_n2, epsilon);
            }
        }
    }
    return val;
}

double integrate_Phi_2DSm_wrt_x(const double r1)
{
    extern const int lenaw_global;
    extern double r1_global, aw_global[lenaw_global];
    extern _parameter prm_global;
    const int Ndiv = 100, Mdiv = 50;
    double integ_result[2] = { 0.0, 0.0 }, integ_err[2] = { 0.0, 0.0 }, r0_range[2] = { r1, 0.5 * prm_global.lambda };
    double solution, val = nan("");
    bool integrable;
    r1_global = r1;
    // n2 = 1 の解があるかどうか探す
    solution = solution_n2_eq_1_2DSm(prm_global, r1, 1.e-13, Ndiv, Mdiv);
    if (isnan(solution)) // n2 = 1 の解がない場合、通常の積分で問題ない
    {
        // r0 積分の範囲は r0_range[0] から r0_range[1]。(prm_global.lambda)
        intde(Phi_2DSm, r0_range[0], r0_range[1], aw_global, &integ_result[0], &integ_err[0]);
        val = integ_result[0];
    }
    else // n2 = 1 の解 solution がある場合、r0 積分領域を solution で分割して和を取る
    {
        intde(Phi_2DSm, r0_range[0], solution - 1.0e-8, aw_global, &integ_result[0], &integ_err[0]);
        intde(Phi_2DSm, solution + 1.0e-8, r0_range[1], aw_global, &integ_result[1], &integ_err[1]);
        if (integ_err[0] < 0 || integ_err[1] < 0 || integ_err[0] > 1.0 || integ_err[1] > 1.0) // 分割した積分が異常値の場合
        {
            integrable = r0range_n2_eq_1_2D(prm_global, r1, r0_range);
            // n2 ~= 1 となる領域のすぐ外側で Phi の各項が十分に 0 に近いとき、その領域からの積分の寄与は 0 とみなし、n2 != 1 の領域だけで積分
            if (integrable)
            {
                intde(Phi_2DSm, r0_range[0], r0_range[1], aw_global, &integ_result[0], &integ_err[0]);
                val = integ_result[0];
            }
            else
            {
                val = nan("");
            }
        }
        else // 分割した積分が正常値の場合
        {
            val = integ_result[0] + integ_result[1];
        }
    }
    return val;
}

double Fex_density_2DSm(const struct _parameter prm, const double rel_error_req)
{
    extern const int lenaw_global;
    extern double aw_global[lenaw_global];
    extern _parameter prm_global;
    const double tiny = 1.e-307;
    double integ_result, integ_err, r1_range[2] = { 0.0, 0.5 * prm.lambda };
    prm_global = prm;
    intdeini(lenaw_global, tiny, rel_error_req, aw_global);
    intde(integrate_Phi_2DSm_wrt_x, r1_range[0], r1_range[1], aw_global, &integ_result, &integ_err);
    // ここまでで、integ_result は単位胞の 1/8 あたりの余剰自由エネルギー
    return 8 * integ_result / pow(prm.lambda, 2);//単位体積あたりの余剰自由エネルギーに変換
}

double Fid_density_approxSm(const int dim, const struct _parameter prm)
{
    double Fid_density = nan("");
    if (dim == 2)
    {
        Fid_density = prm.eta * pow(prm.lambda, 2.0) * (log(prm.eta * prm.lambda * sqrt(prm.gamma / M_PI)) - 3.0 / 2.0);
    }
    else if (dim == 3)
    {
        Fid_density = 0.5 * prm.eta * pow(prm.lambda, 3) * (2 * log(prm.eta * pow(prm.lambda, 3) * pow(prm.gamma / M_PI, 1.5)) - 5);
    }//ここまでで単位胞あたりの理想自由エネルギー
    return Fid_density / pow(prm.lambda, dim);//単位体積あたりの値に変換
}

void F_density_dDSm(const int dim, const struct _parameter prm, const double rel_error_req, double F_density_set[])
{
    //理想自由エネルギー
    /*if (prm.cutoff_r_for_rho < 0.5)
    {
        F_density_set[0] = Fid_density_approxSm(2, prm);//近似
    }
    else
    {*/
    F_density_set[0] = Fid_density_2DSm(prm, 100);
    //}
    //余剰自由エネルギー
    F_density_set[1] = Fex_density_2DSm(prm, rel_error_req);
    //id+ex
    F_density_set[2] = F_density_set[0] + F_density_set[1];
}

// 2次元の単位体積あたりの余剰自由エネルギー。
// Phi_2D の特異点がないものと仮定して台形公式で求める。
double Fex_density_2DSm_naive(const struct _parameter prm, const int n)
{
    double rdiv[n], r[2], integrand_x[n], integrand_y[n], delta, val;
    divide_unit_cell(prm, n, rdiv);
    delta = rdiv[1] - rdiv[0];
    for (int i0 = 0; i0 < n; i0++)
    {
        r[0] = rdiv[i0];
        for (int i1 = 0; i1 < n; i1++)
        {
            r[1] = rdiv[i1];
            integrand_y[i1] = Phi_2DSm(prm, r);
        }
        integrand_x[i0] = trapezoidal_integration(n, integrand_y, delta);
    }
    val = trapezoidal_integration(n, integrand_x, delta);//ここまでで単位胞あたりの理想自由エネルギー
    return val / pow(prm.lambda, 2);// lambda^2 で割って単位体積当たりの値に変換
}
