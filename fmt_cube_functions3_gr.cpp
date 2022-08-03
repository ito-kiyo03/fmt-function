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
double r1_global, r2_global, x1_global, x2_global, x3_global, y1_global, y2_global, y3_global, z1_global, z2_global, z3_global, aw_global[lenaw_global];
double r_x, r_y, r_z, x1_upperlimit, x2_upperlimit, x3_upperlimit, y1_upperlimit, y2_upperlimit, y3_upperlimit, z1_upperlimit, z2_upperlimit, z3_upperlimit;
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


/*double n_alpha_3D_naive(const struct _parameter prm, const int n)
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
                integrand_z[i2] = n_3(prm, r, rho);
            }
            integrand_y[i1] = trapezoidal_integration(n, integrand_z, delta);
        }
        integrand_x[i0] = trapezoidal_integration(n, integrand_y, delta);
    }
    val = trapezoidal_integration(n, integrand_x, delta);//ここまでで単位胞あたりの理想自由エネルギー
    return val / pow(prm.lambda, 3);// lambda^3 で割って単位体積当たりの値に変換
}*/

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

void n_alpha_divide_unit_cell(const struct _parameter prm, const int n, double rdiv[])
{
    for (int i = 0; i < n; i++)
    {
        rdiv[i] = -1.0 + double(i) / (n - 1.0);
    }
}

double rho_calculate(const struct _parameter prm, const int Nr, const double r_integ[], const double rho[][Nr][Nr])
{
    //cout << "caluculate rho  ";
    if (abs(r_integ[0]) <= 1.0 && abs(r_integ[1]) <= 1.0 && abs(r_integ[2]) <= 1.0)
    {
        return 0;
    }
    else
    {
        int ceilx = ceil(r_integ[0] / (2 * prm.lambda)), ceily = ceil(r_integ[1] / (2 * prm.lambda)), ceilz = ceil(r_integ[2] / (2 * prm.lambda));
        return rho[ceilx][ceily][ceilz];
    }
}

/*double n_alpha_rho(const struct _parameter prm, const double r_integrange[], const int Nr, const double r[], const double rho[][Nr][Nr])
{
    if (abs(r_integrange[0]) <= 1.0 && abs(r_integrange[1]) <= 1.0 && abs(r_integrange[2]) <= 1.0)
    {
        return 0;
    }
    else
    {
        int floorx = floor(r_integrange[0] / (2 * prm.lambda)), floory = floor(r_integrange[1] / (2 * prm.lambda)), floorz = floor(r_integrange[2] / (2 * prm.lambda));
        int ceilx = ceil(r_integrange[0] / (2 * prm.lambda)), ceily = ceil(r_integrange[1] / (2 * prm.lambda)), ceilz = ceil(r_integrange[2] / (2 * prm.lambda));
        double delta = sqrt(pow(floorx - ceilx, 2) + pow(floory - ceily, 2) + pow(floorz - ceilz, 2));
        return (rho[floorx][floory][floorz] - rho[ceilx][ceily][ceilz]) * delta / 2;
    }
}*/

double n_alpha_rho(const struct _parameter prm, const double r_integrange[], const int Nr, const double r[], const double rho[][Nr][Nr], const int n)
{
    cout << "caluculate n_alpha  ";
    if (abs(r_integrange[0]) <= 1.0 && abs(r_integrange[1]) <= 1.0 && abs(r_integrange[2]) <= 1.0)
    {
        return 0;
    }
    else
    {
        double r_integ[3], pmsigma[2], sigma[3], r_sigma[3], rdiv[2], delta;
        double integrand_x[n], integrand_y[n], integrand_z[n];
        pmsigma[0] = 1.0 / 2.0; pmsigma[1] = -1.0 / 2.0;//sigma=1.0として計算
        n_alpha_divide_unit_cell(prm, n, rdiv);
        delta = rdiv[1] - rdiv[0];
        for (int i = 0; i < 3; i++)
        {
            sigma[i] = pmsigma[0];
        }

        for (int integ_x = 0; integ_x < n; integ_x++)
        {
            r_integ[0] = r[0] - 0.5 + integ_x / (n - 1);
            for (int integ_y = 0; integ_y < n; integ_x++)
            {
                r_integ[1] = r[1] - 0.5 + integ_y / (n - 1);
                for (int integ_z = 0; integ_z < n; integ_z++)
                {
                    r_integ[2] = r[2] - 0.5 + integ_z / (n - 1);
                    integrand_z[integ_z] = rho_calculate(prm, Nr, r_integ, rho);
                }
                integrand_y[integ_y] = trapezoidal_integration(n, integrand_z, delta);
            }
            integrand_x[integ_x] = trapezoidal_integration(n, integrand_y, delta);
        }
        double val = trapezoidal_integration(n, integrand_x, delta);

        cout << endl;
        return val;
    }
}

double n_0(const struct _parameter prm, const double r_sigma[], const double sigma[], const int Nr, const double r[], const double rho[][Nr][Nr], const int n)
{
    cout << "caluculate n0  ";
    double n0 = 0.0, r_integrange[3];
    for (int i0 = 0; i0 < 2; i0++)
    {
        if (i0 == 0)
        {
            r_integrange[0] = r[0] + sigma[0];
        }
        else
        {
            r_integrange[0] = r[0] - sigma[0];
        }
        for (int i1 = 0; i1 < 2; i1++)
        {
            if (i1 == 0)
            {
                r_integrange[1] = r[1] + sigma[1];
            }
            else
            {
                r_integrange[1] = r[1] - sigma[1];
            }
            for (int i2 = 0; i2 < 2; i2++)
            {
                if (i2 == 0)
                {
                    r_integrange[2] = r[2] + sigma[2];
                }
                else
                {
                    r_integrange[2] = r[2] - sigma[2];
                }
                n0 += n_alpha_rho(prm, r_integrange, Nr, r, rho, n);
            }
        }
    }

    return n0 / 8;
}

double n_1x(const struct _parameter prm, const double r_sigma[], const double sigma[], const int Nr, const double r[], const double rho[][Nr][Nr], const int n)
{
    cout << "caluculate n1x  ";
    double rdiv[n], r_integrange[3], integrand_x[n], integrand_y[n], integrand_z[n], delta, val;
    n_alpha_divide_unit_cell(prm, n, rdiv);
    delta = rdiv[1] - rdiv[0];

    for (int i1 = 0; i1 < 2; i1++)
    {
        if (i1 == 0)
        {
            r_integrange[1] = r[1] + sigma[1];
        }
        else
        {
            r_integrange[1] = r[1] - sigma[1];
        }
        for (int i2 = 0; i2 < n; i2++)
        {
            if (i2 == 0)
            {
                r_integrange[2] = r[2] + sigma[2];
            }
            else
            {
                r_integrange[2] = r[2] - sigma[2];
            }
            for (int i0 = 0; i0 < n; i0++)
            {
                if (sigma[0] > 0)
                {
                    r_integrange[0] = r_sigma[0] + rdiv[i0];
                }
                else
                {
                    r_integrange[0] = r_sigma[0] + 0.5 - rdiv[i0];
                }
                integrand_z[i0] = n_alpha_rho(prm, r_integrange, Nr, r, rho, n);
            }
            integrand_y[i2] = trapezoidal_integration(n, integrand_z, delta);
        }
        integrand_x[i1] = trapezoidal_integration(n, integrand_y, delta);
    }
    val = trapezoidal_integration(n, integrand_x, delta);

    return val / 4;
}

double n_1y(const struct _parameter prm, const double r_sigma[], const double sigma[], const int Nr, const double r[], const double rho[][Nr][Nr], const int n)
{
    cout << "caluculate n1y  ";
    double rdiv[n], r_integrange[3], integrand_x[n], integrand_y[n], integrand_z[n], delta, val;
    n_alpha_divide_unit_cell(prm, n, rdiv);
    delta = rdiv[1] - rdiv[0];

    for (int i0 = 0; i0 < 2; i0++)
    {
        if (i0 == 0)
        {
            r_integrange[0] = r[0] + sigma[0];
        }
        else
        {
            r_integrange[0] = r[0] - sigma[0];
        }
        for (int i2 = 0; i2 < n; i2++)
        {
            if (i2 == 0)
            {
                r_integrange[2] = r[2] + sigma[2];
            }
            else
            {
                r_integrange[2] = r[2] - sigma[2];
            }
            for (int i1 = 0; i1 < n; i1++)
            {
                if (sigma[1] > 0)
                {
                    r_integrange[1] = r_sigma[1] + rdiv[i1];
                }
                else
                {
                    r_integrange[1] = r_sigma[1] + 0.5 - rdiv[i1];
                }
                integrand_y[i1] = n_alpha_rho(prm, r_integrange, Nr, r, rho, n);
            }
            integrand_z[i2] = trapezoidal_integration(n, integrand_z, delta);
        }
        integrand_x[i0] = trapezoidal_integration(n, integrand_y, delta);
    }
    val = trapezoidal_integration(n, integrand_x, delta);

    return val / 4;
}

double n_1z(const struct _parameter prm, const double r_sigma[], const double sigma[], const int Nr, const double r[], const double rho[][Nr][Nr], const int n)
{
    cout << "caluculate n1z  ";
    double rdiv[n], r_integrange[3], integrand_x[n], integrand_y[n], integrand_z[n], delta, val;
    n_alpha_divide_unit_cell(prm, n, rdiv);
    delta = rdiv[1] - rdiv[0];

    for (int i0 = 0; i0 < 2; i0++)
    {
        if (i0 == 0)
        {
            r_integrange[0] = r[0] + sigma[0];
        }
        else
        {
            r_integrange[0] = r[0] - sigma[0];
        }
        for (int i1 = 0; i1 < n; i1++)
        {
            if (sigma[1] > 0)
            {
                r_integrange[1] = r[1] + rdiv[i1];
            }
            else
            {
                r_integrange[1] = r[1]  + 0.5 - rdiv[i1];
            }
            for (int i2 = 0; i2 < n; i2++)
            {
                if (sigma[2] > 0)
                {
                    r_integrange[2] = r_sigma[2] + rdiv[i2];
                }
                else
                {
                    r_integrange[2] = r_sigma[2] + 0.5 - rdiv[i2];
                }
                integrand_z[i2] = n_alpha_rho(prm, r_integrange, Nr, r, rho, n);
            }
            integrand_y[i1] = trapezoidal_integration(n, integrand_z, delta);
        }
        integrand_x[i0] = trapezoidal_integration(n, integrand_y, delta);
    }
    val = trapezoidal_integration(n, integrand_x, delta);

    return val / 4;
}

double n_2x(const struct _parameter prm, const double r_sigma[], const double sigma[], const int Nr, const double r[], const double rho[][Nr][Nr], const int n)
{
    cout << "caluculate n2x  ";
    double rdiv[n], r_integrange[3], integrand_x[n], integrand_y[n], integrand_z[n], delta, val;
    n_alpha_divide_unit_cell(prm, n, rdiv);
    delta = rdiv[1] - rdiv[0];

    for (int i0 = 0; i0 < 2;i0++)
    {
        if (i0 == 0)
        {
            r_integrange[0] = r[0] + sigma[0];
        }
        else
        {
            r_integrange[0] = r[0] - sigma[0];
        }
        for (int i1 = 0; i1 < n; i1++)
        {
            if (sigma[1] > 0)
            {
                r_integrange[1] = r_sigma[1] + rdiv[i1];
            }
            else
            {
                r_integrange[1] = r_sigma[1] - rdiv[i1];
            }
            for (int i2 = 0; i2 < n; i2++)
            {
                if (sigma[2] > 0)
                {
                    r_integrange[2] = r_sigma[2] + rdiv[i2];
                }
                else
                {
                    r_integrange[2] = r_sigma[2] - rdiv[i2];
                }
                integrand_z[i2] = n_alpha_rho(prm, r_integrange, Nr, r, rho, n);
            }
            integrand_y[i1] = trapezoidal_integration(n, integrand_z, delta);
        }
        integrand_x[i0] = trapezoidal_integration(n, integrand_y, delta);
    }
    val = trapezoidal_integration(n, integrand_x, delta);
    
    return val / 2;
}

double n_2y(const struct _parameter prm, const double r_sigma[], const double sigma[], const int Nr, const double r[], const double rho[][Nr][Nr], const int n)
{
    cout << "caluculate n2y  ";
    double rdiv[n], r_integrange[3], integrand_x[n], integrand_y[n], integrand_z[n], delta, val;
    n_alpha_divide_unit_cell(prm, n, rdiv);
    delta = rdiv[1] - rdiv[0];

    for (int i1 = 0; i1 < 2; i1++)
    {
        if (i1 == 0)
        {
            r_integrange[1] = r[1] + sigma[1];
        }
        else
        {
            r_integrange[1] = r[1] - sigma[1];
        }
        for (int i0 = 0; i0 < n; i0++)
        {
            if (sigma[1] > 0)
            {
                r_integrange[0] = r_sigma[0] + rdiv[i0];
            }
            else
            {
                r_integrange[0] = r_sigma[0] - rdiv[i0];
            }
            for (int i2 = 0; i2 < n; i2++)
            {
                if (sigma[2] > 0)
                {
                    r_integrange[2] = r_sigma[2] + rdiv[i2];
                }
                else
                {
                    r_integrange[2] = r_sigma[2] - rdiv[i2];
                }
                integrand_z[i2] = n_alpha_rho(prm, r_integrange, Nr, r, rho, n);
            }
            integrand_x[i0] = trapezoidal_integration(n, integrand_z, delta);
        }
        integrand_y[i1] = trapezoidal_integration(n, integrand_x, delta);
    }
    val = trapezoidal_integration(n, integrand_y, delta);

    return val / 2;
}

double n_2z(const struct _parameter prm, const double r_sigma[], const double sigma[], const int Nr, const double r[], const double rho[][Nr][Nr], const int n)
{
    cout << "caluculate n2z  ";
    double rdiv[n], r_integrange[3], integrand_x[n], integrand_y[n], integrand_z[n], delta, val;
    n_alpha_divide_unit_cell(prm, n, rdiv);
    delta = rdiv[1] - rdiv[0];

    for (int i2 = 0; i2 < 2; i2++)
    {
        if (i2 == 0)
        {
            r_integrange[2] = r[2] + sigma[2];
        }
        else
        {
            r_integrange[2] = r[2] - sigma[2];
        }
        for (int i0 = 0; i0 < n; i0++)
        {
            if (sigma[0] > 0)
            {
                r_integrange[0] = r_sigma[0] + rdiv[i0];
            }
            else
            {
                r_integrange[0] = r_sigma[0] - rdiv[i0];
            }
            for (int i1 = 0; i1 < n; i1++)
            {
                if (sigma[1] > 0)
                {
                    r_integrange[1] = r_sigma[1] + rdiv[i1];
                }
                else
                {
                    r_integrange[1] = r_sigma[1] - rdiv[i1];
                }
                integrand_y[i1] = n_alpha_rho(prm, r_integrange, Nr, r, rho, n);
            }
            integrand_x[i0] = trapezoidal_integration(n, integrand_y, delta);
        }
        integrand_z[i2] = trapezoidal_integration(n, integrand_x, delta);
    }
    val = trapezoidal_integration(n, integrand_z, delta);
    
    return val / 2;
}

double n_3(const struct _parameter prm, const double r_sigma[], const double sigma[], const int Nr, const double r[], const double rho[][Nr][Nr], const int n)
{
    cout << "caluculate n3  ";
    double rdiv[n], r_integrange[3], integrand_x[n], integrand_y[n], integrand_z[n], delta, val;
    n_alpha_divide_unit_cell(prm, n, rdiv);
    delta = rdiv[1] - rdiv[0];
    for (int i0 = 0; i0 < n; i0++)
    {
        if (sigma[0] > 0)
        {
            r_integrange[0] = r_sigma[0] + rdiv[i0];
        }
        else
        {
            r_integrange[0] = r_sigma[0] + 0.5 - rdiv[i0];
        }
        for (int i1 = 0; i1 < n; i1++)
        {
            if (sigma[1] > 0)
            {
                r_integrange[1] = r_sigma[1] + rdiv[i1];
            }
            else
            {
                r_integrange[1] = r_sigma[1] + 0.5 -rdiv[i1];
            }
            for (int i2 = 0; i2 < n; i2++)
            {
                if (sigma[2] > 0)
                {
                    r_integrange[2] = r_sigma[2] + rdiv[i2];
                }
                else
                {
                    r_integrange[2] = r_sigma[2] + 0.5 - rdiv[i2];
                }
                integrand_z[i2] = n_alpha_rho(prm, r_integrange, Nr, r, rho, n);
            }
            integrand_y[i1] = trapezoidal_integration(n, integrand_z, delta);
        }
        integrand_x[i0] = trapezoidal_integration(n, integrand_y, delta);
    }
    val = trapezoidal_integration(n, integrand_x, delta);

    return val / pow(prm.lambda, 3);// lambda^3 で割って単位体積当たりの値に変換
}

double testPhi_n0_derivative(const struct _parameter prm, const int Nr, const double r[], const double rho[][Nr][Nr], const int n)
{
    cout << "caluculate n0_derivative  ";
    double n0_derivative = 0.0, pmsigma[2], sigma[3], r_sigma[3];
    pmsigma[0] = 1.0 / 2.0; pmsigma[1] = -1.0 / 2.0;//sigma=1.0として計算

    for (int i = 0; i < 2; i++)
    {
        r_sigma[0] = r[0] + pmsigma[i];
        sigma[0] = pmsigma[i];
        for (int j = 0; j < 2; j++)
        {
            r_sigma[1] = r[1] + pmsigma[j];
            sigma[1] = pmsigma[j];
            for (int k = 0; k < 2; k++)
            {
                r_sigma[2] = r[2] + pmsigma[k];
                sigma[2] = pmsigma[k];
                n0_derivative += log(1.0 - n_3(prm, r_sigma, sigma, Nr, r, rho, n));
            }
        }
    }

    return -n0_derivative / 8;
}

double testPhi_n1_derivative(const struct _parameter prm, const int Nr, const double r[], const double rho[][Nr][Nr], const int n)
{
    cout << "caluculate n1_derivative  ";
    double n1x_derivative = 0.0, n1y_derivative = 0.0, n1z_derivative = 0.0;
    double r_integ[3], pmsigma[2], sigma[3], r_sigma[3], rdiv[2], delta;
    double integrand_x[n], integrand_y[n], integrand_z[n];
    pmsigma[0] = 1.0 / 2.0; pmsigma[1] = -1.0 / 2.0;//sigma=1.0として計算
    n_alpha_divide_unit_cell(prm, n, rdiv);
    delta = rdiv[1] - rdiv[0];
    for (int i = 0; i < 3; i++)
    {
        sigma[i] = pmsigma[0];
    }

    for (int integ = 0; integ < n; integ++)
    {
        r_integ[0] = r[0] - 0.5 + integ / (n - 1);
        for (int i = 0; i < 2; i++)
        {
            r_sigma[1] = r[1] + pmsigma[i];
            sigma[1] = pmsigma[i];
            for (int j = 0; j < 2; j++)
            {
                r_sigma[2] = r[2] + pmsigma[j];
                sigma[2] = pmsigma[j];
                n1x_derivative += n_2x(prm, r_sigma, sigma, Nr, r, rho, n) / abs(1.0 - n_3(prm, r_sigma, sigma, Nr, r, rho, n));
            }
        }
        integrand_x[integ] = n1x_derivative;
    }
    double val_x = trapezoidal_integration(n, integrand_x, delta);

    for (int integ = 0; integ < n; integ++)
    {
        r_integ[1] = r[1] - 0.5 + integ / (n - 1);
        for (int i = 0; i < 2; i++)
        {
            r_sigma[0] = r[0] + pmsigma[i];
            sigma[0] = pmsigma[i];
            for (int j = 0; j < 2; j++)
            {
                r_sigma[2] = r[2] + pmsigma[j];
                sigma[2] = pmsigma[j];
                n1y_derivative += n_2y(prm, r_sigma, sigma, Nr, r, rho, n) / abs(1.0 - n_3(prm, r_sigma, sigma, Nr, r, rho, n));
            }
        }
        integrand_y[integ] = n1y_derivative;
    }
    double val_y = trapezoidal_integration(n, integrand_y, delta);

    for (int integ = 0; integ < n; integ++)
    {
        r_integ[2] = r[2] - 0.5 + integ / (n - 1);
        for (int i = 0; i < 2; i++)
        {
            r_sigma[0] = r[0] + pmsigma[i];
            sigma[0] = pmsigma[i];
            for (int j = 0; j < 2; j++)
            {
                r_sigma[1] = r[1] + pmsigma[j];
                sigma[1] = pmsigma[j];
                n1z_derivative += n_2z(prm, r_sigma, sigma, Nr, r, rho, n) / abs(1.0 - n_3(prm, r_sigma, sigma, Nr, r, rho, n));
            }
        }
        integrand_z[integ] = n1z_derivative;
    }
    double val_z = trapezoidal_integration(n, integrand_z, delta);

    return (val_x + val_y + val_z) / 4;
}

double testPhi_n2_derivative(const struct _parameter prm, const int Nr, const double r[], const double rho[][Nr][Nr], const int n)
{
    cout << "caluculate n2_derivative  ";
    double n2x_derivative = 0.0, n2y_derivative = 0.0, n2z_derivative = 0.0;
    double r_integ[3], pmsigma[2], sigma[3], r_sigma[3], rdiv[2], delta;
    double integrand_x[n], integrand_y[n], integrand_z[n];
    pmsigma[0] = 1.0 / 2.0; pmsigma[1] = -1.0 / 2.0;//sigma=1.0として計算
    n_alpha_divide_unit_cell(prm, n, rdiv);
    delta = rdiv[1] - rdiv[0];
    for (int i = 0; i < 3; i++)
    {
        sigma[i] = pmsigma[0];
    }

    for (int integ_y = 0; integ_y < n; integ_y++)
    {
        r_integ[1] = r[1] - 0.5 + integ_y / (n - 1);
        for (int integ_z = 0; integ_z < n; integ_z++)
        {
            r_integ[2] = r[2] - 0.5 + integ_z / (n - 1);
            for (int j = 0; j < 2; j++)
            {
                r_sigma[0] = r[0] + pmsigma[j];
                sigma[0] = pmsigma[j];
                n2x_derivative += n_1x(prm, r_sigma, sigma, Nr, r, rho, n) / abs(1.0 - n_3(prm, r_sigma, sigma, Nr, r, rho, n));
                n2x_derivative += n_2y(prm, r_sigma, sigma, Nr, r, rho, n) * n_2z(prm, r_sigma, sigma, Nr, r, rho, n) / pow(1.0 - n_3(prm, r_sigma, sigma, Nr, r, rho, n), 2);
            }
            integrand_z[integ_z] = n2x_derivative;
        }
        integrand_y[integ_y] = trapezoidal_integration(n, integrand_z, delta);
    }
    double val_x = trapezoidal_integration(n, integrand_y, delta);

    for (int integ_x = 0; integ_x < n; integ_x++)
    {
        r_integ[0] = r[0] - 0.5 + integ_x / (n - 1);
        for (int integ_z = 0; integ_z < n; integ_z++)
        {
            r_integ[2] = r[2] - 0.5 + integ_z / (n - 1);
            for (int j = 0; j < 2; j++)
            {
                r_sigma[1] = r[1] + pmsigma[j];
                sigma[1] = pmsigma[j];
                n2y_derivative += n_1y(prm, r_sigma, sigma, Nr, r, rho, n) / abs(1.0 - n_3(prm, r_sigma, sigma, Nr, r, rho, n));
                n2y_derivative += n_2x(prm, r_sigma, sigma, Nr, r, rho, n) * n_2z(prm, r_sigma, sigma, Nr, r, rho, n) / pow(1.0 - n_3(prm, r_sigma, sigma, Nr, r, rho, n), 2);
            }
            integrand_z[integ_z] = n2y_derivative;
        }
        integrand_x[integ_x] = trapezoidal_integration(n, integrand_z, delta);
    }
    double val_y = trapezoidal_integration(n, integrand_x, delta);

    for (int integ_x = 0; integ_x < n; integ_x++)
    {
        r_integ[0] = r[0] - 0.5 + integ_x / (n - 1);
        for (int integ_y = 0; integ_y < n; integ_y++)
        {
            r_integ[1] = r[1] - 0.5 + integ_y / (n - 1);
            for (int j = 0; j < 2; j++)
            {
                r_sigma[2] = r[2] + pmsigma[j];
                sigma[2] = pmsigma[j];
                n2z_derivative += n_1z(prm, r_sigma, sigma, Nr, r, rho, n) / abs(1.0 - n_3(prm, r_sigma, sigma, Nr, r, rho, n));
                n2z_derivative += n_2x(prm, r_sigma, sigma, Nr, r, rho, n) * n_2y(prm, r_sigma, sigma, Nr, r, rho, n) / pow(1.0 - n_3(prm, r_sigma, sigma, Nr, r, rho, n), 2);
            }
            integrand_y[integ_y] = n2z_derivative;
        }
        integrand_x[integ_x] = trapezoidal_integration(n, integrand_y, delta);
    }
    double val_z = trapezoidal_integration(n, integrand_x, delta);

    return (val_x + val_y + val_z) / 2;
}

double testPhi_n3_derivative(const struct _parameter prm, const int Nr, const double r[], const double rho[][Nr][Nr], const int n)
{
    cout << "caluculate n3_derivative  ";
    double n3_derivative = 0.0;
    double r_integ[3], pmsigma[2], sigma[3], r_sigma[3], rdiv[2], delta;
    double integrand_x[n], integrand_y[n], integrand_z[n];
    pmsigma[0] = 1.0 / 2.0; pmsigma[1] = -1.0 / 2.0;//sigma=1.0として計算
    n_alpha_divide_unit_cell(prm, n, rdiv);
    delta = rdiv[1] - rdiv[0];
    for (int i = 0; i < 3; i++)
    {
        sigma[i] = pmsigma[0];
    }

    for (int integ_x = 0; integ_x < n; integ_x++)
    {
        r_integ[0] = r[0] - 0.5 + integ_x / (n - 1);
        for (int integ_y = 0; integ_y < n; integ_x++)
        {
            r_integ[1] = r[1] - 0.5 + integ_y / (n - 1);
            for (int integ_z = 0; integ_z < n; integ_z++)
            {
                r_integ[2] = r[2] - 0.5 + integ_z / (n - 1);
                n3_derivative += n_0(prm, r_sigma, sigma, Nr, r, rho, n) / (1 - n_3(prm, r_sigma, sigma, Nr, r, rho, n));
                n3_derivative += (n_1x(prm, r_sigma, sigma, Nr, r, rho, n) * n_2x(prm, r_sigma, sigma, Nr, r, rho, n) + n_1y(prm, r_sigma, sigma, Nr, r, rho, n) * n_2y(prm, r_sigma, sigma, Nr, r, rho, n) + n_1z(prm, r_sigma, sigma, Nr, r, rho, n) * n_2z(prm, r_sigma, sigma, Nr, r, rho, n)) / pow(1.0 - n_3(prm, r_sigma, sigma, Nr, r, rho, n), 2);
                n3_derivative += n_2x(prm, r_sigma, sigma, Nr, r, rho, n) * n_2y(prm, r_sigma, sigma, Nr, r, rho, n) * n_2z(prm, r_sigma, sigma, Nr, r, rho, n) / pow(1.0 - n_3(prm, r_sigma, sigma, Nr, r, rho, n), 3);
                integrand_z[integ_z] = n3_derivative;
            }
            integrand_y[integ_y] = trapezoidal_integration(n, integrand_z, delta);
        }
        integrand_x[integ_x] = trapezoidal_integration(n, integrand_y, delta);
    }
    double val = trapezoidal_integration(n, integrand_x, delta);

    return val;
}

double testFex_density_3D_rho_derivative(const struct _parameter prm, const int Nr, const double r[], const double rho[][Nr][Nr], double betamu)
{
    extern _parameter prm_global;

    if (r[0] <= 1.0 && r[1] <= 1.0 && r[2] <= 1.0)
    {
        return 0;
    }
    else
    {
        double n0_derivative = testPhi_n0_derivative(prm, Nr, r, rho, 2);
        double n1_derivative = testPhi_n1_derivative(prm, Nr, r, rho, 2);
        double n2_derivative = testPhi_n2_derivative(prm, Nr, r, rho, 2);
        double n3_derivative = testPhi_n3_derivative(prm, Nr, r, rho, 2);

        double Fex_derivative = n0_derivative + n1_derivative + n2_derivative + n3_derivative;
        return exp(-Fex_derivative + betamu);
    }
}