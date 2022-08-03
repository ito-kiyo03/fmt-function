#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
//#include <string>
#include "fmt_cube2.h"
//#include "intde2.h"
using namespace::std;

const int Nr = 5;

int main()
{
    const int dim = 3, Neta = 1, Ngamma = 1, Nvacancy = 10, Ndiv = 21;//いくつ出力するか
    const double eps = 1.0e-8, eta_set[2] = { 0.5, 0.1 }, gamma_set[2] = { 1.00, 1.0 }, vacancy_set[2] = { 0.000001, 0.1 };//{初期値 ,刻み}

    bool err;
    int Fmin_at[Neta][2];
    double r[3], F_density_set[3], testF_density_set[3], Fmin[Neta];//0:id 1:ex 2:id+ex
    struct _parameter prm;

    ofstream f_gr("gr,fmt_cube3,3D.dat");
    //ofstream f_profile("gr,fmt_cube3_profile,3D.dat");
    //ofstream f_setting("gr,fmt_cube3_setting,3D.dat");
    int d = 2;//L:0,Sm:1,X:2
    double rho[Nr][Nr][Nr], G[2][Nr][Nr][Nr], alpha = 0.9, betamu;

    //write_setting(f_setting, dim, Neta, Nvacancy, Ngamma, eps, eta_set, vacancy_set, gamma_set);
    for (int i_eta = 0; i_eta < Neta; i_eta++)//eta
    {
        double eta = eta_set[0] + eta_set[1] * i_eta;
        betamu = log(eta * (1.0 - eta)) - 7 * eta / (1.0 - eta) - 6 * pow(eta, 2) / pow(1.0 - eta, 2) - 2 * pow(eta, 2) / pow(1.0 - eta, 3);
        for (int i_vacancy = 0; i_vacancy < Nvacancy; i_vacancy++)//x_vacancy
        {
            double vacancy = vacancy_set[0] + vacancy_set[1] * i_vacancy;
            for (int i_gamma = 0; i_gamma < Ngamma; i_gamma++)//gamma
            {
                double gamma = gamma_set[0] + gamma_set[1] * i_gamma;
                err = set_parameter(dim, eta, gamma, vacancy, eps, prm);
                cout << eta << " " << prm.lambda << endl;
                double r_set[2] = { 0.0, prm.lambda / Nr };
                if (err)// lambda < 1.0 のときは飛ばす
                {
                    continue;
                }

                for (int rho_i = 0; rho_i < Nr; rho_i++)
                {
                    for (int i_rx = 0; i_rx < Nr; i_rx++)
                    {
                        r[0] = r_set[0] + r_set[1] * i_rx;
                        for (int i_ry = 0; i_ry < Nr; i_ry++)
                        {
                            r[1] = r_set[0] + r_set[1] * i_ry;
                            for (int i_rz = 0; i_rz < Nr; i_rz++)//gamma
                            {
                                r[2] = r_set[0] + r_set[1] * i_rz;
                                if (rho_i == 0)
                                {
                                    if (r[0] <= 1.0 && r[1] <= 1.0 && r[2] <= 1.0)
                                    {
                                        rho[i_rx][i_ry][i_rz] = 0.0;
                                    }
                                    else
                                    {
                                        rho[i_rx][i_ry][i_rz] = eta;
                                    }
                                }
                                /*else if (rho_i == 1)
                                {
                                    G[0][i_rx][i_ry][i_rz] = Fex_density_3D_rho_derivative(prm, r, betamu, eps);
                                    rho[i_rx][i_ry][i_rz] = G[0][i_rx][i_ry][i_rz];
                                }
                                else if (rho_i > 1)
                                {
                                    rho[i_rx][i_ry][i_rz] = alpha * rho[i_rx][i_ry][i_rz] + (1.0 - alpha) * G[0][i_rx][i_ry][i_rz];
                                    G[1][i_rx][i_ry][i_rz] = Fex_density_3D_rho_derivative(prm, r, betamu, eps);
                                    G[0][i_rx][i_ry][i_rz] = G[1][i_rx][i_ry][i_rz];
                                }*/
                                else if (rho_i == 1)
                                {
                                    G[0][i_rx][i_ry][i_rz] = testFex_density_3D_rho_derivative(prm, Nr, r, rho, betamu);
                                    rho[i_rx][i_ry][i_rz] = G[0][i_rx][i_ry][i_rz];
                                }
                                else if (rho_i > 1)
                                {
                                    rho[i_rx][i_ry][i_rz] = alpha * rho[i_rx][i_ry][i_rz] + (1.0 - alpha) * G[0][i_rx][i_ry][i_rz];
                                    G[1][i_rx][i_ry][i_rz] = testFex_density_3D_rho_derivative(prm, Nr, r, rho, betamu);
                                    G[0][i_rx][i_ry][i_rz] = G[1][i_rx][i_ry][i_rz];
                                }
                                
                                f_gr << fixed;
                                f_gr << setprecision(3) << r[0] << " " << setprecision(3) << r[1] << " " << setprecision(3) << r[2] << " " << setprecision(15) << rho[i_rx][i_ry][i_rz] << endl;
                            }
                        }
                    }
                    f_gr << endl;
                }
            }//i_gamma
            f_gr << endl;
        }//i_lambda
        f_gr << endl;
    }//i_eta
    f_gr.close();
    // 自由エネルギー最小の場合の密度プロファイル、荷重密度、Phi を描く
    /*for (int i_eta = 0; i_eta < Neta; i_eta++)
    {
        double eta = eta_set[0] + eta_set[1] * i_eta;
        double vacancy = vacancy_set[0] + vacancy_set[1] * Fmin_at[i_eta][0];
        double gamma = gamma_set[0] + gamma_set[1] * Fmin_at[i_eta][1];
        err = set_parameter(dim, eta, gamma, vacancy, eps, prm);
        bool write_title = false;
        if (i_eta == 0)
        {
            write_title = true;
            f_setting << endl << "When F is minimum:" << endl;
        }
        write_parameter(f_setting, prm, write_title, "");
        write_profile_3D(f_profile, prm, Ndiv, true, true, "#");
    }*/
    //f_setting.close();
    //f_profile.close();
    return 0;
}
