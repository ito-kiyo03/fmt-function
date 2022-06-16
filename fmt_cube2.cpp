#include <cmath>
#include <iostream>
#include <fstream>
//#include <string>
#include "fmt_cube2.h"
//#include "intde2.h"
using namespace::std;
int main()
{
    const int dim = 3, Neta = 1, Ngamma = 3, Nvacancy = 1, Ndiv = 21;//いくつ出力するか
    const double eps = 1.0e-8, eta_set[2] = { 0.99, 0.01 }, gamma_set[2] = { 942.0, 1.0 }, vacancy_set[2] = { 0.000001, 0.001 };//{初期値 ,刻み}

    bool err;
    int Fmin_at[Neta][2];
    double F_density_set[3], Fmin[Neta];//0:id 1:ex 2:id+ex
    struct _parameter prm;

    ofstream f_F("Xfmt_cube2_F,3D,0.99.dat");
    ofstream F_min("Xfmt_cube2_F_min,3D,0.99.dat");
    ofstream f_profile("Xfmt_cube2_profile,3D,0.99.dat");
    ofstream f_setting("Xfmt_cube2_setting,3D,0.99.dat");
    //ofstream rho3D("rho,3D,0.9.dat");
    int d = 2;//L:0,Sm:1,X:2

    write_setting(f_setting, dim, Neta, Nvacancy, Ngamma, eps, eta_set, vacancy_set, gamma_set);
    write_setting(F_min, dim, Neta, Nvacancy, Ngamma, eps, eta_set, vacancy_set, gamma_set);
    F_min << endl;
    for (int i_eta = 0; i_eta < Neta; i_eta++)//eta
    {
        double eta = eta_set[0] + eta_set[1] * i_eta;
        for (int i_vacancy = 0; i_vacancy < Nvacancy; i_vacancy++)//x_vacancy
        {
            double vacancy = vacancy_set[0] + vacancy_set[1] * i_vacancy;
            for (int i_gamma = 0; i_gamma < Ngamma; i_gamma++)//gamma
            {
                double gamma = gamma_set[0] + gamma_set[1] * i_gamma;
                err = set_parameter(dim, eta, gamma, vacancy, eps, prm);
                if (err)// lambda < 1.0 のときは飛ばす
                {
                    continue;
                }

                if (d == 1)
                {//スメスティック相
                    F_density_dDSm(dim, prm, eps, F_density_set);// 自由エネルギーを求める
                }
                else
                {//液相、固相
                    int F_id_approx = F_density_dD(dim, prm, eps, F_density_set);// 自由エネルギーを求める
                    cout << F_id_approx << endl;
                }

                f_F << prm.eta << " " << prm.gamma << " " << prm.lambda << " " << prm.vacancy << " " << F_density_set[0] << " " << F_density_set[1] << " " << F_density_set[2] << endl;

                if ((i_vacancy == 0) && (i_gamma == 0))
                {
                    refresh_F_minimum(i_vacancy, i_gamma, F_density_set[2], Fmin_at[i_eta], Fmin[i_eta]);
                }
                else
                {
                    renew_F_minimum(i_vacancy, i_gamma, F_density_set[2], Fmin_at[i_eta], Fmin[i_eta]);
                }
            }//i_gamma
            f_F << endl;
        }//i_lambda
        f_F << endl;
    }//i_eta
    f_F.close();
    // 自由エネルギー最小の場合の密度プロファイル、荷重密度、Phi を描く
    for (int i_eta = 0; i_eta < Neta; i_eta++)
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
        F_min << prm.eta << " " << prm.gamma << " " << prm.vacancy << " " << prm.lambda << " " << Fmin[i_eta] << endl;
        //write_profile_3Drho(rho3D, prm, Ndiv, write_title, false, "#");
    }
    F_min.close();
    f_setting.close();
    f_profile.close();
    //rho3D.close();
    return 0;
}
