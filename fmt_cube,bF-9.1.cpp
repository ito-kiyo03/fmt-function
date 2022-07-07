#include <cmath>
#include <iostream>
#include <fstream>
//#include <string>
#include "fmt_cube2.h"
//#include "intde2.h"
using namespace::std;
int main()
{
    const int dim = 3, Neta = 1, Ngamma = 10, Nvacancy = 3, Ndiv = 21;//いくつ出力するか
    const double eps = 1.0e-8, eta_set[2] = { 0.97, 0.01 }, gamma_set[2] = { 10090.0, 1.0 }, vacancy_set[2] = { 0.000001, 0.001 };//{初期値 ,刻み}

    bool err;
    int Fmin_at[Neta][2];
    double F_density_set[3], Fmin[Neta];//0:id 1:ex 2:id+ex
    struct _parameter prm;

    //ofstream fig1("X-3.9eta,3D,0.3.txt");
    ofstream test1("n0n1n2n3_min,3Dz,0.9.dat",ios::app);
    ofstream test2("1-n3_min,3Dz,0.9.dat",ios::app);
    ofstream test3("rho_min,3Dz,0.9.dat", ios::app);
    int d = 2;//L:0,Sm:1,X:2

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
                {//スメクティック相
                    F_density_dDSm(dim, prm, eps, F_density_set);// 自由エネルギーを求める
                }
                else
                {//液相、固相
                    F_density_dD(dim, prm, eps, F_density_set);// 自由エネルギーを求める
                }

                /*bool write_title = false;
                if (i_eta == 0)
                {
                    write_title = true;
                }
                write_profile_3Dna(test1, prm, Ndiv, write_title, false, "#");
                write_profile_3DPhi(test2, prm, Ndiv, write_title, true, "#");
                write_profile_3Drho(test3, prm, Ndiv, write_title, false, "#");*/

                //write_profile_DSm(test1, prm, Ndiv, write_title, false, "#");
                //write_profile_2DSmn2(test2, prm, Ndiv, write_title, true, "#");
                //write_profile_2DrhoSm(test3, prm, Ndiv, write_title, false, "#");
                if ((i_vacancy == 0) && (i_gamma == 0))
                {
                    refresh_F_minimum(i_vacancy, i_gamma, F_density_set[2], Fmin_at[i_eta], Fmin[i_eta]);
                }
                else
                {
                    renew_F_minimum(i_vacancy, i_gamma, F_density_set[2], Fmin_at[i_eta], Fmin[i_eta]);
                }
            }//i_gamma
        }//i_lambda
    }//i_eta

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
        }
        write_profile_3Dna(test1, prm, Ndiv, write_title, false, "#");
        write_profile_3DPhi(test2, prm, Ndiv, write_title, true, "#");
        write_profile_3Drho(test3, prm, Ndiv, write_title, false, "#");
        //fig1 << prm.eta << " " << prm.gamma << " " << prm.vacancy << " " << prm.lambda << " " << Fmin[i_eta] - 3.9 * prm.eta << endl;
    }
    //fig1.close();
    test1.close();
    test2.close();
    test3.close();
    return 0;
}
