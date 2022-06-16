#ifndef _FMT_CUBE_H_
#define _FMT_CUBE_H_
#include <fstream>
#include <string>
using namespace::std;
struct _parameter {
	double eta, gamma, vacancy, lambda, coeff_rho, coeff_rhoSm, coeff_t, sqrt_gamma, cutoff_r_for_rho;
  int cutoff_m_for_rho;
};
bool set_parameter(const int dim, const double eta, const double gamma, const double vacancy, const double eps, _parameter& p);
double Fid_density_2D(const struct _parameter p, const int n);
double Fid_density_approx(const int dim, const struct _parameter p);
double Fex_density_2D(const struct _parameter prm, const double rel_error_req);
int F_density_dD(const int dim, const struct _parameter prm, const double rel_error_req, double F_density_set[]);
void refresh_F_minimum(const int i_vacancy, const int i_gamma, const double newF, int Fmin_at[2], double& Fmin);
void renew_F_minimum(const int i_vacancy, const int i_gamma, const double newF, int Fmin_at[2], double& Fmin);
void write_setting(ofstream &ofs, const int dim, const int Neta, const int Nvacancy, const int Ngamma, const double eps, const double eta_set[], const double vacancy_set[], const double gamma_set[]);
void write_parameter(ofstream &ofs, const _parameter prm, const bool write_title, const string pre_string = "");
void write_profile_2D(ofstream &ofs, const _parameter prm, const int Ndiv, const bool write_title, const bool quarter = true, const string pre_string = "");

double Fid_density_3D(const struct _parameter p, const int n);
double Fex_density_3D(const struct _parameter prm, const double rel_error_req);
void write_profile_3D(ofstream& ofs, const _parameter prm, const int Ndiv, const bool write_title, const bool quarter = true, const string pre_string = "");
void write_profile_3Dna(ofstream& ofs, const _parameter prm, const int Ndiv, const bool write_title, const bool quarter = true, const string pre_string = "");
void write_profile_3DPhi(ofstream& ofs, const _parameter prm, const int Ndiv, const bool write_title, const bool quarter = true, const string pre_string = "");
void write_profile_3Drho(ofstream& ofs, const _parameter prm, const int Ndiv, const bool write_title, const bool quarter = true, const string pre_string = "");

double Fid_density_2DSm(const struct _parameter p, const int n);
double Fid_density_approxSm(const int dim, const struct _parameter p);
double Fex_density_2DSm(const struct _parameter prm, const double rel_error_req);
void F_density_dDSm(const int dim, const struct _parameter prm, const double rel_error_req, double F_density_set[]);
void write_profile_2DSm(ofstream& ofs, const _parameter prm, const int Ndiv, const bool write_title, const bool quarter = true, const string pre_string = "");
void write_profile_2DSmn2(ofstream& ofs, const _parameter prm, const int Ndiv, const bool write_title, const bool quarter = true, const string pre_string = "");
void write_profile_2DrhoSm(ofstream& ofs, const _parameter prm, const int Ndiv, const bool write_title, const bool quarter = true, const string pre_string = "");
#endif
