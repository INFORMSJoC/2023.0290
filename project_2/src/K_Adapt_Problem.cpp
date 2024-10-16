#include "K_Adapt_Problem.h"

bool ddid_version; // if True, problem with Decision Dependent Information Discovery, otherwise, classical K-Adaptability problem.


bool K_Adapt_Problem::setInstance()
{
	bool check = true;
	
	check = this->define_sizes();
	check = this->define_maps();


	check = this->define_uncertainty_set();


	check = this->define_C_matrix();
	check = this->define_c_vector(); 

	check = this->define_D_matrix();
	check = this->define_d_vector();

	check = this->define_Q_matrix();
	check = this->define_q_vector();

	
	check = this->define_lbq_vector();
	check = this->define_lbc_vector();
	
	
	check = this->define_X_matrix();
	check = this->define_bx_vector();

	check = this->define_W_matrix();
	check = this->define_bw_vector();

	check = this->define_Y_matrix();
	check = this->define_by_vector();


	check = define_T_matrix();
	check = define_V_matrix();
	check = define_P_matrix();
	check = define_Txi_matrix();
	check = define_Pxi_matrix();
	check = define_Vxi_matrix();
	check = define_H_matrix();
	check = define_h_vector();

	check = define_t_vector();
	check = define_p_vector();
	check = define_v_vector();
	check = define_hmax_vector();

	check = define_Ex_matrix();
	check = define_Ew_matrix();
	check = define_Ey_matrix();
	check = define_g_vector();

	check = define_ub_xi();
	check = define_lb_xi();


	return check;
}

unsigned short K_Adapt_Problem::getNx()
{
	return this->Nx;
}

unsigned short K_Adapt_Problem::getNy()
{
	return this->Ny;
}

unsigned short K_Adapt_Problem::getNw()
{
	return this->Nw;
}

unsigned short K_Adapt_Problem::getNxi()
{
	return this->Nxi;
}

unsigned short K_Adapt_Problem::getNzeta()
{
	return this->Nzeta;
}

unsigned short K_Adapt_Problem::getK()
{
	return this->K;
}

unsigned short K_Adapt_Problem::getL()
{
	return this->L;
}

unsigned short K_Adapt_Problem::getLxi()
{

	return this->Lxi;
}

unsigned short K_Adapt_Problem::getLx()
{

	return this->Lx;
}

unsigned short K_Adapt_Problem::getLy()
{
	return this->Ly;
}

unsigned short K_Adapt_Problem::getLw()
{
	return this->Lw;
}

unsigned short K_Adapt_Problem::getR()
{
	return this->R;
}

Matrix2D<double> K_Adapt_Problem::get1S_uobj_matr_C()
{
	return this->C;
}

Matrix2D<double> K_Adapt_Problem::getDV_uobj_matr_D()
{
	return this->D;
}

Matrix2D<double> K_Adapt_Problem::get2S_uobj_matr_Q()
{
	
	return this->Q;
}

Matrix2D<double> K_Adapt_Problem::get1S_dcon_matr_X()
{
	return this->X;
}

std::vector<double> K_Adapt_Problem::get1S_drhs_vect_bx()
{
	return this->b_x;
}

Matrix2D<double> K_Adapt_Problem::getDV_dcon_matr_W()
{
	return this->W;
}

std::vector<double> K_Adapt_Problem::getDV_drhs_vect_bw()
{
	return this->b_w;
}

Matrix2D<double> K_Adapt_Problem::get2S_dcon_matr_Y()
{
	return this->Y;
}
std::vector<double> K_Adapt_Problem::get2S_drhs_vect_by()
{
	return this->b_y;
}
std::vector<double> K_Adapt_Problem::get1S_dobj_vect_c()
{
	return this->c;
}
std::vector<double> K_Adapt_Problem::get1S_dobj_vect_lbc()
{
	return this->lbc;
}
std::vector<double> K_Adapt_Problem::getDV_dobj_vect_d()
{
	return this->d;
}
std::vector<double> K_Adapt_Problem::get2S_dobj_vect_q()
{
	return this->q;
}
std::vector<double> K_Adapt_Problem::get2S_dobj_vect_lbq()
{
	return this->lbq;
}
Matrix2D<double> K_Adapt_Problem::getUC_xicon_matr_A()
{
	return this->A;
}
Matrix2D<double> K_Adapt_Problem::getUC_zetacon_matr_G()
{
	return this->G;
}
std::vector<double> K_Adapt_Problem::getUC_rhs_vect_b()
{
	return this->b;
}
Matrix2D<double> K_Adapt_Problem::getLink_dcon_matr_Ex()
{
	return this->Ex;
}
Matrix2D<double> K_Adapt_Problem::getLink_dcon_matr_Ew()
{
	return this->Ew;
}
Matrix2D<double> K_Adapt_Problem::getLink_dcon_matr_Ey()
{
	return this->Ey;
}
std::vector<double> K_Adapt_Problem::getLink_dcon_vect_g()
{
	return this->g;
}
Matrix2D<double> K_Adapt_Problem::get1S_ucon_matr_T()
{
	return this->T;
}
Matrix2D<double> K_Adapt_Problem::getDV_ucon_matr_V()
{
	return this->V;
}
Matrix2D<double> K_Adapt_Problem::get2S_ucon_matr_P()
{
	return this->P;
}
Matrix2D<double> K_Adapt_Problem::getUC_ucon_matr_H()
{
	return this->H;
}
std::vector<double> K_Adapt_Problem::get_rhs_ucon_vect_h()
{
	return this->h;
}
Matrix3D<double> K_Adapt_Problem::get1S_ucon_matr_Txi()
{
	return this->Txi;
}
Matrix3D<double> K_Adapt_Problem::getDV_ucon_matr_Vxi()
{
	return Vxi;
}
Matrix3D<double> K_Adapt_Problem::get2S_ucon_matr_Pxi()
{
	return this->Pxi;
}
std::vector<double> K_Adapt_Problem::get1S_ducon_vect_t()
{
	return this->t;
}
std::vector<double> K_Adapt_Problem::getDV_ducon_vect_v()
{
	return this->v;
}
std::vector<double> K_Adapt_Problem::get2S_ducon_vect_p()
{
	return this->p;
}
std::vector<double> K_Adapt_Problem::get_rhs_ducon_vect_hmax()
{
	return this->hmax;
}
std::vector<double> K_Adapt_Problem::get_ub_x()
{
	return this->ub_x;
}
std::vector<double> K_Adapt_Problem::get_ub_w()
{
	return this->ub_w;
}
std::vector<double> K_Adapt_Problem::get_ub_y()
{
	return this->ub_y;
}
std::vector<double> K_Adapt_Problem::get_ub_xi()
{
	return this->ub_xi;
}
std::vector<double> K_Adapt_Problem::get_ub_zeta()
{
	return this->ub_zeta;
}
std::vector<double> K_Adapt_Problem::get_lb_x()
{
	return this->lb_x;
}
std::vector<double> K_Adapt_Problem::get_lb_w()
{
	return this->lb_w;
}
std::vector<double> K_Adapt_Problem::get_lb_y()
{
	return this->lb_y;
}
std::vector<double> K_Adapt_Problem::get_lb_xi()
{
	return this->lb_xi;
}
std::vector<double> K_Adapt_Problem::get_lb_zeta()
{
	return this->lb_zeta;
}
std::string K_Adapt_Problem::getProbName()
{
	return this->problem_name;
}

bool K_Adapt_Problem::define_uncertainty_set()
{
	bool check = true;

	check = define_A_matrix();
	if (!check)
	{

		return false;
	}
	check = define_b_vector();
	if (!check)
	{
		return false;
	}
	return check;
}
