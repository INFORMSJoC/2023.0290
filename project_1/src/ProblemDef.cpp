#include "ProblemDef.h"

bool ProblemDef::setInstance()
{


	/*here we set the configuration*/

	this->decision_dependent_discovery = this->cfg->getValueOfKey<bool>("DDID_VERSION");


	this->max_number_discovery_fraction = this->cfg->getValueOfKey<double>("MAX_NUM_DISCOVERY_FRACTION");

	stringstream ss;
	ss << cfg->getValueOfKey<string>("RESULTS_FOLDER");
	snprintf(results_folder, sizeof(results_folder), "%s", ss.str().c_str());


	bool check = true;

	check = this->define_sizes();
	if (!check)
	{
		std::cout << "Problem in bool ProblemDef::setInstance()" << std::endl;
		getchar();
		return check;
	}
	check = this->define_maps();
	if (!check)
	{
		std::cout << "Problem in bool ProblemDef::setInstance()" << std::endl;
		getchar();
		return check;
	}

	check = this->define_uncertainty_set();
	if (!check)
	{
		std::cout << "Problem in bool ProblemDef::setInstance()" << std::endl;
		getchar();
		return check;
	}

	check = this->define_C_matrix();
	if (!check)
	{
		std::cout << "Problem in bool ProblemDef::setInstance()" << std::endl;
		getchar();
		return check;
	}

	check = this->define_c_vector();
	if (!check)
	{
		std::cout << "Problem in bool ProblemDef::setInstance()" << std::endl;
		getchar();
		return check;
	}

	check = this->define_D_matrix();
	if (!check)
	{
		std::cout << "Problem in bool ProblemDef::setInstance()" << std::endl;
		getchar();
		return check;
	}


	check = this->define_d_vector();
	if (!check)
	{
		std::cout << "Problem in bool ProblemDef::setInstance()" << std::endl;
		getchar();
		return check;
	}


	check = this->define_Q_matrix();
	if (!check)
	{
		std::cout << "Problem in bool ProblemDef::setInstance()" << std::endl;
		getchar();
		return check;
	}


	check = this->define_q_vector();
	if (!check)
	{
		std::cout << "Problem in bool ProblemDef::setInstance()" << std::endl;
		getchar();
		return check;
	}

	check = this->define_lbq_vector();
	if (!check)
	{
		std::cout << "Problem in bool ProblemDef::setInstance()" << std::endl;
		getchar();
		return check;
	}
	check = this->define_lbc_vector();
	if (!check)
	{
		std::cout << "Problem in bool ProblemDef::setInstance()" << std::endl;
		getchar();
		return check;
	}

	check = this->define_X_matrix();
	if (!check)
	{
		std::cout << "Problem in bool ProblemDef::setInstance()" << std::endl;
		getchar();
		return check;
	}
	check = this->define_bx_vector();
	if (!check)
	{
		std::cout << "Problem in bool ProblemDef::setInstance()" << std::endl;
		getchar();
		return check;
	}

	check = this->define_W_matrix();
	if (!check)
	{
		std::cout << "Problem in bool ProblemDef::setInstance()" << std::endl;
		getchar();
		return check;
	}
	check = this->define_bw_vector();
	if (!check)
	{
		std::cout << "Problem in bool ProblemDef::setInstance()" << std::endl;
		getchar();
		return check;
	}

	check = this->define_Y_matrix();
	if (!check)
	{
		std::cout << "Problem in bool ProblemDef::setInstance()" << std::endl;
		getchar();
		return check;
	}
	check = this->define_by_vector();
	if (!check)
	{
		std::cout << "Problem in bool ProblemDef::setInstance()" << std::endl;
		getchar();
		return check;
	}
	check = define_T_matrix();
	if (!check)
	{
		std::cout << "Problem in bool ProblemDef::setInstance()" << std::endl;
		getchar();
		return check;
	}
	check = define_V_matrix();
	if (!check)
	{
		std::cout << "Problem in bool ProblemDef::setInstance()" << std::endl;
		getchar();
		return check;
	}
	check = define_P_matrix();
	if (!check)
	{
		std::cout << "Problem in bool ProblemDef::setInstance()" << std::endl;
		getchar();
		return check;
	}
	check = define_Txi_matrix();
	if (!check)
	{
		std::cout << "Problem in bool ProblemDef::setInstance()" << std::endl;
		getchar();
		return check;
	}
	check = define_Pxi_matrix();
	if (!check)
	{
		std::cout << "Problem in bool ProblemDef::setInstance()" << std::endl;
		getchar();
		return check;
	}
	check = define_Vxi_matrix();
	if (!check)
	{
		std::cout << "Problem in bool ProblemDef::setInstance()" << std::endl;
		getchar();
		return check;
	}
	check = define_H_matrix();
	if (!check)
	{
		std::cout << "Problem in bool ProblemDef::setInstance()" << std::endl;
		getchar();
		return check;
	}
	check = define_h_vector();
	if (!check)
	{
		std::cout << "Problem in bool ProblemDef::setInstance()" << std::endl;
		getchar();
		return check;
	}

	check = define_t_vector();
	if (!check)
	{
		std::cout << "Problem in bool ProblemDef::setInstance()" << std::endl;
		getchar();
		return check;
	}
	check = define_p_vector();
	if (!check)
	{
		std::cout << "Problem in bool ProblemDef::setInstance()" << std::endl;
		getchar();
		return check;
	}
	check = define_v_vector();
	if (!check)
	{
		std::cout << "Problem in bool ProblemDef::setInstance()" << std::endl;
		getchar();
		return check;
	}
	check = define_hmax_vector();
	if (!check)
	{
		std::cout << "Problem in bool ProblemDef::setInstance()" << std::endl;
		getchar();
		return check;
	}

	check = define_Ex_matrix();
	if (!check)
	{
		std::cout << "Problem in bool ProblemDef::setInstance()" << std::endl;
		getchar();
		return check;
	}
	check = define_Ew_matrix();
	if (!check)
	{
		std::cout << "Problem in bool ProblemDef::setInstance()" << std::endl;
		getchar();
		return check;
	}
	check = define_Ey_matrix();
	if (!check)
	{
		std::cout << "Problem in bool ProblemDef::setInstance()" << std::endl;
		getchar();
		return check;
	}
	check = define_g_vector();
	if (!check)
	{
		std::cout << "Problem in bool ProblemDef::setInstance()" << std::endl;
		getchar();
		return check;
	}
	check = define_ub_xi();
	if (!check)
	{
		std::cout << "Problem in bool ProblemDef::setInstance()" << std::endl;
		getchar();
		return check;
	}
	check = define_lb_xi();
	if (!check)
	{
		std::cout << "Problem in bool ProblemDef::setInstance()" << std::endl;
		getchar();
		return check;
	}

	check = define_valid_bounds();
	if (!check)
	{
		std::cout << "Problem in bool ProblemDef::setInstance()" << std::endl;
		getchar();
		return check;
	}


	return check;
}

void ProblemDef::setConfigFile(const char * file)
{
	 this->cfg = new Configuration::ConfigFile(file);/*use it as argument*/
}

void ProblemDef::freeMemory()
{
	delete this->cfg;
}

unsigned short ProblemDef::getNx()
{
	return this->Nx;
}

unsigned short ProblemDef::getNy()
{
	return this->Ny;
}

unsigned short ProblemDef::getNw()
{
	return this->Nw;
}

 unsigned short ProblemDef::getNxi()
{
	return this->Nxi;
}

inline unsigned short ProblemDef::getNzeta()
{
	return this->Nzeta;
}

unsigned short ProblemDef::getK()
{
	return this->K;
}

unsigned short ProblemDef::getL()
{
	return this->L;
}

unsigned short ProblemDef::getLxi()
{

	return this->Lxi;
}

unsigned short ProblemDef::getLx()
{

	return this->Lx;
}

unsigned short ProblemDef::getLy()
{
	return this->Ly;
}

unsigned short ProblemDef::getLw()
{
	return this->Lw;
}

unsigned short ProblemDef::getR()
{
	return this->R;
}

double ProblemDef::getMaxNumDiscoveryFraction()
{
	return this->max_number_discovery_fraction;
}

void ProblemDef::setMaxNumDiscoveryFraction(double mmdf)
{
	this->max_number_discovery_fraction = mmdf;
}

bool ProblemDef::getDDIDversion()
{
	return this->decision_dependent_discovery;
}

void ProblemDef::setDDIversion(bool decision_dependent_discovery)
{
	this->decision_dependent_discovery = decision_dependent_discovery;
}

Matrix2D<double> ProblemDef::get1S_uobj_matr_C()
{
	return this->C;
}

Matrix2D<double> ProblemDef::getDV_uobj_matr_D()
{
	return this->D;
}

Matrix2D<double> ProblemDef::get2S_uobj_matr_Q()
{

	return this->Q;
}

Matrix2D<double> ProblemDef::get1S_dcon_matr_X()
{
	return this->X;
}

std::vector<double> ProblemDef::get1S_drhs_vect_bx()
{
	return this->b_x;
}

Matrix2D<double> ProblemDef::getDV_dcon_matr_W()
{
	return this->W;
}

std::vector<double> ProblemDef::getDV_drhs_vect_bw()
{
	return this->b_w;
}

Matrix2D<double> ProblemDef::get2S_dcon_matr_Y()
{
	return this->Y;
}
std::vector<double> ProblemDef::get2S_drhs_vect_by()
{
	return this->b_y;
}
std::vector<double> ProblemDef::get1S_dobj_vect_c()
{
	return this->c;
}
std::vector<double> ProblemDef::get1S_dobj_vect_lbc()
{
	return this->lbc;
}
std::vector<double> ProblemDef::getDV_dobj_vect_d()
{
	return this->d;
}
std::vector<double> ProblemDef::get2S_dobj_vect_q()
{
	return this->q;
}
std::vector<double> ProblemDef::get2S_dobj_vect_lbq()
{
	return this->lbq;
}
Matrix2D<double> ProblemDef::getUC_xicon_matr_A()
{
	return this->A;
}
Matrix2D<double> ProblemDef::getUC_zetacon_matr_G()
{
	return this->G;
}
std::vector<double> ProblemDef::getUC_rhs_vect_b()
{
	return this->b;
}
Matrix2D<double> ProblemDef::getLink_dcon_matr_Ex()
{
	return this->Ex;
}
Matrix2D<double> ProblemDef::getLink_dcon_matr_Ew()
{
	return this->Ew;
}
Matrix2D<double> ProblemDef::getLink_dcon_matr_Ey()
{
	return this->Ey;
}
std::vector<double> ProblemDef::getLink_dcon_vect_g()
{
	return this->g;
}
Matrix2D<double> ProblemDef::get1S_ucon_matr_T()
{
	return this->T;
}
Matrix2D<double> ProblemDef::getDV_ucon_matr_V()
{
	return this->V;
}
Matrix2D<double> ProblemDef::get2S_ucon_matr_P()
{
	return this->P;
}
Matrix2D<double> ProblemDef::getUC_ucon_matr_H()
{
	return this->H;
}
std::vector<double> ProblemDef::get_rhs_ucon_vect_h()
{
	return this->h;
}
Matrix3D<double> ProblemDef::get1S_ucon_matr_Txi()
{
	return this->Txi;
}
Matrix3D<double> ProblemDef::getDV_ucon_matr_Vxi()
{
	return Vxi;
}
Matrix3D<double> ProblemDef::get2S_ucon_matr_Pxi()
{
	return this->Pxi;
}
std::vector<double> ProblemDef::get1S_ducon_vect_t()
{
	return this->t;
}
std::vector<double> ProblemDef::getDV_ducon_vect_v()
{
	return this->v;
}
std::vector<double> ProblemDef::get2S_ducon_vect_p()
{
	return this->p;
}
std::vector<double> ProblemDef::get_rhs_ducon_vect_hmax()
{
	return this->hmax;
}
std::vector<double> ProblemDef::get_ub_x()
{
	return this->ub_x;
}
std::vector<double> ProblemDef::get_ub_w()
{
	return this->ub_w;
}
std::vector<double> ProblemDef::get_ub_y()
{
	return this->ub_y;
}
std::vector<double> ProblemDef::get_ub_xi()
{
	return this->ub_xi;
}
std::vector<double> ProblemDef::get_ub_zeta()
{
	return this->ub_zeta;
}
std::vector<double> ProblemDef::get_lb_x()
{
	return this->lb_x;
}
std::vector<double> ProblemDef::get_lb_w()
{
	return this->lb_w;
}
std::vector<double> ProblemDef::get_lb_y()
{
	return this->lb_y;
}
std::vector<double> ProblemDef::get_lb_xi()
{
	return this->lb_xi;
}
std::vector<double> ProblemDef::get_lb_zeta()
{
	return this->lb_zeta;
}
std::string ProblemDef::getProbName()
{
	return this->problem_name;
}

double ProblemDef::getValidUB()
{
	return this->valid_ub;
}

double ProblemDef::getValidLB()
{
	return this->valid_lb;
}

bool ProblemDef::define_uncertainty_set()
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