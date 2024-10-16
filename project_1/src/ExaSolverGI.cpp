#include "ExaSolverGI.h"

bool ExaSolverGI::init_cplex(CPXENVptr* cpx_env, CPXLPptr* cpx_lp)
{
	int status = 0;

	// Open env.
	*cpx_env = CPXopenCPLEX(&status);
	if (this->checkCPXstatus(status, cpx_env, cpx_lp)) return false;
	*cpx_lp = CPXcreateprob(*cpx_env, &status, "ssp");
	if (this->checkCPXstatus(status, cpx_env, cpx_lp)) return false;

	return true;
}

int ExaSolverGI::checkCPXstatus(int status, CPXENVptr* cpx_env, CPXLPptr* cpx_lp)
{
	char errmsg[CPXMESSAGEBUFSIZE];
	if (status == 0) return 0;

	CPXgeterrorstring(*cpx_env, status, errmsg);
	printf(" %s \n", errmsg);
	return status;
}

bool ExaSolverGI::free_cplex(CPXENVptr* cpx_env, CPXLPptr* cpx_lp)
{
	int status = 0;
	status = CPXfreeprob(*cpx_env, cpx_lp);
	if (checkCPXstatus(status, cpx_env, cpx_lp)) return false;
	status = CPXcloseCPLEX(cpx_env);
	if (checkCPXstatus(status, cpx_env, cpx_lp)) return false;

	return true;
}


void ExaSolverGI::bm_reset_branc_and_cut_indicators()
{

	bm_num_iter_curr_node = 0;
	bm_num_cut_curr_iter = 0;
	bm_prev_node = 0;
	bm_curr_node = -1;
	bm_lp_prev_iter = 0;

}

bool ExaSolverGI::store_xbar()
{
	bool equal = true;
	if (num_ccg_xib >= ccg_set_xib.getNumRows())
	{
		//std::cout << "Max number of collected scenarios reached" << std::endl;
		//getchar();
		return false;
	}

	int c = 0;
	while (c < num_ccg_xib)
	 {
		equal = true;
	//	for (int i = this->first_pnode; i <= this->last_pnode; ++i)
		for(int i = 0; i < this->problem->getNxi(); ++i)
		{
			/*if different check the next c*/
			if (fabs(this->mlp_curr_sol[this->mlp_v_pr_xi[i]] - ccg_set_xib(c, i)) > EXA_EPS)
			{
				++c;
				equal = false;
				break;
			}
		}
		if (equal)
		{
			//std::cout << "Xi bar already existsed" << std::endl;
			return false;
		}
	} 
	
	/*We take the xi bar now*/
	//for (int i = this->first_pnode; i <= this->last_pnode; ++i)
	for (int i = 0; i < this->problem->getNxi(); ++i)
	{
	this->ccg_set_xib(num_ccg_xib, i) =	this->mlp_curr_sol[this->mlp_v_pr_xi[i]];
	}
	++num_ccg_xib;

	/*Also the xiY*/
	for (int y = 0; y < card_Y; ++y)
	{
		//for (int i = this->first_pnode; i <= this->last_pnode; ++i)
		for (int i = 0; i < this->problem->getNxi(); ++i)
		{
			this->ccg_set_xib(num_ccg_xib, i) = this->mlp_curr_sol[this->mlp_v_xi_y(y, i)];
		}
		++num_ccg_xib;
	}
	return true;
}

bool ExaSolverGI::bm_add_ccg_lb()
{



	if (!this->bm_add_ccg_columns())
	{
		std::cout << "Problems in bm_add_ccg_lb() after bm_add_ccg_columns()" << std::endl;
		return false;
	}

	if (!this->bm_add_ccg_constraints())
	{
		std::cout << "Problems in bm_add_ccg_lb() after bm_add_ccg_constraints()" << std::endl;
		return false;
	}

	return true;
}

bool ExaSolverGI::bm_add_ccg_columns()
{
	 int num_cols = CPXgetnumcols(this->cpx_env_bm, this->cpx_lp_bm);

	int status = 0;
	int vars_ind = 0;
	double lb[1], ub[1];
	double obj[1];
	char* varsname[1];
	char elname[1024];
	varsname[0] = elname;
	char vartype[1];
	double single_node_dur_1;
	double single_node_dur_2;
	int cnt;
	status = 0;

	/*
	 Block of variables for each \xi: y(xi), x(xi), sigma(xi), lambda(xi)
	*/
	for (int xi = 0; xi < this->num_ccg_xib; ++xi)
	{
		
		/*y(\xi, i)*/
		//for (int i = first_pnode; i <= last_pnode; ++i)
		for(int i = 0; i < this->problem->getNy(); ++i)
		{
			this->bm_v_y_xib_i(xi, i) = num_cols;
			obj[0] = 0;

				vartype[0] = 'B';

	sprintf(varsname[0], "y_xi%d_i%d", xi,i);
	lb[0] = 0;
	ub[0] = 1;




	status = CPXnewcols(this->cpx_env_bm, this->cpx_lp_bm, 1, obj, lb, ub, vartype, varsname);
	if (checkCPXstatus(status, &cpx_env_bm, &cpx_lp_bm)) goto TERMINATE;
	++num_cols;
		}



		///*lambda (\xi, r) r, num rows of matrix A*/
		//for (unsigned short r = 0; r < this->num_us_constr; ++r)
		for(int r = 0; r < this->problem->getLxi(); ++r)
		{
			this->bm_v_lamb_xib_i(xi, r) = num_cols;
			obj[0] = 0; 
			vartype[0] = 'C';

			sprintf(varsname[0], "lambda_xi%d_r%d", xi, r);

			lb[0] = 0;
			ub[0] = CPX_INFBOUND;

		status = CPXnewcols(this->cpx_env_bm, this->cpx_lp_bm, 1, obj, lb, ub, vartype, varsname);
		if (checkCPXstatus(status, &cpx_env_bm, &cpx_lp_bm)) goto TERMINATE;		++num_cols;
		}

		/*sigma (\xi, n) where n is the lenght of vector w*/
		//for (unsigned short n = this->first_pnode; n <= this->last_pnode; ++n)
		for(int n = 0; n < this->problem->getNw(); ++n)
		{

			this->bm_v_sig_xib_i(xi, n) = num_cols;
			obj[0] = 0;

			sprintf(varsname[0], "sigma_xi%d_r%d", xi, n);

			lb[0] = - CPX_INFBOUND;
			ub[0] =   CPX_INFBOUND;

			status = CPXnewcols(this->cpx_env_bm, this->cpx_lp_bm, 1, obj, lb, ub, vartype, varsname);
			if (checkCPXstatus(status, &cpx_env_bm, &cpx_lp_bm)) goto TERMINATE;		++num_cols;
		}

	}


TERMINATE:
	if (status)
		return false;
	else
		return true;
}

bool ExaSolverGI::bm_add_ccg_constraints()
{
	double obj[1], rhs[1];
	int status, matbeg[1], nzc;
	char sense[1], vartype[1];
	char* cnstrname[1];
	char elname[1024];
	char ss[1024];
	char locname[1024];
	cnstrname[0] = elname;
	double lb[1];
	double ub[1];
	matbeg[0] = 0;
	double eps_coeff = 0.0001;

	status = 0;

	if (matind.size() < CPXgetnumcols(cpx_env_bm, cpx_lp_bm))
	{
		matind.resize(CPXgetnumcols(cpx_env_bm, cpx_lp_bm));
		matval.resize(CPXgetnumcols(cpx_env_bm, cpx_lp_bm));

	}


	/*CONSTRAINTS OF THE SPP!!!! Y*/
		/*Primal constraints Y*/
	if (this->problem->getLy() > 0)
	{

		vector_loc = this->problem->get2S_drhs_vect_by();
		matr_loc = this->problem->get2S_dcon_matr_Y();
		
		for (int xi = 0; xi < this->num_ccg_xib; ++xi)
		{
			for (unsigned short row = 0; row < this->problem->getLy(); ++row)
			{
				rhs[0] = vector_loc[row];
				sense[0] = 'L';
				nzc = 0;
				sprintf(cnstrname[0], "Yxi_%dy_i%d", xi,row);
				// Now we go trought the "columns", coefficient of the matrices 

				for (unsigned short cols = 0; cols < this->problem->getNy(); ++cols)
					if (fabs(matr_loc(row, cols)) > eps_coeff)
					{
						this->matind[nzc] = this->bm_v_y_xib_i(xi,cols);
						//this->v_y_k(k, cols);
						this->matval[nzc] = matr_loc(row, cols);
						++nzc;
					}
				status = CPXaddrows(cpx_env_bm, cpx_lp_bm, 0, 1, nzc, rhs, sense, matbeg, matind.data(), matval.data(), NULL, cnstrname);
				if (checkCPXstatus(status, &cpx_env_bm, &cpx_lp_bm)) goto TERMINATE;
			}
		}
	}




	matr_loc = this->problem->get2S_uobj_matr_Q();
	// Degree constraints on profitable nodes
	for (int xi = 0; xi < this->num_ccg_xib; ++xi)
	{
		nzc = 0;
		rhs[0] = 0;
		sense[0] = 'E';


		/*Dual constraints*/
		nzc = 0;
		rhs[0] = 0;
		sense[0] = 'E';
	//	for (int i = first_pnode; i <= last_pnode; ++i)
		for(int i = 0; i < this->problem->getNxi(); ++i)
		{
			sprintf(cnstrname[0], "Balance_i%d_xi%d", i, xi);

			nzc = 0;
			rhs[0] = 0;
			sense[0] = 'E';

			matind[nzc] = this->bm_v_y_xib_i(xi, i);
			matval[nzc] = 0;
			//matval[nzc] = -ppxi[i];
			//matr_loc = this->problem->get2S_uobj_matr_Q();
			
			for (int j = 0; j < this->problem->getNy(); ++j)
			{
				matval[nzc] +=  (-matr_loc(i, j));
			}
			++nzc;

	
			/*gamma*/
			matind[nzc] = this->bm_v_sig_xib_i(xi, i);
			matval[nzc] = 1;
			++nzc;

			/*beta*/

			matr_loc = this->problem->getUC_xicon_matr_A();
			for(int l = 0; l < this->problem->getLxi(); ++l)
			if (fabs(matr_loc(l, i)) > eps_coeff)
				{
					matind[nzc] = this->bm_v_lamb_xib_i(xi, l);
					matval[nzc] = matr_loc(l, i);
					++nzc;
				}


			status = CPXaddrows(cpx_env_bm, cpx_lp_bm, 0, 1, nzc, rhs, sense, matbeg, matind.data(), matval.data(), NULL, cnstrname);
			if (checkCPXstatus(status, &cpx_env_bm, &cpx_lp_bm)) goto TERMINATE;
		}

		/*Epigraphic constraints*/
		sprintf(cnstrname[0], "Epigraphic_Constr");

		nzc = 0;
		rhs[0] = 0;
		sense[0] = 'G';

		matind[nzc] = this->bm_v_phi;
		matval[nzc] = 1;
		++nzc;

		this->vector_loc = this->problem->get2S_dobj_vect_q();

		//for (int i = first_pnode; i <= last_pnode; ++i)
		for(int i = 0; i < this->problem->getNy(); ++i)
		{
			matind[nzc] = this->bm_v_y_xib_i(xi, i);
			matval[nzc] = -vector_loc[i];
			++nzc;

			matind[nzc] = this->bm_v_sig_xib_i(xi, i);
			matval[nzc] = - ccg_set_xib(xi, i);
			++nzc;
		}

		vector_loc = this->problem->getUC_rhs_vect_b();

		//for (int l = 0; l < this->num_us_constr; ++l)
		for(int l = 0; l < this->problem->getLxi(); ++l)
			{
				matind[nzc] =   this->bm_v_lamb_xib_i(xi, l);
				matval[nzc] = - vector_loc[l];
				++nzc;
			}
		status = CPXaddrows(cpx_env_bm, cpx_lp_bm, 0, 1, nzc, rhs, sense, matbeg, matind.data(), matval.data(), NULL, cnstrname);
		if (checkCPXstatus(status, &cpx_env_bm, &cpx_lp_bm)) goto TERMINATE;


		/*Product w-sigma constraints*/
		//for (int i = first_pnode; i <= last_pnode; ++i)
		for(int i = 0; i < this->problem->getNw(); ++i)
		{
			nzc = 0;
			rhs[0] = 0;
			sense[0] = 'G';
			sprintf(cnstrname[0], "Linearize_sigmaW_GQ_%d", i);

			matind[nzc] = this->bm_v_sig_xib_i(xi, i);
			matval[nzc] = 1;
			++nzc;

			matind[nzc] = this->bm_v_w_i[i];
			matval[nzc] = BIG_M_SIGMA;
			++nzc;


			status = CPXaddrows(cpx_env_bm, cpx_lp_bm, 0, 1, nzc, rhs, sense, matbeg, matind.data(), matval.data(), NULL, cnstrname);
			if (checkCPXstatus(status, &cpx_env_bm, &cpx_lp_bm)) goto TERMINATE;

			nzc = 0;
			rhs[0] = 0;
			sense[0] = 'L';
			sprintf(cnstrname[0], "Linearize_sigmaW_LQ_%d", i);

			matind[nzc] = this->bm_v_sig_xib_i(xi, i);
			matval[nzc] = 1;
			++nzc;

			matind[nzc] = this->bm_v_w_i[i];
			matval[nzc] = -BIG_M_SIGMA;
			++nzc;
			status = CPXaddrows(cpx_env_bm, cpx_lp_bm, 0, 1, nzc, rhs, sense, matbeg, matind.data(), matval.data(), NULL, cnstrname);
			if (checkCPXstatus(status, &cpx_env_bm, &cpx_lp_bm)) goto TERMINATE;
		}

	}


TERMINATE:
	if (status)
		return false;
	else
		return true;
}

inline bool ExaSolverGI::store_set_of_tight_y()
{
	double slack[1];
	int tau_y_ind;
	int status = 0;

	/*Becase it's a new solution*/
	this->num_tight_y = 0;

	for (int y = 0; y < this->mlp_num_y; ++y)
	{
		tau_y_ind = this->mlp_cnstr_tau_y[y];

		status = CPXgetslack(this->cpx_env_mlp, cpx_lp_mlp, slack, tau_y_ind, tau_y_ind);
		if (checkCPXstatus(status, &cpx_env_mlp, &cpx_lp_mlp)) goto TERMINATE;

		if (std::fabs(slack[0]) < 1e-8)
		{
			this->store_tight_y(y);
		}

	}


TERMINATE:
	if (status != 0)
	{
		std::cout << "Error while storing the tight y in bool ExaSolverGI::store_tight_y()" << std::endl;
		getchar();
		return false;
	}
	else
		return true;
	
}

inline void ExaSolverGI::store_tight_y(int y)
{
	if (num_tight_y >= tight_y.size())
	{
		tight_y.resize(2 * num_tight_y);
		tight_xibar.resize(2 * num_tight_y);
		tight_xi_y.resizeMatrix2D(2 * num_tight_y, this->problem->getNy());

	}


	this->tight_y[this->num_tight_y] = y;


	//for (int i = this->first_pnode; i <= this->last_pnode; ++i)
	for(int i = 0; i < this->problem->getNy(); ++i)
	{
		this->tight_xi_y(num_tight_y, i) = this->mlp_curr_sol[this->mlp_v_xi_y(y,i)];
		this->tight_xibar[i] = this->mlp_curr_sol[this->mlp_v_pr_xi[i]];
	}

	++num_tight_y;
}


bool ExaSolverGI::setConfigFile(const char* conf_file)
{
	//cfg = new Configuration::ConfigFile("config.cfg");

	cfg = new Configuration::ConfigFile(conf_file);

	this->use_xi_bar_scenarios = cfg->getValueOfKey<bool>("EXACT_GI_USE_SCENARIOS");

	
	return true;
}

bool ExaSolverGI::init_data_structures()
{

	//snprintf(results_folder, sizeof(results_folder), "%s", RESULTS_FOLDER);

	//cfg = new Configuration::ConfigFile("config.cfg");

	//stringstream ss;
	//string out_dir;
	//ss << cfg->getValueOfKey<string>("RESULTS_FOLDER");
	//out_dir = ss.str();

	//if (!std::filesystem::create_directory(ss.str()))
	//	cout << "Warning: " << strerror(errno) << endl;

	//snprintf(results_folder, sizeof(results_folder), "%s", ss.str().c_str());




	int num_vrs_op = 0;
	int num_vrs_lp = 0;

	this->rop_v_y_n.resize(this->problem->getNy());
	num_vrs_op += rop_v_y_n.size();


	this->v_beta.resize(this->problem->getLxi());
	num_vrs_op += v_beta.size();

	this->v_gamma.resize(this->problem->getNw());
	num_vrs_op += v_gamma.size();





	/*Now we set for the LP*/
	this->mlp_num_y = 0;
	card_Y = 0;
	Y.resizeMatrix2D(this->problem->getNy(), this->problem->getNy()); // Each row is a solution of the OP.
	y_in_mlp.assign(this->problem->getNy(), false);
	this->mlp_v_pr_xi.resize(this->problem->getNxi());

	num_vrs_lp += mlp_v_pr_xi.size();
	this->mlp_v_xi_y.resizeMatrix2D(this->problem->getNxi(), this->problem->getNy()); /*We start with a size of one route for each node.*/
	mlp_cnstr_tau_y.resize(this->problem->getNy());

	num_vrs_lp += mlp_v_xi_y.getNumElem();
	
	/*Master Problem: assume number of scenarios euqal to number of nodes*/
	bm_v_w_i.resize(this->problem->getNw());

	bm_v_y_xib_i.resizeMatrix2D(MAX_NUM_SCENARIOS, this->problem->getNy());
	bm_v_lamb_xib_i.resizeMatrix2D(MAX_NUM_SCENARIOS,this->problem->getLxi());
	bm_v_sig_xib_i.resizeMatrix2D(MAX_NUM_SCENARIOS, this->problem->getNw());
	
	


	this->num_tight_y = 0;
	tight_y.resize(this->problem->getNy());
	this->tight_xi_y.resizeMatrix2D(this->problem->getNxi(), this->problem->getNy());
	this->tight_xibar.resize(this->problem->getNxi());

	/*matind and matval*/
	matind.resize(std::max(num_vrs_op, num_vrs_lp));
	matval.resize(std::max(num_vrs_op, num_vrs_lp));

	rop_curr_sol.resize(std::max(num_vrs_op, num_vrs_lp));
	mlp_curr_sol.resize(std::max(num_vrs_op, num_vrs_lp));

	bm_curr_sol.resize(10* this->problem->getNw() + 1);

	ccg_set_xib.resizeMatrix2D(MAX_NUM_SCENARIOS, this->problem->getNy());





	/*
	cut counters to 0
	*/
	this->bm_reset_branc_and_cut_indicators();


	num_ccg_xib = 0;

this->num_logic_benders_cuts =0;

bm_num_stored_logic_cuts = 0;
bm_logic_cuts_pool_set.resize(static_cast<std::vector<uint64_t, std::allocator<uint64_t>>::size_type>(10) * this->problem->getNw());
bm_logic_cuts_pool_phi_val.resize( static_cast<std::vector<double, std::allocator<double>>::size_type>(10)		* this->problem->getNw());





	return true;
}


bool ExaSolverGI::clean_data_structure()
{
	delete cfg;

	/*Clean cplex op*/
	int status = 0;
	if (this->cpx_env_rop != NULL &&
		this->cpx_lp_rop != NULL)
	{
		status = CPXfreeprob(this->cpx_env_rop, &this->cpx_lp_rop);
		if (checkCPXstatus(status, &this->cpx_env_rop, &this->cpx_lp_rop)) return false;
		status = CPXcloseCPLEX(&this->cpx_env_rop);
		if (checkCPXstatus(status, &this->cpx_env_rop, &this->cpx_lp_rop)) return false;
	}

	/*Clean mlp*/
	if (this->cpx_env_mlp != NULL &&
		this->cpx_lp_mlp != NULL)
	{
		status = CPXfreeprob(this->cpx_env_mlp, &this->cpx_lp_mlp);
		if (checkCPXstatus(status, &this->cpx_env_mlp, &this->cpx_lp_mlp)) return false;
		status = CPXcloseCPLEX(&this->cpx_env_mlp);
		if (checkCPXstatus(status, &this->cpx_env_mlp, &this->cpx_lp_mlp)) return false;
	}

	/*Clean bm*/
	if (this->cpx_env_bm != NULL &&
		this->cpx_lp_bm != NULL)
	{
		status = CPXfreeprob(this->cpx_env_bm, &this->cpx_lp_bm);
		if (checkCPXstatus(status, &this->cpx_env_bm, &this->cpx_lp_bm)) return false;
		status = CPXcloseCPLEX(&this->cpx_env_bm);
		if (checkCPXstatus(status, &this->cpx_env_bm, &this->cpx_lp_bm)) return false;
	}

	return true;
}

void ExaSolverGI::set_problem(ProblemDef* problem)
{
	this->problem = problem;
}


bool ExaSolverGI::load_and_init(std::vector<bool> w_val)
{
	bool is_ok = true;

	this->w_val = w_val;


	this->bm_cut_lb = this->problem->getValidLB();
	this->bm_best_lb = bm_cut_lb;

	/*without policy ccg is + inf*/
	this->mlp_obj_val =   DBL_MAX;

	this->rop_obj_val = - DBL_MAX;

	rop_is_feasible = true;

	is_ok = init_data_structures();
	if (!is_ok)
	{
		goto TERMINATE;
	}

	is_ok = mlp_build_model();
	if (!is_ok)
	{
		goto TERMINATE;
	}

	is_ok = mlp_update_model();
	if (!is_ok)
	{
		goto TERMINATE;
	}


	is_ok = mlp_solve_model();
	if (!is_ok)
	{
		goto TERMINATE;
	}


	is_ok = rop_build_model();
	if (!is_ok)
	{
		goto TERMINATE;
	}




TERMINATE:
	return is_ok;
}

bool ExaSolverGI::run_sub_problem_solver(double* obj_val)
{
	bool is_ok = true;

	is_ok = mlp_solve_model();
	if (!is_ok)
	{
		goto TERMINATE;
	}


	do {

		is_ok = rop_update_scenario_parameter();
		if (!is_ok)
		{
			goto TERMINATE;
		}

		is_ok = rop_solve_model();
		if (!is_ok)
		{


			goto TERMINATE;
		}

		if (this->mlp_obj_val > this->rop_obj_val + EXA_EPS)
		{
			is_ok = mlp_update_model();
			if (!is_ok)
			{
				goto TERMINATE;
			}
			is_ok = mlp_solve_model();
			if (!is_ok)
			{
				goto TERMINATE;
			}
		}
		else
		{
			if (fabs(this->mlp_obj_val - this->rop_obj_val) > (EXA_EPS))
			{
			
				std::cout << "Check eps in the Column-And-Constraints-Generation" << std::endl;
				getchar();
				getchar();
				getchar();
			}
			break;

		}

	


	} while (true);
	//while (fabs(this->mlp_obj_val - this->rop_obj_val) > EXA_EPS);


	*obj_val = mlp_obj_val;

TERMINATE:
	return is_ok;
}

bool ExaSolverGI::rop_build_model()
{
	bool is_ok = true;
	int status;
	int mip_status;
	CPXLONG      contextmask = 0;
	char errbuf[CPXMESSAGEBUFSIZE];
	int check_mip_stat;

	double start = 0;

	double end = 0;

	is_ok = init_cplex(&this->cpx_env_rop, &this->cpx_lp_rop);
	if (!is_ok)
	{
		goto TERMINATE;
	}

	is_ok = rop_set_cplex();
	if (!is_ok)
	{
		goto TERMINATE;
	}

	is_ok = rop_add_vars();
	if (!is_ok)
	{
		goto TERMINATE;
	}

	is_ok = rop_add_cnstr();
	if (!is_ok)
	{
		goto TERMINATE;
	}





	

TERMINATE:
	return is_ok;
}

bool ExaSolverGI::rop_set_cplex()
{
	bool mod_stat = true;

	int status = 0;

	char errbuf[CPXMESSAGEBUFSIZE];


	status = CPXsetintparam(this->cpx_env_rop, CPX_PARAM_SCRIND, CPX_OFF);
	if (checkCPXstatus(status, &cpx_env_rop, &cpx_lp_rop)) return false;
	status = CPXsetintparam(this->cpx_env_rop, CPX_PARAM_THREADS, 1);
	if (checkCPXstatus(status, &cpx_env_rop, &cpx_lp_rop)) return false;
	status = CPXchgobjsen(this->cpx_env_rop, this->cpx_lp_rop, CPX_MIN);
	if (checkCPXstatus(status, &cpx_env_rop, &cpx_lp_rop)) return false;
	status = CPXsetintparam(this->cpx_env_rop, CPX_PARAM_NUMERICALEMPHASIS, CPX_ON);
	if (checkCPXstatus(status, &cpx_env_rop, &cpx_lp_rop)) return false;
	status = CPXsetdblparam(this->cpx_env_rop, CPXPARAM_Simplex_Tolerances_Feasibility, 1e-9);
	if (checkCPXstatus(status, &cpx_env_rop, &cpx_lp_rop)) return false;
	//status = CPXsetintparam(cpx_env, CPX_PARAM_PREIND, CPX_OFF);
	//if (checkCPXstatus(status))return false;

	/*To be extra numerical accurate*/
	status = CPXsetdblparam(this->cpx_env_rop, CPXPARAM_Simplex_Tolerances_Optimality, 1e-9);
	if (checkCPXstatus(status, &cpx_env_rop, &cpx_lp_rop)) return false;

	status = CPXsetdblparam(this->cpx_env_rop, CPXPARAM_Simplex_Tolerances_Feasibility, 1e-9);
	if (checkCPXstatus(status, &cpx_env_rop, &cpx_lp_rop)) return false;

	status = CPXsetdblparam(this->cpx_env_rop, CPXPARAM_MIP_Tolerances_MIPGap, 1e-9);
	if (checkCPXstatus(status, &cpx_env_rop, &cpx_lp_rop)) return false;
	status = CPXsetdblparam(this->cpx_env_rop, CPXPARAM_MIP_Tolerances_Integrality, 0);
	if (checkCPXstatus(status, &cpx_env_rop, &cpx_lp_rop)) return false;

	double time_limit = this->cfg->getValueOfKey<double>("TIME_LIMIT");

	status = CPXsetdblparam(cpx_env_rop, CPXPARAM_TimeLimit, time_limit);
	if (checkCPXstatus(status, &cpx_env_rop, &cpx_lp_rop)) return false;

	/*
	status = CPXsetintparam(cpx_env, CPXPARAM_MIP_Strategy_VariableSelect, CPX_VARSEL_STRONG);
	if (checkCPXstatus(status)) return false;
	*/


	return true;
}

bool ExaSolverGI::rop_add_vars()
{
	int num_cols = CPXgetnumcols(this->cpx_env_rop, this->cpx_lp_rop);

	int status = 0;
	int vars_ind = 0;
	double lb[1], ub[1];
	double obj[1];
	char* varsname[1];
	char elname[1024];
	varsname[0] = elname;
	char vartype[1];
	double single_node_dur_1;
	double single_node_dur_2;
	int cnt;


	
	double rounded_to_precision = 0;


	/*
	 Routing vatriables, y_kn and x_kij.
	*/
	/*y_kn*/

	if (this->problem->getNy() > 0)
	{
		this->vector_loc = problem->get2S_dobj_vect_q();

		//this->rop_v_y_n.resize(this->problem->getNy());

			for (int i = 0; i < this->problem->getNy(); ++i)
			{


				this->rop_v_y_n[i] = num_cols;

				obj[0] = vector_loc[i];

				//obj[0] = 0;

				vartype[0] = 'B';

				sprintf(varsname[0], "y_i%d", i);
				
				
				lb[0] = 0;
				ub[0] = 1;




				status = CPXnewcols(this->cpx_env_rop, this->cpx_lp_rop, 1, obj, lb, ub, vartype, varsname);
				if (checkCPXstatus(status, &cpx_env_rop, &cpx_lp_rop)) goto TERMINATE;
				++num_cols;
			}
		
	}


	/*Gamma variables*/
	if (this->problem->getNw() > 0)
	{
		//for (unsigned short n = this->first_pnode; n <= this->last_pnode; ++n)
		for(int n = 0; n < problem->getNw();++n)
		{
			v_gamma[n] = num_cols;

			if (
				fabs(this->mlp_curr_sol[mlp_v_pr_xi[n]]) >= 1e-9
				)
			{

				obj[0] = this->mlp_curr_sol[mlp_v_pr_xi[n]];
			}
			else
				obj[0] = 0;

			vartype[0] = 'C';
			sprintf(varsname[0], "gamma_n%d", n);


			if (this->w_val[n])
			{
				lb[0] = -CPX_INFBOUND;
				ub[0] = CPX_INFBOUND;
			}
			else
			{
				lb[0] = 0;
				ub[0] = 0;
			}

			status = CPXnewcols(this->cpx_env_rop, this->cpx_lp_rop, 1, obj, lb, ub, vartype, varsname);
			if (checkCPXstatus(status, &cpx_env_rop, &cpx_lp_rop)) goto TERMINATE;		++num_cols;
		}

	}

	

	/*Beta variables*/
	if (this->problem->getLxi() > 0)
	{
		this->vector_loc = this->problem->getUC_rhs_vect_b();

		for (int r = 0; r < this->problem->getLxi(); ++r)
		{
			v_beta[r] = num_cols;
			obj[0] = this->vector_loc[r];

			vartype[0] = 'C';

			sprintf(varsname[0], "beta_r%d", r);

			lb[0] = 0;
			ub[0] = CPX_INFBOUND;

			status = CPXnewcols(this->cpx_env_rop, this->cpx_lp_rop, 1, obj, lb, ub, vartype, varsname);
			if (checkCPXstatus(status, &cpx_env_rop, &cpx_lp_rop)) goto TERMINATE;		++num_cols;

		}
	}



TERMINATE:

	if (status)
	{
		return false;
	}

	return true;
}

bool ExaSolverGI::rop_add_cnstr()
{
	double obj[1], rhs[1];
	int status, matbeg[1], nzc;
	char sense[1], vartype[1];
	char* cnstrname[1];
	char elname[1024];
	char ss[1024];
	char locname[1024];
	cnstrname[0] = elname;
	double lb[1];
	double ub[1];
	matbeg[0] = 0;
	double coeff_eps = 0.0001;


	 

	/*Primal constraints Y*/
	if (this->problem->getLy() > 0)
	{

		vector_loc = this->problem->get2S_drhs_vect_by();
		matr_loc = this->problem->get2S_dcon_matr_Y();

		for (unsigned short row = 0; row < this->problem->getLy(); ++row)
		{
			rhs[0] = vector_loc[row];
			sense[0] = 'L';
			nzc = 0;
			sprintf(cnstrname[0], "Yy_i%d", row);
			// Now we go trought the "columns", coefficient of the matrices 

			for (unsigned short cols = 0; cols < this->problem->getNy(); ++cols)
				if (fabs(matr_loc(row, cols)) > coeff_eps)
				{
					this->matind[nzc] = this->rop_v_y_n[cols];
						//this->v_y_k(k, cols);
					this->matval[nzc] = matr_loc(row, cols);
					++nzc;
				}
			status = CPXaddrows(cpx_env_rop, cpx_lp_rop, 0, 1, nzc, rhs, sense, matbeg, matind.data(), matval.data(), NULL, cnstrname);
			if (checkCPXstatus(status, &cpx_env_rop, &cpx_lp_rop)) goto TERMINATE;
		}
	}

	/*Dual constraints*/
	/*HEEEREEEEEE we need to write the constraints*/
	if (this->problem->getNxi() > 0)
	{
	


		for (int i = 0; i < this->problem->getNxi(); ++i)
		{
			rhs[0] = 0;
			sense[0] = 'E';
			sprintf(cnstrname[0], "BalanceK_r%d", i);
			nzc = 0;


			/* \gamma \times w*/
			if (this->w_val[i])
			{
				this->matind[nzc] = this->v_gamma[i];
				this->matval[nzc] = 1;
				++nzc;
			}
			else {
				this->matind[nzc] = this->v_gamma[i];
				this->matval[nzc] = 0;
				++nzc;
			}



			/*A^(T)Beta */
			vector_loc = this->problem->getUC_rhs_vect_b();
			matr_loc = this->problem->getUC_xicon_matr_A();


			for (int r = 0; r < this->problem->getLxi(); ++r)
			{
				this->matind[nzc] = this->v_beta[r];
				this->matval[nzc] = matr_loc(r, i);
				++nzc;

			}

			/*now we get Q*/
			matr_loc = this->problem->get2S_uobj_matr_Q();
			for (int col = 0; col < this->problem->getNy(); ++col)
			{
				this->matind[nzc] =   this->rop_v_y_n[col];
				this->matval[nzc] = - matr_loc(i, col);
				++nzc;

			}


			status = CPXaddrows(cpx_env_rop, cpx_lp_rop, 0, 1, nzc, rhs, sense, matbeg, matind.data(), matval.data(), NULL, cnstrname);
			if (checkCPXstatus(status, &cpx_env_rop, &cpx_lp_rop)) goto TERMINATE;
		}


	}






		/*Dual constraints*/
		//nzc = 0;
		//rhs[0] = 0;
		//sense[0] = 'E';

		//for (int i = first_pnode; i <= last_pnode; ++i)
		//{
		//	sprintf(cnstrname[0], "Balance_i%d", i);

		//	nzc = 0;
		//	rhs[0] = 0;
		//	sense[0] = 'E';

		//	matind[nzc] = rop_v_y_n[i];
		//	matval[nzc] = -ppxi[i];
		//	++nzc;


		//	matind[nzc] = v_gamma[i];
		//	matval[nzc] = 1; // Becase we cnage the bound of gamma when w = 0.
		//	//if (w_val[i])
		//	//{
		//	//	//matval[nzc] = -1;
		//	//	matval[nzc] = 1;
		//	//}
		//	//else
		//	//{
		//	//	matval[nzc] = 0;
		//	//}
		//	++nzc;

		//	for (int l = 0; l < this->num_us_constr; ++l)
		//		if (fabs(A(l, i)) > eps_coeff)
		//		{
		//			matind[nzc] = v_beta[l];
		//			matval[nzc] = A(l, i);
		//			++nzc;
		//		}

		//	status = CPXaddrows(cpx_env_op, cpx_lp_op, 0, 1, nzc, rhs, sense, matbeg, matind.data(), matval.data(), NULL, cnstrname);
		//	if (checkCPXstatus(status, &cpx_env_op, &cpx_lp_op)) goto TERMINATE;
		//}



TERMINATE:
	if (status)
		return false;
	else
		return true;
}

bool ExaSolverGI::rop_update_scenario_parameter()
{
	bool mod_stat = true;

	int status;
	int cnt = 0;
	double rounded_to_precision = 0;


	for (int n = 0; n < this->problem->getNxi(); ++n)
	{
		this->matind[cnt] = v_gamma[n];
		if (w_val[n])
		{

			//rounded_to_precision = std::round(this->mlp_curr_sol[this->v_pr_xi[n]] * 1e9) / 1e9;


			this->matval[cnt] = this->mlp_curr_sol[this->mlp_v_pr_xi[n]];

			//	this->matval[cnt] = rounded_to_precision;
		}
		else
		{
			this->matval[cnt] = 0;
		}
		++cnt;
	}


	//for (unsigned short n = this->first_pnode; n <= this->last_pnode; ++n)
	//{

	//	this->matind[cnt] = v_gamma[n];
	//	if (w_val[n])
	//	{

	//		//rounded_to_precision = std::round(this->mlp_curr_sol[this->v_pr_xi[n]] * 1e9) / 1e9;


	//		this->matval[cnt] = this->mlp_curr_sol[this->mlp_v_pr_xi[n]];

	//	//	this->matval[cnt] = rounded_to_precision;
	//	}
	//	else
	//	{
	//		this->matval[cnt] = 0;
	//	}
	//	++cnt;
	//}

	status = CPXchgobj(cpx_env_rop, cpx_lp_rop, cnt, matind.data(), matval.data());
	if (checkCPXstatus(status, &cpx_env_rop, &cpx_lp_rop)) goto TERMINATE;




TERMINATE:
	if (status)
		return false;
	else
		return true;
}

bool ExaSolverGI::rop_solve_model()
{
	bool is_ok = true;
	int status;
	int mip_status;
	CPXLONG      contextmask = 0;
	char errbuf[CPXMESSAGEBUFSIZE];
	int check_mip_stat;

	//status = CPXwriteprob(this->cpx_env_rop, this->cpx_lp_rop, "model_rop.lp", "lp");
	//if (checkCPXstatus(status, &cpx_env_rop, &cpx_lp_rop)) { is_ok = false; goto TERMINATE; };


	status = CPXmipopt(this->cpx_env_rop, this->cpx_lp_rop);
	if (checkCPXstatus(status, &cpx_env_rop, &cpx_lp_rop)) goto TERMINATE;



	check_mip_stat = CPXgetstat(cpx_env_rop, cpx_lp_rop);


	if (check_mip_stat == CPXMIP_INFEASIBLE)
	{
		this->rop_is_feasible = false;

		this->bm_best_lb = -1;
		this->bm_best_ub = -1;




		std::cout << "Nominal problem (Y) is infeasible" << std::endl;
		is_ok = false;
		goto TERMINATE;

	}


	status = CPXgetobjval(cpx_env_rop, cpx_lp_rop, &this->rop_obj_val);
	if (checkCPXstatus(status, &cpx_env_rop, &cpx_lp_rop)) goto TERMINATE;


	if (this->mlp_obj_val > this->rop_obj_val + EXA_EPS)
	{
		status = CPXgetx(cpx_env_rop, cpx_lp_rop, rop_curr_sol.data(), 0, CPXgetnumcols(cpx_env_rop, cpx_lp_rop) - 1);
		if (checkCPXstatus(status, &cpx_env_rop, &cpx_lp_rop)) goto TERMINATE;


		if (Y.getNumRows() <= card_Y)
		{
			Y.resizeMatrix2D(card_Y * 2, 
				this->problem->getNy()
			);
		}

		//for (int i = first_pnode; i <= last_pnode; ++i)
		for(int i = 0; i < this->problem->getNy(); ++i)
		{
			Y(card_Y, i) = round(rop_curr_sol[this->rop_v_y_n[i]]);
		}
		++card_Y;
	}



	/*if (check_mip_stat != CPXMIP_TIME_LIM_INFEAS)
	{
		status = CPXgetobjval(cpx_env, cpx_lp, &this->best_ub);
		if (checkCPXstatus(status)) goto TERMINATE;
	}
	else
	{
		this->best_ub = -1;

	}*/

TERMINATE:
	return is_ok;
}

bool ExaSolverGI::mlp_build_model()
{
	bool is_ok = true;

	
	is_ok = init_cplex(&this->cpx_env_mlp, &this->cpx_lp_mlp);
	if (!is_ok)
	{
		goto TERMINATE;
	}

	is_ok = mlp_set_cplex();
	if (!is_ok)
	{
		goto TERMINATE;
	}



	is_ok = mlp_add_vars();
	if (!is_ok)
	{
		goto TERMINATE;
	}

	is_ok = mlp_add_cnstr();
	if (!is_ok)
	{
		goto TERMINATE;
	}



TERMINATE:
	return is_ok;
}

bool ExaSolverGI::mlp_set_cplex()
{
	bool mod_stat = true;

	int status = 0;

	char errbuf[CPXMESSAGEBUFSIZE];


	status = CPXsetintparam(this->cpx_env_mlp, CPX_PARAM_SCRIND, CPX_OFF);
	if (checkCPXstatus(status, &cpx_env_mlp, &cpx_lp_mlp)) return false;
	status = CPXsetintparam(this->cpx_env_mlp, CPX_PARAM_THREADS, 1);
	if (checkCPXstatus(status, &cpx_env_mlp, &cpx_lp_mlp)) return false;
	status = CPXchgobjsen(this->cpx_env_mlp, this->cpx_lp_mlp, CPX_MAX);
	if (checkCPXstatus(status, &cpx_env_mlp, &cpx_lp_mlp)) return false;
	status = CPXsetintparam(this->cpx_env_mlp, CPX_PARAM_NUMERICALEMPHASIS, CPX_ON);
	if (checkCPXstatus(status, &cpx_env_mlp, &cpx_lp_mlp)) return false;
	status = CPXsetdblparam(this->cpx_env_mlp, CPXPARAM_Simplex_Tolerances_Feasibility, 1e-9);
	if (checkCPXstatus(status, &cpx_env_mlp, &cpx_lp_mlp)) return false;


	/*To be extra numerical accurate*/
	status = CPXsetdblparam(this->cpx_env_mlp, CPXPARAM_Simplex_Tolerances_Optimality, 1e-9);
	if (checkCPXstatus(status, &cpx_env_mlp, &cpx_lp_mlp)) return false;




	return true;
}

bool ExaSolverGI::mlp_add_vars()
{
	 int num_cols = CPXgetnumcols(this->cpx_env_mlp, this->cpx_lp_mlp);

	int status = 0;
	int vars_ind = 0;
	double lb[1], ub[1];
	double obj[1];
	char* varsname[1];
	char elname[1024];
	varsname[0] = elname;
	char vartype[1];

	int cnt;
	double eps = mEPS;

	/*tau*/
	this->mlp_v_tau = num_cols;
	obj[0] = 1;
	vartype[0] = 'C';
	sprintf(varsname[0], "tau");
	lb[0] = -CPX_INFBOUND;
	ub[0] = problem->getValidUB();
	//ub[0] = CPX_INFBOUND - 100;
	//ub[0] = 10000; /*this must be changed to a proper ub*/
	status = CPXnewcols(this->cpx_env_mlp, this->cpx_lp_mlp, 1, obj, lb, ub, NULL, varsname);
	if (checkCPXstatus(status, &cpx_env_mlp, &cpx_lp_mlp)) goto TERMINATE;
	++num_cols;


	/*ov_xi variables*/

	//for (int i = first_pnode; i <= last_pnode; ++i)
	for (int i = 0; i < this->problem->getNxi(); ++i)
	{


		this->mlp_v_pr_xi[i] = num_cols;

		obj[0] = 0;

		vartype[0] = 'C';

		sprintf(varsname[0], "pr_xi_i%d", i);
		lb[0] = - CPX_INFBOUND;
		ub[0] =   CPX_INFBOUND;


		status = CPXnewcols(this->cpx_env_mlp, this->cpx_lp_mlp, 1, obj, lb, ub, NULL, varsname);
		if (checkCPXstatus(status, &cpx_env_mlp, &cpx_lp_mlp)) goto TERMINATE;
		++num_cols;
	}



TERMINATE:
	if (status)
		return false;
	else
		return true;
}

bool ExaSolverGI::mlp_add_cnstr()
{
	double obj[1], rhs[1];
	int status, matbeg[1], nzc;
	char sense[1], vartype[1];
	char* cnstrname[1];
	char elname[1024];
	char ss[1024];
	char locname[1024];
	cnstrname[0] = elname;
	double lb[1];
	double ub[1];
	matbeg[0] = 0;
	double eps_coeff = 0.0001;
	status = 0;


	nzc = 0;
	rhs[0] = 0;
	sense[0] = 'L';

	vector_loc = this->problem->getUC_rhs_vect_b(); 
	matr_loc = this->problem->getUC_xicon_matr_A();

	for (int r = 0; r < this->problem->getLxi(); ++r)
	{
		sprintf(cnstrname[0], "Axi_leq_b_r%d", r);
		nzc = 0;
		rhs[0] = vector_loc[r];



		for(int i = 0; i < this->problem->getNxi(); ++i)
		{

			matind[nzc] = this->mlp_v_pr_xi[i];

			if(fabs(this->matr_loc(r, i)) > eps_coeff)
				matval[nzc] = this->matr_loc(r, i);
			else
				matval[nzc] = 0;

			++nzc;
		}

		status = CPXaddrows(cpx_env_mlp, cpx_lp_mlp, 0, 1, nzc, rhs, sense, matbeg, matind.data(), matval.data(), NULL, cnstrname);
		if (checkCPXstatus(status, &cpx_env_mlp, &cpx_lp_mlp)) goto TERMINATE;
	}


TERMINATE:
	if (status)
		return false;
	else
		return true;
}

bool ExaSolverGI::mlp_update_model()
{
	bool is_ok = true;


	is_ok = add_columns();
	if (!is_ok)
	{
		goto TERMINATE;
	}

	is_ok = add_rows();
	if (!is_ok)
	{
		goto TERMINATE;
	}


TERMINATE:
	return is_ok;
}

bool ExaSolverGI::check_violation_in_Y()
{
	/*double rhs = 0;
	
	int y_star = -1;
	double rhs_star = 0;

	for (int y = 0; y < card_Y; ++y)
		if(!y_in_mlp[y])
	{
			rhs = 0;
		for (int i = first_pnode; i <= this->last_pnode; ++i)
			if(Y(y, i) == 1)
		{
				rhs += this->p[i];
				rhs += this->mlp_curr_sol[]
		}
	}

	return false;*/
	return false;
}

bool ExaSolverGI::mlp_add_column_and_rows(int y)
{
	bool is_ok = true;


	is_ok = add_columns_y(y);
	if (!is_ok)
	{
		goto TERMINATE;
	}

	is_ok = add_rows_y(y);
	if (!is_ok)
	{
		goto TERMINATE;
	}


TERMINATE:
	return is_ok;
}

inline bool ExaSolverGI::add_columns_y(int y)
{
	 int num_cols = CPXgetnumcols(this->cpx_env_mlp, this->cpx_lp_mlp);
	
	this->y_in_mlp[y] = true;

	int status = 0;
	int vars_ind = 0;
	double lb[1], ub[1];
	double obj[1];
	char* varsname[1];
	char elname[1024];
	varsname[0] = elname;
	char vartype[1];

	int cnt;

	if (mlp_v_xi_y.getNumRows() <= y)
	{
		mlp_v_xi_y.resizeMatrix2D(2 * card_Y, this->problem->getNy()
		//	this->num_nodes
		);

		this->mlp_cnstr_tau_y.resize(2 * card_Y);
	}


	/*ov_xi variables*/
		//for (int i = first_pnode; i <= last_pnode; ++i)
	for(int i = 0; i < this->problem->getNy(); ++i)
		{


			this->mlp_v_xi_y(y, i) = num_cols;

			obj[0] = 0;

			vartype[0] = 'C';

			sprintf(varsname[0], "xi_y%d_i%d", y, i);
			lb[0] = - CPX_INFBOUND;
			ub[0] =   CPX_INFBOUND;


			status = CPXnewcols(this->cpx_env_mlp, this->cpx_lp_mlp, 1, obj, NULL, NULL, NULL, varsname);
			if (checkCPXstatus(status, &cpx_env_mlp, &cpx_lp_mlp)) goto TERMINATE;
			++num_cols;
		}

	if (num_cols >= mlp_curr_sol.size())
	{
		mlp_curr_sol.resize(2 * num_cols);
	}

TERMINATE:
	if (status)
		return false;
	else
		return true;
}

inline bool ExaSolverGI::add_rows_y(int y)
{
	double obj[1], rhs[1];
	int status, matbeg[1], nzc;
	char sense[1], vartype[1];
	char* cnstrname[1];
	char elname[1024];
	char ss[1024];
	char locname[1024];
	cnstrname[0] = elname;
	double lb[1];
	double ub[1];
	matbeg[0] = 0;
	double eps_coeff = 0.0001;
	status = 0;
	

	matr_loc   = this->problem->get2S_uobj_matr_Q();
	vector_loc = this->problem->get2S_dobj_vect_q();
	
		nzc = 0;
		rhs[0] = 0;
		sense[0] = 'L';
		sprintf(cnstrname[0], "tau_leq_Y%d", y);
		matind[nzc] = this->mlp_v_tau;
		matval[nzc] = 1;
		++nzc;
		
		
		//for (int i = this->first_pnode; i <= this->last_pnode; ++i)
		for(int i = 0; i < this->problem->getNxi(); ++i)
		{

			matind[nzc] = this->mlp_v_xi_y(y, i);
			matval[nzc] = 0;
			for (int j = 0; j < this->problem->getNy(); ++j)
			{
				matval[nzc] += -(matr_loc(i, j) * Y(y, j));
			}
			++nzc;

			
		}
		
		for (int j = 0; j < this->problem->getNy(); ++j)
			rhs[0] += Y(y, j) * vector_loc[j];


		status = CPXaddrows(cpx_env_mlp, cpx_lp_mlp, 0, 1, nzc, rhs, sense, matbeg, matind.data(), matval.data(), NULL, cnstrname);
		if (checkCPXstatus(status, &cpx_env_mlp, &cpx_lp_mlp)) goto TERMINATE;
	


		matr_loc = this->problem->getUC_xicon_matr_A();
		vector_loc = this->problem->getUC_rhs_vect_b();

	nzc = 0;
	rhs[0] = 0;
	sense[0] = 'L';
		//for (int r = 0; r < this->num_us_constr; ++r)
	for(int r = 0; r < this->problem->getLxi(); ++r)
		{
			nzc = 0;
			//rhs[0] = this->b[r];
			rhs[0] = vector_loc[r];
			sprintf(cnstrname[0], "A_leq_b_Y%d_r%d", y, r);

			//for (int i = this->first_pnode; i <= this->last_pnode; ++i)
			for(int i = 0; i < this->problem->getNxi(); ++i)
			{

				matind[nzc] = this->mlp_v_xi_y(y, i);

				//if (fabs(this->matr_loc(r, i)) > eps_coeff)
					matval[nzc] = this->matr_loc(r, i);
				//else
					matval[nzc] = 0;

				++nzc;
			}

			status = CPXaddrows(cpx_env_mlp, cpx_lp_mlp, 0, 1, nzc, rhs, sense, matbeg, matind.data(), matval.data(), NULL, cnstrname);
			if (checkCPXstatus(status, &cpx_env_mlp, &cpx_lp_mlp)) goto TERMINATE;
		}


	sense[0] = 'E';

	//	for (int i = this->first_pnode; i <= this->last_pnode; ++i)
	for(int i = 0; i < this->problem->getNw(); ++i)
			if (this->w_val[i])
			{
				sprintf(cnstrname[0], "non_ant_Y%d_i%d", y, i);

				nzc = 0;
				rhs[0] = 0;

				matind[nzc] = this->mlp_v_pr_xi[i];
				matval[nzc] = 1;
				++nzc;

				matind[nzc] = this->mlp_v_xi_y(y, i);
				matval[nzc] = -1;
				++nzc;

				status = CPXaddrows(cpx_env_mlp, cpx_lp_mlp, 0, 1, nzc, rhs, sense, matbeg, matind.data(), matval.data(), NULL, cnstrname);
				if (checkCPXstatus(status, &cpx_env_mlp, &cpx_lp_mlp)) goto TERMINATE;
			}

	++mlp_num_y;


TERMINATE:
	if (status)
		return false;
	else
		return true;
}

inline bool ExaSolverGI::add_columns()
{
	 int num_cols = CPXgetnumcols(this->cpx_env_mlp, this->cpx_lp_mlp);

	int status = 0;
	int vars_ind = 0;
	double lb[1], ub[1];
	double obj[1];
	char* varsname[1];
	char elname[1024];
	varsname[0] = elname;
	char vartype[1];

	int cnt;

	if (mlp_v_xi_y.getNumRows() <= card_Y)
	{
		mlp_v_xi_y.resizeMatrix2D(2 * card_Y, this->problem->getNxi());

		mlp_cnstr_tau_y.resize(2 * card_Y);
	}


	/*ov_xi variables*/
	for(int y = mlp_num_y; y < card_Y; ++y)
	//for (int i = first_pnode; i <= last_pnode; ++i)
		for(int i = 0; i < this->problem->getNxi(); ++i)
	{


		this->mlp_v_xi_y(y,i) = num_cols;

		obj[0] = 0;

		vartype[0] = 'C';

		sprintf(varsname[0], "xi_y%d_i%d", y,i);
		lb[0] = - CPX_INFBOUND;
		ub[0] =   CPX_INFBOUND;


		status = CPXnewcols(this->cpx_env_mlp, this->cpx_lp_mlp, 1, obj, lb, ub, NULL, varsname);
		if (checkCPXstatus(status, &cpx_env_mlp, &cpx_lp_mlp)) goto TERMINATE;
		++num_cols;
	}

	if (num_cols >= mlp_curr_sol.size())
	{
		mlp_curr_sol.resize(2 * num_cols);
	}

TERMINATE:
	if (status)
		return false;
	else
		return true;
}

inline bool ExaSolverGI::add_rows()
{
	double obj[1], rhs[1];
	int status, matbeg[1], nzc;
	char sense[1], vartype[1];
	char* cnstrname[1];
	char elname[1024];
	char ss[1024];
	char locname[1024];
	cnstrname[0] = elname;
	double lb[1];
	double ub[1];
	matbeg[0] = 0;
	double eps_coeff = 0.0001;
	status = 0;

	unsigned int num_rows = 0;

	num_rows = CPXgetnumrows(cpx_env_mlp, cpx_lp_mlp);

	this->matr_loc = this->problem->get2S_uobj_matr_Q();
	this->vector_loc = this->problem->get2S_dobj_vect_q();

	nzc = 0;
	rhs[0] = 0;
	sense[0] = 'L';
	for (int y = mlp_num_y; y < card_Y; ++y)
		{
		nzc = 0;
		rhs[0] = 0;
		sprintf(cnstrname[0], "tau_leq_Y%d", y);
		matind[nzc] = this->mlp_v_tau;
		matval[nzc] = 1;
		++nzc;

		//for (int i = this->first_pnode; i <= this->last_pnode; ++i)
		for(int i = 0; i < this->problem->getNxi(); ++i)
		{

			matind[nzc] = this->mlp_v_xi_y(y, i);
			matval[nzc] = 0;/*reset*/
			for (int j = 0; j < this->problem->getNy(); ++j)
			{
				matval[nzc] += (- matr_loc(i,j) * Y(y, j));
			}
			//matval[nzc] = -this->ppxi[i] * Y(y, i);
			++nzc;

			rhs[0] += Y(y, i) * vector_loc[i];
			//rhs[0] += Y(y, i) * this->pp[i];
		}
		/*You need this later to check tighness*/
		this->mlp_cnstr_tau_y[y] = num_rows;

		status = CPXaddrows(cpx_env_mlp, cpx_lp_mlp, 0, 1, nzc, rhs, sense, matbeg, matind.data(), matval.data(), NULL, cnstrname);
		if (checkCPXstatus(status, &cpx_env_mlp, &cpx_lp_mlp)) goto TERMINATE;
		++num_rows;
		}



	nzc = 0;
	rhs[0] = 0;
	sense[0] = 'L';

	matr_loc   = this->problem->getUC_xicon_matr_A(); 
	vector_loc = this->problem->getUC_rhs_vect_b();

for (int y = mlp_num_y; y < card_Y; ++y)
	//for (int r = 0; r < this->num_us_constr; ++r)
	for(int r = 0; r < this->problem->getLxi(); ++r)
	{
		nzc = 0;
		//rhs[0] = this->b[r];
		rhs[0] = vector_loc[r];
		sprintf(cnstrname[0], "balance_Y%d_r%d", y,r);

		//for (int i = this->first_pnode; i <= this->last_pnode; ++i)
		for(int i = 0; i < this->problem->getNxi(); ++i)
		{

			matind[nzc] = this->mlp_v_xi_y(y, i);
			matval[nzc] = matr_loc(r, i);


		/*	if (fabs(this->A(r, i)) > eps_coeff)
				matval[nzc] = this->A(r, i);
			else
				matval[nzc] = 0;*/

			++nzc;
		}

		status = CPXaddrows(cpx_env_mlp, cpx_lp_mlp, 0, 1, nzc, rhs, sense, matbeg, matind.data(), matval.data(), NULL, cnstrname);
		if (checkCPXstatus(status, &cpx_env_mlp, &cpx_lp_mlp)) goto TERMINATE;
	}


sense[0] = 'E';
for (int y = mlp_num_y; y < card_Y; ++y)
	//for (int i = this->first_pnode; i <= this->last_pnode; ++i)
	for(int i = 0; i < this->problem->getNxi(); ++i)
		if (this->w_val[i])
		{
			sprintf(cnstrname[0], "non_ant_Y%d_i%d", y, i);

			nzc = 0;
			rhs[0] = 0;
			
			matind[nzc] = this->mlp_v_pr_xi[i];
			matval[nzc] = 1;
			++nzc;

			matind[nzc] = this->mlp_v_xi_y(y, i);
			matval[nzc] = -1;
			++nzc;

			status = CPXaddrows(cpx_env_mlp, cpx_lp_mlp, 0, 1, nzc, rhs, sense, matbeg, matind.data(), matval.data(), NULL, cnstrname);
			if (checkCPXstatus(status, &cpx_env_mlp, &cpx_lp_mlp)) goto TERMINATE;
		}

	mlp_num_y = card_Y;

TERMINATE:
	if (status)
		return false;
	else
		return true;
}

bool ExaSolverGI::mlp_solve_model()
{
	bool is_ok = true;
	int status;
	int mip_status;
	CPXLONG      contextmask = 0;
	char errbuf[CPXMESSAGEBUFSIZE];
	int check_prob_stat;

	//status = CPXwriteprob(this->cpx_env_mlp, this->cpx_lp_mlp, "model_mlp.lp", "lp");
	//if (checkCPXstatus(status, &cpx_env_mlp, &cpx_lp_mlp)) { is_ok = false; goto TERMINATE; };

	
	status = CPXlpopt(this->cpx_env_mlp, this->cpx_lp_mlp);
	if (checkCPXstatus(status, &cpx_env_mlp, &cpx_lp_mlp)) goto TERMINATE;

	//check_prob_stat = CPXgetstat(this->cpx_env_mlp, this->cpx_lp_mlp);
	//
	//if (check_prob_stat == CPX_STAT_UNBOUNDED
	//	||
	//	check_prob_stat ==  CPX_STAT_INForUNBD
	//	)
	//{
	//	status = CPXgetray(cpx_env_mlp, cpx_lp_mlp, mlp_curr_sol.data());
	//	
	//	this->mlp_obj_val = DBL_MAX;
	//}else{

	status = CPXgetx(cpx_env_mlp, cpx_lp_mlp, mlp_curr_sol.data(), 0, CPXgetnumcols(cpx_env_mlp, cpx_lp_mlp) - 1);
	if (checkCPXstatus(status, &cpx_env_mlp, &cpx_lp_mlp)) goto TERMINATE;

	status = CPXgetobjval(cpx_env_mlp, cpx_lp_mlp, &this->mlp_obj_val);
	if (checkCPXstatus(status, &cpx_env_mlp, &cpx_lp_mlp)) goto TERMINATE;

	//if (fabs(this->mlp_obj_val - 1.00) < mEPS)
	//{
	//	std::cout << "Check here" << std::endl;
	//}
	//}



TERMINATE:
	return is_ok;
}

bool ExaSolverGI::subproblem_solver_update_discovery_parameters(const std::vector<bool>* w_val)
{
	bool is_ok = true;
	this->w_val = *w_val;

	is_ok = this->mlp_update_discovery_parameters();
	if (!is_ok)
	{
		goto TERMINATE;
	}
	is_ok = rop_update_discovery_parameters();
	if (!is_ok)
	{
		goto TERMINATE;
	}
	TERMINATE:
	return is_ok;
}

bool ExaSolverGI::run_exact_method(double* objective_value)
{
	bool is_ok = true;

	/*Initilizing the sep problem*/
	w_val.assign(this->problem->getNw(),
		true);
	
	is_ok = this->load_and_init(w_val);
	if (!is_ok)
	{
		goto TERMINATE;
	}

is_ok = this->bm_build_model();
if (!is_ok)
{
	goto TERMINATE;
}

is_ok = bm_solve_model();
if (!is_ok)
{
	goto TERMINATE;
}
*objective_value = this->bm_best_ub;

TERMINATE:
return is_ok;

}


inline bool ExaSolverGI::mlp_update_discovery_parameters()
{
	double rhs[1];
	int status, matbeg[1], nzc;
	char sense[1], vartype[1];
	char* cnstrname[1];
	char elname[1024];
	char ss[1024];
	char locname[1024];
	cnstrname[0] = elname;
	double lb[1];
	double ub[1];
	matbeg[0] = 0;
	double eps_coeff = 0.0001;
	status = 0;

	bool check = true;


	status = CPXdelrows(this->cpx_env_mlp, cpx_lp_mlp, 0, CPXgetnumrows(cpx_env_mlp, cpx_lp_mlp) - 1);
	if (checkCPXstatus(status, &cpx_env_mlp, &cpx_lp_mlp)) goto TERMINATE;

	/*Now we add all the constraints */
	/* first set the number of y to 0.*/
	mlp_num_y = 0;
	check = mlp_add_cnstr();
	if (!check)
	{
		std::cout << "Error in ExaSolverGI::mlp_update_discovery_parameters()" << std::endl;
		goto TERMINATE;
	}

	check = add_rows();
	if (!check)
	{
		std::cout << "Error in ExaSolverGI::mlp_update_discovery_parameters()" << std::endl;
		goto TERMINATE;
	}


	
	
TERMINATE:
	if (status || !check)
	{
		std::cout << "Error in ExaSolverGI::mlp_update_discovery_parameters()" << std::endl;
		return false;
	}
	else
		return true;
}

bool ExaSolverGI::rop_update_discovery_parameters()
{
	char lu[2];
	double bd[2];
	int cnt = 0;
	int status = 0;

	for(int i = 0; i < this->problem->getNw(); ++i)
	{
		if (this->w_val[i])
		{
			cnt = 0;

			this->matind[cnt] = this->v_gamma[i];
			lu[cnt] = 'L';
			bd[cnt] = - CPX_INFBOUND;
			++cnt;
			
			this->matind[cnt] = this->v_gamma[i];
			lu[cnt] = 'U';
			bd[cnt] =   CPX_INFBOUND;
			++cnt;
			
			status = CPXchgbds(this->cpx_env_rop, this->cpx_lp_rop, cnt, this->matind.data(), lu, bd);
			if (checkCPXstatus(status, &cpx_env_rop, &cpx_lp_rop)) goto TERMINATE;

		}
		else
		{
			cnt = 0;
			this->matind[cnt] = this->v_gamma[i];
			lu[cnt] = 'B';
			bd[cnt] = 0;
			++cnt;

			status = CPXchgbds(this->cpx_env_rop, this->cpx_lp_rop, cnt, this->matind.data(), lu, bd);
			if (checkCPXstatus(status, &cpx_env_rop, &cpx_lp_rop)) goto TERMINATE;

		}
	}


TERMINATE:
	if (status)
		return false;
	else
		return true;
}

bool ExaSolverGI::bm_build_model()
{
	bool is_ok = true;
	int status;
	int mip_status;
	CPXLONG      contextmask = 0;
	char errbuf[CPXMESSAGEBUFSIZE];
	int check_mip_stat;


	is_ok = init_cplex(&this->cpx_env_bm, &this->cpx_lp_bm);
	if (!is_ok)
	{
		goto TERMINATE;
	}

	is_ok = bm_set_cplex();
	if (!is_ok)
	{
		goto TERMINATE;
	}

	is_ok = bm_add_vars();
	if (!is_ok)
	{
		goto TERMINATE;
	}

	is_ok = bm_add_cnstr();
	if (!is_ok)
	{
		goto TERMINATE;
	}

	is_ok = bm_add_starting_cuts();
	if (!is_ok)
	{
		goto TERMINATE;
	}

	is_ok = bm_add_ccg_lb();
	if (!is_ok)
	{
		goto TERMINATE;
	}

	is_ok = set_branching_score();
	if (!is_ok)
	{
		goto TERMINATE;
	}


	if (bm_curr_sol.size() < CPXgetnumcols(this->cpx_env_bm, this->cpx_lp_bm))
	{
		bm_curr_sol.resize(CPXgetnumcols(this->cpx_env_bm, this->cpx_lp_bm));
	}



	contextmask |= CPX_CALLBACKCONTEXT_CANDIDATE;
	contextmask |= CPX_CALLBACKCONTEXT_RELAXATION;
	if (contextmask != 0) {
		//	 We are done and now we register our callback function. 
		status = CPXcallbacksetfunc(cpx_env_bm, cpx_lp_bm, contextmask, bm_general_callback, this);
		if (status != 0) {
			fprintf(stderr, "Failed to add callback: %s\n",
				CPXgeterrorstring(cpx_env_bm, status, errbuf));
			goto TERMINATE;
		}
	}
TERMINATE:
	return is_ok;
}

bool ExaSolverGI::bm_set_cplex()
{
	bool mod_stat = true;

	int status = 0;

	char errbuf[CPXMESSAGEBUFSIZE];


	status = CPXsetintparam(this->cpx_env_bm, CPX_PARAM_SCRIND, CPX_ON);
	if (checkCPXstatus(status, &cpx_env_bm, &cpx_lp_bm)) return false;
	status = CPXsetintparam(this->cpx_env_bm, CPX_PARAM_THREADS, 1);
	if (checkCPXstatus(status, &cpx_env_bm, &cpx_lp_bm)) return false;
	status = CPXchgobjsen(this->cpx_env_bm, this->cpx_lp_bm, CPX_MIN);
	if (checkCPXstatus(status, &cpx_env_bm, &cpx_lp_bm)) return false;
	status = CPXsetintparam(this->cpx_env_bm, CPX_PARAM_NUMERICALEMPHASIS, CPX_ON);
	if (checkCPXstatus(status, &cpx_env_bm, &cpx_lp_bm)) return false;
	status = CPXsetdblparam(this->cpx_env_bm, CPXPARAM_Simplex_Tolerances_Feasibility, 1e-9);
	if (checkCPXstatus(status, &cpx_env_bm, &cpx_lp_bm)) return false;
	//status = CPXsetintparam(cpx_env, CPX_PARAM_PREIND, CPX_OFF);
	//if (checkCPXstatus(status))return false;

	/*To be extra numerical accurate*/
	status = CPXsetdblparam(this->cpx_env_bm, CPXPARAM_Simplex_Tolerances_Optimality, 1e-9);
	if (checkCPXstatus(status, &cpx_env_bm, &cpx_lp_bm)) return false;

	status = CPXsetdblparam(this->cpx_env_bm, CPXPARAM_Simplex_Tolerances_Feasibility, 1e-9);
	if (checkCPXstatus(status, &cpx_env_bm, &cpx_lp_bm)) return false;

	//status = CPXsetdblparam(this->cpx_env_bm, CPXPARAM_MIP_Tolerances_MIPGap, 1e-9);
	//if (checkCPXstatus(status, &cpx_env_bm, &cpx_lp_bm)) return false;

	status = CPXsetdblparam(this->cpx_env_bm, CPXPARAM_MIP_Tolerances_Integrality, 0);
	if (checkCPXstatus(status, &cpx_env_bm, &cpx_lp_bm)) return false;

	double time_limit = this->cfg->getValueOfKey<double>("TIME_LIMIT");

	status = CPXsetdblparam(cpx_env_bm, CPXPARAM_TimeLimit, time_limit);
	if (checkCPXstatus(status, &cpx_env_bm, &cpx_lp_bm)) return false;

	/*
	status = CPXsetintparam(cpx_env, CPXPARAM_MIP_Strategy_VariableSelect, CPX_VARSEL_STRONG);
	if (checkCPXstatus(status)) return false;
	*/


	return true;
}

bool ExaSolverGI::bm_add_vars()
{

	 int num_cols = CPXgetnumcols(this->cpx_env_bm, this->cpx_lp_bm);

	int status = 0;
	int vars_ind = 0;
	double lb[1], ub[1];
	double obj[1];
	char* varsname[1];
	char elname[1024];
	varsname[0] = elname;
	char vartype[1];
	double single_node_dur_1;
	double single_node_dur_2;
	int cnt;


	/*
	Phi variable
	*/
	this->bm_v_phi = num_cols;
	obj[0] = 1;
	vartype[0] = 'C';
	sprintf(varsname[0], "Phi");
	//this->bm_best_lb;
	//lb[0] = -CPX_INFBOUND;
	lb[0] = this->problem->getValidLB();
	ub[0] = this->problem->getValidUB();
		//CPX_INFBOUND; 
	status = CPXnewcols(this->cpx_env_bm, this->cpx_lp_bm, 1, obj, lb, ub, vartype, varsname);
	if (checkCPXstatus(status, &cpx_env_bm, &cpx_lp_bm)) goto TERMINATE;
	++num_cols;

	/*
	First stage variables W: discovery variables
	*/
	/*Only profitable nodes.*/
	//r (int i = first_pnode; i <= last_pnode; ++i)
	if (this->problem->getNw() > 0)
	{
		for (int i = 0; i < this->problem->getNw(); ++i) {
			this->bm_v_w_i[i] = num_cols;
			obj[0] = this->problem->getDV_dobj_vect_d()[i];
			
			//this->disc_cost_n[i]; /*This can stay into objective.*/

			vartype[0] = 'B';
			sprintf(varsname[0], "w_n%d", i);
			if (this->problem->getDDIDversion())
			{
				lb[0] = 0;
				ub[0] = 1;

			}
			else {
				lb[0] = 1;
				ub[0] = 1;
			}


			status = CPXnewcols(this->cpx_env_bm, this->cpx_lp_bm, 1, obj, lb, ub, vartype, varsname);
			if (checkCPXstatus(status, &cpx_env_bm, &cpx_lp_bm)) goto TERMINATE;

			++num_cols;
		}
	}

	

TERMINATE:
	if (status)
		return false;
	else
		return true;
}

bool ExaSolverGI::bm_add_cnstr()
{
	double obj[1], rhs[1];
	int status, matbeg[1], nzc;
	char sense[1], vartype[1];
	char* cnstrname[1];
	char elname[1024];
	char ss[1024];
	char locname[1024];
	cnstrname[0] = elname;
	double lb[1];
	double ub[1];
	matbeg[0] = 0;
	double coeff_eps = mEPS;

	if (matind.size() < CPXgetnumcols(cpx_env_bm, cpx_lp_bm))
	{
		matind.resize(CPXgetnumcols(cpx_env_bm, cpx_lp_bm));
		matval.resize(CPXgetnumcols(cpx_env_bm, cpx_lp_bm));

	}

	// W matrices constraints
	this->vector_loc = this->problem->getDV_drhs_vect_bw();
	this->matr_loc = this->problem->getDV_dcon_matr_W();


	for (unsigned short row = 0; row < this->problem->getLw(); ++row)
	{
		rhs[0] = vector_loc[row];
		sense[0] = 'L';
		nzc = 0;
		sprintf(cnstrname[0], "Ww_%d", row);
		// Now we go trought the "columns", coefficient of the matrices 

		for (unsigned short cols = 0; cols < this->problem->getNw(); ++cols)
			if (fabs(matr_loc(row, cols)) > coeff_eps)
			{
				this->matind[nzc] = this->bm_v_w_i[cols];
				this->matval[nzc] = matr_loc(row, cols);
				++nzc;
			}
		status = CPXaddrows(cpx_env_bm, cpx_lp_bm, 0, 1, nzc, rhs, sense, matbeg, matind.data(), matval.data(), NULL, cnstrname);
		if (checkCPXstatus(status, &cpx_env_bm, &cpx_lp_bm)) goto TERMINATE;

	}




TERMINATE:
	if (status)
		return false;
	else
		return true;
}

bool ExaSolverGI::set_branching_score()
{

	//CPXPARAM_MIP_Strategy_Order

	bm_branch_priority_score.assign(CPXgetnumcols(cpx_env_bm, cpx_lp_bm), 0);
	int status = 0;
	int	cnt = 0;
	for (int i = 0; i < this->problem->getNw(); ++i)
			//for (int i = first_pnode; i <= last_pnode; ++i)
			{
				matind[cnt] = bm_v_w_i[i];
				bm_branch_priority_score[cnt] = 1;
				++cnt;

	
			}

		status = CPXcopyorder(this->cpx_env_bm, this->cpx_lp_bm, cnt, matind.data(), this->bm_branch_priority_score.data(),
			CPX_BRANCH_GLOBAL);
		if (checkCPXstatus(status, &cpx_env_bm, &cpx_lp_bm)) goto TERMINATE;

	TERMINATE:
		if (status)
			return false;
		else
			return true;
	
}

bool ExaSolverGI::bm_solve_model()
{
	bool is_ok = true;
	int status;
	int mip_status;
	CPXLONG      contextmask = 0;
	char errbuf[CPXMESSAGEBUFSIZE];
	int check_mip_stat;

	double start = 0;

	double end = 0;

	//status = CPXwriteprob(this->cpx_env_bm, this->cpx_lp_bm, "masterSPP.lp", "lp");
	//if (checkCPXstatus(status, &cpx_env_bm, &cpx_lp_bm)) { is_ok = false; goto TERMINATE; };


	start = clock();
	status = CPXmipopt(this->cpx_env_bm, this->cpx_lp_bm);
	if (checkCPXstatus(status, &cpx_env_bm, &cpx_lp_bm)) goto TERMINATE;
	end = clock();
	this->bm_branch_and_cut_time = (double)(end - start) / (double)CLK_TCK;


	check_mip_stat = CPXgetstat(cpx_env_bm, cpx_lp_bm);

	if (bm_curr_sol.size() < CPXgetnumcols(cpx_env_bm, cpx_lp_bm))
	{
		bm_curr_sol.resize(CPXgetnumcols(cpx_env_bm, cpx_lp_bm));
	}

	if (check_mip_stat != CPXMIP_TIME_LIM_INFEAS)
	{
		status = CPXgetobjval(cpx_env_bm, cpx_lp_bm, &this->bm_best_ub);
		if (checkCPXstatus(status, &cpx_env_bm, &cpx_lp_bm)) goto TERMINATE;

		status = CPXgetx(cpx_env_bm, cpx_lp_bm, bm_curr_sol.data(), 0, CPXgetnumcols(cpx_env_bm, cpx_lp_bm) - 1);
		if (checkCPXstatus(status, &cpx_env_bm, &cpx_lp_bm)) goto TERMINATE;

	}
	else
	{
		this->bm_best_ub = 0;

	}


	status = CPXgetbestobjval(cpx_env_bm, cpx_lp_bm, &this->bm_best_lb);
	if (checkCPXstatus(status, &cpx_env_bm, &cpx_lp_bm)) goto TERMINATE;



	status = CPXgetmiprelgap(cpx_env_bm, cpx_lp_bm, &this->cplex_gap);
	if (checkCPXstatus(status, &cpx_env_bm, &cpx_lp_bm)) goto TERMINATE;

	this->bandc_nodes = CPXgetnodecnt(cpx_env_bm, cpx_lp_bm);

	/*
	if (bandc_nodes == 0 )
	{
		this->root_node_lb = bm_best_ub;
	}
	*/

TERMINATE:
	return is_ok;
}

bool ExaSolverGI::bm_add_starting_cuts()
{
	double rhs[1];
	int status, matbeg[1], nzc;
	char sense[1], vartype[1];
	char* cnstrname[1];
	char elname[1024];
	char ss[1024];
	char locname[1024];
	cnstrname[0] = elname;
	double lb[1];
	double ub[1];
	matbeg[0] = 0;
	bool is_ok;
	double phi_val = 0;
	//uint64_t curr_cut = 0;

	bitset<MAX_NUM_CUTS_ELEM> curr_cut (0);
	
	double start = 0;
	double end = 0;

	status = false;

	start = clock();

		/*Full discovery cut*/
		//curr_cut = 0;
		curr_cut.reset();
		//for (int i = first_pnode; i <= last_pnode; ++i)
		for(int i = 0; i < this->problem->getNw(); ++i)
		{
			this->w_val[i] = true;


			if (this->w_val[i])
			{
				this->bm_add_elem_to_logic_cut(&curr_cut, i);
			}
		}


		is_ok = this->subproblem_solver_update_discovery_parameters(&w_val);
		if (!is_ok)
		{
			goto TERMINATE;
		}

		is_ok = this->run_sub_problem_solver(&phi_val);
		if (!is_ok)
		{
			goto TERMINATE;
		}
		
		this->bm_add_logic_cut_to_pool(&curr_cut, &phi_val);
		curr_cut.reset();
		this->bm_cut_lb = phi_val;
		this->bm_best_lb = phi_val;


		/*Now we save xi*/
		if (this->use_xi_bar_scenarios)
			//(ccg_set_xib.getNumRows() > 0)
		{
			store_xbar();
		}
		//else {
			nzc = 0;
			sense[0] = 'G';
			rhs[0] = phi_val;

			this->matind[nzc] = this->bm_v_phi;
			this->matval[nzc] = 1;
			++nzc;
			sprintf(cnstrname[0], "Full_discovery_cut");

			//for (int i = first_pnode; i <= last_pnode; ++i)
			for (int i = 0; i < this->problem->getNw(); ++i)
				if (!w_val[i])
				{
					this->matind[nzc] = this->bm_v_w_i[i];
					this->matval[nzc] = (phi_val - this->bm_cut_lb);
					++nzc;

				}

			if (!this->bm_use_info_cuts)
				for (int i = 0; i < this->problem->getNw(); ++i)
					if (w_val[i])
					{
						rhs[0] += this->bm_cut_lb - phi_val;


						this->matind[nzc] = this->bm_v_w_i[i];
						this->matval[nzc] = this->bm_cut_lb - phi_val;
						++nzc;
					}

			


			status = CPXaddrows(cpx_env_bm, cpx_lp_bm, 0, 1, nzc, rhs, sense, matbeg, matind.data(), matval.data(), NULL, cnstrname);
			if (checkCPXstatus(status, &cpx_env_bm, &cpx_lp_bm)) goto TERMINATE;
		//}




		if (this->problem->getDDIDversion()
			) {
			//curr_cut = 0;
			curr_cut.reset();

			/*First 0 sensor*/
			//for (int i = first_pnode; i <= last_pnode; ++i)
			for (int i = 0; i < this->problem->getNw(); ++i)
			{
				this->w_val[i] = false;

				if (this->w_val[i])
				{
					this->bm_add_elem_to_logic_cut(&curr_cut, i);
				}


			}
			is_ok = this->subproblem_solver_update_discovery_parameters(&w_val);
			if (!is_ok)
			{
				goto TERMINATE;
			}

			is_ok = this->run_sub_problem_solver(&phi_val);
			if (!is_ok)
			{
				goto TERMINATE;
			}


			this->bm_add_logic_cut_to_pool(&curr_cut, &phi_val);

			/*Now we save xi*/
			if (this->use_xi_bar_scenarios)
			{
				store_xbar();
			}
			else {
/*If we don't warm start wih the scenarios, we use warm start with a cut*/
			nzc = 0;
			sense[0] = 'G';
			rhs[0] = phi_val;


			this->matind[nzc] = this->bm_v_phi;
			this->matval[nzc] = 1;
			++nzc;
			sprintf(cnstrname[0], "No_discovery_cut");

			//for (int i = first_pnode; i <= last_pnode; ++i)
			for (int i = 0; i < this->problem->getNw(); ++i)
				if (!w_val[i])
				{

					this->matind[nzc] = this->bm_v_w_i[i];
					this->matval[nzc] = (phi_val - this->bm_cut_lb);
					++nzc;


				}

			if (!this->bm_use_info_cuts)
				for (int i = 0; i < this->problem->getNw(); ++i)
					if (w_val[i])
					{
						rhs[0] += this->bm_cut_lb - phi_val;


						this->matind[nzc] = this->bm_v_w_i[i];
						this->matval[nzc] = this->bm_cut_lb - phi_val;
						++nzc;
					}


			status = CPXaddrows(cpx_env_bm, cpx_lp_bm, 0, 1, nzc, rhs, sense, matbeg, matind.data(), matval.data(), NULL, cnstrname);
			if (checkCPXstatus(status, &cpx_env_bm, &cpx_lp_bm)) goto TERMINATE;
		}


			/*Cardinality Max -1*/
			/*These cuts didn't help much*/
			if (false) {
				for (int istar = 0; istar < this->problem->getNw(); ++istar)
				{
					//curr_cut = 0;
					curr_cut.reset();
					//for (int i = first_pnode; i <= last_pnode; ++i)
					for (int i = 0; i < this->problem->getNw(); ++i)
					{
						if (i == istar)
						{
							this->w_val[i] = false;
						}
						else {
							this->w_val[i] = true;
						}

						if (this->w_val[i])
						{
							this->bm_add_elem_to_logic_cut(&curr_cut, i);
						}

					}

					is_ok = this->subproblem_solver_update_discovery_parameters(&w_val);
					if (!is_ok)
					{
						goto TERMINATE;
					}

					is_ok = this->run_sub_problem_solver(&phi_val);
					if (!is_ok)
					{
						goto TERMINATE;
					}

					/*Now we save xi*/
					if (this->use_xi_bar_scenarios)
						//if (ccg_set_xib.getNumRows() > 0)
					{
						store_xbar();
					}

					this->bm_add_logic_cut_to_pool(&curr_cut, &phi_val);
					//curr_cut = 0;
					curr_cut.reset();

					nzc = 0;
					sense[0] = 'G';
					rhs[0] = phi_val;

					this->matind[nzc] = this->bm_v_phi;
					this->matval[nzc] = 1;
					++nzc;
					sprintf(cnstrname[0], "CrdMaxminusOne_discovery_cut_i%d", istar);

					//for (int i = first_pnode; i <= last_pnode; ++i)
					for (int i = 0; i < this->problem->getNw(); ++i)
						if (!w_val[i])
						{

							this->matind[nzc] = this->bm_v_w_i[i];
							//this->matval[nzc] = -phi_val;
							this->matval[nzc] = (phi_val - this->bm_cut_lb);
							++nzc;


							/*	this->matind[nzc] = this->bm_v_w_i[i];
								this->matval[nzc] = -phi_val;
								++nzc;*/
						}
					status = CPXaddrows(cpx_env_bm, cpx_lp_bm, 0, 1, nzc, rhs, sense, matbeg, matind.data(), matval.data(), NULL, cnstrname);
					if (checkCPXstatus(status, &cpx_env_bm, &cpx_lp_bm)) goto TERMINATE;
				}

				/*Cardinality 1*/
				for (int istar = 0; istar < this->problem->getNw(); ++istar)
				{
					//curr_cut = 0;
					curr_cut.reset();

					for (int i = 0; i < this->problem->getNw(); ++i)
					{
						if (i == istar)
						{
							this->w_val[i] = true;
						}
						else {
							this->w_val[i] = false;
						}


						if (this->w_val[i])
						{
							this->bm_add_elem_to_logic_cut(&curr_cut, i);
						}

					}

					is_ok = this->subproblem_solver_update_discovery_parameters(&w_val);
					if (!is_ok)
					{
						goto TERMINATE;
					}

					is_ok = this->run_sub_problem_solver(&phi_val);
					if (!is_ok)
					{
						goto TERMINATE;
					}

					/*Now we save xi*/
					if (this->use_xi_bar_scenarios)
					{
						store_xbar();
					}


					this->bm_add_logic_cut_to_pool(&curr_cut, &phi_val);


					nzc = 0;
					sense[0] = 'G';
					rhs[0] = phi_val;

					this->matind[nzc] = this->bm_v_phi;
					this->matval[nzc] = 1;
					++nzc;
					sprintf(cnstrname[0], "CrdMaxmOne_discovery_cut_i%d", istar);

					for (int i = 0; i < this->problem->getNw(); ++i)
						if (!w_val[i])
						{

							this->matind[nzc] = this->bm_v_w_i[i];
							this->matval[nzc] = (phi_val - this->bm_cut_lb);
							++nzc;

						
						}
					status = CPXaddrows(cpx_env_bm, cpx_lp_bm, 0, 1, nzc, rhs, sense, matbeg, matind.data(), matval.data(), NULL, cnstrname);
					if (checkCPXstatus(status, &cpx_env_bm, &cpx_lp_bm)) goto TERMINATE;
				}
			}

}
		end = clock();
		this->bm_starting_cuts_time = (double)(end - start) / (double)CLK_TCK;


TERMINATE:
	if (status || !is_ok)
		return false;
	else
		return true;
}

inline int CPXPUBLIC ExaSolverGI::bm_general_callback(CPXCALLBACKCONTEXTptr context, CPXLONG contextid, void* userdata)
{
	bool cut_added = false;
	bool w_integer = false;
	int status = 0;
	bool is_ok = true;
	ExaSolverGI* solver = (ExaSolverGI*)userdata;
	double val;
	double treshold;


	status = CPXcallbackgetinfodbl(context, CPXCALLBACKINFO_BEST_BND, &solver->bm_best_lb);
		//(solver->cpx_env_bm, solver->cpx_lp_bm, &solver->bm_best_lb);
	if (status != 0)
	{
		solver->checkCPXstatus(status, &solver->cpx_env_bm, &solver->cpx_lp_bm);

		return status;
	}
	/*Update the bound for the cuts*/
	if(solver->bm_best_lb > solver->bm_cut_lb)
		solver->bm_cut_lb = solver->bm_best_lb;


	if (contextid == CPX_CALLBACKCONTEXT_CANDIDATE)
	{

		status = CPXcallbackgetcandidatepoint(context, solver->bm_curr_sol.data(), 0, CPXgetnumcols(solver->cpx_env_bm, solver->cpx_lp_bm) - 1, &val);
		if (status != 0)
		{
			solver->checkCPXstatus(status, &solver->cpx_env_bm, &solver->cpx_lp_bm);

			return status;
		}


	

		if (!cut_added)
		{
			is_ok = solver->bm_separate_logic_benders_cut(context, contextid, NULL, &cut_added);
		}
	

		


	}
	else
		if (contextid == CPX_CALLBACKCONTEXT_RELAXATION)
		{


			/*Only CPX 12.10*/
			CPXLONG node_id, node_depth;

			CPXcallbackgetinfolong(context, CPXCALLBACKINFO_NODEUID, &node_id);
			if (status != 0)
			{
				solver->checkCPXstatus(status, &solver->cpx_env_bm, &solver->cpx_lp_bm);			return status;
			}

			CPXcallbackgetinfolong(context, CPXCALLBACKINFO_NODEDEPTH, &node_depth);
			if (status != 0)
			{
				solver->checkCPXstatus(status, &solver->cpx_env_bm, &solver->cpx_lp_bm);			return status;
			}



			status = CPXcallbackgetrelaxationpoint(context, solver->bm_curr_sol.data(), 0, CPXgetnumcols(solver->cpx_env_bm, solver->cpx_lp_bm) - 1, &val);
			if (status != 0)
			{
				solver->checkCPXstatus(status, &solver->cpx_env_bm, &solver->cpx_lp_bm);				return status;
			}


			/*if w are integer, we add the benders cut*/
	/*		w_integer = true;
			for(int i = 0; i < solver->problem->getLw(); ++i)
			{
				if (
					(solver->bm_curr_sol[solver->bm_v_w_i[i]] > mEPS)
					&&
					(1.00 - solver->bm_curr_sol[solver->bm_v_w_i[i]] > mEPS)
					)
				{
					w_integer = false;
					break;
				}

			}*/



			/*cut per node if fractional.*/
		/*	if (!w_integer)
			{*/
				if (node_id == solver->bm_curr_node)
				{

					solver->bm_curr_node = node_id;

					if (
						(node_depth <= 1
							&&
							/*		fabs(val - solver->bm_lp_prev_iter) < (EXA_Benders_MIN_COEFF_LP_IMPR * fabs(solver->bm_lp_prev_iter))
									&&*/
							solver->bm_num_iter_curr_node >= EXA_Benders_MAX_NUM_ITER_NODE_DEPTH_1)
						||
						(node_depth > 1
							&&
							solver->bm_num_iter_curr_node >= EXA_Benders_MAX_NUM_ITER_NODE)
						)
					{

						return 0; // Branch.
					}
					++solver->bm_num_iter_curr_node;

				}
				else {

					solver->bm_reset_branc_and_cut_indicators();
					solver->bm_curr_node = node_id;
					solver->bm_lp_prev_iter = val;
				}

				/*If we didn't quit, we add the cut.*/
				is_ok = solver->bm_separate_logic_benders_cut(context, contextid, NULL, &cut_added);
		//}
		//	else {
		//		/*We have integer w. This can happen only if we use scenarios.*/
		//		is_ok = solver->bm_separate_logic_benders_cut(context, contextid, NULL, &cut_added);
		//	}

					


		}


	return status;

}

inline bool ExaSolverGI::bm_separate_logic_benders_cut(CPXCALLBACKCONTEXTptr context, CPXLONG contextid, double var_treshold, bool* cut_added)
{
	bool is_ok = true;
	double obj = 0;

	int status = 0;

	double rhs[1];
	int matbeg[1], nzc;
	char sense[1];
	int force[1];
	int local[1];
	force[0] = CPX_USECUT_FORCE;
	local[0] = 0;
	matbeg[0] = 0;

	bool cut_is_old = false;
	int cut_index;
	//uint64_t curr_cut = 0;
	bitset<MAX_NUM_CUTS_ELEM> curr_cut(0);
	*cut_added = false;


	for (int i = 0; i < this->problem->getNw(); ++i)
		{

			if (this->bm_curr_sol[this->bm_v_w_i[i]] > mEPS)
			{
				this->w_val[i] = true;
			}
			else
			{
				this->w_val[i] = false;
			}


			if (this->w_val[i])
			{
				this->bm_add_elem_to_logic_cut(&curr_cut, i);
			}
		}
		

		cut_is_old = this->bm_is_logic_cut_already_added(&curr_cut, &cut_index);
		if (cut_is_old)
		{
			obj = this->bm_logic_cuts_pool_phi_val[cut_index];
		}
		else {


			is_ok = this->subproblem_solver_update_discovery_parameters(&w_val);
			if (!is_ok)
			{
				goto TERMINATE;
			}

			is_ok = this->run_sub_problem_solver(&obj);
			if (!is_ok)
			{
				goto TERMINATE;
			}

			this->bm_add_logic_cut_to_pool(&curr_cut, &obj);
		}

	if (
		((obj - this->bm_curr_sol[this->bm_v_phi]) / (obj + EXA_EPS)) > EXA_EPS
		//((obj - this->bm_curr_sol[this->bm_v_phi]) / (obj + EXA_EPS)) < EXA_EPS
		//obj > this->bm_curr_sol[this->bm_v_phi] + EXA_EPS
		&&
		///*Second condition to avoid problems with cplex 12.10 rejecting candidate.*/
		!(cut_is_old 
			&& 
			//((obj - this->bm_curr_sol[this->bm_v_phi]) / (obj + EXA_EPS)) < EXA_EPS
			((obj - this->bm_curr_sol[this->bm_v_phi]) / (obj + EXA_EPS)) < EXA_EPS
			&&
			contextid == CPX_CALLBACKCONTEXT_CANDIDATE
			)
		)
	{

		nzc = 0;
		sense[0] = 'G';
		rhs[0] = obj;
		

		this->matind[nzc] = this->bm_v_phi;
		this->matval[nzc] = 1;
		++nzc;

		for (int i = 0; i < this->problem->getNw(); ++i)
			if(!w_val[i])
		{

				this->matind[nzc] = this->bm_v_w_i[i];
				this->matval[nzc] = (obj - this->bm_best_lb); /*Updating bound*/
				++nzc;

		
		
		}

		/*If we don't use Information cuts we add this part as well.*/
		if (!this->bm_use_info_cuts)
			for (int i = 0; i < this->problem->getNw(); ++i)
				if (w_val[i])
				{
					rhs[0] += this->bm_cut_lb - obj;


					this->matind[nzc] = this->bm_v_w_i[i];
					this->matval[nzc] = this->bm_cut_lb - obj;
					++nzc;
				}


		if (contextid == CPX_CALLBACKCONTEXT_RELAXATION)
		{


			status = CPXcallbackaddusercuts(context, 1, nzc, rhs,
				sense, matbeg, matind.data(), matval.data(), force, local);
			if (checkCPXstatus(status, &cpx_env_bm, &cpx_lp_bm)) return false;

			++num_logic_benders_cuts;
			*cut_added = true;
		}
		else
			if (contextid == CPX_CALLBACKCONTEXT_CANDIDATE)
			{
				status = CPXcallbackrejectcandidate(context, 1, nzc, rhs,
					sense, matbeg, matind.data(), matval.data());
				if (checkCPXstatus(status, &cpx_env_bm, &cpx_lp_bm)) return false;
				*cut_added = true;

				++num_logic_benders_cuts;

			}
	}
	else
	{
		if (contextid == CPX_CALLBACKCONTEXT_CANDIDATE)
		{

			

			this->store_set_of_tight_y();
		}

	}

TERMINATE:

	return true;
}

inline void ExaSolverGI::bm_add_elem_to_logic_cut(bitset<MAX_NUM_CUTS_ELEM>* set, unsigned int elem)
{
	*set |= (unsigned int)pow(2.00, (double)(elem - 1));
}

inline bool ExaSolverGI::bm_is_logic_cut_already_added(bitset<MAX_NUM_CUTS_ELEM>* cut, int* index)
{
	for (int i = 0; i < bm_num_stored_logic_cuts; ++i)
		if (*cut == bm_logic_cuts_pool_set[i])
		{
			*index = i;

			return true;
		}

	return false;
}

inline void ExaSolverGI::bm_add_logic_cut_to_pool(bitset<MAX_NUM_CUTS_ELEM>* cut, const double* phi)
{

	if (bm_num_stored_logic_cuts >= bm_logic_cuts_pool_set.size())
	{

		bm_logic_cuts_pool_set.resize(2 * bm_num_stored_logic_cuts);
		bm_logic_cuts_pool_phi_val.resize(2 * bm_num_stored_logic_cuts);

	}
	bm_logic_cuts_pool_set[bm_num_stored_logic_cuts] = *cut;
	bm_logic_cuts_pool_phi_val[bm_num_stored_logic_cuts] = *phi;
	++bm_num_stored_logic_cuts;

}


