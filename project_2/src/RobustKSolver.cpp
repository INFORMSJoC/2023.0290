#include "RobustKSolver.h"


bool RobustKSolver::init_data_structures()
{
	unsigned int num_var = 0;

	num_var = /*PHI*/ 1 + 
		this->problem->getNx() + (this->problem->getNy() * this->problem->getK() + this->problem->getNx()+ 
		this->problem->getNw())* this->problem->getK() +
		this->problem->getNw() + this->problem->getK()/*alpha_k*/ + 
		this->problem->getNxi() + this->problem->getNxi() * this->problem->getK()  /*Beta and beta k*/+
		+ 2 * (this->problem->getNxi() * this->problem->getK()) /*gamma and gamma_k*/;

	/*matind.resize(MAX_NUM_VARS);
	matval.resize(MAX_NUM_VARS);*/

	//matind.resize(num_var);
	//matval.resize(num_var);
	//master_curr_sol.resize(num_var);
	//slave_curr_sol.resize(num_var);

	matind.resize(MAX_NUM_VARS);
	matval.resize(MAX_NUM_VARS);
	master_curr_sol.resize(MAX_NUM_VARS);
	slave_curr_sol.resize(MAX_NUM_VARS);

	
	this->v_dual_set.resize(1);

	v_dual_set[0].z_k.assignMatrix2D(this->problem->getK(), this->problem->getL(), 0);

	//v_dual_set[1].z_k.assignMatrix2D(this->problem->getK(), this->problem->getL(), 1);

	//v_dual_set[2].z_k.assignMatrix2D(this->problem->getK(), this->problem->getL(), 1);
	//v_dual_set[2].z_k(0, 0) = 1;
	//v_dual_set[2].z_k(1, 0) = 0;

	//v_dual_set[3].z_k.assignMatrix2D(this->problem->getK(), this->problem->getL(), 1);
	//v_dual_set[3].z_k(0, 0) = 0;
	//v_dual_set[3].z_k(1, 0) = 1;


	return true;
}

int RobustKSolver::checkCPXstatus(CPXENVptr env, int status)
{
		char errmsg[CPXMESSAGEBUFSIZE];
		if (status == 0) return 0;

		CPXgeterrorstring(env, status, errmsg);
		printf(" %s \n", errmsg);
		return status;
	
}

bool RobustKSolver::free_cplex(CPXENVptr* env, CPXLPptr* lp)
{
	
		int status = 0;
		status = CPXfreeprob(*env, lp);
		if (checkCPXstatus(*env, status)) return false;
		status = CPXcloseCPLEX(env);
		if (checkCPXstatus(*env, status)) return false;

		return true;
	
}

bool RobustKSolver::free_cplex_problem(CPXENVptr* env, CPXLPptr* lp)
{

	int status = 0;
	status = CPXfreeprob(*env, lp);
	if (checkCPXstatus(*env, status)) return false;


	return true;
}

bool RobustKSolver::init_cplex_data_str()
{
	int status = 0;
	
	this->cpx_env_master = CPXopenCPLEX(&status);
	if (this->checkCPXstatus(this->cpx_env_master, status)) return false;
	this->cpx_lp_master = CPXcreateprob(this->cpx_env_master, &status,this->problem->getProbName().c_str());
	if (checkCPXstatus(this->cpx_env_master, status))return false;

	return true;
}

bool RobustKSolver::set_cplex_parameters()
{
	bool mod_stat = true;

	int status = 0;

	char errbuf[CPXMESSAGEBUFSIZE];


	status = CPXsetintparam(this->cpx_env_master, CPX_PARAM_SCRIND, CPX_ON);
	if (checkCPXstatus(this->cpx_env_master, status)) return false;
	status = CPXsetintparam(this->cpx_env_master, CPX_PARAM_THREADS, 1);
	if (checkCPXstatus(this->cpx_env_master, status)) return false;
	status = CPXchgobjsen(this->cpx_env_master, this->cpx_lp_master, CPX_MIN);
	if (checkCPXstatus(this->cpx_env_master, status)) return false;
	status = CPXsetintparam(this->cpx_env_master, CPX_PARAM_NUMERICALEMPHASIS, CPX_ON);
	if (checkCPXstatus(this->cpx_env_master, status)) return false;
	status = CPXsetdblparam(this->cpx_env_master, CPXPARAM_Simplex_Tolerances_Feasibility, 1e-9);
	if (checkCPXstatus(this->cpx_env_master, status)) return false;

	status = CPXsetdblparam(this->cpx_env_master, CPXPARAM_Simplex_Tolerances_Optimality, 1e-9);
	if (checkCPXstatus(this->cpx_env_master, status)) return false;
	status = CPXsetdblparam(this->cpx_env_master, CPXPARAM_MIP_Tolerances_MIPGap, 1e-9);
	if (checkCPXstatus(this->cpx_env_master, status)) return false;
	status = CPXsetdblparam(this->cpx_env_master, CPXPARAM_MIP_Tolerances_Integrality, 0);
	if (checkCPXstatus(this->cpx_env_master, status)) return false;
	

	status = CPXsetdblparam(this->cpx_env_master, CPXPARAM_TimeLimit, TIME_LIMIT);
	if (checkCPXstatus(this->cpx_env_master, status)) return false;


	return true;



}

bool RobustKSolver::set_configuration_file(char* config_file)
{
	
	this->cfg = new Configuration::ConfigFile(config_file);

	return true;
	
}

void RobustKSolver::del_config_file()
{
	delete this->cfg;
}

bool RobustKSolver::column_and_constraints_generation(double* obj_val)
{
	bool check = true;
	int status = 0;

	char errbuf[CPXMESSAGEBUFSIZE];

	int check_mip_stat;

	double start = 0;

	double end = 0;

	double elapsed_time;

	double master_obj_val = 0;
	double master_unc_cost = 0;
	double slave_obj_val = 0;

	check = build_master_model();
	if (!check)
	{
		goto TERMINATE;
	}

	do
	{
		check = master_solve_problem(&master_obj_val, &master_unc_cost);
		if (!check)
		{
			goto TERMINATE;
		}

		check = slave_build_problem();
		if (!check)
		{
			goto TERMINATE;
		}

		check = slave_solve_problem(&slave_obj_val);
		if (!check)
		{
			goto TERMINATE;
		}

		if (slave_obj_val > master_unc_cost)
		{
			check = slave_store_z_values();
			if (!check)
			{
				goto TERMINATE;
			}

			check = master_update_problem(v_dual_set.size() - 1);
			if (!check)
			{
				goto TERMINATE;
			}
		}
		/*Free the slave*/
		check = free_cplex(&cpx_env_slave, &cpx_lp_slave);
		if (!check)
		{
			goto TERMINATE;
		}


	}
	//while (true);
	while (master_unc_cost + COEFF_EPS < slave_obj_val);


	check = free_cplex(&cpx_env_master, &cpx_lp_master);
	if (!check)
	{
		goto TERMINATE;
	}

	*obj_val = master_obj_val;

TERMINATE:
	if (!check)
		return false;
	else
		return true;

}

bool RobustKSolver::solve_problem()
{
	bool check = true;
	int status = 0;

	char errbuf[CPXMESSAGEBUFSIZE];

	int check_mip_stat;

	double start = 0;

	double end = 0;

	double elapsed_time;

	double slave_obj_val = 0;

	check = build_master_model();
	if (!check)
	{
		goto TERMINATE;
	}

	start = clock();
	status = CPXmipopt(this->cpx_env_master, this->cpx_lp_master);
	if (checkCPXstatus(this->cpx_env_master, status)) goto TERMINATE;
	end = clock();
	this->comp_time = (double)(end - start) / (double)CLK_TCK;


	status = CPXgetobjval(this->cpx_env_master, this->cpx_lp_master, &this->bub);
	if (checkCPXstatus(this->cpx_env_master, status)) goto TERMINATE;

	status =CPXgetbestobjval(this->cpx_env_master, this->cpx_lp_master, &this->blb);
	if (checkCPXstatus(this->cpx_env_master, status)) goto TERMINATE;

	status = CPXgetmiprelgap(this->cpx_env_master, this->cpx_lp_master, &this->cplex_gap);
	if (checkCPXstatus(this->cpx_env_master, status)) goto TERMINATE;


	status = CPXgetx(this->cpx_env_master, this->cpx_lp_master, this->master_curr_sol.data(),
		0, CPXgetnumcols(cpx_env_master, cpx_lp_master)-1);
	if (checkCPXstatus(this->cpx_env_master, status)) goto TERMINATE;

	this->free_cplex(&this->cpx_env_master, &this->cpx_lp_master);


TERMINATE:
	if (!check)
		return false;
	else
		return true;
}

bool RobustKSolver::appr_solver(double* obj_val)
{
	bool check = true;
	int status = 0;

	char errbuf[CPXMESSAGEBUFSIZE];

	int check_mip_stat;

	double start = 0;

	double end = 0;

	double elapsed_time;



	check = appr_build_model();
	if (!check)
	{
		goto TERMINATE;
	}

check = appr_solve_problem(obj_val);
		if (!check)
		{
			goto TERMINATE;
		}


	TERMINATE:

		if (status)
		{
			return false;
		}

		return true;
}

void RobustKSolver::set_problem_info(K_Adapt_Problem* problem)
{
	this->problem = problem;
}

bool RobustKSolver::build_master_model()
{
	bool check = true;
	int status;
	bool use_optimistic_cuts = cfg->getValueOfKey<bool>("USE_OPTIMISTIC_CUTS");
	bool use_rlt_cuts = cfg->getValueOfKey<bool>("USE_RLT_CUTS");


	/* Data structures*/
	check = this->init_data_structures();
	if (!check)
	{
		goto TERMINATE;
	}

	check = this->init_cplex_data_str();
	if (!check)
	{
		goto TERMINATE;
	}

	check = this->set_cplex_parameters();
	if (!check)
	{
		goto TERMINATE;
	}

	/*
	Variables
	*/
	check = this->add_primal_variables();
	if (!check)
	{
		goto TERMINATE;
	}

	for (int z = 0; z < v_dual_set.size(); ++z)
	{
		check = this->add_dualization_variables(z);
	if (!check)
	{
		goto TERMINATE;
	}

	check = this->add_bilinear_variables(z);
	if (!check)
	{
		goto TERMINATE;
	}
	}
	
	

	/*
	Constraints
	*/
	check = this->add_deterministic_constraints();
	if (!check)
	{
		goto TERMINATE;
	}



	if (use_optimistic_cuts)
		{

			check = add_primal_valid_inequalities();
			if (!check)
			{
			goto TERMINATE;
			}
		}
	
	
	for (int z = 0; z < v_dual_set.size(); ++z)
	{
		check = this->add_dualization_constraints(z);
		if (!check)
		{
			goto TERMINATE;
		}

		check = this->add_linearization_constraints(z);
		if (!check)
		{
			goto TERMINATE;
		}


		if (z == 0)
		{
			check = this->add_dual_vars_valid_inequalities(z);
			if (!check)
			{
				goto TERMINATE;
			}
		}
		if (use_rlt_cuts)
		{
			check = this->add_RLT_AlphaXWY_valid_inequalities(z);
			if (!check)
			{
				goto TERMINATE;
			}

			check = add_RLT_AlphaExwy_valid_inequalities(z);
			if (!check)
			{
				goto TERMINATE;
			}


			check = add_RLT_Primal_valid_inequalities(z);
			if (!check)
			{
				goto TERMINATE;
			}
		}
	}



	//status = CPXwriteprob(this->cpx_env_master, this->cpx_lp_master, "master.lp", "lp");
	//if (this->checkCPXstatus(this->cpx_env_master, status)) { check = false; goto TERMINATE; };

TERMINATE: 
	if (!check)
		return false;
	else
		return true;
}

bool RobustKSolver::add_primal_variables()
{
	int num_rows = CPXgetnumcols(this->cpx_env_master, this->cpx_lp_master);
	
	int status = 0;
	int vars_ind = 0;
	double lb[1], ub[1];
	double obj[1];
	char* varsname[1];
	char elname[1024];
	varsname[0] = elname;
	char vartype[1];

	/*
	Phi variable
	*/
	this->v_phi = num_rows;
	obj[0] = 1;
	vartype[0] = 'C';
	sprintf(varsname[0], "Phi");
	lb[0] = -CPX_INFBOUND;
	ub[0] = CPX_INFBOUND;
	status = CPXnewcols(this->cpx_env_master, this->cpx_lp_master, 1, obj, lb, ub, vartype, varsname);
	if (checkCPXstatus(this->cpx_env_master, status)) goto TERMINATE;
	++num_rows;

	/*
	First stage variables
	*/
	if (this->problem->getNx() > 0)
	{
		this->v_x.resize(problem->getNx());

		for (int i = 0; i < this->problem->getNx(); ++i)
		{
			this->v_x[i] = num_rows;

			obj[0] = problem->get1S_dobj_vect_c()[i];

			vartype[0] = 'B';
			
			sprintf(varsname[0], "x_i%d", i);
			lb[0] = 0;
			ub[0] = 1;



			status = CPXnewcols(this->cpx_env_master, this->cpx_lp_master, 1, obj, lb, ub, vartype, varsname);
			if (checkCPXstatus(this->cpx_env_master, status)) goto TERMINATE;
			++num_rows;
		}
	}
	/*
	Discovery variables
	*/
	if (this->problem->getNw() > 0)
	{
		this->v_w.resize(this->problem->getNw());

		for (int i = 0; i < this->problem->getNw(); ++i)
		{
			this->v_w[i] = num_rows;

			obj[0] = this->problem->getDV_dobj_vect_d()[i];
			
			vartype[0] = 'B';

			sprintf(varsname[0], "w_i%d", i);
			
			if (!ddid_version)
			{
				lb[0] = 1;
				ub[0] = 1;
			}
			else {
				lb[0] = 0;
				ub[0] = 1;
			}
			status = CPXnewcols(this->cpx_env_master, this->cpx_lp_master, 1, obj, lb, ub, vartype, varsname);
			if (checkCPXstatus(this->cpx_env_master, status)) goto TERMINATE;
			++num_rows;
		}
	}
	/*
	Second stage variables.
	*/
	if (this->problem->getNy() > 0)
	{
	

		this->v_y_k.resizeMatrix2D(this->problem->getK(), this->problem->getNy());

		for(int k = 0; k < this->problem->getK(); ++k)
		{
	

			for (int i = 0; i < this->problem->getNy(); ++i)
			{
				

				this->v_y_k(k, i) = num_rows;

				obj[0] = 0;

				vartype[0] = 'B';

				sprintf(varsname[0], "y_k%d_i%d", k,i);
				lb[0] = 0;
				ub[0] = 1;

	


				status = CPXnewcols(this->cpx_env_master, this->cpx_lp_master, 1, obj, lb, ub, vartype, varsname);
				if (checkCPXstatus(this->cpx_env_master, status)) goto TERMINATE;
				++num_rows;
			}
		}
	}

TERMINATE:
	
	if (status)
	{
		return false;
	}

	return true;
}

bool RobustKSolver::add_dualization_variables(unsigned int z)
{
	int status = 0;
	int vars_ind = 0;
	double lb[1], ub[1];
	double obj[1];
	char* varsname[1];
	char elname[1024];
	varsname[0] = elname;
	char vartype[1];

	double eps = 0.001;
	
	bool use_tighter_mccormick = this->cfg->getValueOfKey<bool>("USE_TIGHTER_MCCORMICK");


	int num_cols = CPXgetnumcols(this->cpx_env_master, this->cpx_lp_master);

	// Alpha variables: size = K.
	v_dual_set[z].v_alpha.resize(this->problem->getK());
	for (unsigned short k = 0; k < this->problem->getK(); ++k)
	{
		v_dual_set[z].v_alpha[k] = num_cols;

		obj[0] = 0;

		vartype[0] = 'C';

		sprintf(varsname[0], "alpha_z%d_k%d", z,k);


		//if (USE_STANDARD_FORMULATION)
		if(!use_tighter_mccormick)
		{
			lb[0] = 0;
			ub[0] = 1;
		}
		else {


			if (k == 0)
			{
				lb[0] = std::max(0.0, (1.00 / (double)this->problem->getK()) - eps);
				ub[0] = 1;
			}
			else
			{
				lb[0] = 0;
				ub[0] = std::min(1.00, (1.00 / (double)(k + 1)) + eps);
			}
		}

		status = CPXnewcols(this->cpx_env_master, this->cpx_lp_master, 1, obj, lb, ub, vartype, varsname);
		if (checkCPXstatus(this->cpx_env_master, status)) goto TERMINATE;
		++num_cols;
	}

	// Beta variables: size = Lxi {Number of constraints in the Uncertainty set}
	this->v_dual_set[z].v_beta.resize(this->problem->getLxi());
	for(unsigned short r = 0; r < this->problem->getLxi(); ++r)
	{
		this->v_dual_set[z].v_beta[r] = num_cols;

		obj[0] = 0;

		vartype[0] = 'C';

		sprintf(varsname[0], "beta_z%d_r%d", z, r);

		lb[0] = 0;
		ub[0] = CPX_INFBOUND;

		status = CPXnewcols(this->cpx_env_master, this->cpx_lp_master, 1, obj, lb, ub, vartype, varsname);
		if (checkCPXstatus(this->cpx_env_master, status)) goto TERMINATE;
		++num_cols;
	}
	// Beta^k variables. Each Beta^k has a size of Lxi.
	//this->v_dual_set[z].v_beta_k.resize(this->problem->getK());
	this->v_dual_set[z].v_beta_k.resizeMatrix2D(this->problem->getK(), this->problem->getLxi());

	for (unsigned short k = 0; k < this->problem->getK(); ++k)
	{
		// Resize the vector for each k
		//v_dual_set[z].v_beta_k[k].resize(this->problem->getLxi());

		for (unsigned short r = 0; r < this->problem->getLxi(); ++r)
		{
			//v_dual_set[z].v_beta_k[k][r] = num_cols;

			v_dual_set[z].v_beta_k(k, r) = num_cols;

			obj[0] = 0;

			vartype[0] = 'C';

			sprintf(varsname[0], "betaK_z%d_k%d_r%d", z, k,r);

			lb[0] = 0;
			ub[0] = CPX_INFBOUND;

			status = CPXnewcols(this->cpx_env_master, this->cpx_lp_master, 1, obj, lb, ub, vartype, varsname);
			if (checkCPXstatus(this->cpx_env_master, status)) goto TERMINATE;
			++num_cols;
		}
	}

	// Gamma^k variables
	//this->v_dual_set[z].v_gamma_k.resize(this->problem->getK());

	this->v_dual_set[z].v_gamma_k.resizeMatrix2D(this->problem->getK(), this->problem->getNxi());

	for (unsigned short k = 0; k < this->problem->getK(); ++k)
	{
		// Resize the vector for each k: equal to number of xi.
		//v_dual_set[z].v_gamma_k[k].resize(this->problem->getNxi());

		for (unsigned short r = 0; r < this->problem->getNxi(); ++r)
		{

			//v_dual_set[z].v_gamma_k[k][r] = num_cols;

			v_dual_set[z].v_gamma_k(k, r) = num_cols;

			obj[0] = 0;

			vartype[0] = 'C';

			sprintf(varsname[0], "t_gamma_z%d_k%d_r%d", z, k, r);

			lb[0] = - CPX_INFBOUND;
			ub[0] = CPX_INFBOUND;

			status = CPXnewcols(this->cpx_env_master, this->cpx_lp_master, 1, obj, lb, ub, vartype, varsname);
			if (checkCPXstatus(this->cpx_env_master, status)) goto TERMINATE;
			++num_cols;
		}
	}

TERMINATE:
	if (status)
	{
		return false;
	}
	return true;
}

bool RobustKSolver::add_bilinear_variables(unsigned int z)
{
	int num_cols = CPXgetnumcols(this->cpx_env_master, this->cpx_lp_master);

	int status = 0;
	int vars_ind = 0;
	double lb[1], ub[1];
	double obj[1];
	char* varsname[1];
	char elname[1024];
	varsname[0] = elname;
	char vartype[1];

	/*
	x times alpha
	*/
	if (this->problem->getNx() > 0 
		&&
		this->problem->get1S_uobj_matr_C().getNumElem() > 0
		)
	{
		this->v_dual_set[z].v_t_x_k.resizeMatrix2D(this->problem->getK(), this->problem->getNx());

		//this->v_dual_set[z].v_t_x_k.resize(this->problem->getK());

		for (int k = 0; k < this->problem->getK(); ++k)
		{
			//this->v_dual_set[z].v_t_x_k[k].resize(this->problem->getNx());

			for (int i = 0; i < this->problem->getNx(); ++i)
			{
				//this->v_dual_set[z].v_t_x_k[k][i] = num_cols;

				this->v_dual_set[z].v_t_x_k(k, i) = num_cols;
				obj[0] = 0;
				vartype[0] = 'C';

				sprintf(varsname[0], "t_x_z%d_k%d_i%d",z,k,i);
				lb[0] = 0;
				ub[0] = 1;

				status = CPXnewcols(this->cpx_env_master, this->cpx_lp_master, 1, obj, lb, ub, vartype, varsname);
				if (checkCPXstatus(this->cpx_env_master, status)) goto TERMINATE;
				++num_cols;
			}

		}
	}
	/*
	w times alpha
	*/
	if (this->problem->getNw() > 0
		&&
		this->problem->getDV_uobj_matr_D().getNumElem() > 0)
	{
		//this->v_dual_set[z].v_t_w_k.resize(this->problem->getK());

		this->v_dual_set[z].v_t_w_k.resizeMatrix2D(this->problem->getK(), this->problem->getNw());

		for (int k = 0; k < this->problem->getK(); ++k)
		{
			//this->v_dual_set[z].v_t_w_k[k].resize(this->problem->getNw());

			for (int i = 0; i < this->problem->getNw(); ++i)
			{
				//this->v_dual_set[z].v_t_w_k[k][i] = num_cols;

				this->v_dual_set[z].v_t_w_k(k, i) = num_cols;

				obj[0] = 0;

				vartype[0] = 'C';

				sprintf(varsname[0], "t_w_z%d_k%d_i%d", z,k, i);

				lb[0] = 0;
				ub[0] = 1;

				status = CPXnewcols(this->cpx_env_master, this->cpx_lp_master, 1, obj, lb, ub, vartype, varsname);
				if (checkCPXstatus(this->cpx_env_master, status)) goto TERMINATE;
				++num_cols;
			}
		}
	}
	/*
	y^k time alpha_k
	*/
	if (this->problem->getNy() > 0)
	{
		
		this->v_dual_set[z].v_t_y_k.resizeMatrix2D(this->problem->getK(), this->problem->getNy());

		for (int k = 0; k < this->problem->getK(); ++k)
		{
			for (int i = 0; i < this->problem->getNy(); ++i)
			{


				this->v_dual_set[z].v_t_y_k(k, i) = num_cols;

				obj[0] = 0;

				vartype[0] = 'C';

				sprintf(varsname[0], "t_y_z%d_k%d_i%d", z, k, i);

				lb[0] = 0;
				ub[0] = 1;

				status = CPXnewcols(this->cpx_env_master, this->cpx_lp_master, 1, obj, lb, ub, vartype, varsname);
				if (checkCPXstatus(this->cpx_env_master, status)) goto TERMINATE;
				++num_cols;
			}
		}
	}

	/*gamma_k times w*/

	this->v_dual_set[z].v_t_gamma_k.resizeMatrix2D(this->problem->getK(), this->problem->getNw());

	for (unsigned short k = 0; k < this->problem->getK(); ++k)
	{
	

		for (unsigned short r = 0; r < this->problem->getNw(); ++r)
		{

			v_dual_set[z].v_t_gamma_k(k, r) = num_cols;

			obj[0] = 0;

			vartype[0] = 'C';

			sprintf(varsname[0], "t_gammaK_z%d_k%d_r%d", z, k, r);

			lb[0] = -CPX_INFBOUND;
			ub[0] = CPX_INFBOUND;

			status = CPXnewcols(this->cpx_env_master, this->cpx_lp_master, 1, obj, lb, ub, vartype, varsname);
			if (checkCPXstatus(this->cpx_env_master, status)) goto TERMINATE;
			++num_cols;
		}
	}
	


TERMINATE:

	if (status)
	{
		return false;
	}

	return true;
}

bool RobustKSolver::add_deterministic_constraints()
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
	status = 0;
	double coeff_eps = COEFF_EPS;

	std::vector<double> b_loc = this->problem->get1S_drhs_vect_bx();
	Matrix2D<double> Matr_loc = this->problem->get1S_dcon_matr_X();



	// X matrix constraints
	for (unsigned short row = 0; row < this->problem->getLx(); ++row)
	{
		rhs[0] = b_loc[row];
		sense[0] = 'L';
		nzc = 0;
		sprintf(cnstrname[0], "Xx_%d", row);
		// Now we go through the "columns", coefficients of the matrices 

		for (unsigned short cols = 0; cols < this->problem->getNx(); ++cols)
			if (fabs(Matr_loc(row,cols)) > coeff_eps)
			{
				this->matind[nzc] = this->v_x[cols];
				this->matval[nzc] = Matr_loc(row,cols);
				++nzc;
			}
		status = CPXaddrows(this->cpx_env_master, this->cpx_lp_master, 0, 1, nzc, rhs, sense, matbeg, this->matind.data(), this->matval.data(), NULL, cnstrname);
		if (checkCPXstatus(this->cpx_env_master, status)) goto TERMINATE;

	}

	// W matrices constraints 
	b_loc = this->problem->getDV_drhs_vect_bw();
	Matr_loc = this->problem->getDV_dcon_matr_W();

	for (unsigned short row = 0; row < this->problem->getLw(); ++row)
	{
		rhs[0] = b_loc[row];
		sense[0] = 'L';
		nzc = 0;
		sprintf(cnstrname[0], "Ww_%d", row);
		// Now we go trought the "columns", coefficient of the matrices 

		for (unsigned short cols = 0; cols < this->problem->getNw(); ++cols)
			if (fabs(Matr_loc(row, cols)) > coeff_eps)
			{
				this->matind[nzc] = this->v_w[cols];
				this->matval[nzc] = Matr_loc(row,cols);
				++nzc;
			}
		status = CPXaddrows(this->cpx_env_master, this->cpx_lp_master, 0, 1, nzc, rhs, sense, matbeg, this->matind.data(), this->matval.data(), NULL, cnstrname);
		if (checkCPXstatus(this->cpx_env_master, status)) goto TERMINATE;

	}


	// Y matrices constraints for each k

	b_loc = this->problem->get2S_drhs_vect_by();
	Matr_loc = this->problem->get2S_dcon_matr_Y();
	
	for (int k = 0; k < this->problem->getK(); ++k)
	{
		for (unsigned short row = 0; row < this->problem->getLy(); ++row)
		{
			rhs[0] = b_loc[row];
			sense[0] = 'L';
			nzc = 0;
			sprintf(cnstrname[0], "Yy_k%d_%d", k,row);
			// Now we go trought the "columns", coefficient of the matrices 

			for (unsigned short cols = 0; cols < this->problem->getNy(); ++cols)
				if (fabs(Matr_loc(row,cols)) > coeff_eps)
				{
					this->matind[nzc] = this->v_y_k(k,cols);
					this->matval[nzc] = Matr_loc(row, cols);
					++nzc;
				}
			status = CPXaddrows(this->cpx_env_master, this->cpx_lp_master, 0, 1, nzc, rhs, sense, matbeg, this->matind.data(), this->matval.data(), NULL, cnstrname);
			if (checkCPXstatus(this->cpx_env_master, status)) goto TERMINATE;

		}
	}
	// Ex(x) + Ew(w) + Ey(y) <=g
	if (this->problem->getR() > 0)
	{
		b_loc = this->problem->getLink_dcon_vect_g();
		

		for (int row = 0; row < this->problem->getR(); ++row)
			for (int k = 0; k < this->problem->getK(); ++k)
			{
		
			sense[0] = 'L';
			nzc = 0;
			rhs[0] = b_loc[row];

			sprintf(cnstrname[0], "LinkXWY_k%d_r%d", k, row);
			
			Matr_loc = this->problem->getLink_dcon_matr_Ex();
			if (Matr_loc.getNumElem() > 0)
			{
				for (int col = 0; col < this->problem->getNx(); ++col)
					if (fabs(Matr_loc(row, col)) > coeff_eps)
					{
						this->matind[nzc] = this->v_x[col];
						this->matval[nzc] = Matr_loc(row, col);
						++nzc;
					}
			}
			Matr_loc = this->problem->getLink_dcon_matr_Ew();
			if (Matr_loc.getNumElem() > 0)
			{
				for (int col = 0; col < this->problem->getNw(); ++col)
					if (fabs(Matr_loc(row, col)) > coeff_eps)
					{
						this->matind[nzc] = this->v_w[col];
						this->matval[nzc] = Matr_loc(row, col);
						++nzc;
					}
			}

			Matr_loc = this->problem->getLink_dcon_matr_Ey();
			if (Matr_loc.getNumElem() > 0)
			{
				for (int col = 0; col < this->problem->getNy(); ++col)
					if (fabs(Matr_loc(row, col)) > coeff_eps)
					{
						this->matind[nzc] = this->v_y_k(k,col);
						this->matval[nzc] = Matr_loc(row, col);
						++nzc;
					}
			}

			status = CPXaddrows(this->cpx_env_master, this->cpx_lp_master, 0, 1, nzc, rhs, sense, matbeg, this->matind.data(), this->matval.data(), NULL, cnstrname);
			if (checkCPXstatus(this->cpx_env_master, status)) goto TERMINATE;
			}
	}


TERMINATE:
	if (status)
	{
		return false;
	}

	return true;
}

bool RobustKSolver::add_dualization_constraints(unsigned int z)
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

	double coeff_eps = COEFF_EPS;
	double int_eps = INT_EPS;
	double big_m = CNSTR_UNC_BIG_M;

	Matrix2D<double> A_loc = this->problem->getUC_xicon_matr_A();
	Matrix2D<double> Q_loc = this->problem->get2S_uobj_matr_Q();
	Matrix2D<double> C_loc = this->problem->get1S_uobj_matr_C();
	Matrix2D<double> D_loc = this->problem->getDV_uobj_matr_D();

	Matrix2D<double> M_loc;
	Matrix3D<double> M3_loc;

	std::vector<double> b_loc;
	double coeff = 0;

	/*\phi >= uncertainty_cost*/
	sense[0] = 'G';
	rhs[0] = 0;
	nzc = 0;
	sprintf(cnstrname[0], "objective_constr_z%d", z);

	// \phi
	this->matind[nzc] = this->v_phi;
	this->matval[nzc] = 1;
	++nzc;
	// \beta and \beta_k
	for (int bb = 0; bb < this->problem->getLxi(); ++bb)
	{
		this->matind[nzc] = this->v_dual_set[z].v_beta[bb];
		this->matval[nzc] = -this->problem->getUC_rhs_vect_b()[bb];
		++nzc;

		if (nzc >= MAX_NUM_VARS)
		{
			getchar();
		}

		for (int k = 0; k < this->problem->getK(); ++k)
		{
			this->matind[nzc] = this->v_dual_set[z].v_beta_k(k, bb);
			this->matval[nzc] = -this->problem->getUC_rhs_vect_b()[bb];
			++nzc;

			if (nzc >= MAX_NUM_VARS)
			{
				getchar();
			}
		}
	}
	// alpha in objective only if CU.
	if (this->problem->getL() > 0)
	{
		b_loc = this->problem->get_rhs_ucon_vect_h();

		for (int k = 0; k < this->problem->getK(); ++k)
		{
			coeff = 0;

			for (int l = 0; l < this->problem->getL(); ++l)
				if (this->v_dual_set[z].z_k(k, l) > 0
					&&
					fabs(b_loc[l]) > coeff_eps
					)
				{
					coeff += b_loc[l];
				}

			if (fabs(coeff) > coeff_eps)
			{
				this->matind[nzc] = this->v_dual_set[z].v_alpha[k];
				this->matval[nzc] = (coeff * big_m);
				//this->matval[nzc] = -(coeff * big_m);
				++nzc;
			}

		}
	}

	// tilde_x
	M_loc = this->problem->get1S_ucon_matr_T();
	for (int k = 0; k < this->problem->getK(); ++k)
		for (int xi = 0; xi < this->problem->getNx(); ++xi)
		{

			//this->matind[nzc] = this->v_x[xi];
			coeff = 0;
			/*Tx part*/

			for (int l = 0; l < this->problem->getL(); ++l)
				if (this->v_dual_set[z].z_k(k, l) > 0
					&&
					fabs(M_loc(l, xi)) > coeff_eps)
				{
					coeff -= M_loc(l, xi);
					//coeff += M_loc(l, xi);
					//this->matval[nzc] += 
				}


			if (fabs(coeff) > coeff_eps)
			{
				this->matind[nzc] = this->v_dual_set[z].v_t_x_k(k, xi);
				this->matval[nzc] = coeff * big_m;
				++nzc;
			}


			if (nzc >= MAX_NUM_VARS)
			{
				getchar();
			}
		}

	//tilde_w
	
		M_loc = this->problem->getDV_ucon_matr_V();
		for (int k = 0; k < this->problem->getK(); ++k)
			for (int wi = 0; wi < this->problem->getNw(); ++wi)
			{

				//this->matind[nzc] = this->v_x[xi];
				coeff = 0;
				/*Vw part*/

				for (int l = 0; l < this->problem->getL(); ++l)
					if (this->v_dual_set[z].z_k(k, l) > 0
						&&
						fabs(M_loc(l, wi)) > coeff_eps)
					{
						coeff -= M_loc(l, wi);
						//coeff += M_loc(l, wi);
						//this->matval[nzc] += M_loc(l, wi);
					}

				if (fabs(coeff) > coeff_eps)
				{
					this->matind[nzc] = this->v_dual_set[z].v_t_w_k(k, wi);
					this->matval[nzc] = coeff * big_m;
					++nzc;
				}

				if (nzc >= MAX_NUM_VARS)
				{
					getchar();
				}
			}
	

	// tilde_y
	
		M_loc = this->problem->get2S_ucon_matr_P();
	for (int k = 0; k < this->problem->getK(); ++k)
	{
		for (int yi = 0; yi < this->problem->getNy(); ++yi)
		{

			//this->matind[nzc] = this->v_dual_set[z].v_t_y_k(k,yi);
			coeff = -this->problem->get2S_dobj_vect_q()[yi];
			//this->matval[nzc] = 
			//coeff = 0;

			/*Py part*/
			for (int l = 0; l < this->problem->getL(); ++l)
				if (this->v_dual_set[z].z_k(k, l) > 0
					&&
					fabs(M_loc(l, yi)) > coeff_eps)
				{
					coeff -= (M_loc(l, yi) * big_m);
					//coeff += M_loc(l, yi) * big_m;
					//this->matval[nzc] += M_loc(l, yi) * big_m;
				}


			if (fabs(coeff) > coeff_eps)
			{
				this->matind[nzc] = this->v_dual_set[z].v_t_y_k(k, yi);
				this->matval[nzc] = coeff;
				++nzc;
			}

			if (nzc >= MAX_NUM_VARS)
			{
				getchar();
			}
		}

	}


	if (nzc > 0)
	{
		status = CPXaddrows(this->cpx_env_master, this->cpx_lp_master, 0, 1, nzc, rhs, sense, matbeg, this->matind.data(), this->matval.data(), NULL, cnstrname);
		if (checkCPXstatus(this->cpx_env_master, status)) goto TERMINATE;
	}


	/*Sum over alpha = 1*/
	sense[0] = 'E';
	rhs[0] = 1;
	nzc = 0;

	sprintf(cnstrname[0], "Sum_k_Alpha_%d = 1",z);

	for (int k = 0; k < this->problem->getK(); ++k)
	{
		this->matind[nzc] = this->v_dual_set[z].v_alpha[k];
		this->matval[nzc] = 1;
		++nzc;

	}
	status = CPXaddrows(this->cpx_env_master, this->cpx_lp_master, 0, 1, nzc, rhs, sense, matbeg, this->matind.data(), this->matval.data(), NULL, cnstrname);
	if (checkCPXstatus(this->cpx_env_master, status)) goto TERMINATE;
	

	/*At Beta + \sum_k w \circ gamma_k* , \foraeach Nxi (one constr for each component of the \xi vector  */
	/*
	k is the policy index.
	row is the constraints end (equivalenlty) the component of vector \xi
	*/

	for (int row = 0; row < this->problem->getNxi(); ++row)
	{
		rhs[0] = 0;
		sense[0] = 'E';
		sprintf(cnstrname[0], "Balance_z%d_n%d", z, row);
		nzc = 0;

		if (A_loc.getNumElem() > 0)
		{
			// A is transpose. THerefore row and col are inverted.
			for (int col = 0; col < this->problem->getLxi(); ++col)
				if (fabs(
					A_loc(col, row)
					//this->problem->getUC_xicon_matr_A()(col,row)
				)
			> coeff_eps)
				{
					this->matind[nzc] = this->v_dual_set[z].v_beta[col];
					this->matval[nzc] = A_loc(col, row);
					++nzc;
				}
		}

		for (int k = 0; k < this->problem->getK(); ++k)
		{
			this->matind[nzc] = this->v_dual_set[z].v_t_gamma_k(k, row);
			this->matval[nzc] = 1;
			++nzc;
		}

		
		status = CPXaddrows(this->cpx_env_master, this->cpx_lp_master, 0, 1, nzc, rhs, sense, matbeg, this->matind.data(), this->matval.data(), NULL, cnstrname);
		if (checkCPXstatus(this->cpx_env_master, status)) goto TERMINATE;
	
	}


	
	for (int k = 0; k < this->problem->getK(); ++k)
		for (int row = 0; row < this->problem->getNxi(); ++row)
		{
			rhs[0] = 0;
			sense[0] = 'E';
			sprintf(cnstrname[0], "BalanceK_z%d_k%d_r%d", z, k, row);
			nzc = 0;

			/*
			At Beta^k
			*/
			if (A_loc.getNumElem() > 0)
			{
				// A is transpose. Therefore row and col are inverted.

				for (int col = 0; col < this->problem->getLxi(); ++col)
					if (
						fabs(
							//this->problem->getUC_xicon_matr_A()(col,row)
							A_loc(col, row)
						) > coeff_eps)
					{
						this->matind[nzc] = this->v_dual_set[z].v_beta_k(k, col);
						this->matval[nzc] = A_loc(col, row);
						++nzc;
					}
			}

			/*
			 -(w \circ gamma)
			*/
			this->matind[nzc] = this->v_dual_set[z].v_t_gamma_k(k, row);
			this->matval[nzc] = -1;
			++nzc;


			//Cx \times alpha_k + T_n x \alpha 
			
				M3_loc = this->problem->get1S_ucon_matr_Txi();
			for (int col = 0; col < this->problem->getNx(); ++col)
			{
				coeff = 0;

				if (fabs(C_loc(row, col)) > coeff_eps)
				{

					coeff = -C_loc(row, col);
					//++nzc;

				}

				for (int l = 0; l < this->problem->getL(); ++l)
					if (this->v_dual_set[z].z_k(k, l) > 0
						&&
						fabs(M3_loc(row, l, col)) > coeff_eps)
					{
						coeff += -(M3_loc(row, l, col) * big_m);
					}

				if (fabs(coeff) > coeff_eps)
				{
					this->matind[nzc] = this->v_dual_set[z].v_t_x_k(k, col);
					this->matval[nzc] = coeff;
					++nzc;

				}
			}
		

			//Dw\times alpha_k + V_n w \alpha 
			if (this->problem->getDV_ucon_matr_Vxi().getNumElem() > 0 || D_loc.getNumElem() > 0) 
			{
				M3_loc = this->problem->getDV_ucon_matr_Vxi();
				for (int col = 0; col < this->problem->getNw(); ++col)
				{
					coeff = 0;

					if (fabs(D_loc(row, col)) > coeff_eps)
					{

						coeff = -D_loc(row, col);
						//++nzc;

					}
					for (int l = 0; l < this->problem->getL(); ++l)
						if (this->v_dual_set[z].z_k(k, l) > coeff_eps
							&&
							fabs(M3_loc(row, l, col)) > coeff_eps)
						{
							coeff += -(M3_loc(row, l, col) * big_m);
						}

					if (fabs(coeff) > coeff_eps)
					{
						this->matind[nzc] = this->v_dual_set[z].v_t_w_k(k, col);
						this->matval[nzc] = coeff;
						++nzc;
					}
				}
			}

			/*if (this->problem->get2S_ucon_matr_Pxi().getNumElem() > 0)
			{*/
				//Qy\times alpha_k + P_n w \alpha 
				M3_loc = this->problem->get2S_ucon_matr_Pxi();
				for (int col = 0; col < this->problem->getNy(); ++col)
				{
					coeff = 0;

					if (fabs(Q_loc(row, col)) > coeff_eps)
					{

						coeff = -Q_loc(row, col);

					}

					for (int l = 0; l < this->problem->getL(); ++l)
						if (this->v_dual_set[z].z_k(k, l) > 0
							&&
							fabs(M3_loc(row, l, col)) > coeff_eps)
						{
							coeff += -(M3_loc(row, l, col) * big_m);
						}

					if (fabs(coeff) > coeff_eps)
					{
						this->matind[nzc] = this->v_dual_set[z].v_t_y_k(k, col);
						this->matval[nzc] = coeff;
						++nzc;
					}
				}
			//}

			/*alpha_k times H */
				if (this->problem->getUC_ucon_matr_H().getNumElem() > 0)
				{
					M_loc = this->problem->getUC_ucon_matr_H();
					coeff = 0;
					for (int l = 0; l < this->problem->getL(); ++l)
						if (this->v_dual_set[z].z_k(k, l) > 0
							&&
							fabs(M_loc(l, row)) > coeff_eps)
						{
							coeff += M_loc(l, row) * big_m;
						}
					if (fabs(coeff) > coeff_eps)
					{
						this->matind[nzc] = this->v_dual_set[z].v_alpha[k];
						this->matval[nzc] = coeff;
						++nzc;
					}



					////Dw \times alpha_k
					//if (D_loc.getNumElem() > 0)
					//{
					//	for (int col = 0; col < this->problem->getNw(); ++col)
					//		if(fabs(this->problem->getDV_uobj_matr_D()(row,col)) > coeff_eps)
					//	{
					//		this->matind[nzc] = this->v_dual_set[z].v_t_w_k(k,col);
					//		this->matval[nzc] = - this->problem->getDV_uobj_matr_D()(row,col);
					//		++nzc;
					//	}

					//}

					////Qy^k \times alpha_k
					//if (Q_loc.getNumElem() > 0)
					//{
					//	for (int col = 0; col < this->problem->getNy(); ++col)
					//		if (fabs(Q_loc(row,col)) > coeff_eps)
					//	{
					//		this->matind[nzc] = this->v_dual_set[z].v_t_y_k(k,col);
					//		this->matval[nzc] = - Q_loc(row,col);
					//		++nzc;
					//	}
					//}
					/*if (z != 0)
					{*/
				}
			status = CPXaddrows(this->cpx_env_master, this->cpx_lp_master, 0, 1, nzc, rhs, sense, matbeg, this->matind.data(), this->matval.data(), NULL, cnstrname);
			if (checkCPXstatus(this->cpx_env_master, status)) goto TERMINATE;
			//}
		
		}

	/*
	* Constraints for \zeta and \zeta^k
	*/
	if (this->problem->getNzeta() > 0 &&
		this->problem->getUC_zetacon_matr_G().getNumElem() > 0
		)
	{

		/*
		\zeta
		*/
		rhs[0] = 0;
		sense[0] = 'E';
		nzc = 0;

		for (int row = 0; row < this->problem->getNzeta(); ++row)
		{
			
			sprintf(cnstrname[0], "BalanceZeta_r%d", row);
			nzc = 0;

			for (int col = 0; col < this->problem->getNzeta(); ++col)
				if(fabs(this->problem->getUC_zetacon_matr_G()(col,row)) > coeff_eps)
			{
					this->matind[nzc] = this->v_dual_set[z].v_beta[col];
					this->matval[nzc] = this->problem->getUC_zetacon_matr_G()(col,row);
					++nzc;
			}
			status = CPXaddrows(this->cpx_env_master, this->cpx_lp_master, 0, 1, nzc, rhs, sense, matbeg, this->matind.data(), this->matval.data(), NULL, cnstrname);
			if (checkCPXstatus(this->cpx_env_master, status)) goto TERMINATE;
		}
		/*
		\zeta^k
		*/
		rhs[0] = 0;
		sense[0] = 'E';
		nzc = 0;

		for(int k = 0; k < this->problem->getK(); ++k)
			for (int row = 0; row < this->problem->getNzeta(); ++row)
			{

			sprintf(cnstrname[0], "BalanceZetaK_k%d_r%d", k,row);
			nzc = 0;

			for (int col = 0; col < this->problem->getNzeta(); ++col)
				if (fabs(this->problem->getUC_zetacon_matr_G()(col,row)) > coeff_eps)
				{
					this->matind[nzc] = this->v_dual_set[z].v_beta_k(k,row);
					this->matval[nzc] = this->problem->getUC_zetacon_matr_G()(col,row);
					++nzc;
				}
			status = CPXaddrows(this->cpx_env_master, this->cpx_lp_master, 0, 1, nzc, rhs, sense, matbeg, this->matind.data(), this->matval.data(), NULL, cnstrname);
			if (checkCPXstatus(this->cpx_env_master, status)) goto TERMINATE;
			}


	}


TERMINATE:
	if (status)
	{
		return false;
	}

	return true;
}

bool RobustKSolver::add_linearization_constraints(unsigned int z)
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

	double coeff_eps = COEFF_EPS;
	double big_m = BIG_M;
	status = 0;
	/*
	 x * alpha_k
	*/
	if (this->problem->getNx() > 0
		&&
		this->problem->get1S_uobj_matr_C().getNumElem() > 0)
	{
		for (int k = 0; k < this->problem->getK(); ++k)
			for (int xi = 0; xi < this->problem->getNx(); ++xi)
			{
				/*
					* t_x_k <= x
					*/
				rhs[0] = 0;
				sense[0] = 'L';
				sprintf(cnstrname[0], "Linear_t_x_(1)_k%d_i%d", k, xi);
				nzc = 0;

				this->matind[nzc] = this->v_dual_set[z].v_t_x_k(k,xi);
				this->matval[nzc] = 1;
				++nzc;

				this->matind[nzc] = this->v_x[xi];
				this->matval[nzc] = -1;
				++nzc;

				status = CPXaddrows(this->cpx_env_master, this->cpx_lp_master, 0, 1, nzc, rhs, sense, matbeg, this->matind.data(), this->matval.data(), NULL, cnstrname);
				if (checkCPXstatus(this->cpx_env_master, status)) goto TERMINATE;

				/*
				* t_x_k <= alpha_k
				*/
				rhs[0] = 0;
				sense[0] = 'L';
				sprintf(cnstrname[0], "Linear_t_x_(2)_k%d_i%d", k, xi);
				nzc = 0;

				this->matind[nzc] = this->v_dual_set[z].v_t_x_k(k,xi);
				this->matval[nzc] = 1;
				++nzc;

				this->matind[nzc] = this->v_dual_set[z].v_alpha[k];
				this->matval[nzc] = -1;
				++nzc;

				status = CPXaddrows(this->cpx_env_master, this->cpx_lp_master, 0, 1, nzc, rhs, sense, matbeg, this->matind.data(), this->matval.data(), NULL, cnstrname);
				if (checkCPXstatus(this->cpx_env_master, status)) goto TERMINATE;

				/*
				t_x_k - alpha_k - x >= -1
				*/
				rhs[0] = -1;
				sense[0] = 'G';
				sprintf(cnstrname[0], "Linear_t_x_(3)_k%d_i%d", k, xi);
				nzc = 0;

				this->matind[nzc] = this->v_dual_set[z].v_t_x_k(k,xi);
				this->matval[nzc] = 1;
				++nzc;

				this->matind[nzc] = this->v_dual_set[z].v_alpha[k];
				this->matval[nzc] = -1;
				++nzc;

				this->matind[nzc] = this->v_x[xi];
				this->matval[nzc] = -1;
				++nzc;

				status = CPXaddrows(this->cpx_env_master, this->cpx_lp_master, 0, 1, nzc, rhs, sense, matbeg, this->matind.data(), this->matval.data(), NULL, cnstrname);
				if (checkCPXstatus(this->cpx_env_master, status)) goto TERMINATE;
			}
	}

	/*
	* w * alpha_k
	* 
	*/
	if (this->problem->getNw() > 0
		&&
		this->problem->getDV_uobj_matr_D().getNumElem() > 0)
	{
		for (int k = 0; k < this->problem->getK(); ++k)
			for (int wi = 0; wi < this->problem->getNw(); ++wi)
			{
				/*
					* t_w_k <= w
					*/
				rhs[0] = 0;
				sense[0] = 'L';
				sprintf(cnstrname[0], "Linear_t_w_(1)_k%d_i%d", k, wi);
				nzc = 0;

				this->matind[nzc] = this->v_dual_set[z].v_t_w_k(k,wi);
				this->matval[nzc] = 1;
				++nzc;

				this->matind[nzc] = this->v_w[wi];
				this->matval[nzc] = -1;
				++nzc;

				status = CPXaddrows(this->cpx_env_master, this->cpx_lp_master, 0, 1, nzc, rhs, sense, matbeg, this->matind.data(), this->matval.data(), NULL, cnstrname);
				if (checkCPXstatus(this->cpx_env_master, status)) goto TERMINATE;

				/*
				* t_w_k <= alpha_k
				*/
				rhs[0] = 0;
				sense[0] = 'L';
				sprintf(cnstrname[0], "Linear_t_w_(2)_k%d_i%d", k, wi);
				nzc = 0;

				this->matind[nzc] = this->v_dual_set[z].v_t_w_k(k,wi);
				this->matval[nzc] = 1;
				++nzc;

				this->matind[nzc] = this->v_dual_set[z].v_alpha[k];
				this->matval[nzc] = -1;
				++nzc;

				status = CPXaddrows(this->cpx_env_master, this->cpx_lp_master, 0, 1, nzc, rhs, sense, matbeg, this->matind.data(), this->matval.data(), NULL, cnstrname);
				if (checkCPXstatus(this->cpx_env_master, status)) goto TERMINATE;

				/*
				t_w_k - alpha_k - w(UB alpha) >= -1(UN alpha)
				*/
				rhs[0] = -1;
				sense[0] = 'G';
				sprintf(cnstrname[0], "Linear_t_x_(3)_k%d_i%d", k, wi);
				nzc = 0;

				this->matind[nzc] = this->v_dual_set[z].v_t_w_k(k,wi);
				this->matval[nzc] = 1;
				++nzc;

				this->matind[nzc] = this->v_dual_set[z].v_alpha[k];
				this->matval[nzc] = -1;
				++nzc;

				this->matind[nzc] = this->v_w[wi];
				this->matval[nzc] = -1;
				++nzc;

				status = CPXaddrows(this->cpx_env_master, this->cpx_lp_master, 0, 1, nzc, rhs, sense, matbeg, this->matind.data(), this->matval.data(), NULL, cnstrname);
				if (checkCPXstatus(this->cpx_env_master, status)) goto TERMINATE;
			}
	}



	/*
	* y^k * alpha_k
	*/
	for(int k = 0; k < this->problem->getK(); ++k)
		for (int yi = 0; yi < this->problem->getNy(); ++yi)
		{
			/*
			* t_y_k <= y_k
			*/
			rhs[0] = 0;
			sense[0] = 'L';
			sprintf(cnstrname[0], "Linear_t_y_(1)_k%d_i%d", k, yi);
			nzc = 0;
			
			this->matind[nzc] = this->v_dual_set[z].v_t_y_k(k,yi);
			this->matval[nzc] = 1;
			++nzc;

			this->matind[nzc] = this->v_y_k(k,yi);
			this->matval[nzc] = -1;
			++nzc;

			status = CPXaddrows(this->cpx_env_master, this->cpx_lp_master, 0, 1, nzc, rhs, sense, matbeg, this->matind.data(), this->matval.data(), NULL, cnstrname);
			if (checkCPXstatus(this->cpx_env_master, status)) goto TERMINATE;

			/*
			* t_y_k <= alpha_k
			*/
			rhs[0] = 0;
			sense[0] = 'L';
			sprintf(cnstrname[0], "Linear_t_y_(2)_k%d_i%d", k, yi);
			nzc = 0;

			this->matind[nzc] = this->v_dual_set[z].v_t_y_k(k,yi);
			this->matval[nzc] = 1;
			++nzc;

			this->matind[nzc] = this->v_dual_set[z].v_alpha[k];
			this->matval[nzc] = -1;
			++nzc;

			status = CPXaddrows(this->cpx_env_master, this->cpx_lp_master, 0, 1, nzc, rhs, sense, matbeg, this->matind.data(), this->matval.data(), NULL, cnstrname);
			if (checkCPXstatus(this->cpx_env_master, status)) goto TERMINATE;

			/*
			t_y_k - alpha_k - y_k(UB alpha) >= -1(UB alpha)
			*/
			rhs[0] = -1;
			sense[0] = 'G';
			sprintf(cnstrname[0], "Linear_t_y_(3)_k%d_i%d", k, yi);
			nzc = 0;

			this->matind[nzc] = this->v_dual_set[z].v_t_y_k(k,yi);
			this->matval[nzc] = 1;
			++nzc;

			this->matind[nzc] = this->v_dual_set[z].v_alpha[k];
			this->matval[nzc] = -1;
			++nzc;

			this->matind[nzc] = this->v_y_k(k,yi);
			this->matval[nzc] = -1;
			++nzc;

			status = CPXaddrows(this->cpx_env_master, this->cpx_lp_master, 0, 1, nzc, rhs, sense, matbeg, this->matind.data(), this->matval.data(), NULL, cnstrname);
			if (checkCPXstatus(this->cpx_env_master, status)) goto TERMINATE;
		}

	/*
	w * gamma_k
	*/
	if (this->problem->getNw() > 0)
	{


		for (int k = 0; k < this->problem->getK(); ++k)
			for (int wi = 0; wi < this->problem->getNw(); ++wi)
			{
				/*
				* t_gamma_k <= w M
				*/
				rhs[0] = 0;
				sense[0] = 'L';
				sprintf(cnstrname[0], "Linear_t_gamma_(1)_k%d_i%d", k, wi);
				nzc = 0;

				this->matind[nzc] = this->v_dual_set[z].v_t_gamma_k(k,wi);
				this->matval[nzc] = 1;
				++nzc;

				this->matind[nzc] = this->v_w[wi];
				this->matval[nzc] = -big_m;
				++nzc;

				status = CPXaddrows(this->cpx_env_master, this->cpx_lp_master, 0, 1, nzc, rhs, sense, matbeg, this->matind.data(), this->matval.data(), NULL, cnstrname);
				if (checkCPXstatus(this->cpx_env_master, status)) goto TERMINATE;

				/*
				* t_gamma_k >= -w M
				*/
				rhs[0] = 0;
				sense[0] = 'G';
				sprintf(cnstrname[0], "Linear_t_gamma_(1)_k%d_i%d", k, wi);
				nzc = 0;

				this->matind[nzc] = this->v_dual_set[z].v_t_gamma_k(k,wi);
				this->matval[nzc] = 1;
				++nzc;

				this->matind[nzc] = this->v_w[wi];
				this->matval[nzc] = big_m;
				++nzc;

				status = CPXaddrows(this->cpx_env_master, this->cpx_lp_master, 0, 1, nzc, rhs, sense, matbeg, this->matind.data(), this->matval.data(), NULL, cnstrname);
				if (checkCPXstatus(this->cpx_env_master, status)) goto TERMINATE;


				/*
				t_gamma_k - gamma_k + Mw <= M
				*/
				rhs[0] = big_m;
				sense[0] = 'L';
				sprintf(cnstrname[0], "Linear_t_gamma_(3)_k%d_i%d", k, wi);
				nzc = 0;

				this->matind[nzc] = this->v_dual_set[z].v_t_gamma_k(k,wi);
				this->matval[nzc] = 1;
				++nzc;

				this->matind[nzc] = this->v_dual_set[z].v_gamma_k(k,wi);
				this->matval[nzc] = -1;
				++nzc;

				this->matind[nzc] = this->v_w[wi];
				this->matval[nzc] = big_m;
				++nzc;

				status = CPXaddrows(this->cpx_env_master, this->cpx_lp_master, 0, 1, nzc, rhs, sense, matbeg, this->matind.data(), this->matval.data(), NULL, cnstrname);
				if (checkCPXstatus(this->cpx_env_master, status)) goto TERMINATE;




				/*
				t_gamma_k - gamma_k - Mw >= -M
				*/
				rhs[0] = -big_m;
				sense[0] = 'G';
				sprintf(cnstrname[0], "Linear_t_gamma_(3)_k%d_i%d", k, wi);
				nzc = 0;

				this->matind[nzc] = this->v_dual_set[z].v_t_gamma_k(k,wi);
				this->matval[nzc] = 1;
				++nzc;

				this->matind[nzc] = this->v_dual_set[z].v_gamma_k(k,wi);
				this->matval[nzc] = -1;
				++nzc;

				this->matind[nzc] = this->v_w[wi];
				this->matval[nzc] = -big_m;
				++nzc;

				status = CPXaddrows(this->cpx_env_master, this->cpx_lp_master, 0, 1, nzc, rhs, sense, matbeg, this->matind.data(), this->matval.data(), NULL, cnstrname);
				if (checkCPXstatus(this->cpx_env_master, status)) goto TERMINATE;


			}
	}
TERMINATE:
	if (status)
	{
		return false;
	}

	return true;
}

bool RobustKSolver::add_primal_valid_inequalities()
{
	/*
	For ^k
	*/
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
	status = 0;
	double coeff_eps = COEFF_EPS;

	std::vector<double> y_coeff_loc = this->problem->get2S_dobj_vect_lbq();
	std::vector<double> x_coeff_loc = this->problem->get1S_dobj_vect_lbc();
	std::vector<double> c_loc = this->problem->get1S_dobj_vect_c();


	
		for (int k = 0; k < this->problem->getK(); ++k)
		{
			rhs[0] = 0;
			sense[0] = 'G';
			nzc = 0;

			sprintf(cnstrname[0], "Det_Obj_cut_k%d", k);

			this->matind[nzc] = this->v_phi;
			this->matval[nzc] = 1;
			++nzc;


			for (unsigned short cols = 0; cols < this->problem->getNx(); ++cols)
				if (fabs(x_coeff_loc[cols]) > coeff_eps)
				{
					this->matind[nzc] = this->v_x[cols];
					this->matval[nzc] = -x_coeff_loc[cols];
					++nzc;
				}



			for (unsigned short cols = 0; cols < this->problem->getNy(); ++cols)
				if (fabs(y_coeff_loc[cols]) > coeff_eps)
				{
					this->matind[nzc] = this->v_y_k(k, cols);
					this->matval[nzc] = -y_coeff_loc[cols];
					++nzc;
				}

			status = CPXaddrows(this->cpx_env_master, this->cpx_lp_master, 0, 1, nzc, rhs, sense, matbeg, this->matind.data(), this->matval.data(), NULL, cnstrname);
			if (checkCPXstatus(this->cpx_env_master, status)) goto TERMINATE;
		}


		/*
		Here we add the approx constr cut.
		*/

		/*This is only for constraint uncertainty.*/
		y_coeff_loc = this->problem->get2S_ducon_vect_p();
		x_coeff_loc = this->problem->get1S_ducon_vect_t();

		for (int k = 0; k < this->problem->getK(); ++k)
			for (int l = 0; l < this->problem->getL(); ++l)
			{
				rhs[0] = this->problem->get_rhs_ducon_vect_hmax()[l];
				sense[0] = 'L';
				nzc = 0;

				sprintf(cnstrname[0], "Det_CU_cut_k%d_l%d", k, l);


				for (unsigned short col = 0; col < this->problem->getNx(); ++col)
					if (fabs(y_coeff_loc[col]) > coeff_eps)
					{
						this->matind[nzc] = this->v_x[col];
						this->matval[nzc] = x_coeff_loc[col];
						++nzc;
					}

				for (unsigned short col = 0; col < this->problem->getNy(); ++col)
					if (fabs(y_coeff_loc[col]) > coeff_eps)
					{
						this->matind[nzc] = this->v_y_k(k, col);
						this->matval[nzc] = y_coeff_loc[col];
						++nzc;
					}

				status = CPXaddrows(this->cpx_env_master, this->cpx_lp_master, 0, 1, nzc, rhs, sense, matbeg, this->matind.data(), this->matval.data(), NULL, cnstrname);
				if (checkCPXstatus(this->cpx_env_master, status)) goto TERMINATE;

			}
	


TERMINATE:

	if (status)
	{
		return false;
	}

	return true;
}

bool RobustKSolver::add_dual_vars_valid_inequalities(unsigned int z)
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

	double coeff_eps = COEFF_EPS;
	double big_m = BIG_M;
	double alpha_bound_eps = alphaBoundEps;
	status = 0;

	int cnt = 0;
	char lu[2];
	double bd[2];
	unsigned int num_k;

	bool use_alpha_symmetry_breaking =this->cfg->getValueOfKey<bool>("USE_ALPHA_SYMMETRY_BREAKING");
	bool use_tighter_mccormick = this->cfg->getValueOfKey<bool>("USE_TIGHTER_MCCORMICK");
	
	// alpha_k >= alpha_(k+1)
	sense[0] = 'G';
	rhs[0] = 0;
	nzc = 0;

	if (use_alpha_symmetry_breaking) {
		for (int k = 0; k < this->problem->getK() - 1; ++k)
		{
			nzc = 0;

			sprintf(cnstrname[0], "Alpha_k%d_geq_alpha_kk%d ", k, k + 1);

			this->matind[nzc] = this->v_dual_set[z].v_alpha[k];
			this->matval[nzc] = 1;
			++nzc;

			this->matind[nzc] = this->v_dual_set[z].v_alpha[k + 1];
			this->matval[nzc] = -1;
			++nzc;

			if (nzc > 0)
			{
				status = CPXaddrows(this->cpx_env_master, this->cpx_lp_master, 0, 1, nzc, rhs, sense, matbeg, this->matind.data(), this->matval.data(), NULL, cnstrname);
				if (checkCPXstatus(this->cpx_env_master, status)) goto TERMINATE;
			}
		}
	}
	/*
	Tight alpha bounds
	*/
	//num_k = this->problem->getK();
	//for (int k = 0; k < this->problem->getK(); ++k)
	//	if(false)
	//{
	//	cnt = 0;
	//	this->matind[cnt] = this->v_dual_set[z].v_alpha[k];

	//	if (k == 0)
	//	{
	//		bd[cnt] = (1.00 / (double)num_k) - alpha_bound_eps;
	//		lu[cnt] = 'L';
	//		++cnt;
	//		status = CPXchgbds(this->cpx_env_master, this->cpx_lp_master, cnt, this->matind.data(), lu, bd);
	//		if (checkCPXstatus(this->cpx_env_master, status)) goto TERMINATE;

	//		cnt = 0;
	//		bd[cnt] = 1;
	//		lu[cnt] = 'U';
	//		++cnt;
	//		status = CPXchgbds(this->cpx_env_master, this->cpx_lp_master, cnt, this->matind.data(), lu, bd);
	//		if (checkCPXstatus(this->cpx_env_master, status)) goto TERMINATE;
	//		
	//	}
	//	else
	//	{
	//		bd[cnt] = 0;
	//		lu[cnt] = 'L';
	//		++cnt;
	//		status = CPXchgbds(this->cpx_env_master, this->cpx_lp_master, cnt, this->matind.data(), lu, bd);
	//		if (checkCPXstatus(this->cpx_env_master, status)) goto TERMINATE;

	//		cnt = 0;
	//		bd[cnt] = (1.00 / (double)(k+1)) + alpha_bound_eps;
	//		lu[cnt] = 'U';
	//		++cnt;
	//		status = CPXchgbds(this->cpx_env_master, this->cpx_lp_master, cnt, this->matind.data(), lu, bd);
	//		if (checkCPXstatus(this->cpx_env_master, status)) goto TERMINATE;
	//	}
	//}
	
	
	/*
	McCormik: like the old ones but with new UB of alpha
	*/
	/*
	 x * alpha_k
	*/
	if (use_tighter_mccormick)
	{
		if (this->problem->getNx() > 0
			&&
			this->problem->get1S_uobj_matr_C().getNumElem() > 0)
		{
			for (int k = 0; k < this->problem->getK(); ++k)
				for (int xi = 0; xi < this->problem->getNx(); ++xi)
					if (k > 0) // For k = 0 the UB is 1.
					{


						/*
						t_x_k - alpha_k - x >= -1
						*/

						rhs[0] = -((1.00 / (double)(k + 1)) + alpha_bound_eps);

						sense[0] = 'G';
						sprintf(cnstrname[0], "Linear_t_x_(3)_k%d_i%d", k, xi);
						nzc = 0;

						this->matind[nzc] = this->v_dual_set[z].v_t_x_k(k, xi);
						this->matval[nzc] = 1;
						++nzc;

						this->matind[nzc] = this->v_dual_set[z].v_alpha[k];
						this->matval[nzc] = -1;
						++nzc;

						this->matind[nzc] = this->v_x[xi];
						this->matval[nzc] = -((1.00 / (double)(k + 1)) + alpha_bound_eps);
						++nzc;

						status = CPXaddrows(this->cpx_env_master, this->cpx_lp_master, 0, 1, nzc, rhs, sense, matbeg, this->matind.data(), this->matval.data(), NULL, cnstrname);
						if (checkCPXstatus(this->cpx_env_master, status)) goto TERMINATE;
					}
		}

		/*
		* w * alpha_k
		*
		*/
		if (this->problem->getNw() > 0
			&&
			this->problem->getDV_uobj_matr_D().getNumElem() > 0)
		{
			for (int k = 0; k < this->problem->getK(); ++k)
				for (int wi = 0; wi < this->problem->getNw(); ++wi)
					if (k > 0)
					{


						/*
						t_w_k - alpha_k - w(UB alpha) >= -1(UN alpha)
						*/
						rhs[0] = -((1.00 / (double)(k + 1)) + alpha_bound_eps);
						sense[0] = 'G';
						sprintf(cnstrname[0], "Linear_t_x_(3)_k%d_i%d", k, wi);
						nzc = 0;

						this->matind[nzc] = this->v_dual_set[z].v_t_w_k(k, wi);
						this->matval[nzc] = 1;
						++nzc;

						this->matind[nzc] = this->v_dual_set[z].v_alpha[k];
						this->matval[nzc] = -1;
						++nzc;

						this->matind[nzc] = this->v_w[wi];
						this->matval[nzc] = -((1.00 / (double)(k + 1)) + alpha_bound_eps);
						++nzc;

						status = CPXaddrows(this->cpx_env_master, this->cpx_lp_master, 0, 1, nzc, rhs, sense, matbeg, this->matind.data(), this->matval.data(), NULL, cnstrname);
						if (checkCPXstatus(this->cpx_env_master, status)) goto TERMINATE;
					}
		}



		/*
		* y^k * alpha_k
		*/
		for (int k = 0; k < this->problem->getK(); ++k)
			for (int yi = 0; yi < this->problem->getNy(); ++yi)
				if (k > 0)
				{



					/*
					t_y_k - alpha_k - y_k(UB alpha) >= -1(UB alpha)
					*/
					rhs[0] = -((1.00 / (double)(k + 1)) + alpha_bound_eps);
					sense[0] = 'G';
					sprintf(cnstrname[0], "Linear_t_y_(3)_k%d_i%d", k, yi);
					nzc = 0;

					this->matind[nzc] = this->v_dual_set[z].v_t_y_k(k, yi);
					this->matval[nzc] = 1;
					++nzc;

					this->matind[nzc] = this->v_dual_set[z].v_alpha[k];
					this->matval[nzc] = -1;
					++nzc;

					this->matind[nzc] = this->v_y_k(k, yi);
					this->matval[nzc] = -((1.00 / (double)(k + 1)) + alpha_bound_eps);
					++nzc;

					status = CPXaddrows(this->cpx_env_master, this->cpx_lp_master, 0, 1, nzc, rhs, sense, matbeg, this->matind.data(), this->matval.data(), NULL, cnstrname);
					if (checkCPXstatus(this->cpx_env_master, status)) goto TERMINATE;
				}
	}
TERMINATE:
	if (status)
	{
		return false;
	}

	return true;
}

bool RobustKSolver::add_RLT_AlphaXWY_valid_inequalities(unsigned int z)
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
	status = 0;
	double coeff_eps = COEFF_EPS;

	std::vector<double> b_loc;
	Matrix2D<double> Matr_loc;



	// X matrix constraints
	if (this->problem->getLx() > 0
		&&
		this->problem->get1S_uobj_matr_C().getNumElem() > 0)
	{

		b_loc = this->problem->get1S_drhs_vect_bx();
		Matr_loc = this->problem->get1S_dcon_matr_X();

		for (int k = 0; k < this->problem->getK(); ++k)
		{

			for (unsigned short row = 0; row < this->problem->getLx(); ++row)
			{
				sense[0] = 'L';
				nzc = 0;
				sprintf(cnstrname[0], "t_Xx_k%d_%d", k, row);

				// alpha_var: if the rhs differenc than 0 it multiplies for alpha.
				if (fabs(b_loc[row]) > 0)
				{
					rhs[0] = 0;
					this->matind[nzc] = this->v_dual_set[z].v_alpha[k];
					this->matval[nzc] = -b_loc[row];
					++nzc;
				}
				else {

					rhs[0] = b_loc[row];
				}

				// Now we go through the "columns", coefficients of the matrices 

				for (unsigned short cols = 0; cols < this->problem->getNx(); ++cols)
					if (fabs(Matr_loc(row, cols)) > coeff_eps)
					{
						this->matind[nzc] = this->v_dual_set[z].v_t_x_k(k, cols);
						this->matval[nzc] = Matr_loc(row, cols);
						++nzc;
					}
				status = CPXaddrows(this->cpx_env_master, this->cpx_lp_master, 0, 1, nzc, rhs, sense, matbeg, this->matind.data(), this->matval.data(), NULL, cnstrname);
				if (checkCPXstatus(this->cpx_env_master, status)) goto TERMINATE;
			}

		}
	}
	// W matrices constraints 
	if (this->problem->getLw() > 0
		&&
		this->problem->getDV_uobj_matr_D().getNumElem() > 0)
	{
		b_loc = this->problem->getDV_drhs_vect_bw();
		Matr_loc = this->problem->getDV_dcon_matr_W();

		for (int k = 0; k < this->problem->getK(); ++k)
		{
			for (unsigned short row = 0; row < this->problem->getLw(); ++row)
			{
				sense[0] = 'L';
				nzc = 0;
				sprintf(cnstrname[0], "t_Ww_%dk_%d", k, row);
				// Now we go trought the "columns", coefficient of the matrices 

					// alpha_var: if the rhs differenc than 0 it multiplies for alpha.
				if (fabs(b_loc[row]) > 0)
				{
					rhs[0] = 0;
					this->matind[nzc] = this->v_dual_set[z].v_alpha[k];
					this->matval[nzc] = -b_loc[row];
					++nzc;
				}
				else {

					rhs[0] = b_loc[row];
				}


				for (unsigned short col = 0; col < this->problem->getNw(); ++col)
					if (fabs(Matr_loc(row, col)) > coeff_eps)
					{
						this->matind[nzc] = this->v_dual_set[z].v_t_w_k(k, col);
						this->matval[nzc] = Matr_loc(row, col);
						++nzc;
					}

				status = CPXaddrows(this->cpx_env_master, this->cpx_lp_master, 0, 1, nzc, rhs, sense, matbeg, this->matind.data(), this->matval.data(), NULL, cnstrname);
				if (checkCPXstatus(this->cpx_env_master, status)) goto TERMINATE;
			}
		}
	}

	// Y matrices constraints for each k

	b_loc = this->problem->get2S_drhs_vect_by();
	Matr_loc = this->problem->get2S_dcon_matr_Y();

	if (this->problem->getLy() > 0
		&&
		this->problem->get2S_uobj_matr_Q().getNumElem())
	{
	for (int k = 0; k < this->problem->getK(); ++k)
	{
		for (unsigned short row = 0; row < this->problem->getLy(); ++row)
		{

			sense[0] = 'L';
			nzc = 0;
			sprintf(cnstrname[0], "t_Yy_k%d_%d", k, row);

			if (fabs(b_loc[row]) > 0)
			{
				rhs[0] = 0;
				this->matind[nzc] = this->v_dual_set[z].v_alpha[k];
				this->matval[nzc] = -b_loc[row];
				++nzc;
			}
			else {

				rhs[0] = b_loc[row];
			}
			
			
			// Now we go trought the "columns", coefficient of the matrices 

			for (unsigned short cols = 0; cols < this->problem->getNy(); ++cols)
				if (fabs(Matr_loc(row, cols)) > coeff_eps)
				{
					this->matind[nzc] = this->v_dual_set[z].v_t_y_k(k, cols);
					this->matval[nzc] = Matr_loc(row, cols);
					++nzc;
				}
			status = CPXaddrows(this->cpx_env_master, this->cpx_lp_master, 0, 1, nzc, rhs, sense, matbeg, this->matind.data(), this->matval.data(), NULL, cnstrname);
			if (checkCPXstatus(this->cpx_env_master, status)) goto TERMINATE;

		}
	}
	}
TERMINATE:
	if (status)
	{
		return false;
	}

	return true;

	
}

bool RobustKSolver::add_RLT_AlphaExwy_valid_inequalities(unsigned int z)
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
	status = 0;
	double coeff_eps = COEFF_EPS;

	std::vector<double> b_loc = this->problem->get1S_drhs_vect_bx();
	Matrix2D<double> Matr_loc = this->problem->get1S_dcon_matr_X();




	if (this->problem->getR() > 0)
	{
		b_loc = this->problem->getLink_dcon_vect_g();


		for (int row = 0; row < this->problem->getR(); ++row)
			for (int k = 0; k < this->problem->getK(); ++k)
			{

				sense[0] = 'L';
				nzc = 0;
				rhs[0] = 0;

				if (fabs(b_loc[row]) > coeff_eps)
				{
					this->matind[nzc] = this->v_dual_set[z].v_alpha[k];
					this->matval[nzc] = -b_loc[row];
					++nzc;
				}

				sprintf(cnstrname[0], "LinkXWY_k%d_r%d", k, row);

				Matr_loc = this->problem->getLink_dcon_matr_Ex();
				if (Matr_loc.getNumElem() > 0)
				{
					for (int col = 0; col < this->problem->getNx(); ++col)
						if (fabs(Matr_loc(row, col)) > coeff_eps)
						{
							//this->matind[nzc] = this->v_x[col];
							this->matind[nzc] = this->v_dual_set[z].v_t_x_k(k, col);
							this->matval[nzc] = Matr_loc(row, col);
							++nzc;
						}
				}
				Matr_loc = this->problem->getLink_dcon_matr_Ew();
				if (Matr_loc.getNumElem() > 0)
				{
					for (int col = 0; col < this->problem->getNw(); ++col)
						if (fabs(Matr_loc(row, col)) > coeff_eps)
						{
							//this->matind[nzc] = this->v_w[col];
							this->matind[nzc] = this->v_dual_set[z].v_t_w_k(k, col);
							this->matval[nzc] = Matr_loc(row, col);
							++nzc;
						}
				}

				Matr_loc = this->problem->getLink_dcon_matr_Ey();
				if (Matr_loc.getNumElem() > 0)
				{
					for (int col = 0; col < this->problem->getNy(); ++col)
						if (fabs(Matr_loc(row, col)) > coeff_eps)
						{
							//this->matind[nzc] = this->v_y_k(k, col);
							this->matind[nzc] = this->v_dual_set[z].v_t_y_k(k, col);
							this->matval[nzc] = Matr_loc(row, col);
							++nzc;
						}
				}

				status = CPXaddrows(this->cpx_env_master, this->cpx_lp_master, 0, 1, nzc, rhs, sense, matbeg, this->matind.data(), this->matval.data(), NULL, cnstrname);
				if (checkCPXstatus(this->cpx_env_master, status)) goto TERMINATE;
			}
	}


TERMINATE:
	if (status)
		return false;
	else
		return true;

}

bool RobustKSolver::add_RLT_Primal_valid_inequalities(unsigned int z)
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
	status = 0;
	double coeff_eps = COEFF_EPS;

	std::vector<double> y_coeff_loc;
	std::vector<double> x_coeff_loc;


	y_coeff_loc = this->problem->get2S_ducon_vect_p();
	x_coeff_loc = this->problem->get1S_ducon_vect_t();

	for (int k = 0; k < this->problem->getK(); ++k)
		for (int l = 0; l < this->problem->getL(); ++l)
		{

			rhs[0] = 0;
			sense[0] = 'L';
			nzc = 0;

			sprintf(cnstrname[0], "RLT_Det_CU_cut_k%d_l%d", k, l);

			if (fabs(this->problem->get_rhs_ducon_vect_hmax()[l]) > coeff_eps)
			{
				this->matind[nzc] = this->v_dual_set[z].v_alpha[k];
				this->matval[nzc] = -this->problem->get_rhs_ducon_vect_hmax()[l];
				++nzc;
			}
				
		

			for (unsigned short col = 0; col < this->problem->getNx(); ++col)
				if (fabs(y_coeff_loc[col]) > coeff_eps)
				{
					this->matind[nzc] = this->v_dual_set[z].v_t_x_k(k, col);
					this->matval[nzc] = x_coeff_loc[col];
					++nzc;
				}

			for (unsigned short col = 0; col < this->problem->getNy(); ++col)
				if (fabs(y_coeff_loc[col]) > coeff_eps)
				{
					this->matind[nzc] = this->v_dual_set[z].v_t_y_k(k, col);
					this->matval[nzc] = y_coeff_loc[col];
					++nzc;
				}

			status = CPXaddrows(this->cpx_env_master, this->cpx_lp_master, 0, 1, nzc, rhs, sense, matbeg, this->matind.data(), this->matval.data(), NULL, cnstrname);
			if (checkCPXstatus(this->cpx_env_master, status)) goto TERMINATE;

		}


TERMINATE:
	if (status)
		return false;
	else
		return true;
}

bool RobustKSolver::add_prev_iter_Phi_LB(unsigned int z)
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
	double coeff_eps = COEFF_EPS;

	std::vector<double> x_coeff_loc = this->problem->get1S_dobj_vect_c();


	rhs[0] = this->master_curr_blb;
	sense[0] = 'G';
	nzc = 0;

	sprintf(cnstrname[0], "Phi_lb_iter_%d", z-1);

	this->matind[nzc] = this->v_phi;
	this->matval[nzc] = 1;
	++nzc;


	for (unsigned short cols = 0; cols < this->problem->getNx(); ++cols)
		if (fabs(x_coeff_loc[cols]) > coeff_eps)
		{
			this->matind[nzc] = this->v_x[cols];
			this->matval[nzc] = x_coeff_loc[cols];
			++nzc;
		}



	status = CPXaddrows(this->cpx_env_master, this->cpx_lp_master, 0, 1, nzc, rhs, sense, matbeg, this->matind.data(), this->matval.data(), NULL, cnstrname);
	if (checkCPXstatus(this->cpx_env_master, status)) goto TERMINATE;

TERMINATE:
	if (status)
		return false;
	else
		return true;
}

bool RobustKSolver::master_update_problem(unsigned int z)
{
	bool check = true;
	int status;



	/*
	Variables
	*/
	check = this->add_dualization_variables(v_dual_set.size()-1);
	if (!check)
	{
		goto TERMINATE;
	}

	check = this->add_bilinear_variables(v_dual_set.size() - 1);
	if (!check)
	{
		goto TERMINATE;
	}


	/*
	Constraints
	*/
	check = this->add_dualization_constraints(v_dual_set.size() - 1);
	if (!check)
	{
		goto TERMINATE;
	}
	check = this->add_linearization_constraints(v_dual_set.size() - 1);
	if (!check)
	{
		goto TERMINATE;
	}

	if (v_dual_set.size() - 1 == 0)
	{
		check = this->add_dual_vars_valid_inequalities(v_dual_set.size() - 1);
		if (!check)
		{
			goto TERMINATE;
		}
	}

	check = this->add_RLT_AlphaXWY_valid_inequalities(v_dual_set.size() - 1);
	if (!check)
	{
		goto TERMINATE;
	}

	check = add_RLT_AlphaExwy_valid_inequalities(v_dual_set.size() - 1);
	if (!check)
	{
		goto TERMINATE;
	}

	check = add_prev_iter_Phi_LB(v_dual_set.size() - 1);
	if (!check)
	{
		goto TERMINATE;
	}


	status = CPXwriteprob(this->cpx_env_master, this->cpx_lp_master, "master.lp", "lp");
	if (this->checkCPXstatus(this->cpx_env_master, status)) { check = false; goto TERMINATE; };

TERMINATE:
	if (!check)
		return false;
	else
		return true;
}

bool RobustKSolver::master_solve_problem(double* obj_val, double* unc_cost)
{

	bool check = true;
	int status;
	int mip_status;

	status = CPXmipopt(this->cpx_env_master, this->cpx_lp_master);
	if (checkCPXstatus(this->cpx_env_master, status)) goto TERMINATE;






	mip_status = CPXgetstat(this->cpx_env_master, this->cpx_lp_master);

	status = CPXgetx(this->cpx_env_master, this->cpx_lp_master, this->master_curr_sol.data(),
		0, CPXgetnumcols(cpx_env_master, cpx_lp_master) - 1);
	if (checkCPXstatus(this->cpx_env_master, status)) goto TERMINATE;

	*unc_cost = master_curr_sol[this->v_phi];


	double best_lb;

	status = CPXgetbestobjval(this->cpx_env_master, this->cpx_lp_master,&this->master_curr_blb);
	if (checkCPXstatus(this->cpx_env_master, status)) goto TERMINATE;

	status = CPXgetobjval(this->cpx_env_master, this->cpx_lp_master, obj_val);
	if (checkCPXstatus(this->cpx_env_master, status)) goto TERMINATE;

	this->master_curr_bub = *obj_val;


TERMINATE:
	if (!check)
		return false;
	else
		return true;
}

bool RobustKSolver::slave_init_cplex_data_str()
{
	int status = 0;

	this->cpx_env_slave = CPXopenCPLEX(&status);
	if (this->checkCPXstatus(this->cpx_env_slave, status)) return false;
	this->cpx_lp_slave = CPXcreateprob(this->cpx_env_slave, &status, this->problem->getProbName().c_str());
	if (checkCPXstatus(this->cpx_env_slave, status))return false;

	return true;
}

bool RobustKSolver::slave_set_cplex_parameters()
{
	bool mod_stat = true;

	int status = 0;

	char errbuf[CPXMESSAGEBUFSIZE];


	status = CPXsetintparam(this->cpx_env_slave, CPX_PARAM_SCRIND, CPX_ON);
	if (checkCPXstatus(this->cpx_env_slave, status)) return false;
	status = CPXsetintparam(this->cpx_env_slave, CPX_PARAM_THREADS, 1);
	if (checkCPXstatus(this->cpx_env_slave, status)) return false;
	status = CPXchgobjsen(this->cpx_env_slave, this->cpx_lp_slave, CPX_MAX);
	if (checkCPXstatus(this->cpx_env_slave, status)) return false;
	status = CPXsetintparam(this->cpx_env_slave, CPX_PARAM_NUMERICALEMPHASIS, CPX_ON);
	if (checkCPXstatus(this->cpx_env_slave, status)) return false;
	status = CPXsetdblparam(this->cpx_env_slave, CPXPARAM_Simplex_Tolerances_Feasibility, 1e-9);
	if (checkCPXstatus(this->cpx_env_slave, status)) return false;
	status = CPXsetdblparam(this->cpx_env_slave, CPXPARAM_Simplex_Tolerances_Optimality, 1e-9);
	if (checkCPXstatus(this->cpx_env_slave, status)) return false;
	status = CPXsetdblparam(this->cpx_env_slave, CPXPARAM_MIP_Tolerances_MIPGap, 1e-9);
	if (checkCPXstatus(this->cpx_env_slave, status)) return false;
	status = CPXsetdblparam(this->cpx_env_slave, CPXPARAM_MIP_Tolerances_Integrality, 0);
	if (checkCPXstatus(this->cpx_env_slave, status)) return false;



	//status = CPXsetdblparam(this->cpx_env_slave, CPX_PARAM_EPRHS, 1e-9);
	//if (checkCPXstatus(this->cpx_env_slave, status)) return false;
	//status = CPXsetdblparam(this->cpx_env_slave, CPX_PARAM_EPOPT, 1e-9);
	//if (checkCPXstatus(this->cpx_env_slave, status)) return false;
	//status = CPXsetdblparam(this->cpx_env_slave, CPX_PARAM_EPGAP, 1e-9);
	//if (checkCPXstatus(this->cpx_env_slave, status)) return false;
	//status = CPXsetdblparam(this->cpx_env_slave, CPXPARAM_MIP_Tolerances_Integrality, 1e-12);
	//if (checkCPXstatus(this->cpx_env_slave, status)) return false;
	

	//status = CPXsetdblparam(g_spp_milp.cpx_env, CPXPARAM_TimeLimit, TIME_LIMIT);
	//if (checkStatus(g_spp_milp.cpx_env, status)) goto TERMINATE;


	return true;
}

bool RobustKSolver::slave_build_problem()
{
	bool check = true;
	int status;



	check = this->slave_init_cplex_data_str();
	if (!check)
	{
		goto TERMINATE;
	}

	check = this->slave_set_cplex_parameters();
	if (!check)
	{
		goto TERMINATE;
	}


	/*
	Variables
	*/
	check = this->slave_add_variables();
	if (!check)
	{
		goto TERMINATE;
	}

	
	/*
	Constraints
	*/
	check = this->slave_add_constraints();
	if (!check)
	{
		goto TERMINATE;
	}
	check = this->slave_add_linearization_constraints();
	if (!check)
	{
		goto TERMINATE;
	}

	status = CPXwriteprob(this->cpx_env_slave, this->cpx_lp_slave, "slave.lp", "lp");
	if (this->checkCPXstatus(this->cpx_env_slave, status)) { check = false; goto TERMINATE; };

TERMINATE:
	if (!check)
		return false;
	else
		return true;
}

bool RobustKSolver::slave_add_variables()
{
	unsigned short num_cols = CPXgetnumcols(this->cpx_env_slave, this->cpx_lp_slave);

	int status = 0;
	int vars_ind = 0;
	double lb[1], ub[1];
	double obj[1];
	char* varsname[1];
	char elname[1024];
	varsname[0] = elname;
	char vartype[1];

	/*
	Tau variable
	*/
	this->vs_tau = num_cols;
	obj[0] = 1;
	vartype[0] = 'C';
	sprintf(varsname[0], "Tau");
	lb[0] = -CPX_INFBOUND;
	ub[0] = CPX_INFBOUND;
	status = CPXnewcols(this->cpx_env_slave, this->cpx_lp_slave, 1, obj, lb, ub, NULL, varsname);
	if (checkCPXstatus(this->cpx_env_slave, status)) goto TERMINATE;
	++num_cols;

	/*ov_Xi variables*/
	this->vs_ov_xi.resize(this->problem->getNxi());
	for (unsigned int n = 0; n < this->problem->getNxi(); ++n)
	{
		this->vs_ov_xi[n] = num_cols;
		obj[0] = 0;
		vartype[0] = 'C';
		sprintf(varsname[0], "ov_xi_n%d", n);
		lb[0] = -CPX_INFBOUND;
		ub[0] = CPX_INFBOUND;
		status = CPXnewcols(this->cpx_env_slave, this->cpx_lp_slave, 1, obj, lb, ub, NULL, varsname);
		if (checkCPXstatus(this->cpx_env_slave, status)) goto TERMINATE;
		++num_cols;

	}
	/*ov_zeta */
	this->vs_ov_zeta.resize(this->problem->getNzeta());
	for (unsigned int n = 0; n < this->problem->getNzeta(); ++n)
	{
		this->vs_ov_zeta[n] = num_cols;
		obj[0] = 0;
		vartype[0] = 'C';
		sprintf(varsname[0], "ov_xi_n%d", n);
		lb[0] = -CPX_INFBOUND;
		ub[0] = CPX_INFBOUND;
		status = CPXnewcols(this->cpx_env_slave, this->cpx_lp_slave, 1, obj, lb, ub, NULL, varsname);
		if (checkCPXstatus(this->cpx_env_slave, status)) goto TERMINATE;
		++num_cols;
	}
	/*xi_k */
	this->vs_xi_k.resizeMatrix2D(this->problem->getK(), this->problem->getNxi());
	for (int k = 0; k < this->problem->getK(); ++k)
		for (unsigned int n = 0; n < this->problem->getNxi(); ++n)
		{
			this->vs_xi_k(k, n) = num_cols;
			obj[0] = 0;
			vartype[0] = 'C';
			sprintf(varsname[0], "xi_k%d_n%d", k, n);

			lb[0] = -CPX_INFBOUND;
			ub[0] = CPX_INFBOUND;

			status = CPXnewcols(this->cpx_env_slave, this->cpx_lp_slave, 1, obj, lb, ub, NULL, varsname);
			if (checkCPXstatus(this->cpx_env_slave, status)) goto TERMINATE;
			++num_cols;
		}
	/*zeta_k */
	if (this->problem->getNzeta() > 0)
	{
	this->vs_zeta_k.resizeMatrix2D(this->problem->getK(), this->problem->getNzeta());

	for (int k = 0; k < this->problem->getK(); ++k)
		for (unsigned int n = 0; n < this->problem->getNzeta(); ++n)
		{
			this->vs_zeta_k(k, n) = num_cols;
			obj[0] = 0;
			vartype[0] = 'C';
			sprintf(varsname[0], "zeta_k%d_n%d", k, n);

			lb[0] = -CPX_INFBOUND;
			ub[0] = CPX_INFBOUND;

			status = CPXnewcols(this->cpx_env_slave, this->cpx_lp_slave, 1, obj, lb, ub, NULL, varsname);
			if (checkCPXstatus(this->cpx_env_slave, status)) goto TERMINATE;
			++num_cols;
		}
	}

	/*z_k*/
	this->vs_z_k.resizeMatrix2D(this->problem->getK(), this->problem->getL());

	for (int k = 0; k < this->problem->getK(); ++k)
		for (unsigned int l = 0; l < this->problem->getL(); ++l)
		{
			this->vs_z_k(k, l) = num_cols;
			obj[0] = 0;
			vartype[0] = 'B';
			sprintf(varsname[0], "z_k%d_l%d", k, l);

			lb[0] = 0;
			ub[0] = 1;

			status = CPXnewcols(this->cpx_env_slave, this->cpx_lp_slave, 1, obj, lb, ub, NULL, varsname);
			if (checkCPXstatus(this->cpx_env_slave, status)) goto TERMINATE;
			++num_cols;
		}

	/*chi^{kl}_n*/
	this->vs_chi_kl.resizeMatrix3D(this->problem->getK(), this->problem->getL(), this->problem->getNxi());

	for (unsigned short k = 0; k < this->problem->getK(); ++k)
		for (unsigned short l = 0; l < this->problem->getL(); ++l)
			for (unsigned short n = 0; n < this->problem->getNxi(); ++n)
			{
				this->vs_chi_kl(k, l, n) = num_cols;
				obj[0] = 0;
				vartype[0] = 'C';
				sprintf(varsname[0], "chi_k%d_l%d_n%d", k, l,n);

				lb[0] = -CPX_INFBOUND;
				ub[0] =  CPX_INFBOUND;

				status = CPXnewcols(this->cpx_env_slave, this->cpx_lp_slave, 1, obj, lb, ub, NULL, varsname);
				if (checkCPXstatus(this->cpx_env_slave, status)) goto TERMINATE;
				++num_cols;
			}

TERMINATE:

	if (status)
	{
		return false;
	}

	return true;
}

bool RobustKSolver::slave_add_constraints()
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
	status = 0;
	double coeff_eps = COEFF_EPS;
	double int_eps = INT_EPS;
	double big_m = CNSTR_UNC_BIG_M;
	double coeff = 0;

	std::vector<double> b_loc;
	Matrix2D<double> Matr2D_loc;
	Matrix3D<double> Matr3D_loc;



	/*Tau constraints: at the end*/
	for (unsigned int k = 0; k < this->problem->getK(); ++k)
	{
		rhs[0] = 0;
		sense[0] = 'L';
		nzc = 0;

		sprintf(cnstrname[0], "ConstrTau_k%d", k);


		b_loc = this->problem->get2S_dobj_vect_q();

		for (unsigned int i = 0; i < this->problem->getNy(); ++i)
			if (master_curr_sol[this->v_y_k(k, i)] > int_eps)
			{ 

				rhs[0] += (b_loc[i] * master_curr_sol[this->v_y_k(k, i)]);
				//rhs[0] += (b_loc[i]);
				//* master_curr_sol[this->v_y_k(k, i)]);
			}

		

		b_loc = this->problem->get_rhs_ucon_vect_h();

		this->matind[nzc] = this->vs_tau;
		this->matval[nzc] = 1;
		++nzc;


		for (unsigned int n = 0; n < this->problem->getNxi(); ++n)
		{
			//this->matind[nzc] = this->vs_xi_k(k, n);
			coeff = 0;
			//this->matval[nzc] = 0;

		/*\xiCx*/
			Matr2D_loc = this->problem->get1S_uobj_matr_C();
			for (unsigned int i = 0; i < this->problem->getNx(); ++i)
				if (
					fabs(Matr2D_loc(n, i)) > coeff_eps
					&&
					fabs(this->master_curr_sol[this->v_x[i]]) > int_eps
					)
				{
						//coeff += -Matr2D_loc(n, i); /*With this we do not multiply by the var val, becase it shoud be one.*/

						coeff += (-(Matr2D_loc(n, i) * master_curr_sol[this->v_x[i]]));
					
				}

			/*\xiDw*/
			Matr2D_loc = this->problem->getDV_uobj_matr_D();
			for (unsigned int i = 0; i < this->problem->getNw(); ++i)
				if (
					fabs(Matr2D_loc(n, i)) > coeff_eps
					&&
					fabs(this->master_curr_sol[this->v_w[i]]) > int_eps
				
					)
				{
					//coeff += -Matr2D_loc(n, i);

					coeff += (-(Matr2D_loc(n, i) * master_curr_sol[this->v_w[i]]));
					
					
				}


			/*\xiQy*/

			Matr2D_loc = this->problem->get2S_uobj_matr_Q();
			for (unsigned int i = 0; i < this->problem->getNy(); ++i)
				if (
					fabs(Matr2D_loc(n, i)) > coeff_eps
					&&
					fabs(this->master_curr_sol[this->v_y_k(k, i)]) > int_eps
					)
				{

					//coeff += -Matr2D_loc(n, i);
					
					coeff += (-(Matr2D_loc(n, i) * master_curr_sol[this->v_y_k(k, i)]));
					
					
				}

			if (fabs(coeff) > coeff_eps)
			{
				this->matind[nzc] = this->vs_xi_k(k, n);
				this->matval[nzc] = coeff;
				++nzc;
			}
		}


		/*zk \top (Tx + Vw + Pyk -h)M*/
		for (unsigned int l = 0; l < this->problem->getL(); ++l)
		{
			//this->matind[nzc] = this->vs_z_k(k, l);
			coeff = 0;
			/*z_k_l h_l*/
			coeff += (b_loc[l]*big_m);

			/*z_k_l * Tx */
			Matr2D_loc = this->problem->get1S_ucon_matr_T();
			for (unsigned int i = 0; i < this->problem->getNx(); ++i)
				if (fabs(Matr2D_loc(l, i)) > coeff_eps
					&&
					fabs(this->master_curr_sol[this->v_x[i]]) > int_eps
					)
			{
				coeff += (-(Matr2D_loc(l, i) * master_curr_sol[this->v_x[i]]) * big_m);
			}

			/*z_k_l * Vw */
			Matr2D_loc = this->problem->getDV_ucon_matr_V();
			for (unsigned int i = 0; i < this->problem->getNw(); ++i)
				if (fabs(Matr2D_loc(l, i)) > coeff_eps
					&&
					fabs(this->master_curr_sol[this->v_w[i]]) > int_eps
					)
				{
					coeff += (-(Matr2D_loc(l, i) * master_curr_sol[this->v_w[i]]) * big_m);

				}

			/*z_k_l * Py */
			Matr2D_loc = this->problem->get2S_ucon_matr_P();
			for (unsigned int i = 0; i < this->problem->getNy(); ++i)
				if (fabs(Matr2D_loc(l, i)) > coeff_eps
					&&
					fabs(this->master_curr_sol[this->v_y_k(k,i)]) > int_eps
					)
				{
					coeff += (-(Matr2D_loc(l, i) * master_curr_sol[this->v_y_k(k, i)] * big_m));

				}

			if (fabs(coeff) > coeff_eps)
			{
				this->matind[nzc] = this->vs_z_k(k, l);
				this->matval[nzc] = coeff;
				//this->matval[nzc] = coeff * big_m;
				++nzc;
			}

		}
		/*\chi times the mess: \chi_(kln)*/
		//Matr3D_loc
		Matr2D_loc = this->problem->getUC_ucon_matr_H();

		for(unsigned int l = 0; l < this->problem->getL(); ++l)
			for (unsigned int n = 0; n < this->problem->getNxi(); ++n)
			{
				//this->matind[nzc] = this->vs_chi_kl(k, l, n);
				coeff = 0;

				/*Hln * \chi*/
				coeff += (Matr2D_loc(l, n) * big_m);


				Matr3D_loc = this->problem->get1S_ucon_matr_Txi();
				for (unsigned int i = 0; i < this->problem->getNx(); ++i)
					if (
						fabs(Matr3D_loc(n,l,i)) > coeff_eps
						&&
						fabs(this->master_curr_sol[this->v_x[i]]) > int_eps)
					{

						coeff += (-(Matr3D_loc(n, l, i) * master_curr_sol[this->v_x[i]]) * big_m);
					}

				Matr3D_loc = this->problem->getDV_ucon_matr_Vxi();
				for (unsigned int i = 0; i < this->problem->getNw(); ++i)
					if (
						fabs(Matr3D_loc(n, l, i)) > coeff_eps
						&&
						fabs(this->master_curr_sol[this->v_w[i]]) > int_eps)
					{

						coeff += (-(Matr3D_loc(n, l, i) * master_curr_sol[this->v_w[i]]) * big_m);
					}

				Matr3D_loc = this->problem->get2S_ucon_matr_Pxi();
				for (unsigned int i = 0; i < this->problem->getNy(); ++i)
					if (
						fabs(Matr3D_loc(n, l, i)) > coeff_eps
						&&
						fabs(this->master_curr_sol[this->v_y_k(k,i)]) > int_eps)
					{

						coeff += (-(Matr3D_loc(n, l, i) * master_curr_sol[this->v_y_k(k,i)]) * big_m);
					}


				if (fabs(coeff) > coeff_eps)
				{
					this->matind[nzc] = this->vs_chi_kl(k, l, n);
					this->matval[nzc] = coeff;
					//this->matval[nzc] = coeff * big_m;
					++nzc;
				}


			}

		status = CPXaddrows(this->cpx_env_slave, this->cpx_lp_slave, 0, 1, nzc, rhs, sense, matbeg, this->matind.data(), this->matval.data(), NULL, cnstrname);
		if (checkCPXstatus(this->cpx_env_slave, status)) goto TERMINATE;
	}

	/*Uncertainty set constraints*/
	b_loc = this->problem->getUC_rhs_vect_b();
	Matr2D_loc = this->problem->getUC_xicon_matr_A();
	for (unsigned int row = 0; row < this->problem->getLxi(); ++row)
	{
		rhs[0] = b_loc[row];
		sense[0] = 'L';
		nzc = 0;
		sprintf(cnstrname[0], "UNCset_%d", row);

		for (unsigned int col = 0; col < this->problem->getNxi(); ++col)
			if(fabs(Matr2D_loc(row,col)) > coeff_eps)
		{
				this->matind[nzc] = this->vs_ov_xi[col];
				this->matval[nzc] = Matr2D_loc(row, col);
				++nzc;
		}
		/*
		Add G and zeta here
		*/

		status = CPXaddrows(this->cpx_env_slave, this->cpx_lp_slave, 0, 1, nzc, rhs, sense, matbeg, this->matind.data(), this->matval.data(), NULL, cnstrname);
		if (checkCPXstatus(this->cpx_env_slave, status)) goto TERMINATE;
	}

	/*Uncertainty set constraints for each k*/
	for (unsigned int k = 0; k < this->problem->getK(); ++k)
		for (unsigned int row = 0; row < this->problem->getLxi(); ++row)
	{
		rhs[0] = b_loc[row];
		sense[0] = 'L';
		nzc = 0;
		sprintf(cnstrname[0], "UNCset_k%d_%d", k,row);

		for (unsigned int col = 0; col < this->problem->getNxi(); ++col)
			if (fabs(Matr2D_loc(row, col)) > coeff_eps)
			{
				this->matind[nzc] = this->vs_xi_k(k, col);
				this->matval[nzc] = Matr2D_loc(row, col);
				++nzc;
			}
		/*
		Add G and zeta here
		*/

		status = CPXaddrows(this->cpx_env_slave, this->cpx_lp_slave, 0, 1, nzc, rhs, sense, matbeg, this->matind.data(), this->matval.data(), NULL, cnstrname);
		if (checkCPXstatus(this->cpx_env_slave, status)) goto TERMINATE;
	}

	/* w\xi = w\xi_k*/
	for (unsigned int k = 0; k < this->problem->getK(); ++k)
		for (unsigned int n = 0; n < this->problem->getNxi(); ++n)
			if(this->master_curr_sol[this->v_w[n]] > 0)
		{
				rhs[0] = 0;
				sense[0] = 'E';
				nzc = 0;

				sprintf(cnstrname[0], "DiscoveryConstr_k%d_%d", k, n);


				this->matind[nzc] = this->vs_ov_xi[n];
				this->matval[nzc] = 1;
				++nzc;

				this->matind[nzc] = this->vs_xi_k(k,n);
				this->matval[nzc] = -1;
				++nzc;

				status = CPXaddrows(this->cpx_env_slave, this->cpx_lp_slave, 0, 1, nzc, rhs, sense, matbeg, this->matind.data(), this->matval.data(), NULL, cnstrname);
				if (checkCPXstatus(this->cpx_env_slave, status)) goto TERMINATE;

		}

	/*zk <= 1*/
	for (unsigned int k = 0; k < this->problem->getK(); ++k)
	{
		nzc = 0;
		rhs[0] = 1;
		sense[0] = 'L';
		sprintf(cnstrname[0], "Z_leq_1_k%d", k);

		for(unsigned int l = 0; l < this->problem->getL(); ++l)
		{
			this->matind[nzc] = this->vs_z_k(k, l);
			this->matval[nzc] = 1;
			++nzc;
		}
		status = CPXaddrows(this->cpx_env_slave, this->cpx_lp_slave, 0, 1, nzc, rhs, sense, matbeg, this->matind.data(), this->matval.data(), NULL, cnstrname);
		if (checkCPXstatus(this->cpx_env_slave, status)) goto TERMINATE;
	}


TERMINATE:

	if (status)
	{
		return false;
	}

	return true;
}

bool RobustKSolver::slave_add_linearization_constraints()
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
	status = 0;
	double coeff_eps = COEFF_EPS;
	double int_eps = INT_EPS;
	double big_m = CNSTR_UNC_BIG_M;
	double coeff = 0;

	std::vector<double> b_loc;
	Matrix2D<double> Matr2D_loc;
	Matrix3D<double> Matr3D_loc;

	b_loc = this->problem->get_ub_xi();


	/*\chi_{kln} <= z_{kl}UB_\xi  */
	b_loc = this->problem->get_ub_xi();
	for(unsigned int k = 0; k < this->problem->getK(); ++k)
		for (unsigned int l = 0; l < this->problem->getL(); ++l)
			for(unsigned int n = 0; n < this->problem->getNxi(); ++n)
		{
			rhs[0] = 0;
			nzc = 0;
			sense[0] = 'L';
			sprintf(cnstrname[0], "bil(1)_chi_geq_zUBxi_k%d_l%d_n%d", k,l,n);

			this->matind[nzc] = this->vs_chi_kl(k, l, n);
			this->matval[nzc] = 1;
			++nzc;

			this->matind[nzc] = this->vs_z_k(k, l);
			this->matval[nzc] = -b_loc[n];
			++nzc;

			status = CPXaddrows(this->cpx_env_slave, this->cpx_lp_slave, 0, 1, nzc, rhs, sense, matbeg, this->matind.data(), this->matval.data(), NULL, cnstrname);
			if (checkCPXstatus(this->cpx_env_slave, status)) goto TERMINATE;
		}


	/*\chi_{kln} >= z_{kl}LB_\xi  */
	b_loc = this->problem->get_lb_xi();
	for (unsigned int k = 0; k < this->problem->getK(); ++k)
		for (unsigned int l = 0; l < this->problem->getL(); ++l)
			for (unsigned int n = 0; n < this->problem->getNxi(); ++n)
			{
				rhs[0] = 0;
				nzc = 0;
				sense[0] = 'G';
				sprintf(cnstrname[0], "bil(2)_chi_leq_zLBxi_k%d_l%d_n%d", k, l, n);

				this->matind[nzc] = this->vs_chi_kl(k, l, n);
				this->matval[nzc] = 1;
				++nzc;

				this->matind[nzc] = this->vs_z_k(k, l);
				this->matval[nzc] = -b_loc[n];
				++nzc;

				status = CPXaddrows(this->cpx_env_slave, this->cpx_lp_slave, 0, 1, nzc, rhs, sense, matbeg, this->matind.data(), this->matval.data(), NULL, cnstrname);
				if (checkCPXstatus(this->cpx_env_slave, status)) goto TERMINATE;
			}

	/*\chi_{kln} <= z_{kl}LB_\xi + \xi - LB  */
	b_loc = this->problem->get_lb_xi();
	for (unsigned int k = 0; k < this->problem->getK(); ++k)
		for (unsigned int l = 0; l < this->problem->getL(); ++l)
			for (unsigned int n = 0; n < this->problem->getNxi(); ++n)
			{
				rhs[0] = -b_loc[n];
				nzc = 0;
				sense[0] = 'L';
				sprintf(cnstrname[0], "bil(3)_chi_leq_zLBxi_xi_-LB_k%d_l%d_n%d", k, l, n);

				this->matind[nzc] = this->vs_chi_kl(k, l, n);
				this->matval[nzc] = 1;
				++nzc;

				this->matind[nzc] = this->vs_z_k(k, l);
				this->matval[nzc] = -b_loc[n];
				++nzc;

				this->matind[nzc] = this->vs_xi_k(k, n);
				this->matval[nzc] = -1;
				++nzc;


				status = CPXaddrows(this->cpx_env_slave, this->cpx_lp_slave, 0, 1, nzc, rhs, sense, matbeg, this->matind.data(), this->matval.data(), NULL, cnstrname);
				if (checkCPXstatus(this->cpx_env_slave, status)) goto TERMINATE;
			}

	/*\chi_{kln} >= z_{kl}UB_\xi + \xi - UB  */
	b_loc = this->problem->get_ub_xi();
	for (unsigned int k = 0; k < this->problem->getK(); ++k)
		for (unsigned int l = 0; l < this->problem->getL(); ++l)
			for (unsigned int n = 0; n < this->problem->getNxi(); ++n)
			{
				rhs[0] = -b_loc[n];
				nzc = 0;
				sense[0] = 'G';
				sprintf(cnstrname[0], "bil(4)_chi_geq_zLBxi_xi_-UB_k%d_l%d_n%d", k, l, n);

				this->matind[nzc] = this->vs_chi_kl(k, l, n);
				this->matval[nzc] = 1;
				++nzc;

				this->matind[nzc] = this->vs_z_k(k, l);
				this->matval[nzc] = -b_loc[n];
				++nzc;

				this->matind[nzc] = this->vs_xi_k(k, n);
				this->matval[nzc] = -1;
				++nzc;


				status = CPXaddrows(this->cpx_env_slave, this->cpx_lp_slave, 0, 1, nzc, rhs, sense, matbeg, this->matind.data(), this->matval.data(), NULL, cnstrname);
				if (checkCPXstatus(this->cpx_env_slave, status)) goto TERMINATE;
			}


TERMINATE:

	if (status)
	{
		return false;
	}

	return true;
}

bool RobustKSolver::slave_solve_problem(double* obj_val)
{

	bool check = true;
	int status;
	
	status = CPXlpopt(this->cpx_env_slave, this->cpx_lp_slave);
if (checkCPXstatus(this->cpx_env_slave, status)) goto TERMINATE;

	//status = CPXmipopt(this->cpx_env_slave, this->cpx_lp_slave);
	//if (checkCPXstatus(this->cpx_env_slave, status)) goto TERMINATE;

	status = CPXgetx(this->cpx_env_slave, this->cpx_lp_slave, this->slave_curr_sol.data(),
		0, CPXgetnumcols(cpx_env_slave, cpx_lp_slave) - 1);
	if (checkCPXstatus(this->cpx_env_slave, status)) goto TERMINATE;

	status = CPXgetobjval(cpx_env_slave, cpx_lp_slave, obj_val);
	if (checkCPXstatus(this->cpx_env_slave, status)) goto TERMINATE;

	
TERMINATE:
	if (status)
		return false;
	else
		return true;
}

bool RobustKSolver::slave_store_z_values()
{
	bool check = true;
	double int_eps = INT_EPS;
	
	v_dual_set.push_back(dual_vars_t());
	int zz = v_dual_set.size() - 1;

	v_dual_set[zz].z_k.assignMatrix2D(this->problem->getK(), this->problem->getL(), 0);

	for (int k = 0; k < this->problem->getK(); ++k)
	{

		for(int l = 0; l < this->problem->getL(); ++l)
			if (this->slave_curr_sol[vs_z_k(k, l)] > int_eps)
			{
				v_dual_set[zz].z_k(k, l) = 1;
				break;
			}

	}


		return true;
}

bool RobustKSolver::appr_init_cplex_data_str()
{
	int status = 0;

	this->cpx_env_milp = CPXopenCPLEX(&status);
	if (this->checkCPXstatus(this->cpx_env_milp, status)) return false;
	this->cpx_lp_milp = CPXcreateprob(this->cpx_env_milp, &status, this->problem->getProbName().c_str());
	if (checkCPXstatus(this->cpx_env_milp, status))return false;

	return true;
}

bool RobustKSolver::appr_set_cplex_parameters()
{
	bool mod_stat = true;

	int status = 0;

	char errbuf[CPXMESSAGEBUFSIZE];


	status = CPXsetintparam(this->cpx_env_milp, CPX_PARAM_SCRIND, CPX_ON);
	if (checkCPXstatus(this->cpx_env_milp, status)) return false;
	status = CPXsetintparam(this->cpx_env_milp, CPX_PARAM_THREADS, 1);
	if (checkCPXstatus(this->cpx_env_milp, status)) return false;
	status = CPXchgobjsen(this->cpx_env_milp, this->cpx_lp_milp, CPX_MIN);
	if (checkCPXstatus(this->cpx_env_milp, status)) return false;
	status = CPXsetintparam(this->cpx_env_milp, CPX_PARAM_NUMERICALEMPHASIS, CPX_ON);
	if (checkCPXstatus(this->cpx_env_milp, status)) return false;
	status = CPXsetdblparam(this->cpx_env_milp, CPXPARAM_Simplex_Tolerances_Feasibility, 1e-9);
	if (checkCPXstatus(this->cpx_env_milp, status)) return false;

	status = CPXsetdblparam(this->cpx_env_milp, CPXPARAM_Simplex_Tolerances_Optimality, 1e-9);
	if (checkCPXstatus(this->cpx_env_milp, status)) return false;
	status = CPXsetdblparam(this->cpx_env_milp, CPXPARAM_MIP_Tolerances_MIPGap, 1e-9);
	if (checkCPXstatus(this->cpx_env_milp, status)) return false;
	status = CPXsetdblparam(this->cpx_env_milp, CPXPARAM_MIP_Tolerances_Integrality, 0);
	if (checkCPXstatus(this->cpx_env_milp, status)) return false;


	//status = CPXsetdblparam(g_spp_milp.cpx_env, CPXPARAM_TimeLimit, TIME_LIMIT);
	//if (checkStatus(g_spp_milp.cpx_env, status)) goto TERMINATE;


	return true;


}

bool RobustKSolver::appr_build_model()
{
	bool check = true;
	int status;


	/* Data structures*/
	check = this->init_data_structures();
	if (!check)
	{
		goto TERMINATE;
	}

	check = this->appr_init_cplex_data_str();
	if (!check)
	{
		goto TERMINATE;
	}

	check = this->appr_set_cplex_parameters();
	if (!check)
	{
		goto TERMINATE;
	}

	/*
	Variables
	*/
	check = this->appr_add_primal_variables();
	if (!check)
	{
		goto TERMINATE;
	}

	
	check = this->appr_add_dualization_variables();
	if (!check)
	{
			goto TERMINATE;
	}

	check = this->appr_add_bilinear_variables();
	if (!check)
	{
			goto TERMINATE;
	}
	



	/*
	Constraints
	*/
	check = this->appr_add_deterministic_constraints();
	if (!check)
	{
		goto TERMINATE;
	}

	//check = appr_add_primal_valid_inequalities();
	//if (!check)
	//{
	//	goto TERMINATE;
	//}
	////

	
		check = this->appr_add_dualization_constraints();
		if (!check)
		{
			goto TERMINATE;
		}

		check = this->appr_add_linearization_constraints();
		if (!check)
		{
			goto TERMINATE;
		}

		check = this->appr_add_RLT_AlphaXWY_valid_inequalities();
		if (!check)
		{
			goto TERMINATE;
		}

		check = appr_add_RLT_AlphaExwy_valid_inequalities();
		if (!check)
		{
			goto TERMINATE;
		}


		check = appr_add_RLT_Primal_valid_inequalities();
		if (!check)
		{
			goto TERMINATE;
		}

	


	status = CPXwriteprob(this->cpx_env_milp, this->cpx_lp_milp, "approximation.lp", "lp");
	if (this->checkCPXstatus(this->cpx_env_milp, status)) { check = false; goto TERMINATE; };

TERMINATE:
	if (!check)
		return false;
	else
		return true;
}

bool RobustKSolver::appr_solve_problem(double* obj_val)
{
	bool check = true;
	int status;
	int mip_status;

	status = CPXmipopt(this->cpx_env_milp, this->cpx_lp_milp);
	if (checkCPXstatus(this->cpx_env_milp, status)) goto TERMINATE;


	mip_status = CPXgetstat(this->cpx_env_milp, this->cpx_lp_milp);

	status = CPXgetx(this->cpx_env_milp, this->cpx_lp_milp, this->master_curr_sol.data(),
		0, CPXgetnumcols(cpx_env_milp, cpx_lp_milp) - 1);
	if (checkCPXstatus(this->cpx_env_milp, status)) goto TERMINATE;

	
	double best_lb;

	status = CPXgetbestobjval(this->cpx_env_milp, this->cpx_lp_milp, &this->master_curr_blb);
	if (checkCPXstatus(this->cpx_env_milp, status)) goto TERMINATE;

	status = CPXgetobjval(this->cpx_env_milp, this->cpx_lp_milp, obj_val);
	if (checkCPXstatus(this->cpx_env_milp, status)) goto TERMINATE;

	this->master_curr_bub = *obj_val;


TERMINATE:
	if (!check)
		return false;
	else
		return true;
}

bool RobustKSolver::appr_add_primal_variables()
{
	unsigned short num_rows = CPXgetnumcols(this->cpx_env_milp, this->cpx_lp_milp);

	int status = 0;
	int vars_ind = 0;
	double lb[1], ub[1];
	double obj[1];
	char* varsname[1];
	char elname[1024];
	varsname[0] = elname;
	char vartype[1];

	/*
	Phi variable
	*/
	this->v_phi = num_rows;
	obj[0] = 1;
	vartype[0] = 'C';
	sprintf(varsname[0], "Phi");
	lb[0] = -CPX_INFBOUND;
	ub[0] = CPX_INFBOUND;
	status = CPXnewcols(this->cpx_env_milp, this->cpx_lp_milp, 1, obj, lb, ub, vartype, varsname);
	if (checkCPXstatus(this->cpx_env_milp, status)) goto TERMINATE;
	++num_rows;

	/*
	First stage variables
	*/
	if (this->problem->getNx() > 0)
	{
		this->v_x.resize(problem->getNx());

		for (int i = 0; i < this->problem->getNx(); ++i)
		{
			this->v_x[i] = num_rows;

			obj[0] = problem->get1S_dobj_vect_c()[i];

			vartype[0] = 'B';

			sprintf(varsname[0], "x_i%d", i);
			lb[0] = 0;
			ub[0] = 1;

		

			status = CPXnewcols(this->cpx_env_milp, this->cpx_lp_milp, 1, obj, lb, ub, vartype, varsname);
			if (checkCPXstatus(this->cpx_env_milp, status)) goto TERMINATE;
			++num_rows;
		}
	}
	/*
	Discovery variables
	*/
	if (this->problem->getNw() > 0)
	{
		this->v_w.resize(this->problem->getNw());

		for (int i = 0; i < this->problem->getNw(); ++i)
		{
			this->v_w[i] = num_rows;

			obj[0] = this->problem->getDV_dobj_vect_d()[i];

			vartype[0] = 'B';

			sprintf(varsname[0], "w_i%d", i);

			lb[0] = 0;
			ub[0] = 1;

			status = CPXnewcols(this->cpx_env_milp, this->cpx_lp_milp, 1, obj, lb, ub, vartype, varsname);
			if (checkCPXstatus(this->cpx_env_milp, status)) goto TERMINATE;
			++num_rows;
		}
	}
	/*
	Second stage variables.
	*/
	if (this->problem->getNy() > 0)
	{


		this->v_y_k.resizeMatrix2D(this->problem->getK(), this->problem->getNy());

		for (int k = 0; k < this->problem->getK(); ++k)
		{


			for (int i = 0; i < this->problem->getNy(); ++i)
			{


				this->v_y_k(k, i) = num_rows;

				obj[0] = 0;

				vartype[0] = 'B';

				sprintf(varsname[0], "y_k%d_i%d", k, i);
				lb[0] = 0;
				ub[0] = 1;

				


				status = CPXnewcols(this->cpx_env_milp, this->cpx_lp_milp, 1, obj, lb, ub, vartype, varsname);
				if (checkCPXstatus(this->cpx_env_milp, status)) goto TERMINATE;
				++num_rows;
			}
		}
	}

TERMINATE:

	if (status)
	{
		return false;
	}

	return true;
}

bool RobustKSolver::appr_add_dualization_variables()
{
	int status = 0;
	int vars_ind = 0;
	double lb[1], ub[1];
	double obj[1];
	char* varsname[1];
	char elname[1024];
	varsname[0] = elname;
	char vartype[1];


	unsigned short num_cols = CPXgetnumcols(this->cpx_env_milp, this->cpx_lp_milp);

	// Alpha variables: size = K.
	//appr_dual_vars
	appr_dual_vars.v_alpha.resize(this->problem->getK());
	for (unsigned short k = 0; k < this->problem->getK(); ++k)
	{
		appr_dual_vars.v_alpha[k] = num_cols;

		obj[0] = 0;

		vartype[0] = 'C';

		sprintf(varsname[0], "alpha_k%d", k);

		lb[0] = 0;
		ub[0] = 1;


		status = CPXnewcols(this->cpx_env_milp, this->cpx_lp_milp, 1, obj, lb, ub, vartype, varsname);
		if (checkCPXstatus(this->cpx_env_milp, status)) goto TERMINATE;
		++num_cols;
	}

	// Beta variables: size = Lxi {Number of constraints in the Uncertainty set}
	this->appr_dual_vars.v_beta.resize(this->problem->getLxi());
	for (unsigned short r = 0; r < this->problem->getLxi(); ++r)
	{
		this->appr_dual_vars.v_beta[r] = num_cols;

		obj[0] = 0;

		vartype[0] = 'C';

		sprintf(varsname[0], "beta_r%d", r);

		lb[0] = 0;
		ub[0] = CPX_INFBOUND;

		status = CPXnewcols(this->cpx_env_milp, this->cpx_lp_milp, 1, obj, lb, ub, vartype, varsname);
		if (checkCPXstatus(this->cpx_env_milp, status)) goto TERMINATE;
		++num_cols;
	}
	// Beta^k variables. Each Beta^k has a size of Lxi.
	//this->v_dual_set[z].v_beta_k.resize(this->problem->getK());
	this->appr_dual_vars.v_beta_k.resizeMatrix2D(this->problem->getK(), this->problem->getLxi());

	for (unsigned short k = 0; k < this->problem->getK(); ++k)
	{
		// Resize the vector for each k
	

		for (unsigned short r = 0; r < this->problem->getLxi(); ++r)
		{
			

			appr_dual_vars.v_beta_k(k, r) = num_cols;

			obj[0] = 0;

			vartype[0] = 'C';

			sprintf(varsname[0], "betaK_k%d_r%d", k, r);

			lb[0] = 0;
			ub[0] = CPX_INFBOUND;

			status = CPXnewcols(this->cpx_env_milp, this->cpx_lp_milp, 1, obj, lb, ub, vartype, varsname);
			if (checkCPXstatus(this->cpx_env_milp, status)) goto TERMINATE;
			++num_cols;
		}
	}

	// Gamma^k variables
	this->appr_dual_vars.v_gamma_k.resizeMatrix2D(this->problem->getK(), this->problem->getNxi());

	for (unsigned short k = 0; k < this->problem->getK(); ++k)
	{
	

		for (unsigned short r = 0; r < this->problem->getNxi(); ++r)
		{


			appr_dual_vars.v_gamma_k(k, r) = num_cols;

			obj[0] = 0;

			vartype[0] = 'C';

			sprintf(varsname[0], "t_gamma_k%d_r%d", k, r);

			lb[0] = -CPX_INFBOUND;
			ub[0] = CPX_INFBOUND;

			status = CPXnewcols(this->cpx_env_milp, this->cpx_lp_milp, 1, obj, lb, ub, vartype, varsname);
			if (checkCPXstatus(this->cpx_env_milp, status)) goto TERMINATE;
			++num_cols;
		}
	}

	// The other variables
	//rho
	this->appr_dual_vars.v_rho.resize(this->problem->getK());
	for(int k = 0; k < this->problem->getK(); ++k)
	{
		appr_dual_vars.v_rho[k] = num_cols;
		obj[0] = 0;

		vartype[0] = 'C';

		sprintf(varsname[0], "rho_k%d", k);

		lb[0] = 0;
		ub[0] = CPX_INFBOUND;

		status = CPXnewcols(this->cpx_env_milp, this->cpx_lp_milp, 1, obj, lb, ub, vartype, varsname);
		if (checkCPXstatus(this->cpx_env_milp, status)) goto TERMINATE;
		++num_cols;
	}

	this->appr_dual_vars.v_ov_mu_kl.resizeMatrix3D(this->problem->getK(), this->problem->getL(), this->problem->getNxi());
	for (int k = 0; k < this->problem->getK(); ++k)
		for (int l = 0; l < this->problem->getL(); ++l)
			for (int xi = 0; xi < this->problem->getNxi(); ++xi)
			{
				appr_dual_vars.v_ov_mu_kl(k,l,xi) = num_cols;
				obj[0] = 0;

				vartype[0] = 'C';

				sprintf(varsname[0], "ov_mu_k%d_l%d_xi%d", k,l,xi);

				lb[0] = 0;
				ub[0] = CPX_INFBOUND;

				status = CPXnewcols(this->cpx_env_milp, this->cpx_lp_milp, 1, obj, lb, ub, vartype, varsname);
				if (checkCPXstatus(this->cpx_env_milp, status)) goto TERMINATE;
				++num_cols;
			}

	this->appr_dual_vars.v_un_mu_kl.resizeMatrix3D(this->problem->getK(), this->problem->getL(), this->problem->getNxi());
	for (int k = 0; k < this->problem->getK(); ++k)
		for (int l = 0; l < this->problem->getL(); ++l)
			for (int xi = 0; xi < this->problem->getNxi(); ++xi)
			{
				appr_dual_vars.v_un_mu_kl(k, l, xi) = num_cols;
				obj[0] = 0;

				vartype[0] = 'C';

				sprintf(varsname[0], "un_mu_k%d_l%d_xi%d", k, l, xi);

				lb[0] = -CPX_INFBOUND;
				ub[0] = 0;

				status = CPXnewcols(this->cpx_env_milp, this->cpx_lp_milp, 1, obj, lb, ub, vartype, varsname);
				if (checkCPXstatus(this->cpx_env_milp, status)) goto TERMINATE;
				++num_cols;
			}



	this->appr_dual_vars.v_ov_pi_kl.resizeMatrix3D(this->problem->getK(), this->problem->getL(), this->problem->getNxi());
	for (int k = 0; k < this->problem->getK(); ++k)
		for (int l = 0; l < this->problem->getL(); ++l)
			for (int xi = 0; xi < this->problem->getNxi(); ++xi)
			{
				appr_dual_vars.v_ov_pi_kl(k, l, xi) = num_cols;
				
				obj[0] = 0;
					//this->problem->get_lb_xi()[xi];

				vartype[0] = 'C';

				sprintf(varsname[0], "ov_pi_k%d_l%d_xi%d", k, l, xi);

				lb[0] = 0;
				ub[0] = CPX_INFBOUND;

				status = CPXnewcols(this->cpx_env_milp, this->cpx_lp_milp, 1, obj, lb, ub, vartype, varsname);
				if (checkCPXstatus(this->cpx_env_milp, status)) goto TERMINATE;
				++num_cols;
			}


	this->appr_dual_vars.v_un_pi_kl.resizeMatrix3D(this->problem->getK(), this->problem->getL(), this->problem->getNxi());
	for (int k = 0; k < this->problem->getK(); ++k)
		for (int l = 0; l < this->problem->getL(); ++l)
			for (int xi = 0; xi < this->problem->getNxi(); ++xi)
			{
				appr_dual_vars.v_un_pi_kl(k, l, xi) = num_cols;

				obj[0] = 0;
					//this->problem->get_lb_xi()[xi];

				vartype[0] = 'C';

				sprintf(varsname[0], "un_pi_k%d_l%d_xi%d", k, l, xi);

				lb[0] = -CPX_INFBOUND;
				ub[0] = 0;

				status = CPXnewcols(this->cpx_env_milp, this->cpx_lp_milp, 1, obj, lb, ub, vartype, varsname);
				if (checkCPXstatus(this->cpx_env_milp, status)) goto TERMINATE;
				++num_cols;
			}



TERMINATE:
	if (status)
	{
		return false;
	}
	return true;
}

bool RobustKSolver::appr_add_bilinear_variables()
{
	unsigned short num_cols = CPXgetnumcols(this->cpx_env_milp, this->cpx_lp_milp);

	int status = 0;
	int vars_ind = 0;
	double lb[1], ub[1];
	double obj[1];
	char* varsname[1];
	char elname[1024];
	varsname[0] = elname;
	char vartype[1];

	/*
	x times alpha
	*/
	if (this->problem->getNx() > 0
		&&
		this->problem->get1S_uobj_matr_C().getNumElem() > 0
		)
	{
		appr_dual_vars.v_t_x_k.resizeMatrix2D(this->problem->getK(), this->problem->getNx());


		for (int k = 0; k < this->problem->getK(); ++k)
		{

			for (int i = 0; i < this->problem->getNx(); ++i)
			{

				this->appr_dual_vars.v_t_x_k(k, i) = num_cols;
				obj[0] = 0;
				vartype[0] = 'C';

				sprintf(varsname[0], "t_x_k%d_i%d", k, i);
				lb[0] = 0;
				ub[0] = 1;

				status = CPXnewcols(this->cpx_env_milp, this->cpx_lp_milp, 1, obj, lb, ub, vartype, varsname);
				if (checkCPXstatus(this->cpx_env_milp, status)) goto TERMINATE;
				++num_cols;
			}

		}
	}
	/*
	w times alpha
	*/
	if (this->problem->getNw() > 0
		&&
		this->problem->getDV_uobj_matr_D().getNumElem() > 0)
	{

		this->appr_dual_vars.v_t_w_k.resizeMatrix2D(this->problem->getK(), this->problem->getNw());

		for (int k = 0; k < this->problem->getK(); ++k)
		{

			for (int i = 0; i < this->problem->getNw(); ++i)
			{
				

				this->appr_dual_vars.v_t_w_k(k, i) = num_cols;

				obj[0] = 0;

				vartype[0] = 'C';

				sprintf(varsname[0], "t_w_k%d_i%d", k, i);

				lb[0] = 0;
				ub[0] = 1;

				status = CPXnewcols(this->cpx_env_milp, this->cpx_lp_milp, 1, obj, lb, ub, vartype, varsname);
				if (checkCPXstatus(this->cpx_env_milp, status)) goto TERMINATE;
				++num_cols;
			}
		}
	}
	/*
	y^k time alpha_k
	*/
	if (this->problem->getNy() > 0)
	{

		this->appr_dual_vars.v_t_y_k.resizeMatrix2D(this->problem->getK(), this->problem->getNy());

		for (int k = 0; k < this->problem->getK(); ++k)
		{
			for (int i = 0; i < this->problem->getNy(); ++i)
			{


				this->appr_dual_vars.v_t_y_k(k, i) = num_cols;

				obj[0] = 0;

				vartype[0] = 'C';

				sprintf(varsname[0], "t_y_k%d_i%d", k, i);

				lb[0] = 0;
				ub[0] = 1;

				status = CPXnewcols(this->cpx_env_milp, this->cpx_lp_milp, 1, obj, lb, ub, vartype, varsname);
				if (checkCPXstatus(this->cpx_env_milp, status)) goto TERMINATE;
				++num_cols;
			}
		}
	}

	/*gamma_k times w*/

	this->appr_dual_vars.v_t_gamma_k.resizeMatrix2D(this->problem->getK(), this->problem->getNw());

	for (unsigned short k = 0; k < this->problem->getK(); ++k)
	{


		for (unsigned short r = 0; r < this->problem->getNw(); ++r)
		{

			appr_dual_vars.v_t_gamma_k(k, r) = num_cols;

			obj[0] = 0;

			vartype[0] = 'C';

			sprintf(varsname[0], "t_gammaK_z%d_k%d_r%d", k, r);

			lb[0] = -CPX_INFBOUND;
			ub[0] = CPX_INFBOUND;

			status = CPXnewcols(this->cpx_env_milp, this->cpx_lp_milp, 1, obj, lb, ub, vartype, varsname);
			if (checkCPXstatus(this->cpx_env_milp, status)) goto TERMINATE;
			++num_cols;
		}
	}

TERMINATE:

	if (status)
	{
		return false;
	}

	return true;
}

bool RobustKSolver::appr_add_deterministic_constraints()
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
	status = 0;
	double coeff_eps = COEFF_EPS;

	std::vector<double> b_loc = this->problem->get1S_drhs_vect_bx();
	Matrix2D<double> Matr_loc = this->problem->get1S_dcon_matr_X();



	// X matrix constraints
	for (unsigned short row = 0; row < this->problem->getLx(); ++row)
	{
		rhs[0] = b_loc[row];
		sense[0] = 'L';
		nzc = 0;
		sprintf(cnstrname[0], "Xx_%d", row);
		// Now we go through the "columns", coefficients of the matrices 

		for (unsigned short cols = 0; cols < this->problem->getNx(); ++cols)
			if (fabs(Matr_loc(row, cols)) > coeff_eps)
			{
				this->matind[nzc] = this->v_x[cols];
				this->matval[nzc] = Matr_loc(row, cols);
				++nzc;
			}
		status = CPXaddrows(this->cpx_env_milp, this->cpx_lp_milp, 0, 1, nzc, rhs, sense, matbeg, this->matind.data(), this->matval.data(), NULL, cnstrname);
		if (checkCPXstatus(this->cpx_env_milp, status)) goto TERMINATE;

	}

	// W matrices constraints 
	b_loc = this->problem->getDV_drhs_vect_bw();
	Matr_loc = this->problem->getDV_dcon_matr_W();

	for (unsigned short row = 0; row < this->problem->getLw(); ++row)
	{
		rhs[0] = b_loc[row];
		sense[0] = 'L';
		nzc = 0;
		sprintf(cnstrname[0], "Ww_%d", row);
		// Now we go trought the "columns", coefficient of the matrices 

		for (unsigned short cols = 0; cols < this->problem->getNw(); ++cols)
			if (fabs(Matr_loc(row, cols)) > coeff_eps)
			{
				this->matind[nzc] = this->v_w[cols];
				this->matval[nzc] = Matr_loc(row, cols);
				++nzc;
			}
		status = CPXaddrows(this->cpx_env_milp, this->cpx_lp_milp, 0, 1, nzc, rhs, sense, matbeg, this->matind.data(), this->matval.data(), NULL, cnstrname);
		if (checkCPXstatus(this->cpx_env_milp, status)) goto TERMINATE;

	}


	// Y matrices constraints for each k

	b_loc = this->problem->get2S_drhs_vect_by();
	Matr_loc = this->problem->get2S_dcon_matr_Y();

	for (int k = 0; k < this->problem->getK(); ++k)
	{
		for (unsigned short row = 0; row < this->problem->getLy(); ++row)
		{
			rhs[0] = b_loc[row];
			sense[0] = 'L';
			nzc = 0;
			sprintf(cnstrname[0], "Yy_k%d_%d", k, row);
			// Now we go trought the "columns", coefficient of the matrices 

			for (unsigned short cols = 0; cols < this->problem->getNy(); ++cols)
				if (fabs(Matr_loc(row, cols)) > coeff_eps)
				{
					this->matind[nzc] = this->v_y_k(k, cols);
					this->matval[nzc] = Matr_loc(row, cols);
					++nzc;
				}
			status = CPXaddrows(this->cpx_env_milp, this->cpx_lp_milp, 0, 1, nzc, rhs, sense, matbeg, this->matind.data(), this->matval.data(), NULL, cnstrname);
			if (checkCPXstatus(this->cpx_env_milp, status)) goto TERMINATE;

		}
	}
	// Ex(x) + Ew(w) + Ey(y) <=g
	if (this->problem->getR() > 0)
	{
		b_loc = this->problem->getLink_dcon_vect_g();


		for (int row = 0; row < this->problem->getR(); ++row)
			for (int k = 0; k < this->problem->getK(); ++k)
			{

				sense[0] = 'L';
				nzc = 0;
				rhs[0] = b_loc[row];

				sprintf(cnstrname[0], "LinkXWY_k%d_r%d", k, row);

				Matr_loc = this->problem->getLink_dcon_matr_Ex();
				if (Matr_loc.getNumElem() > 0)
				{
					for (int col = 0; col < this->problem->getNx(); ++col)
						if (fabs(Matr_loc(row, col)) > coeff_eps)
						{
							this->matind[nzc] = this->v_x[col];
							this->matval[nzc] = Matr_loc(row, col);
							++nzc;
						}
				}
				Matr_loc = this->problem->getLink_dcon_matr_Ew();
				if (Matr_loc.getNumElem() > 0)
				{
					for (int col = 0; col < this->problem->getNw(); ++col)
						if (fabs(Matr_loc(row, col)) > coeff_eps)
						{
							this->matind[nzc] = this->v_w[col];
							this->matval[nzc] = Matr_loc(row, col);
							++nzc;
						}
				}

				Matr_loc = this->problem->getLink_dcon_matr_Ey();
				if (Matr_loc.getNumElem() > 0)
				{
					for (int col = 0; col < this->problem->getNy(); ++col)
						if (fabs(Matr_loc(row, col)) > coeff_eps)
						{
							this->matind[nzc] = this->v_y_k(k, col);
							this->matval[nzc] = Matr_loc(row, col);
							++nzc;
						}
				}

				status = CPXaddrows(this->cpx_env_milp, this->cpx_lp_milp, 0, 1, nzc, rhs, sense, matbeg, this->matind.data(), this->matval.data(), NULL, cnstrname);
				if (checkCPXstatus(this->cpx_env_milp, status)) goto TERMINATE;
			}
	}


TERMINATE:
	if (status)
	{
		return false;
	}

	return true;
}

bool RobustKSolver::appr_add_primal_valid_inequalities()
{
	/*
	For ^k
	*/
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
	status = 0;
	double coeff_eps = COEFF_EPS;

	std::vector<double> y_coeff_loc = this->problem->get2S_dobj_vect_lbq();
	std::vector<double> x_coeff_loc = this->problem->get1S_dobj_vect_lbc();
	std::vector<double> c_loc = this->problem->get1S_dobj_vect_c();

	for (int k = 0; k < this->problem->getK(); ++k)
	{
		rhs[0] = 0;
		sense[0] = 'G';
		nzc = 0;

		sprintf(cnstrname[0], "Det_Obj_cut_k%d", k);

		this->matind[nzc] = this->v_phi;
		this->matval[nzc] = 1;
		++nzc;


		for (unsigned short cols = 0; cols < this->problem->getNx(); ++cols)
			if (fabs(x_coeff_loc[cols]) > coeff_eps)
			{
				this->matind[nzc] = this->v_x[cols];
				this->matval[nzc] = -x_coeff_loc[cols];
				++nzc;
			}



		for (unsigned short cols = 0; cols < this->problem->getNy(); ++cols)
			if (fabs(y_coeff_loc[cols]) > coeff_eps)
			{
				this->matind[nzc] = this->v_y_k(k, cols);
				this->matval[nzc] = -y_coeff_loc[cols];
				++nzc;
			}

		status = CPXaddrows(this->cpx_env_milp, this->cpx_lp_milp, 0, 1, nzc, rhs, sense, matbeg, this->matind.data(), this->matval.data(), NULL, cnstrname);
		if (checkCPXstatus(this->cpx_env_milp, status)) goto TERMINATE;
	}


	/*
	Here we add the approx constr cut.
	*/
	y_coeff_loc = this->problem->get2S_ducon_vect_p();
	x_coeff_loc = this->problem->get1S_ducon_vect_t();

	for (int k = 0; k < this->problem->getK(); ++k)
		for (int l = 0; l < this->problem->getL(); ++l)
		{
			rhs[0] = this->problem->get_rhs_ducon_vect_hmax()[l];
			sense[0] = 'L';
			nzc = 0;

			sprintf(cnstrname[0], "Det_CU_cut_k%d_l%d", k,l);


			for (unsigned short col = 0; col < this->problem->getNx(); ++col)
				if (fabs(y_coeff_loc[col]) > coeff_eps)
				{
					this->matind[nzc] = this->v_x[col];
					this->matval[nzc] = x_coeff_loc[col];
					++nzc;
				}

			for (unsigned short col = 0; col < this->problem->getNy(); ++col)
				if (fabs(y_coeff_loc[col]) > coeff_eps)
				{
					this->matind[nzc] = this->v_y_k(k, col);
					this->matval[nzc] = y_coeff_loc[col];
					++nzc;
				}

			status = CPXaddrows(this->cpx_env_milp, this->cpx_lp_milp, 0, 1, nzc, rhs, sense, matbeg, this->matind.data(), this->matval.data(), NULL, cnstrname);
			if (checkCPXstatus(this->cpx_env_milp, status)) goto TERMINATE;

		}

TERMINATE:
	if (status)
	{
		return false;
	}

	return true;

}

bool RobustKSolver::appr_add_dualization_constraints()
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

	double coeff_eps = COEFF_EPS;
	double int_eps = INT_EPS;
	double big_m = CNSTR_UNC_BIG_M;

	Matrix2D<double> A_loc = this->problem->getUC_xicon_matr_A();
	Matrix2D<double> Q_loc = this->problem->get2S_uobj_matr_Q();
	Matrix2D<double> C_loc = this->problem->get1S_uobj_matr_C();
	Matrix2D<double> D_loc = this->problem->getDV_uobj_matr_D();

	Matrix2D<double> M_loc;
	Matrix3D<double> M3_loc;

	std::vector<double> b_loc;
	
	double coeff = 0;

	/*\phi >= uncertainty_cost*/
	sense[0] = 'G';
	rhs[0] = 0;
	nzc = 0;
	sprintf(cnstrname[0], "objective_constr");

	// \phi
	this->matind[nzc] = this->v_phi;
	this->matval[nzc] = 1;
	++nzc;
	// \beta and \beta_k
	for (int bb = 0; bb < this->problem->getLxi(); ++bb)
	{
		this->matind[nzc] = this->appr_dual_vars.v_beta[bb];
		this->matval[nzc] = -this->problem->getUC_rhs_vect_b()[bb];
		++nzc;

		if (nzc >= MAX_NUM_VARS)
		{
			getchar();
		}

		for (int k = 0; k < this->problem->getK(); ++k)
		{
			this->matind[nzc] = this->appr_dual_vars.v_beta_k(k, bb);
			this->matval[nzc] = -this->problem->getUC_rhs_vect_b()[bb];
			++nzc;

			if (nzc >= MAX_NUM_VARS)
			{
				getchar();
			}
		}
	}
	

	
	// tilde_y
	M_loc = this->problem->get2S_ucon_matr_P();
	for (int k = 0; k < this->problem->getK(); ++k)
	{
		for (int yi = 0; yi < this->problem->getNy(); ++yi)
		{

			coeff = -this->problem->get2S_dobj_vect_q()[yi];
			
			if (fabs(coeff) > coeff_eps)
			{
				this->matind[nzc] = this->appr_dual_vars.v_t_y_k(k, yi);
				this->matval[nzc] = coeff;
				++nzc;
			}

			if (nzc >= MAX_NUM_VARS)
			{
				getchar();
			}
		}

	}

	// rho
	for (int k = 0; k < this->problem->getK(); ++k)
	{
		this->matind[nzc] = this->appr_dual_vars.v_rho[k];
		this->matval[nzc] = -1;
		++nzc;

		if (nzc >= MAX_NUM_VARS)
		{
			getchar();
		}
	}


	// ov_pi and un_pi
	for (int k = 0; k < this->problem->getK(); ++k)
		for(int l = 0; l < this->problem->getL(); ++l)
			for (int xi = 0; xi < this->problem->getNxi(); ++xi)
			{
				coeff = this->problem->get_lb_xi()[xi];

				if (fabs(coeff) > coeff_eps)
				{
					this->matind[nzc] = this->appr_dual_vars.v_ov_pi_kl(k, l, xi);
					this->matval[nzc] = coeff;
					++nzc;
				}

				if (nzc >= MAX_NUM_VARS)
				{
					getchar();
				}

				coeff = this->problem->get_ub_xi()[xi];

				if (fabs(coeff) > coeff_eps)
				{
					this->matind[nzc] = this->appr_dual_vars.v_un_pi_kl(k, l, xi);
					this->matval[nzc] = coeff;
					++nzc;
				}

				if (nzc >= MAX_NUM_VARS)
				{
					getchar();
				}


			}


	status = CPXaddrows(this->cpx_env_milp, this->cpx_lp_milp, 0, 1, nzc, rhs, sense, matbeg, this->matind.data(), this->matval.data(), NULL, cnstrname);
	if (checkCPXstatus(this->cpx_env_milp, status)) goto TERMINATE;



	/*Sum over alpha = 1*/
	sense[0] = 'E';
	rhs[0] = 1;
	nzc = 0;

	sprintf(cnstrname[0], "Sum_k_Alpha = 1");

	for (int k = 0; k < this->problem->getK(); ++k)
	{
		this->matind[nzc] = this->appr_dual_vars.v_alpha[k];
		this->matval[nzc] = 1;
		++nzc;

	}
	status = CPXaddrows(this->cpx_env_milp, this->cpx_lp_milp, 0, 1, nzc, rhs, sense, matbeg, this->matind.data(), this->matval.data(), NULL, cnstrname);
	if (checkCPXstatus(this->cpx_env_milp, status)) goto TERMINATE;


	/*At Beta + \sum_k w \circ gamma_k* , \foraeach Nxi (one constr for each component of the \xi vector  */
	/*
	k is the policy index.
	row is the constraints end (equivalenlty) the component of vector \xi
	*/

	for (int row = 0; row < this->problem->getNxi(); ++row)
	{
		rhs[0] = 0;
		sense[0] = 'E';
		sprintf(cnstrname[0], "Balance_n%d", row);
		nzc = 0;

		if (A_loc.getNumElem() > 0)
		{
			// A is transpose. THerefore row and col are inverted.
			for (int col = 0; col < this->problem->getLxi(); ++col)
				if (fabs(
					A_loc(col, row)
					//this->problem->getUC_xicon_matr_A()(col,row)
				)
			> coeff_eps)
				{
					this->matind[nzc] = this->appr_dual_vars.v_beta[col];
					this->matval[nzc] = A_loc(col, row);
					++nzc;
				}
		}

		for (int k = 0; k < this->problem->getK(); ++k)
		{
			this->matind[nzc] = this->appr_dual_vars.v_t_gamma_k(k, row);
			this->matval[nzc] = 1;
			++nzc;
		}


		status = CPXaddrows(this->cpx_env_milp, this->cpx_lp_milp, 0, 1, nzc, rhs, sense, matbeg, this->matind.data(), this->matval.data(), NULL, cnstrname);
		if (checkCPXstatus(this->cpx_env_milp, status)) goto TERMINATE;

	}



	for (int k = 0; k < this->problem->getK(); ++k)
		for (int row = 0; row < this->problem->getNxi(); ++row)
		{
			rhs[0] = 0;
			sense[0] = 'E';
			sprintf(cnstrname[0], "BalanceK_k%d_r%d", k, row);
			nzc = 0;

			/*
			At Beta^k
			*/
			if (A_loc.getNumElem() > 0)
			{
				// A is transpose. Therefore row and col are inverted.

				for (int col = 0; col < this->problem->getLxi(); ++col)
					if (
						fabs(
							//this->problem->getUC_xicon_matr_A()(col,row)
							A_loc(col, row)
						) > coeff_eps)
					{
						this->matind[nzc] = this->appr_dual_vars.v_beta_k(k, col);
						this->matval[nzc] = A_loc(col, row);
						++nzc;
					}
			}

			/*
			 -(w \circ gamma)
			*/
			this->matind[nzc] = this->appr_dual_vars.v_t_gamma_k(k, row);
			this->matval[nzc] = -1;
			++nzc;


			//Cx \times alpha_k + T_n x \alpha 
			M3_loc = this->problem->get1S_ucon_matr_Txi();
			for (int col = 0; col < this->problem->getNx(); ++col)
			{
				coeff = 0;

				if (fabs(C_loc(row, col)) > coeff_eps)
				{

					coeff = -C_loc(row, col);
					//++nzc;

				}

		/*		for (int l = 0; l < this->problem->getL(); ++l)
					if (this->appr_dual_vars.z_k(k, l) > 0
						&&
						fabs(M3_loc(row, l, col)) > coeff_eps)
					{
						coeff += -(M3_loc(row, l, col) * big_m);
					}*/

				if (fabs(coeff) > coeff_eps)
				{
					this->matind[nzc] = this->appr_dual_vars.v_t_x_k(k, col);
					this->matval[nzc] = coeff;
					++nzc;

				}
			}

			//Dw\times alpha_k + V_n w \alpha 
			M3_loc = this->problem->getDV_ucon_matr_Vxi();
			for (int col = 0; col < this->problem->getNw(); ++col)
			{
				coeff = 0;

				if (fabs(D_loc(row, col)) > coeff_eps)
				{

					coeff = -D_loc(row, col);
					//++nzc;

				}
			/*	for (int l = 0; l < this->problem->getL(); ++l)
					if (this->appr_dual_vars.z_k(k, l) > coeff_eps
						&&
						fabs(M3_loc(row, l, col)) > coeff_eps)
					{
						coeff += -(M3_loc(row, l, col) * big_m);
					}*/

				if (fabs(coeff) > coeff_eps)
				{
					this->matind[nzc] = this->appr_dual_vars.v_t_w_k(k, col);
					this->matval[nzc] = coeff;
					++nzc;
				}
			}

			//Qy\times alpha_k + P_n w \alpha 
			M3_loc = this->problem->get2S_ucon_matr_Pxi();
			for (int col = 0; col < this->problem->getNy(); ++col)
			{
				coeff = 0;

				if (fabs(Q_loc(row, col)) > coeff_eps)
				{

					coeff = -Q_loc(row, col);

				}

				/*for (int l = 0; l < this->problem->getL(); ++l)
					if (this->appr_dual_vars.z_k(k, l) > 0
						&&
						fabs(M3_loc(row, l, col)) > coeff_eps)
					{
						coeff += -(M3_loc(row, l, col) * big_m);
					}*/

				if (fabs(coeff) > coeff_eps)
				{
					this->matind[nzc] = this->appr_dual_vars.v_t_y_k(k, col);
					this->matval[nzc] = coeff;
					++nzc;
				}
			}

			/*alpha_k times H */
	/*		M_loc = this->problem->getUC_ucon_matr_H();
			coeff = 0;
			for (int l = 0; l < this->problem->getL(); ++l)
				if (this->appr_dual_vars.z_k(k, l) > 0
					&&
					fabs(M_loc(l, row)) > coeff_eps)
				{
					coeff += M_loc(l, row) * big_m;
				}
			if (fabs(coeff) > coeff_eps)
			{
				this->matind[nzc] = this->appr_dual_vars.v_alpha[k];
				this->matval[nzc] = coeff;
				++nzc;
			}*/

			for (int l = 0; l < this->problem->getL(); ++l)
			{
				this->matind[nzc] = this->appr_dual_vars.v_ov_pi_kl(k, l, row);
				this->matval[nzc] = -1;
				++nzc;

				this->matind[nzc] = this->appr_dual_vars.v_un_pi_kl(k, l, row);
				this->matval[nzc] = -1;
				++nzc;
			}

		
			status = CPXaddrows(this->cpx_env_milp, this->cpx_lp_milp, 0, 1, nzc, rhs, sense, matbeg, this->matind.data(), this->matval.data(), NULL, cnstrname);
			if (checkCPXstatus(this->cpx_env_milp, status)) goto TERMINATE;
			
		}

		// Dual for z^k_l variables

		b_loc = this->problem->get_rhs_ucon_vect_h();
		for (int k = 0; k < this->problem->getK(); ++k)
			for (int l = 0; l < this->problem->getL(); ++l)
			{
				rhs[0] = 0;
				sense[0] = 'L';
				sprintf(cnstrname[0], "BalanceZvar_k%d_l%d", k, l);
				nzc = 0;

				/*rho_k*/
				this->matind[nzc] = this->appr_dual_vars.v_rho[k];
				this->matval[nzc] = 1;
				++nzc;

				/*alpha*/
				if (fabs(b_loc[l]) > coeff_eps)
				{
					this->matind[nzc] = this->appr_dual_vars.v_alpha[k];
					this->matval[nzc] = (b_loc[l]*big_m);
					++nzc;
				}


				/*t_x_k*/
				M_loc = this->problem->get1S_ucon_matr_T();
				for (int x = 0; x < this->problem->getNx(); ++x)
				{
					coeff = 0;
					

					coeff -= M_loc(l, x);

					/*for (int l = 0; l < this->problem->getL(); ++l)
						if (
							fabs(M_loc(l, x)) > coeff_eps)
						{
							
						}*/

					if (fabs(coeff) > coeff_eps)
					{
						this->matind[nzc] = this->appr_dual_vars.v_t_x_k(k, x);
						this->matval[nzc] = coeff * big_m;
						++nzc;
					}
					if (nzc >= MAX_NUM_VARS)
					{
						getchar();
					}
				}
				/*t_w_k*/
				M_loc = this->problem->getDV_ucon_matr_V();
				for (int w = 0; w < this->problem->getNw(); ++w)
				{
					coeff = 0;


				/*	for (int l = 0; l < this->problem->getL(); ++l)
						if (
							fabs(M_loc(l, w)) > coeff_eps)
						{
							
						}*/

					coeff -= M_loc(l, w);

					if (fabs(coeff) > coeff_eps)
					{
						this->matind[nzc] = this->appr_dual_vars.v_t_w_k(k, w);
						this->matval[nzc] = coeff * big_m;
						++nzc;
					}
					if (nzc >= MAX_NUM_VARS)
					{
						getchar();
					}
				}
				/*t_y_k*/
				M_loc = this->problem->get2S_ucon_matr_P();
				for (int y = 0; y < this->problem->getNy(); ++y)
				{
					coeff = 0;


					/*for (int l = 0; l < this->problem->getL(); ++l)
						if (
							fabs(M_loc(l, y)) > coeff_eps)
						{
							
						}*/
					coeff -= M_loc(l, y);
					if (fabs(coeff) > coeff_eps)
					{
						this->matind[nzc] = this->appr_dual_vars.v_t_y_k(k, y);
						this->matval[nzc] = coeff * big_m;
						++nzc;
					}
					if (nzc >= MAX_NUM_VARS)
					{
						getchar();
					}
				}


				/*ov_mu,  un_pi*/
				b_loc = this->problem->get_ub_xi();
				for (int xi = 0; xi < this->problem->getNxi(); ++xi)
					if(fabs(b_loc[xi] > coeff_eps))
					{
					
						this->matind[nzc] = this->appr_dual_vars.v_ov_mu_kl(k, l, xi);
						this->matval[nzc] = - b_loc[xi];
						++nzc;

						this->matind[nzc] = this->appr_dual_vars.v_un_pi_kl(k, l, xi);
						this->matval[nzc] = -b_loc[xi];
						++nzc;
					}

				/*un_mu, ov_pi,*/
				b_loc = this->problem->get_lb_xi();
				for (int xi = 0; xi < this->problem->getNxi(); ++xi)
					if (fabs(b_loc[xi] > coeff_eps))
					{

						this->matind[nzc] = this->appr_dual_vars.v_un_mu_kl(k, l, xi);
						this->matval[nzc] = -b_loc[xi];
						++nzc;

						this->matind[nzc] = this->appr_dual_vars.v_ov_pi_kl(k, l, xi);
						this->matval[nzc] = -b_loc[xi];
						++nzc;
					}



				status = CPXaddrows(this->cpx_env_milp, this->cpx_lp_milp, 0, 1, nzc, rhs, sense, matbeg, this->matind.data(), this->matval.data(), NULL, cnstrname);
				if (checkCPXstatus(this->cpx_env_milp, status)) goto TERMINATE;
			}

			/*Chi Balance constraints*/
			M_loc = this->problem->getUC_ucon_matr_H();
			for (int k = 0; k < this->problem->getK(); ++k)
				for (int l = 0; l < this->problem->getL(); ++l)
					for (int xi = 0; xi < this->problem->getNxi(); ++xi)
					{
						rhs[0] = 0;
						sense[0] = 'E';
						sprintf(cnstrname[0], "BalanceZvar_k%d_l%d", k, l);
						nzc = 0;

						/*alpha_k times h */
		if (fabs(M_loc(l,xi)) > coeff_eps)
		{
			this->matind[nzc] = this->appr_dual_vars.v_alpha[k];
			this->matval[nzc] = (M_loc(l, xi)*big_m);
			++nzc;
		}

		M3_loc = this->problem->get1S_ucon_matr_Txi();
		for (int x = 0; x < this->problem->getNx(); ++x)
			if(fabs(M3_loc(xi,l,x)) > coeff_eps)
			{

				this->matind[nzc] = this->appr_dual_vars.v_t_x_k(k, x);
				this->matval[nzc] = -(M3_loc(xi, l, x) * big_m);
				++nzc;

			}

		M3_loc = this->problem->getDV_ucon_matr_Vxi();
		for (int w = 0; w < this->problem->getNw(); ++w)
			if (fabs(M3_loc(xi, l, w)) > coeff_eps)
			{

				this->matind[nzc] = this->appr_dual_vars.v_t_w_k(k, w);
				this->matval[nzc] = -(M3_loc(xi, l, w) * big_m);
				++nzc;
			}

		M3_loc = this->problem->get2S_ucon_matr_Pxi();
		for (int y = 0; y < this->problem->getNy(); ++y)
			if (fabs(M3_loc(xi, l, y)) > coeff_eps)
			{

				this->matind[nzc] = this->appr_dual_vars.v_t_y_k(k, y);
				this->matval[nzc] = -(M3_loc(xi, l, y) * big_m);
				++nzc;
			}


		this->matind[nzc] = this->appr_dual_vars.v_ov_mu_kl(k, l, xi);
		this->matval[nzc] = 1;
		++nzc;

		this->matind[nzc] = this->appr_dual_vars.v_un_mu_kl(k, l, xi);
		this->matval[nzc] = 1;
		++nzc;

		this->matind[nzc] = this->appr_dual_vars.v_ov_pi_kl(k, l, xi);
		this->matval[nzc] = 1;
		++nzc;

		this->matind[nzc] = this->appr_dual_vars.v_un_pi_kl(k, l, xi);
		this->matval[nzc] = 1;
		++nzc;



		status = CPXaddrows(this->cpx_env_milp, this->cpx_lp_milp, 0, 1, nzc, rhs, sense, matbeg, this->matind.data(), this->matval.data(), NULL, cnstrname);
		if (checkCPXstatus(this->cpx_env_milp, status)) goto TERMINATE;

		}




	/*
	* Constraints for \zeta and \zeta^k
	*/
	if (this->problem->getNzeta() > 0 &&
		this->problem->getUC_zetacon_matr_G().getNumElem() > 0
		)
	{

		/*
		\zeta
		*/
		rhs[0] = 0;
		sense[0] = 'E';
		nzc = 0;

		for (int row = 0; row < this->problem->getNzeta(); ++row)
		{

			sprintf(cnstrname[0], "BalanceZeta_r%d", row);
			nzc = 0;

			for (int col = 0; col < this->problem->getNzeta(); ++col)
				if (fabs(this->problem->getUC_zetacon_matr_G()(col, row)) > coeff_eps)
				{
					this->matind[nzc] = this->appr_dual_vars.v_beta[col];
					this->matval[nzc] = this->problem->getUC_zetacon_matr_G()(col, row);
					++nzc;
				}
			status = CPXaddrows(this->cpx_env_milp, this->cpx_lp_milp, 0, 1, nzc, rhs, sense, matbeg, this->matind.data(), this->matval.data(), NULL, cnstrname);
			if (checkCPXstatus(this->cpx_env_milp, status)) goto TERMINATE;
		}
		/*
		\zeta^k
		*/
		rhs[0] = 0;
		sense[0] = 'E';
		nzc = 0;

		for (int k = 0; k < this->problem->getK(); ++k)
			for (int row = 0; row < this->problem->getNzeta(); ++row)
			{

				sprintf(cnstrname[0], "BalanceZetaK_k%d_r%d", k, row);
				nzc = 0;

				for (int col = 0; col < this->problem->getNzeta(); ++col)
					if (fabs(this->problem->getUC_zetacon_matr_G()(col, row)) > coeff_eps)
					{
						this->matind[nzc] = this->appr_dual_vars.v_beta_k(k, row);
						this->matval[nzc] = this->problem->getUC_zetacon_matr_G()(col, row);
						++nzc;
					}
				status = CPXaddrows(this->cpx_env_milp, this->cpx_lp_milp, 0, 1, nzc, rhs, sense, matbeg, this->matind.data(), this->matval.data(), NULL, cnstrname);
				if (checkCPXstatus(this->cpx_env_milp, status)) goto TERMINATE;
			}


	}


TERMINATE:
	if (status)
	{
		return false;
	}

	return true;
}




bool RobustKSolver::appr_add_linearization_constraints()
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

	double coeff_eps = COEFF_EPS;
	double big_m = BIG_M;
	status = 0;
	/*
	 x * alpha_k
	*/
	if (this->problem->getNx() > 0
		&&
		this->problem->get1S_uobj_matr_C().getNumElem() > 0)
	{
		for (int k = 0; k < this->problem->getK(); ++k)
			for (int xi = 0; xi < this->problem->getNx(); ++xi)
			{
				/*
					* t_x_k <= x
					*/
				rhs[0] = 0;
				sense[0] = 'L';
				sprintf(cnstrname[0], "Linear_t_x_(1)_k%d_i%d", k, xi);
				nzc = 0;
				
				this->matind[nzc] = this->appr_dual_vars.v_t_x_k(k, xi);
				this->matval[nzc] = 1;
				++nzc;

				this->matind[nzc] = this->v_x[xi];
				this->matval[nzc] = -1;
				++nzc;

				status = CPXaddrows(this->cpx_env_milp, this->cpx_lp_milp, 0, 1, nzc, rhs, sense, matbeg, this->matind.data(), this->matval.data(), NULL, cnstrname);
				if (checkCPXstatus(this->cpx_env_milp, status)) goto TERMINATE;

				/*
				* t_x_k <= alpha_k
				*/
				rhs[0] = 0;
				sense[0] = 'L';
				sprintf(cnstrname[0], "Linear_t_x_(2)_k%d_i%d", k, xi);
				nzc = 0;

				this->matind[nzc] = this->appr_dual_vars.v_t_x_k(k, xi);
				this->matval[nzc] = 1;
				++nzc;

				this->matind[nzc] = this->appr_dual_vars.v_alpha[k];
				this->matval[nzc] = -1;
				++nzc;

				status = CPXaddrows(this->cpx_env_milp, this->cpx_lp_milp, 0, 1, nzc, rhs, sense, matbeg, this->matind.data(), this->matval.data(), NULL, cnstrname);
				if (checkCPXstatus(this->cpx_env_milp, status)) goto TERMINATE;

				/*
				t_x_k - alpha_k - x >= -1
				*/
				rhs[0] = -1;
				sense[0] = 'G';
				sprintf(cnstrname[0], "Linear_t_x_(3)_k%d_i%d", k, xi);
				nzc = 0;

				this->matind[nzc] = this->appr_dual_vars.v_t_x_k(k, xi);
				this->matval[nzc] = 1;
				++nzc;

				this->matind[nzc] = this->appr_dual_vars.v_alpha[k];
				this->matval[nzc] = -1;
				++nzc;

				this->matind[nzc] = this->v_x[xi];
				this->matval[nzc] = -1;
				++nzc;

				status = CPXaddrows(this->cpx_env_milp, this->cpx_lp_milp, 0, 1, nzc, rhs, sense, matbeg, this->matind.data(), this->matval.data(), NULL, cnstrname);
				if (checkCPXstatus(this->cpx_env_milp, status)) goto TERMINATE;
			}
	}

	/*
	* w * alpha_k
	*
	*/
	if (this->problem->getNw() > 0
		&&
		this->problem->getDV_uobj_matr_D().getNumElem() > 0)
	{
		for (int k = 0; k < this->problem->getK(); ++k)
			for (int wi = 0; wi < this->problem->getNw(); ++wi)
			{
				/*
					* t_w_k <= w
					*/
				rhs[0] = 0;
				sense[0] = 'L';
				sprintf(cnstrname[0], "Linear_t_w_(1)_k%d_i%d", k, wi);
				nzc = 0;

				this->matind[nzc] = this->appr_dual_vars.v_t_w_k(k, wi);
				this->matval[nzc] = 1;
				++nzc;

				this->matind[nzc] = this->v_w[wi];
				this->matval[nzc] = -1;
				++nzc;

				status = CPXaddrows(this->cpx_env_milp, this->cpx_lp_milp, 0, 1, nzc, rhs, sense, matbeg, this->matind.data(), this->matval.data(), NULL, cnstrname);
				if (checkCPXstatus(this->cpx_env_milp, status)) goto TERMINATE;

				/*
				* t_w_k <= alpha_k
				*/
				rhs[0] = 0;
				sense[0] = 'L';
				sprintf(cnstrname[0], "Linear_t_w_(2)_k%d_i%d", k, wi);
				nzc = 0;

				this->matind[nzc] = this->appr_dual_vars.v_t_w_k(k, wi);
				this->matval[nzc] = 1;
				++nzc;

				this->matind[nzc] = this->appr_dual_vars.v_alpha[k];
				this->matval[nzc] = -1;
				++nzc;

				status = CPXaddrows(this->cpx_env_milp, this->cpx_lp_milp, 0, 1, nzc, rhs, sense, matbeg, this->matind.data(), this->matval.data(), NULL, cnstrname);
				if (checkCPXstatus(this->cpx_env_milp, status)) goto TERMINATE;

				/*
				t_w_k - alpha_k - w(UB alpha) >= -1(UN alpha)
				*/
				rhs[0] = -1;
				sense[0] = 'G';
				sprintf(cnstrname[0], "Linear_t_x_(3)_k%d_i%d", k, wi);
				nzc = 0;

				this->matind[nzc] = this->appr_dual_vars.v_t_w_k(k, wi);
				this->matval[nzc] = 1;
				++nzc;

				this->matind[nzc] = this->appr_dual_vars.v_alpha[k];
				this->matval[nzc] = -1;
				++nzc;

				this->matind[nzc] = this->v_w[wi];
				this->matval[nzc] = -1;
				++nzc;

				status = CPXaddrows(this->cpx_env_milp, this->cpx_lp_milp, 0, 1, nzc, rhs, sense, matbeg, this->matind.data(), this->matval.data(), NULL, cnstrname);
				if (checkCPXstatus(this->cpx_env_milp, status)) goto TERMINATE;
			}
	}



	/*
	* y^k * alpha_k
	*/
	for (int k = 0; k < this->problem->getK(); ++k)
		for (int yi = 0; yi < this->problem->getNy(); ++yi)
		{
			/*
			* t_y_k <= y_k
			*/
			rhs[0] = 0;
			sense[0] = 'L';
			sprintf(cnstrname[0], "Linear_t_y_(1)_k%d_i%d", k, yi);
			nzc = 0;

			this->matind[nzc] = this->appr_dual_vars.v_t_y_k(k, yi);
			this->matval[nzc] = 1;
			++nzc;

			this->matind[nzc] = this->v_y_k(k, yi);
			this->matval[nzc] = -1;
			++nzc;

			status = CPXaddrows(this->cpx_env_milp, this->cpx_lp_milp, 0, 1, nzc, rhs, sense, matbeg, this->matind.data(), this->matval.data(), NULL, cnstrname);
			if (checkCPXstatus(this->cpx_env_milp, status)) goto TERMINATE;

			/*
			* t_y_k <= alpha_k
			*/
			rhs[0] = 0;
			sense[0] = 'L';
			sprintf(cnstrname[0], "Linear_t_y_(2)_k%d_i%d", k, yi);
			nzc = 0;

			this->matind[nzc] = this->appr_dual_vars.v_t_y_k(k, yi);
			this->matval[nzc] = 1;
			++nzc;

			this->matind[nzc] = this->appr_dual_vars.v_alpha[k];
			this->matval[nzc] = -1;
			++nzc;

			status = CPXaddrows(this->cpx_env_milp, this->cpx_lp_milp, 0, 1, nzc, rhs, sense, matbeg, this->matind.data(), this->matval.data(), NULL, cnstrname);
			if (checkCPXstatus(this->cpx_env_milp, status)) goto TERMINATE;

			/*
			t_y_k - alpha_k - y_k(UB alpha) >= -1(UB alpha)
			*/
			rhs[0] = -1;
			sense[0] = 'G';
			sprintf(cnstrname[0], "Linear_t_y_(3)_k%d_i%d", k, yi);
			nzc = 0;

			this->matind[nzc] = this->appr_dual_vars.v_t_y_k(k, yi);
			this->matval[nzc] = 1;
			++nzc;

			this->matind[nzc] = this->appr_dual_vars.v_alpha[k];
			this->matval[nzc] = -1;
			++nzc;

			this->matind[nzc] = this->v_y_k(k, yi);
			this->matval[nzc] = -1;
			++nzc;

			status = CPXaddrows(this->cpx_env_milp, this->cpx_lp_milp, 0, 1, nzc, rhs, sense, matbeg, this->matind.data(), this->matval.data(), NULL, cnstrname);
			if (checkCPXstatus(this->cpx_env_milp, status)) goto TERMINATE;
		}

	/*
	w * gamma_k
	*/
	if (this->problem->getNw() > 0)
	{
		/*
		Check that Nw = Nxi
		*/

		for (int k = 0; k < this->problem->getK(); ++k)
			for (int wi = 0; wi < this->problem->getNw(); ++wi)
			{
				/*
				* t_gamma_k <= w M
				*/
				rhs[0] = 0;
				sense[0] = 'L';
				sprintf(cnstrname[0], "Linear_t_gamma_(1)_k%d_i%d", k, wi);
				nzc = 0;

				this->matind[nzc] = this->appr_dual_vars.v_t_gamma_k(k, wi);
				this->matval[nzc] = 1;
				++nzc;

				this->matind[nzc] = this->v_w[wi];
				this->matval[nzc] = -big_m;
				++nzc;

				status = CPXaddrows(this->cpx_env_milp, this->cpx_lp_milp, 0, 1, nzc, rhs, sense, matbeg, this->matind.data(), this->matval.data(), NULL, cnstrname);
				if (checkCPXstatus(this->cpx_env_milp, status)) goto TERMINATE;

				/*
				* t_gamma_k >= -w M
				*/
				rhs[0] = 0;
				sense[0] = 'G';
				sprintf(cnstrname[0], "Linear_t_gamma_(1)_k%d_i%d", k, wi);
				nzc = 0;

				this->matind[nzc] = this->appr_dual_vars.v_t_gamma_k(k, wi);
				this->matval[nzc] = 1;
				++nzc;

				this->matind[nzc] = this->v_w[wi];
				this->matval[nzc] = big_m;
				++nzc;

				status = CPXaddrows(this->cpx_env_milp, this->cpx_lp_milp, 0, 1, nzc, rhs, sense, matbeg, this->matind.data(), this->matval.data(), NULL, cnstrname);
				if (checkCPXstatus(this->cpx_env_milp, status)) goto TERMINATE;


				/*
				t_gamma_k - gamma_k + Mw <= M
				*/
				rhs[0] = big_m;
				sense[0] = 'L';
				sprintf(cnstrname[0], "Linear_t_gamma_(3)_k%d_i%d", k, wi);
				nzc = 0;

				this->matind[nzc] = this->appr_dual_vars.v_t_gamma_k(k, wi);
				this->matval[nzc] = 1;
				++nzc;

				this->matind[nzc] = this->appr_dual_vars.v_gamma_k(k, wi);
				this->matval[nzc] = -1;
				++nzc;

				this->matind[nzc] = this->v_w[wi];
				this->matval[nzc] = big_m;
				++nzc;

				status = CPXaddrows(this->cpx_env_milp, this->cpx_lp_milp, 0, 1, nzc, rhs, sense, matbeg, this->matind.data(), this->matval.data(), NULL, cnstrname);
				if (checkCPXstatus(this->cpx_env_milp, status)) goto TERMINATE;




				/*
				t_gamma_k - gamma_k - Mw >= -M
				*/
				rhs[0] = -big_m;
				sense[0] = 'G';
				sprintf(cnstrname[0], "Linear_t_gamma_(3)_k%d_i%d", k, wi);
				nzc = 0;

				this->matind[nzc] = this->appr_dual_vars.v_t_gamma_k(k, wi);
				this->matval[nzc] = 1;
				++nzc;

				this->matind[nzc] = this->appr_dual_vars.v_gamma_k(k, wi);
				this->matval[nzc] = -1;
				++nzc;

				this->matind[nzc] = this->v_w[wi];
				this->matval[nzc] = -big_m;
				++nzc;

				status = CPXaddrows(this->cpx_env_milp, this->cpx_lp_milp, 0, 1, nzc, rhs, sense, matbeg, this->matind.data(), this->matval.data(), NULL, cnstrname);
				if (checkCPXstatus(this->cpx_env_milp, status)) goto TERMINATE;


			}
	}
TERMINATE:
	if (status)
	{
		return false;
	}

	return true;
}

bool RobustKSolver::appr_add_RLT_AlphaXWY_valid_inequalities()
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
	status = 0;
	double coeff_eps = COEFF_EPS;

	std::vector<double> b_loc;
	Matrix2D<double> Matr_loc;

	

	// X matrix constraints
	if (this->problem->getLx() > 0
		&&
		this->problem->get1S_uobj_matr_C().getNumElem() > 0)
	{

		b_loc = this->problem->get1S_drhs_vect_bx();
		Matr_loc = this->problem->get1S_dcon_matr_X();

		for (int k = 0; k < this->problem->getK(); ++k)
		{

			for (unsigned short row = 0; row < this->problem->getLx(); ++row)
			{
				sense[0] = 'L';
				nzc = 0;
				sprintf(cnstrname[0], "t_Xx_k%d_%d", k, row);

				// alpha_var: if the rhs differenc than 0 it multiplies for alpha.
				if (fabs(b_loc[row]) > 0)
				{
					rhs[0] = 0;
					this->matind[nzc] = this->appr_dual_vars.v_alpha[k];
					this->matval[nzc] = -b_loc[row];
					++nzc;
				}
				else {

					rhs[0] = b_loc[row];
				}

				// Now we go through the "columns", coefficients of the matrices 

				for (unsigned short cols = 0; cols < this->problem->getNx(); ++cols)
					if (fabs(Matr_loc(row, cols)) > coeff_eps)
					{
						this->matind[nzc] = this->appr_dual_vars.v_t_x_k(k, cols);
						this->matval[nzc] = Matr_loc(row, cols);
						++nzc;
					}
				status = CPXaddrows(this->cpx_env_milp, this->cpx_lp_milp, 0, 1, nzc, rhs, sense, matbeg, this->matind.data(), this->matval.data(), NULL, cnstrname);
				if (checkCPXstatus(this->cpx_env_milp, status)) goto TERMINATE;
			}

		}
	}
	// W matrices constraints 
	if (this->problem->getLw() > 0
		&&
		this->problem->getDV_uobj_matr_D().getNumElem() > 0)
	{
		b_loc = this->problem->getDV_drhs_vect_bw();
		Matr_loc = this->problem->getDV_dcon_matr_W();

		for (int k = 0; k < this->problem->getK(); ++k)
		{
			for (unsigned short row = 0; row < this->problem->getLw(); ++row)
			{
				sense[0] = 'L';
				nzc = 0;
				sprintf(cnstrname[0], "t_Ww_%dk_%d", k, row);
				// Now we go trought the "columns", coefficient of the matrices 

					// alpha_var: if the rhs differenc than 0 it multiplies for alpha.
				if (fabs(b_loc[row]) > 0)
				{
					rhs[0] = 0;
					this->matind[nzc] = this->appr_dual_vars.v_alpha[k];
					this->matval[nzc] = -b_loc[row];
					++nzc;
				}
				else {

					rhs[0] = b_loc[row];
				}


				for (unsigned short col = 0; col < this->problem->getNw(); ++col)
					if (fabs(Matr_loc(row, col)) > coeff_eps)
					{
						this->matind[nzc] = this->appr_dual_vars.v_t_w_k(k, col);
						this->matval[nzc] = Matr_loc(row, col);
						++nzc;
					}

				status = CPXaddrows(this->cpx_env_milp, this->cpx_lp_milp, 0, 1, nzc, rhs, sense, matbeg, this->matind.data(), this->matval.data(), NULL, cnstrname);
				if (checkCPXstatus(this->cpx_env_milp, status)) goto TERMINATE;
			}
		}
	}

	// Y matrices constraints for each k

	b_loc = this->problem->get2S_drhs_vect_by();
	Matr_loc = this->problem->get2S_dcon_matr_Y();

	if (this->problem->getLy() > 0
		&&
		this->problem->get2S_uobj_matr_Q().getNumElem())
	{
		for (int k = 0; k < this->problem->getK(); ++k)
		{
			for (unsigned short row = 0; row < this->problem->getLy(); ++row)
			{

				sense[0] = 'L';
				nzc = 0;
				sprintf(cnstrname[0], "t_Yy_k%d_%d", k, row);

				if (fabs(b_loc[row]) > 0)
				{
					rhs[0] = 0;
					this->matind[nzc] = this->appr_dual_vars.v_alpha[k];
					this->matval[nzc] = -b_loc[row];
					++nzc;
				}
				else {

					rhs[0] = b_loc[row];
				}


				// Now we go trought the "columns", coefficient of the matrices 

				for (unsigned short cols = 0; cols < this->problem->getNy(); ++cols)
					if (fabs(Matr_loc(row, cols)) > coeff_eps)
					{
						this->matind[nzc] = this->appr_dual_vars.v_t_y_k(k, cols);
						this->matval[nzc] = Matr_loc(row, cols);
						++nzc;
					}
				status = CPXaddrows(this->cpx_env_milp, this->cpx_lp_milp, 0, 1, nzc, rhs, sense, matbeg, this->matind.data(), this->matval.data(), NULL, cnstrname);
				if (checkCPXstatus(this->cpx_env_milp, status)) goto TERMINATE;

			}
		}
	}

TERMINATE:
	if (status)
	{
		return false;
	}
	return true;

}

bool RobustKSolver::appr_add_RLT_AlphaExwy_valid_inequalities()
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
	status = 0;
	double coeff_eps = COEFF_EPS;

	std::vector<double> b_loc = this->problem->get1S_drhs_vect_bx();
	Matrix2D<double> Matr_loc = this->problem->get1S_dcon_matr_X();



	if (this->problem->getR() > 0)
	{
		b_loc = this->problem->getLink_dcon_vect_g();


		for (int row = 0; row < this->problem->getR(); ++row)
			for (int k = 0; k < this->problem->getK(); ++k)
			{

				sense[0] = 'L';
				nzc = 0;
				rhs[0] = 0;

				if (fabs(b_loc[row]) > coeff_eps)
				{
					this->matind[nzc] = this->appr_dual_vars.v_alpha[k];
					this->matval[nzc] = -b_loc[row];
					++nzc;
				}

				sprintf(cnstrname[0], "LinkXWY_k%d_r%d", k, row);

				Matr_loc = this->problem->getLink_dcon_matr_Ex();
				if (Matr_loc.getNumElem() > 0)
				{
					for (int col = 0; col < this->problem->getNx(); ++col)
						if (fabs(Matr_loc(row, col)) > coeff_eps)
						{
							this->matind[nzc] = this->appr_dual_vars.v_t_x_k(k, col);
							this->matval[nzc] = Matr_loc(row, col);
							++nzc;
						}
				}
				Matr_loc = this->problem->getLink_dcon_matr_Ew();
				if (Matr_loc.getNumElem() > 0)
				{
					for (int col = 0; col < this->problem->getNw(); ++col)
						if (fabs(Matr_loc(row, col)) > coeff_eps)
						{
							this->matind[nzc] = this->appr_dual_vars.v_t_w_k(k, col);
							this->matval[nzc] = Matr_loc(row, col);
							++nzc;
						}
				}

				Matr_loc = this->problem->getLink_dcon_matr_Ey();
				if (Matr_loc.getNumElem() > 0)
				{
					for (int col = 0; col < this->problem->getNy(); ++col)
						if (fabs(Matr_loc(row, col)) > coeff_eps)
						{
							this->matind[nzc] = this->appr_dual_vars.v_t_y_k(k, col);
							this->matval[nzc] = Matr_loc(row, col);
							++nzc;
						}
				}

				status = CPXaddrows(this->cpx_env_milp, this->cpx_lp_milp, 0, 1, nzc, rhs, sense, matbeg, this->matind.data(), this->matval.data(), NULL, cnstrname);
				if (checkCPXstatus(this->cpx_env_milp, status)) goto TERMINATE;
			}
	}


TERMINATE:
	if (status)
		return false;
	else
		return true;
}

bool RobustKSolver::appr_add_RLT_Primal_valid_inequalities()
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
	status = 0;
	double coeff_eps = COEFF_EPS;

	std::vector<double> y_coeff_loc;
	std::vector<double> x_coeff_loc;


	y_coeff_loc = this->problem->get2S_ducon_vect_p();
	x_coeff_loc = this->problem->get1S_ducon_vect_t();

	for (int k = 0; k < this->problem->getK(); ++k)
		for (int l = 0; l < this->problem->getL(); ++l)
		{

			rhs[0] = 0;
			sense[0] = 'L';
			nzc = 0;

			sprintf(cnstrname[0], "RLT_Det_CU_cut_k%d_l%d", k, l);

			if (fabs(this->problem->get_rhs_ducon_vect_hmax()[l]) > coeff_eps)
			{
				this->matind[nzc] = this->appr_dual_vars.v_alpha[k];
				this->matval[nzc] = -this->problem->get_rhs_ducon_vect_hmax()[l];
				++nzc;
			}



			for (unsigned short col = 0; col < this->problem->getNx(); ++col)
				if (fabs(y_coeff_loc[col]) > coeff_eps)
				{
					this->matind[nzc] = this->appr_dual_vars.v_t_x_k(k, col);
					this->matval[nzc] = x_coeff_loc[col];
					++nzc;
				}

			for (unsigned short col = 0; col < this->problem->getNy(); ++col)
				if (fabs(y_coeff_loc[col]) > coeff_eps)
				{
					this->matind[nzc] = this->appr_dual_vars.v_t_y_k(k, col);
					this->matval[nzc] = y_coeff_loc[col];
					++nzc;
				}

			status = CPXaddrows(this->cpx_env_milp, this->cpx_lp_milp, 0, 1, nzc, rhs, sense, matbeg, this->matind.data(), this->matval.data(), NULL, cnstrname);
			if (checkCPXstatus(this->cpx_env_milp, status)) goto TERMINATE;

		}


TERMINATE:
	if (status)
		return false;
	else
		return true;
}


