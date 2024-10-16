#include "ROPEUsolver.h"

inline void ROPEUsolver::reset_branc_and_cut_indicators()
{
	num_iter_curr_node = 0;
	num_cut_curr_iter = 0;
	prev_node = 0;
	curr_node = -1;
	lp_prev_iter = 0;

}

inline int ROPEUsolver::checkCPXstatus(int status)
{
	char errmsg[CPXMESSAGEBUFSIZE];
	if (status == 0) return 0;

	CPXgeterrorstring(this->cpx_env, status, errmsg);
	printf(" %s \n", errmsg);
	return status;
}

bool ROPEUsolver::free_cplex()
{
	int status = 0;
	status = CPXfreeprob(this->cpx_env, &this->cpx_lp);
	if (checkCPXstatus(status)) return false;
	status = CPXcloseCPLEX(&this->cpx_env);
	if (checkCPXstatus(status)) return false;

	return true;
}

bool ROPEUsolver::init_data_structures()
{

	//snprintf(results_folder, sizeof(results_folder), "%s", RESULTS_FOLDER);

	//cfg = new Configuration::ConfigFile("config.cfg");

	stringstream ss;
	string out_dir; 
	ss << cfg->getValueOfKey<string>("RESULTS_FOLDER");
	out_dir = ss.str();

	if(!std::filesystem::create_directory(ss.str()))
		cout << "Warning: " << strerror(errno) << endl;
	
	snprintf(results_folder, sizeof(results_folder), "%s", ss.str().c_str());

	/*Algo Configuration*/
	this->alpha_symmetry_breaking = this->cfg->getValueOfKey<bool>("ALPHA_SYMMETRY_BREAKING");
	this->alpha_bound_reduction = this->cfg->getValueOfKey<bool>("ALPHA_BOUNDS_REDUCTION");
	this->xtilde_cuts = this->cfg->getValueOfKey<bool>("X_TILDE_CUTS");
	this->best_scenario_cuts = this->cfg->getValueOfKey<bool>("BEST_SCENARIO_CUTS");

	int num_vrs = 0;

	big_m_gamma = BIG_M_GAMMA;

	this->v_w_n.resize(this->num_nodes);
	num_vrs += this->num_nodes;

	this->v_y_kn.resizeMatrix2D(this->K, this->num_nodes);
	num_vrs += v_y_kn.getNumElem();

	this->v_x_kij.resizeMatrix3D(this->K, this->num_nodes, this->num_nodes);
	num_vrs += v_x_kij.getNumElem();

	this->v_t_y_kn.resizeMatrix2D(this->K, this->num_nodes);
	num_vrs += v_t_y_kn.getNumElem();

	this->v_t_x_kij.resizeMatrix3D(this->K, this->num_nodes, this->num_nodes);
	num_vrs += v_t_x_kij.getNumElem();


	v_alpha_k.resize(this->K); 
	num_vrs += this->K;


	v_gamma_kn.resizeMatrix2D(this->K,this->num_nodes); 
	num_vrs += v_gamma_kn.getNumElem();

	/*Bilinear*/
	v_t_gamma_kn.resizeMatrix2D(this->K, this->num_nodes); 
	num_vrs += v_t_gamma_kn.getNumElem();

	v_t_y_kn.resizeMatrix2D(this->K, this->num_nodes);
	num_vrs += v_t_y_kn.getNumElem();

	v_t_x_kij.resizeMatrix3D(this->K, this->num_nodes, this->num_nodes);
	num_vrs += v_t_x_kij.getNumElem();


	/*uncertainty set data structure: for now I assume the uncertainty set is defined by num_nodes*2 constr + 1 (one for each node plus budget.)*/

	A.assignMatrix2D((4 * this->num_nodes) + 1, this->num_nodes,0.00);

	v_beta_l.resize(this->num_nodes*2 + 1); 
	num_vrs += v_beta_l.size();

	v_beta_kl.resizeMatrix2D(this->K, this->num_nodes*2 + 1); // query with K, num_con.
	num_vrs += v_beta_kl.getNumElem();

	b.resize(this->num_nodes*2 + 1);

	best_scenario_profit_i.resize(this->num_nodes);

	/*matind and matval*/
	matind.resize(num_vrs);
	matval.resize(num_vrs);
	branch_priority_score.resize(num_vrs);

	curr_sol.resize(num_vrs);

	/*Alpha Bounds Objects*/
	this->lb_alpha_k.resize(K);
	this->ub_alpha_k.resize(K);


	/*Separations Objects*/

	sep_mcut_set_1.assign(this->num_nodes, 0);
	sep_mcut_set_2.assign(this->num_nodes, 0);

	sep_queue.elem.assign(MAX_NUM_ELE_IN_Q, 0);

	sep_un_graph.edge.resize(MAX_NUM_EDGES);
	sep_un_graph.nodes_map.resize(this->num_nodes);
	sep_un_graph.conn_comp.resize(this->num_nodes);

	sep_di_graph.ver_label.resize(this->num_nodes);
	sep_di_graph.level.resize(this->num_nodes);
	sep_di_graph.adj_archs.resizeMatrix2D(this->num_nodes,2 * this->num_nodes);
	sep_di_graph.adj_size.resize(this->num_nodes);

	node_subtour_subset.resize(this->num_nodes + 1);

	t_y_obj_coeff_val.resize(this->num_nodes);
	t_y_uc_coeff_val.resize(this->num_nodes);

	/*
	cut counters to 0
	*/
	this->num_dvar_cuts = 0;
	this->num_tildevar_cuts = 0;
	this->num_dyna_mccormik = 0;
	this->root_node_processed = false;
	this->reset_branc_and_cut_indicators();


	alrijne_node_collecting_time = cfg->getValueOfKey<double>("ALRIJNE_COLL_NODE_TIME");



	return true;
}

bool ROPEUsolver::init_out_file()
{
	if (out_file != NULL)
		return true;


	if (out_file == NULL)
	{
		char locname[1024];
		char ss[1024];

		snprintf(locname, sizeof(locname), "%s%s%s%s", results_folder,"\\","ResKapp_", inst_name);
		snprintf(ss, sizeof(ss), "%s%s", locname, ".log");
		out_file = fopen(ss, "w");

		if (out_file == NULL)
		{
			printf("Error creating output file.\n");
			return false;
		}
	}
	return true;
}

bool ROPEUsolver::clean_data_structure()
{
	/*Clean cplex*/
	int status = 0;
	status = CPXfreeprob(this->cpx_env, &this->cpx_lp);
	if (checkCPXstatus(status)) return false;
	status = CPXcloseCPLEX(&this->cpx_env);
	if (checkCPXstatus(status)) return false;


	if (out_file != NULL)
		fclose(out_file);


	delete cfg;


	return true;
}

bool ROPEUsolver::define_uncertainty_set()
{
	switch (this->uncertainty_set_type)
	{
	case 1:
		this->define_uncertainty_set_type_1();
		break;
	case 2:
		this->define_uncertainty_set_type_2();
		break;
	case 3:
		this->define_uncertainty_set_type_3();
		break;
	case 4:
		this->define_uncertainty_set_type_4();
		break;
	default:
		std::cout << "Invalid code for US type!" << std::endl;
		getchar;
		return false;
		break;
	}

	return true;
}

bool ROPEUsolver::init_cplex()
{
	int status = 0;

	// Open env.
	this->cpx_env= CPXopenCPLEX(&status);
	if (this->checkCPXstatus(status)) return false;
	this->cpx_lp = CPXcreateprob(cpx_env, &status, inst_name);
	if (this->checkCPXstatus(status)) return false;

//Z:
//	if (status)
//	{
//		return false;
//	}

	return true;
}

bool ROPEUsolver::set_cplex()
{
	bool mod_stat = true;

	int status = 0;

	char errbuf[CPXMESSAGEBUFSIZE];


	status = CPXsetintparam(this->cpx_env, CPX_PARAM_SCRIND, CPX_ON);
	if (checkCPXstatus(status)) return false;
	status = CPXsetintparam(this->cpx_env, CPX_PARAM_THREADS, 1);
	if (checkCPXstatus( status)) return false;
	status = CPXchgobjsen(this->cpx_env, this->cpx_lp, CPX_MIN);
	if (checkCPXstatus(status)) return false;
	status = CPXsetintparam(this->cpx_env, CPX_PARAM_NUMERICALEMPHASIS, CPX_ON);
	if (checkCPXstatus(status)) return false;
	status = CPXsetdblparam(this->cpx_env, CPXPARAM_Simplex_Tolerances_Feasibility, 1e-9);
	if (checkCPXstatus(status)) return false;
	//status = CPXsetintparam(cpx_env, CPX_PARAM_PREIND, CPX_OFF);
	//if (checkCPXstatus(status))return false;

	/*To be extra numerical accurate*/
	status = CPXsetdblparam(this->cpx_env, CPXPARAM_Simplex_Tolerances_Optimality, 1e-9);
	if (checkCPXstatus(status)) return false;
	status = CPXsetdblparam(this->cpx_env, CPXPARAM_MIP_Tolerances_MIPGap, 1e-9);
	if (checkCPXstatus(status)) return false;
	status = CPXsetdblparam(this->cpx_env, CPXPARAM_MIP_Tolerances_Integrality, 0);
	if (checkCPXstatus(status)) return false;

	double time_limit = this->cfg->getValueOfKey<double>("TIME_LIMIT");

	status = CPXsetdblparam(cpx_env, CPXPARAM_TimeLimit, time_limit);
	if (checkCPXstatus(status)) return false;

	/*
	status = CPXsetintparam(cpx_env, CPXPARAM_MIP_Strategy_VariableSelect, CPX_VARSEL_STRONG);
	if (checkCPXstatus(status)) return false;
	*/

	/*status = CPXsetintparam(cpx_env, CPXPARAM_MIP_Strategy_Search, CPX_MIPSEARCH_TRADITIONAL);
	if (checkCPXstatus(status)) return false;
	*/

	return true;
}

bool ROPEUsolver::define_alpha_bounds()
{
	//Tight alpha bounds
	double alpha_bound_eps = mEPS;
		
	for (int k = 0; k < K; ++k)
	{
		/*Heuristic*/
	/*	this->lb_alpha_k[k] = std::max(0.000, (1.00 / (double)K) - alpha_bound_eps);
		this->ub_alpha_k[k] = std::min(1.000, (1.00 / (double)K) + alpha_bound_eps);*/

		/*Exact*/
		if (alpha_bound_reduction)
		{
			if (k == 0)
			{
				this->lb_alpha_k[k] = std::max(0.0, (1.00 / (double)K) - alpha_bound_eps);
				this->ub_alpha_k[k] = 1;
			}
			else
			{

				this->lb_alpha_k[k] = 0;
				this->ub_alpha_k[k] = std::min(1.00, (1.00 / (double)(k + 1)) + alpha_bound_eps);

			}
	}
		else
		{
			this->lb_alpha_k[k] = 0;
			this->ub_alpha_k[k] = 1;

		}
	}
	return true;
}

bool ROPEUsolver::add_deterministic_variables()
{
	unsigned short num_cols = CPXgetnumcols(this->cpx_env, this->cpx_lp);

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
	this->v_phi = num_cols;
	obj[0] = 1;
	vartype[0] = 'C';
	sprintf(varsname[0], "Phi");
	lb[0] = -CPX_INFBOUND;
	ub[0] = CPX_INFBOUND;
	status = CPXnewcols(this->cpx_env, this->cpx_lp, 1, obj, lb, ub, vartype, varsname);
	if (checkCPXstatus(status)) goto TERMINATE;
	++num_cols;

	/*
	First stage variables W: discovery variables
	*/
	/*Only profitable nodes.*/
		for (int i = first_pnode; i <= last_pnode; ++i)
		{
			this->v_w_n[i] = num_cols;

			obj[0] = this->disc_cost_n[i]; /*This can stay into objective.*/

			vartype[0] = 'B';
			sprintf(varsname[0], "w_n%d", i);
			if (decision_depende_discovery)
			{
				lb[0] = 0;
				ub[0] = 1;
			}
			else {
				lb[0] = 1;
				ub[0] = 1;
			}

		/*	if ( 
				(i <= 8 && i != 4) || i == 13)
			{
				lb[0] = 1;
				ub[0] = 1;
			}*/


			status = CPXnewcols(this->cpx_env, this->cpx_lp, 1, obj, lb, ub, vartype, varsname);
			if (checkCPXstatus(status)) goto TERMINATE;
			++num_cols;
		}
	

	/*
	Second stage variables--> Routing vatriables, y_kn and x_kij.
	*/
		/*y_kn*/
		for (int k = 0; k < this->K; ++k)
			for (int i = first_pnode; i <= last_pnode; ++i)
			{


				this->v_y_kn(k, i) = num_cols;

				obj[0] = 0;

				vartype[0] = 'B';

				sprintf(varsname[0], "y_k%d_i%d", k, i);
				lb[0] = 0;
				ub[0] = 1;




				status = CPXnewcols(this->cpx_env, this->cpx_lp, 1, obj, lb, ub, vartype, varsname);
				if (checkCPXstatus(status)) goto TERMINATE;
				++num_cols;
			}

		/*x_kij*/
		for (int k = 0; k < this->K; ++k)
			for (int i = start_node; i <= end_node; ++i)
				for (int j = start_node; j <= end_node; ++j)
					if(this->edge_exists(i,j))
					{
					single_node_dur_1 = 0;
					single_node_dur_2 = 0;

					this->v_x_kij(k, i,j) = num_cols;

					obj[0] = 0;

					vartype[0] = 'B';

					sprintf(varsname[0], "x_k%d_i%d_j%d", k, i,j);
					lb[0] = 0;

			

					if (i != start_node && j != end_node)
					{
						/*Two possible orders.*/
						single_node_dur_1 += this->t_ij(start_node, i);
						single_node_dur_1 += this->t_ij(i, j);
						single_node_dur_1 += this->t_ij(j, end_node);


						single_node_dur_2 += this->t_ij(start_node, j);
						single_node_dur_2 += this->t_ij(i, j);
						single_node_dur_2 += this->t_ij(i, end_node);


					}

					if (i == start_node && j != end_node)
					{
						/*In this case one possible order.*/
						single_node_dur_1 += this->t_ij(start_node, j);
						single_node_dur_1 += this->t_ij(j, end_node);

						single_node_dur_2 = single_node_dur_1;

					}

					if (i != start_node && j == end_node)
					{
						/*In this case one possible order.*/
						single_node_dur_1 += this->t_ij(start_node, i);
						single_node_dur_1 += this->t_ij(i, end_node);

						single_node_dur_2 = single_node_dur_1;

					}







					if (std::min(single_node_dur_1, single_node_dur_2) > max_dur + EPS)
					{

						ub[0] = 0;

					}
					else
					{
						ub[0] = 1;
					}

					status = CPXnewcols(this->cpx_env, this->cpx_lp, 1, obj, lb, ub, vartype, varsname);
					if (checkCPXstatus(status)) goto TERMINATE;
					++num_cols;
					}

		if (BandC_SET_BRANCHING_PRIORITY)
		{
			cnt = 0;
			for (int k = 0; k < this->K; ++k)
				for (int i = first_pnode; i <= last_pnode; ++i)
				{
					matind[cnt] = v_y_kn(k, i);
					branch_priority_score[cnt] = 1;
					++cnt;

					matind[cnt] = this->v_w_n[i];
					branch_priority_score[cnt] = 1;
					++cnt;
				}



			status = CPXcopyorder(this->cpx_env, this->cpx_lp, cnt, matind.data(), this->branch_priority_score.data(),
				CPX_BRANCH_GLOBAL);
		}


TERMINATE:

	if (status)
	{
		return false;
	}

	return true;
}

bool ROPEUsolver::add_dualization_variables()
{
	int status = 0;
	int vars_ind = 0;
	double lb[1], ub[1];
	double obj[1];
	char* varsname[1];
	char elname[1024];
	varsname[0] = elname;
	char vartype[1];



	unsigned short num_cols = CPXgetnumcols(this->cpx_env, this->cpx_lp);

	// Alpha variables: size = K.
	
	for (unsigned short k = 0; k < this->K; ++k)
	{
		v_alpha_k[k] = num_cols;

		obj[0] = 0;

		vartype[0] = 'C';

		sprintf(varsname[0], "alpha_k%d", k);

		lb[0] = this->lb_alpha_k[k];
		ub[0] = this->ub_alpha_k[k];


		status = CPXnewcols(this->cpx_env, this->cpx_lp, 1, obj, lb, ub, vartype, varsname);
		if (checkCPXstatus(status)) goto TERMINATE;
		++num_cols;
	}

	// Beta variables: size = Lxi {Number of constraints in the Uncertainty set}
	for (unsigned short r = 0; r < this->num_us_constr; ++r)
	{
		
		v_beta_l[r] = num_cols;

		obj[0] = 0;

		vartype[0] = 'C';

		sprintf(varsname[0], "beta_r%d", r);

		lb[0] = 0;
		ub[0] = CPX_INFBOUND;

		//if (r == this->num_us_constr - 1)
		//	ub[0] = 0;

		status = CPXnewcols(this->cpx_env, this->cpx_lp, 1, obj, lb, ub, vartype, varsname);
		if (checkCPXstatus(status)) goto TERMINATE;
		++num_cols;
	}

	// Beta^k variables. Each Beta^k has a size of Lxi.
	for (unsigned short k = 0; k < this->K; ++k)
		for (unsigned short r = 0; r < this->num_us_constr; ++r)
		{

			v_beta_kl(k, r) = num_cols;

			obj[0] = 0;

			vartype[0] = 'C';

			sprintf(varsname[0], "betaK_k%d_r%d", k, r);

			lb[0] = 0;
			ub[0] = CPX_INFBOUND;

			/*Used to try big M*/
			/*if (r == this->num_us_constr - 1)
				ub[0] = 0;*/

			status = CPXnewcols(this->cpx_env, this->cpx_lp, 1, obj, lb, ub, vartype, varsname);
			if (checkCPXstatus(status)) goto TERMINATE;
			++num_cols;
		}
	

	// Gamma^k variables
	for (unsigned short k = 0; k < this->K; ++k)
		for (unsigned short n = this->first_pnode; n <= this->last_pnode; ++n)
		{

			v_gamma_kn(k, n) = num_cols;

			obj[0] = 0;

			vartype[0] = 'C';

			sprintf(varsname[0], "gamma_k%d_n%d", k, n);

		//	lb[0] = 0;
		//	ub[0] = 25;

			lb[0] = -CPX_INFBOUND;
			ub[0] = CPX_INFBOUND;

			//lb[0] = -1;
			//ub[0] = 0;

			status = CPXnewcols(this->cpx_env, this->cpx_lp, 1, obj, lb, ub, vartype, varsname);
			if (checkCPXstatus(status)) goto TERMINATE;
			++num_cols;
		}
	

TERMINATE:
	if (status)
	{
		return false;
	}
	return true;
}

bool ROPEUsolver::add_linearization_variables()
{
	unsigned short num_cols = CPXgetnumcols(this->cpx_env, this->cpx_lp);

	int status = 0;
	int vars_ind = 0;
	double lb[1], ub[1];
	double obj[1];
	char* varsname[1];
	char elname[1024];
	varsname[0] = elname;
	char vartype[1];
	double single_node_dur_1 = 0;
	double single_node_dur_2 = 0;


	/*
	y^k time alpha_k
	*/
	/*y_kn*/
	for (int k = 0; k < this->K; ++k)
		for (int i = first_pnode; i <= last_pnode; ++i)
		{


			this->v_t_y_kn(k, i) = num_cols;

			obj[0] = 0;

			vartype[0] = 'C';

			sprintf(varsname[0], "t_y_k%d_i%d", k, i);
			lb[0] = 0;
			//ub[0] = 1;
			ub[0] = this->ub_alpha_k[k];

			status = CPXnewcols(this->cpx_env, this->cpx_lp, 1, obj, lb, ub, vartype, varsname);
			if (checkCPXstatus(status)) goto TERMINATE;
			++num_cols;
		}

	/*x_t_kij*/
	for (int k = 0; k < this->K; ++k)
		for (int i = start_node; i <= end_node; ++i)
			for (int j = start_node; j <= end_node; ++j)
				if (this->edge_exists(i, j))
				{
					single_node_dur_1 = 0;
					single_node_dur_2 = 0;
					this->v_t_x_kij(k, i, j) = num_cols;

					obj[0] = 0;

					vartype[0] = 'C';

					sprintf(varsname[0], "t_x_k%d_i%d_j%d", k, i, j);
					//lb[0] = 0;
					//ub[0] = 1;
					lb[0] = 0;


					if (i != start_node && j != end_node)
					{
						single_node_dur_1 += this->t_ij(start_node, i);
						single_node_dur_1 += this->t_ij(i, j);
						single_node_dur_1 += this->t_ij(j, end_node);

						single_node_dur_2 += this->t_ij(start_node, j);
						single_node_dur_2 += this->t_ij(i, j);
						single_node_dur_2 += this->t_ij(i, end_node);
					}


					if (i == start_node && j != end_node)
					{
						/*In this case one possible order.*/
						single_node_dur_1 += this->t_ij(start_node, j);
						single_node_dur_1 += this->t_ij(j, end_node);

						single_node_dur_2 = single_node_dur_1;

					}

					if (i != start_node && j == end_node)
					{
						/*In this case one possible order.*/
						single_node_dur_1 += this->t_ij(start_node, i);
						single_node_dur_1 += this->t_ij(i, end_node);

						single_node_dur_2 = single_node_dur_1;

					}




					if (std::min(single_node_dur_1, single_node_dur_2) > max_dur + EPS)
					{

						ub[0] = 0;

					}
					else
					{
						ub[0] = this->ub_alpha_k[k];
					}

					status = CPXnewcols(this->cpx_env, this->cpx_lp, 1, obj, lb, ub, vartype, varsname);
					if (checkCPXstatus(status)) goto TERMINATE;
					++num_cols;
				}

	// t_gamma^k variables
	//for (unsigned short k = 0; k < this->K; ++k)
	//	for (unsigned short n = this->first_pnode; n <= this->last_pnode; ++n)
	//	{

	//		v_t_gamma_kn(k, n) = num_cols;

	//		obj[0] = 0;

	//		vartype[0] = 'C';

	//		sprintf(varsname[0], "t_gamma_k%d_n%d", k, n);

	//		lb[0] = -CPX_INFBOUND;
	//		ub[0] = CPX_INFBOUND;

	//		status = CPXnewcols(this->cpx_env, this->cpx_lp, 1, obj, lb, ub, vartype, varsname);
	//		if (checkCPXstatus(status)) goto TERMINATE;
	//		++num_cols;
	//	}


TERMINATE:
	if (status)
	{
		return false;
	}
	return true;
}

bool ROPEUsolver::add_deterministic_constraints()
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

	// Max number of discoveries
	rhs[0] = this->max_num_disc;
	sense[0] = 'L';
	nzc = 0;
	sprintf(cnstrname[0], "Max_num_discoveries");

	for (int i = this->first_pnode; i <= this->last_pnode; ++i)
	{
		matind[nzc] = this->v_w_n[i];
		matval[nzc] = 1;

		++nzc;
	}
	status = CPXaddrows(cpx_env, cpx_lp, 0, 1, nzc, rhs, sense, matbeg, matind.data(), matval.data(), NULL, cnstrname);
	if (checkCPXstatus(status)) goto TERMINATE;


	// Degree constraints on profitable nodes
	nzc = 0;
	//rhs[0] = 2;
	rhs[0] = 0;
	sense[0] = 'E';
	for (int k = 0; k < K; ++k)
		for (int i = this->first_pnode; i <= this->last_pnode; ++i)
		{
			sprintf(cnstrname[0], "Degrees_k%d_i%d", k, i);
			nzc = 0;
			for (int ii = i + 1; ii <= end_node; ++ii)
				if(this->edge_exists(i,ii))
			{
				matind[nzc] = this->v_x_kij(k, i, ii);
				matval[nzc] = 1;
				++nzc;
			}
			for (int ii = start_node; ii < i; ++ii)
				if (this->edge_exists(ii, i))
			{
				matind[nzc] = this->v_x_kij(k, ii, i);
				matval[nzc] = 1;
				++nzc;
			}

		/*	matind[nzc] = v_y_kn(k, i);
			matval[nzc] = 2;
			++nzc;*/
			matind[nzc] = v_y_kn(k, i);
			matval[nzc] = -2;
			++nzc; 

			status = CPXaddrows(cpx_env, cpx_lp, 0, 1, nzc, rhs, sense, matbeg, matind.data(), matval.data(), NULL, cnstrname);
			if (checkCPXstatus(status)) goto TERMINATE;
		}

	// Degree constraint  on start node 
	rhs[0] = 1;
	sense[0] = 'E';
	nzc = 0;
	for (int k = 0; k < K; ++k)
	{
		nzc = 0;
		sprintf(cnstrname[0], "DegreeDepotOut_s%d", k);

		for (int j = this->start_node; j <= this->end_node; ++j)
			if(this->edge_exists(start_node,j))
		{
			//matind[nzc] = g_v_xk_kij[s][g_start_node][j];
			matind[nzc] = v_x_kij(k, this->start_node, j);
			matval[nzc] = 1;
			++nzc;

		}
		status = CPXaddrows(cpx_env, cpx_lp, 0, 1, nzc, rhs, sense, matbeg, matind.data(), matval.data(), NULL, cnstrname);
		if (checkCPXstatus(status)) goto TERMINATE;

	}
	// Degree constraint  on end node 
	rhs[0] = 1;
	sense[0] = 'E';
	nzc = 0;
	for (int k = 0; k < K; ++k)
	{
		nzc = 0;
		sprintf(cnstrname[0], "DegreeDepotIn_s%d", k);

		for (int j = start_node; j <= end_node; ++j)
			if (this->edge_exists(j, end_node))
		{
			//g_matind[nzc] = g_v_xk_kij[s][j][g_end_node];
			matind[nzc] = v_x_kij(k, j, end_node);
			matval[nzc] = 1;
			++nzc;

		}
		status = CPXaddrows(cpx_env, cpx_lp, 0, 1, nzc, rhs, sense, matbeg, matind.data(), matval.data(), NULL, cnstrname);
		if (checkCPXstatus(status)) goto TERMINATE;
	}

	nzc = 0;
	rhs[0] = max_dur;
	sense[0] = 'L';
	for (int k = 0; k < K; ++k)
	{
		nzc = 0;
		sprintf(cnstrname[0], "Max_Dur_s%d", k);
		for (int i = start_node; i < end_node; ++i)
			for (int j =start_node+1; j <= end_node; ++j)
				if(this->edge_exists(i,j))
			{
			//	g_matind[nzc] = g_v_xk_kij[s][i][j];
				matind[nzc] = v_x_kij(k, i, j);
				matval[nzc] = t_ij(i,j);
				++nzc;
			}
		status = CPXaddrows(cpx_env, cpx_lp, 0, 1, nzc, rhs, sense, matbeg, matind.data(), matval.data(), NULL, cnstrname);
		if (checkCPXstatus(status)) goto TERMINATE;
	}

	//Logic
	//for (int k = 0; k < this->K; ++k)
	//	for (int i = start_node; i < end_node; ++i)
	//		for (int j = i + 1; j <= end_node; ++j)
	//		{
	//			/* Separate for integer and add also for tile if violeted */
	//			if (i != start_node)
	//			{
	//				rhs[0] = 0;
	//				sense[0] = 'L';
	//				nzc = 0;

	//				matind[nzc] = v_x_kij(k, i, j);
	//				matval[nzc] = 1;
	//				++nzc;

	//				matind[nzc] = v_y_kn(k, i);
	//				matval[nzc] = -1;
	//				++nzc;

	//				status = CPXaddrows(cpx_env, cpx_lp, 0, 1, nzc, rhs, sense, matbeg, matind.data(), matval.data(), NULL, cnstrname);
	//				if (checkCPXstatus(status)) goto TERMINATE;

	//				if (xtilde_cuts)
	//				{
	//					rhs[0] = 0;
	//					sense[0] = 'L';
	//					nzc = 0;

	//					matind[nzc] = v_t_x_kij(k, i, j);
	//					matval[nzc] = 1;
	//					++nzc;

	//					matind[nzc] = v_t_y_kn(k, i);
	//					matval[nzc] = -1;
	//					++nzc;

	//					status = CPXaddrows(cpx_env, cpx_lp, 0, 1, nzc, rhs, sense, matbeg, matind.data(), matval.data(), NULL, cnstrname);
	//					if (checkCPXstatus(status)) goto TERMINATE;

	//				}
	//			}



	//			if (j != end_node)
	//			{


	//				rhs[0] = 0;
	//				sense[0] = 'L';
	//				nzc = 0;

	//				matind[nzc] = v_x_kij(k, i, j);
	//				matval[nzc] = 1;
	//				++nzc;


	//				matind[nzc] = v_y_kn(k, j);
	//				matval[nzc] = -1;
	//				++nzc;

	//				status = CPXaddrows(cpx_env, cpx_lp, 0, 1, nzc, rhs, sense, matbeg, matind.data(), matval.data(), NULL, cnstrname);
	//				if (checkCPXstatus(status)) goto TERMINATE;

	//				if (xtilde_cuts)
	//				{
	//					rhs[0] = 0;
	//					sense[0] = 'L';
	//					nzc = 0;
	//			
	//					matind[nzc] = v_t_x_kij(k, i, j);
	//					matval[nzc] = 1;
	//					++nzc;


	//					matind[nzc] = v_t_y_kn(k, j);
	//					matval[nzc] = -1;
	//					++nzc;

	//					status = CPXaddrows(cpx_env, cpx_lp, 0, 1, nzc, rhs, sense, matbeg, matind.data(), matval.data(), NULL, cnstrname);
	//					if (checkCPXstatus(status)) goto TERMINATE;
	//				}



	//			}
	//		}




TERMINATE:
	if (status)
		return false;
	else
		return true;
}

bool ROPEUsolver::add_objective_constraint()
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

	/*\phi >= uncertainty_cost*/
	sense[0] = 'G';
	rhs[0] = 0;
	nzc = 0;
	sprintf(cnstrname[0], "Objective_constr");

	// \phi
	this->matind[nzc] = this->v_phi;
	this->matval[nzc] = 1;
	++nzc;
	// \beta and \beta_k
	for (int l = 0; l < this->num_us_constr; ++l)
	{
		this->matind[nzc] = this->v_beta_l[l];
		this->matval[nzc] = -this->b[l];
		++nzc;

		for (int k = 0; k < this->K; ++k)
		{
			this->matind[nzc] = this->v_beta_kl(k, l);
			this->matval[nzc] = -this->b[l];
			++nzc;
		}
	}


	// tilde_y
	for(int k = 0; k < this->K; ++k)
		for (int i = first_pnode; i <= last_pnode; ++i)
		{
			this->matind[nzc] = v_t_y_kn(k, i);
			this->matval[nzc] = -this->t_y_obj_coeff_val[i];
			//this->matval[nzc] = this->node_det_profit[i];
		//	this->matval[nzc] = - this->node_det_profit[i];
			++nzc;
		}



	status = CPXaddrows(this->cpx_env, this->cpx_lp, 0, 1, nzc, rhs, sense, matbeg, this->matind.data(), this->matval.data(), NULL, cnstrname);
	if (checkCPXstatus(status)) goto TERMINATE;

TERMINATE:
	if (status)
		return false;
	else return true;
}

bool ROPEUsolver::add_valid_inequalities_objective_lower_bound_constraint()
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

	for (int k = 0; k < this->K; ++k)
	{
		/*\phi >= - y * best_scenario_i*/
	sense[0] = 'G';
	rhs[0] = 0;
	nzc = 0;
	sprintf(cnstrname[0], "Objective_constr");

	// \phi
	this->matind[nzc] = this->v_phi;
	this->matval[nzc] = 1;
	++nzc;

	// y 

	for (int i = first_pnode; i <= last_pnode; ++i)
	{
		this->matind[nzc] = v_y_kn(k, i);
		this->matval[nzc] = this->best_scenario_profit_i[i];
			//this->node_det_profit[i];
			//this->max_prof[i];
		//this->matval[nzc] = -this->min_prof[i];
		++nzc;
	}
	status = CPXaddrows(this->cpx_env, this->cpx_lp, 0, 1, nzc, rhs, sense, matbeg, this->matind.data(), this->matval.data(), NULL, cnstrname);
	if (checkCPXstatus(status)) goto TERMINATE;
}


TERMINATE:
	if (status)
		return false;
	else return true;
}

bool ROPEUsolver::add_dualization_constraints()
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

	// Constraint on alpha
	nzc = 0;
	rhs[0] = 1;
	sense[0] = 'E';
	sprintf(cnstrname[0], "Alpha_constr");
	for (int k = 0; k < K; ++k)
	{
		matind[nzc] = v_alpha_k[k];
		matval[nzc] = 1;
		++nzc;
	}
	status = CPXaddrows(cpx_env, cpx_lp, 0, 1, nzc, rhs, sense, matbeg, matind.data(), matval.data(), NULL, cnstrname);
	if (checkCPXstatus(status)) goto TERMINATE;
	////

		// Uncertainty balance for each k and xi variable
		for (int i = first_pnode; i <= last_pnode; ++i)
		{
			sprintf(cnstrname[0], "Balance_i%d", i);

			nzc = 0;
			rhs[0] = 0;
			sense[0] = 'E';

	
			for (int k = 0; k < this->K; ++k)
			{
				//matind[nzc] = v_t_gamma_kn(k, i);
				matind[nzc] = v_gamma_kn(k, i);
				matval[nzc] = 1;
				++nzc;
			}

			for (int l = 0; l < this->num_us_constr; ++l)
				if (fabs(A(l, i)) > eps_coeff)
				{
					matind[nzc] = v_beta_l[l];
					matval[nzc] = A(l, i);
					++nzc;
				}

			status = CPXaddrows(cpx_env, cpx_lp, 0, 1, nzc, rhs, sense, matbeg, matind.data(), matval.data(), NULL, cnstrname);
			if (checkCPXstatus(status)) goto TERMINATE;


			
		}



	// Uncertainty balance for each k and xi variable
	for (int k = 0; k < K; ++k)
		for (int i = first_pnode; i <= last_pnode; ++i)
		{
			sprintf(cnstrname[0], "Balance_k%d_i%d", k, i);

			nzc = 0;
			rhs[0] = 0;
			sense[0] = 'E';

			matind[nzc] = v_t_y_kn(k, i);
			matval[nzc] = - this->t_y_uc_coeff_val[i];
			//matval[nzc] = 0.5 * this->node_det_profit[i];
			//matval[nzc] = - 0.5 * this->node_det_profit[i];
			++nzc;


			//matind[nzc] = v_t_gamma_kn(k,i);
			matind[nzc] = v_gamma_kn(k, i);
			matval[nzc] = -1;
			++nzc;

			for (int l = 0; l < this->num_us_constr; ++l)
				if(fabs(A(l,i)) > eps_coeff)
			{
					matind[nzc] = v_beta_kl(k, l);
					matval[nzc] = A(l, i);
					++nzc;
			}

			status = CPXaddrows(cpx_env, cpx_lp, 0, 1, nzc, rhs, sense, matbeg, matind.data(), matval.data(), NULL, cnstrname);
			if (checkCPXstatus(status)) goto TERMINATE;



			/// Additional constraints to check big M
			
			if (
				(i <= 8 && i != 4) || i == 13)
						{

			nzc = 0;
			rhs[0] = 0;
			sense[0] = 'E';


				for (int l = 0; l < this->num_us_constr; ++l)
				if(fabs(A(l,i)) > eps_coeff)
				{
					matind[nzc] = v_beta_kl(k, l);
					matval[nzc] = A(l, i);
					++nzc;
				}

			//status = CPXaddrows(cpx_env, cpx_lp, 0, 1, nzc, rhs, sense, matbeg, matind.data(), matval.data(), NULL, cnstrname);
			//if (checkCPXstatus(status)) goto TERMINATE;

			}
			


		}


TERMINATE:
	if (status)
		return false;
	else
		return true;

}

bool ROPEUsolver::add_linearization_constraints_old()
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

	// Linking gamt , gam and w (1)
	for (int k = 0; k < K; ++k)
		for (int i = first_pnode; i <=last_pnode; ++i)
		{
			nzc = 0;
			rhs[0] = big_m_gamma;
			sense[0] = 'L';
			sprintf(cnstrname[0], "Link_gamt_gam_w_1_k%d_i%d", k, i);

			matind[nzc] = v_t_gamma_kn(k, i);
				//g_vt_gammak_ki[k][i];
			matval[nzc] = 1;
			++nzc;

			matind[nzc] = v_gamma_kn(k, i);
				//g_v_gammak_ki[k][i];
			matval[nzc] = -1;
			++nzc;

			matind[nzc] = v_w_n[i];
			matval[nzc] = 
				big_m_gamma;
			++nzc;


			status = CPXaddrows(cpx_env, cpx_lp, 0, 1, nzc, rhs, sense, matbeg, matind.data(), matval.data(), NULL, cnstrname);
			if (checkCPXstatus(status)) goto TERMINATE;
		}

	// Linking gamt , gam and w (2)
	for (int k = 0; k < K; ++k)
		for (int i = first_pnode; i <= last_pnode; ++i)
		{
			nzc = 0;
			rhs[0] = -big_m_gamma;
			sense[0] = 'G';
			sprintf(cnstrname[0], "Link_gamt_gam_w_2_k%d_i%d", k, i);

			matind[nzc] = v_t_gamma_kn(k, i);
				//[k][i];
			matval[nzc] = 1;
			++nzc;

			matind[nzc] = v_gamma_kn(k,i);
			matval[nzc] = -1;
			++nzc;

			matind[nzc] = v_w_n[i];
			matval[nzc] = -big_m_gamma;
			++nzc;

			status = CPXaddrows(cpx_env, cpx_lp, 0, 1, nzc, rhs, sense, matbeg, matind.data(), matval.data(), NULL, cnstrname);
			if (checkCPXstatus(status)) goto TERMINATE;
		}

	// Linking gamt  and w (1)
	for (int k = 0; k < K; ++k)
		for (int i = first_pnode; i <= last_pnode; ++i)
		{
			nzc = 0;
			rhs[0] = 0;
			sense[0] = 'L';
			sprintf(cnstrname[0], "Link_gamt_w_1_k%d_i%d", k, i);


			matind[nzc] = v_t_gamma_kn(k,i);
			matval[nzc] = 1;
			++nzc;

			matind[nzc] = v_w_n[i];
			matval[nzc] = -big_m_gamma;
			++nzc;

			status = CPXaddrows(cpx_env, cpx_lp, 0, 1, nzc, rhs, sense, matbeg, matind.data(), matval.data(), NULL, cnstrname);
			if (checkCPXstatus(status)) goto TERMINATE;
		}

	// Linking gamt  and w (2)
	for (int k = 0; k < K; ++k)
		for (int i = first_pnode; i <= last_pnode; ++i)
		{
			nzc = 0;
			rhs[0] = 0;
			sense[0] = 'G';
			sprintf(cnstrname[0], "Link_gamt_w_2_k%d_i%d", k, i);


			matind[nzc] = v_t_gamma_kn(k,i);
			matval[nzc] = 1;
			++nzc;

			matind[nzc] = v_w_n[i];
			matval[nzc] = big_m_gamma;
			++nzc;

			status = CPXaddrows(cpx_env, cpx_lp, 0, 1, nzc, rhs, sense, matbeg, matind.data(), matval.data(), NULL, cnstrname);
			if (checkCPXstatus(status)) goto TERMINATE;
		}

	// Linking y and yk
	for (int k = 0; k < K; ++k)
		for (int i = first_pnode; i <= last_pnode; ++i)
		{
			nzc = 0;
			rhs[0] = 0;
			sense[0] = 'L';
			sprintf(cnstrname[0], "Link_yt_y_k%d_i%d", k, i);

			matind[nzc] = v_t_y_kn(k,i);
			matval[nzc] = 1;
			++nzc;

			matind[nzc] = v_y_kn(k,i);
			matval[nzc] = -1;
			++nzc;

			status = CPXaddrows(cpx_env, cpx_lp, 0, 1, nzc, rhs, sense, matbeg, matind.data(), matval.data(), NULL, cnstrname);
			if (checkCPXstatus(status)) goto TERMINATE;
		}

	// Linking yt and alphak
	for (int k = 0; k < K; ++k)
		for (int i = first_pnode; i <= last_pnode; ++i)
		{
			nzc = 0;
			rhs[0] = 0;
			sense[0] = 'L';
			sprintf(cnstrname[0], "Link_yt_alpha_k%d_i%d", k, i);

			matind[nzc] = v_t_y_kn(k,i);
			matval[nzc] = 1;
			++nzc;

			matind[nzc] = v_alpha_k[k];
			matval[nzc] = -1;
			++nzc;

			status = CPXaddrows(cpx_env, cpx_lp, 0, 1, nzc, rhs, sense, matbeg, matind.data(), matval.data(), NULL, cnstrname);
			if (checkCPXstatus(status)) goto TERMINATE;
		}

	// Link yt, y and alpha
	for (int k = 0; k < K; ++k)
		for (int i = first_pnode; i <= last_pnode; ++i)
		{

			nzc = 0;
			rhs[0] = -1;
			sense[0] = 'G';
			sprintf(cnstrname[0], "Link_yt_y_alpha_k%d_i%d", k, i);

			matind[nzc] = v_t_y_kn(k,i);
			matval[nzc] = 1;
			++nzc;

			matind[nzc] = v_alpha_k[k];
			matval[nzc] = -1;
			++nzc;

			matind[nzc] = v_y_kn(k,i);
			matval[nzc] = -1;
			++nzc;


			status = CPXaddrows(cpx_env, cpx_lp, 0, 1, nzc, rhs, sense, matbeg, matind.data(), matval.data(), NULL, cnstrname);
			if (checkCPXstatus(status)) goto TERMINATE;
		}


	///// X////

	// Linking t_x_kij and x_kij
	for (int k = 0; k < K; ++k)
		for (int i = start_node; i <= end_node; ++i)
			for (int j = start_node; j <= end_node; ++j)
				if(this->edge_exists(i,j))
		{
			nzc = 0;
			rhs[0] = 0;
			sense[0] = 'L';
			sprintf(cnstrname[0], "Link_tx_x_k%d_i%d_j%d", k, i,j);

			matind[nzc] = v_t_x_kij(k, i, j);
				//v_t_y_kn(k, i);
			matval[nzc] = 1;
			++nzc;

			matind[nzc] = v_x_kij(k, i, j);
				//v_y_kn(k, i);
			matval[nzc] = -1;
			++nzc;

			status = CPXaddrows(cpx_env, cpx_lp, 0, 1, nzc, rhs, sense, matbeg, matind.data(), matval.data(), NULL, cnstrname);
			if (checkCPXstatus(status)) goto TERMINATE;
		}

	// Linking xt and alphak
	for (int k = 0; k < K; ++k)
		for (int i = start_node; i <= end_node; ++i)
			for (int j = start_node; j <= end_node; ++j)
				if (this->edge_exists(i, j))
		{
			nzc = 0;
			rhs[0] = 0;
			sense[0] = 'L';
			sprintf(cnstrname[0], "Link_tx_alpha_k%d_i%d_j%d", k, i,j);

			matind[nzc] = v_t_x_kij(k, i, j);
				//v_t_y_kn(k, i);
			matval[nzc] = 1;
			++nzc;

			matind[nzc] = v_alpha_k[k];
			matval[nzc] = -1;
			++nzc;

			status = CPXaddrows(cpx_env, cpx_lp, 0, 1, nzc, rhs, sense, matbeg, matind.data(), matval.data(), NULL, cnstrname);
			if (checkCPXstatus(status)) goto TERMINATE;
		}

	// Link t_x, x and alpha
	for (int k = 0; k < K; ++k)
		for (int i = start_node; i <= end_node; ++i)
			for (int j = start_node; j <= end_node; ++j)
				if (this->edge_exists(i, j))
		{

			nzc = 0;
			rhs[0] = -1;
			sense[0] = 'G';
			sprintf(cnstrname[0], "Link_t_x_y_alpha_k%d_i%d_j%d", k, i,j);

			matind[nzc] = v_t_x_kij(k,i,j);
				//v_t_y_kn(k, i);
			matval[nzc] = 1;
			++nzc;

			matind[nzc] = v_alpha_k[k];
			matval[nzc] = -1;
			++nzc;

			matind[nzc] = v_x_kij(k, i, j);
				//v_y_kn(k, i);
			matval[nzc] = -1;
			++nzc;


			status = CPXaddrows(cpx_env, cpx_lp, 0, 1, nzc, rhs, sense, matbeg, matind.data(), matval.data(), NULL, cnstrname);
			if (checkCPXstatus(status)) goto TERMINATE;
		}



TERMINATE:
	if (status)
		return false;
	else
		return true;


}

bool ROPEUsolver::add_linearization_constraints()
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

	
	// Linking gamt  and w (1)
	for (int k = 0; k < K; ++k)
		for (int i = first_pnode; i <= last_pnode; ++i)
		{
			nzc = 0;
			rhs[0] = 0;
			sense[0] = 'L';
			sprintf(cnstrname[0], "Link_gamt_w_1_k%d_i%d", k, i);


			//matind[nzc] = v_t_gamma_kn(k, i);
			matind[nzc] = v_gamma_kn(k, i);
			matval[nzc] = 1;
			++nzc;

			matind[nzc] = v_w_n[i];
			matval[nzc] = -big_m_gamma;
			++nzc;

			status = CPXaddrows(cpx_env, cpx_lp, 0, 1, nzc, rhs, sense, matbeg, matind.data(), matval.data(), NULL, cnstrname);
			if (checkCPXstatus(status)) goto TERMINATE;
		}

	// Linking gamt  and w (2)
	for (int k = 0; k < K; ++k)
		for (int i = first_pnode; i <= last_pnode; ++i)
		{
			nzc = 0;
			rhs[0] = 0;
			sense[0] = 'G';
			sprintf(cnstrname[0], "Link_gamt_w_2_k%d_i%d", k, i);


			//matind[nzc] = v_t_gamma_kn(k, i);
			matind[nzc] = v_gamma_kn(k, i);
			matval[nzc] = 1;
			++nzc;

			matind[nzc] = v_w_n[i];
			matval[nzc] = big_m_gamma;
			++nzc;

			status = CPXaddrows(cpx_env, cpx_lp, 0, 1, nzc, rhs, sense, matbeg, matind.data(), matval.data(), NULL, cnstrname);
			if (checkCPXstatus(status)) goto TERMINATE;
		}

	// Linking y and yk
	for (int k = 0; k < K; ++k)
		for (int i = first_pnode; i <= last_pnode; ++i)
		{
			nzc = 0;
			rhs[0] = 0;
			sense[0] = 'L';
			sprintf(cnstrname[0], "Link_yt_y_k%d_i%d", k, i);

			matind[nzc] = v_t_y_kn(k, i);
			matval[nzc] = 1;
			++nzc;

			matind[nzc] = v_y_kn(k, i);
			matval[nzc] = -1;
			++nzc;

			status = CPXaddrows(cpx_env, cpx_lp, 0, 1, nzc, rhs, sense, matbeg, matind.data(), matval.data(), NULL, cnstrname);
			if (checkCPXstatus(status)) goto TERMINATE;
		}

	// Linking yt and alphak
	for (int k = 0; k < K; ++k)
		for (int i = first_pnode; i <= last_pnode; ++i)
		{
			nzc = 0;
			rhs[0] = 0;
			sense[0] = 'L';
			sprintf(cnstrname[0], "Link_yt_alpha_k%d_i%d", k, i);

			matind[nzc] = v_t_y_kn(k, i);
			matval[nzc] = 1;
			++nzc;

			matind[nzc] = v_alpha_k[k];
			matval[nzc] = -1;
			++nzc;

			status = CPXaddrows(cpx_env, cpx_lp, 0, 1, nzc, rhs, sense, matbeg, matind.data(), matval.data(), NULL, cnstrname);
			if (checkCPXstatus(status)) goto TERMINATE;
		}

	// Link yt, y and alpha
	for (int k = 0; k < K; ++k)
		for (int i = first_pnode; i <= last_pnode; ++i)
		{

			nzc = 0;
			rhs[0] = -1;
			sense[0] = 'G';
			sprintf(cnstrname[0], "Link_yt_y_alpha_k%d_i%d", k, i);

			matind[nzc] = v_t_y_kn(k, i);
			matval[nzc] = 1;
			++nzc;

			matind[nzc] = v_alpha_k[k];
			matval[nzc] = -1;
			++nzc;

			matind[nzc] = v_y_kn(k, i);
			matval[nzc] = -1;
			++nzc;


			status = CPXaddrows(cpx_env, cpx_lp, 0, 1, nzc, rhs, sense, matbeg, matind.data(), matval.data(), NULL, cnstrname);
			if (checkCPXstatus(status)) goto TERMINATE;
		}


TERMINATE:
	if (status)
		return false;
	else
		return true;
}

bool ROPEUsolver::add_strengthened_linearization_constraints_old()
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

	// Linking gamt , gam and w (1)
	for (int k = 0; k < K; ++k)
		for (int i = first_pnode; i <= last_pnode; ++i)
		{
			nzc = 0;
			rhs[0] = big_m_gamma;
			sense[0] = 'L';
			sprintf(cnstrname[0], "Link_gamt_gam_w_1_k%d_i%d", k, i);

			matind[nzc] = v_t_gamma_kn(k, i);
			//g_vt_gammak_ki[k][i];
			matval[nzc] = 1;
			++nzc;

			matind[nzc] = v_gamma_kn(k, i);
			//g_v_gammak_ki[k][i];
			matval[nzc] = -1;
			++nzc;

			matind[nzc] = v_w_n[i];
			matval[nzc] = big_m_gamma;
			++nzc;


			status = CPXaddrows(cpx_env, cpx_lp, 0, 1, nzc, rhs, sense, matbeg, matind.data(), matval.data(), NULL, cnstrname);
			if (checkCPXstatus(status)) goto TERMINATE;
		}

	// Linking gamt , gam and w (2)
	for (int k = 0; k < K; ++k)
		for (int i = first_pnode; i <= last_pnode; ++i)
		{
			nzc = 0;
			rhs[0] = -big_m_gamma;
			sense[0] = 'G';
			sprintf(cnstrname[0], "Link_gamt_gam_w_2_k%d_i%d", k, i);

			matind[nzc] = v_t_gamma_kn(k, i);
			//[k][i];
			matval[nzc] = 1;
			++nzc;

			matind[nzc] = v_gamma_kn(k, i);
			matval[nzc] = -1;
			++nzc;

			matind[nzc] = v_w_n[i];
			matval[nzc] = -big_m_gamma;
			++nzc;

			status = CPXaddrows(cpx_env, cpx_lp, 0, 1, nzc, rhs, sense, matbeg, matind.data(), matval.data(), NULL, cnstrname);
			if (checkCPXstatus(status)) goto TERMINATE;
		}

	// Linking gamt  and w (1)
	for (int k = 0; k < K; ++k)
		for (int i = first_pnode; i <= last_pnode; ++i)
		{
			nzc = 0;
			rhs[0] = 0;
			sense[0] = 'L';
			sprintf(cnstrname[0], "Link_gamt_w_1_k%d_i%d", k, i);


			matind[nzc] = v_t_gamma_kn(k, i);
			matval[nzc] = 1;
			++nzc;

			matind[nzc] = v_w_n[i];
			matval[nzc] = -big_m_gamma;
			++nzc;

			status = CPXaddrows(cpx_env, cpx_lp, 0, 1, nzc, rhs, sense, matbeg, matind.data(), matval.data(), NULL, cnstrname);
			if (checkCPXstatus(status)) goto TERMINATE;
		}

	// Linking gamt  and w (2)
	for (int k = 0; k < K; ++k)
		for (int i = first_pnode; i <= last_pnode; ++i)
		{
			nzc = 0;
			rhs[0] = 0;
			sense[0] = 'G';
			sprintf(cnstrname[0], "Link_gamt_w_2_k%d_i%d", k, i);


			matind[nzc] = v_t_gamma_kn(k, i);
			matval[nzc] = 1;
			++nzc;

			matind[nzc] = v_w_n[i];
			matval[nzc] = big_m_gamma;
			++nzc;

			status = CPXaddrows(cpx_env, cpx_lp, 0, 1, nzc, rhs, sense, matbeg, matind.data(), matval.data(), NULL, cnstrname);
			if (checkCPXstatus(status)) goto TERMINATE;
		}


	// Linking y and yk
	/*t_y <= y * UB alpha OK */
	for (int k = 0; k < K; ++k)
		for (int i = first_pnode; i <= last_pnode; ++i)
		{
			nzc = 0;
			rhs[0] = 0;
			sense[0] = 'L';
			sprintf(cnstrname[0], "Link_yt_y_k%d_i%d", k, i);

			matind[nzc] = v_t_y_kn(k, i);
			matval[nzc] = 1;
			++nzc;

			matind[nzc] = v_y_kn(k, i);
			matval[nzc] = -this->ub_alpha_k[k];
			++nzc;

			status = CPXaddrows(cpx_env, cpx_lp, 0, 1, nzc, rhs, sense, matbeg, matind.data(), matval.data(), NULL, cnstrname);
			if (checkCPXstatus(status)) goto TERMINATE;
		}

	// Linking yt and alphak
	/*t_y <= alpha * UB_y - LB_alpha OK  */
	for (int k = 0; k < K; ++k)
		for (int i = first_pnode; i <= last_pnode; ++i)
		{
			nzc = 0;
			rhs[0] = -this->lb_alpha_k[k];
			sense[0] = 'L';
			sprintf(cnstrname[0], "Link_yt_alpha_k%d_i%d", k, i);

			matind[nzc] = v_t_y_kn(k, i);
			matval[nzc] = 1;
			++nzc;

			matind[nzc] = v_alpha_k[k];
			matval[nzc] = -1;
			++nzc;

			matind[nzc] = v_y_kn(k, i);
			matval[nzc] = -this->lb_alpha_k[k];
			++nzc;

			status = CPXaddrows(cpx_env, cpx_lp, 0, 1, nzc, rhs, sense, matbeg, matind.data(), matval.data(), NULL, cnstrname);
			if (checkCPXstatus(status)) goto TERMINATE;
		}

	// Link yt, y and alpha
	/*t_y - alpha_k - y(UB_alpha) >= -UB_alpha *UB_y OK  */
	for (int k = 0; k < K; ++k)
		for (int i = first_pnode; i <= last_pnode; ++i)
		{

			nzc = 0;
			rhs[0] = -this->ub_alpha_k[k];
			sense[0] = 'G';
			sprintf(cnstrname[0], "Link_yt_y_alpha_k%d_i%d", k, i);

			matind[nzc] = v_t_y_kn(k, i);
			matval[nzc] = 1;
			++nzc;

			matind[nzc] = v_alpha_k[k];
			matval[nzc] = -1;
			++nzc;

			matind[nzc] = v_y_kn(k, i);
			matval[nzc] = -this->ub_alpha_k[k];
			++nzc;


			status = CPXaddrows(cpx_env, cpx_lp, 0, 1, nzc, rhs, sense, matbeg, matind.data(), matval.data(), NULL, cnstrname);
			if (checkCPXstatus(status)) goto TERMINATE;
		}


	// Link yt, y and alpha
/*t_y  - y(LB_alpha) >=0 *UB_y OK  */
	for (int k = 0; k < K; ++k)
		for (int i = first_pnode; i <= last_pnode; ++i)
		{

			nzc = 0;
			rhs[0] = 0;
			sense[0] = 'G';
			sprintf(cnstrname[0], "Link_yt_y_alpha_k%d_i%d", k, i);

			matind[nzc] = v_t_y_kn(k, i);
			matval[nzc] = 1;
			++nzc;

			matind[nzc] = v_y_kn(k, i);
			matval[nzc] = -this->lb_alpha_k[k];
			++nzc;


			status = CPXaddrows(cpx_env, cpx_lp, 0, 1, nzc, rhs, sense, matbeg, matind.data(), matval.data(), NULL, cnstrname);
			if (checkCPXstatus(status)) goto TERMINATE;
		}
	


	///// X////

	// Linking t_x_kij and x_kij
	/*t_x <= x * UB alpha */
	for (int k = 0; k < K; ++k)
		for (int i = start_node; i <= end_node; ++i)
			for (int j = start_node; j <= end_node; ++j)
				if (this->edge_exists(i, j))
				{
					nzc = 0;
					rhs[0] = 0;
					sense[0] = 'L';
					sprintf(cnstrname[0], "Link_tx_x_k%d_i%d_j%d", k, i, j);

					matind[nzc] = v_t_x_kij(k, i, j);
					//v_t_y_kn(k, i);
					matval[nzc] = 1;
					++nzc;

					matind[nzc] = v_x_kij(k, i, j);
					//v_y_kn(k, i);
					matval[nzc] = -this->ub_alpha_k[k];
					++nzc;

					/*	status = CPXaddrows(cpx_env, cpx_lp, 0, 1, nzc, rhs, sense, matbeg, matind.data(), matval.data(), NULL, cnstrname);
						if (checkCPXstatus(status)) goto TERMINATE;*/
				}


	// Linking xt and alphak
		/*t_x <= alpha * UB_x - LB_alpha */
	for (int k = 0; k < K; ++k)
		for (int i = start_node; i <= end_node; ++i)
			for (int j = start_node; j <= end_node; ++j)
				if (this->edge_exists(i, j))
				{
					nzc = 0;
					rhs[0] = -this->lb_alpha_k[k];
					sense[0] = 'L';
					sprintf(cnstrname[0], "Link_tx_alpha_k%d_i%d_j%d", k, i, j);

					matind[nzc] = v_t_x_kij(k, i, j);
					//v_t_y_kn(k, i);
					matval[nzc] = 1;
					++nzc;

					matind[nzc] = v_x_kij(k, i, j);
					//v_t_y_kn(k, i);
					matval[nzc] = -this->lb_alpha_k[k];
					++nzc;

					matind[nzc] = v_alpha_k[k];
					matval[nzc] = -1;
					++nzc;

				/*	status = CPXaddrows(cpx_env, cpx_lp, 0, 1, nzc, rhs, sense, matbeg, matind.data(), matval.data(), NULL, cnstrname);
					if (checkCPXstatus(status)) goto TERMINATE;*/
				}

	// Link t_x, x and alpha
		/*t_x - alpha_k - x(UB_alpha) >= -UB_alpha *UB_y   */
	for (int k = 0; k < K; ++k)
		for (int i = start_node; i <= end_node; ++i)
			for (int j = start_node; j <= end_node; ++j)
				if (this->edge_exists(i, j))
				{

					nzc = 0;
					rhs[0] = - this->ub_alpha_k[k];
					sense[0] = 'G';
					sprintf(cnstrname[0], "Link_t_x_y_alpha_k%d_i%d_j%d", k, i, j);

					matind[nzc] = v_t_x_kij(k, i, j);
					//v_t_y_kn(k, i);
					matval[nzc] = 1;
					++nzc;

					matind[nzc] = v_alpha_k[k];
					matval[nzc] = -1;
					++nzc;

					matind[nzc] = v_x_kij(k, i, j);
					//v_y_kn(k, i);
					matval[nzc] = -this->ub_alpha_k[k];
					++nzc;


				/*status = CPXaddrows(cpx_env, cpx_lp, 0, 1, nzc, rhs, sense, matbeg, matind.data(), matval.data(), NULL, cnstrname);
				if (checkCPXstatus(status)) goto TERMINATE;*/
				}

	// Link xt, x
/*t_x  - x(LB_alpha) >=0 *UB_y OK  */
	for (int k = 0; k < K; ++k)
		for (int i = start_node; i <= end_node; ++i)
			for (int j = start_node; j <= end_node; ++j)
				if (this->edge_exists(i, j))
		{

			nzc = 0;
			rhs[0] = 0;
			sense[0] = 'G';
			sprintf(cnstrname[0], "Link_xt_x_alpha_k%d_i%d", k, i);

			matind[nzc] = v_t_x_kij(k, i,j);
			matval[nzc] = 1;
			++nzc;

			matind[nzc] = v_x_kij(k, i,j);
			matval[nzc] = -this->lb_alpha_k[k];
			++nzc;


			/*status = CPXaddrows(cpx_env, cpx_lp, 0, 1, nzc, rhs, sense, matbeg, matind.data(), matval.data(), NULL, cnstrname);
			if (checkCPXstatus(status)) goto TERMINATE;*/
		}


TERMINATE:
	if (status)
		return false;
	else
		return true;
}

bool ROPEUsolver::add_strengthened_linearization_constraints()
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



	// Linking gamt  and w (1)
	for (int k = 0; k < K; ++k)
		for (int i = first_pnode; i <= last_pnode; ++i)
		{
			nzc = 0;
			rhs[0] = 0;
			sense[0] = 'L';
			sprintf(cnstrname[0], "Link_gamt_w_1_k%d_i%d", k, i);



			matind[nzc] = v_gamma_kn(k, i);
			matval[nzc] = 1;
			++nzc;

			matind[nzc] = v_w_n[i];
			matval[nzc] = -big_m_gamma; // Original line of code
			++nzc;

			status = CPXaddrows(cpx_env, cpx_lp, 0, 1, nzc, rhs, sense, matbeg, matind.data(), matval.data(), NULL, cnstrname);
			if (checkCPXstatus(status)) goto TERMINATE;
		}

	// Linking gamt  and w (2)
	for (int k = 0; k < K; ++k)
		for (int i = first_pnode; i <= last_pnode; ++i)
		{
			nzc = 0;
			rhs[0] = 0;
			sense[0] = 'G';
			sprintf(cnstrname[0], "Link_gamt_w_2_k%d_i%d", k, i);


			matind[nzc] = v_gamma_kn(k, i);
			matval[nzc] = 1;
			++nzc;

			matind[nzc] = v_w_n[i];
			matval[nzc] = big_m_gamma; //Original line of code
			++nzc;

			status = CPXaddrows(cpx_env, cpx_lp, 0, 1, nzc, rhs, sense, matbeg, matind.data(), matval.data(), NULL, cnstrname);
			if (checkCPXstatus(status)) goto TERMINATE;
		}


	//CPXaddindcontr()
	//CPXaddindconstr()
	rhs[0] = 0;
	sense[0] = 'E';
	for (int k = 0; k < K; ++k)
		for (int i = first_pnode; i <= last_pnode; ++i)
		{
			sprintf(cnstrname[0], "Indicator_w_gamma_k%d_i%d", k, i);
			nzc = 0;

			matind[nzc] = v_gamma_kn(k, i);
			matval[nzc] = 1;
			++nzc;

			//status = CPXaddindconstr(cpx_env, cpx_lp, this->v_w_n[i], 1, nzc, rhs[0], sense[0], matind.data(), matval.data(), cnstrname[0]);
			//if (checkCPXstatus(status)) goto TERMINATE;

		}

	// Linking y and yk
	/*t_y <= y * UB alpha OK */
	for (int k = 0; k < K; ++k)
		for (int i = first_pnode; i <= last_pnode; ++i)
		{
			nzc = 0;
			rhs[0] = 0;
			sense[0] = 'L';
			sprintf(cnstrname[0], "Link_(1)_yt_y_k%d_i%d", k, i);

			matind[nzc] = v_t_y_kn(k, i);
			matval[nzc] = 1;
			++nzc;

			matind[nzc] = v_y_kn(k, i);
			matval[nzc] = -this->ub_alpha_k[k];
			++nzc;

			status = CPXaddrows(cpx_env, cpx_lp, 0, 1, nzc, rhs, sense, matbeg, matind.data(), matval.data(), NULL, cnstrname);
			if (checkCPXstatus(status)) goto TERMINATE;
		}

	// Linking yt and alphak
	/*t_y <= alpha * UB_y - LB_alpha OK  */
	for (int k = 0; k < K; ++k)
		for (int i = first_pnode; i <= last_pnode; ++i)
		{
			nzc = 0;
			rhs[0] = -this->lb_alpha_k[k];
			sense[0] = 'L';
			sprintf(cnstrname[0], "Link_(2)_yt_alpha_k%d_i%d", k, i);

			matind[nzc] = v_t_y_kn(k, i);
			matval[nzc] = 1;
			++nzc;

			matind[nzc] = v_alpha_k[k];
			matval[nzc] = -1;
			++nzc;

			matind[nzc] = v_y_kn(k, i);
			matval[nzc] = -this->lb_alpha_k[k];
			++nzc;

			status = CPXaddrows(cpx_env, cpx_lp, 0, 1, nzc, rhs, sense, matbeg, matind.data(), matval.data(), NULL, cnstrname);
			if (checkCPXstatus(status)) goto TERMINATE;
		}

	// Link yt, y and alpha
	/*t_y - alpha_k - y(UB_alpha) >= -UB_alpha *UB_y OK  */
	for (int k = 0; k < K; ++k)
		for (int i = first_pnode; i <= last_pnode; ++i)
		{

			nzc = 0;
			rhs[0] = -this->ub_alpha_k[k];
			sense[0] = 'G';
			sprintf(cnstrname[0], "Link_(3)_yt_y_alpha_k%d_i%d", k, i);

			matind[nzc] = v_t_y_kn(k, i);
			matval[nzc] = 1;
			++nzc;

			matind[nzc] = v_alpha_k[k];
			matval[nzc] = -1;
			++nzc;

			matind[nzc] = v_y_kn(k, i);
			matval[nzc] = -this->ub_alpha_k[k];
			++nzc;


			status = CPXaddrows(cpx_env, cpx_lp, 0, 1, nzc, rhs, sense, matbeg, matind.data(), matval.data(), NULL, cnstrname);
			if (checkCPXstatus(status)) goto TERMINATE;
		}


	// Link yt, y and alpha
/*t_y  - y(LB_alpha) >=0 *UB_y OK  */
	for (int k = 0; k < K; ++k)
		for (int i = first_pnode; i <= last_pnode; ++i)
		{

			nzc = 0;
			rhs[0] = 0;
			sense[0] = 'G';
			sprintf(cnstrname[0], "Link_(4)_yt_y_alpha_k%d_i%d", k, i);

			matind[nzc] = v_t_y_kn(k, i);
			matval[nzc] = 1;
			++nzc;

			matind[nzc] = v_y_kn(k, i);
			matval[nzc] = -this->lb_alpha_k[k];
			++nzc;


			status = CPXaddrows(cpx_env, cpx_lp, 0, 1, nzc, rhs, sense, matbeg, matind.data(), matval.data(), NULL, cnstrname);
			if (checkCPXstatus(status)) goto TERMINATE;
		}


TERMINATE:
	if (status)
		return false;
	else
		return true;
}

bool ROPEUsolver::add_valid_inequalities_alpha_symmetries()
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

	double coeff_eps =mEPS;
	double big_m = BIG_M_GAMMA;
	double alpha_bound_eps = mEPS;
	status = 0;

	int cnt = 0;
	char lu[2];
	double bd[2];
	unsigned int num_k;

	// alpha_k >= alpha_(k+1)
	sense[0] = 'G';
	rhs[0] = 0;
	nzc = 0;

	for (int k = 0; k < K - 1; ++k)
	{
		nzc = 0;

		sprintf(cnstrname[0], "Alpha_k%d_geq_alpha_kk%d ", k, k + 1);

		this->matind[nzc] = v_alpha_k[k];
		this->matval[nzc] = 1;
		++nzc;

		this->matind[nzc] = v_alpha_k[k + 1];
		this->matval[nzc] = -1;
		++nzc;

		if (nzc > 0)
		{
			status = CPXaddrows(this->cpx_env, this->cpx_lp, 0, 1, nzc, rhs, sense, matbeg, this->matind.data(), this->matval.data(), NULL, cnstrname);
			if (checkCPXstatus(status)) goto TERMINATE;
		}
	}
	/*
	Tight alpha bounds: already defined
	*/
	//num_k = K;
	//for (int k = 0; k < K; ++k)
	//	{
	//		cnt = 0;
	//		this->matind[cnt] = this->v_alpha_k[k];

	//		if (k == 0)
	//		{

	//			this->lb_alpha_k[k] = (1.00 / (double)num_k) - alpha_bound_eps;

	//			bd[cnt] = lb_alpha_k[k];
	//			lu[cnt] = 'L';
	//			++cnt;
	//			status = CPXchgbds(this->cpx_env, this->cpx_lp, cnt, this->matind.data(), lu, bd);
	//			if (checkCPXstatus(status)) goto TERMINATE;


	//			cnt = 0;
	//			this->ub_alpha_k[k] = 1;
	//			bd[cnt] = 1;
	//			lu[cnt] = 'U';
	//			++cnt;
	//			status = CPXchgbds(this->cpx_env, this->cpx_lp, cnt, this->matind.data(), lu, bd);
	//			if (checkCPXstatus( status)) goto TERMINATE;

	//		}
	//		else
	//		{
	//			this->lb_alpha_k[k] = 0;

	//			bd[cnt] = lb_alpha_k[k];
	//			lu[cnt] = 'L';
	//			++cnt;
	//			status = CPXchgbds(this->cpx_env, this->cpx_lp, cnt, this->matind.data(), lu, bd);
	//			if (checkCPXstatus(status)) goto TERMINATE;

	//			cnt = 0;
	//			this->ub_alpha_k[k] = (1.00 / (double)(k + 1)) + alpha_bound_eps;
	//			bd[cnt] = ub_alpha_k[k];
	//			lu[cnt] = 'U';
	//			++cnt;
	//			status = CPXchgbds(this->cpx_env, this->cpx_lp, cnt, this->matind.data(), lu, bd);
	//			if (checkCPXstatus(status)) goto TERMINATE;
	//		}
	//	}


TERMINATE:
	if (status)
	{
		return false;
	}

	return true;
}

bool ROPEUsolver::add_RLT_valid_inequalties()
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



	// Degree constraints on profitable nodes
	nzc = 0;
	//rhs[0] = 2;
	rhs[0] = 0;
	sense[0] = 'E';
	for (int k = 0; k < K; ++k)
		for (int i = this->first_pnode; i <= this->last_pnode; ++i)
		{
			sprintf(cnstrname[0], "Degrees_k%d_i%d", k, i);
			nzc = 0;

		/*	matind[nzc] = v_alpha_k[k];
			matval[nzc] = -2;
			++nzc;*/
		

			for (int ii = i + 1; ii <= end_node; ++ii)
				if (this->edge_exists(i, ii))
				{
					matind[nzc] = this->v_t_x_kij(k, i, ii);
					matval[nzc] = 1;
					++nzc;
				}
			for (int ii = start_node; ii < i; ++ii)
				if (this->edge_exists(ii, i))
				{
					matind[nzc] = this->v_t_x_kij(k, ii, i);
					matval[nzc] = 1;
					++nzc;
				}


			matind[nzc] = v_t_y_kn(k, i);
			matval[nzc] = -2;
			++nzc;

			/*matind[nzc] = v_t_y_kn(k, i);
			matval[nzc] = 2;
			++nzc;*/

			status = CPXaddrows(cpx_env, cpx_lp, 0, 1, nzc, rhs, sense, matbeg, matind.data(), matval.data(), NULL, cnstrname);
			if (checkCPXstatus(status)) goto TERMINATE;
		}

	// Degree constraint  on start node 
	//rhs[0] = 1;
	rhs[0] = 0;
	sense[0] = 'E';
	nzc = 0;
	for (int k = 0; k < K; ++k)
	{
		nzc = 0;

		matind[nzc] = v_alpha_k[k];
		matval[nzc] = -1;
		++nzc;

		sprintf(cnstrname[0], "DegreeDepotIn_s%d", k);

		for (int j = this->start_node; j <= this->end_node; ++j)
			if (this->edge_exists(start_node, j))
			{
				
				matind[nzc] = v_t_x_kij(k, this->start_node, j);
				matval[nzc] = 1;
				++nzc;

			}
		status = CPXaddrows(cpx_env, cpx_lp, 0, 1, nzc, rhs, sense, matbeg, matind.data(), matval.data(), NULL, cnstrname);
		if (checkCPXstatus(status)) goto TERMINATE;

	}
	// Degree constraint  on end node tilde
	//rhs[0] = 1;
	rhs[0] = 0;
	sense[0] = 'E';
	nzc = 0;
	for (int k = 0; k < K; ++k)
	{
		nzc = 0;

		matind[nzc] = v_alpha_k[k];
		matval[nzc] = -1;
		++nzc;

		sprintf(cnstrname[0], "DegreeDepotIn_s%d", k);

		for (int j = start_node; j <= end_node; ++j)
			if (this->edge_exists(j, end_node))
			{
				
				matind[nzc] = v_t_x_kij(k, j, end_node);
				matval[nzc] = 1;
				++nzc;

			}
		status = CPXaddrows(cpx_env, cpx_lp, 0, 1, nzc, rhs, sense, matbeg, matind.data(), matval.data(), NULL, cnstrname);
		if (checkCPXstatus(status)) goto TERMINATE;
	}

	nzc = 0;
	rhs[0] = max_dur;
	rhs[0] = 0;
	sense[0] = 'L';
	for (int k = 0; k < K; ++k)
	{
		nzc = 0;

		matind[nzc] = v_alpha_k[k];
		matval[nzc] = - max_dur;
		++nzc;

		sprintf(cnstrname[0], "Max_Dur_s%d", k);
		for (int i = start_node; i <= end_node; ++i)
			for (int j = start_node; j <= end_node; ++j)
				if (this->edge_exists(i, j))
				{
				
					matind[nzc] = v_t_x_kij(k, i, j);
					matval[nzc] = t_ij(i, j);
					++nzc;
				}
		status = CPXaddrows(cpx_env, cpx_lp, 0, 1, nzc, rhs, sense, matbeg, matind.data(), matval.data(), NULL, cnstrname);
		if (checkCPXstatus(status)) goto TERMINATE;
	}

	//Mccormic x alpha

	//separate_McCormick_x_alpha
	
		//for (int k = 0; k < this->K; ++k)
		//	for (int i = start_node; i < end_node; ++i)
		//		for (int j = i + 1; j <= end_node; ++j)
		//		{

		//			//Mc1
		//			nzc = 0;
		//			rhs[0] = 0;
		//			sense[0] = 'L';


		//			matind[nzc] = v_t_x_kij(k, i, j);
		//			matval[nzc] = 1;
		//			++nzc;

		//			matind[nzc] = v_x_kij(k, i, j);
		//			matval[nzc] = -this->ub_alpha_k[k];
		//			++nzc;

		//			status = CPXaddrows(cpx_env, cpx_lp, 0, 1, nzc, rhs, sense, matbeg, matind.data(), matval.data(), NULL, cnstrname);
		//			if (checkCPXstatus(status)) goto TERMINATE;

		//			//Mc2
		//			nzc = 0;
		//			rhs[0] = -this->lb_alpha_k[k];
		//			sense[0] = 'L';

		//			matind[nzc] = v_t_x_kij(k, i, j);
		//			matval[nzc] = 1;
		//			++nzc;

		//			matind[nzc] = v_x_kij(k, i, j);
		//			matval[nzc] = -this->lb_alpha_k[k];
		//			++nzc;

		//			matind[nzc] = v_alpha_k[k];
		//			matval[nzc] = -1;
		//			++nzc;
		//			status = CPXaddrows(cpx_env, cpx_lp, 0, 1, nzc, rhs, sense, matbeg, matind.data(), matval.data(), NULL, cnstrname);
		//			if (checkCPXstatus(status)) goto TERMINATE;

		//			//Mc3
		//			nzc = 0;
		//			rhs[0] = -this->ub_alpha_k[k];
		//			sense[0] = 'G';

		//			matind[nzc] = v_t_x_kij(k, i, j);
		//			matval[nzc] = 1;
		//			++nzc;

		//			matind[nzc] = v_alpha_k[k];
		//			matval[nzc] = -1;
		//			++nzc;

		//			matind[nzc] = v_x_kij(k, i, j);
		//			matval[nzc] = -this->ub_alpha_k[k];
		//			++nzc;
		//			status = CPXaddrows(cpx_env, cpx_lp, 0, 1, nzc, rhs, sense, matbeg, matind.data(), matval.data(), NULL, cnstrname);
		//			if (checkCPXstatus(status)) goto TERMINATE;

		//			//Mc4
		//			nzc = 0;
		//			rhs[0] = 0;
		//			sense[0] = 'G';

		//			matind[nzc] = v_t_x_kij(k, i, j);
		//			matval[nzc] = 1;
		//			++nzc;

		//			matind[nzc] = v_x_kij(k, i, j);
		//			matval[nzc] = -this->lb_alpha_k[k];
		//			++nzc;
		//			status = CPXaddrows(cpx_env, cpx_lp, 0, 1, nzc, rhs, sense, matbeg, matind.data(), matval.data(), NULL, cnstrname);
		//			if (checkCPXstatus(status)) goto TERMINATE;
		//		}
	



TERMINATE:
	if (status)
		return false;
	else
		return true;
}

inline int CPXPUBLIC ROPEUsolver::general_callback(CPXCALLBACKCONTEXTptr context, CPXLONG contextid, void* userdata)
{
	bool cut_added = false;
	int status = 0;
	bool is_ok = true;
	ROPEUsolver* solver = (ROPEUsolver*)userdata;
	double val;
	double treshold;


	if (contextid == CPX_CALLBACKCONTEXT_CANDIDATE)
	{

		status = CPXcallbackgetcandidatepoint(context, solver->curr_sol.data(), 0, CPXgetnumcols(solver->cpx_env, solver->cpx_lp) - 1, &val);
		if (status != 0)
		{
			solver->checkCPXstatus(status);
			return status;
		}

		if (!cut_added)
			is_ok = solver->separate_heuristic_subtours(context, contextid, 0, &cut_added);

		if (!cut_added)
			if (solver->xtilde_cuts)
			{
				is_ok = solver->separate_McCormick_x_alpha(context, contextid, &cut_added);
			}
		


	}

	if (contextid == CPX_CALLBACKCONTEXT_RELAXATION)
	{
		/*Only CPX 12.10*/
		CPXLONG node_id, node_depth;

		CPXcallbackgetinfolong(context, CPXCALLBACKINFO_NODEUID, &node_id);
		if (status != 0)
		{
			solver->checkCPXstatus(status);
			return status;
		}

		CPXcallbackgetinfolong(context, CPXCALLBACKINFO_NODEDEPTH, &node_depth);
		if (status != 0)
		{
			solver->checkCPXstatus(status);
			return status;
		}




		/*	CPXLONG node_count;
			status = CPXcallbackgetinfolong(context, CPXCALLBACKINFO_NODECOUNT, &node_count);
			if (status != 0)
			{
				solver->checkCPXstatus(status);
				return status;
			}*/

		status = CPXcallbackgetrelaxationpoint(context, solver->curr_sol.data(), 0, CPXgetnumcols(solver->cpx_env, solver->cpx_lp) - 1, &val);
		if (status != 0)
		{
			solver->checkCPXstatus(status);
			return status;
		}


		if (node_id > 1) /*Do this after root node.*/
		{
			if (node_id == solver->curr_node)
			{
				++solver->num_iter_curr_node;
				//num_cut_curr_iter = 0;
				//prev_node = 0;
				solver->curr_node = node_id;

				//if (fabs(val - solver->lp_prev_iter) < fabs(BandC_MIN_COEFF_LP_IMPR * solver->lp_prev_iter)
				//	&&
				//	solver->num_iter_curr_node >= BandC_MAX_NUM_ITER_NODE
				//	)
				//{
				//	return 0; // Branch.
				//}

				//if (
				///*	val - solver->lp_prev_iter < (BandC_MIN_COEFF_LP_IMPR * solver->lp_prev_iter)
				//	&& */
				//	solver->num_cut_curr_iter >= BandC_NUM_CUT_PER_NODE)
				//{
				//	return 0; // Branch.
				//}

			}
			else {
				solver->reset_branc_and_cut_indicators();
				solver->curr_node = node_id;
				solver->lp_prev_iter = val;
				solver->num_cut_curr_iter = 0;

			}
		}






		if (solver->xtilde_cuts)
		{
			is_ok = solver->separate_McCormick_x_alpha(context, contextid, &cut_added);
		}

		
		is_ok = solver->separate_logic_cuts(context, contextid, &cut_added);


		if (!cut_added ||
			node_id <= 1
			)
		{
			treshold = 0;

			while (!cut_added
				||
				node_depth <= 1
				)
			{
				is_ok = solver->separate_heuristic_subtours(context, contextid, treshold, &cut_added);

				if (solver->xtilde_cuts)
				{
					if (!cut_added
						||
						node_depth <= 1
						)
						is_ok = solver->separate_tilde_heuristic_subtours(context, contextid, treshold, &cut_added);
				}

				treshold += TRESHOLD_STEP;

				if (treshold > TRESHOLD_MAX)
					break;
			}
		}


		if (node_depth <= 1
			&&
			!solver->root_node_processed)
		{
			if (!cut_added || node_id <= 1)
			is_ok = solver->separate_exact_subtuors(context, contextid, &cut_added);
			
			if (solver->xtilde_cuts)
			{
				if (!cut_added || node_id <= 1)
				is_ok = solver->separate_tilde_exact_subtours(context, contextid, &cut_added);
			}

			solver->root_node_lb = val;
			if (!cut_added)
			{
				solver->root_node_processed = true;
			}

		}

		//if (solver->xtilde_cuts)
		//{
		//	if (!cut_added ||
		//		node_id <= 1
		//		)
		//	{
		//		treshold = 0;

		//		while (!cut_added
		//			||
		//			node_depth <= 1
		//			)
		//		{
		//			//is_ok = solver->separate_heuristic_subtours(context, contextid, treshold, &cut_added);

		//			if (solver->xtilde_cuts)
		//			{
		//				if (!cut_added
		//					||
		//					node_depth <= 1
		//					)
		//					is_ok = solver->separate_tilde_heuristic_subtours(context, contextid, treshold, &cut_added);
		//			}

		//			treshold += TRESHOLD_STEP;

		//			if (treshold > TRESHOLD_MAX)
		//				break;
		//		}
		//	}

		//	/*Next set*/
		//	if (node_depth <= 1
		//		&&
		//		!solver->root_node_processed)
		//	{
		//		//if (!cut_added || node_id <= 1)
		//		//is_ok = solver->separate_exact_subtuors(context, contextid, &cut_added);
		//		if (solver->xtilde_cuts)
		//		{
		//			//if (!cut_added || node_id <= 1)
		//			is_ok = solver->separate_tilde_exact_subtours(context, contextid, &cut_added);
		//		}

		//		solver->root_node_lb = val;

		//		if (!cut_added)
		//		{
		//			solver->root_node_processed = true;
		//		}

		//	}

		//}

	}

	return 0;
		}
	

inline bool ROPEUsolver::separate_heuristic_subtours(CPXCALLBACKCONTEXTptr context, CPXLONG contextid, double var_treshold ,bool* cut_added)
{
	int status = 0;
	double val;

	//std::vector<short> k_ind;
	//k_ind.assign(this->K, 0);
	//double score = 0;
	//double max_score = 0;
	//int k = 0;
	//bool keep_going = true;
	
	
	
	//while(keep_going)
	//{
	//	k = -1;
	//for (int k_cut = 0; k_cut < this->K; ++k_cut)
	//	if (k_ind[k_cut] == 0
	//		&& 
	//		(contextid == CPX_CALLBACKCONTEXT_CANDIDATE
	//					||
	//					curr_sol[v_alpha_k[k]] > cEPS)
	//		)
	//	{
	//		score = min(fabs(this->ub_alpha_k[k] - curr_sol[this->v_alpha_k[k]]), fabs(this->lb_alpha_k[k] - curr_sol[this->v_alpha_k[k]]));
	//		if (score > max_score)
	//		{
	//			 k = k_cut;
	//			max_score = score;
	//		}
	//	}
	//if(k != -1)
	//	k_ind[k] = 1;
	//else
	//{
	//	keep_going = false;
	//	break;
	//}

	for (int k = 0; k < this->K; ++k)
		if(contextid == CPX_CALLBACKCONTEXT_CANDIDATE
			||
			curr_sol[v_alpha_k[k]] > cEPS
			)
		{ 
	
		build_sol_graph(curr_sol.data(), k, var_treshold);
		connected_components();

		if (sep_un_graph.num_cc > 1)
		{
			double obj[1], rhs[1];
			int matbeg[1], nzc;
			char sense[1];
			sense[0] = 'L';
			matbeg[0] = 0;
			rhs[0] = 0;
			double lhs_v = 0;
			double rhs_v = 0;

			int force[1];
			int local[1];
			force[0] = CPX_USECUT_FORCE;
			//force[0] = CPX_USECUT_FILTER;
			local[0] = 0;
			int yi_max = 0;
			double yi_max_val = 0;


			for (int cc = 0; cc < sep_un_graph.num_cc; ++cc)
				if (cc != sep_un_graph.conn_comp[0]
					&&
					cc != sep_un_graph.conn_comp[sep_un_graph.num_nodes - 1])
				{
					nzc = 0;
					rhs[0] = 0;
					lhs_v = 0;
					// Saving connected component
				//	short connected_comp[MAX_NUM_NODES];
					node_subtour_subset[0] = 0;
					bool first = true;


					
					for (int v = 1; v < sep_un_graph.num_nodes; ++v)
						if (sep_un_graph.conn_comp[v] == cc)
						{
							++node_subtour_subset[0];
							node_subtour_subset[node_subtour_subset[0]] = sep_un_graph.nodes_map[v];

						}

					if (node_subtour_subset[0] < 2)
					{
						continue;
					}


					for (int kk = 0; kk < this->K; ++kk)
					{

						nzc = 0;
						rhs[0] = 0;
						lhs_v = 0;
						yi_max_val = 0;
						yi_max = 1;

						for (int i = 1; i <= node_subtour_subset[0]; ++i)
						{
							if (this->curr_sol[v_y_kn(kk, node_subtour_subset[i])] > yi_max_val)
							{
								yi_max = i;
								yi_max_val = this->curr_sol[v_y_kn(kk, node_subtour_subset[i])];
							}
						}

						
							for (int i = 1; i <= node_subtour_subset[0]; ++i)
							{
								if (i != yi_max)
								{

									matind[nzc] = this->v_y_kn(kk, node_subtour_subset[i]);
									matval[nzc] = -1;

									++nzc;

									lhs_v -= this->curr_sol[v_y_kn(kk, node_subtour_subset[i])];

								}
							}
						for (int i = 1; i < node_subtour_subset[0]; ++i)
							for (int j = i + 1; j <= node_subtour_subset[0]; ++j)
							{
								matind[nzc] = v_x_kij(kk, node_subtour_subset[i], node_subtour_subset[j]);
								matval[nzc] = 1;
								++nzc;


								lhs_v += this->curr_sol[v_x_kij(kk, node_subtour_subset[i], node_subtour_subset[j])];
							}
						if (nzc > 0
							&&
							lhs_v > rhs[0] + cEPS
							)
						{
							if (contextid == CPX_CALLBACKCONTEXT_CANDIDATE)
							{

								status = CPXcallbackrejectcandidate(context, 1, nzc, rhs,
									sense, matbeg, matind.data(), matval.data());
								if (status)checkCPXstatus(status);
							}
							else
								if (contextid == CPX_CALLBACKCONTEXT_RELAXATION)
								{

									++this->num_cut_curr_iter;


									//std::cout << "Cut" << std::endl;
									status = CPXcallbackaddusercuts(context, 1, nzc, rhs,
										sense, matbeg, matind.data(), matval.data(), force, local);
									if (status)checkCPXstatus(status);
								}
							++this->num_dvar_cuts;



							/*Tilde cut*/
							if (xtilde_cuts
									&&
									(
										curr_sol[v_alpha_k[k]] > this->lb_alpha_k[k] + mEPS
										&&
										curr_sol[v_alpha_k[k]] < this->ub_alpha_k[k] - mEPS
										)
								)
							{
								nzc = 0;
								rhs[0] = 0;
								lhs_v = 0;
								yi_max = 1;
								yi_max_val = 0;
								for (int i = 1; i <= node_subtour_subset[0]; ++i)
								{
									if (this->curr_sol[v_t_y_kn(kk, node_subtour_subset[i])] > yi_max_val)
									{
										yi_max = i;
										yi_max_val = this->curr_sol[v_t_y_kn(kk, node_subtour_subset[i])];
									}
								}

								for (int i = 1; i <= node_subtour_subset[0]; ++i)
								{
									if (i != yi_max)
									{
										matind[nzc] = this->v_t_y_kn(kk, node_subtour_subset[i]);
										matval[nzc] = -1;

										++nzc;


										lhs_v -= this->curr_sol[v_t_y_kn(kk, node_subtour_subset[i])];


									}
								}

								for (int i = 1; i < node_subtour_subset[0]; ++i)
									for (int j = i + 1; j <= node_subtour_subset[0]; ++j)
									{
										matind[nzc] = v_t_x_kij(kk, node_subtour_subset[i], node_subtour_subset[j]);
										matval[nzc] = 1;
										++nzc;
										lhs_v += this->curr_sol[v_t_x_kij(kk, node_subtour_subset[i], node_subtour_subset[j])];
									}



								if (nzc > 0
									&&
									lhs_v > rhs[0]
									+ cEPS
									)
								{



									if (contextid == CPX_CALLBACKCONTEXT_CANDIDATE)
									{
										status = CPXcallbackrejectcandidate(context, 1, nzc, rhs,
											sense, matbeg, matind.data(), matval.data());
										if (status)checkCPXstatus(status);
									}
									else
										if (contextid == CPX_CALLBACKCONTEXT_RELAXATION)
										{

										//	++this->num_cut_curr_iter;


											status = CPXcallbackaddusercuts(context, 1, nzc, rhs,
												sense, matbeg, matind.data(), matval.data(), force, local);
											if (status)checkCPXstatus(status);
										}

									++this->num_tildevar_cuts;
								}
								(*cut_added) = true;
								if (status != 0)
								{

									return false;
								}
							}

							//return true;
						}
					
					}
				}
		}
	}

	return true;
}

inline bool ROPEUsolver::separate_tilde_heuristic_subtours(CPXCALLBACKCONTEXTptr context, CPXLONG contextid, double var_treshold,bool* cut_added)
{
	int status = 0;
	double val;

	/*std::vector<short> k_ind;
	k_ind.assign(this->K, 0);
	double score = 0;
	double max_score = 0;
	int k = 0;
	bool keep_going = true;*/



	/*while (keep_going)
	{
		k = -1;
		for (int k_cut = 0; k_cut < this->K; ++k_cut)
			if (k_ind[k_cut] == 0
				&&
				(contextid == CPX_CALLBACKCONTEXT_CANDIDATE
					||
					curr_sol[v_alpha_k[k]] > cEPS)
				)
			{
				score = min(fabs(this->ub_alpha_k[k] - curr_sol[this->v_alpha_k[k]]), fabs(this->lb_alpha_k[k] - curr_sol[this->v_alpha_k[k]]));
				if (score > max_score)
				{
					k = k_cut;
					max_score = score;
				}
			}
		if (k != -1)
			k_ind[k] = 1;
		else
		{
			keep_going = false;
			break;
		}*/

		for (int k = 0; k < this->K; ++k)
		if (contextid == CPX_CALLBACKCONTEXT_CANDIDATE
			||
			(curr_sol[v_alpha_k[k]] > this->lb_alpha_k[k] +  mEPS
				&&
				curr_sol[v_alpha_k[k]] < this->ub_alpha_k[k] - mEPS
				)
			)
		{ 

		
			build_tilde_sol_graph(curr_sol.data(), k, var_treshold * curr_sol[v_alpha_k[k]]);
			connected_components();

			if (sep_un_graph.num_cc > 1)
			{
				double obj[1], rhs[1];
				int matbeg[1], nzc;
				char sense[1];
				sense[0] = 'L';
				matbeg[0] = 0;
				rhs[0] = 0;
				double lhs_v = 0;
				double rhs_v = 0;

				int force[1];
				int local[1];
				force[0] = CPX_USECUT_FORCE;
				//force[0] = CPX_USECUT_FILTER;
				local[0] = 0;
				int yi_max = 0;
				double yi_max_val = 0;


				for (int cc = 0; cc < sep_un_graph.num_cc; ++cc)
					if (cc != sep_un_graph.conn_comp[0]
						&&
						cc != sep_un_graph.conn_comp[sep_un_graph.num_nodes - 1])
					{
						nzc = 0;
						rhs[0] = 0;
						lhs_v = 0;
						// Saving connected component for the Tilde cut
						//short connected_comp[MAX_NUM_NODES];
						node_subtour_subset[0] = 0;



						for (int v = 1; v < sep_un_graph.num_nodes; ++v)
							if (sep_un_graph.conn_comp[v] == cc)
							{
								++node_subtour_subset[0];
								node_subtour_subset[node_subtour_subset[0]] = sep_un_graph.nodes_map[v];

							}

						if (node_subtour_subset[0] < 2)
						{
							continue;
						}

						/* Now we first have tilde cut*/
						for (int kk = 0; kk < this->K; ++kk)
						{

							nzc = 0;
							rhs[0] = 0;
							lhs_v = 0;
							yi_max_val = 0;
							yi_max = 1;

							for (int i = 1; i <= node_subtour_subset[0]; ++i)
							{
								if (this->curr_sol[v_y_kn(kk, node_subtour_subset[i])] > yi_max_val)
								{
									yi_max = i;
									yi_max_val = this->curr_sol[v_t_y_kn(kk, node_subtour_subset[i])];
								}
							}


							for (int i = 1; i <= node_subtour_subset[0]; ++i)
							{
								if (i != yi_max)
								{

									matind[nzc] = this->v_t_y_kn(kk, node_subtour_subset[i]);
									matval[nzc] = -1;

									++nzc;

									lhs_v -= this->curr_sol[v_t_y_kn(kk, node_subtour_subset[i])];

								}
							}
							for (int i = 1; i < node_subtour_subset[0]; ++i)
								for (int j = i + 1; j <= node_subtour_subset[0]; ++j)
								{
									matind[nzc] = v_t_x_kij(kk, node_subtour_subset[i], node_subtour_subset[j]);
									matval[nzc] = 1;
									++nzc;


									lhs_v += this->curr_sol[v_t_x_kij(kk, node_subtour_subset[i], node_subtour_subset[j])];
								}
							if (nzc > 0
								&&
								lhs_v > rhs[0] + cEPS
								)
							{
								
								if (contextid == CPX_CALLBACKCONTEXT_CANDIDATE)
								{
									status = CPXcallbackrejectcandidate(context, 1, nzc, rhs,
										sense, matbeg, matind.data(), matval.data());
									if (status)checkCPXstatus(status);
								}
								else
									if (contextid == CPX_CALLBACKCONTEXT_RELAXATION)
									{
										status = CPXcallbackaddusercuts(context, 1, nzc, rhs,
											sense, matbeg, matind.data(), matval.data(), force, local);
										if (status)checkCPXstatus(status);

										++this->num_cut_curr_iter;

									}
								++this->num_tildevar_cuts;

								/*Integer variables*/
								nzc = 0;
								rhs[0] = 0;
								lhs_v = 0;
								yi_max = 1;
								yi_max_val = 0;
								for (int i = 1; i <= node_subtour_subset[0]; ++i)
								{
									if (this->curr_sol[v_y_kn(k, node_subtour_subset[i])] > yi_max_val)
									{
										yi_max = i;
										yi_max_val = this->curr_sol[v_y_kn(k, node_subtour_subset[i])];
									}
								}

								for (int i = 1; i <= node_subtour_subset[0]; ++i)
								{
									if (i != yi_max)
									{
										matind[nzc] = this->v_y_kn(kk, node_subtour_subset[i]);
										matval[nzc] = -1;

										++nzc;


										lhs_v -= this->curr_sol[v_y_kn(kk, node_subtour_subset[i])];


									}
								}

								for (int i = 1; i < node_subtour_subset[0]; ++i)
									for (int j = i + 1; j <= node_subtour_subset[0]; ++j)
									{
										matind[nzc] = v_x_kij(kk, node_subtour_subset[i], node_subtour_subset[j]);
										matval[nzc] = 1;
										++nzc;
										lhs_v += this->curr_sol[v_x_kij(kk, node_subtour_subset[i], node_subtour_subset[j])];
									}



								if (nzc > 0
									&&
									lhs_v > rhs[0]
									+ cEPS
									)
								{
									if (contextid == CPX_CALLBACKCONTEXT_CANDIDATE)
									{
										status = CPXcallbackrejectcandidate(context, 1, nzc, rhs,
											sense, matbeg, matind.data(), matval.data());
										if (status)checkCPXstatus(status);
									}
									else
										if (contextid == CPX_CALLBACKCONTEXT_RELAXATION)
										{
											status = CPXcallbackaddusercuts(context, 1, nzc, rhs,
												sense, matbeg, matind.data(), matval.data(), force, local);
											if (status)checkCPXstatus(status);

											++this->num_cut_curr_iter;

										}

									++this->num_dvar_cuts;

								}
								(*cut_added) = true;
								if (status != 0)
								{

									return false;
								}

								//return true;
							}
						}
					}
			}
		}

	return true;
}

inline bool ROPEUsolver::separate_exact_subtuors(CPXCALLBACKCONTEXTptr context, CPXLONG contextid,bool* cut_added)
{

	int status = 0;
	double val = 0;
	double obj[1], rhs[1];
	int matbeg[1], nzc;
	char sense[1];
	int force[1];
	int local[1];
	//force[0] = CPX_USECUT_FILTER;
	force[0] = CPX_USECUT_FORCE;
	local[0] = 0;
	/*std::vector<short> k_ind;
	k_ind.assign(this->K, 0);
	double score = 0;
	double max_score = 0;
	int k = 0;
	bool keep_going = true;*/

	//short connected_comp[MAX_NUM_NODES];
	double lhs_val = 0;

	//std::cout << "Here!" << std::endl;

//	for (int k = 0; k < this->K; ++k)
//{
//	for (int i = this->first_pnode; i <= this->last_pnode; ++i)
//		if (curr_sol[this->v_y_kn(k,i)] > mEPS)
//		{
//			std::cout << "y_k" << k << "_i" << i << " = "<< curr_sol[this->v_y_kn(k, i)] << std::endl;
//		}
//
//	for(int i = start_node; i <= end_node; ++i)
//		for (int j = start_node; j <= end_node; ++j)
//			if (this->edge_exists(i, j)
//				&&
//				curr_sol[this->v_x_kij(k, i, j)] > mEPS
//				)
//			{
//				std::cout << "x_k" << k << "_i" << i <<"_j"<<j <<" = "<< curr_sol[this->v_x_kij(k, i, j)] << std::endl;
//
//			}
//
//}
	/*while (keep_going)
	{
		k = -1;
		for (int k_cut = 0; k_cut < this->K; ++k_cut)
			if (k_ind[k_cut] == 0
				&&
				(contextid == CPX_CALLBACKCONTEXT_CANDIDATE
					||
					curr_sol[v_alpha_k[k]] > cEPS)
				)
			{
				score = min(fabs(this->ub_alpha_k[k] - curr_sol[this->v_alpha_k[k]]), fabs(this->lb_alpha_k[k] - curr_sol[this->v_alpha_k[k]]));
				if (score > max_score)
				{
					k = k_cut;
					max_score = score;
				}
			}
		if (k != -1)
			k_ind[k] = 1;
		else
		{
			keep_going = false;
			break;
		}*/

	for(int k = 0; k < this->K; ++k)
		if (curr_sol[v_alpha_k[k]] > cEPS)
		{
			build_sol_directed_graph(curr_sol.data(), k);
			// GSEC bsed on max flow
			node_subtour_subset[0] = 0;
			double yi_max_val = 0;
			int yi_max = 0;
			for (int s = 0; s < sep_di_graph.num_ver - 1; ++s)
				for (int t = s + 1; t < sep_di_graph.num_ver; ++t)
				{

					node_subtour_subset[0] = 0;

					int flow = dinic_max_flow(s, t);
					if (flow < (2 * scaleMXflow))
					{
						sense[0] = 'L';
						matbeg[0] = 0;
						rhs[0] = 0;

						for (int v = start_node + 1;
							//v < sep_di_graph.ver_label[sep_di_graph.num_ver - 1];
							v < end_node;
							++v)
							if (sep_mcut_set_1[v] == 1)
							{
								++node_subtour_subset[0];
								node_subtour_subset[node_subtour_subset[0]] = v;
							}


						if (node_subtour_subset[0] < 2)
						{
							continue;
						}


						for (int kk = 0; kk < K; ++kk)
						{
							nzc = 0;
							lhs_val = 0;
							//	node_subtour_subset[0] = 0;
							yi_max_val = 0;
							yi_max = 1;


							for (int i = 1; i <= node_subtour_subset[0]; ++i)
							{
								if (this->curr_sol[v_y_kn(kk, node_subtour_subset[i])] > yi_max_val)
								{
									yi_max = i;
									yi_max_val = this->curr_sol[v_y_kn(kk, node_subtour_subset[i])];
								}
							}


							for (int i = 1; i <= node_subtour_subset[0]; ++i)
							{
								if (i != yi_max)
								{

									matind[nzc] = v_y_kn(kk, node_subtour_subset[i]);
									matval[nzc] = -1;



									lhs_val -= curr_sol[v_y_kn(kk, node_subtour_subset[i])];
									++nzc;


								}
							}

							for (int i = 1; i < node_subtour_subset[0]; ++i)
								for (int j = i + 1; j <= node_subtour_subset[0]; ++j)
								{
									matind[nzc] = v_x_kij(kk, node_subtour_subset[i], node_subtour_subset[j]);
									matval[nzc] = 1;
									lhs_val += curr_sol[v_x_kij(kk, node_subtour_subset[i], node_subtour_subset[j])];

									++nzc;

								}

							if (lhs_val > rhs[0] + cEPS)
							{
								/*std::cout << "Cut!" << std::endl;
								getchar();*/

								status = CPXcallbackaddusercuts(context, 1, nzc, rhs,
									sense, matbeg, matind.data(), matval.data(), force, local);
								if (status)checkCPXstatus(status);

								++this->num_dvar_cuts;

								++this->num_cut_curr_iter;


							}



							if (status != 0)
							{
								return false;
							}

							/*Tilde cut*/
							if(this->xtilde_cuts
								&&
							(curr_sol[v_alpha_k[k]] > this->lb_alpha_k[k] + mEPS
								&&
								curr_sol[v_alpha_k[k]] < this->ub_alpha_k[k] - mEPS)
								)
							{
							nzc = 0;
							lhs_val = 0;
							//node_subtour_subset[0] = 0;
							yi_max_val = 0;
							yi_max = 1;


							for (int i = 1; i <= node_subtour_subset[0]; ++i)
							{
								if (this->curr_sol[v_t_y_kn(kk, node_subtour_subset[i])] > yi_max_val)
								{
									yi_max = i;
									yi_max_val = this->curr_sol[v_t_y_kn(kk, node_subtour_subset[i])];
								}
							}
							for (int i = 1; i <= node_subtour_subset[0]; ++i)
							{
								if (i != yi_max)
								{

									matind[nzc] = v_t_y_kn(kk, node_subtour_subset[i]);
									matval[nzc] = -1;
									lhs_val -= curr_sol[v_t_y_kn(kk, node_subtour_subset[i])];
									++nzc;


								}
							}

							for (int i = 1; i < node_subtour_subset[0]; ++i)
								for (int j = i + 1; j <= node_subtour_subset[0]; ++j)
								{
									matind[nzc] = v_t_x_kij(kk, node_subtour_subset[i], node_subtour_subset[j]);
									matval[nzc] = 1;
									lhs_val += curr_sol[v_t_x_kij(kk, node_subtour_subset[i], node_subtour_subset[j])];

									++nzc;

								}

							if (lhs_val > rhs[0] + cEPS)
							{
								/*	std::cout << "Cut!" << std::endl;
									getchar();*/
								status = CPXcallbackaddusercuts(context, 1, nzc, rhs,
									sense, matbeg, matind.data(), matval.data(), force, local);
								if (status)checkCPXstatus(status);

								++this->num_tildevar_cuts;

							//	++this->num_cut_curr_iter;


							}
							if (status != 0)
							{
								return false;
							}
						}
						}
					}
				}
		}
	
	return true;
}

inline bool ROPEUsolver::separate_tilde_exact_subtours(CPXCALLBACKCONTEXTptr context, CPXLONG contextid, bool* cut_added)
{
	int status = 0;
	double val = 0;
	double obj[1], rhs[1];
	int matbeg[1], nzc;
	char sense[1];
	int force[1];
	int local[1];
	force[0] = CPX_USECUT_FORCE;
	//force[0] - CPX_USECUT_FILTER;
	local[0] = 0;

	//short connected_comp[MAX_NUM_NODES];
	double lhs_val = 0;
	//std::vector<short> k_ind;
	//k_ind.assign(this->K, 0);
	//double score = 0;
	//double max_score = 0;
	//int k = 0;
	//bool keep_going = true;

	//std::cout << "Here!" << std::endl;

	for (int k = 0; k < this->K; ++k)
		if (contextid == CPX_CALLBACKCONTEXT_CANDIDATE
			||
			(curr_sol[v_alpha_k[k]] > this->lb_alpha_k[k] + mEPS
				&&
				curr_sol[v_alpha_k[k]] < this->ub_alpha_k[k] - mEPS
				))
		if (curr_sol[v_alpha_k[k]] > cEPS)
		{





	/*while (keep_going)
	{
		k = -1;
		for (int k_cut = 0; k_cut < this->K; ++k_cut)
			if (k_ind[k_cut] == 0
				&&
				(contextid == CPX_CALLBACKCONTEXT_CANDIDATE
					||
					curr_sol[v_alpha_k[k]] > cEPS)
				)
			{
				score = min(fabs(this->ub_alpha_k[k] - curr_sol[this->v_alpha_k[k]]), fabs(this->lb_alpha_k[k] - curr_sol[this->v_alpha_k[k]]));
				if (score > max_score)
				{
					k = k_cut;
					max_score = score;
				}
			}
		if (k != -1)
			k_ind[k] = 1;
		else
		{
			keep_going = false;
			break;
		}*/
			build_tilde_sol_directed_graph(curr_sol.data(), k);
			// GSEC bsed on max flow
			node_subtour_subset[0] = 0;
			double yi_max_val = 0;
			int yi_max = 0;
			for (int s = 0; s < sep_di_graph.num_ver - 1; ++s)
				for (int t = s + 1; t < sep_di_graph.num_ver; ++t)
				{

					node_subtour_subset[0] = 0;

					int flow = dinic_max_flow(s, t);
					if (flow < (2 * scaleMXflow))
					{
						sense[0] = 'L';
						matbeg[0] = 0;
						rhs[0] = 0;

						for (int v = start_node + 1; 
							//v < sep_di_graph.ver_label[sep_di_graph.num_ver - 1];
							v < end_node;
							++v)
							if (sep_mcut_set_1[v] == 1)
							{
								++node_subtour_subset[0];
								node_subtour_subset[node_subtour_subset[0]] = v;
							}





						if (node_subtour_subset[0] < 2)
						{
							continue;
						}




						for (int kk = 0; kk < K; ++kk)
						{
							nzc = 0;
							lhs_val = 0;
							node_subtour_subset[0] = 0;
							yi_max_val = 0;
							yi_max = 1;


							for (int i = 1; i <= node_subtour_subset[0]; ++i)
							{
								if (this->curr_sol[v_t_y_kn(kk, node_subtour_subset[i])] > yi_max_val)
								{
									yi_max = i;
									yi_max_val = this->curr_sol[v_t_y_kn(kk, node_subtour_subset[i])];
								}
							}


							for (int i = 1; i <= node_subtour_subset[0]; ++i)
							{
								if (i != yi_max)
								{

									matind[nzc] = v_t_y_kn(kk, node_subtour_subset[i]);
									matval[nzc] = -1;



									lhs_val -= curr_sol[v_t_y_kn(kk, node_subtour_subset[i])];
									++nzc;


								}
							}

							for (int i = 1; i < node_subtour_subset[0]; ++i)
								for (int j = i + 1; j <= node_subtour_subset[0]; ++j)
								{
									matind[nzc] = v_t_x_kij(kk, node_subtour_subset[i], node_subtour_subset[j]);
									matval[nzc] = 1;
									lhs_val += curr_sol[v_t_x_kij(kk, node_subtour_subset[i], node_subtour_subset[j])];

									++nzc;

								}

							if (lhs_val > rhs[0] + cEPS)
							{
								*cut_added = true;
			/*					std::cout << "Cut!" << std::endl;
								getchar();*/

								status = CPXcallbackaddusercuts(context, 1, nzc, rhs,
									sense, matbeg, matind.data(), matval.data(), force, local);
								if (status)checkCPXstatus(status);

								++this->num_tildevar_cuts;

							}
							if (status != 0)
							{
								return false;
							}

							/*Cut for integer variables*/

							nzc = 0;
							lhs_val = 0;
							node_subtour_subset[0] = 0;
							yi_max_val = 0;
							yi_max = 1;


							for (int i = 1; i <= node_subtour_subset[0]; ++i)
							{
								if (this->curr_sol[v_y_kn(kk, node_subtour_subset[i])] > yi_max_val)
								{
									yi_max = i;
									yi_max_val = this->curr_sol[v_y_kn(kk, node_subtour_subset[i])];
								}
							}
							for (int i = 1; i <= node_subtour_subset[0]; ++i)
							{
								if (i != yi_max)
								{

									matind[nzc] = v_y_kn(kk, node_subtour_subset[i]);
									matval[nzc] = -1;
									lhs_val -= curr_sol[v_y_kn(kk, node_subtour_subset[i])];
									++nzc;


								}
							}

							for (int i = 1; i < node_subtour_subset[0]; ++i)
								for (int j = i + 1; j <= node_subtour_subset[0]; ++j)
								{
									matind[nzc] = v_x_kij(kk, node_subtour_subset[i], node_subtour_subset[j]);
									matval[nzc] = 1;
									lhs_val += curr_sol[v_x_kij(kk, node_subtour_subset[i], node_subtour_subset[j])];

									++nzc;

								}

							if (lhs_val > rhs[0] + cEPS)
							{
								*cut_added = true;

							
								status = CPXcallbackaddusercuts(context, 1, nzc, rhs,
									sense, matbeg, matind.data(), matval.data(), force, local);
								if (status)checkCPXstatus(status);

								++this->num_dvar_cuts;

								++this->num_cut_curr_iter;

							}
							if (status != 0)
							{
								return false;
							}
						}
					}
				}
		}

	return true;
}

inline bool ROPEUsolver::separate_logic_cuts(CPXCALLBACKCONTEXTptr context, CPXLONG contextid, bool* cut_added)
{
	int status = 0;
	double val = 0;
	double obj[1], rhs[1];
	int matbeg[1], nzc;
	char sense[1];
	int force[1];
	int local[1];
	force[0] = CPX_USECUT_FORCE;
	local[0] = 0;
	matbeg[0] = 0;
	double lhs_v = 0;

	if(false){
	int i_star, j_star, k_star, y_star;
	bool tilde_var = false;
	double max_viol = 0;
	for (int k = 0; k < this->K; ++k)
		if (curr_sol[v_alpha_k[k]] > cEPS)
			for (int i = start_node; i < end_node; ++i)
				for (int j = i + 1; j <= end_node; ++j)
				{
					if (i != start_node)
					{
						if (curr_sol[v_x_kij(k, i, j)] - curr_sol[v_y_kn(k, i)] >= max_viol)
						{
							max_viol = curr_sol[v_x_kij(k, i, j)] - curr_sol[v_y_kn(k, i)];
							tilde_var = false;
							i_star = i;
							j_star = j;
							k_star = k;
							y_star = i;

						}

						if (curr_sol[v_t_x_kij(k, i, j)] - curr_sol[v_t_y_kn(k, i)] >= max_viol)
						{
							max_viol = curr_sol[v_t_x_kij(k, i, j)] - curr_sol[v_t_y_kn(k, i)];
							tilde_var = true;
							i_star = i;
							j_star = j;
							k_star = k;
							y_star = i;

						}

					}
					if (j != end_node)
					{
						if (curr_sol[v_x_kij(k, i, j)] - curr_sol[v_y_kn(k, j)] >= max_viol)
						{
							max_viol = curr_sol[v_x_kij(k, i, j)] - curr_sol[v_y_kn(k, j)];
							tilde_var = false;
							i_star = i;
							j_star = j;
							k_star = k;
							y_star = j;
						}

						if (curr_sol[v_t_x_kij(k, i, j)] - curr_sol[v_t_y_kn(k, j)] >= max_viol)
						{
							max_viol = curr_sol[v_t_x_kij(k, i, j)] - curr_sol[v_t_y_kn(k, j)];
							tilde_var = true;
							i_star = i;
							j_star = j;
							k_star = k;
							y_star = j;
						}



					}
				}

	if (max_viol > cEPS)
	{

		switch (tilde_var)
		{
		case true:

			rhs[0] = 0;
			sense[0] = 'L';
			nzc = 0;
			matind[nzc] = v_t_x_kij(k_star, i_star, j_star);
			matval[nzc] = 1;
			++nzc;
			matind[nzc] = v_t_y_kn(k_star, y_star);
			matval[nzc] = -1;
			++nzc;
			++this->num_tildevar_cuts;
			break;


		case false:
			rhs[0] = 0;
			sense[0] = 'L';
			nzc = 0;
			matind[nzc] = v_x_kij(k_star, i_star, j_star);
			matval[nzc] = 1;
			++nzc;
			matind[nzc] = v_y_kn(k_star, y_star);
			matval[nzc] = -1;
			++nzc;
			++this->num_dvar_cuts;
			break;

		default:
			return false;


		}



		status = CPXcallbackaddusercuts(context, 1, nzc, rhs,
			sense, matbeg, matind.data(), matval.data(), force, local);

		++this->num_cut_curr_iter;

		if (status)
		{
			checkCPXstatus(status);
			//fprintf(stdout, "Cut Logic!\n")

			return false;
		}
	}

}

	//if (true){
		// Logic VI
		for (int k = 0; k < this->K; ++k)
			if (curr_sol[v_alpha_k[k]] > cEPS)
				for (int i = start_node; i < end_node; ++i)
					for (int j = i + 1; j <= end_node; ++j)
					{

						/* Separate for integer and add also for tile if violeted */
						if (i != start_node
							&&
							curr_sol[v_x_kij(k, i, j)] > curr_sol[v_y_kn(k, i)] + cEPS)
						{
							*cut_added = true;

							// Add the cuts for i on y and x. and for each k.
							for (int kk = 0; kk < this->K; ++kk)
								//	for (int kk = k; kk <= k; ++kk)
							{
								rhs[0] = 0;
								sense[0] = 'L';
								nzc = 0;
								lhs_v = 0;

								matind[nzc] = v_x_kij(kk, i, j);
								matval[nzc] = 1;
								lhs_v += curr_sol[v_x_kij(kk, i, j)];
								++nzc;

								matind[nzc] = v_y_kn(kk, i);
								matval[nzc] = -1;
								lhs_v -= curr_sol[v_y_kn(kk, i)];
								++nzc;

								if (lhs_v > rhs[0] + cEPS)
								{
									status = CPXcallbackaddusercuts(context, 1, nzc, rhs,
										sense, matbeg, matind.data(), matval.data(), force, local);

									++this->num_cut_curr_iter;
								//	return true;

									if (status)
									{
										checkCPXstatus(status);
										//fprintf(stdout, "Cut Logic!\n")

										return false;
									}
									++this->num_dvar_cuts;

								}


							}
						}

						if (j != end_node
							&&
							curr_sol[v_x_kij(k, i, j)] > curr_sol[v_y_kn(k, j)] + cEPS)
						{
							*cut_added = true;

							for (int kk = 0; kk < this->K; ++kk)
								//for (int kk = k; kk <= k; ++kk)
							{
								// Add the cuts for j, x and y and for each k
								rhs[0] = 0;
								sense[0] = 'L';
								nzc = 0;
								lhs_v = 0;

								matind[nzc] = v_x_kij(kk, i, j);
								matval[nzc] = 1;
								++nzc;
								lhs_v += curr_sol[v_x_kij(kk, i, j)];

								matind[nzc] = v_y_kn(kk, j);
								matval[nzc] = -1;
								lhs_v -= curr_sol[v_y_kn(kk, j)];
								++nzc;

								if (lhs_v > rhs[0] + mEPS)
								{

									++this->num_cut_curr_iter;


									status = CPXcallbackaddusercuts(context, 1, nzc, rhs,
										sense, matbeg, matind.data(), matval.data(), force, local);
									if (status)
									{
										checkCPXstatus(status);
										//fprintf(stdout, "Cut Logic!\n")
										return false;
									}
								//	return true;

									++this->num_dvar_cuts;
								}
							}
						}

						/*Nowe separate for Tilde and check also Intger*/
						if (this->xtilde_cuts && false)
						{
							if (i != start_node
								&&
								curr_sol[v_t_x_kij(k, i, j)] > curr_sol[v_t_y_kn(k, i)] + cEPS)
							{
								*cut_added = true;
								// Add the cuts for i on y and x. and for each k.
								for (int kk = 0; kk < this->K; ++kk)
									for (int kk = k; kk <= k; ++kk)

									{

										rhs[0] = 0;
										sense[0] = 'L';
										nzc = 0;
										lhs_v = 0;


										matind[nzc] = v_t_x_kij(kk, i, j);
										matval[nzc] = 1;
										lhs_v += curr_sol[v_t_x_kij(kk, i, j)];
										++nzc;

										matind[nzc] = v_t_y_kn(kk, i);
										matval[nzc] = -1;
										lhs_v -= curr_sol[v_t_y_kn(kk, i)];
										++nzc;


										if (lhs_v > rhs[0] + cEPS)
										{

											++this->num_cut_curr_iter;

											status = CPXcallbackaddusercuts(context, 1, nzc, rhs,
												sense, matbeg, matind.data(), matval.data(), force, local);
											if (status)
											{
												checkCPXstatus(status);
												//fprintf(stdout, "Cut Logic!\n")
												return false;
											}
											++this->num_tildevar_cuts;
											//return true;


										}


									}
							}

							if (j != end_node
								&&
								curr_sol[v_t_x_kij(k, i, j)] > curr_sol[v_t_y_kn(k, j)] + cEPS)
							{
								*cut_added = true;

								for (int kk = 0; kk < this->K; ++kk)
									//		for (int kk = k; kk <= k; ++kk)

								{
									// Add the cuts for j, x and y and for each k
									rhs[0] = 0;
									sense[0] = 'L';
									nzc = 0;
									lhs_v = 0;

									matind[nzc] = v_t_x_kij(kk, i, j);
									matval[nzc] = 1;
									lhs_v += curr_sol[v_t_x_kij(kk, i, j)];
									++nzc;

									matind[nzc] = v_t_y_kn(kk, j);
									matval[nzc] = -1;
									lhs_v -= curr_sol[v_t_y_kn(kk, j)];
									++nzc;


									if (lhs_v > rhs[0] + cEPS)
									{

									//	++this->num_cut_curr_iter;


										status = CPXcallbackaddusercuts(context, 1, nzc, rhs,
											sense, matbeg, matind.data(), matval.data(), force, local);
										if (status)
										{
											checkCPXstatus(status);
											//fprintf(stdout, "Cut Logic!\n")
											return false;
										}
										++this->num_tildevar_cuts;
									//	return true;

									}
								}
							}
						}

						/*	if (*cut_added)
							{
								return true;
							}*/

					}
	//}
	return true;
}

bool ROPEUsolver::separate_McCormick_x_alpha(CPXCALLBACKCONTEXTptr context, CPXLONG contextid, bool* cut_added)
{
	int status = 0;
	double val = 0;
	double obj[1], rhs[1];
	int matbeg[1], nzc;
	char sense[1];
	int force[1];
	int local[1];
	force[0] = CPX_USECUT_FORCE;
	local[0] = 0;
	matbeg[0] = 0;
	double lhs_v = 0;


	int i_star, j_star, k_star, mc_star;
	double max_viol = 0;
	/*We check for he most violated*/
	for (int k = 0; k < this->K; ++k)
		for (int i = start_node; i < end_node; ++i)
			for (int j = i + 1; j <= end_node; ++j)
				if (
					fabs(curr_sol[v_alpha_k[k]] * curr_sol[v_x_kij(k, i, j)] - curr_sol[v_t_x_kij(k, i, j)]) > mEPS
					)
				{
					/*Mc 1*/
					lhs_v = 0;
					rhs[0] = 0;
					lhs_v += curr_sol[v_t_x_kij(k, i, j)];
					lhs_v -= curr_sol[v_x_kij(k, i, j)] * this->ub_alpha_k[k];

					if (lhs_v - rhs[0] >= max_viol)
					{
						max_viol = lhs_v - rhs[0];
						i_star = i;
						j_star = j;
						k_star = k;
						mc_star = 1;
					}

					/*Mc 2*/
					rhs[0] = -this->lb_alpha_k[k];
					lhs_v = 0;
					lhs_v += curr_sol[v_t_x_kij(k, i, j)];
					lhs_v -= curr_sol[v_x_kij(k, i, j)] * this->lb_alpha_k[k];
					lhs_v -= curr_sol[v_alpha_k[k]];

					if (lhs_v - rhs[0] >= max_viol)
					{
						max_viol = lhs_v - rhs[0];
						i_star = i;
						j_star = j;
						k_star = k;
						mc_star = 2;
					}

					/*Mc 3*/
					rhs[0] = -this->ub_alpha_k[k];
					lhs_v = 0;

					lhs_v += curr_sol[v_t_x_kij(k, i, j)];
					lhs_v -= curr_sol[v_alpha_k[k]];
					lhs_v -= curr_sol[v_x_kij(k, i, j)] * ub_alpha_k[k];
					//	if (lhs_v < rhs[0] - mEPS)
					if (rhs[0] - lhs_v >= max_viol)
					{
						max_viol = rhs[0] - lhs_v;
						i_star = i;
						j_star = j;
						k_star = k;
						mc_star = 3;
					}

					/*Mc 4*/
					rhs[0] = 0;
					lhs_v = 0;

					lhs_v += curr_sol[v_t_x_kij(k, i, j)];
					lhs_v -= curr_sol[v_x_kij(k, i, j)] * this->lb_alpha_k[k];

					//if (lhs_v < rhs[0] - mEPS)
					if (rhs[0] - lhs_v >= max_viol)
					{
						max_viol = rhs[0] - lhs_v;
						i_star = i;
						j_star = j;
						k_star = k;
						mc_star = 4;
					}

				}

	if (max_viol > cEPS)
	{
		*cut_added = true;
		switch (mc_star)
		{
		case 1:
			nzc = 0;
			rhs[0] = 0;
			sense[0] = 'L';
			lhs_v = 0;

			matind[nzc] = v_t_x_kij(k_star, i_star, j_star);
			matval[nzc] = 1;
			++nzc;
			lhs_v += curr_sol[v_t_x_kij(k_star, i_star, j_star)];

			matind[nzc] = v_x_kij(k_star, i_star, j_star);
			matval[nzc] = -this->ub_alpha_k[k_star];
			++nzc;
			lhs_v -= curr_sol[v_x_kij(k_star, i_star, j_star)] * this->ub_alpha_k[k_star];
			break;

		case 2:
			nzc = 0;
			rhs[0] = -this->lb_alpha_k[k_star];
			sense[0] = 'L';

			matind[nzc] = v_t_x_kij(k_star, i_star, j_star);
			matval[nzc] = 1;
			++nzc;

			matind[nzc] = v_x_kij(k_star, i_star, j_star);
			matval[nzc] = -this->lb_alpha_k[k_star];
			++nzc;

			matind[nzc] = v_alpha_k[k_star];
			matval[nzc] = -1;
			++nzc;
			break;

		case 3:
			nzc = 0;
			rhs[0] = -this->ub_alpha_k[k_star];
			sense[0] = 'G';

			matind[nzc] = v_t_x_kij(k_star, i_star, j_star);
			matval[nzc] = 1;
			++nzc;

			matind[nzc] = v_alpha_k[k_star];
			matval[nzc] = -1;
			++nzc;

			matind[nzc] = v_x_kij(k_star, i_star, j_star);
			matval[nzc] = -this->ub_alpha_k[k_star];
			++nzc;
			break;

		case 4:
			nzc = 0;
			rhs[0] = 0;
			sense[0] = 'G';

			matind[nzc] = v_t_x_kij(k_star, i_star, j_star);
			matval[nzc] = 1;
			++nzc;

			matind[nzc] = v_x_kij(k_star, i_star, j_star);
			matval[nzc] = -this->lb_alpha_k[k_star];
			++nzc;
			break;

		default:
			break;
		}

		if (contextid == CPX_CALLBACKCONTEXT_RELAXATION)
		{

			status = CPXcallbackaddusercuts(context, 1, nzc, rhs,
				sense, matbeg, matind.data(), matval.data(), force, local);
			if (status)checkCPXstatus(status);

			++this->num_cut_curr_iter;

		}
		else
			if (contextid == CPX_CALLBACKCONTEXT_CANDIDATE)
			{

				status = CPXcallbackrejectcandidate(context, 1, nzc, rhs,
					sense, matbeg, matind.data(), matval.data());
				if (status)checkCPXstatus(status);
			}
		++this->num_cut_curr_iter;
		++this->num_dyna_mccormik;

	}



	//{
		//	if (false)
		//	{
		//		for (int k = 0; k < this->K; ++k)
		//			for (int i = start_node; i < end_node; ++i)
		//				for (int j = i + 1; j <= end_node; ++j)
		//					if (
		//						fabs(curr_sol[v_alpha_k[k]] * curr_sol[v_x_kij(k, i, j)] - curr_sol[v_t_x_kij(k, i, j)]) > mEPS
		//						)
		//					{
		//
		//
		//						/*Here we check for violation*/
		//
		//
		//						/*McCormick (1)*/
		//						nzc = 0;
		//						rhs[0] = 0;
		//						sense[0] = 'L';
		//						lhs_v = 0;
		//
		//						matind[nzc] = v_t_x_kij(k, i, j);
		//						matval[nzc] = 1;
		//						++nzc;
		//						lhs_v += curr_sol[v_t_x_kij(k, i, j)];
		//
		//						matind[nzc] = v_x_kij(k, i, j);
		//						matval[nzc] = -this->ub_alpha_k[k];
		//						++nzc;
		//						lhs_v -= curr_sol[v_x_kij(k, i, j)] * this->ub_alpha_k[k];
		//
		//						if (lhs_v > rhs[0] + mEPS
		//							&&
		//							nzc > 0)
		//						{
		//							*cut_added = true;
		//
		//
		//							if (contextid == CPX_CALLBACKCONTEXT_RELAXATION)
		//							{
		//
		//								status = CPXcallbackaddusercuts(context, 1, nzc, rhs,
		//									sense, matbeg, matind.data(), matval.data(), force, local);
		//								if (status)checkCPXstatus(status);
		//
		//								++this->num_cut_curr_iter;
		//							}
		//							else
		//								if (contextid == CPX_CALLBACKCONTEXT_CANDIDATE)
		//								{
		//
		//									status = CPXcallbackrejectcandidate(context, 1, nzc, rhs,
		//										sense, matbeg, matind.data(), matval.data());
		//									if (status)checkCPXstatus(status);
		//								}
		//
		//
		//
		//							/*	status = CPXcallbackaddusercuts(context, 1, nzc, rhs,
		//									sense, matbeg, matind.data(), matval.data(), force, local);*/
		//
		//							++this->num_dyna_mccormik;
		//
		//
		//							if (status)
		//							{
		//								checkCPXstatus(status);
		//								return false;
		//							}
		//						}
		//
		//						/*McCormick (2)*/
		//						nzc = 0;
		//						rhs[0] = -this->lb_alpha_k[k];
		//						sense[0] = 'L';
		//						lhs_v = 0;
		//
		//						matind[nzc] = v_t_x_kij(k, i, j);
		//						matval[nzc] = 1;
		//						++nzc;
		//						lhs_v += curr_sol[v_t_x_kij(k, i, j)];
		//
		//						matind[nzc] = v_x_kij(k, i, j);
		//						matval[nzc] = -this->lb_alpha_k[k];
		//						++nzc;
		//						lhs_v -= curr_sol[v_x_kij(k, i, j)] * this->lb_alpha_k[k];
		//
		//						matind[nzc] = v_alpha_k[k];
		//						matval[nzc] = -1;
		//						++nzc;
		//						lhs_v -= curr_sol[v_alpha_k[k]];
		//
		//						if (lhs_v > rhs[0] + mEPS
		//							&&
		//							nzc > 0)
		//						{
		//							*cut_added = true;
		//
		//
		//
		//
		//							if (contextid == CPX_CALLBACKCONTEXT_RELAXATION)
		//							{
		//
		//								status = CPXcallbackaddusercuts(context, 1, nzc, rhs,
		//									sense, matbeg, matind.data(), matval.data(), force, local);
		//								if (status)checkCPXstatus(status);
		//
		//								++this->num_cut_curr_iter;
		//
		//							}
		//							else
		//								if (contextid == CPX_CALLBACKCONTEXT_CANDIDATE)
		//								{
		//
		//									status = CPXcallbackrejectcandidate(context, 1, nzc, rhs,
		//										sense, matbeg, matind.data(), matval.data());
		//									if (status)checkCPXstatus(status);
		//								}
		//
		//
		//
		//
		//
		//							//status = CPXcallbackaddusercuts(context, 1, nzc, rhs,
		//							//	sense, matbeg, matind.data(), matval.data(), force, local);
		//
		//							++this->num_dyna_mccormik;
		//
		//
		//							if (status)
		//							{
		//								checkCPXstatus(status);
		//								return false;
		//							}
		//						}
		//
		//						/*McCormick (3)*/
		//						nzc = 0;
		//						rhs[0] = -this->ub_alpha_k[k];
		//						sense[0] = 'G';
		//						lhs_v = 0;
		//
		//						matind[nzc] = v_t_x_kij(k, i, j);
		//						matval[nzc] = 1;
		//						++nzc;
		//						lhs_v += curr_sol[v_t_x_kij(k, i, j)];
		//
		//						matind[nzc] = v_alpha_k[k];
		//						matval[nzc] = -1;
		//						++nzc;
		//						lhs_v -= curr_sol[v_alpha_k[k]];
		//
		//						matind[nzc] = v_x_kij(k, i, j);
		//						matval[nzc] = -this->ub_alpha_k[k];
		//						++nzc;
		//						lhs_v -= curr_sol[v_x_kij(k, i, j)] * ub_alpha_k[k];
		//
		//						if (lhs_v < rhs[0] - mEPS
		//							&&
		//							nzc > 0)
		//						{
		//							*cut_added = true;
		//
		//
		//
		//							if (contextid == CPX_CALLBACKCONTEXT_RELAXATION)
		//							{
		//
		//								status = CPXcallbackaddusercuts(context, 1, nzc, rhs,
		//									sense, matbeg, matind.data(), matval.data(), force, local);
		//								if (status)checkCPXstatus(status);
		//
		//								++this->num_cut_curr_iter;
		//
		//							}
		//							else
		//								if (contextid == CPX_CALLBACKCONTEXT_CANDIDATE)
		//								{
		//
		//									status = CPXcallbackrejectcandidate(context, 1, nzc, rhs,
		//										sense, matbeg, matind.data(), matval.data());
		//									if (status)checkCPXstatus(status);
		//								}
		//
		//
		//
		//							/*	status = CPXcallbackaddusercuts(context, 1, nzc, rhs,
		//									sense, matbeg, matind.data(), matval.data(), force, local);*/
		//
		//							++this->num_dyna_mccormik;
		//
		//
		//							if (status)
		//							{
		//								checkCPXstatus(status);
		//								return false;
		//							}
		//						}
		//
		//						/*McCormick (4)*/
		//						nzc = 0;
		//						rhs[0] = 0;
		//						sense[0] = 'G';
		//						lhs_v = 0;
		//
		//						matind[nzc] = v_t_x_kij(k, i, j);
		//						matval[nzc] = 1;
		//						++nzc;
		//						lhs_v += curr_sol[v_t_x_kij(k, i, j)];
		//
		//						matind[nzc] = v_x_kij(k, i, j);
		//						matval[nzc] = -this->lb_alpha_k[k];
		//						++nzc;
		//						lhs_v -= curr_sol[v_x_kij(k, i, j)] * this->lb_alpha_k[k];
		//
		//						if (lhs_v < rhs[0] - mEPS
		//							&&
		//							nzc > 0)
		//						{
		//							*cut_added = true;
		//
		//
		//
		//							if (contextid == CPX_CALLBACKCONTEXT_RELAXATION)
		//							{
		//
		//								status = CPXcallbackaddusercuts(context, 1, nzc, rhs,
		//									sense, matbeg, matind.data(), matval.data(), force, local);
		//								if (status)checkCPXstatus(status);
		//
		//								++this->num_cut_curr_iter;
		//
		//							}
		//							else
		//								if (contextid == CPX_CALLBACKCONTEXT_CANDIDATE)
		//								{
		//
		//									status = CPXcallbackrejectcandidate(context, 1, nzc, rhs,
		//										sense, matbeg, matind.data(), matval.data());
		//									if (status)checkCPXstatus(status);
		//								}
		//
		//
		//							/*	status = CPXcallbackaddusercuts(context, 1, nzc, rhs,
		//									sense, matbeg, matind.data(), matval.data(), force, local);*/
		//
		//							++this->num_dyna_mccormik;
		//
		//
		//							if (status)
		//							{
		//								checkCPXstatus(status);
		//								return false;
		//							}
		//						}
		//						//return true;
		//
		//					}
		//
		//}
	//}
					return true;
}



bool ROPEUsolver::define_uncertainty_set_type_1()
{
	int node, col, row;

	row = 0;

		/*\xi in [0,1] and \sum_i \xi <= B */
			 // First the upper bound xi <= 1: therefore we stop at Nxi;
		for (node = this->start_node; node <= this->end_node; ++node)
			if (this->node_is_profitable[node])
			{

				this->A(row, node) = 1;
				this->b[row] = 1;
				++row;


			}

	for (node = this->start_node; node <= this->end_node; ++node)
		if (this->node_is_profitable[node])
		{

			this->A(row, node) = -1;
			this->b[row] = 0;
			++row;


		}

	/*Last constraints about the budget.*/
	double tot_det_profit = 0;
	//
	for (node = this->start_node; node <= this->end_node; ++node)
		if (this->node_is_profitable[node])
		{
			tot_det_profit += this->node_det_profit[node];

			this->A(row, node) = 1;
		}
	this->b[row] = this->uncertainty_parameter;
	++row;
	this->num_us_constr = row;


	/*Now we define the US cost*/
	//- y(p - p*0.5*\xi) ->  y(-p + p*0.5*\xi)
	for (int i = start_node; i <= this->end_node; ++i)
	{
		this->best_scenario_profit_i[i] = this->node_det_profit[i]; /*\xi = 0 for each entry*/

		this->t_y_obj_coeff_val[i] =- this->node_det_profit[i]; // Coefficent without xi

		this->t_y_uc_coeff_val[i] = this->node_det_profit[i] * 
			0.5
			; // Coefficient with xi: the proft can be reduced at most by 50%

	}



	return true;
}

bool ROPEUsolver::define_uncertainty_set_type_2()
{
	int node, col, row;
	double fraction = this->uncertainty_parameter;

	row = 0;

		/*\xi in [0,1] and \sum_i \xi = 1 */
	/*y(P\xi)*/
			 // First the upper bound xi <= 1: therefore we stop at Nxi;
		for (node = this->start_node; node <= this->end_node; ++node)
			if (this->node_is_profitable[node])
			{

				this->A(row, node) = 1;
				//this->b[row] = 1;
				this->b[row] =fraction;
				++row;


			}

		for (node = this->start_node; node <= this->end_node; ++node)
			if (this->node_is_profitable[node])
			{

				this->A(row, node) = -1;
				this->b[row] = 0;
				++row;


			}


		/*Last 2 constraint about the sum of the total is equal to 1.*/
		//
		for (node = this->start_node; node <= this->end_node; ++node)
			if (this->node_is_profitable[node])
			{

				this->A(row, node) = 1;
			}

		this->b[row] = (1 + 1e-7);
		//this->b[row] = 1;
		++row;

		for (node = this->start_node; node <= this->end_node; ++node)
			if (this->node_is_profitable[node])
			{

				this->A(row, node) = -1;
			}
		//this->b[row] = -1;
		this->b[row] = - (1 - 1e-7);
		++row;


		this->num_us_constr = row;
		/*Now we define the US cost*/
		// -y(P\xi) where P is the total profit -> y(-P\xi)
		for (int i = start_node; i <= this->end_node; ++i)
		{
			//this->best_scenario_profit_i[i] = this->tot_deterministic_profit; /*\xi = 1 for each entry: then eveyrone has the max profit*/
			this->best_scenario_profit_i[i] = 1;
			
			this->t_y_obj_coeff_val[i] =- 0; // Coefficent without xi

			//this->t_y_uc_coeff_val[i] =- this->tot_deterministic_profit; // Coefficient with xi
			this->t_y_uc_coeff_val[i] = -1; /*100%*/
		}

	

	return true;
}

bool ROPEUsolver::define_uncertainty_set_type_3()
{
	int node, col, row;
	double fraction = this->uncertainty_parameter;

	row = 0;

	/*\xi in [-1,1] and \sum_i \xi =0 */
	//	//- y(p - p*0.5*\xi) --> y(-p + p*0.5*\xi)

		 // First the upper bound xi <= 1: therefore we stop at Nxi;
	for (node = this->start_node; node <= this->end_node; ++node)
		if (this->node_is_profitable[node])
		{

			this->A(row, node) = 1;
			this->b[row] = 1;
			++row;


		}

	for (node = this->start_node; node <= this->end_node; ++node)
		if (this->node_is_profitable[node])
		{

			this->A(row, node) = -1;
			this->b[row] = 1;
			++row;


		}

	/*Last 2 constraints about the budget.*/
	//
	for (node = this->start_node; node <= this->end_node; ++node)
		if (this->node_is_profitable[node])
		{

			this->A(row, node) = 1;
		}
	//this->b[row] = this->uncertainty_parameter;
	this->b[row] =0;
	++row;

	/*Last constraints about the budget.*/
	//
	for (node = this->start_node; node <= this->end_node; ++node)
		if (this->node_is_profitable[node])
		{

			this->A(row, node) = -1;
		}
	//this->b[row] = -this->uncertainty_parameter;
	this->b[row] =0;
	++row;



	this->num_us_constr = row;


	/*Now we define the US cost*/
	// y(p - p*0.5*\xi)
	for (int i = start_node; i <= this->end_node; ++i)
	{
		this->best_scenario_profit_i[i] = this->node_det_profit[i] + (this->node_det_profit[i]* fraction); /*\xi = -1 for each entry*/

		this->t_y_obj_coeff_val[i] = - this->node_det_profit[i]; // Coefficent without xi

		this->t_y_uc_coeff_val[i] = this->node_det_profit[i] * fraction;

	//	this->t_y_uc_coeff_val[i] = this->node_det_profit[i] * 0.5; // Coefficient with xi: the proft can be reduced at most by 50%

	}



	return true;
}

bool ROPEUsolver::define_uncertainty_set_type_4()
{
	int node, col, row;


	/*\xi in [(1-teta) p, (1+teta)p] and \sum_i \xi = tot_profit*/
	//	//- y(*\xi)
	double theta = this->uncertainty_parameter;
	
	/*Now we define the US cost*/
// y(\xi)
	for (int i = start_node; i <= this->end_node; ++i)
	{
		this->best_scenario_profit_i[i] = this->node_det_profit[i]  + this->node_det_profit[i] * theta; /*\xi at max for each node*/

		this->t_y_obj_coeff_val[i] = 0; // Coefficent without xi

		this->t_y_uc_coeff_val[i] = -1; // Coefficient with xi: the profit is exaclyt xi now

	}
	
	row = 0;
		 // First the upper bound xi <= (1+theta)*p : therefore we stop at Nxi;
	for (node = this->start_node; node <= this->end_node; ++node)
		if (this->node_is_profitable[node])
		{

			this->A(row, node) = 1;
			this->b[row] = (1.00 + theta) * this->node_det_profit[node];
			++row;


		}
	// First the upper bound xi >= (1-theta)*p : therefore we stop at Nxi;

	for (node = this->start_node; node <= this->end_node; ++node)
		if (this->node_is_profitable[node])
		{

			this->A(row, node) = -1;
			this->b[row] =- (1 - theta) * this->node_det_profit[node];
			++row;


		}

	/*Last 2 constraints about the budget.*/
	//
	for (node = this->start_node; node <= this->end_node; ++node)
		if (this->node_is_profitable[node])
		{

			this->A(row, node) = 1;
		}
	this->b[row] = this->tot_deterministic_profit;
	++row;

	/*Last constraints about the budget.*/
	//
	for (node = this->start_node; node <= this->end_node; ++node)
		if (this->node_is_profitable[node])
		{

			this->A(row, node) = -1;
		}
	this->b[row] = -this->tot_deterministic_profit;
	++row;
	this->num_us_constr = row;



	return true;
}

bool ROPEUsolver::print_debug_curr_sol()
{

	/*Printing interesting information*/
	std::cout << "Instance_Name" << "	" << "Max_Tour_Duration" << "	" << "num_K" <<
		"	" << "Best_Upper_Bound" << "	" << "Best_Lower_Bound" << "	" << "Cplex_Opt_Gap" << "	" << "Comp_Time" <<"	"<< "Exact_W_Evaluation" <<
		std::endl;
	std::cout << this->inst_name << "	" << this->max_dur << "	" << this->K <<
		"	" << this->best_ub << "	" << this->best_lb << "	" << this->cplex_gap << "	" << this->comput_time <<"	" <<this->exact_w_evaluation <<std::endl;



	//for (int k = 0; k < this->K; ++k)
	//{
	//	for (int i = this->first_pnode; i <= this->last_pnode; ++i)
	//		if (curr_sol[this->v_y_kn(k,i)] > mEPS)
	//		{
	//			std::cout << "y_k" << k << "_i" << i << " = 1" << std::endl;
	//		}

	//	for(int i = start_node; i <= end_node; ++i)
	//		for (int j = start_node; j <= end_node; ++j)
	//			if (this->edge_exists(i, j)
	//				&&
	//				curr_sol[this->v_x_kij(k, i, j)] > mEPS
	//				)
	//			{
	//				std::cout << "x_k" << k << "_i" << i <<"_j"<<j <<" = 1" << std::endl;

	//			}

	//}


	return true;
}

bool ROPEUsolver::print_solution_file()
{

	unsigned short num_discovery = 0;
	double discovery_cost = 0;
	double minus_profi = 0;

	if (init_out_file())
	{
		if (out_file == NULL)
		{

			printf("Error writing solution file!\n");
			return false;
		}

		fprintf(out_file, "++++ Instance Info ++++\n");
		fprintf(out_file, "Instance name: %s\n", inst_name);
		fprintf(out_file, "Network name: %s\n", network_name);

		fprintf(out_file, "Number of nodes: %d\n", num_nodes);
		fprintf(out_file, "Number of profitable nodes: %d\n", num_prof_nodes);
		fprintf(out_file, "Max number of discovery: %d\n", max_num_disc);
		fprintf(out_file, "Max duration: %d\n", max_dur);
		fprintf(out_file, "Alrijne Collecting Time: %lf\n", alrijne_node_collecting_time);

		fprintf(out_file, "++++ Uncertainty  ++++\n");
		fprintf(out_file, "Number of K-policy: %d\n", K);
		fprintf(out_file, "Uncertainty Set Type: %d\n", this->uncertainty_set_type);
		fprintf(out_file, "Uncertainty Set Parameter: %lf\n", this->uncertainty_parameter);
		fprintf(out_file, "Total Deterministic Profit in the Network: %lf\n", this->tot_deterministic_profit);



		fprintf(out_file, "++++ Solution of the Model ++++\n");
		fprintf(out_file, "-----------------------------------------------");
		fprintf(out_file, "\nObjective:\n");
		fprintf(out_file, "Best Lower Bound: %lf\n", best_lb);
		fprintf(out_file, "Best upper bound: %lf\n", best_ub);
		fprintf(out_file, "Elapsed time: %lf\n", this->comput_time);
		fprintf(out_file, "Cplex Gap: %lf\n", this->cplex_gap);
		fprintf(out_file, "Exact W evaluation: %lf\n", this->exact_w_evaluation);

		
		fprintf(out_file, "B&C nodes: %d\n", bandc_nodes);

		fprintf(out_file, "Number Added VI on integer variables: %d\n", this->num_dvar_cuts);
		fprintf(out_file, "Number Added VI on tilde variables: %d\n", this->num_tildevar_cuts);
		fprintf(out_file, "Number Added McCormick VI: %d\n", this->num_dyna_mccormik);




		fprintf(out_file, "\nVariables:\n");

		fprintf(out_file, "\nDiscoveries: w_i (i is a profitable node). \n");
		for (int i = start_node + 1; i < end_node; ++i)
			if (curr_sol[v_w_n[i]] > mEPS)
			{
				fprintf(out_file, " w_[node:%d] = %lf\n", i, curr_sol[v_w_n[i]]);
				++num_discovery;
				discovery_cost += this->disc_cost_n[i];
			}

		fprintf(out_file, "Beta Variables. \n");
		for (unsigned short r = 0; r < this->num_us_constr; ++r)
			if (std::abs(curr_sol[v_beta_l[r]]) > mEPS)
		{
			fprintf(out_file, "beta_l[%d] = %lf\n", r, curr_sol[v_beta_l[r]]);

		}

		fprintf(out_file, "K-policies. \n");
		for (int k = 0; k < K; ++k)
		{
			fprintf(out_file, "Visited nodes:\n");
			for (int i = start_node + 1; i < end_node; ++i)
				if (curr_sol[v_y_kn(k, i)] > mEPS)
				{
					fprintf(out_file, "y_k[%d]_i[%d] = %.0f\n", k, i, curr_sol[v_y_kn(k, i)]);
					fprintf(out_file, "tilde_y_k[%d]_i[%d] = %.2f \n", k, i, curr_sol[v_t_y_kn(k, i)]);

							//fprintf(out_file, "ov_betak_[%d]_[%d] = %.lf:\n", k,i, curr_sol[this->v_beta_l]);
						//	fprintf(out_file, "gamma_[%d]_[%d] = %.lf:\n", k, i, curr_sol[g_v_gammak_ki[k][i]]);
							//fprintf(out_file, "tilde_gamma_[%d]_[%d] = %.lf:\n", k, i, curr_sol[g_vt_gammak_ki[k][i]]);
				}

			fprintf(out_file, "Duals:\n");
			fprintf(out_file, "alpha_k[%d] = %lf\n", k, curr_sol[v_alpha_k[k]]);
				for (unsigned short n = this->first_pnode; n <= this->last_pnode; ++n)
					if(std::abs(curr_sol[v_gamma_kn(k, n)]) > mEPS)
				{
					fprintf(out_file, "gamma_k[%d]_n[%d] = %lf\n", k,n ,curr_sol[v_gamma_kn(k, n)]);

				}
				for (unsigned short r = 0; r < this->num_us_constr; ++r)
					if (std::abs(curr_sol[v_beta_kl(k, r)]) > mEPS)
				{
					fprintf(out_file, "ov_betak_[%d]_[%d] = %lf\n", k, r, curr_sol[v_beta_kl(k,r)]);

				}



			fprintf(out_file, "Policy-%d description\n", k);
			int curr_node = start_node;
			int prev_node = -1;
			double route_length = 0;
			fprintf(out_file, "<%d", curr_node);

			while (curr_node != end_node)
			{
				for (int j = start_node; j <= end_node; ++j)
					if (j != curr_node
						&&
						j != prev_node
						&&
						curr_sol[v_x_kij(k, std::min(curr_node, j), std::max(curr_node, j))] > mEPS)
					{
						route_length += t_ij(std::min(curr_node, j), std::max(curr_node, j));
						//[imin(curr_node, j)] [imax(curr_node, j)] ;

						prev_node = curr_node;
						curr_node = j;

						break;
					}
				fprintf(out_file, "-%d", curr_node);
			}
			fprintf(out_file, "> : Length = %lf\n", route_length);
			fprintf(out_file, "_________________________________________________\n");
		}


			fprintf(out_file, "Legenda_Model:");
			fprintf(out_file, "Instance_Name\t");
			fprintf(out_file, "Network_Name\t");

			fprintf(out_file, "Num_K\t");

			fprintf(out_file, "Num_Nodes\t");
			fprintf(out_file, "Num_Prof_Nodes\t");
			fprintf(out_file, "Max_Tour_Duration\t");


			fprintf(out_file, "Max_Num_Discoveries\t");
			fprintf(out_file, "Max_Discoveries_Fraction\t");
			fprintf(out_file, "Best_LB\t");
			fprintf(out_file, "Best_Ub\t");
			fprintf(out_file, "Cplex_Gap\t");
			fprintf(out_file, "Root_Node_LB\t");
			fprintf(out_file, "Elapsed_Time\t");
			fprintf(out_file, "Num_VI_int_var\t");
			fprintf(out_file, "Num_VI_tilde_var\t");
			fprintf(out_file, "Num_Dyn_McCormick\t");
			fprintf(out_file, "Num_BandC_Nodes\t");
			fprintf(out_file, "Number_Discoveries\t");
			fprintf(out_file, "Profit_cost\t");
			fprintf(out_file, "Discovery_cost\t");
			fprintf(out_file, "US_Type\t");
			fprintf(out_file, "US_Param\t");
			fprintf(out_file, "US_Tot_Det_Profit\t");
			fprintf(out_file, "Exact_W_Evaluation\n");



			fprintf(out_file, "Table_Model:");
			fprintf(out_file, "%s\t", inst_name);
			fprintf(out_file, "%s\t", network_name);
			fprintf(out_file, "%d\t", K);
			fprintf(out_file, "%d\t", this->num_nodes);
			fprintf(out_file, "%d\t", this->num_prof_nodes);
			fprintf(out_file, "%d\t", this->max_dur);
			fprintf(out_file, "%d\t", this->max_num_disc);
			fprintf(out_file, "%lf\t", this->max_number_discovery_fraction);
			fprintf(out_file, "%lf\t", this->best_lb);
			fprintf(out_file, "%lf\t", this->best_ub);
			fprintf(out_file, "%lf\t", this->cplex_gap);
			fprintf(out_file, "%lf\t", this->root_node_lb);
			fprintf(out_file, "%lf\t", this->comput_time);
			fprintf(out_file, "%d\t", this->num_dvar_cuts);
			fprintf(out_file, "%d\t", this->num_tildevar_cuts);
			fprintf(out_file, "%d\t", this->num_dyna_mccormik);
			fprintf(out_file, "%d\t", this->bandc_nodes);
			fprintf(out_file, "%d\t", num_discovery);
			fprintf(out_file, "%lf\t", this->best_ub - discovery_cost);
			fprintf(out_file, "%lf\t", discovery_cost);
			fprintf(out_file, "%d\t", uncertainty_set_type);
			fprintf(out_file, "%lf\t",uncertainty_parameter);
			fprintf(out_file, "%lf\t", this->tot_deterministic_profit);
			fprintf(out_file, "%lf\n", this->exact_w_evaluation);


		
	}
	return true;
}

bool ROPEUsolver::print_solution_latex()
{
	std::ofstream latex_out_file;

	char locname[1024];
	char ss[1024];
	snprintf(locname, sizeof(locname), "%s%s%s", results_folder, "\\Sol_Latex_Kapp_", inst_name);
	snprintf(ss, sizeof(ss), "%s%s", locname, ".txt");
	
	latex_out_file.open(ss);
	//latex_out_file = fopen(ss, "w");

	if (!latex_out_file.is_open())
	{
		printf("Error creating output file.\n");
		return false;
	}
	latex_out_file << "\\begin{tikzpicture}" << std::endl;

	latex_out_file << "\\begin{axis} [legend pos=south east," << std::endl;
	latex_out_file << "scatter/classes={%" << std::endl;
	latex_out_file << "	start={mark=triangle*}," << std::endl;
	latex_out_file << "	end={mark=triangle},"<< std::endl;
	latex_out_file << "discovery_on={mark=square*},"<< std::endl;
	latex_out_file << "discovery_off={mark=o}}]" << std::endl;
	
	latex_out_file << "\\addplot[scatter, only marks,  scatter src = explicit symbolic, forget plot] table[meta = label]{"<< std::endl;
	latex_out_file << "x     y      label" << std::endl;

	/*
	Here the nodes
	*/

	/*First start and end*/
	latex_out_file <<std::to_string(this->node_x_coord[start_node]) <<"	"<< std::to_string(this->node_y_coord[start_node])<<"	"<<"start"<<std::endl;
	latex_out_file << std::to_string(this->node_x_coord[end_node]) << "	" << std::to_string(this->node_y_coord[end_node]) << "	" << "end" << std::endl;
	/*Profitable nodes*/
	for (int n = first_pnode; n <= last_pnode; ++n)
	{
		if (curr_sol[this->v_w_n[n]] > EPS)
		{
			latex_out_file << std::to_string(this->node_x_coord[n]) << "	" << std::to_string(this->node_y_coord[n]) << "	" << "discovery_on" << std::endl;

		}else
		{
			latex_out_file << std::to_string(this->node_x_coord[n]) << "	" << std::to_string(this->node_y_coord[n]) << "	" << "discovery_off" << std::endl;

		}

	}


	latex_out_file << "};"<< std::endl;


	
	/*Now all the policy*/
	for (int k = 0; k < this->K; ++k)
	{
		if (k == 0)
		{
			latex_out_file << "	\\addplot[thick, dashed, mark=none] coordinates {" << std::endl;
		}
		if (k == 1)
		{
			latex_out_file << "	\\addplot[thick,dash dot,mark=none, color=blue] coordinates {" << std::endl;
		}
		if (k == 2)
		{
			latex_out_file << "	\\addplot[mark=none, color=green] coordinates {" << std::endl;
		}
		if (k == 3)
		{
			latex_out_file << "	\\addplot[thick,dotted,color=red,mark=none] coordinates {" << std::endl;
		}
		if (k > 3)
		{
			latex_out_file << "	\\addplot[thick,dotted,color=red,mark=none] coordinates {" << std::endl;
		}

		/*
		Here used edges
		*/

		int curr_node = start_node;
		int prev_node = -1;

		while (curr_node != end_node)
		{
			for (int j = start_node; j <= end_node; ++j)
				if (j != curr_node
								&&
								j != prev_node
								&&
								curr_sol[v_x_kij(k, std::min(curr_node, j), std::max(curr_node, j))] > EPS)
			{
					prev_node = curr_node;
					curr_node = j;
					break;
			}
			latex_out_file << "(" << this->node_x_coord[prev_node] << ", " << this->node_y_coord[prev_node] << ")" << std::endl;
			
			if (curr_node == end_node)
			{
				latex_out_file << "(" << this->node_x_coord[curr_node] << ", " << this->node_y_coord[curr_node] << ")" << std::endl;
			}
		}

		latex_out_file << "};" << std::endl;

	}
	latex_out_file << "\\addlegendentry{{\\scriptsize k = 1}}" << std::endl;
	latex_out_file << "\\addlegendentry{{\\scriptsize k = 2}}" << std::endl;
	latex_out_file << "\\addlegendentry{{\\scriptsize k = 3}}" << std::endl;
	latex_out_file << "\\addlegendentry{{\\scriptsize k = 4}}" << std::endl;
	latex_out_file << "\\end{axis}" << std::endl;
	latex_out_file << "\\end{tikzpicture}" << std::endl;

	latex_out_file.close();

	return true;
}

bool ROPEUsolver::build_and_solve_model()
{
	bool is_ok = true;
	int status;
	int mip_status;
	CPXLONG      contextmask = 0;
	char errbuf[CPXMESSAGEBUFSIZE];
	int check_mip_stat;

	double start = 0;

	double end = 0;

	is_ok = init_cplex();
	if (!is_ok)
	{
		goto TERMINATE;
	}

	is_ok = set_cplex();
	if (!is_ok)
	{
		goto TERMINATE;
	}

	is_ok = define_alpha_bounds();
	if (!is_ok)
	{
		goto TERMINATE;
	}

	is_ok = add_deterministic_variables();
	if (!is_ok)
	{
		goto TERMINATE;
	}

	is_ok = add_dualization_variables();
	if (!is_ok)
	{
		goto TERMINATE;
	}

	is_ok = add_linearization_variables();
	if (!is_ok)
	{
		goto TERMINATE;
	}

	is_ok = add_deterministic_constraints();
	if (!is_ok)
	{
		goto TERMINATE;
	}

	is_ok = add_dualization_constraints();
	if (!is_ok)
	{
		goto TERMINATE;
	}

	if (this->alpha_symmetry_breaking)
	{
		is_ok = add_valid_inequalities_alpha_symmetries();
		if (!is_ok)
		{
			goto TERMINATE;
		}
	}
	is_ok = add_strengthened_linearization_constraints();
	//is_ok = add_linearization_constraints();
	if (!is_ok)
	{
		goto TERMINATE;
	}

	is_ok = add_objective_constraint();
	if (!is_ok)
	{
		goto TERMINATE;
	}

	if (best_scenario_cuts)
	{
		is_ok = add_valid_inequalities_objective_lower_bound_constraint();
		if (!is_ok)
		{
			goto TERMINATE;
		}
	}

	if (xtilde_cuts)
	{
		is_ok = add_RLT_valid_inequalties();
		if (!is_ok)
		{
			goto TERMINATE;
		}
	}

	contextmask |= CPX_CALLBACKCONTEXT_CANDIDATE;
	contextmask|= CPX_CALLBACKCONTEXT_RELAXATION;
	if (contextmask != 0) {
		//	 We are done and now we register our callback function. 
		status = CPXcallbacksetfunc(cpx_env, cpx_lp, contextmask, general_callback, this);
		if (status != 0) {
			fprintf(stderr, "Failed to add callback: %s\n",
				CPXgeterrorstring(cpx_env, status, errbuf));
			goto TERMINATE;
		}
	}



	status = CPXwriteprob(this->cpx_env, this->cpx_lp, "model.lp", "lp");
if (this->checkCPXstatus(status)) { is_ok = false; goto TERMINATE; };

	start = clock();
	status = CPXmipopt(this->cpx_env, this->cpx_lp);
	if (checkCPXstatus(status)) goto TERMINATE;
	end = clock();
	this->comput_time =(double)(end - start) / (double)CLK_TCK;


	check_mip_stat = CPXgetstat(cpx_env, cpx_lp);

		status = CPXgetx(cpx_env, cpx_lp, curr_sol.data(), 0, CPXgetnumcols(cpx_env, cpx_lp) - 1);
	if (checkCPXstatus(status)) goto TERMINATE;




	if (check_mip_stat != CPXMIP_TIME_LIM_INFEAS)
	{
		status = CPXgetobjval(cpx_env, cpx_lp, &this->best_ub);
		if (checkCPXstatus(status)) goto TERMINATE;
	}
	else
	{
		this->best_ub = -1;

	}


	status = CPXgetbestobjval(cpx_env, cpx_lp, &this->best_lb);
	if (checkCPXstatus(status)) goto TERMINATE;



	status = CPXgetmiprelgap(cpx_env, cpx_lp, &this->cplex_gap);
	if (checkCPXstatus(status)) goto TERMINATE;

	bandc_nodes = CPXgetnodecnt(cpx_env, cpx_lp);



	TERMINATE:
	return is_ok;
}

bool ROPEUsolver::get_w_values(std::vector<bool>& w_val)
{
	bool is_ok = true;

	w_val.assign(this->num_nodes, false);

	for (int i = start_node + 1; i < end_node; ++i)
		if (curr_sol[v_w_n[i]] > EPS)
		{
			w_val[i] = true;
		}
	


	return is_ok;
}

bool ROPEUsolver::read_distance_matrix(FILE* inp)
{

	//double coll_time_ij = ALRIJNE_COLL_TIME * 100.00;

//double coll_time_dep_i = ALRIJNE_COLL_TIME * 100.00 * 0.5;




	double coll_time_ij = alrijne_node_collecting_time * 100.00;

	double coll_time_dep_i = alrijne_node_collecting_time * 100.00 * 0.5;

	t_ij(this->start_node, this->end_node) = 0;
	this->edge_exists(this->start_node, this->end_node) = true;


	int i, j;




	for (i = 0; i < this->num_nodes - 1; i++)
	{


		for (j = 0; j < i; j++)
		{
			//Use lf format specifier, %c is for character
			if (!fscanf(inp, "%lf", &t_ij(j, i)))
			{
				std::cout << "Error while reading matrix in EXAsolver::read_distance_matrix(FILE* inp)" << std::endl;
				getchar();
				return false;
				//break;
			}
			else
			{
				this->edge_exists(j, i) = true;

				t_ij(j, i) = (int)(t_ij(j, i) * 100.00);
				t_ij(i, j) = t_ij(i, j);


				if (i != start_node && i != end_node && j != start_node && j != end_node)
				{
					t_ij(j, i) = t_ij(j, i) + coll_time_ij;
					t_ij(i, j) = t_ij(i, j) + coll_time_ij;
				}
				else {
					t_ij(j, i) = t_ij(j, i) + coll_time_dep_i;
					t_ij(i, j) = t_ij(i, j) + coll_time_dep_i;
				}



				if (j == 0)
				{
					t_ij(i, this->num_nodes - 1) = t_ij(j, i);
					t_ij(this->num_nodes - 1, i) = t_ij(j, i);

					this->edge_exists(i, this->num_nodes - 1) = true;
				}
			}
			// mat[i][j] -= '0'; 
			//printf("%lf\n", mat[i][j]); //Use lf format specifier, \n is for new line
		}


	}

	//Now profit
	this->node_det_profit[start_node] = 0;
	this->node_det_profit[end_node] = 0;
	for (int i = first_pnode; i <= last_pnode; ++i)
	{
		if (!fscanf(inp, "%lf", &this->node_det_profit[i]))
		{
			std::cout << "Error while reading profit in EXAsolver::read_distance_matrix(FILE* inp)" << std::endl;
			getchar();
			return false;
			//break;
		}
	}

	return true;
}

bool ROPEUsolver::read_instance(char* file, char* config_file)
{
	FILE* inp;
	errno_t err;
	char line[1024];
	char* ss;
	char cline[1024];

	bool read_matrix = false;


	if (config_file == NULL)
	{
		cfg = new Configuration::ConfigFile("config.cfg");
		std::cout << "Standard configuration file config.cfg loaded" << std::endl;
	}
	else
		cfg = new Configuration::ConfigFile(config_file);

	//cfg = new Configuration::ConfigFile("config.cfg");

	this->K = this->cfg->getValueOfKey<int>("NUM_K");
	this->decision_depende_discovery = this->cfg->getValueOfKey<bool>("DDID_VERSION");
	this->max_number_discovery_fraction = this->cfg->getValueOfKey<double>("MAX_NUM_DISCOVERY_FRACTION");
	this->uncertainty_set_type = this->cfg->getValueOfKey<int>("US_TYPE");
	this->uncertainty_parameter = this->cfg->getValueOfKey<double>("US_PARAM");

	read_matrix = this->cfg->getValueOfKey<bool>("READ_D_MATRIX");

	alrijne_node_collecting_time = cfg->getValueOfKey<double>("ALRIJNE_COLL_NODE_TIME");


	memset(inst_name, 0, sizeof(inst_name));

	char ignore[1024];
	if ((err = fopen_s(&inp, file, "r")) != 0) {
		printf(" Cannot open the input file \n");
		exit(1);
		return false;
	}
	int pos_l = -1;
	int pos_r = -1;
	for (int i = strlen(file); i >= 0; --i)
	{
		if (file[i] == '.' && pos_r == -1)
			pos_r = i;
		if (file[i] == '\\' && pos_l == -1)
			pos_l = i + 1;

		if (pos_l != -1 && pos_r != -1)
			break;


	}
	if (pos_l == -1)
		pos_l = 0;


	for (int i = 0; i <= pos_r - pos_l; ++i)
	{
		inst_name[i] = file[pos_l + i];
	}
	inst_name[strlen(inst_name) - 1] = 0;

	snprintf(network_name, sizeof(network_name), "%s", inst_name);
	int fract = (max_number_discovery_fraction * 100.00);
	int us_par = 0;
	if (this->uncertainty_parameter - ((int)uncertainty_parameter) < 0.001)
		us_par = uncertainty_parameter;
	else
		us_par = (uncertainty_parameter * 100.00);


	snprintf(inst_name, sizeof(inst_name), "%s_%s%d_%s%d_%s%d_%s%d_%s%d", inst_name, "K", K, "D", decision_depende_discovery, "md0", fract, "US", uncertainty_set_type, "USpar", us_par);

	printf("\n%s\n", network_name);
	printf("\n%s\n", inst_name);






	// Tmax
	fgets(line, sizeof(line), inp);
	ss = strtok(line, " ");
	//max_dur = round(atof(ss) * 100.00);
	double dur = round(atof(ss));
	max_dur = (dur * 100.00);
	//max_dur *= 100;


	/*fgets(line, sizeof(line), inp);
	ss = strtok(line, " ");
	ss = strtok(NULL, " ");
	g_num_nodes = atoi(ss);*/

	if(!read_matrix)
	{
	num_nodes = 0;
	start_node = 0;

	node_x_coord.push_back(0.0);
	node_y_coord.push_back(0.0);
	node_det_profit.push_back(0.0);
	fgets(line, sizeof(line), inp);
	sscanf(line, "%lf\t%lf\t%lf\n", &node_x_coord[start_node], &node_y_coord[start_node], &node_det_profit[start_node]);

	//disc_cost_n[start_node] = 0;
	disc_cost_n.push_back(0.0);


	//fgets(line, sizeof(line), inp);
	//sscanf(line, "%lf\t%lf\t%d\n", &node_x_coord[num_nodes - 1], &node_x_coord[num_nodes - 1], &node_x_coord[num_nodes - 1]);
	//g_sens_cost[MAX_NUM_NODES - 1] = 0;

	double temp_x_coord, temp_y_coord, temp_node_det_prof;

	fgets(line, sizeof(line), inp);

	node_x_coord.push_back(0.0);
	node_y_coord.push_back(0.0);
	node_det_profit.push_back(0.0);
	sscanf(line, "%lf\t%lf\t%lf\n", &temp_x_coord, &temp_y_coord, &temp_node_det_prof);




	this->tot_deterministic_profit = 0;
	// All the other nodes
	int node_i = 1;
	while (fgets(line, sizeof(line), inp) != NULL)
	{



		node_x_coord.push_back(0.0);
		node_y_coord.push_back(0.0);
		node_det_profit.push_back(0.0);
		sscanf(line, "%lf\t%lf\t%lf\n", &node_x_coord[node_i], &node_y_coord[node_i], &node_det_profit[node_i]);

		tot_deterministic_profit += node_det_profit[node_i];

		// Sensors cost
		disc_cost_n.push_back(0.0);
		//disc_cost_n[node_i] = 0;
		++node_i;
	}

	node_x_coord[node_i] = temp_x_coord;
	node_y_coord[node_i] = temp_y_coord;
	node_det_profit[node_i] = temp_node_det_prof;


	end_node = node_i;
	num_nodes = end_node + 1;
	num_prof_nodes = num_nodes - 2;

	first_pnode = start_node + 1;
	last_pnode = end_node - 1;

	node_is_profitable.assign(num_nodes, true);

	node_is_profitable[start_node] = false;
	node_is_profitable[end_node] = false;


	// num sesnors 
	if (decision_depende_discovery)
	{
		max_num_disc = ceil(num_prof_nodes * this->max_number_discovery_fraction);

	}
	else {
		max_num_disc = num_prof_nodes;
	}

	// Travel times
	this->t_ij.assignMatrix2D(this->num_nodes, this->num_nodes, 0);

	// Node existence
	this->edge_exists.assignMatrix2D(this->num_nodes, this->num_nodes, false);
	double dt_ij = 0;
	// Euclidean, profitable nodes and existing edges.
	for (int i = start_node; i < end_node; ++i)
		for (int j = i + 1; j <= end_node; ++j)
		{
			/*t_ij(i,j) = round(
				100.0 *
				sqrt(pow(node_x_coord[i] - node_x_coord[j], 2.0) + pow(node_y_coord[i] - node_y_coord[j], 2.0)));*/
				//((floorf(original * (double)INSTANCE_DISCR))) / (double)INSTANCE_DISCR;

				/*	dt_ij = std::sqrt(std::pow(node_x_coord[i] - node_x_coord[j], 2.0) + std::pow(node_y_coord[i] - node_y_coord[j], 2.0));
				dt_ij = floorf(100.0 * dt_ij)/100.0;
				t_ij(i,j) = dt_ij;*/

			t_ij(i, j) = //std::round(
				(int)(100.00 * std::sqrt(std::pow(node_x_coord[i] - node_x_coord[j], 2.0) + std::pow(node_y_coord[i] - node_y_coord[j], 2.0)))
				//)
				;
			//t_ij(i,j) = t_ij(i,j);

			this->edge_exists(i, j) = true;
		}
}
else
	{
	/*Num nodes*/
	fgets(line, sizeof(line), inp);
	ss = strtok(line, " ");
	num_nodes = (atoi(ss));

	this->start_node = 0;
	this->end_node = this->num_nodes - 1;

	this->first_pnode = start_node + 1;
	this->last_pnode = end_node - 1;

	num_prof_nodes = this->num_nodes - 2;

	node_is_profitable.assign(num_nodes, true);

	node_is_profitable[start_node] = false;
	node_is_profitable[end_node] = false;



	// num sesnors 
	if (decision_depende_discovery)
	{
		max_num_disc = ceil(num_prof_nodes * this->max_number_discovery_fraction);

	}
	else
	{
		max_num_disc = num_prof_nodes;
	}


	this->node_x_coord.assign(this->num_nodes, 0);
	this->node_y_coord.assign(this->num_nodes, 0);



	/*Fake profit*/
	node_det_profit.assign(this->num_nodes, 1);
	node_det_profit[start_node] = 0;
	node_det_profit[end_node] = 0;
	tot_deterministic_profit = 1;

	//tot_deterministic_profit = this->num_prof_nodes; /* 1 time num prof nodes */
	disc_cost_n.assign(this->num_nodes, 0);

	// Travel times
	this->t_ij.assignMatrix2D(this->num_nodes, this->num_nodes, 0);

	// Node existence
	this->edge_exists.assignMatrix2D(this->num_nodes, this->num_nodes, false);

	this->read_distance_matrix(inp);
	}


	return true;
}


/*Methods for Graph*/
inline bool ROPEUsolver::init_queue()
{
	sep_queue.num_el = 0;
	sep_queue.next = 0;
	return true;
}


inline bool ROPEUsolver::queue_is_empty()
{
	if (sep_queue.num_el == 0
		||
		sep_queue.num_el == sep_queue.next)
		return true;
	else
		return false;
}

inline bool ROPEUsolver::push_back(int el)
{
	sep_queue.elem[sep_queue.num_el] = el;
	++sep_queue.num_el;
	return true;
}

inline int ROPEUsolver::pop_front()
{
	int elem = sep_queue.elem[sep_queue.next];
	++sep_queue.next;
	return elem;
}

inline int ROPEUsolver::build_sol_graph(const double* curr_sol, int k, double var_treshold = 0)
{
	sep_un_graph.num_nodes = 0;
	sep_un_graph.num_edges = 0;

	sep_un_graph.num_nodes = 1;
	sep_un_graph.nodes_map[sep_un_graph.num_nodes - 1] = start_node;
		//g_start_node;

	for (int i = start_node + 1; i < end_node; ++i)
		if (curr_sol[v_y_kn(k,i)] > mEPS)
		{

			sep_un_graph.nodes_map[sep_un_graph.num_nodes] = i;
			++sep_un_graph.num_nodes;
		}

	sep_un_graph.nodes_map[sep_un_graph.num_nodes] = end_node;
	++sep_un_graph.num_nodes;

	for (int i = 0; i < sep_un_graph.num_nodes - 1; ++i)
		for (int j = i + 1; j < sep_un_graph.num_nodes; ++j)
			if (curr_sol[v_x_kij(k,sep_un_graph.nodes_map[i],sep_un_graph.nodes_map[j])] > mEPS + var_treshold)
			//if (curr_sol[g_v_xk_kij[k][sep_un_graph.nodes_map[i]][sep_un_graph.nodes_map[j]]] > mEPS)
			{

				// Add edge.
				sep_un_graph.edge[sep_un_graph.num_edges].first = i;
				sep_un_graph.edge[sep_un_graph.num_edges].second = j;
				++sep_un_graph.num_edges;
			}
	return 0;
}

inline int ROPEUsolver::build_tilde_sol_graph(const double* curr_sol, int k, double val_threshold = 0)
{
	sep_un_graph.num_nodes = 0;
	sep_un_graph.num_edges = 0;

	sep_un_graph.num_nodes = 1;
	sep_un_graph.nodes_map[sep_un_graph.num_nodes - 1] = start_node;

	for (int i = start_node + 1; i < end_node; ++i)
		if (curr_sol[v_t_y_kn(k, i)] > mEPS )
		{

			sep_un_graph.nodes_map[sep_un_graph.num_nodes] = i;
			++sep_un_graph.num_nodes;
		}

	sep_un_graph.nodes_map[sep_un_graph.num_nodes] = end_node;
	++sep_un_graph.num_nodes;

	for (int i = 0; i < sep_un_graph.num_nodes - 1; ++i)
		for (int j = i + 1; j < sep_un_graph.num_nodes; ++j)
			if (curr_sol[v_t_x_kij(k, sep_un_graph.nodes_map[i], sep_un_graph.nodes_map[j])] > mEPS + val_threshold)
			{
				sep_un_graph.edge[sep_un_graph.num_edges].first = i;
				sep_un_graph.edge[sep_un_graph.num_edges].second = j;
				++sep_un_graph.num_edges;
			}
	return 0;
}


unsigned int ROPEUsolver::connected_components()
{
	unsigned int i;
	unsigned int component = 0;

	for (int i = 0; i < sep_un_graph.num_nodes; ++i)
		sep_un_graph.conn_comp[i] = -1;



	for (i = 0; i < sep_un_graph.num_nodes; i++) {
		if (sep_un_graph.conn_comp[i] == -1) {
			connected_components_recursive(i, component);
			component++;
		}
	}
	sep_un_graph.num_cc = component;
	return component;
}

inline void ROPEUsolver::connected_components_recursive(unsigned int vertex, unsigned int component)
{
	unsigned int i;
	/* Put this vertex in the current component */
	sep_un_graph.conn_comp[vertex] = component;
	//COMPONENTS[vertex] = component;
	for (i = 0; i < sep_un_graph.num_edges; i++) {
		if (sep_un_graph.edge[i].first == vertex || sep_un_graph.edge[i].second == vertex) {
			/* Adjacent */
			const unsigned int neighbour = sep_un_graph.edge[i].first == vertex ?
				sep_un_graph.edge[i].second : sep_un_graph.edge[i].first;
			if (sep_un_graph.conn_comp[neighbour] == -1) {
				/* Not yet visited */
				connected_components_recursive(neighbour, component);
			}
		}
	}
}


inline bool ROPEUsolver::build_sol_directed_graph(const double* curr_sol, int k)
{


	sep_di_graph.num_ver = 0;



	sep_di_graph.num_ver = 1;
	sep_di_graph.ver_label[sep_di_graph.num_ver - 1] = start_node;




	for (int i = start_node + 1; i < end_node; ++i)
		if (curr_sol[v_y_kn(k,i)] > mEPS)
		{
			sep_di_graph.ver_label[sep_di_graph.num_ver] = i;
			sep_di_graph.adj_size[sep_di_graph.num_ver] = 0;
			++sep_di_graph.num_ver;

		}

	// Init number of archs to 0.
	for (int i = 0; i < sep_di_graph.num_ver; ++i)
		sep_di_graph.adj_size[i] = 0;


	// First arcs not related with depot
	for (int i = 1; i < sep_di_graph.num_ver - 1; ++i)
		for (int j = i + 1; j < sep_di_graph.num_ver; ++j)
			if (curr_sol[v_x_kij(k,sep_di_graph.ver_label[i],sep_di_graph.ver_label[j])] > mEPS)
		//	if (curr_sol[g_cg_price.v_x_ij[sep_di_graph.ver_label[i]][sep_di_graph.ver_label[j]]] > mmEPS)
			{
			

				double dq =
					(curr_sol[
						v_x_kij(k,
						sep_di_graph.ver_label[i]
						,sep_di_graph.ver_label[j])
					]) 
					* scaleMXflow;
					//(curr_sol[g_cg_price.v_x_ij[sep_di_graph.ver_label[i]][sep_di_graph.ver_label[j]]]) * scaleMXflow;


				int q = (int)dq;
				add_archs(i, j, q);
			}

	// Now the depot (the capacity is given by the sum)
	for (int i = 1; i < sep_di_graph.num_ver; ++i)
		if (curr_sol[v_x_kij(k,start_node,sep_di_graph.ver_label[i])] > mEPS
			||
			curr_sol[v_x_kij(k,sep_di_graph.ver_label[i],end_node)] > mEPS)
	/*	if (curr_sol[g_cg_price.v_x_ij[g_start_node][sep_di_graph.ver_label[i]]] > mEPS
			||
			curr_sol[g_cg_price.v_x_ij[sep_di_graph.ver_label[i]][g_end_node]] > mEPS)*/
		{
			double dq =
				(curr_sol[v_x_kij(k,start_node,sep_di_graph.ver_label[i])] +
					curr_sol[v_x_kij(k,sep_di_graph.ver_label[i],end_node)]) * scaleMXflow;
				/*(curr_sol[v_x_kij[g_start_node][sep_di_graph.ver_label[i]]] +
					curr_sol[g_cg_price.v_x_ij[sep_di_graph.ver_label[i]][g_end_node]]) * scaleMXflow;*/
			int q = (int)dq;

			add_archs(0, i, q);
		}


	return 0;
}

inline bool ROPEUsolver::build_tilde_sol_directed_graph(const double* curr_sol, int k)
{
	sep_di_graph.num_ver = 0;



	sep_di_graph.num_ver = 1;
	sep_di_graph.ver_label[sep_di_graph.num_ver - 1] = start_node;




	for (int i = start_node + 1; i < end_node; ++i)
		if (curr_sol[v_t_y_kn(k, i)] > mEPS)
		{
			sep_di_graph.ver_label[sep_di_graph.num_ver] = i;
			sep_di_graph.adj_size[sep_di_graph.num_ver] = 0;
			++sep_di_graph.num_ver;

		}

	// Init number of archs to 0.
	for (int i = 0; i < sep_di_graph.num_ver; ++i)
		sep_di_graph.adj_size[i] = 0;


	// First arcs not related with depot
	for (int i = 1; i < sep_di_graph.num_ver - 1; ++i)
		for (int j = i + 1; j < sep_di_graph.num_ver; ++j)
			if (curr_sol[v_t_x_kij(k, sep_di_graph.ver_label[i], sep_di_graph.ver_label[j])] > mEPS)
			{


				double dq =
					(curr_sol[
						v_t_x_kij(k,
							sep_di_graph.ver_label[i]
							, sep_di_graph.ver_label[j])
					])
					* scaleMXflow;


						int q = (int)dq;
						add_archs(i, j, q);
			}

	// Now the depot (the capacity is given by the sum)
	for (int i = 1; i < sep_di_graph.num_ver; ++i)
		if (curr_sol[v_t_x_kij(k, start_node, sep_di_graph.ver_label[i])] > mEPS
			||
			curr_sol[v_t_x_kij(k, sep_di_graph.ver_label[i], end_node)] > mEPS)
	
		{
			double dq =
				(curr_sol[v_t_x_kij(k, start_node, sep_di_graph.ver_label[i])] +
					curr_sol[v_t_x_kij(k, sep_di_graph.ver_label[i], end_node)]) * scaleMXflow;
			
			int q = (int)dq;

			add_archs(0, i, q);
		}


	return 0;
}

inline void ROPEUsolver::add_archs(int a1, int a2, int q)
{
	int pos_a1 = sep_di_graph.adj_size[a1];
	int pos_a2 = sep_di_graph.adj_size[a2];

	sep_di_graph.adj_archs(a1,pos_a1).a1 = a1;
	sep_di_graph.adj_archs(a1,pos_a1).a2 = a2;

	sep_di_graph.adj_archs(a1,pos_a1).q = q;
	sep_di_graph.adj_archs(a1,pos_a1).rev = pos_a2;
	sep_di_graph.adj_archs(a1,pos_a1).flow = 0; // init flow to 0

	// Now the reverse
	sep_di_graph.adj_archs(a2,pos_a2).a1 = a2;
	sep_di_graph.adj_archs(a2,pos_a2).a2 = a1;

	sep_di_graph.adj_archs(a2,pos_a2).q = q;
	sep_di_graph.adj_archs(a2,pos_a2).rev = pos_a1;
	sep_di_graph.adj_archs(a2,pos_a2).flow = 0; // init flow to 0


	++sep_di_graph.adj_size[a1];
	++sep_di_graph.adj_size[a2];
}

inline bool ROPEUsolver::bfs_di_graph(int s, int t)
{
	int level;
	int flow;

	for (int i = 0; i < sep_di_graph.num_ver; ++i)
	{
		sep_di_graph.level[i] = -1;
	}
	sep_di_graph.level[s] = 0;

	// init queue
	init_queue();
	push_back(s);

	while (!queue_is_empty())
	{
		int u = pop_front();

		for (int a = 0; a < sep_di_graph.adj_size[u]; ++a)
		{
			level = sep_di_graph.level[sep_di_graph.adj_archs(u,a).a2];
			flow = sep_di_graph.adj_archs(u,a).flow;

			if (level < 0 && flow < sep_di_graph.adj_archs(u,a).q)
			{
				sep_di_graph.level[sep_di_graph.adj_archs(u,a).a2] = sep_di_graph.level[u] + 1;

				push_back(sep_di_graph.adj_archs(u,a).a2);

			}
		}

	}

	return sep_di_graph.level[t] < 0 ? false : true;
}

inline int ROPEUsolver::send_flow_di_graph(int s, int flow, int t, short* start)
{
	// Sink reached 
	if (s == t)
		return flow;

	// Traverse adiacent edges
	for (; start[s] < sep_di_graph.adj_size[s]; ++start[s])
	{
		arch_t* e = &(sep_di_graph.adj_archs(s,start[s]));

		if (sep_di_graph.level[e->a2] == sep_di_graph.level[s] + 1
			&&
			e->flow < e->q
			)
		{
			int curr_flow = std::min(flow, e->q - e->flow);
			int temp_flow = send_flow_di_graph(e->a2, curr_flow, t, start);

			if (temp_flow > 0)
			{
				// flow to current edge
				e->flow += temp_flow;

				// subtract flow from reverse edge of curr edge
				int rev_arc = e->rev;
				sep_di_graph.adj_archs(e->a2,e->rev).flow -= temp_flow;
				return temp_flow;

			}

		}

	}


	return 0;
}

inline int ROPEUsolver::dinic_max_flow(int s, int t)
{
	// Corner case 
	if (s == t)
		return -1;

	// store how many edges are visited 
		// from V { 0 to V } 
	short start[MAX_NUM_VERTICES];
	for (int i = 0; i < sep_di_graph.num_ver; ++i)
		start[i] = 0;


	int total = 0;  // Initialize result 

	// Augment the flow while there is path 
	// from source to sink 
	while (bfs_di_graph(s, t) == true)
	{

		for (int i = 0; i < sep_di_graph.num_ver; ++i)
			start[i] = 0;

		int flow;
		do {
			flow = send_flow_di_graph(s, INT_MAX, t, start);
			total += flow;
		} while (flow > 0);


	}

	for (int i = 0; i < end_node; ++i)
	{
		sep_mcut_set_1[i] = 0;
		//sep_mcut_sets[mCUT_S1][i] = 0;
	}

	int qi;
	bool start_is_in_q = false;
	for (int q = 0; q < sep_queue.num_el; ++q)
	{
		if (sep_queue.elem[q] == 0)
		{
			start_is_in_q = true;
			break;
		}
	}

	if (start_is_in_q)
	{
		for (int i = 0; i < sep_di_graph.num_ver; ++i)
		{
			sep_mcut_set_1[sep_di_graph.ver_label[i]] = 1;
			//sep_mcut_sets[mCUT_S1][sep_di_graph.ver_label[i]] = 1;
		}

		for (int q = 0; q < sep_queue.num_el; ++q)
		{
			qi = sep_di_graph.ver_label[sep_queue.elem[q]];

			sep_mcut_set_1[qi] = 0;
			//sep_mcut_sets[mCUT_S1][qi] = 0;
		}
	}
	else
	{
		for (int q = 0; q < sep_queue.num_el; ++q)
		{
			qi = sep_di_graph.ver_label[sep_queue.elem[q]];

			sep_mcut_set_1[qi] = 1;
			//sep_mcut_sets[mCUT_S1][qi] = 1;
		}
	}


	// return maximum flow 
	return total;
}