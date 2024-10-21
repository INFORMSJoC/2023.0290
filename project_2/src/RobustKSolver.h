#pragma once
#include "K_Adapt_Problem.h"
#include "ConfigFile.h"
#include "ilcplex/cplexx.h"
#include "ilcplex/cplex.h"

#define MAX_NUM_VARS 1000000
//#define USE_STANDARD_FORMULATION true // Set to true to use the old K-adaptability formulation, false to use what we propose in the paper.
#define TIME_LIMIT 7200

typedef struct dual_vars
{
	 //This is used to check for which constraints we force infeas. In the approximation model, is used to store vars.
	Matrix2D<unsigned short> z_k;
	std::vector<unsigned short> v_alpha; // Alpha variables
	std::vector<unsigned short> v_beta; 
	Matrix2D<unsigned short> v_beta_k;
	Matrix2D<unsigned short> v_gamma_k; //

	// Bilinear variables
	Matrix2D<unsigned short> v_t_gamma_k; // gamma times w.
	Matrix2D<unsigned short> v_t_y_k; // alpha times y. Components of y_k are moltiplied only for alpha k. You can then query as a y variables. 
	Matrix2D<unsigned short> v_t_w_k; // alpha times w. Since we have one alpha for each k, we query with [k][w_ind]
	Matrix2D<unsigned short> v_t_x_k; // alpha times x. Since we have one alpha for each k, we query with [k][x_ind]

	// Dual for continuos approximation of CU. approximation
	std::vector<unsigned short> v_rho; // for each k.
	Matrix3D<unsigned short> v_ov_mu_kl; // for each kl a vector.
	Matrix3D<unsigned short> v_un_mu_kl; // for each kl a vector.
	Matrix3D<unsigned short> v_ov_pi_kl; // for each kl a vector.
	Matrix3D<unsigned short> v_un_pi_kl; // for each kl a vector.

}dual_vars_t;





class RobustKSolver
{
private:

	K_Adapt_Problem* problem;

	/*Cplex structures for the master problem.*/
	CPXENVptr cpx_env_master;
	CPXLPptr  cpx_lp_master;

	/*Cplex structure for the speration problem in CU*/
	CPXENVptr cpx_env_slave;
	CPXLPptr  cpx_lp_slave;

	/*Cpex structures for approximation problem. Note that is equivalent to OU for L = 0*/
	CPXENVptr cpx_env_milp;
	CPXLPptr  cpx_lp_milp;


	Configuration::ConfigFile* cfg;


	// Matind and Matval, and Solution
	std::vector<int> matind;
	std::vector<double> matval;
	std::vector<double> master_curr_sol;
	double master_curr_blb;
	double master_curr_bub;
	std::vector<double> slave_curr_sol;



	// Phi var
	unsigned short v_phi;
	// Now we have continer for deterministic variables 
	std::vector<unsigned short> v_x; 
	std::vector<unsigned short> v_w;
	Matrix2D<unsigned short> v_y_k;
	// Container for dual variables
	std::vector<dual_vars_t> v_dual_set;
	/*Set of dual variables for approximated model*/
	dual_vars_t appr_dual_vars;


	//General methods
	bool init_data_structures();
	int  checkCPXstatus(CPXENVptr env, int status);
	bool free_cplex(CPXENVptr* env, CPXLPptr* lp);
	bool free_cplex_problem(CPXENVptr* env,CPXLPptr* lp);
	


	/*Master Problem*/
	bool init_cplex_data_str();
	bool set_cplex_parameters();
	bool build_master_model();
	bool add_primal_variables();
	bool add_dualization_variables(unsigned int z);
	bool add_bilinear_variables(unsigned int z);
	bool add_deterministic_constraints();// Xx <= b_x, Ww <=b_w, Yy <=b_y, Linking <=g
	bool add_primal_valid_inequalities(); // Optimistic cuts
	bool add_dualization_constraints(unsigned int z);
	bool add_linearization_constraints(unsigned int z);
	bool add_dual_vars_valid_inequalities(unsigned int z);
	bool add_RLT_AlphaXWY_valid_inequalities(unsigned int z);
	bool add_RLT_AlphaExwy_valid_inequalities(unsigned int z);
	bool add_RLT_Primal_valid_inequalities(unsigned int z);
	bool add_prev_iter_Phi_LB(unsigned int z);
	bool master_update_problem(unsigned int z);
	bool master_solve_problem(double* obj_val, double* unc_cost); // Solve the master and save the solution

	
	/*Separation Problem*/
	unsigned short vs_tau;
	std::vector<unsigned short> vs_ov_xi, vs_ov_zeta;
	Matrix2D<unsigned short> vs_xi_k, vs_zeta_k, vs_z_k;
	Matrix3D<unsigned short> vs_chi_kl;

	// The solution is stored in curr_sol, then you can build the model.
	bool slave_init_cplex_data_str();
	bool slave_set_cplex_parameters();
	bool slave_build_problem();
	bool slave_add_variables();
	bool slave_add_constraints();
	bool slave_add_linearization_constraints();
	bool slave_solve_problem(double* obj_val);
	bool slave_store_z_values();


	/*Approximation: deprecated*/
	bool appr_init_cplex_data_str();
	bool appr_set_cplex_parameters();
	bool appr_build_model();
	bool appr_solve_problem(double* obj_val);
	bool appr_add_primal_variables();
	bool appr_add_dualization_variables();
	bool appr_add_bilinear_variables();
	bool appr_add_deterministic_constraints();// Xx <= b_x, Ww <=b_w, Yy <=b_y, Linking <=g
	bool appr_add_primal_valid_inequalities();
	bool appr_add_dualization_constraints(); // Later
	bool appr_add_linearization_constraints();
	bool appr_add_dual_vars_valid_inequalities();
	bool appr_add_RLT_AlphaXWY_valid_inequalities();
	bool appr_add_RLT_AlphaExwy_valid_inequalities();
	bool appr_add_RLT_Primal_valid_inequalities();



public:
	//Set configuration file
	bool set_configuration_file(char* config_file);
	void del_config_file();
	/*C&C G*/
	bool column_and_constraints_generation(double* obj_val);
	bool solve_problem();
	bool appr_solver(double* obj_val);
	void set_problem_info(K_Adapt_Problem* problem);

	double cplex_gap;
	double bub;
	double blb;
	double comp_time;

	std::vector<unsigned short> get_v_w_ind()
	{
		return this->v_w;
	}

	Matrix2D<unsigned short> get_v_y_k_ind()
	{
		return this->v_y_k; 
	}

	std::vector<double> getCurrSol() {
		return this->master_curr_sol;
	}

};

