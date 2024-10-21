#pragma once
#include "ProblemDef.h"
#include "Utilities.h"
#include "ilcplex/cplex.h"
#include "ConfigFile.h"
#include <bitset>


#define EPS 0.1
#define mEPS 0.0001 /*Small eps*/
#define cEPS 0.005 /* cut Small eps*/


/*Treshold for cut EXACT algo*/
#define EXA_EPS 1e-4
/*Teshodol stopping creteria subproblem (difference between Master and Sub in CG)*/
//#define EXA_EPS_conv 1e-5

/*Master Benders*/
#define EXA_Benders_MIN_COEFF_LP_IMPR 0.01
#define EXA_Benders_MAX_NUM_ITER_NODE 1
#define EXA_Benders_MAX_NUM_ITER_NODE_DEPTH_1 1


#define BIG_M_SIGMA 100
#define MAX_NUM_SCENARIOS 10000


/*B&C*/
#define BandC_NUM_CUT_PER_NODE 500 
#define BandC_MIN_COEFF_LP_IMPR 0.01
#define BandC_MAX_NUM_ITER_NODE 5
#define BandC_SET_BRANCHING_PRIORITY true


#define MAX_NUM_CUTS_ELEM 1000 /*this is the constant to init bitset*/

class ExaSolverGI
{

//	char inst_name[1024];
//	char network_name[1024];
//	FILE* out_file;
	//char results_folder[1024];

	ProblemDef* problem;

//	bool decision_depende_discovery; /*true for decision dependent discovery version. False for w = e.*/
//	double max_number_discovery_fraction; /* \in [0,1]  ceil(max_number_discovery_fraction * num_obeservable_elements) = max_num_disc*/
	/********************/


	bool init_cplex(CPXENVptr* cpx_env,
		CPXLPptr* cpx_lp);
	int  checkCPXstatus(int status, CPXENVptr* cpx_env,
		CPXLPptr* cpx_lp);
	bool free_cplex(CPXENVptr* cpx_env_op,
		CPXLPptr* cpx_lp_op);



	/*Discovery*/
	std::vector<bool> w_val; /*True if component i is observed 0 otherwise. */


	Configuration::ConfigFile* cfg;






	/*To write constraints and store solutions*/
	std::vector<int> matind;
	std::vector<double> matval;
	std::vector<double> rop_curr_sol;
	std::vector<double> mlp_curr_sol;
	std::vector<double> bm_curr_sol;
	std::vector<double> vector_loc;
	Matrix2D<double> matr_loc;

	double mlp_obj_val;
	double rop_obj_val;
	bool rop_is_feasible; /*This check if the nominal problem Y is feasible.*/
	double bm_obj_val;

	double bm_cut_lb; // The full discovery bounds used for the cuts. Updated every time.



	unsigned int card_Y;
	Matrix2D<unsigned int> Y; // Each row is a solution of the OP.
	std::vector<bool> y_in_mlp; // With respect of the columns of Y, y_in_mlp[y] tells if y has been added to the mlp.


	/*Cplex structures for OP problem.*/
	CPXENVptr cpx_env_rop;
	CPXLPptr  cpx_lp_rop;

	/*Cplex structures for LP */
	CPXENVptr cpx_env_mlp;
	CPXLPptr  cpx_lp_mlp;


	/*OP problem */
	std::vector<unsigned short> rop_v_y_n; /*node*/
	std::vector<unsigned short> v_beta; /*dual Uncertainty set constraints*/
	std::vector<unsigned short> v_gamma; /*dual non-anticipativity constraints*/


	/*LP master*/
	unsigned int mlp_num_y;
	short mlp_v_tau;
	std::vector<unsigned int> mlp_v_pr_xi; /*Prime xi: it does not increase.*/
	Matrix2D<unsigned int> mlp_v_xi_y; /* (y, n): one set of variable for each solution of the OP*/
	std::vector<unsigned int>mlp_cnstr_tau_y; /*Index for policy contraints*/


	/*Master benders*/
	long int bm_num_iter_curr_node;
	long int bm_num_cut_curr_iter;
	long int bm_prev_node;
	long int bm_curr_node;
	double bm_lp_prev_iter;
	void bm_reset_branc_and_cut_indicators();






	/*Sub-sub problem*/
	bool rop_build_model();
	bool rop_set_cplex();
	bool rop_add_vars();
	bool rop_add_cnstr();
	bool rop_update_scenario_parameter(); /*update xi_prime*/
	bool rop_solve_model();

	/*LP problem*/
	bool mlp_build_model();
	bool mlp_set_cplex();
	bool mlp_add_vars();
	bool mlp_add_cnstr();
	bool mlp_update_model(); /*add columns and rows*/
	bool add_columns();
	bool add_rows();
	bool mlp_solve_model();

	/*Methods to clean the Y. Deprecated for now*/
	bool check_violation_in_Y();
	bool mlp_add_column_and_rows(int y); // Add a specified y solution from set Y.
	bool add_columns_y(int y); //Add columns for a  y solution from set Y.
	bool add_rows_y(int y); //Add rows for a  y solution from set Y.

	/*Update subproblem models*/

	//bool subproblem_solver_update_discovery_parameters(const std::vector<bool>* w_val);
	bool mlp_update_discovery_parameters();
	bool rop_update_discovery_parameters();

	/*Bender Master Problem*/
	CPXENVptr cpx_env_bm;
	CPXLPptr  cpx_lp_bm;
	unsigned int bm_v_phi;
	std::vector<unsigned int> bm_v_w_i;
	bool bm_build_model();
	bool bm_set_cplex();
	bool bm_add_vars();
	bool bm_add_cnstr();
	bool set_branching_score();
	bool bm_solve_model();
	bool bm_add_starting_cuts();
	static int CPXPUBLIC bm_general_callback(CPXCALLBACKCONTEXTptr context, CPXLONG contextid, void* userdata);
	bool bm_separate_logic_benders_cut(CPXCALLBACKCONTEXTptr context, CPXLONG contextid, double var_treshold, bool* cut_added);

	/*Storing cuts*/
	uint64_t bm_num_stored_logic_cuts;
	//std::vector<uint64_t> bm_logic_cuts_pool_set; /*Query adding a bit set*/
	std::vector<double> bm_logic_cuts_pool_phi_val;
	/**/
	std::vector<bitset<MAX_NUM_CUTS_ELEM>> bm_logic_cuts_pool_set;



	void bm_add_elem_to_logic_cut(bitset<MAX_NUM_CUTS_ELEM>* set, unsigned int elem);
	bool bm_is_logic_cut_already_added(bitset<MAX_NUM_CUTS_ELEM>* cut, int* index);
	void bm_add_logic_cut_to_pool(bitset<MAX_NUM_CUTS_ELEM>* cut, const double* phi);

	std::vector<int> bm_branch_priority_score;


	/*Data structure for the additional variables CGG to have a starting Lower Bound.*/
	/*Second stage problem part*/
	Matrix2D<unsigned int>bm_v_y_xib_i; /*Second stage*/
	Matrix2D<unsigned int> bm_v_lamb_xib_i; /* Dual of matrix A. One vector lambda for each xi_bar and for each constraint in A */
	Matrix2D<unsigned int> bm_v_sig_xib_i; /*Dual of xi_bar = xi*/

	/*Ise information cuts or not*/
	bool bm_use_info_cuts;
	/*Data structure to store the scenarios*/
	bool use_xi_bar_scenarios;
	unsigned int num_ccg_xib;
	Matrix2D<double> ccg_set_xib;
	std::vector<unsigned int> cnstr_ccg_xi; /*index of the constraints. To check which xi is tight.*/
	bool store_xbar();
	bool bm_add_ccg_lb();/*This add variables and constraints for each stored xi*/
	bool bm_add_ccg_columns();
	bool bm_add_ccg_constraints();


	/*Stat*/
	double bm_best_lb;
	double bm_best_ub;
	double bm_branch_and_cut_time;
	double bm_starting_cuts_time;
	double cplex_gap;
	unsigned int bandc_nodes;
	unsigned int num_logic_benders_cuts;
	double root_node_lb;

	/*Stores tight policies*/
	bool store_set_of_tight_y();
	void store_tight_y(int y);
	unsigned int num_tight_y;
	std::vector<unsigned int> tight_y;
	Matrix2D<double> tight_xi_y;
	std::vector<double> tight_xibar;


public:
	//bool read_instance(char* file);

	bool setConfigFile(const char* conf_file);

	bool init_data_structures();

	bool clean_data_structure();

	void set_problem(ProblemDef* problem); 

	bool load_and_init(std::vector<bool> w_val);

	bool run_sub_problem_solver(double* obj_val);

	bool subproblem_solver_update_discovery_parameters(const std::vector<bool>* w_val);

	bool run_exact_method(double* objective_value);

	//bool print_solution_file();

	//bool print_solution_latex();

	double getGap()
	{
		return this->cplex_gap;
	}
	double getBestLB()
	{
		return this->bm_best_lb;
	}
	double getBestUB()
	{
		return this->bm_best_ub;
	}
	double getBandCTime()
	{
		return this->bm_branch_and_cut_time;
	}
	double getStartTime()
	{
		return this->bm_starting_cuts_time;
	}

	unsigned int getNum_ccg_xib()
	{
		return this->num_ccg_xib;
	}

	bool getUseXiBarScenarios()
	{
		return this->use_xi_bar_scenarios;
	}

	bool getRop_is_feasible()
	{
		return this->rop_is_feasible;
	}

	std::vector<double> getBM_curr_sol() {
		return this->bm_curr_sol;
	}

	std::vector<unsigned int> getIndW() {
		return this->bm_v_w_i;
	}
};

