#pragma once
#include "Utilities.h"
#include "ilcplex/cplex.h"
#include "ConfigFile.h"
#include <filesystem>

#define EPS 0.1
#define mEPS 0.0001 /*Small eps*/
#define cEPS 0.005 /* cut Small eps*/


/*Treshold for cut EXACT algo*/
#define EXA_EPS 1e-5


/*Master Benders*/
#define EXA_Benders_MIN_COEFF_LP_IMPR 0.01
#define EXA_Benders_MAX_NUM_ITER_NODE 1
#define EXA_Benders_MAX_NUM_ITER_NODE_DEPTH_1 100

#define BIG_M_SIGMA 100
#define MAX_NUM_SCENARIOS 10000 


/*B&C*/
#define BandC_NUM_CUT_PER_NODE 200 
#define BandC_MIN_COEFF_LP_IMPR 0.01
#define BandC_MAX_NUM_ITER_NODE 5
#define BandC_SET_BRANCHING_PRIORITY true

#define TRESHOLD_MAX 0.5 /*Subtour sepation*/
#define TRESHOLD_STEP 0.25 /*Subtour sepation*/






class EXAsolver
{

private:

	/*
Instance info
*/
	char inst_name[1024];
	char network_name[1024];
	FILE* out_file = NULL;
	char results_folder[1024];

	/*From configuration*/
	unsigned short uncertainty_set_type; /*Uncertainty set types: 1 or 2 for now*/
	double uncertainty_parameter; /*Uncertainty set parameter. For unceertainty set type 1 is the budget. For US4 is the theta, for other US not used*/

	bool decision_depende_discovery; /*true for decision dependent discovery version. False for w = e.*/
	double max_number_discovery_fraction; /* \in [0,1]  ceil(max_number_discovery_fraction * num_obeservable_elements) = max_num_disc*/
	/********************/


	bool init_cplex(CPXENVptr* cpx_env,
		CPXLPptr* cpx_lp);
	int  checkCPXstatus(int status, CPXENVptr* cpx_env,
	CPXLPptr*  cpx_lp);
	bool free_cplex(CPXENVptr* cpx_env_op,
		CPXLPptr* cpx_lp_op);

	bool init_out_file();


	

	/*Master Benders Problem*/


	/*OP data*/

	/*Discovery*/
	std::vector<bool> w_val; /*True if component i is observed 0 otherwise. */


	/*Nodes*/
	unsigned short num_nodes;
	unsigned short num_prof_nodes; // Number of profitable node.
	unsigned short start_node;
	unsigned short end_node;
	unsigned short first_pnode;
	unsigned short last_pnode;

	/*Alrijne collecting time*/
	double alrijne_node_collecting_time;
	/*Rhs constraints*/
	unsigned short max_num_disc;
	long int max_dur; // Max duration 

	/*Costs and profits*/
	std::vector<double>  disc_cost_n;

	std::vector<double> node_det_profit; //  Deterministic instance profit.
	double tot_deterministic_profit;
	std::vector<double> node_x_coord;
	std::vector<double> node_y_coord;

	Matrix2D<double> t_ij; // Travel time.
	Matrix2D<unsigned short> edge_exists; /*Used with 0 or 1.*/
	std::vector<bool> node_is_profitable;


	Configuration::ConfigFile* cfg;


	/*To write constraints*/
	std::vector<int> matind;
	std::vector<double> matval;
	std::vector<double> rop_curr_sol;
	std::vector<double> mlp_curr_sol;
	std::vector<double> bm_curr_sol;
	double mlp_obj_val;
	double rop_obj_val;
	double bm_obj_val;

	double bm_cut_lb; // Initialized as the full discovery bound. Updated at each call of the callback.


	unsigned short num_us_constr;
	Matrix2D<double> A; // Uncertainty set for xi vars.
	std::vector<double> b; // Rhs for uncertainty set.
	/*Obejective function object*/
	std::vector<double>ppxi; // Function of \xi in the OP.
	std::vector<double> pp; // part not function of P.


	unsigned int card_Y;
	Matrix2D<unsigned int> Y; // Each row is a solution of the OP.
	std::vector<bool> y_in_mlp; // With respect of the columns of Y, y_in_mlp[y] tells if y has been added to the mlp.


	/*Cplex structures for OP problem.*/
	CPXENVptr cpx_env_op;
	CPXLPptr  cpx_lp_op;

	/*Cplex structures for LP */
	CPXENVptr cpx_env_mlp;
	CPXLPptr  cpx_lp_mlp;


	/*OP problem */
	std::vector<unsigned short> rop_v_y_n; /*node*/
	Matrix2D<unsigned short> v_x_ij; /*node i, node j*/
	std::vector<unsigned short> v_beta; /*dual Uncertainty set constraints*/
	std::vector<unsigned short> v_gamma; /*dual non-anticipativity constraints*/





	/*LP master*/
	unsigned int mlp_num_y;
	short mlp_v_tau;
	std::vector<unsigned int> mlp_v_pr_xi; /*Prime xi: it does not increase.*/
	Matrix2D<unsigned int> mlp_v_xi_y; /* (y, n): one set of variable for each solution of the OP*/
	std::vector<unsigned int>mlp_cnstr_tau_y ; /*Index for policy contraints*/



	/*Callback */
	
	/*
	Callbackand cutting plane methods
		*/

	static int CPXPUBLIC rop_general_callback(CPXCALLBACKCONTEXTptr context, CPXLONG contextid, void* userdata);


	bool separate_heuristic_subtours(CPXCALLBACKCONTEXTptr context, CPXLONG contextid, double var_treshold, bool* cut_added);
	bool separate_exact_subtuors(CPXCALLBACKCONTEXTptr context, CPXLONG contextid, bool* cut_added);
	bool separate_logic_cuts(CPXCALLBACKCONTEXTptr context, CPXLONG contextid, bool* cut_added);

	long int num_iter_curr_node;
	long int num_cut_curr_iter;
	long int prev_node;
	long int curr_node;
	double lp_prev_iter;
	void reset_branc_and_cut_indicators();

	/*Master benders*/
	long int bm_num_iter_curr_node;
	long int bm_num_cut_curr_iter;
	long int bm_prev_node;
	long int bm_curr_node;
	double bm_lp_prev_iter;
	void bm_reset_branc_and_cut_indicators();


	/*Separations methods*/
	un_graph_t sep_un_graph;
	di_graph_t sep_di_graph;



	queue_t sep_queue;
	bool init_queue();
	bool queue_is_empty();
	bool push_back(int el);
	int pop_front();


	// Connected components
	int rop_build_sol_graph(const double* curr_sol,  double val_treshold);
	unsigned int connected_components();
	void connected_components_recursive(unsigned int vertex,
		unsigned int component);

	int bm_build_sol_graph(const double* curr_sol, double val_treshold, int xi_bar);

	// Max-Flow with Dinic's Alg.
	bool build_sol_directed_graph(const double* curr_sol);// create the directed graph from curr sol.
	bool bm_build_sol_directed_graph(const double* curr_sol, int xi);// create the directed graph from curr sol.
	void add_archs(int a1, int a2, int q);
	bool bfs_di_graph(int s, int t); // breadth first search
	int send_flow_di_graph(int s, int flow, int t, short* start);
	int dinic_max_flow(int s, int t);
	// min cut: rectrive the min cut from the max_flow.
	/*Min cuts set*/

	std::vector<bool> sep_mcut_set_1;
	std::vector<bool> sep_mcut_set_2;
	std::vector<unsigned short> node_subtour_subset;



	/*ROP problem*/
	bool rop_build_model();
	bool rop_set_cplex();
	bool rop_add_vars();
	bool rop_add_cnstr();
	bool rop_update_scenario_parameter_model(); /*update xi_prime*/
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

	uint64_t bm_num_stored_logic_cuts;
	std::vector<uint64_t> bm_logic_cuts_pool_set; /*Query adding a bit set*/
	std::vector<double> bm_logic_cuts_pool_phi_val;

	void bm_add_elem_to_logic_cut(uint64_t* set, unsigned int elem);
	bool bm_is_logic_cut_already_added(uint64_t* cut, int* index);
	void bm_add_logic_cut_to_pool(uint64_t* cut, const double* phi);

	std::vector<int> bm_branch_priority_score;


	/*Data structure for the additional variables CGG to have a starting Lower Bound.*/
	/*Second stage problem part*/
	Matrix2D<unsigned int>bm_v_y_xib_i; /*Second stage*/
	Matrix3D<unsigned int>bm_v_x_xib_ij; /*Routing*/
	Matrix2D<unsigned int> bm_v_lamb_xib_i; /* Dual of matrix A. One vector lambda for each xi_bar and for each constraint in A */
	Matrix2D<unsigned int> bm_v_sig_xib_i; /*Dual of xi_bar = xi*/
	
	/*Normal cuts vs INFO cuts*/
	bool bm_use_info_cuts;
    /*Data structure to store the scenarios*/
	bool bm_use_scenarios_ccg_xib;
	unsigned int num_ccg_xib;
	Matrix2D<double> ccg_set_xib;
	std::vector<unsigned int> cnstr_ccg_xi; /*index of the constraints. To check which xi is tight.*/
	bool store_xbar();
	bool bm_add_ccg_lb();/*This add variables and constraints for each stored xi*/
	bool bm_add_ccg_columns();
	bool bm_add_ccg_constraints();
	bool bm_separate_heuristic_subtours(CPXCALLBACKCONTEXTptr context, CPXLONG contextid, double var_treshold, bool* cut_added, int xi_bar);
	bool bm_separate_exact_subtuors(CPXCALLBACKCONTEXTptr context, CPXLONG contextid, bool* cut_added, int xi_bar);


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

	/*Reading Alrijne Case matrix*/

	bool read_distance_matrix(FILE* inp);

	public: 
		bool read_instance(char* file, char* config_file);

		bool init_data_structures();

		bool clean_data_structure();


		bool define_uncertainty_set();

		bool define_uncertainty_set_type_1();

		bool define_uncertainty_set_type_2();

		bool define_uncertainty_set_type_3();

		bool define_uncertainty_set_type_4();

		bool load_and_init(std::vector<bool> w_val);

		bool run_sub_problem_solver(double* obj_val);

		bool subproblem_solver_update_discovery_parameters(const std::vector<bool>* w_val);

		bool run_exact_method(double* objective_value);

		bool print_solution_file();

		bool print_solution_latex();


};

