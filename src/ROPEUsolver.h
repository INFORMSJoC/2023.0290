#pragma once
#include "Utilities.h"
#include "ilcplex/cplex.h"
#include "ConfigFile.h"
#include <filesystem>

/****************************************/
//#define MAX_NUM_DISCOVERY_FRACTION 0.5 -->set using config file
//#define DDID_VERSION true -->set using config file
//#define NUM_K 3 --> Set using config file
//#define US_TYPE 1 -->set using config file
//#define US_BUDGET 3 -->set using config file
/***************************************/

#define BIG_M_GAMMA 1000/*big M for gamma linearization*/
#define EPS 0.1
#define mEPS 0.0001 /*Small eps*/
#define cEPS 0.005 /* cut Small eps*/



/*B&C*/
#define BandC_NUM_CUT_PER_NODE 500
#define BandC_MIN_COEFF_LP_IMPR 0.01
#define BandC_MAX_NUM_ITER_NODE 5
#define BandC_SET_BRANCHING_PRIORITY true

#define TRESHOLD_MAX 0.5 /*Subtour sepation*/
#define TRESHOLD_STEP 0.25 /*Subtour sepation*/

class ROPEUsolver
{
private:

	/*
	Instance info
	*/
	char inst_name[1024];
	char network_name[1024];
	FILE* out_file = NULL;
	char results_folder[1024];
	
	Configuration::ConfigFile* cfg;

	/*From configuration*/
	unsigned short K; /*Number of policies*/
	unsigned short uncertainty_set_type; /*Uncertainty set types: 1 or 2 for now*/
	double uncertainty_parameter; /*Uncertainty set parameter. For unceertainty set type 1 is the budget. For US4 is the theta, for other US not used*/

	bool decision_depende_discovery; /*true for decision dependent discovery version. False for w = e.*/
	double max_number_discovery_fraction; /* \in [0,1]  ceil(max_number_discovery_fraction * num_obeservable_elements) = max_num_disc*/
	/********************/

	/***Algo configuration******/
	bool alpha_bound_reduction;
	bool alpha_symmetry_breaking;
	bool xtilde_cuts;
	bool best_scenario_cuts;

	/*****************************/
	double big_m_gamma;



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


	/*Cplex structures for the master problem.*/
	CPXENVptr cpx_env;
	CPXLPptr  cpx_lp;

	/*To write constraints*/
	std::vector<int> matind;
	std::vector<int> branch_priority_score;
	std::vector<double> matval;
	std::vector<double> curr_sol;

	/*Variables container*/
	// Phi var
	unsigned short v_phi;

	std::vector<unsigned short> v_w_n; /*Node*/
	Matrix2D<unsigned short> v_y_kn; /*K-pol, node*/
	Matrix3D<unsigned short> v_x_kij; /*K-pol, node i, node j*/
	
	/*Duals*/
	std::vector<unsigned short> v_alpha_k; // query with k.
	std::vector<unsigned short> v_beta_l; // query with num_con.
	Matrix2D<unsigned short> v_beta_kl;  // query with K, num_con.
	Matrix2D<unsigned short> v_gamma_kn; // quesry with k, nu, nodes.

	/*Bilinear*/
	Matrix2D<unsigned short> v_t_gamma_kn; // gamma times w.
	Matrix2D<unsigned short> v_t_y_kn; // alpha times y. Components of y_k are moltiplied only for alpha k. You can then query as a y variables. 
	Matrix3D<unsigned short> v_t_x_kij; // alpha times x.


	/*Uncertanty set */
	Matrix2D<double> A; // Uncertainty set for xi vars.
	Matrix2D<double> G; // Uncertainty set for zeta vars size L_xi times Nzeta.
	std::vector<double> b; // Rhs for uncertainty set.
	std::vector<double> t_y_obj_coeff_val; // Coefficient in the objective constraints.
	std::vector<double> t_y_uc_coeff_val; // Coefficient in the balance constraints for tilde_y.


	std::vector<double> best_scenario_profit_i; // Minimum profit according to the used US. Yhis is used for the LB.
	unsigned short num_us_constr;


	/*B&C tools*/
	long int num_iter_curr_node;
	long int num_cut_curr_iter;
	long int prev_node;
	long int curr_node;
	double lp_prev_iter;
	void reset_branc_and_cut_indicators();

	/*
	Alpha bounds for VI
	*/
	std::vector<double> ub_alpha_k; // UB alpha_k
	std::vector<double> lb_alpha_k; // LB alpha_k

	/*Generic methods*/
	int  checkCPXstatus(int status);
	bool free_cplex();


	bool init_out_file();

	bool init_cplex();

	bool set_cplex();

	bool define_alpha_bounds();

	bool add_deterministic_variables();

	bool add_dualization_variables();

	bool add_linearization_variables();

	bool add_deterministic_constraints();

	bool add_objective_constraint();

	bool add_valid_inequalities_objective_lower_bound_constraint();

	bool add_dualization_constraints();

	bool add_linearization_constraints_old();

	bool add_linearization_constraints();

	bool add_strengthened_linearization_constraints_old();

	bool add_strengthened_linearization_constraints();

	bool add_valid_inequalities_alpha_symmetries();

	bool add_RLT_valid_inequalties();

	

	/*
	Callback and cutting plane methods
	*/

	static int CPXPUBLIC general_callback(CPXCALLBACKCONTEXTptr context, CPXLONG contextid, void* userdata);


	bool separate_heuristic_subtours(CPXCALLBACKCONTEXTptr context, CPXLONG contextid,double var_treshold,bool* cut_added);
	bool separate_tilde_heuristic_subtours(CPXCALLBACKCONTEXTptr context, CPXLONG contextid, double var_treshold, bool* cut_added);
	bool separate_exact_subtuors(CPXCALLBACKCONTEXTptr context, CPXLONG contextid,bool* cut_added);
	bool separate_tilde_exact_subtours(CPXCALLBACKCONTEXTptr context, CPXLONG contextid, bool* cut_added);
	bool separate_logic_cuts(CPXCALLBACKCONTEXTptr context, CPXLONG contextid, bool* cut_added);
	bool separate_McCormick_x_alpha(CPXCALLBACKCONTEXTptr context, CPXLONG contextid, bool* cut_added);




	/*Separations methods*/
	un_graph_t sep_un_graph;
	di_graph_t sep_di_graph;



	queue_t sep_queue;
	bool init_queue();
	bool queue_is_empty();
	bool push_back(int el);
	int pop_front();


	// Connected components
	 int build_sol_graph(const double* curr_sol, int k_pol, double val_treshold);
	 int build_tilde_sol_graph(const double* curr_sol, int k_pol, double val_treshold);


	unsigned int connected_components();

	void connected_components_recursive(unsigned int vertex,
		unsigned int component);

	// Max-Flow with Dinic's Alg.
	bool build_sol_directed_graph(const double* curr_sol, int k);// create the directed graph from curr sol.
	bool build_tilde_sol_directed_graph(const double* curr_sol, int k);// create the directed graph from curr sol.
	void add_archs(int a1, int a2, int q);
	bool bfs_di_graph(int s, int t); // breadth first search
	int send_flow_di_graph(int s, int flow, int t, short* start);
	int dinic_max_flow(int s, int t);
	/*Min cuts set*/
	 std::vector<bool> sep_mcut_set_1;
	 std::vector<bool> sep_mcut_set_2;
	 std::vector<unsigned short> node_subtour_subset;


	 /*Callback*/
	 bool root_node_processed;
	 

	 bool read_distance_matrix(FILE* inp);

public:
	

	// Orienteering
	bool read_instance(char* file, char* config_file);

	bool init_data_structures();

	bool clean_data_structure();

	bool define_uncertainty_set();

	bool define_uncertainty_set_type_1();

	bool define_uncertainty_set_type_2();

	bool define_uncertainty_set_type_3();

	bool define_uncertainty_set_type_4(); 


	bool print_debug_curr_sol();

	bool print_solution_file();

	bool print_solution_latex();

	bool build_and_solve_model();

	bool get_w_values(std::vector<bool>& w_val);


	double exact_w_evaluation;

	double best_ub;
	double best_lb;
	double cplex_gap;
	double comput_time;
	double root_node_lb;
	long unsigned int bandc_nodes;

	long unsigned int num_dvar_cuts; /*cuts on det variables*/
	long unsigned int num_tildevar_cuts; /*cuts on tilde variables, not mccormic*/
	long unsigned int num_dyna_mccormik; /*Dynamically added Mccormic*/

};

