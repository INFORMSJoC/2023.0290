#pragma once
#include "Utilities.h"
#include "ilcplex/cplex.h"
#include "ConfigFile.h"
#include <filesystem>

#define EPS 0.1
#define mEPS 0.0001 /*Small eps*/
#define cEPS 0.005 /* cut Small eps*/

#define STARTING_POLICY_SET_SIZE 100
#define MAX_NUM_POLICY 100000 /*Not used now*/

#define MAX_NUM_CNSTR_PER_C_US 10

#define TRESHOLD_MAX 0.5 /*Subtour sepation*/
#define TRESHOLD_STEP 0.25 /*Subtour sepation*/

#define BIG_M_SIGMA 100

/*Treshold for cut EXACT algo*/
#define EXA_EPS 1e-6
/*Teshodol stopping creteria subproblem (difference between Master and Sub in CG)*/
#define EXA_EPS_conv 1e-5

class ExaDCCG
{

	typedef struct master_problem {

		CPXENVptr cpx_env;
		CPXLPptr  cpx_lp;


		/*Parameters*/
		unsigned short card_Y; /*cardinality of Y*/
		Matrix2D<unsigned int> Y; // Each row is a solution y. You query solution element i of y Y(y,i).



		std::vector<unsigned short> v_w_i; /*one w for each i in N_w*/
		std::vector<unsigned short> v_alpha_y; /*One alpha for each policy y in Y*/
		Matrix2D<unsigned int> v_beta_y_l; /*One beta for each policy y in Y and for each constraints of US*/
		Matrix2D<unsigned int> v_gamma_y_i; /*One gamma for each policy y in Y and for each i in N_w*/
		std::vector<unsigned int> v_Beta_l;/*One Beta for each constraint in US. this is \pi in the picture*/

		std::vector<unsigned int> c_dual_xibar; /*Index of constraints associated to xi bar in order to get the */

		/*Solution*/
		std::vector<double> curr_sol;
		std::vector<double> curr_duals;
		std::vector<double> curr_xibar_i; /*Current dual variable xibar*/
		double curr_obj_val;
		double curr_solver_gap;


	}t_master_problem;

	typedef struct slave_problem {

		CPXENVptr cpx_env;
		CPXLPptr  cpx_lp;

		std::vector<double>  xibar_i; /*Parameters it gets from the master. One for each i \in N_w*/


		std::vector<unsigned short> v_w_i; /*one for each observable item*/
		std::vector<unsigned short> v_y_i; /*one foe each i in N_w*/

		Matrix2D<unsigned short> v_x_ij; /*node i, node j, OP specific*/
		std::vector<unsigned short> v_beta; /*dual Uncertainty set constraints*/
		std::vector<unsigned short> v_gamma; /*dual non-anticipativity constraints*/

		/*Solution*/
		std::vector<double> curr_sol;
		double curr_obj_val;

	}t_slave_problem;


private:
	/*To write models*/
	std::vector<int> matind;
	std::vector<double> matval;
	std::vector<char> matchar;

	/*Master problem*/
	master_problem mp_mod;
	bool mp_build_model();
	bool mp_set_cplex();
	bool mp_add_vars();
	bool mp_add_cnstr();
	bool mp_solve_model();
	bool mp_obtain_linear_relaxation();
	bool mp_solve_linear_relaxation();
	bool mp_get_xibar(); /*get the xibar by solving the linear relaxation with fixed w.*/
	bool mp_solver_master_problem(); /*it wraps everything for the master.*/

	slave_problem sp_mod;
	bool sp_build_model();
	bool sp_set_cplex();
	bool sp_add_vars();
	bool sp_add_cnstr();
	bool sp_update_model(); /*update xi_prime*/
	bool sp_solve_model();
	static int CPXPUBLIC sp_general_callback(CPXCALLBACKCONTEXTptr context, CPXLONG contextid, void* userdata);
	bool separate_heuristic_subtours(CPXCALLBACKCONTEXTptr context, CPXLONG contextid, double var_treshold, bool* cut_added);
	bool separate_exact_subtuors(CPXCALLBACKCONTEXTptr context, CPXLONG contextid, bool* cut_added);
	bool separate_routing_logic_cuts(CPXCALLBACKCONTEXTptr context, CPXLONG contextid, bool* cut_added);
	


	/*Dual CCG algorithm*/
	bool init_CCG_algorithm(); /*we init the parameters here and we do the first iteraion*/
	bool generate_initial_xibar(); /*each xibar random between lb and ub. This is specific for the US we are using.*/
	bool generate_initial_y(); /*specific for op.*/

	/*stats*/
	double ccg_time;
	double ccg_gap;
	double ccg_iter;



	/*
Instance info
*/
	Configuration::ConfigFile* cfg;

	char inst_name[1024];
	char network_name[1024];
	FILE* out_file;
	char results_folder[1024];



	/*From configuration*/
	unsigned short uncertainty_set_type; /*Uncertainty set types: 1 or 2 for now*/
	double uncertainty_parameter; /*Uncertainty set parameter. For unceertainty set type 1 is the budget. For US4 is the theta, for other US not used*/

	bool decision_depende_discovery; /*true for decision dependent discovery version. False for w = e.*/
	double max_number_discovery_fraction; /* \in [0,1]  ceil(max_number_discovery_fraction * num_obeservable_elements) = max_num_disc*/
	/********************/


	/*Data containers*/
	unsigned short num_us_constr;
	Matrix2D<double> A; // Uncertainty set for xi vars.
	std::vector<double> b; // Rhs for uncertainty set.
	/*Obejective function object*/
	std::vector<double>ppxi; // Function of \xi in the OP.
	std::vector<double> pp; // part not function of P.

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





	/*Methods*/

	bool init_cplex(CPXENVptr* cpx_env,
		CPXLPptr* cpx_lp);
	int  checkCPXstatus(int status, CPXENVptr* cpx_env,
		CPXLPptr* cpx_lp);
	bool free_cplex(CPXENVptr* cpx_env_op,
		CPXLPptr* cpx_lp_op);

	bool init_out_file();


	/*Reading Alrijne Case matrix*/

	bool read_distance_matrix(FILE* inp);



	/*Separations methods*/
	un_graph_t sep_un_graph;
	di_graph_t sep_di_graph;



	queue_t sep_queue;
	bool init_queue();
	bool queue_is_empty();
	bool push_back(int el);
	int pop_front();
	// Connected components
	int build_sol_graph(const double* curr_sol, double val_treshold);
	unsigned int connected_components();
	void connected_components_recursive(unsigned int vertex,
		unsigned int component);
	// Max-Flow with Dinic's Alg.
	bool build_sol_directed_graph(const double* curr_sol, double val_treshold);// create the directed graph from curr sol.
	void add_archs(int a1, int a2, int q);
	bool bfs_di_graph(int s, int t); // breadth first search
	int send_flow_di_graph(int s, int flow, int t, short* start);
	int dinic_max_flow(int s, int t);
	// min cut: rectrive the min cut from the max_flow.
	/*Min cuts set*/
	std::vector<bool> sep_mcut_set_1;
	std::vector<bool> sep_mcut_set_2;
	std::vector<unsigned short> node_subtour_subset;


public:
	bool read_instance(char* file);

	bool init_data_structures();

	bool clean_data_structure();

	/*run algo*/
	bool run_dual_ccg_algorithm();




	bool define_uncertainty_set();

	bool define_uncertainty_set_type_1();

	bool define_uncertainty_set_type_2();

	bool define_uncertainty_set_type_3();

	bool define_uncertainty_set_type_4();

};

