#include "ExaDCCG.h"

bool ExaDCCG::init_cplex(CPXENVptr* cpx_env, CPXLPptr* cpx_lp)
{
	int status = 0;

	// Open env.
	*cpx_env = CPXopenCPLEX(&status);
	if (this->checkCPXstatus(status, cpx_env, cpx_lp)) return false;
	*cpx_lp = CPXcreateprob(*cpx_env, &status, inst_name);
	if (this->checkCPXstatus(status, cpx_env, cpx_lp)) return false;

	return true;
}

int ExaDCCG::checkCPXstatus(int status, CPXENVptr* cpx_env, CPXLPptr* cpx_lp)
{
	char errmsg[CPXMESSAGEBUFSIZE];
	if (status == 0) return 0;

	CPXgeterrorstring(*cpx_env, status, errmsg);
	printf(" %s \n", errmsg);
	return status;
}


bool ExaDCCG::free_cplex(CPXENVptr* cpx_env, CPXLPptr* cpx_lp)
{
	int status = 0;
	status = CPXfreeprob(*cpx_env, cpx_lp);
	if (checkCPXstatus(status, cpx_env, cpx_lp)) return false;
	status = CPXcloseCPLEX(cpx_env);
	if (checkCPXstatus(status, cpx_env, cpx_lp)) return false;

	return true;
}



/*Subtours elimination constraints*/

inline bool ExaDCCG::init_queue()
{
	sep_queue.num_el = 0;
	sep_queue.next = 0;
	return true;
}

inline bool ExaDCCG::queue_is_empty()
{
	if (sep_queue.num_el == 0
		||
		sep_queue.num_el == sep_queue.next)
		return true;
	else
		return false;
}

inline bool ExaDCCG::push_back(int el)
{
	sep_queue.elem[sep_queue.num_el] = el;
	++sep_queue.num_el;
	return true;
}

inline int ExaDCCG::pop_front()
{
	int elem = sep_queue.elem[sep_queue.next];
	++sep_queue.next;
	return elem;
}

bool ExaDCCG::build_sol_directed_graph(const double* curr_sol, double val_treshold)
{
	sep_un_graph.num_nodes = 0;
	sep_un_graph.num_edges = 0;

	sep_un_graph.num_nodes = 1;
	sep_un_graph.nodes_map[sep_un_graph.num_nodes - 1] = start_node;
	//g_start_node;

	for (int i = start_node + 1; i < end_node; ++i)
		if (curr_sol[sp_mod.v_y_i[i]] > mEPS)
		{

			sep_un_graph.nodes_map[sep_un_graph.num_nodes] = i;
			++sep_un_graph.num_nodes;
		}

	sep_un_graph.nodes_map[sep_un_graph.num_nodes] = end_node;
	++sep_un_graph.num_nodes;

	for (int i = 0; i < sep_un_graph.num_nodes - 1; ++i)
		for (int j = i + 1; j < sep_un_graph.num_nodes; ++j)
			if (curr_sol[sp_mod.v_x_ij(sep_un_graph.nodes_map[i], sep_un_graph.nodes_map[j])] > mEPS + val_treshold)
				//if (curr_sol[g_v_xk_kij[k][sep_un_graph.nodes_map[i]][sep_un_graph.nodes_map[j]]] > mEPS)
			{

				// Add edge.
				sep_un_graph.edge[sep_un_graph.num_edges].first = i;
				sep_un_graph.edge[sep_un_graph.num_edges].second = j;
				++sep_un_graph.num_edges;
			}
	return 0;
}




inline void ExaDCCG::add_archs(int a1, int a2, int q)
{
	int pos_a1 = sep_di_graph.adj_size[a1];
	int pos_a2 = sep_di_graph.adj_size[a2];

	sep_di_graph.adj_archs(a1, pos_a1).a1 = a1;
	sep_di_graph.adj_archs(a1, pos_a1).a2 = a2;

	sep_di_graph.adj_archs(a1, pos_a1).q = q;
	sep_di_graph.adj_archs(a1, pos_a1).rev = pos_a2;
	sep_di_graph.adj_archs(a1, pos_a1).flow = 0; // init flow to 0

	// Now the reverse
	sep_di_graph.adj_archs(a2, pos_a2).a1 = a2;
	sep_di_graph.adj_archs(a2, pos_a2).a2 = a1;

	sep_di_graph.adj_archs(a2, pos_a2).q = q;
	sep_di_graph.adj_archs(a2, pos_a2).rev = pos_a1;
	sep_di_graph.adj_archs(a2, pos_a2).flow = 0; // init flow to 0


	++sep_di_graph.adj_size[a1];
	++sep_di_graph.adj_size[a2];
}

inline bool ExaDCCG::bfs_di_graph(int s, int t)
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
			level = sep_di_graph.level[sep_di_graph.adj_archs(u, a).a2];
			flow = sep_di_graph.adj_archs(u, a).flow;

			if (level < 0 && flow < sep_di_graph.adj_archs(u, a).q)
			{
				sep_di_graph.level[sep_di_graph.adj_archs(u, a).a2] = sep_di_graph.level[u] + 1;

				push_back(sep_di_graph.adj_archs(u, a).a2);

			}
		}

	}

	return sep_di_graph.level[t] < 0 ? false : true;
}

inline int ExaDCCG::send_flow_di_graph(int s, int flow, int t, short* start)
{
	// Sink reached 
	if (s == t)
		return flow;

	// Traverse adiacent edges
	for (; start[s] < sep_di_graph.adj_size[s]; ++start[s])
	{
		arch_t* e = &(sep_di_graph.adj_archs(s, start[s]));

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
				sep_di_graph.adj_archs(e->a2, e->rev).flow -= temp_flow;
				return temp_flow;

			}

		}

	}


	return 0;
}

inline int ExaDCCG::dinic_max_flow(int s, int t)
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

int ExaDCCG::build_sol_graph(const double* curr_sol, double val_treshold)
{
	sep_un_graph.num_nodes = 0;
	sep_un_graph.num_edges = 0;

	sep_un_graph.num_nodes = 1;
	sep_un_graph.nodes_map[sep_un_graph.num_nodes - 1] = start_node;
	//g_start_node;

	for (int i = start_node + 1; i < end_node; ++i)
		if (curr_sol[this->sp_mod.v_y_i[i]] > mEPS)
		{

			sep_un_graph.nodes_map[sep_un_graph.num_nodes] = i;
			++sep_un_graph.num_nodes;
		}

	sep_un_graph.nodes_map[sep_un_graph.num_nodes] = end_node;
	++sep_un_graph.num_nodes;

	for (int i = 0; i < sep_un_graph.num_nodes - 1; ++i)
		for (int j = i + 1; j < sep_un_graph.num_nodes; ++j)
			if (curr_sol[this->sp_mod.v_x_ij(sep_un_graph.nodes_map[i], sep_un_graph.nodes_map[j])] > mEPS + val_treshold)
				//if (curr_sol[g_v_xk_kij[k][sep_un_graph.nodes_map[i]][sep_un_graph.nodes_map[j]]] > mEPS)
			{

				// Add edge.
				sep_un_graph.edge[sep_un_graph.num_edges].first = i;
				sep_un_graph.edge[sep_un_graph.num_edges].second = j;
				++sep_un_graph.num_edges;
			}
	return 0;
}

unsigned int ExaDCCG::connected_components()
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

inline void ExaDCCG::connected_components_recursive(unsigned int vertex, unsigned int component)
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

/*Reading instance and initialize it*/
bool ExaDCCG::init_out_file()
{
	if (out_file != NULL)
		return true;


	if (out_file == NULL)
	{
		char locname[1024];
		char ss[1024];

		snprintf(locname, sizeof(locname), "%s%s%s%s", results_folder, "\\", "Res_ExaDCCG_", inst_name);
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

bool ExaDCCG::read_distance_matrix(FILE* inp)
{



	/*double coll_time_ij = alrijne_node_collecting_time * 100.00;

	double coll_time_dep_i = alrijne_node_collecting_time * 100.00 * 0.5;*/

	double coll_time_ij = 0;

	double coll_time_dep_i = 0;

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
	///*Print tij here: */
	//std::filebuf fb;
	//fb.open("alrijne_matrix.txt", std::ios::out);
	//std::ostream os(&fb);
	//

	//for (i = 0; i < this->num_nodes - 1; i++)
	//{
	//	for (j = 0; j < i; j++)
	//	{
	//		os << t_ij(j, i) /100.0 << "\t";
	//	}
	//	os << "\n";
	//}
	//fb.close();
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
		tot_deterministic_profit += node_det_profit[i];
	}

	return true;
}

bool ExaDCCG::read_instance(char* file)
{
	FILE* inp;
	errno_t err;
	char line[1024];
	char* ss;
	char cline[1024];

	bool read_matrix = false;

	cfg = new Configuration::ConfigFile("config.cfg");

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


	snprintf(inst_name, sizeof(inst_name), "%s_%s%d_%s%d_%s%d_%s%d", inst_name, "D", decision_depende_discovery, "md0", fract, "US", uncertainty_set_type, "USpar", us_par);

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


	if (!read_matrix)
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
				t_ij(i, j) =
					(int)(100.00 * std::sqrt(std::pow(node_x_coord[i] - node_x_coord[j], 2.0) + std::pow(node_y_coord[i] - node_y_coord[j], 2.0)));

				t_ij(j, i) = t_ij(i, j);

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



		/*Fake profit: this is updated in the matrix*/
		node_det_profit.assign(this->num_nodes, 1);
		node_det_profit[start_node] = 0;
		node_det_profit[end_node] = 0;

		//	tot_deterministic_profit = 1; /*For Alrijne ww consider 100% the toal that is constant.*/
		tot_deterministic_profit = 0; // we now hae profit
		//this->num_prof_nodes; /* 1 time num prof nodes */
		disc_cost_n.assign(this->num_nodes, 0);

		// Travel times
		this->t_ij.assignMatrix2D(this->num_nodes, this->num_nodes, 0);

		// Node existence
		this->edge_exists.assignMatrix2D(this->num_nodes, this->num_nodes, false);

		this->read_distance_matrix(inp);

	}



	return true;
}

bool ExaDCCG::init_data_structures()
{

	//snprintf(results_folder, sizeof(results_folder), "%s", RESULTS_FOLDER);

	//cfg = new Configuration::ConfigFile("config.cfg");

	stringstream ss;
	string out_dir;
	ss << cfg->getValueOfKey<string>("RESULTS_FOLDER");
	out_dir = ss.str();

	if (!std::filesystem::create_directory(ss.str()))
		cout << "Warning: " << strerror(errno) << endl;

	snprintf(results_folder, sizeof(results_folder), "%s", ss.str().c_str());




	/*uncertainty set data structure: for now I assume the uncertainty set is defined by num_nodes*2 constr + 1 (one for each node plus budget.)*/

	A.assignMatrix2D((MAX_NUM_CNSTR_PER_C_US * this->num_nodes) + 1, this->num_nodes, 0.00);

	b.resize((MAX_NUM_CNSTR_PER_C_US * this->num_nodes) + 1);

	ppxi.resize(this->num_nodes);

	pp.resize(this->num_nodes); // part not function of P.




	unsigned int num_vrs_slave = 0;
	unsigned int num_vrs_master = 0;



	this->sp_mod.v_y_i.resize(this->num_nodes);
	num_vrs_slave += sp_mod.v_y_i.size();

	this->sp_mod.v_x_ij.resizeMatrix2D(this->num_nodes, this->num_nodes);
	num_vrs_slave += sp_mod.v_x_ij.getNumElem();

	/*Exchange 4 with the max number of constrasints*/
	this->sp_mod.v_beta.resize(b.size() * this->num_nodes);
	num_vrs_slave += sp_mod.v_beta.size();

	this->sp_mod.v_gamma.resize(this->num_nodes);
	num_vrs_slave += sp_mod.v_gamma.size();

	this->sp_mod.v_w_i.resize(this->num_nodes);
	num_vrs_slave += sp_mod.v_w_i.size();

	/*solution*/
	sp_mod.curr_sol.resize(num_vrs_slave);



	/*Master Problem*/
	mp_mod.v_w_i.resize(this->num_nodes);
	num_vrs_master += mp_mod.v_w_i.size();
	mp_mod.card_Y = 0;
	mp_mod.v_alpha_y.resize(STARTING_POLICY_SET_SIZE);
	mp_mod.Y.resizeMatrix2D(STARTING_POLICY_SET_SIZE, this->num_nodes);
	num_vrs_master += mp_mod.Y.getNumElem();
	mp_mod.v_beta_y_l.resizeMatrix2D(STARTING_POLICY_SET_SIZE, b.size());
	num_vrs_master += mp_mod.v_beta_y_l.getNumElem();
	mp_mod.v_gamma_y_i.resizeMatrix2D(STARTING_POLICY_SET_SIZE, this->num_nodes);
	num_vrs_master += mp_mod.v_gamma_y_i.getNumElem();
	mp_mod.v_Beta_l.resize(b.size());
	num_vrs_master += mp_mod.v_Beta_l.size();

	mp_mod.curr_sol.resize(num_vrs_master);
	mp_mod.curr_duals.resize(num_vrs_master);
	mp_mod.c_dual_xibar.resize(num_vrs_master);


	mp_mod.curr_xibar_i.resize(this->num_nodes);

	/*matind and matval*/
	matind.resize(std::max(num_vrs_slave, num_vrs_master));
	matval.resize(std::max(num_vrs_slave, num_vrs_master));
	matchar.resize(std::max(num_vrs_slave, num_vrs_master));



	/*Separations Objects*/
	sep_mcut_set_1.assign(this->num_nodes, 0);
	sep_mcut_set_2.assign(this->num_nodes, 0);
	sep_queue.elem.assign(MAX_NUM_ELE_IN_Q, 0);
	sep_un_graph.edge.resize(MAX_NUM_EDGES);
	sep_un_graph.nodes_map.resize(this->num_nodes);
	sep_un_graph.conn_comp.resize(this->num_nodes);
	sep_di_graph.ver_label.resize(this->num_nodes);
	sep_di_graph.level.resize(this->num_nodes);
	sep_di_graph.adj_archs.resizeMatrix2D(this->num_nodes, 2 * this->num_nodes);
	sep_di_graph.adj_size.resize(this->num_nodes);
	node_subtour_subset.resize(this->num_nodes + 1);

	
	if (!this->define_uncertainty_set())
	{
		return false;
	}

	return true;
}

bool ExaDCCG::clean_data_structure()
{
	delete cfg;

	/*Clean cplex op*/
	int status = 0;
	//status = CPXfreeprob(this->cpx_env_op, &this->cpx_lp_op);
	//if (checkCPXstatus(status, &this->cpx_env_op, &this->cpx_lp_op)) return false;
	//status = CPXcloseCPLEX(&this->cpx_env_op);
	//if (checkCPXstatus(status, &this->cpx_env_op, &this->cpx_lp_op)) return false;

	///*Clean mlp*/
	//status = CPXfreeprob(this->cpx_env_mlp, &this->cpx_lp_mlp);
	//if (checkCPXstatus(status, &this->cpx_env_mlp, &this->cpx_lp_mlp)) return false;
	//status = CPXcloseCPLEX(&this->cpx_env_mlp);
	//if (checkCPXstatus(status, &this->cpx_env_mlp, &this->cpx_lp_mlp)) return false;


	///*Clean bm*/
	//status = CPXfreeprob(this->cpx_env_bm, &this->cpx_lp_bm);
	//if (checkCPXstatus(status, &this->cpx_env_bm, &this->cpx_lp_bm)) return false;
	//status = CPXcloseCPLEX(&this->cpx_env_bm);
	//if (checkCPXstatus(status, &this->cpx_env_bm, &this->cpx_lp_bm)) return false;

	if (out_file != NULL)
		fclose(out_file);





	return true;
}

bool ExaDCCG::run_dual_ccg_algorithm()
{
	bool check = true;

	unsigned int iter = 0;
	double start = 0;
	double end = 0;
	double curr_time;
	
	double prova;


	check = this->init_data_structures();
	if (!check)
	{
		goto TERMINATE;
	}


	start = clock();
	check = this->init_CCG_algorithm();
	if (!check)
	{
		goto TERMINATE;
	}
	++iter;



		while (
			fabs(sp_mod.curr_obj_val - mp_mod.curr_obj_val) > mEPS
			)
	{
	
		check =  mp_solver_master_problem();
		if (!check)
		{
			goto TERMINATE;
		}

		check = sp_update_model();
		if (!check)
		{
			goto TERMINATE;
		}

		check = sp_solve_model();
		if (!check)
		{
			goto TERMINATE;
		}
		++iter;
		end = clock();

		curr_time = (double)(end - start) / (double)CLK_TCK;

		std::cout << "Iter " << iter << "| UB = " << mp_mod.curr_obj_val << " | LB =  " << sp_mod.curr_obj_val << " | Time = " << curr_time<<std::endl;
	}
	end = clock();
	ccg_time = (double)(end - start) / (double)CLK_TCK;
	ccg_iter = iter;

TERMINATE:
	return check;

}

bool ExaDCCG::define_uncertainty_set()
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

bool ExaDCCG::define_uncertainty_set_type_1()
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


	///*Now we define the US cost*/
	////- y(p - p*0.5*\xi) ->  y(-p + p*0.5*\xi)
	for (int i = start_node; i <= this->end_node; ++i)
	{
		this->pp[i] = -this->node_det_profit[i];

		this->ppxi[i] = this->node_det_profit[i] * 0.5;


	}



	return true;
}

bool ExaDCCG::define_uncertainty_set_type_2()
{
	int node, col, row;
	double fraction = this->uncertainty_parameter;
	row = 0;
	/*US parama percentage of the total amount*/
	/*\xi in [0, fraction * B] and \sum_i \xi = B */
/*y(\xi)   B = 100% bijvoorbeeld*/
		 // First the upper bound xi <= 1: therefore we stop at Nxi;
	for (node = this->start_node; node <= this->end_node; ++node)
		if (this->node_is_profitable[node])
		{

			this->A(row, node) = 1;
			this->b[row] = fraction;
			//this->b[row] = fraction * this->tot_deterministic_profit;
			++row;


		}

	for (node = this->start_node; node <= this->end_node; ++node)
		if (this->node_is_profitable[node])
		{

			this->A(row, node) = -1;
			this->b[row] = 0;
			++row;


		}
	/*Last 2 constraint about the sum of the total is equal to 100% of the profit.*/
	//
	for (node = this->start_node; node <= this->end_node; ++node)
		if (this->node_is_profitable[node])
		{

			this->A(row, node) = 1;
		}
	//this->b[row] = this->tot_deterministic_profit;
	this->b[row] = 1 + 1e-7;
	++row;

	for (node = this->start_node; node <= this->end_node; ++node)
		if (this->node_is_profitable[node])
		{

			this->A(row, node) = -1;
		}
	//this->b[row] = -this->tot_deterministic_profit;
	this->b[row] = -(1 - 1e-7);
	++row;
	this->num_us_constr = row;


	/*Now we define the US cost*/
	// -y(\xi) where P is the total profit -> y(-\xi)
	for (int i = start_node; i <= this->end_node; ++i)
	{
		//	this->best_scenario_profit_i[i] = this->tot_deterministic_profit; /*\xi = 1 for each entry: then eveyrone has the max profit*/

		//	this->t_y_obj_coeff_val[i] = -0; // Coefficent without xi
		this->pp[i] = 0;
		//	this->t_y_uc_coeff_val[i] = -this->tot_deterministic_profit; // Coefficient with xi

		this->ppxi[i] = -1;
		//this->ppxi[i] = -this->tot_deterministic_profit;
	}


	return true;
}

bool ExaDCCG::define_uncertainty_set_type_3()
{
	int node, col, row;

	row = 0;

	double fraction = this->uncertainty_parameter;

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
	this->b[row] = 0;
	++row;

	/*Last constraints about the budget.*/
	//
	for (node = this->start_node; node <= this->end_node; ++node)
		if (this->node_is_profitable[node])
		{

			this->A(row, node) = -1;
		}
	//this->b[row] = -this->uncertainty_parameter;
	this->b[row] = 0;
	++row;



	this->num_us_constr = row;


	/*Now we define the US cost*/
	// -y(p - p*0.5*\xi) --> -yp + yp0.5\xi 
	for (int i = start_node; i <= this->end_node; ++i)
	{
		this->pp[i] = -this->node_det_profit[i];
		//this->t_y_obj_coeff_val[i] = -this->node_det_profit[i]; // Coefficent without xi

		//this->t_y_uc_coeff_val[i] = this->node_det_profit[i] * 0.5; // Coefficient with xi: the proft can be reduced at most by 50%
		this->ppxi[i] = this->node_det_profit[i] * fraction;

	}



	return true;
}

bool ExaDCCG::define_uncertainty_set_type_4()
{
	int node, col, row;


	/*\xi in [(1-teta) p, (1+teta)p] and \sum_i \xi = tot_profit*/
	//	//- y(*\xi)
	double theta = this->uncertainty_parameter;

	/*Now we define the US cost*/
//- y(\xi)
	for (int i = start_node; i <= this->end_node; ++i)
	{
		this->pp[i] = 0;

		this->ppxi[i] = -1;

		//	this->t_y_obj_coeff_val[i] = 0; // Coefficent without xi

		//	this->t_y_uc_coeff_val[i] = -1; // Coefficient with xi: the profit is exaclyt xi now

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
			this->b[row] = -(1 - theta) * this->node_det_profit[node];
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

bool ExaDCCG::mp_build_model()
{
	bool is_ok = true;

	is_ok = init_cplex(&this->mp_mod.cpx_env, &this->mp_mod.cpx_lp);
	if (!is_ok)
	{
		goto TERMINATE;
	}

	is_ok = mp_set_cplex();
	if (!is_ok)
	{
		goto TERMINATE;
	}

	is_ok = mp_add_vars();
	if (!is_ok)
	{
		goto TERMINATE;
	}

	is_ok = mp_add_cnstr();
	if (!is_ok)
	{
		goto TERMINATE;
	}


TERMINATE:
	return is_ok;
}

bool ExaDCCG::mp_set_cplex()
{
	bool mod_stat = true;

	int status = 0;

	char errbuf[CPXMESSAGEBUFSIZE];


	status = CPXsetintparam(this->mp_mod.cpx_env, CPX_PARAM_SCRIND, CPX_ON);
	if (checkCPXstatus(status, &mp_mod.cpx_env, &this->mp_mod.cpx_lp)) return false;
	status = CPXsetintparam(this->mp_mod.cpx_env, CPX_PARAM_THREADS, 1);
	if (checkCPXstatus(status, &mp_mod.cpx_env, &this->mp_mod.cpx_lp)) return false;
	status = CPXchgobjsen(this->mp_mod.cpx_env, this->mp_mod.cpx_lp, CPX_MIN);
	if (checkCPXstatus(status, &mp_mod.cpx_env, &mp_mod.cpx_lp)) return false;
	status = CPXsetintparam(this->mp_mod.cpx_env, CPX_PARAM_NUMERICALEMPHASIS, CPX_ON);
	if (checkCPXstatus(status, &mp_mod.cpx_env, &mp_mod.cpx_lp)) return false;
	status = CPXsetdblparam(this->mp_mod.cpx_env, CPXPARAM_Simplex_Tolerances_Feasibility, 1e-9);
	if (checkCPXstatus(status, &mp_mod.cpx_env, &mp_mod.cpx_lp)) return false;
	//status = CPXsetintparam(cpx_env, CPX_PARAM_PREIND, CPX_OFF);
	//if (checkCPXstatus(status))return false;

	/*To be extra numerical accurate*/
	status = CPXsetdblparam(this->mp_mod.cpx_env, CPXPARAM_Simplex_Tolerances_Optimality, 1e-9);
	if (checkCPXstatus(status, &mp_mod.cpx_env, &mp_mod.cpx_lp)) return false;

	status = CPXsetdblparam(this->mp_mod.cpx_env, CPXPARAM_Simplex_Tolerances_Feasibility, 1e-9);
	if (checkCPXstatus(status, &mp_mod.cpx_env, &mp_mod.cpx_lp)) return false;

	status = CPXsetdblparam(this->mp_mod.cpx_env, CPXPARAM_MIP_Tolerances_MIPGap, 1e-9);
	if (checkCPXstatus(status, &mp_mod.cpx_env, &mp_mod.cpx_lp)) return false;
	status = CPXsetdblparam(this->mp_mod.cpx_env, CPXPARAM_MIP_Tolerances_Integrality, 0);
	if (checkCPXstatus(status, &mp_mod.cpx_env, &mp_mod.cpx_lp)) return false;

	double time_limit = this->cfg->getValueOfKey<double>("TIME_LIMIT");

	status = CPXsetdblparam(mp_mod.cpx_env, CPXPARAM_TimeLimit, time_limit);
	if (checkCPXstatus(status, &mp_mod.cpx_env, &mp_mod.cpx_lp)) return false;

	
	return true;
}

bool ExaDCCG::mp_add_vars()
{
	unsigned int num_cols = CPXgetnumcols(this->mp_mod.cpx_env, this->mp_mod.cpx_lp);

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
	First stage variables W: discovery variables
	*/
	/*Only profitable nodes.*/
	for (int i = first_pnode; i <= last_pnode; ++i)
	{
		this->mp_mod.v_w_i[i] = num_cols;
		obj[0] = this->disc_cost_n[i];

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


		status = CPXnewcols(mp_mod.cpx_env, mp_mod.cpx_lp, 1, obj, lb, ub, vartype, varsname);
		if (checkCPXstatus(status, &mp_mod.cpx_env, &mp_mod.cpx_lp)) goto TERMINATE;
		++num_cols;
	}
	/*alpha variables*/
	for (int y = 0; y < this->mp_mod.card_Y; ++y)
	{
		this->mp_mod.v_alpha_y[y] = num_cols;
		obj[0] = 0;

		for (int i = first_pnode; i <= last_pnode; ++i)
		{
			obj[0] += this->pp[i] * this->mp_mod.Y(y, i);
		}
		
		vartype[0] = 'C';
		sprintf(varsname[0], "alpha_y%d", y);
		lb[0] = 0;
		ub[0] = 1;

		status = CPXnewcols(mp_mod.cpx_env, mp_mod.cpx_lp, 1, obj, lb, ub, NULL, varsname);
		if (checkCPXstatus(status, &mp_mod.cpx_env, &mp_mod.cpx_lp)) goto TERMINATE;
		++num_cols;
	}

	/*gamma variables*/
	for (int y = 0; y < this->mp_mod.card_Y; ++y)
		for (int i = first_pnode; i <= last_pnode; ++i)
		{
			this->mp_mod.v_gamma_y_i(y,i) = num_cols;
			obj[0] = 0;


			vartype[0] = 'C';
			sprintf(varsname[0], "gamma_y%d_i%d", y,i);

			lb[0] = - BIG_M_SIGMA;
			ub[0] =   BIG_M_SIGMA;

			status = CPXnewcols(mp_mod.cpx_env, mp_mod.cpx_lp, 1, obj, lb, ub, NULL, varsname);
			if (checkCPXstatus(status, &mp_mod.cpx_env, &mp_mod.cpx_lp)) goto TERMINATE;
			++num_cols;
		}

	/*beta variables*/
	for (int y = 0; y < this->mp_mod.card_Y; ++y)
		for (int l = 0; l < this->num_us_constr; ++l)
		{
			this->mp_mod.v_beta_y_l(y, l) = num_cols;
			obj[0] = this->b[l];


			vartype[0] = 'C';
			sprintf(varsname[0], "beta_y%d_i%d", y, l);

			lb[0] = 0;
			ub[0] = BIG_M_SIGMA;



			status = CPXnewcols(mp_mod.cpx_env, mp_mod.cpx_lp, 1, obj, lb, ub, NULL, varsname);
			if (checkCPXstatus(status, &mp_mod.cpx_env, &mp_mod.cpx_lp)) goto TERMINATE;
			++num_cols;
		}

	for (int l = 0; l < this->num_us_constr; ++l)
	{
		this->mp_mod.v_Beta_l[l] = num_cols;
		obj[0] = this->b[l];


		vartype[0] = 'C';
		sprintf(varsname[0], "Beta_l%d", l);

		lb[0] = 0;
		ub[0] = BIG_M_SIGMA;


		status = CPXnewcols(mp_mod.cpx_env, mp_mod.cpx_lp, 1, obj, lb, ub, NULL, varsname);
		if (checkCPXstatus(status, &mp_mod.cpx_env, &mp_mod.cpx_lp)) goto TERMINATE;
		++num_cols;
	}


TERMINATE:
	if (status)
		return false;
	else
		return true;
}

bool ExaDCCG::mp_add_cnstr()
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


	unsigned int num_rows = CPXgetnumrows(this->mp_mod.cpx_env, this->mp_mod.cpx_lp);


	
	if (matind.size() < CPXgetnumcols(mp_mod.cpx_env, mp_mod.cpx_lp))
	{
		matind.resize(CPXgetnumcols(mp_mod.cpx_env, mp_mod.cpx_lp));
		matval.resize(CPXgetnumcols(mp_mod.cpx_env, mp_mod.cpx_lp));

	}


	// Max number of discoveries: w \in W
	rhs[0] = this->max_num_disc;
	sense[0] = 'L';
	nzc = 0;
	sprintf(cnstrname[0], "Max_num_discoveries");

	for (int i = this->first_pnode; i <= this->last_pnode; ++i)
	{
		matind[nzc] = this->mp_mod.v_w_i[i];
		matval[nzc] = 1;

		++nzc;
	}
	status = CPXaddrows(mp_mod.cpx_env, mp_mod.cpx_lp, 0, 1, nzc, rhs, sense, matbeg, matind.data(), matval.data(), NULL, cnstrname);
	if (checkCPXstatus(status, &mp_mod.cpx_env, &mp_mod.cpx_lp)) goto TERMINATE;
	++num_rows;

	/*\sum_y \in Y alpha_y = 1 */
	rhs[0] = 1;
	sense[0] = 'E';
	nzc = 0;
	sprintf(cnstrname[0], "sum_alphs");

	for (int y = 0; y < this->mp_mod.card_Y; ++y)
	{
		matind[nzc] = this->mp_mod.v_alpha_y[y];
		matval[nzc] = 1;
		++nzc;
	}
	status = CPXaddrows(mp_mod.cpx_env, mp_mod.cpx_lp, 0, 1, nzc, rhs, sense, matbeg, matind.data(), matval.data(), NULL, cnstrname);
	if (checkCPXstatus(status, &mp_mod.cpx_env, &mp_mod.cpx_lp)) goto TERMINATE;
	++num_rows;


	/*
	\forall y \in Y
	\forall i \in N_w
	 \sum_l 1:NumRows  A_(l,i) * \beta^y_l = \gamma_y_i + \alpha_y * y_i 
	*/
	rhs[0] = 0;
	sense[0] = 'E';
	nzc = 0;
	for (int y = 0; y < this->mp_mod.card_Y; ++y)
		for (int i = this->first_pnode; i <= this->last_pnode; ++i)
		{
			nzc = 0;
			sprintf(cnstrname[0], "BalanceXi_y%d_i%d", y, i);


			matind[nzc] =   this->mp_mod.v_alpha_y[y];
			matval[nzc] = - (this->mp_mod.Y(y,i) * this->ppxi[i]);
			++nzc;

			matind[nzc] =   this->mp_mod.v_gamma_y_i(y, i);
			matval[nzc] = - 1;
			++nzc;


			for (int l = 0; l < this->num_us_constr; ++l)
			{
				matind[nzc] = this->mp_mod.v_beta_y_l(y, l);
				matval[nzc] = this->A(l, i);
				++nzc;
			}

			status = CPXaddrows(mp_mod.cpx_env, mp_mod.cpx_lp, 0, 1, nzc, rhs, sense, matbeg, matind.data(), matval.data(), NULL, cnstrname);
			if (checkCPXstatus(status, &mp_mod.cpx_env, &mp_mod.cpx_lp)) goto TERMINATE;
			++num_rows;
		}


	rhs[0] = 0;
	sense[0] = 'E';
	nzc = 0;
	for (int i = this->first_pnode; i <= this->last_pnode; ++i)
	{
		nzc = 0;
		sprintf(cnstrname[0], "BalanceXiBar_i%d", i);


		for (int y = 0; y < this->mp_mod.card_Y; ++y)
		{
			matind[nzc] = this->mp_mod.v_gamma_y_i(y, i);
			matval[nzc] = 1;
			++nzc;
		}

		for (int l = 0; l < this->num_us_constr; ++l)
		{
			matind[nzc] = this->mp_mod.v_Beta_l[l];
			matval[nzc] = this->A(l, i);
			++nzc;
		}

		status = CPXaddrows(mp_mod.cpx_env, mp_mod.cpx_lp, 0, 1, nzc, rhs, sense, matbeg, matind.data(), matval.data(), NULL, cnstrname);
		if (checkCPXstatus(status, &mp_mod.cpx_env, &mp_mod.cpx_lp)) goto TERMINATE;
		
		this->mp_mod.c_dual_xibar[i] = num_rows;
		++num_rows;
	}

	/*Linearization gamma - w*/
	for (int i = first_pnode; i <= last_pnode; ++i)
		for(int y = 0; y < mp_mod.card_Y; ++y)
	{
		nzc = 0;
		rhs[0] = 0;
		sense[0] = 'G';
		sprintf(cnstrname[0], "Linearize_gammaW_GQ_y%d_i%d", y,i);

		matind[nzc] = this->mp_mod.v_gamma_y_i(y, i);
		matval[nzc] = 1;
		++nzc;

		matind[nzc] = this->mp_mod.v_w_i[i];
		matval[nzc] = BIG_M_SIGMA;
		++nzc;


		status = CPXaddrows(mp_mod.cpx_env, mp_mod.cpx_lp, 0, 1, nzc, rhs, sense, matbeg, matind.data(), matval.data(), NULL, cnstrname);
		if (checkCPXstatus(status, &mp_mod.cpx_env, &mp_mod.cpx_lp)) goto TERMINATE;

		nzc = 0;
		rhs[0] = 0;
		sense[0] = 'L';
		sprintf(cnstrname[0], "Linearize_gammaW_LQ_y%d_i%d", y,i);

		matind[nzc] = this->mp_mod.v_gamma_y_i(y, i);
		matval[nzc] = 1;
		++nzc;

		matind[nzc] = this->mp_mod.v_w_i[i];
		matval[nzc] = -BIG_M_SIGMA;
		++nzc;
		status = CPXaddrows(mp_mod.cpx_env, mp_mod.cpx_lp, 0, 1, nzc, rhs, sense, matbeg, matind.data(), matval.data(), NULL, cnstrname);
		if (checkCPXstatus(status, &mp_mod.cpx_env, &mp_mod.cpx_lp)) goto TERMINATE;
	}

TERMINATE:
	if (status)
		return false;
	else
		return true;
}

bool ExaDCCG::mp_solve_model()
{
	bool is_ok = true;
	int status;
	int mip_status;
	CPXLONG      contextmask = 0;
	char errbuf[CPXMESSAGEBUFSIZE];
	int check_mip_stat;

	double start = 0;

	double end = 0;

	status = CPXwriteprob(this->mp_mod.cpx_env, this->mp_mod.cpx_lp, "pene.lp", "lp");
	if (checkCPXstatus(status, &mp_mod.cpx_env, &mp_mod.cpx_lp)) { is_ok = false; goto TERMINATE; };

	
	start = clock();
	status = CPXmipopt(this->mp_mod.cpx_env, this->mp_mod.cpx_lp);
	if (checkCPXstatus(status, &mp_mod.cpx_env, &mp_mod.cpx_lp)) goto TERMINATE;
	end = clock();
	//this->bm_branch_and_cut_time = (double)(end - start) / (double)CLK_TCK;


	check_mip_stat = CPXgetstat(mp_mod.cpx_env, mp_mod.cpx_lp);

	if (mp_mod.curr_sol.size() < CPXgetnumcols(mp_mod.cpx_env, mp_mod.cpx_lp))
	{
		mp_mod.curr_sol.resize(CPXgetnumcols(mp_mod.cpx_env, mp_mod.cpx_lp));
	}

	if (check_mip_stat != CPXMIP_TIME_LIM_INFEAS)
	{
		status = CPXgetobjval(mp_mod.cpx_env, mp_mod.cpx_lp, &this->mp_mod.curr_obj_val);
		if (checkCPXstatus(status, &mp_mod.cpx_env, &mp_mod.cpx_lp)) goto TERMINATE;


		/*Do we need this all the time? Yes, at least for the w.*/
		status = CPXgetx(mp_mod.cpx_env, mp_mod.cpx_lp, mp_mod.curr_sol.data(), 0, CPXgetnumcols(mp_mod.cpx_env, mp_mod.cpx_lp) - 1);
		if (checkCPXstatus(status, &mp_mod.cpx_env, &mp_mod.cpx_lp)) goto TERMINATE;

	}
	else
	{
		this->mp_mod.curr_obj_val = DBL_MAX;

	}


	//status = CPXgetbestobjval(cpx_env_bm, cpx_lp_bm, &this->bm_best_lb);
	//if (checkCPXstatus(status, &cpx_env_bm, &cpx_lp_bm)) goto TERMINATE;



	status = CPXgetmiprelgap(mp_mod.cpx_env, mp_mod.cpx_lp, &this->mp_mod.curr_solver_gap);
	if (checkCPXstatus(status, &mp_mod.cpx_env, &mp_mod.cpx_lp)) goto TERMINATE;


	//this->bandc_nodes = CPXgetnodecnt(cpx_env_bm, cpx_lp_bm);

	/*
	if (bandc_nodes == 0 )
	{
		this->root_node_lb = bm_best_ub;
	}
	*/

TERMINATE:
	return is_ok;
}

bool ExaDCCG::mp_obtain_linear_relaxation()
{
	/*
	1: change variables type for w.
	2: fix their bounds.
	*/

	bool mod_stat = true;
	int status;
	int cnt = 0;

	for (unsigned short n = this->first_pnode; n <= this->last_pnode; ++n)
	{

		this->matind[cnt] = mp_mod.v_w_i[n];
		matchar[cnt] = 'C';
		//this->matval[cnt] = this->mp_mod.curr_xibar_i[n];
		++cnt;
	}
	status = CPXchgctype(mp_mod.cpx_env, mp_mod.cpx_lp, cnt, matind.data(), matchar.data());
	if (checkCPXstatus(status, &sp_mod.cpx_env, &sp_mod.cpx_lp)) goto TERMINATE;


	/*Fix the w bounds*/
	cnt = 0;
	for (unsigned short n = this->first_pnode; n <= this->last_pnode; ++n)
	{

		this->matind[cnt] = mp_mod.v_w_i[n];
		matchar[cnt] = 'B';
		this->matval[cnt] = std::round(mp_mod.curr_sol[mp_mod.v_w_i[n]]);
		++cnt;
	}
	status = CPXchgbds(mp_mod.cpx_env, mp_mod.cpx_lp, cnt, matind.data(), matchar.data(), matval.data());
	if (checkCPXstatus(status, &sp_mod.cpx_env, &sp_mod.cpx_lp)) goto TERMINATE;


	status = CPXchgprobtype(mp_mod.cpx_env, mp_mod.cpx_lp, CPXPROB_LP);
	if (checkCPXstatus(status, &sp_mod.cpx_env, &sp_mod.cpx_lp)) goto TERMINATE;




TERMINATE:
	if (status)
		return false;
	else
		return true;
}

bool ExaDCCG::mp_solve_linear_relaxation()
{
	bool is_ok = true;
	int status;
	char errbuf[CPXMESSAGEBUFSIZE];
	int check_stat;

	/*status = CPXsetdefaults(this->mp_mod.cpx_env);
	if (checkCPXstatus(status, &mp_mod.cpx_env, &mp_mod.cpx_lp)) { is_ok = false; goto TERMINATE; };*/


	status = CPXwriteprob(this->mp_mod.cpx_env, this->mp_mod.cpx_lp, "peneLP.lp", "lp");
	if (checkCPXstatus(status, &mp_mod.cpx_env, &mp_mod.cpx_lp)) { is_ok = false; goto TERMINATE; };


	status = CPXlpopt(this->mp_mod.cpx_env, this->mp_mod.cpx_lp);
	if (checkCPXstatus(status, &mp_mod.cpx_env, &mp_mod.cpx_lp)) { is_ok = false; goto TERMINATE; }
	//this->bm_branch_and_cut_time = (double)(end - start) / (double)CLK_TCK;


	check_stat = CPXgetstat(mp_mod.cpx_env, mp_mod.cpx_lp);

	if (mp_mod.curr_duals.size() < CPXgetnumrows(mp_mod.cpx_env, mp_mod.cpx_lp))
	{
		mp_mod.curr_duals.resize(CPXgetnumrows(mp_mod.cpx_env, mp_mod.cpx_lp));
	}


	status = CPXgetpi(mp_mod.cpx_env, mp_mod.cpx_lp, mp_mod.curr_duals.data(), 0, CPXgetnumrows(mp_mod.cpx_env, mp_mod.cpx_lp) - 1);
	if (checkCPXstatus(status, &mp_mod.cpx_env, &mp_mod.cpx_lp)) { is_ok = false; goto TERMINATE; }

	/*now save the xibar*/
	for (int i = first_pnode; i <= last_pnode; ++i)
	{
		mp_mod.curr_xibar_i[i] = mp_mod.curr_duals[mp_mod.c_dual_xibar[i]];
	}



TERMINATE:
	return is_ok;
}

bool ExaDCCG::mp_get_xibar()
{
	bool is_ok = true;

	is_ok = this->mp_obtain_linear_relaxation();
	if (!is_ok)
	{
		goto TERMINATE;
	}

	is_ok = mp_solve_linear_relaxation();
	if (!is_ok)
	{
		goto TERMINATE;
	}

TERMINATE:
	return is_ok;
}

bool ExaDCCG::mp_solver_master_problem()
{

	bool check = true;
	int status;
	int mip_status;
	char errbuf[CPXMESSAGEBUFSIZE];
	int check_mip_stat;

	std::cout << "MASTER PROBLEM" << std::endl;

	check = this->mp_build_model();
	if (!check)
	{
		std::cout << "Problem in mp_solver_master_proble(),  mp_build_model();" << std::endl;
		goto TERMINATE;
	}

	check = this->mp_solve_model();
	if (!check)
	{
		std::cout << "Problem in mp_solver_master_proble(),  this->mp_solve_model();;" << std::endl;

		goto TERMINATE;
	}

	check = this->mp_get_xibar();
	if (!check)
	{
		std::cout << "Problem in mp_solver_master_proble(), this->mp_get_xibar();" << std::endl;

		goto TERMINATE;
	}

	/*now we free the problem*/
	this->free_cplex(&mp_mod.cpx_env, &mp_mod.cpx_lp);
	if (!check)
	{
		std::cout << "Problem in mp_solver_master_proble(),  this->free_cplex(&mp_mod.cpx_env, &mp_mod.cpx_lp);" << std::endl;

		goto TERMINATE;
	}

TERMINATE: 
	return check;
}



bool ExaDCCG::sp_build_model()
{
	bool is_ok = true;
	int status;
	int mip_status;
	CPXLONG      contextmask = 0;
	char errbuf[CPXMESSAGEBUFSIZE];
	int check_mip_stat;

	double start = 0;

	double end = 0;

	is_ok = init_cplex(&this->sp_mod.cpx_env, &this->sp_mod.cpx_lp);
	if (!is_ok)
	{
		goto TERMINATE;
	}

	is_ok = sp_set_cplex();
	if (!is_ok)
	{
		goto TERMINATE;
	}

	is_ok = sp_add_vars();
	if (!is_ok)
	{
		goto TERMINATE;
	}

	is_ok = sp_add_cnstr();
	if (!is_ok)
	{
		goto TERMINATE;
	}

	contextmask |= CPX_CALLBACKCONTEXT_CANDIDATE;
	contextmask |= CPX_CALLBACKCONTEXT_RELAXATION;
	if (contextmask != 0) {
		//	 We are done and now we register our callback function. 
		status = CPXcallbacksetfunc(sp_mod.cpx_env, sp_mod.cpx_lp, contextmask, sp_general_callback, this);
		if (status != 0) {
			fprintf(stderr, "Failed to add callback: %s\n",
				CPXgeterrorstring(sp_mod.cpx_env, status, errbuf));
			goto TERMINATE;
		}
	}


TERMINATE:
	return is_ok;
}

bool ExaDCCG::sp_set_cplex()
{
	bool mod_stat = true;

	int status = 0;

	char errbuf[CPXMESSAGEBUFSIZE];

	
	status = CPXsetintparam(this->sp_mod.cpx_env, CPX_PARAM_SCRIND, CPX_ON);
	if (checkCPXstatus(status, &sp_mod.cpx_env, &sp_mod.cpx_lp)) return false;
	status = CPXsetintparam(this->sp_mod.cpx_env, CPX_PARAM_THREADS, 1);
	if (checkCPXstatus(status, &sp_mod.cpx_env, &sp_mod.cpx_lp)) return false;
	status = CPXchgobjsen(this->sp_mod.cpx_env, this->sp_mod.cpx_lp, CPX_MIN);
	if (checkCPXstatus(status, &sp_mod.cpx_env, &sp_mod.cpx_lp)) return false;
	status = CPXsetintparam(this->sp_mod.cpx_env, CPX_PARAM_NUMERICALEMPHASIS, CPX_ON);
	if (checkCPXstatus(status, &sp_mod.cpx_env, &sp_mod.cpx_lp)) return false;
	status = CPXsetdblparam(this->sp_mod.cpx_env, CPXPARAM_Simplex_Tolerances_Feasibility, 1e-9);
	if (checkCPXstatus(status, &sp_mod.cpx_env, &sp_mod.cpx_lp)) return false;
	//status = CPXsetintparam(cpx_env, CPX_PARAM_PREIND, CPX_OFF);
	//if (checkCPXstatus(status))return false;

	/*To be extra numerical accurate*/
	status = CPXsetdblparam(this->sp_mod.cpx_env, CPXPARAM_Simplex_Tolerances_Optimality, 1e-9);
	if (checkCPXstatus(status, &sp_mod.cpx_env, &sp_mod.cpx_lp)) return false;

	status = CPXsetdblparam(this->sp_mod.cpx_env, CPXPARAM_Simplex_Tolerances_Feasibility, 1e-9);
	if (checkCPXstatus(status, &sp_mod.cpx_env, &sp_mod.cpx_lp)) return false;

	status = CPXsetdblparam(this->sp_mod.cpx_env, CPXPARAM_MIP_Tolerances_MIPGap, 1e-9);
	if (checkCPXstatus(status, &sp_mod.cpx_env, &sp_mod.cpx_lp)) return false;
	status = CPXsetdblparam(this->sp_mod.cpx_env, CPXPARAM_MIP_Tolerances_Integrality, 0);
	if (checkCPXstatus(status, &sp_mod.cpx_env, &sp_mod.cpx_lp)) return false;

	double time_limit = this->cfg->getValueOfKey<double>("TIME_LIMIT");

	status = CPXsetdblparam(sp_mod.cpx_env, CPXPARAM_TimeLimit, time_limit);
	if (checkCPXstatus(status, &sp_mod.cpx_env, &sp_mod.cpx_lp)) return false;

	/*
	status = CPXsetintparam(cpx_env, CPXPARAM_MIP_Strategy_VariableSelect, CPX_VARSEL_STRONG);
	if (checkCPXstatus(status)) return false;
	*/


	return true;
}

bool ExaDCCG::sp_add_vars()
{
	unsigned short num_cols = CPXgetnumcols(this->sp_mod.cpx_env, this->sp_mod.cpx_lp);

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

	/*Discovery variables*/
	for (int i = first_pnode; i <= last_pnode; ++i)
	{


		this->sp_mod.v_w_i[i] = num_cols;
		obj[0] = this->disc_cost_n[i];

		vartype[0] = 'B';

		sprintf(varsname[0], "w_i%d", i);
		lb[0] = 0;
		ub[0] = 1;

		status = CPXnewcols(this->sp_mod.cpx_env, this->sp_mod.cpx_lp, 1, obj, lb, ub, vartype, varsname);
		if (checkCPXstatus(status, &sp_mod.cpx_env, &sp_mod.cpx_lp)) goto TERMINATE;
		++num_cols;
	}



	/*
	 Routing vatriables, y_kn and x_kij.
	*/
	/*y_kn*/

	for (int i = first_pnode; i <= last_pnode; ++i)
	{


		this->sp_mod.v_y_i[i] = num_cols;

		obj[0] = this->pp[i];

		vartype[0] = 'B';

		sprintf(varsname[0], "y_i%d", i);
		lb[0] = 0;
		ub[0] = 1;




		status = CPXnewcols(this->sp_mod.cpx_env, this->sp_mod.cpx_lp, 1, obj, lb, ub, vartype, varsname);
		if (checkCPXstatus(status, &sp_mod.cpx_env, &sp_mod.cpx_lp)) goto TERMINATE;
		++num_cols;
	}

	/*x_ij*/

	for (int i = start_node; i <= end_node; ++i)
		for (int j = start_node; j <= end_node; ++j)
			if (this->edge_exists(i, j))
			{
				single_node_dur_1 = 0;
				single_node_dur_2 = 0;

				this->sp_mod.v_x_ij(i, j) = num_cols;

				obj[0] = 0;

				vartype[0] = 'B';

				sprintf(varsname[0], "x_i%d_j%d", i, j);
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

				status = CPXnewcols(this->sp_mod.cpx_env, this->sp_mod.cpx_lp, 1, obj, lb, ub, vartype, varsname);
				if (checkCPXstatus(status, &sp_mod.cpx_env, &sp_mod.cpx_lp)) goto TERMINATE;
				++num_cols;
			}

	/*Dual variables */
	for (unsigned short r = 0; r < this->num_us_constr; ++r)
	{


		this->sp_mod.v_beta[r] = num_cols;

		//v_beta[r] = num_cols;

		obj[0] = this->b[r];

		vartype[0] = 'C';

		sprintf(varsname[0], "beta_r%d", r);

		lb[0] = 0;
		ub[0] = CPX_INFBOUND;

		status = CPXnewcols(this->sp_mod.cpx_env, this->sp_mod.cpx_lp, 1, obj, lb, ub, vartype, varsname);
		if (checkCPXstatus(status, &sp_mod.cpx_env, &sp_mod.cpx_lp)) goto TERMINATE;		++num_cols;
	}


	for (unsigned short n = this->first_pnode; n <= this->last_pnode; ++n)
	{
		this->sp_mod.v_gamma[n] = num_cols;


		if (/*this->w_val[n]
			&&*/
			fabs(this->mp_mod.curr_xibar_i[n]) > 1e-9
			)
		{

			obj[0] = this->mp_mod.curr_xibar_i[n];
		}
		else
			obj[0] = 0;

		vartype[0] = 'C';

		sprintf(varsname[0], "gamma_n%d", n);

		
		lb[0] = -CPX_INFBOUND;
		ub[0] = CPX_INFBOUND;
		
		

		status = CPXnewcols(this->sp_mod.cpx_env, this->sp_mod.cpx_lp, 1, obj, lb, ub, vartype, varsname);
		if (checkCPXstatus(status, &sp_mod.cpx_env, &sp_mod.cpx_lp)) goto TERMINATE;		++num_cols;
	}


TERMINATE:

	if (status)
	{
		return false;
	}

	return true;
}

bool ExaDCCG::sp_add_cnstr()
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


	// Constraints on the discovery budget

	nzc = 0;
	rhs[0] = 
		this->max_num_disc;
	sense[0] = 'L';
	sprintf(cnstrname[0], "Discovery_budget");

	for (int i = first_pnode; i <= last_pnode; ++i)
	{

		matind[nzc] = this->sp_mod.v_w_i[i];
		matval[nzc] = 1;
		++nzc;

	}
	status = CPXaddrows(sp_mod.cpx_env, sp_mod.cpx_lp, 0, 1, nzc, rhs, sense, matbeg, matind.data(), matval.data(), NULL, cnstrname);
	if (checkCPXstatus(status, &sp_mod.cpx_env, &sp_mod.cpx_lp)) goto TERMINATE;

	// Degree constraints on profitable nodes
	nzc = 0;
	rhs[0] = 0;
	sense[0] = 'E';

	for (int i = this->first_pnode; i <= this->last_pnode; ++i)
	{
		sprintf(cnstrname[0], "Degrees_i%d", i);
		nzc = 0;
		for (int ii = i + 1; ii <= end_node; ++ii)
			if (this->edge_exists(i, ii))
			{
				matind[nzc] = this->sp_mod.v_x_ij(i, ii);
				matval[nzc] = 1;
				++nzc;
			}
		for (int ii = start_node; ii < i; ++ii)
			if (this->edge_exists(ii, i))
			{
				matind[nzc] = this->sp_mod.v_x_ij(ii, i);
				matval[nzc] = 1;
				++nzc;
			}


		matind[nzc] = sp_mod.v_y_i[i];
		matval[nzc] = -2;
		++nzc;

		status = CPXaddrows(sp_mod.cpx_env, sp_mod.cpx_lp, 0, 1, nzc, rhs, sense, matbeg, matind.data(), matval.data(), NULL, cnstrname);
		if (checkCPXstatus(status, &sp_mod.cpx_env, &sp_mod.cpx_lp)) goto TERMINATE;
	}

	// Degree constraint  on start node 
	rhs[0] = 1;
	sense[0] = 'E';
	nzc = 0;

	sprintf(cnstrname[0], "DegreeDepotOut");
	for (int j = this->start_node; j <= this->end_node; ++j)
		if (this->edge_exists(start_node, j))
		{
			matind[nzc] = sp_mod.v_x_ij(this->start_node, j);
			matval[nzc] = 1;
			++nzc;

		}
	status = CPXaddrows(sp_mod.cpx_env, sp_mod.cpx_lp, 0, 1, nzc, rhs, sense, matbeg, matind.data(), matval.data(), NULL, cnstrname);
	if (checkCPXstatus(status, &sp_mod.cpx_env, &sp_mod.cpx_lp)) goto TERMINATE;


	// Degree constraint  on end node 
	rhs[0] = 1;
	sense[0] = 'E';
	nzc = 0;
	sprintf(cnstrname[0], "DegreeDepotIn");

	for (int j = start_node; j <= end_node; ++j)
		if (this->edge_exists(j, end_node))
		{
			matind[nzc] = sp_mod.v_x_ij(j, end_node);
			matval[nzc] = 1;
			++nzc;

		}
	status = CPXaddrows(sp_mod.cpx_env, sp_mod.cpx_lp, 0, 1, nzc, rhs, sense, matbeg, matind.data(), matval.data(), NULL, cnstrname);
	if (checkCPXstatus(status, &sp_mod.cpx_env, &sp_mod.cpx_lp)) goto TERMINATE;


	nzc = 0;
	rhs[0] = max_dur;
	sense[0] = 'L';
	sprintf(cnstrname[0], "Max_Dur");
	for (int i = start_node; i < end_node; ++i)
		for (int j = start_node + 1; j <= end_node; ++j)
			if (this->edge_exists(i, j))
			{
			
				matind[nzc] = sp_mod.v_x_ij(i, j);
				matval[nzc] = t_ij(i, j);
				++nzc;
			}
	status = CPXaddrows(sp_mod.cpx_env, sp_mod.cpx_lp, 0, 1, nzc, rhs, sense, matbeg, matind.data(), matval.data(), NULL, cnstrname);
	if (checkCPXstatus(status, &sp_mod.cpx_env, &sp_mod.cpx_lp)) goto TERMINATE;


	/*Dual constraints*/
	nzc = 0;
	rhs[0] = 0;
	sense[0] = 'E';

	for (int i = first_pnode; i <= last_pnode; ++i)
	{
		sprintf(cnstrname[0], "Balance_i%d", i);

		nzc = 0;
		rhs[0] = 0;
		sense[0] = 'E';

		matind[nzc] = sp_mod.v_y_i[i];
		matval[nzc] = -ppxi[i];
		++nzc;


		matind[nzc] = sp_mod.v_gamma[i];
		matval[nzc] = 1; // Becase we cnage the bound of gamma when w = 0.

		++nzc;

		for (int l = 0; l < this->num_us_constr; ++l)
			if (fabs(A(l, i)) > eps_coeff)
			{
				matind[nzc] = sp_mod.v_beta[l];
				matval[nzc] = A(l, i);
				++nzc;
			}

		status = CPXaddrows(sp_mod.cpx_env, sp_mod.cpx_lp, 0, 1, nzc, rhs, sense, matbeg, matind.data(), matval.data(), NULL, cnstrname);
		if (checkCPXstatus(status, &sp_mod.cpx_env, &sp_mod.cpx_lp)) goto TERMINATE;
	}

	/*Big-m gamma w*/
			/*Product w-sigma constraints*/
	for (int i = first_pnode; i <= last_pnode; ++i)
	{
		nzc = 0;
		rhs[0] = 0;
		sense[0] = 'G';
		sprintf(cnstrname[0], "Linearize_gammaW_GQ_%d", i);

		matind[nzc] = this->sp_mod.v_gamma[i];
		matval[nzc] = 1;
		++nzc;

		matind[nzc] = this->sp_mod.v_w_i[i];
		matval[nzc] = BIG_M_SIGMA;
		++nzc;


		status = CPXaddrows(sp_mod.cpx_env, sp_mod.cpx_lp, 0, 1, nzc, rhs, sense, matbeg, matind.data(), matval.data(), NULL, cnstrname);
		if (checkCPXstatus(status, &sp_mod.cpx_env, &sp_mod.cpx_lp)) goto TERMINATE;

		nzc = 0;
		rhs[0] = 0;
		sense[0] = 'L';
		sprintf(cnstrname[0], "Linearize_gammaW_LQ_%d", i);

		matind[nzc] = this->sp_mod.v_gamma[i];
		matval[nzc] = 1;
		++nzc;

		matind[nzc] = this->sp_mod.v_w_i[i];
		matval[nzc] = -BIG_M_SIGMA;
		++nzc;
		status = CPXaddrows(sp_mod.cpx_env, sp_mod.cpx_lp, 0, 1, nzc, rhs, sense, matbeg, matind.data(), matval.data(), NULL, cnstrname);
		if (checkCPXstatus(status, &sp_mod.cpx_env, &sp_mod.cpx_lp)) goto TERMINATE;
	}


TERMINATE:
	if (status)
		return false;
	else
		return true;
}

bool ExaDCCG::sp_update_model()
{
	bool mod_stat = true;

	int status;
	int cnt = 0;
	double rounded_to_precision = 0;

	for (unsigned short n = this->first_pnode; n <= this->last_pnode; ++n)
	{

		this->matind[cnt] = sp_mod.v_gamma[n];
			//rounded_to_precision = std::round(this->mlp_curr_sol[this->v_pr_xi[n]] * 1e9) / 1e9;
		this->matval[cnt] = this->mp_mod.curr_xibar_i[n];
		++cnt;
	}

	status = CPXchgobj(sp_mod.cpx_env, sp_mod.cpx_lp, cnt, matind.data(), matval.data());
	if (checkCPXstatus(status, &sp_mod.cpx_env, &sp_mod.cpx_lp)) goto TERMINATE;



	/*try fixing w*/
	//cnt = 0;
	//for (unsigned short n = this->first_pnode; n <= this->last_pnode; ++n)
	//{

	//	this->matind[cnt] = sp_mod.v_w_i[n];
	//		//mp_mod.v_w_i[n];
	//	matchar[cnt] = 'B';
	//	this->matval[cnt] = std::round(mp_mod.curr_sol[mp_mod.v_w_i[n]]);
	//	++cnt;
	//}
	//status = CPXchgbds(sp_mod.cpx_env, sp_mod.cpx_lp, cnt, matind.data(), matchar.data(), matval.data());
	//if (checkCPXstatus(status, &sp_mod.cpx_env, &sp_mod.cpx_lp)) goto TERMINATE;


TERMINATE:
	if (status)
		return false;
	else
		return true;

}

bool ExaDCCG::sp_solve_model()
{
	bool is_ok = true;
	int status;
	int mip_status;
	CPXLONG      contextmask = 0;
	char errbuf[CPXMESSAGEBUFSIZE];
	int check_mip_stat;

		status = CPXwriteprob(this->sp_mod.cpx_env, this->sp_mod.cpx_lp, "slave.lp", "lp");
		if (checkCPXstatus(status, &sp_mod.cpx_env, &sp_mod.cpx_lp)) { is_ok = false; goto TERMINATE; };


	status = CPXmipopt(this->sp_mod.cpx_env, this->sp_mod.cpx_lp);
	if (checkCPXstatus(status, &sp_mod.cpx_env, &sp_mod.cpx_lp)) goto TERMINATE;



	check_mip_stat = CPXgetstat(sp_mod.cpx_env, sp_mod.cpx_lp);

	status = CPXgetobjval(sp_mod.cpx_env, sp_mod.cpx_lp, &this->sp_mod.curr_obj_val);
	if (checkCPXstatus(status, &sp_mod.cpx_env, &sp_mod.cpx_lp)) goto TERMINATE;


	if (this->mp_mod.curr_obj_val > this->sp_mod.curr_obj_val + EXA_EPS)
	{
		status = CPXgetx(sp_mod.cpx_env, sp_mod.cpx_lp, sp_mod.curr_sol.data(), 0, CPXgetnumcols(sp_mod.cpx_env, sp_mod.cpx_lp) - 1);
		if (checkCPXstatus(status, &sp_mod.cpx_env, &sp_mod.cpx_lp)) goto TERMINATE;


		if (mp_mod.Y.getNumRows() <= mp_mod.card_Y)
		{
			mp_mod.Y.resizeMatrix2D(mp_mod.card_Y * 2, num_nodes);
			mp_mod.v_alpha_y.resize(mp_mod.card_Y * 2);
			mp_mod.v_beta_y_l.resizeMatrix2D(mp_mod.card_Y * 2, this->num_us_constr);
			mp_mod.v_gamma_y_i.resizeMatrix2D(mp_mod.card_Y * 2, num_nodes);

		}

		for (int i = first_pnode; i <= last_pnode; ++i)
		{
			mp_mod.Y(mp_mod.card_Y, i) = round(sp_mod.curr_sol[this->sp_mod.v_y_i[i]]);
		}
		++mp_mod.card_Y;
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


int CPXPUBLIC ExaDCCG::sp_general_callback(CPXCALLBACKCONTEXTptr context, CPXLONG contextid, void* userdata)
{
	bool cut_added = false;
	int status = 0;
	bool is_ok = true;
	ExaDCCG* solver = (ExaDCCG*)userdata;
	double val;
	double treshold;


	if (contextid == CPX_CALLBACKCONTEXT_CANDIDATE)
	{

		status = CPXcallbackgetcandidatepoint(context, solver->sp_mod.curr_sol.data(), 0, CPXgetnumcols(solver->sp_mod.cpx_env, solver->sp_mod.cpx_lp) - 1, &val);
		if (status != 0)
		{
			solver->checkCPXstatus(status, &solver->sp_mod.cpx_env, &solver->sp_mod.cpx_lp);

			return status;
		}

		if (!cut_added)
			is_ok = solver->separate_heuristic_subtours(context, contextid, 0, &cut_added);


	}

	if (contextid == CPX_CALLBACKCONTEXT_RELAXATION)
	{
		/*Only if  CPX >= 12.10*/
		CPXLONG node_id, node_depth;




		CPXcallbackgetinfolong(context, CPXCALLBACKINFO_NODEUID, &node_id);
		if (status != 0)
		{
			solver->checkCPXstatus(status, &solver->sp_mod.cpx_env, &solver->sp_mod.cpx_lp);			return status;
		}

		CPXcallbackgetinfolong(context, CPXCALLBACKINFO_NODEDEPTH, &node_depth);
		if (status != 0)
		{
			solver->checkCPXstatus(status, &solver->sp_mod.cpx_env, &solver->sp_mod.cpx_lp);			return status;
		}





		status = CPXcallbackgetrelaxationpoint(context, solver->sp_mod.curr_sol.data(), 0, CPXgetnumcols(solver->sp_mod.cpx_env, solver->sp_mod.cpx_lp) - 1, &val);
		if (status != 0)
		{
			solver->checkCPXstatus(status, &solver->sp_mod.cpx_env, &solver->sp_mod.cpx_lp);				return status;
		}






		is_ok = solver->separate_routing_logic_cuts(context, contextid, &cut_added);

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



				treshold += TRESHOLD_STEP;

				if (treshold > TRESHOLD_MAX)
					break;
			}
		}

		if (node_depth <= 1)
		{
			is_ok = solver->separate_exact_subtuors(context, contextid, &cut_added);
		}


	}

	return 0;
}

bool ExaDCCG::separate_heuristic_subtours(CPXCALLBACKCONTEXTptr context, CPXLONG contextid, double var_treshold, bool* cut_added)
{
	int status = 0;
	double val;

	
	build_sol_directed_graph(sp_mod.curr_sol.data(), var_treshold);
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
				nzc = 0;
				rhs[0] = 0;
				lhs_v = 0;
				yi_max_val = 0;
				yi_max = 1;

				for (int i = 1; i <= node_subtour_subset[0]; ++i)
				{
					if (this->sp_mod.curr_sol[sp_mod.v_y_i[node_subtour_subset[i]]] > yi_max_val)
					{
						yi_max = i;
						yi_max_val = this->sp_mod.curr_sol[sp_mod.v_y_i[node_subtour_subset[i]]];
					}
				}


				for (int i = 1; i <= node_subtour_subset[0]; ++i)
				{
					if (i != yi_max)
					{
						
						matind[nzc] = sp_mod.v_y_i[node_subtour_subset[i]];
						matval[nzc] = -1;

						++nzc;

						lhs_v -= sp_mod.curr_sol[sp_mod.v_y_i[node_subtour_subset[i]]];

					}
				}
				for (int i = 1; i < node_subtour_subset[0]; ++i)
					for (int j = i + 1; j <= node_subtour_subset[0]; ++j)
					{
						matind[nzc] = 
							sp_mod.v_x_ij(node_subtour_subset[i], node_subtour_subset[j]);
						matval[nzc] = 1;
						++nzc;


						lhs_v += this->sp_mod.curr_sol[sp_mod.v_x_ij(node_subtour_subset[i], node_subtour_subset[j])];
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
						if (status)checkCPXstatus(status, &sp_mod.cpx_env, &sp_mod.cpx_lp);
					}
					else
						if (contextid == CPX_CALLBACKCONTEXT_RELAXATION)
						{
							//std::cout << "Cut" << std::endl;
							status = CPXcallbackaddusercuts(context, 1, nzc, rhs,
								sense, matbeg, matind.data(), matval.data(), force, local);
							if (status)checkCPXstatus(status, &sp_mod.cpx_env, &sp_mod.cpx_lp);
						}
				}
			}
	}
	return true;
}

bool ExaDCCG::separate_exact_subtuors(CPXCALLBACKCONTEXTptr context, CPXLONG contextid, bool* cut_added)
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


	double lhs_val = 0;


	build_sol_directed_graph(sp_mod.curr_sol.data(),0.0);
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

				nzc = 0;
				lhs_val = 0;
				//	node_subtour_subset[0] = 0;
				yi_max_val = 0;
				yi_max = 1;


				for (int i = 1; i <= node_subtour_subset[0]; ++i)
				{
					if (this->sp_mod.curr_sol[sp_mod.v_y_i[node_subtour_subset[i]]] > yi_max_val)
					{
						yi_max = i;
						yi_max_val = this->sp_mod.curr_sol[sp_mod.v_y_i[node_subtour_subset[i]]];
					}
				}


				for (int i = 1; i <= node_subtour_subset[0]; ++i)
				{
					if (i != yi_max)
					{

						matind[nzc] = sp_mod.v_y_i[node_subtour_subset[i]];
						matval[nzc] = -1;



						lhs_val -= sp_mod.curr_sol[sp_mod.v_y_i[node_subtour_subset[i]]];
						++nzc;


					}
				}

				for (int i = 1; i < node_subtour_subset[0]; ++i)
					for (int j = i + 1; j <= node_subtour_subset[0]; ++j)
					{
						matind[nzc] = sp_mod.v_x_ij(node_subtour_subset[i], node_subtour_subset[j]);
						matval[nzc] = 1;
						lhs_val += sp_mod.curr_sol[sp_mod.v_x_ij(node_subtour_subset[i], node_subtour_subset[j])];

						++nzc;

					}

				if (lhs_val > rhs[0] + cEPS)
				{
					status = CPXcallbackaddusercuts(context, 1, nzc, rhs,
						sense, matbeg, matind.data(), matval.data(), force, local);
					if (status)checkCPXstatus(status,&sp_mod.cpx_env,&sp_mod.cpx_lp);
				}
				if (status != 0)
				{
					return false;
				}
			}
		}
	return true;
}

bool ExaDCCG::separate_routing_logic_cuts(CPXCALLBACKCONTEXTptr context, CPXLONG contextid, bool* cut_added)
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

	// Logic VI
	for (int i = start_node; i < end_node; ++i)
		for (int j = i + 1; j <= end_node; ++j)
		{

			/* Separate for integer and add also for tile if violeted */
			if (i != start_node
				&&
				sp_mod.curr_sol[sp_mod.v_x_ij(i, j)] > sp_mod.curr_sol[sp_mod.v_y_i[i]] + cEPS)
			{
				*cut_added = true;

				// Add the cuts for i on y and x. 

				rhs[0] = 0;
				sense[0] = 'L';
				nzc = 0;
				lhs_v = 0;

				matind[nzc] = sp_mod.v_x_ij(i, j);
				matval[nzc] = 1;
				lhs_v += sp_mod.curr_sol[sp_mod.v_x_ij(i, j)];
				++nzc;

				matind[nzc] = sp_mod.v_y_i[i];
				matval[nzc] = -1;
				lhs_v -= sp_mod.curr_sol[sp_mod.v_y_i[i]];
				++nzc;

				if (lhs_v > rhs[0] + cEPS)
				{
					status = CPXcallbackaddusercuts(context, 1, nzc, rhs,
						sense, matbeg, matind.data(), matval.data(), force, local);

					if (status)
					{
						if (status)checkCPXstatus(status, &sp_mod.cpx_env, &sp_mod.cpx_lp);

						return false;
					}

				}



			}

			if (j != end_node
				&&
				sp_mod.curr_sol[sp_mod.v_x_ij(i, j)] > sp_mod.curr_sol[sp_mod.v_y_i[j]] + cEPS)
			{
				*cut_added = true;


				// Add the cuts for j, x and y and for each k
				rhs[0] = 0;
				sense[0] = 'L';
				nzc = 0;
				lhs_v = 0;

				matind[nzc] = sp_mod.v_x_ij(i, j);
				matval[nzc] = 1;
				++nzc;
				lhs_v += sp_mod.curr_sol[sp_mod.v_x_ij(i, j)];

				matind[nzc] = sp_mod.v_y_i[j];
				matval[nzc] = -1;
				lhs_v -= sp_mod.curr_sol[sp_mod.v_y_i[j]];
				++nzc;

				if (lhs_v > rhs[0] + mEPS)
				{
					status = CPXcallbackaddusercuts(context, 1, nzc, rhs,
						sense, matbeg, matind.data(), matval.data(), force, local);
					if (status)
					{
						if (status)checkCPXstatus(status, &sp_mod.cpx_env, &sp_mod.cpx_lp);
						//fprintf(stdout, "Cut Logic!\n")
						return false;
					}
				}

			}
		}

	return true;
}

bool ExaDCCG::init_CCG_algorithm()
{
	bool check = true;
	int status;
	int mip_status;
	char errbuf[CPXMESSAGEBUFSIZE];
	int check_mip_stat;


	double start = 0;
	double end = 0;

	/*
	1) Reset data structure
	2) populate Y
	3) Init and solve master problem
	4) Get xi bar
	5) init and update slave
	6) Go back to 3 an repreat until convergence
	*/

	check = this->generate_initial_y();
	if (!check)
	{
		goto TERMINATE;
	}

	check =  this->mp_solver_master_problem();
	if (!check)
	{
		goto TERMINATE;
	}

	check = this->sp_build_model(); 
	if (!check)
	{
		goto TERMINATE;
	}

	check = this->sp_solve_model();
	if (!check)
	{
		goto TERMINATE;
	}

TERMINATE:
	return check;
}

bool ExaDCCG::generate_initial_y()
{
	mp_mod.card_Y = 0;

	/*minimum extension all the time*/
	unsigned short i;
	unsigned int T = 0;
	i = 0;

	for (unsigned short j = first_pnode; j <= last_pnode; ++j)
		if(
			//true
			T + this->t_ij(i,j) + this->t_ij(j,0) <= this->max_dur
			)
	{
			mp_mod.Y(mp_mod.card_Y, j) = 1;
			T = T + this->t_ij(i, j);
			i = j;
	}
		else
		{
			break;
		}
	++i;
	for (; i <= last_pnode; ++i)
	{
		mp_mod.Y(mp_mod.card_Y, i) = 0;

	}

	++mp_mod.card_Y;

	return true;
}
