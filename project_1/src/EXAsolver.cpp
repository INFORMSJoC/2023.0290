#include "EXAsolver.h"

inline bool EXAsolver::init_queue()
{
	sep_queue.num_el = 0;
	sep_queue.next = 0;
	return true;
}

inline bool EXAsolver::queue_is_empty()
{
	if (sep_queue.num_el == 0
		||
		sep_queue.num_el == sep_queue.next)
		return true;
	else
		return false;
}

inline bool EXAsolver::push_back(int el)
{
	sep_queue.elem[sep_queue.num_el] = el;
	++sep_queue.num_el;
	return true;
}

inline int EXAsolver::pop_front()
{
	int elem = sep_queue.elem[sep_queue.next];
	++sep_queue.next;
	return elem;
}




int EXAsolver::rop_build_sol_graph(const double* curr_sol, double val_treshold = 0)
{
	sep_un_graph.num_nodes = 0;
	sep_un_graph.num_edges = 0;

	sep_un_graph.num_nodes = 1;
	sep_un_graph.nodes_map[sep_un_graph.num_nodes - 1] = start_node;
	//g_start_node;

	for (int i = start_node + 1; i < end_node; ++i)
		if (curr_sol[rop_v_y_n[i]] > mEPS)
		{

			sep_un_graph.nodes_map[sep_un_graph.num_nodes] = i;
			++sep_un_graph.num_nodes;
		}

	sep_un_graph.nodes_map[sep_un_graph.num_nodes] = end_node;
	++sep_un_graph.num_nodes;

	for (int i = 0; i < sep_un_graph.num_nodes - 1; ++i)
		for (int j = i + 1; j < sep_un_graph.num_nodes; ++j)
			if (curr_sol[v_x_ij(sep_un_graph.nodes_map[i], sep_un_graph.nodes_map[j])] > mEPS + val_treshold)
				//if (curr_sol[g_v_xk_kij[k][sep_un_graph.nodes_map[i]][sep_un_graph.nodes_map[j]]] > mEPS)
			{

				// Add edge.
				sep_un_graph.edge[sep_un_graph.num_edges].first = i;
				sep_un_graph.edge[sep_un_graph.num_edges].second = j;
				++sep_un_graph.num_edges;
			}
	return 0;
}




unsigned int EXAsolver::connected_components()
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

inline void EXAsolver::connected_components_recursive(unsigned int vertex, unsigned int component)
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

int EXAsolver::bm_build_sol_graph(const double* curr_sol, double val_treshold, int xi_bar)
{
	sep_un_graph.num_nodes = 0;
	sep_un_graph.num_edges = 0;

	sep_un_graph.num_nodes = 1;
	sep_un_graph.nodes_map[sep_un_graph.num_nodes - 1] = start_node;
	//g_start_node;

	for (int i = start_node + 1; i < end_node; ++i)
		if(curr_sol[bm_v_y_xib_i(xi_bar, i)] > mEPS)
		//if (curr_sol[rop_v_y_n[i]] > mEPS)
		{

			sep_un_graph.nodes_map[sep_un_graph.num_nodes] = i;
			++sep_un_graph.num_nodes;
		}

	sep_un_graph.nodes_map[sep_un_graph.num_nodes] = end_node;
	++sep_un_graph.num_nodes;

	for (int i = 0; i < sep_un_graph.num_nodes - 1; ++i)
		for (int j = i + 1; j < sep_un_graph.num_nodes; ++j)
			if (curr_sol[bm_v_x_xib_ij(xi_bar,sep_un_graph.nodes_map[i], sep_un_graph.nodes_map[j])] > mEPS + val_treshold)
			//if (curr_sol[v_x_ij(sep_un_graph.nodes_map[i], sep_un_graph.nodes_map[j])] > mEPS + val_treshold)
			{

				// Add edge.
				sep_un_graph.edge[sep_un_graph.num_edges].first = i;
				sep_un_graph.edge[sep_un_graph.num_edges].second = j;
				++sep_un_graph.num_edges;
			}
	return 0;
}


inline bool EXAsolver::build_sol_directed_graph(const double* curr_sol)
{


	sep_di_graph.num_ver = 0;



	sep_di_graph.num_ver = 1;
	sep_di_graph.ver_label[sep_di_graph.num_ver - 1] = start_node;




	for (int i = start_node + 1; i < end_node; ++i)
		if (curr_sol[rop_v_y_n[i]] > mEPS)
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
			if (curr_sol[v_x_ij(sep_di_graph.ver_label[i], sep_di_graph.ver_label[j])] > mEPS)
				//	if (curr_sol[g_cg_price.v_x_ij[sep_di_graph.ver_label[i]][sep_di_graph.ver_label[j]]] > mmEPS)
			{


				double dq =
					(curr_sol[
						v_x_ij(sep_di_graph.ver_label[i]
							, sep_di_graph.ver_label[j])
					])
					* scaleMXflow;
						//(curr_sol[g_cg_price.v_x_ij[sep_di_graph.ver_label[i]][sep_di_graph.ver_label[j]]]) * scaleMXflow;


						int q = (int)dq;
						add_archs(i, j, q);
			}

	// Now the depot (the capacity is given by the sum)
	for (int i = 1; i < sep_di_graph.num_ver; ++i)
		if (curr_sol[v_x_ij(start_node, sep_di_graph.ver_label[i])] > mEPS
			||
			curr_sol[v_x_ij(sep_di_graph.ver_label[i], end_node)] > mEPS)
			/*	if (curr_sol[g_cg_price.v_x_ij[g_start_node][sep_di_graph.ver_label[i]]] > mEPS
					||
					curr_sol[g_cg_price.v_x_ij[sep_di_graph.ver_label[i]][g_end_node]] > mEPS)*/
		{
			double dq =
				(curr_sol[v_x_ij(start_node, sep_di_graph.ver_label[i])] +
					curr_sol[v_x_ij(sep_di_graph.ver_label[i], end_node)]) * scaleMXflow;
			/*(curr_sol[v_x_kij[g_start_node][sep_di_graph.ver_label[i]]] +
				curr_sol[g_cg_price.v_x_ij[sep_di_graph.ver_label[i]][g_end_node]]) * scaleMXflow;*/
			int q = (int)dq;

			add_archs(0, i, q);
		}


	return 0;
}

bool EXAsolver::bm_build_sol_directed_graph(const double* curr_sol, int xi)
{

	sep_di_graph.num_ver = 0;



	sep_di_graph.num_ver = 1;
	sep_di_graph.ver_label[sep_di_graph.num_ver - 1] = start_node;




	for (int i = start_node + 1; i < end_node; ++i)
		//if (curr_sol[rop_v_y_n[i]] > mEPS)
		if (curr_sol[bm_v_y_xib_i(xi,i)] > mEPS)
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
			if (curr_sol[
				bm_v_x_xib_ij(xi, sep_di_graph.ver_label[i], sep_di_graph.ver_label[j])] > mEPS)
			//if (curr_sol[v_x_ij(sep_di_graph.ver_label[i], sep_di_graph.ver_label[j])] > mEPS)
			{


				double dq =
					(curr_sol[
						bm_v_x_xib_ij(xi, sep_di_graph.ver_label[i], sep_di_graph.ver_label[j])
					])
					* scaleMXflow;
						//(curr_sol[g_cg_price.v_x_ij[sep_di_graph.ver_label[i]][sep_di_graph.ver_label[j]]]) * scaleMXflow;


						int q = (int)dq;
						add_archs(i, j, q);
			}

	// Now the depot (the capacity is given by the sum)
	for (int i = 1; i < sep_di_graph.num_ver; ++i)
		if (curr_sol[bm_v_x_xib_ij(xi,start_node, sep_di_graph.ver_label[i])] > mEPS
			||
			curr_sol[bm_v_x_xib_ij(xi,sep_di_graph.ver_label[i], end_node)] > mEPS)
			/*	if (curr_sol[g_cg_price.v_x_ij[g_start_node][sep_di_graph.ver_label[i]]] > mEPS
					||
					curr_sol[g_cg_price.v_x_ij[sep_di_graph.ver_label[i]][g_end_node]] > mEPS)*/
		{
			double dq =
				(curr_sol[bm_v_x_xib_ij(xi,start_node, sep_di_graph.ver_label[i])] +
					curr_sol[bm_v_x_xib_ij(xi,sep_di_graph.ver_label[i], end_node)]) * scaleMXflow;
			/*(curr_sol[v_x_kij[g_start_node][sep_di_graph.ver_label[i]]] +
				curr_sol[g_cg_price.v_x_ij[sep_di_graph.ver_label[i]][g_end_node]]) * scaleMXflow;*/
			int q = (int)dq;

			add_archs(0, i, q);
		}


	return 0;
}


inline void EXAsolver::add_archs(int a1, int a2, int q)
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

inline bool EXAsolver::bfs_di_graph(int s, int t)
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

inline int EXAsolver::send_flow_di_graph(int s, int flow, int t, short* start)
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

inline int EXAsolver::dinic_max_flow(int s, int t)
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


bool EXAsolver::init_cplex(CPXENVptr* cpx_env, CPXLPptr* cpx_lp)
{
	int status = 0;

	// Open env.
	*cpx_env = CPXopenCPLEX(&status);
	if (this->checkCPXstatus(status, cpx_env, cpx_lp)) return false;
	*cpx_lp = CPXcreateprob(*cpx_env, &status, inst_name);
	if (this->checkCPXstatus(status, cpx_env, cpx_lp)) return false;

	return true;
}

int EXAsolver::checkCPXstatus(int status, CPXENVptr* cpx_env, CPXLPptr* cpx_lp)
{
	char errmsg[CPXMESSAGEBUFSIZE];
	if (status == 0) return 0;

	CPXgeterrorstring(*cpx_env, status, errmsg);
	printf(" %s \n", errmsg);
	return status;
}

bool EXAsolver::free_cplex(CPXENVptr* cpx_env, CPXLPptr* cpx_lp)
{
	int status = 0;
	status = CPXfreeprob(*cpx_env, cpx_lp);
	if (checkCPXstatus(status, cpx_env, cpx_lp)) return false;
	status = CPXcloseCPLEX(cpx_env);
	if (checkCPXstatus(status, cpx_env, cpx_lp)) return false;

	return true;
}

bool EXAsolver::init_out_file()
{
	if (out_file != NULL)
		return true;


	if (out_file == NULL)
	{
		char locname[1024];
		char ss[1024];

		snprintf(locname, sizeof(locname), "%s%s%s%s", results_folder, "\\", "Res_Exa_", inst_name);
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

inline bool EXAsolver::separate_heuristic_subtours(CPXCALLBACKCONTEXTptr context, CPXLONG contextid, double var_treshold, bool* cut_added)
{
	int status = 0;
	double val;

			rop_build_sol_graph(rop_curr_sol.data(),var_treshold);
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
								if (this->rop_curr_sol[rop_v_y_n[node_subtour_subset[i]]] > yi_max_val)
								{
									yi_max = i;
									yi_max_val = this->rop_curr_sol[rop_v_y_n[node_subtour_subset[i]]];
								}
							}


							for (int i = 1; i <= node_subtour_subset[0]; ++i)
							{
								if (i != yi_max)
								{

									matind[nzc] = this->rop_v_y_n[node_subtour_subset[i]];
									matval[nzc] = -1;

									++nzc;

									lhs_v -= this->rop_curr_sol[rop_v_y_n[node_subtour_subset[i]]];

								}
							}
							for (int i = 1; i < node_subtour_subset[0]; ++i)
								for (int j = i + 1; j <= node_subtour_subset[0]; ++j)
								{
									matind[nzc] = v_x_ij(node_subtour_subset[i], node_subtour_subset[j]);
									matval[nzc] = 1;
									++nzc;


									lhs_v += this->rop_curr_sol[v_x_ij(node_subtour_subset[i], node_subtour_subset[j])];
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
									if (status)checkCPXstatus(status, &cpx_env_op, &cpx_lp_op);
								}
								else
									if (contextid == CPX_CALLBACKCONTEXT_RELAXATION)
									{
										//std::cout << "Cut" << std::endl;
										status = CPXcallbackaddusercuts(context, 1, nzc, rhs,
											sense, matbeg, matind.data(), matval.data(), force, local);
										if (status)checkCPXstatus(status, &cpx_env_op, &cpx_lp_op);
									}
							}
					}
			}
	return true;
}



inline bool EXAsolver::separate_logic_cuts(CPXCALLBACKCONTEXTptr context, CPXLONG contextid, bool* cut_added)
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
							rop_curr_sol[v_x_ij(i, j)] > rop_curr_sol[rop_v_y_n[i]] + cEPS)
						{
							*cut_added = true;

							// Add the cuts for i on y and x. 
							
								rhs[0] = 0;
								sense[0] = 'L';
								nzc = 0;
								lhs_v = 0;

								matind[nzc] = v_x_ij( i, j);
								matval[nzc] = 1;
								lhs_v += rop_curr_sol[v_x_ij( i, j)];
								++nzc;

								matind[nzc] = rop_v_y_n[i];
								matval[nzc] = -1;
								lhs_v -= rop_curr_sol[rop_v_y_n[i]];
								++nzc;

								if (lhs_v > rhs[0] + cEPS)
								{
									status = CPXcallbackaddusercuts(context, 1, nzc, rhs,
										sense, matbeg, matind.data(), matval.data(), force, local);

									if (status)
									{
										if (status)checkCPXstatus(status, &cpx_env_op, &cpx_lp_op);
										//fprintf(stdout, "Cut Logic!\n")

										return false;
									}

								}


							
						}
						if (j != end_node
							&&
							rop_curr_sol[v_x_ij(i, j)] > rop_curr_sol[rop_v_y_n[j]] + cEPS)
						{
							*cut_added = true;

						
								// Add the cuts for j, x and y and for each k
								rhs[0] = 0;
								sense[0] = 'L';
								nzc = 0;
								lhs_v = 0;

								matind[nzc] = v_x_ij(i, j);
								matval[nzc] = 1;
								++nzc;
								lhs_v += rop_curr_sol[v_x_ij( i, j)];

								matind[nzc] = rop_v_y_n[j];
								matval[nzc] = -1;
								lhs_v -= rop_curr_sol[rop_v_y_n[j]];
								++nzc;

								if (lhs_v > rhs[0] + mEPS)
								{
									status = CPXcallbackaddusercuts(context, 1, nzc, rhs,
										sense, matbeg, matind.data(), matval.data(), force, local);
									if (status)
									{
										if (status)checkCPXstatus(status, &cpx_env_op, &cpx_lp_op);
										//fprintf(stdout, "Cut Logic!\n")
										return false;
									}
								}
							
						}
					}
	
	return true;
}

inline void EXAsolver::reset_branc_and_cut_indicators()
{
	
		num_iter_curr_node = 0;
		num_cut_curr_iter = 0;
		prev_node = 0;
		curr_node = -1;
		lp_prev_iter = 0;

	
}

void EXAsolver::bm_reset_branc_and_cut_indicators()
{

	bm_num_iter_curr_node = 0;
	bm_num_cut_curr_iter = 0;
	bm_prev_node = 0;
	bm_curr_node = -1;
	bm_lp_prev_iter = 0;

}



inline bool EXAsolver::separate_exact_subtuors(CPXCALLBACKCONTEXTptr context, CPXLONG contextid, bool* cut_added)
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


			build_sol_directed_graph(rop_curr_sol.data());
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
								if (this->rop_curr_sol[rop_v_y_n[node_subtour_subset[i]]] > yi_max_val)
								{
									yi_max = i;
									yi_max_val = this->rop_curr_sol[rop_v_y_n[node_subtour_subset[i]]];
								}
							}


							for (int i = 1; i <= node_subtour_subset[0]; ++i)
							{
								if (i != yi_max)
								{

									matind[nzc] = rop_v_y_n[node_subtour_subset[i]];
									matval[nzc] = -1;



									lhs_val -= rop_curr_sol[rop_v_y_n[node_subtour_subset[i]]];
									++nzc;


								}
							}

							for (int i = 1; i < node_subtour_subset[0]; ++i)
								for (int j = i + 1; j <= node_subtour_subset[0]; ++j)
								{
									matind[nzc] = v_x_ij(node_subtour_subset[i], node_subtour_subset[j]);
									matval[nzc] = 1;
									lhs_val += rop_curr_sol[v_x_ij(node_subtour_subset[i], node_subtour_subset[j])];

									++nzc;

								}

							if (lhs_val > rhs[0] + cEPS)
							{
								status = CPXcallbackaddusercuts(context, 1, nzc, rhs,
									sense, matbeg, matind.data(), matval.data(), force, local);
								if(status)checkCPXstatus(status, &cpx_env_op, &cpx_lp_op);
							}
							if (status != 0)
							{
								return false;
							}
					}
				}
	return true;
}



inline int CPXPUBLIC EXAsolver::rop_general_callback(CPXCALLBACKCONTEXTptr context, CPXLONG contextid, void* userdata)
{
	bool cut_added = false;
	int status = 0;
	bool is_ok = true;
	EXAsolver* solver = (EXAsolver*)userdata;
	double val;
	double treshold;


	if (contextid == CPX_CALLBACKCONTEXT_CANDIDATE)
	{

		status = CPXcallbackgetcandidatepoint(context, solver->rop_curr_sol.data(), 0, CPXgetnumcols(solver->cpx_env_op, solver->cpx_lp_op) - 1, &val);
		if (status != 0)
		{
			solver->checkCPXstatus(status, &solver->cpx_env_op, &solver->cpx_lp_op);

			return status;
		}
	
		if (!cut_added)
			is_ok = solver->separate_heuristic_subtours(context, contextid, 0, &cut_added);


	}

	if (contextid == CPX_CALLBACKCONTEXT_RELAXATION)
	{
		/*Only CPX 12.10*/
		CPXLONG node_id, node_depth;

	
		

		CPXcallbackgetinfolong(context, CPXCALLBACKINFO_NODEUID, &node_id);
		if (status != 0)
		{
			solver->checkCPXstatus(status, &solver->cpx_env_op, &solver->cpx_lp_op);			return status;
		}

		CPXcallbackgetinfolong(context, CPXCALLBACKINFO_NODEDEPTH, &node_depth);
		if (status != 0)
		{
			solver->checkCPXstatus(status, &solver->cpx_env_op, &solver->cpx_lp_op);			return status;
		}





		status = CPXcallbackgetrelaxationpoint(context, solver->rop_curr_sol.data(), 0, CPXgetnumcols(solver->cpx_env_op, solver->cpx_lp_op) - 1, &val);
		if (status != 0)
		{
			solver->checkCPXstatus(status, &solver->cpx_env_op, &solver->cpx_lp_op);				return status;
		}


		if (node_id > 1) /*Do this after root node.*/
		{
			if (node_id == solver->curr_node)
			{
				++solver->num_iter_curr_node;
				//num_cut_curr_iter = 0;
				//prev_node = 0;
				solver->curr_node = node_id;

				if (fabs(val - solver->lp_prev_iter) < fabs(BandC_MIN_COEFF_LP_IMPR * solver->lp_prev_iter)
					&&
					solver->num_iter_curr_node >= BandC_MAX_NUM_ITER_NODE
					)
				{
					return 0; // Branch.
				}

		

			}
			else {
				solver->reset_branc_and_cut_indicators();
				solver->curr_node = node_id;
				solver->lp_prev_iter = val;
				solver->num_cut_curr_iter = 0;

			}
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


bool EXAsolver::store_xbar()
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
		for (int i = this->first_pnode; i <= this->last_pnode; ++i)
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
	for (int i = this->first_pnode; i <= this->last_pnode; ++i)
	{
	this->ccg_set_xib(num_ccg_xib, i) =	this->mlp_curr_sol[this->mlp_v_pr_xi[i]];
	}
	++num_ccg_xib;

	/*Also the xiY*/
	for (int y = 0; y < card_Y; ++y)
	{
		for (int i = this->first_pnode; i <= this->last_pnode; ++i)
		{
			this->ccg_set_xib(num_ccg_xib, i) = this->mlp_curr_sol[this->mlp_v_xi_y(y, i)];
		}
		++num_ccg_xib;
	}
	return true;
}

bool EXAsolver::bm_add_ccg_lb()
{

	/*lets manually add a scenario to check*/
	//num_ccg_xib = 0; /*reset*/
	//for (int i = this->first_pnode; i <= 5; ++i)
	//{
	//	this->ccg_set_xib(num_ccg_xib, i) = 0.2;
	//}
	//this->ccg_set_xib(num_ccg_xib, 1) = 0; 
	//this->ccg_set_xib(num_ccg_xib, 2) = 0.20;
	//this->ccg_set_xib(num_ccg_xib, 3) = 0.20;
	//this->ccg_set_xib(num_ccg_xib, 4) = 0;
	//this->ccg_set_xib(num_ccg_xib, 5) = 0.06667;
	//this->ccg_set_xib(num_ccg_xib, 6) = 0;
	//this->ccg_set_xib(num_ccg_xib, 7) = 0.0667;
	//this->ccg_set_xib(num_ccg_xib, 8) = 0.20;
	//this->ccg_set_xib(num_ccg_xib, 9) = 0.20;
	//this->ccg_set_xib(num_ccg_xib, 10) = 0.0667;
	//++num_ccg_xib;

	//this->ccg_set_xib(num_ccg_xib, 10) = 0;
	//this->ccg_set_xib(num_ccg_xib, 9) = 0.20;
	//this->ccg_set_xib(num_ccg_xib, 8) = 0.20;
	//this->ccg_set_xib(num_ccg_xib, 7) = 0;
	//this->ccg_set_xib(num_ccg_xib, 6) = 0.06667;
	//this->ccg_set_xib(num_ccg_xib, 5) = 0;
	//this->ccg_set_xib(num_ccg_xib, 4) = 0.0667;
	//this->ccg_set_xib(num_ccg_xib, 3) = 0.20;
	//this->ccg_set_xib(num_ccg_xib, 2) = 0.20;
	//this->ccg_set_xib(num_ccg_xib, 1) = 0.0667;
	//++num_ccg_xib;


	//this->ccg_set_xib(num_ccg_xib, 10) = 0.20;
	//this->ccg_set_xib(num_ccg_xib, 9) = 0.20;
	//this->ccg_set_xib(num_ccg_xib, 8) = 0.20;
	//this->ccg_set_xib(num_ccg_xib, 7) = 0.20;
	//this->ccg_set_xib(num_ccg_xib, 6) = 0.20;
	//this->ccg_set_xib(num_ccg_xib, 5) = 0;
	//this->ccg_set_xib(num_ccg_xib, 4) = 0;
	//this->ccg_set_xib(num_ccg_xib, 3) = 0;
	//this->ccg_set_xib(num_ccg_xib, 2) = 0;
	//this->ccg_set_xib(num_ccg_xib, 1) = 0;
	//++num_ccg_xib;

	this->bm_add_ccg_columns();

	this->bm_add_ccg_constraints();

	return true;
}

bool EXAsolver::bm_add_ccg_columns()
{
	unsigned int num_cols = CPXgetnumcols(this->cpx_env_bm, this->cpx_lp_bm);

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
		for (int i = first_pnode; i <= last_pnode; ++i)
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

		/*x( \xi, i ,j)*/
		for (int i = start_node; i <= end_node; ++i)
			for (int j = start_node; j <= end_node; ++j)
				if (this->edge_exists(i, j))
				{
				single_node_dur_1 = 0;
				single_node_dur_2 = 0;

				this->bm_v_x_xib_ij(xi, i, j) = num_cols;

				obj[0] = 0;

				vartype[0] = 'B';

				sprintf(varsname[0], "x_xi%d_i%d_j%d",xi ,i, j);
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

				status = CPXnewcols(this->cpx_env_bm, this->cpx_lp_bm, 1, obj, lb, ub, vartype, varsname);
				if (checkCPXstatus(status, &cpx_env_bm, &cpx_lp_bm)) goto TERMINATE;
				++num_cols;
				}

		///*lambda (\xi, r) r, num rows of matrix A*/
		for (unsigned short r = 0; r < this->num_us_constr; ++r)
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
		for (unsigned short n = this->first_pnode; n <= this->last_pnode; ++n)
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

bool EXAsolver::bm_add_ccg_constraints()
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


	// Degree constraints on profitable nodes
	for (int xi = 0; xi < this->num_ccg_xib; ++xi)
	{
		nzc = 0;
		rhs[0] = 0;
		sense[0] = 'E';

		for (int i = this->first_pnode; i <= this->last_pnode; ++i)
		{
			sprintf(cnstrname[0], "Degrees_i%d_xi%d", i, xi);
			nzc = 0;
			for (int ii = i + 1; ii <= end_node; ++ii)
				if (this->edge_exists(i, ii))
				{
					matind[nzc] = this->bm_v_x_xib_ij(xi,i, ii);
					matval[nzc] = 1;
					++nzc;
				}
			for (int ii = start_node; ii < i; ++ii)
				if (this->edge_exists(ii, i))
				{
					matind[nzc] = this->bm_v_x_xib_ij(xi,ii, i);
					matval[nzc] = 1;
					++nzc;
				}

			matind[nzc] = bm_v_y_xib_i(xi,i);
			matval[nzc] = -2;
			++nzc;

			status = CPXaddrows(cpx_env_bm, cpx_lp_bm, 0, 1, nzc, rhs, sense, matbeg, matind.data(), matval.data(), NULL, cnstrname);
			if (checkCPXstatus(status, &cpx_env_bm, &cpx_lp_bm)) goto TERMINATE;
		}

		// Degree constraint  on start node 
		rhs[0] = 1;
		sense[0] = 'E';
		nzc = 0;

		sprintf(cnstrname[0], "DegreeDepotOut_xi%d", xi);
		for (int j = this->start_node; j <= this->end_node; ++j)
			if (this->edge_exists(start_node, j))
			{
				matind[nzc] = bm_v_x_xib_ij(xi,this->start_node, j);
				matval[nzc] = 1;
				++nzc;

			}

		status = CPXaddrows(cpx_env_bm, cpx_lp_bm, 0, 1, nzc, rhs, sense, matbeg, matind.data(), matval.data(), NULL, cnstrname);
		if (checkCPXstatus(status, &cpx_env_bm, &cpx_lp_bm)) goto TERMINATE;


		// Degree constraint  on end node 
		rhs[0] = 1;
		sense[0] = 'E';
		nzc = 0;
		sprintf(cnstrname[0], "DegreeDepotIn_xi%d", xi);

		for (int j = start_node; j <= end_node; ++j)
			if (this->edge_exists(j, end_node))
			{
				matind[nzc] = bm_v_x_xib_ij(xi,j, end_node);
				matval[nzc] = 1;
				++nzc;

			}

		status = CPXaddrows(cpx_env_bm, cpx_lp_bm, 0, 1, nzc, rhs, sense, matbeg, matind.data(), matval.data(), NULL, cnstrname);
		if (checkCPXstatus(status, &cpx_env_bm, &cpx_lp_bm)) goto TERMINATE;


		nzc = 0;
		rhs[0] = max_dur;
		sense[0] = 'L';
		sprintf(cnstrname[0], "Max_Dur_xi%d", xi);
		for (int i = start_node; i < end_node; ++i)
			for (int j = start_node + 1; j <= end_node; ++j)
				if (this->edge_exists(i, j))
				{
					matind[nzc] = bm_v_x_xib_ij(xi,i, j);
					matval[nzc] = t_ij(i, j);
					++nzc;
				}

		status = CPXaddrows(cpx_env_bm, cpx_lp_bm, 0, 1, nzc, rhs, sense, matbeg, matind.data(), matval.data(), NULL, cnstrname);
		if (checkCPXstatus(status, &cpx_env_bm, &cpx_lp_bm)) goto TERMINATE;


		/*Dual constraints*/
		nzc = 0;
		rhs[0] = 0;
		sense[0] = 'E';
		for (int i = first_pnode; i <= last_pnode; ++i)
		{
			sprintf(cnstrname[0], "Balance_i%d_xi%d", i, xi);

			nzc = 0;
			rhs[0] = 0;
			sense[0] = 'E';

			matind[nzc] = this->bm_v_y_xib_i(xi, i);
			matval[nzc] = -ppxi[i];
			++nzc;

		

			
			matind[nzc] = this->bm_v_sig_xib_i(xi, i);
			matval[nzc] = 1;
			++nzc;
	/*				matind[nzc] = this->bm_v_sig_xib_i(xi, i);
		matval[nzc] = -1;
		++nzc;*/

			for (int l = 0; l < this->num_us_constr; ++l)
				if (fabs(A(l, i)) > eps_coeff)
				{
					matind[nzc] = this->bm_v_lamb_xib_i(xi, l);
					matval[nzc] = A(l, i);
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

		for (int i = first_pnode; i <= last_pnode; ++i)
		{
			matind[nzc] = this->bm_v_y_xib_i(xi, i);
			matval[nzc] = -this->pp[i];
			++nzc;

			matind[nzc] = this->bm_v_sig_xib_i(xi, i);
			matval[nzc] = -ccg_set_xib(xi, i);
			++nzc;
		}

		for (int l = 0; l < this->num_us_constr; ++l)
			{
				matind[nzc] = this->bm_v_lamb_xib_i(xi, l);
				matval[nzc] = -this->b[l];
				++nzc;
			}
		status = CPXaddrows(cpx_env_bm, cpx_lp_bm, 0, 1, nzc, rhs, sense, matbeg, matind.data(), matval.data(), NULL, cnstrname);
		if (checkCPXstatus(status, &cpx_env_bm, &cpx_lp_bm)) goto TERMINATE;


		/*Product w-sigma constraints*/
		for (int i = first_pnode; i <= last_pnode; ++i)
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

bool EXAsolver::bm_separate_heuristic_subtours(CPXCALLBACKCONTEXTptr context, CPXLONG contextid, double var_treshold, bool* cut_added, int xi_bar)
{
	int status = 0;
	double val;

	bm_build_sol_graph(bm_curr_sol.data(), var_treshold, xi_bar);
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


			for (int xi = 0; xi < this->num_ccg_xib; ++xi)
			{

				nzc = 0;
				rhs[0] = 0;
				lhs_v = 0;
				yi_max_val = 0;
				yi_max = 1;


					for (int i = 1; i <= node_subtour_subset[0]; ++i)
					{
						if (this->bm_curr_sol[bm_v_y_xib_i(xi, node_subtour_subset[i])] > yi_max_val)
						{
							yi_max = i;
							yi_max_val = this->bm_curr_sol[bm_v_y_xib_i(xi, node_subtour_subset[i])];
						}
					}


				for (int i = 1; i <= node_subtour_subset[0]; ++i)
				{
					if (i != yi_max)
					{

						matind[nzc] = this->bm_v_y_xib_i(xi, node_subtour_subset[i]);
						//this->rop_v_y_n[node_subtour_subset[i]];
						matval[nzc] = -1;

						++nzc;

						lhs_v -= this->bm_curr_sol[bm_v_y_xib_i(xi, node_subtour_subset[i])];

					}
				}
				for (int i = 1; i < node_subtour_subset[0]; ++i)
					for (int j = i + 1; j <= node_subtour_subset[0]; ++j)
					{
						matind[nzc] = this->bm_v_x_xib_ij(xi, node_subtour_subset[i], node_subtour_subset[j]);
						matval[nzc] = 1;
						++nzc;


						lhs_v += this->bm_curr_sol[bm_v_x_xib_ij(xi, node_subtour_subset[i], node_subtour_subset[j])];
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
						if (status)checkCPXstatus(status, &cpx_env_bm, &cpx_lp_bm);
					}
					else
						if (contextid == CPX_CALLBACKCONTEXT_RELAXATION)
						{
							//std::cout << "Cut" << std::endl;
							status = CPXcallbackaddusercuts(context, 1, nzc, rhs,
								sense, matbeg, matind.data(), matval.data(), force, local);
							if (status)checkCPXstatus(status, &cpx_env_bm, &cpx_lp_bm);
						}

					*cut_added = true;
				}
			}
			}
	}
	return true;
}

bool EXAsolver::bm_separate_exact_subtuors(CPXCALLBACKCONTEXTptr context, CPXLONG contextid, bool* cut_added, int xi_bar)
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


	bm_build_sol_directed_graph(bm_curr_sol.data(), xi_bar);
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
				yi_max_val = 0;
				yi_max = 1;

				for (int xi = 0; xi < this->num_ccg_xib; ++xi)
					{


					nzc = 0;
					lhs_val = 0;
					yi_max_val = 0;
					yi_max = 1;

					for (int i = 1; i <= node_subtour_subset[0]; ++i)
					{
						if (this->bm_curr_sol[bm_v_y_xib_i(xi, node_subtour_subset[i])] > yi_max_val)
						{
							yi_max = i;
							yi_max_val = this->bm_curr_sol[bm_v_y_xib_i(xi, node_subtour_subset[i])];
						}
					}


				for (int i = 1; i <= node_subtour_subset[0]; ++i)
				{
					if (i != yi_max)
					{

						matind[nzc] = bm_v_y_xib_i(xi, node_subtour_subset[i]);
						//rop_v_y_n[node_subtour_subset[i]];
						matval[nzc] = -1;



						lhs_val -= bm_curr_sol[bm_v_y_xib_i(xi, node_subtour_subset[i])];
						//rop_curr_sol[rop_v_y_n[node_subtour_subset[i]]];
						++nzc;


					}
				}

				for (int i = 1; i < node_subtour_subset[0]; ++i)
					for (int j = i + 1; j <= node_subtour_subset[0]; ++j)
					{
						matind[nzc] = bm_v_x_xib_ij(xi, node_subtour_subset[i], node_subtour_subset[j]);
						//v_x_ij(node_subtour_subset[i], node_subtour_subset[j]);
						matval[nzc] = 1;
						lhs_val += bm_curr_sol[bm_v_x_xib_ij(xi, node_subtour_subset[i], node_subtour_subset[j])];

						++nzc;

					}

				if (lhs_val > rhs[0] + cEPS)
				{
					status = CPXcallbackaddusercuts(context, 1, nzc, rhs,
						sense, matbeg, matind.data(), matval.data(), force, local);
					if (status)checkCPXstatus(status, &cpx_env_bm, &cpx_lp_bm);

					*cut_added = true;

				}
				if (status != 0)
				{
					return false;
				}
			}
			}
		}
	return true;
}


inline bool EXAsolver::store_set_of_tight_y()
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
		std::cout << "Error while storing the tight y in bool EXAsolver::store_tight_y()" << std::endl;
		getchar();
		return false;
	}
	else
		return true;
	
}

inline void EXAsolver::store_tight_y(int y)
{
	if (num_tight_y >= tight_y.size())
	{
		tight_y.resize(2 * num_tight_y);
		tight_xibar.resize(2 * num_tight_y);
		tight_xi_y.resizeMatrix2D(2 * num_tight_y, this->num_nodes);

	}


	this->tight_y[this->num_tight_y] = y;


	for (int i = this->first_pnode; i <= this->last_pnode; ++i)
	{
		this->tight_xi_y(num_tight_y, i) = this->mlp_curr_sol[this->mlp_v_xi_y(y,i)];
		this->tight_xibar[i] = this->mlp_curr_sol[this->mlp_v_pr_xi[i]];
	}

	++num_tight_y;
}

bool EXAsolver::read_distance_matrix(FILE* inp)
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

				t_ij(j, i) =(int) (t_ij(j, i) * 100.00);
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

bool EXAsolver::read_instance(char* file, char* config_file)
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
		disc_cost_n.assign(this->num_nodes,0);

		// Travel times
		this->t_ij.assignMatrix2D(this->num_nodes, this->num_nodes, 0);

		// Node existence
		this->edge_exists.assignMatrix2D(this->num_nodes, this->num_nodes, false);

		this->read_distance_matrix(inp);

	}


	
	return true;
}


bool EXAsolver::init_data_structures()
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




	int num_vrs_op = 0;
	int num_vrs_lp = 0;

	this->rop_v_y_n.resize(this->num_nodes);
	num_vrs_op += rop_v_y_n.size();

	this->v_x_ij.resizeMatrix2D(this->num_nodes, this->num_nodes);
	num_vrs_op += v_x_ij.getNumElem();

	this->v_beta.resize(4 * this->num_nodes);
	num_vrs_op += v_beta.size();

	this->v_gamma.resize(this->num_nodes);
	num_vrs_op += v_gamma.size();



	/*uncertainty set data structure: for now I assume the uncertainty set is defined by num_nodes*2 constr + 1 (one for each node plus budget.)*/

	A.assignMatrix2D((4 * this->num_nodes) + 1, this->num_nodes, 0.00);

	b.resize((4 * this->num_nodes) + 1);

	ppxi.resize(this->num_nodes);

	pp.resize(this->num_nodes); // part not function of P.


	/*Now we set for the LP*/
	this->mlp_num_y = 0;
	card_Y = 0;
	Y.resizeMatrix2D(this->num_nodes, this->num_nodes); // Each row is a solution of the OP.
	y_in_mlp.assign(this->num_nodes, false);
	this->mlp_v_pr_xi.resize(this->num_nodes);

	num_vrs_lp += mlp_v_pr_xi.size();
	this->mlp_v_xi_y.resizeMatrix2D(this->num_nodes, this->num_nodes); /*We start with a size of one route for each node.*/
	mlp_cnstr_tau_y.resize(this->num_nodes);

	num_vrs_lp += mlp_v_xi_y.getNumElem();
	
	bm_v_w_i.resize(this->num_nodes);

	this->bm_use_scenarios_ccg_xib =  cfg->getValueOfKey<bool>("EXACT_OP_USE_SCENARIOS");

	this->bm_use_info_cuts = cfg->getValueOfKey<bool>("EXACT_USE_INFO_CUTS");

	if (bm_use_scenarios_ccg_xib)
	{

		ccg_set_xib.resizeMatrix2D(MAX_NUM_SCENARIOS, this->num_nodes);
		bm_v_y_xib_i.resizeMatrix2D(MAX_NUM_SCENARIOS, this->num_nodes);
		bm_v_x_xib_ij.resizeMatrix3D(MAX_NUM_SCENARIOS, this->num_nodes, this->num_nodes);
		bm_v_lamb_xib_i.resizeMatrix2D(MAX_NUM_SCENARIOS, b.size());
		bm_v_sig_xib_i.resizeMatrix2D(MAX_NUM_SCENARIOS, this->num_nodes);
	}
	


	this->num_tight_y = 0;
	tight_y.resize(this->num_nodes);
	this->tight_xi_y.resizeMatrix2D(this->num_nodes, this->num_nodes);
	this->tight_xibar.resize(this->num_nodes);

	/*matind and matval*/
	matind.resize(std::max(num_vrs_op, num_vrs_lp));
	matval.resize(std::max(num_vrs_op, num_vrs_lp));

	rop_curr_sol.resize(std::max(num_vrs_op, num_vrs_lp));
	mlp_curr_sol.resize(std::max(num_vrs_op, num_vrs_lp));

	bm_curr_sol.resize(10*this->num_nodes+1);



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


	/*
	cut counters to 0
	*/
	this->reset_branc_and_cut_indicators();
	this->bm_reset_branc_and_cut_indicators();


	num_ccg_xib = 0;

this->num_logic_benders_cuts =0;

bm_num_stored_logic_cuts = 0;
bm_logic_cuts_pool_set.resize(static_cast<std::vector<uint64_t, std::allocator<uint64_t>>::size_type>(10) * this->num_nodes);
bm_logic_cuts_pool_phi_val.resize( static_cast<std::vector<double, std::allocator<double>>::size_type>(10)		*     this->num_nodes);





	return true;
}


bool EXAsolver::clean_data_structure()
{
	delete cfg;

	/*Clean cplex op*/
	int status = 0;
	status = CPXfreeprob(this->cpx_env_op, &this->cpx_lp_op);
	if (checkCPXstatus(status, &this->cpx_env_op, &this->cpx_lp_op)) return false;
	status = CPXcloseCPLEX(&this->cpx_env_op);
	if (checkCPXstatus(status, &this->cpx_env_op, &this->cpx_lp_op)) return false;

	/*Clean mlp*/
	status = CPXfreeprob(this->cpx_env_mlp, &this->cpx_lp_mlp);
	if (checkCPXstatus(status, &this->cpx_env_mlp, &this->cpx_lp_mlp)) return false;
	status = CPXcloseCPLEX(&this->cpx_env_mlp);
	if (checkCPXstatus(status, &this->cpx_env_mlp, &this->cpx_lp_mlp)) return false;


	/*Clean bm*/
	status = CPXfreeprob(this->cpx_env_bm, &this->cpx_lp_bm);
	if (checkCPXstatus(status, &this->cpx_env_bm, &this->cpx_lp_bm)) return false;
	status = CPXcloseCPLEX(&this->cpx_env_bm);
	if (checkCPXstatus(status, &this->cpx_env_bm, &this->cpx_lp_bm)) return false;

	if (out_file != NULL)
		fclose(out_file);





	return true;
}



bool EXAsolver::define_uncertainty_set()
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


bool EXAsolver::define_uncertainty_set_type_1()
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

bool EXAsolver::define_uncertainty_set_type_2()
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

bool EXAsolver::define_uncertainty_set_type_3()
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

bool EXAsolver::define_uncertainty_set_type_4()
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





bool EXAsolver::load_and_init(std::vector<bool> w_val)
{
	bool is_ok = true;

	this->w_val = w_val;

	this->bm_cut_lb = - DBL_MAX;
	this->bm_best_lb = - DBL_MAX;


	is_ok = init_data_structures();
	if (!is_ok)
	{
		goto TERMINATE;
	}

	is_ok = this->define_uncertainty_set();
	if (!is_ok)
	{
		goto TERMINATE;
	}

	is_ok = mlp_build_model();
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

bool EXAsolver::run_sub_problem_solver(double* obj_val)
{
	bool is_ok = true;

	is_ok = mlp_solve_model();
	if (!is_ok)
	{
		goto TERMINATE;
	}


	do {

		is_ok = rop_update_scenario_parameter_model();
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
			}
			break;

		}

	


	} while (true);
	//while (fabs(this->mlp_obj_val - this->rop_obj_val) > EXA_EPS);


	*obj_val = mlp_obj_val;

TERMINATE:
	return is_ok;
}

bool EXAsolver::rop_build_model()
{
	bool is_ok = true;
	int status;
	int mip_status;
	CPXLONG      contextmask = 0;
	char errbuf[CPXMESSAGEBUFSIZE];
	int check_mip_stat;

	double start = 0;

	double end = 0;

	is_ok = init_cplex(&this->cpx_env_op, &this->cpx_lp_op);
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

	contextmask |= CPX_CALLBACKCONTEXT_CANDIDATE;
	contextmask |= CPX_CALLBACKCONTEXT_RELAXATION;
	if (contextmask != 0) {
		//	 We are done and now we register our callback function. 
		status = CPXcallbacksetfunc(cpx_env_op, cpx_lp_op, contextmask, rop_general_callback, this);
		if (status != 0) {
			fprintf(stderr, "Failed to add callback: %s\n",
				CPXgeterrorstring(cpx_env_op, status, errbuf));
			goto TERMINATE;
		}
	}



	

TERMINATE:
	return is_ok;
}

bool EXAsolver::rop_set_cplex()
{
	bool mod_stat = true;

	int status = 0;

	char errbuf[CPXMESSAGEBUFSIZE];


	status = CPXsetintparam(this->cpx_env_op, CPX_PARAM_SCRIND, CPX_OFF);
	if (checkCPXstatus(status, &cpx_env_op, &cpx_lp_op)) return false;
	status = CPXsetintparam(this->cpx_env_op, CPX_PARAM_THREADS, 1);
	if (checkCPXstatus(status, &cpx_env_op, &cpx_lp_op)) return false;
	status = CPXchgobjsen(this->cpx_env_op, this->cpx_lp_op, CPX_MIN);
	if (checkCPXstatus(status, &cpx_env_op, &cpx_lp_op)) return false;
	status = CPXsetintparam(this->cpx_env_op, CPX_PARAM_NUMERICALEMPHASIS, CPX_ON);
	if (checkCPXstatus(status, &cpx_env_op, &cpx_lp_op)) return false;
	status = CPXsetdblparam(this->cpx_env_op, CPXPARAM_Simplex_Tolerances_Feasibility, 1e-9);
	if (checkCPXstatus(status, &cpx_env_op, &cpx_lp_op)) return false;
	//status = CPXsetintparam(cpx_env, CPX_PARAM_PREIND, CPX_OFF);
	//if (checkCPXstatus(status))return false;

	/*To be extra numerical accurate*/
	status = CPXsetdblparam(this->cpx_env_op, CPXPARAM_Simplex_Tolerances_Optimality, 1e-9);
	if (checkCPXstatus(status, &cpx_env_op, &cpx_lp_op)) return false;

	status = CPXsetdblparam(this->cpx_env_op, CPXPARAM_Simplex_Tolerances_Feasibility, 1e-9);
	if (checkCPXstatus(status, &cpx_env_op, &cpx_lp_op)) return false;

	status = CPXsetdblparam(this->cpx_env_op, CPXPARAM_MIP_Tolerances_MIPGap, 1e-9);
	if (checkCPXstatus(status, &cpx_env_op, &cpx_lp_op)) return false;
	status = CPXsetdblparam(this->cpx_env_op, CPXPARAM_MIP_Tolerances_Integrality, 0);
	if (checkCPXstatus(status, &cpx_env_op, &cpx_lp_op)) return false;

	double time_limit = this->cfg->getValueOfKey<double>("TIME_LIMIT");

	status = CPXsetdblparam(cpx_env_op, CPXPARAM_TimeLimit, time_limit);
	if (checkCPXstatus(status, &cpx_env_op, &cpx_lp_op)) return false;

	/*
	status = CPXsetintparam(cpx_env, CPXPARAM_MIP_Strategy_VariableSelect, CPX_VARSEL_STRONG);
	if (checkCPXstatus(status)) return false;
	*/


	return true;
}

bool EXAsolver::rop_add_vars()
{
	unsigned short num_cols = CPXgetnumcols(this->cpx_env_op, this->cpx_lp_op);

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

		for (int i = first_pnode; i <= last_pnode; ++i)
		{


			this->rop_v_y_n[i] = num_cols;

			obj[0] = this->pp[i];

			vartype[0] = 'B';

			sprintf(varsname[0], "y_i%d", i);
			lb[0] = 0;
			ub[0] = 1;




			status = CPXnewcols(this->cpx_env_op, this->cpx_lp_op, 1, obj, lb, ub, vartype, varsname);
			if (checkCPXstatus(status, &cpx_env_op, &cpx_lp_op)) goto TERMINATE;
			++num_cols;
		}

	/*x_ij*/

		for (int i = start_node; i <= end_node; ++i)
			for (int j = start_node; j <= end_node; ++j)
				if (this->edge_exists(i, j))
				{
					single_node_dur_1 = 0;
					single_node_dur_2 = 0;

					this->v_x_ij(i, j) = num_cols;

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

					status = CPXnewcols(this->cpx_env_op, this->cpx_lp_op, 1, obj, lb, ub, vartype, varsname);
					if (checkCPXstatus(status, &cpx_env_op, &cpx_lp_op)) goto TERMINATE;
					++num_cols;
				}

	/*Dual variables */
	for (unsigned short r = 0; r < this->num_us_constr; ++r)
	{

		v_beta[r] = num_cols;

		obj[0] = this->b[r];

		vartype[0] = 'C';

		sprintf(varsname[0], "beta_r%d", r);

		lb[0] = 0;
		ub[0] = CPX_INFBOUND;

		status = CPXnewcols(this->cpx_env_op, this->cpx_lp_op, 1, obj, lb, ub, vartype, varsname);
		if (checkCPXstatus(status, &cpx_env_op, &cpx_lp_op)) goto TERMINATE;		++num_cols;
	}


		for (unsigned short n = this->first_pnode; n <= this->last_pnode; ++n)
		{

			v_gamma[n] = num_cols;

			if (/*this->w_val[n] 
				&&*/
				fabs(this->mlp_curr_sol[mlp_v_pr_xi[n]]) >= 1e-9
				)
			{

				//rounded_to_precision = std::round(this->mlp_curr_sol[this->v_pr_xi[n]] * 1e9) / 1e9;

				

				obj[0] = this->mlp_curr_sol[mlp_v_pr_xi[n]];
				//obj[0] = rounded_to_precision;
			}
			else
				obj[0] = 0;

			vartype[0] = 'C';

			sprintf(varsname[0], "gamma_n%d",n);

		
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

			status = CPXnewcols(this->cpx_env_op, this->cpx_lp_op, 1, obj, lb, ub, vartype, varsname);
			if (checkCPXstatus(status, &cpx_env_op, &cpx_lp_op)) goto TERMINATE;		++num_cols;
		}


TERMINATE:

	if (status)
	{
		return false;
	}

	return true;
}

bool EXAsolver::rop_add_cnstr()
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
					matind[nzc] = this->v_x_ij(i, ii);
					matval[nzc] = 1;
					++nzc;
				}
			for (int ii = start_node; ii < i; ++ii)
				if (this->edge_exists(ii, i))
				{
					matind[nzc] = this->v_x_ij(ii, i);
					matval[nzc] = 1;
					++nzc;
				}

	
			matind[nzc] = rop_v_y_n[i];
			matval[nzc] = -2;
			++nzc;

			status = CPXaddrows(cpx_env_op, cpx_lp_op, 0, 1, nzc, rhs, sense, matbeg, matind.data(), matval.data(), NULL, cnstrname);
			if (checkCPXstatus(status, &cpx_env_op, &cpx_lp_op)) goto TERMINATE;
		}

	// Degree constraint  on start node 
	rhs[0] = 1;
	sense[0] = 'E';
	nzc = 0;

		sprintf(cnstrname[0], "DegreeDepotOut");
		for (int j = this->start_node; j <= this->end_node; ++j)
			if (this->edge_exists(start_node, j))
			{
				//matind[nzc] = g_v_xk_kij[s][g_start_node][j];
				matind[nzc] = v_x_ij(this->start_node, j);
				matval[nzc] = 1;
				++nzc;

			}
		status = CPXaddrows(cpx_env_op, cpx_lp_op, 0, 1, nzc, rhs, sense, matbeg, matind.data(), matval.data(), NULL, cnstrname);
		if (checkCPXstatus(status, &cpx_env_op, &cpx_lp_op)) goto TERMINATE;

	
	// Degree constraint  on end node 
	rhs[0] = 1;
	sense[0] = 'E';
	nzc = 0;
		sprintf(cnstrname[0], "DegreeDepotIn");

		for (int j = start_node; j <= end_node; ++j)
			if (this->edge_exists(j, end_node))
			{
				//g_matind[nzc] = g_v_xk_kij[s][j][g_end_node];
				matind[nzc] = v_x_ij(j, end_node);
				matval[nzc] = 1;
				++nzc;

			}
		status = CPXaddrows(cpx_env_op, cpx_lp_op, 0, 1, nzc, rhs, sense, matbeg, matind.data(), matval.data(), NULL, cnstrname);
		if (checkCPXstatus(status, &cpx_env_op, &cpx_lp_op)) goto TERMINATE;
	

	nzc = 0;
	rhs[0] = max_dur;
	sense[0] = 'L';
			sprintf(cnstrname[0], "Max_Dur");
		for (int i = start_node; i < end_node; ++i)
			for (int j = start_node + 1; j <= end_node; ++j)
				if (this->edge_exists(i, j))
				{
					//	g_matind[nzc] = g_v_xk_kij[s][i][j];
					matind[nzc] = v_x_ij(i, j);
					matval[nzc] = t_ij(i, j);
					++nzc;
				}
		status = CPXaddrows(cpx_env_op, cpx_lp_op, 0, 1, nzc, rhs, sense, matbeg, matind.data(), matval.data(), NULL, cnstrname);
		if (checkCPXstatus(status, &cpx_env_op, &cpx_lp_op)) goto TERMINATE;


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

			matind[nzc] = rop_v_y_n[i];
			matval[nzc] = -ppxi[i];
			++nzc;


			matind[nzc] = v_gamma[i];
			matval[nzc] = 1; // Becase we cnage the bound of gamma when w = 0.
			//if (w_val[i])
			//{
			//	//matval[nzc] = -1;
			//	matval[nzc] = 1;
			//}
			//else
			//{
			//	matval[nzc] = 0;
			//}
			++nzc;

			for (int l = 0; l < this->num_us_constr; ++l)
				if (fabs(A(l, i)) > eps_coeff)
				{
					matind[nzc] = v_beta[l];
					matval[nzc] = A(l, i);
					++nzc;
				}

			status = CPXaddrows(cpx_env_op, cpx_lp_op, 0, 1, nzc, rhs, sense, matbeg, matind.data(), matval.data(), NULL, cnstrname);
			if (checkCPXstatus(status, &cpx_env_op, &cpx_lp_op)) goto TERMINATE;
		}



TERMINATE:
	if (status)
		return false;
	else
		return true;
}

bool EXAsolver::rop_update_scenario_parameter_model()
{
	bool mod_stat = true;

	int status;
	int cnt = 0;
	double rounded_to_precision = 0;

	for (unsigned short n = this->first_pnode; n <= this->last_pnode; ++n)
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

	status = CPXchgobj(cpx_env_op, cpx_lp_op, cnt, matind.data(), matval.data());
	if (checkCPXstatus(status, &cpx_env_op, &cpx_lp_op)) goto TERMINATE;




TERMINATE:
	if (status)
		return false;
	else
		return true;
}

bool EXAsolver::rop_solve_model()
{
	bool is_ok = true;
	int status;
	int mip_status;
	CPXLONG      contextmask = 0;
	char errbuf[CPXMESSAGEBUFSIZE];
	int check_mip_stat;

//	status = CPXwriteprob(this->cpx_env_op, this->cpx_lp_op, "model_rop.lp", "lp");
//	if (checkCPXstatus(status, &cpx_env_op, &cpx_lp_op)) { is_ok = false; goto TERMINATE; };


	status = CPXmipopt(this->cpx_env_op, this->cpx_lp_op);
	if (checkCPXstatus(status, &cpx_env_op, &cpx_lp_op)) goto TERMINATE;



	check_mip_stat = CPXgetstat(cpx_env_op, cpx_lp_op);

	status = CPXgetobjval(cpx_env_op, cpx_lp_op, &this->rop_obj_val);
	if (checkCPXstatus(status, &cpx_env_op, &cpx_lp_op)) goto TERMINATE;


	if (this->mlp_obj_val > this->rop_obj_val + EXA_EPS)
	{
		status = CPXgetx(cpx_env_op, cpx_lp_op, rop_curr_sol.data(), 0, CPXgetnumcols(cpx_env_op, cpx_lp_op) - 1);
		if (checkCPXstatus(status, &cpx_env_op, &cpx_lp_op)) goto TERMINATE;


		if (Y.getNumRows() <= card_Y)
		{
			Y.resizeMatrix2D(card_Y * 2, num_nodes);
		}

		for (int i = first_pnode; i <= last_pnode; ++i)
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

bool EXAsolver::mlp_build_model()
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

bool EXAsolver::mlp_set_cplex()
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

bool EXAsolver::mlp_add_vars()
{
	unsigned short num_cols = CPXgetnumcols(this->cpx_env_mlp, this->cpx_lp_mlp);

	int status = 0;
	int vars_ind = 0;
	double lb[1], ub[1];
	double obj[1];
	char* varsname[1];
	char elname[1024];
	varsname[0] = elname;
	char vartype[1];

	int cnt;


	/*tau*/
	this->mlp_v_tau = num_cols;
	obj[0] = 1;
	vartype[0] = 'C';
	sprintf(varsname[0], "tau");
	lb[0] = -CPX_INFBOUND;
	//ub[0] = CPX_INFBOUND;
	ub[0] = 0; /*Minimum prize cannt be lower than 0.*/
	status = CPXnewcols(this->cpx_env_mlp, this->cpx_lp_mlp, 1, obj, lb, ub, NULL, varsname);
	if (checkCPXstatus(status, &cpx_env_mlp, &cpx_lp_mlp)) goto TERMINATE;
	++num_cols;


	/*ov_xi variables*/

	for (int i = first_pnode; i <= last_pnode; ++i)
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

bool EXAsolver::mlp_add_cnstr()
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

	for (int r = 0; r < this->num_us_constr; ++r)
	{
		sprintf(cnstrname[0], "Axi_leq_b_r%d", r);
		nzc = 0;
		rhs[0] = this->b[r];

		for (int i = this->first_pnode; i <= this->last_pnode; ++i)
		{

			matind[nzc] = this->mlp_v_pr_xi[i];

			if(fabs(this->A(r, i)) > eps_coeff)
				matval[nzc] = this->A(r, i);
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

bool EXAsolver::mlp_update_model()
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

bool EXAsolver::check_violation_in_Y()
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

bool EXAsolver::mlp_add_column_and_rows(int y)
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

inline bool EXAsolver::add_columns_y(int y)
{
	unsigned short num_cols = CPXgetnumcols(this->cpx_env_mlp, this->cpx_lp_mlp);
	
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
		mlp_v_xi_y.resizeMatrix2D(2 * card_Y, this->num_nodes);

		this->mlp_cnstr_tau_y.resize(2 * card_Y);
	}


	/*ov_xi variables*/
		for (int i = first_pnode; i <= last_pnode; ++i)
		{


			this->mlp_v_xi_y(y, i) = num_cols;

			obj[0] = 0;

			vartype[0] = 'C';

			sprintf(varsname[0], "xi_y%d_i%d", y, i);
			lb[0] = -CPX_INFBOUND;
			ub[0] = CPX_INFBOUND;


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

inline bool EXAsolver::add_rows_y(int y)
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
		sprintf(cnstrname[0], "tau_leq_Y%d", y);
		matind[nzc] = this->mlp_v_tau;
		matval[nzc] = 1;
		++nzc;
		for (int i = this->first_pnode; i <= this->last_pnode; ++i)
		{

			matind[nzc] = this->mlp_v_xi_y(y, i);
			matval[nzc] = -this->ppxi[i] * Y(y, i);
			++nzc;

			rhs[0] += Y(y, i) * this->pp[i];
		}

		status = CPXaddrows(cpx_env_mlp, cpx_lp_mlp, 0, 1, nzc, rhs, sense, matbeg, matind.data(), matval.data(), NULL, cnstrname);
		if (checkCPXstatus(status, &cpx_env_mlp, &cpx_lp_mlp)) goto TERMINATE;
	



	nzc = 0;
	rhs[0] = 0;
	sense[0] = 'L';
		for (int r = 0; r < this->num_us_constr; ++r)
		{
			nzc = 0;
			rhs[0] = this->b[r];
			sprintf(cnstrname[0], "balance_Y%d_r%d", y, r);

			for (int i = this->first_pnode; i <= this->last_pnode; ++i)
			{

				matind[nzc] = this->mlp_v_xi_y(y, i);

				if (fabs(this->A(r, i)) > eps_coeff)
					matval[nzc] = this->A(r, i);
				else
					matval[nzc] = 0;

				++nzc;
			}

			status = CPXaddrows(cpx_env_mlp, cpx_lp_mlp, 0, 1, nzc, rhs, sense, matbeg, matind.data(), matval.data(), NULL, cnstrname);
			if (checkCPXstatus(status, &cpx_env_mlp, &cpx_lp_mlp)) goto TERMINATE;
		}


	sense[0] = 'E';
		for (int i = this->first_pnode; i <= this->last_pnode; ++i)
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

inline bool EXAsolver::add_columns()
{
	unsigned short num_cols = CPXgetnumcols(this->cpx_env_mlp, this->cpx_lp_mlp);

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
		mlp_v_xi_y.resizeMatrix2D(2 * card_Y, this->num_nodes);

		mlp_cnstr_tau_y.resize(2 * card_Y);
	}


	/*ov_xi variables*/
	for(int y = mlp_num_y; y < card_Y; ++y)
	for (int i = first_pnode; i <= last_pnode; ++i)
	{


		this->mlp_v_xi_y(y,i) = num_cols;

		obj[0] = 0;

		vartype[0] = 'C';

		sprintf(varsname[0], "xi_y%d_i%d", y,i);
		lb[0] = -CPX_INFBOUND;
		ub[0] = CPX_INFBOUND;


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

inline bool EXAsolver::add_rows()
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

		for (int i = this->first_pnode; i <= this->last_pnode; ++i)
		{

			matind[nzc] = this->mlp_v_xi_y(y, i);
			matval[nzc] = -this->ppxi[i] * Y(y, i);
			++nzc;

			rhs[0] += Y(y, i) * this->pp[i];
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
	
for (int y = mlp_num_y; y < card_Y; ++y)
	for (int r = 0; r < this->num_us_constr; ++r)
	{
		nzc = 0;
		rhs[0] = this->b[r];
		sprintf(cnstrname[0], "balance_Y%d_r%d", y,r);

		for (int i = this->first_pnode; i <= this->last_pnode; ++i)
		{

			matind[nzc] = this->mlp_v_xi_y(y, i);

			if (fabs(this->A(r, i)) > eps_coeff)
				matval[nzc] = this->A(r, i);
			else
				matval[nzc] = 0;

			++nzc;
		}

		status = CPXaddrows(cpx_env_mlp, cpx_lp_mlp, 0, 1, nzc, rhs, sense, matbeg, matind.data(), matval.data(), NULL, cnstrname);
		if (checkCPXstatus(status, &cpx_env_mlp, &cpx_lp_mlp)) goto TERMINATE;
	}


sense[0] = 'E';
for (int y = mlp_num_y; y < card_Y; ++y)
	for (int i = this->first_pnode; i <= this->last_pnode; ++i)
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

bool EXAsolver::mlp_solve_model()
{
	bool is_ok = true;
	int status;
	int mip_status;
	CPXLONG      contextmask = 0;
	char errbuf[CPXMESSAGEBUFSIZE];
	int check_mip_stat;

	//status = CPXwriteprob(this->cpx_env_mlp, this->cpx_lp_mlp, "model_mlp.lp", "lp");
	//if (checkCPXstatus(status, &cpx_env_mlp, &cpx_lp_mlp)) { is_ok = false; goto TERMINATE; };

	
	status = CPXlpopt(this->cpx_env_mlp, this->cpx_lp_mlp);
	if (checkCPXstatus(status, &cpx_env_mlp, &cpx_lp_mlp)) goto TERMINATE;


	

	status = CPXgetx(cpx_env_mlp, cpx_lp_mlp, mlp_curr_sol.data(), 0, CPXgetnumcols(cpx_env_mlp, cpx_lp_mlp) - 1);
	if (checkCPXstatus(status, &cpx_env_mlp, &cpx_lp_mlp)) goto TERMINATE;

	status = CPXgetobjval(cpx_env_mlp, cpx_lp_mlp, &this->mlp_obj_val);
	if (checkCPXstatus(status, &cpx_env_mlp, &cpx_lp_mlp)) goto TERMINATE;




TERMINATE:
	return is_ok;
}

bool EXAsolver::subproblem_solver_update_discovery_parameters(const std::vector<bool>* w_val)
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

bool EXAsolver::run_exact_method(double* objective_value)
{
	bool is_ok = true;

	/*Initilizing the sep problem*/
	w_val.assign(this->num_nodes, true);
	
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

bool EXAsolver::print_solution_file()
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
		fprintf(out_file, "Uncertainty Set Type: %d\n", this->uncertainty_set_type);
		fprintf(out_file, "Uncertainty Set Parameter: %lf\n", this->uncertainty_parameter);
		fprintf(out_file, "Total Deterministic Profit in the Network: %lf\n", this->tot_deterministic_profit);




		fprintf(out_file, "++++ Solution of the Model ++++\n");
		fprintf(out_file, "-----------------------------------------------");
		fprintf(out_file, "\nObjective:\n");
		fprintf(out_file, "Best Lower Bound: %lf\n", bm_best_lb);
		fprintf(out_file, "Best upper bound: %lf\n", bm_best_ub);
		fprintf(out_file, "Elapsed time Branch-and-Cut: %lf\n", this->bm_branch_and_cut_time);
		fprintf(out_file, "Elapsed time Starting Cuts: %lf\n", this->bm_starting_cuts_time);

		fprintf(out_file, "Cplex Gap: %lf\n", this->cplex_gap);

		fprintf(out_file, "B&C nodes: %d\n", bandc_nodes);

		fprintf(out_file, "Number Logic Benders Cuts: %d\n", this->num_logic_benders_cuts);
		fprintf(out_file, "Number Starting Scenarios: %d\n", this->num_ccg_xib);
		fprintf(out_file, "Cardinality of Y: %d\n", this->card_Y);
		fprintf(out_file, "Number of Tight y in the sub-problem: %d\n", this->num_tight_y);






		fprintf(out_file, "\nVariables Master Problem:\n");
		fprintf(out_file, "\nPhi =%lf \n", bm_curr_sol[bm_v_phi]);

		fprintf(out_file, "\nDiscoveries: w_i (i is a profitable node). \n");
		for (int i = start_node + 1; i < end_node; ++i)
			if (bm_curr_sol[bm_v_w_i[i]] > mEPS)
			{
				fprintf(out_file, " w_[node:%d] = %lf\n", i, bm_curr_sol[bm_v_w_i[i]]);
				++num_discovery;
				discovery_cost += this->disc_cost_n[i];
			}

	
		fprintf(out_file, "\nScenarios and Tight Ys Sub-Problem\n");

		fprintf(out_file, "\nOverline Scenario: ovxi_i (i is a profitable node). \n");
		for (int i = start_node + 1; i < end_node; ++i)
			{
				fprintf(out_file, " ov_xi_[node:%d] = %lf\n", i, this->tight_xibar[i]);
			}

		fprintf(out_file, "\n Tigth Ys: \n");
		for (int ind_y = 0; ind_y < this->num_tight_y; ++ind_y)
		{
			fprintf(out_file, "\n y = %d { \n", tight_y[ind_y]);
			for (int i = start_node + 1; i < end_node; ++i)
			{
				fprintf(out_file, "Node %d: {\n", i);
				fprintf(out_file, "y_[node:%d] = %d\n", i, this->Y(this->tight_y[ind_y], i));
				fprintf(out_file, "xi_y_i_[y:%d][node:%d] = %lf\n",tight_y[ind_y], i, this->mlp_curr_sol[this->mlp_v_xi_y(tight_y[ind_y],i)]);
				fprintf(out_file, "} \n");

			}
			fprintf(out_file, "\n } \n");

		}




			//fprintf(out_file, "Policy-%d description\n", k);
			//int curr_node = start_node;
			//int prev_node = -1;
			//double route_length = 0;
			//fprintf(out_file, "<%d", curr_node);

			//while (curr_node != end_node)
			//{
			//	for (int j = start_node; j <= end_node; ++j)
			//		if (j != curr_node
			//			&&
			//			j != prev_node
			//			&&
			//			curr_sol[v_x_kij(k, std::min(curr_node, j), std::max(curr_node, j))] > mEPS)
			//		{
			//			route_length += t_ij(std::min(curr_node, j), std::max(curr_node, j));
			//			//[imin(curr_node, j)] [imax(curr_node, j)] ;

			//			prev_node = curr_node;
			//			curr_node = j;

			//			break;
			//		}
			//	fprintf(out_file, "-%d", curr_node);
			//}
			//fprintf(out_file, "> : Length = %lf\n", route_length);
			//fprintf(out_file, "_________________________________________________\n");
		


		fprintf(out_file, "Legenda_Model:");
		fprintf(out_file, "Instance_Name\t");
		fprintf(out_file, "Network_Name\t");


		fprintf(out_file, "Num_Nodes\t");
		fprintf(out_file, "Num_Prof_Nodes\t");
		fprintf(out_file, "Max_Tour_Duration\t");


		fprintf(out_file, "Max_Num_Discoveries\t");
		fprintf(out_file, "Max_Discoveries_Fraction\t");
		fprintf(out_file, "Best_LB\t");
		fprintf(out_file, "Best_Ub\t");
		fprintf(out_file, "Cplex_Gap\t");
		fprintf(out_file, "Root_Node_LB\t");
		fprintf(out_file, "Branch_and_Cut_Elapsed_Time\t");
		fprintf(out_file, "Starting_Cuts_Elapsed_Time\t");
		fprintf(out_file, "Tot_Time\t");
		fprintf(out_file, "Num_Logic_Benders_Cuts\t");
		fprintf(out_file, "Number_Starting_Scenarios\t");
		fprintf(out_file, "Num_BandC_Nodes\t");
		fprintf(out_file, "Number_Discoveries\t");
		fprintf(out_file, "Profit_cost\t");
		fprintf(out_file, "Discovery_cost\t");
		fprintf(out_file, "Card_Y\t");
		fprintf(out_file, "Num_Tight_Y\t");
		fprintf(out_file, "US_Type\t");
		fprintf(out_file, "US_Param\t");	
		fprintf(out_file, "US_Tot_Det_Profit\n");



		fprintf(out_file, "Table_Model:");
		fprintf(out_file, "%s\t", inst_name);
		fprintf(out_file, "%s\t", network_name);
		fprintf(out_file, "%d\t", this->num_nodes);
		fprintf(out_file, "%d\t", this->num_prof_nodes);
		fprintf(out_file, "%d\t", this->max_dur);
		fprintf(out_file, "%d\t", this->max_num_disc);
		fprintf(out_file, "%lf\t", this->max_number_discovery_fraction);
		fprintf(out_file, "%lf\t", this->bm_best_lb);
		fprintf(out_file, "%lf\t", this->bm_best_ub);
		fprintf(out_file, "%lf\t", this->cplex_gap);
		fprintf(out_file, "%lf\t", this->root_node_lb);
		fprintf(out_file, "%lf\t", this->bm_branch_and_cut_time);
		fprintf(out_file, "%lf\t", this->bm_starting_cuts_time);
		fprintf(out_file, "%lf\t", this->bm_starting_cuts_time + this->bm_branch_and_cut_time);
		fprintf(out_file, "%d\t", this->num_logic_benders_cuts);
		fprintf(out_file, "%d\t", this->num_ccg_xib);
		fprintf(out_file, "%d\t", this->bandc_nodes);
		fprintf(out_file, "%d\t", num_discovery);
		fprintf(out_file, "%lf\t", this->bm_best_ub - discovery_cost);
		fprintf(out_file, "%lf\t", discovery_cost);
		fprintf(out_file, "%d\t", this->card_Y);
		fprintf(out_file, "%d\t", this->num_tight_y);
		fprintf(out_file, "%d\t", uncertainty_set_type);
		fprintf(out_file, "%lf\t", uncertainty_parameter);
		fprintf(out_file, "%lf\n", this->tot_deterministic_profit);


	}
	return true;
}

bool EXAsolver::print_solution_latex()
{
	std::ofstream latex_out_file;

	char locname[1024];
	char ss[1024];
	snprintf(locname, sizeof(locname), "%s%s%s", results_folder, "\\Sol_Latex_Exa_", inst_name);
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
	latex_out_file << "	end={mark=triangle}," << std::endl;
	latex_out_file << "discovery_on={mark=square*}," << std::endl;
	latex_out_file << "discovery_off={mark=o}}]" << std::endl;

	latex_out_file << "\\addplot[scatter, only marks,  scatter src = explicit symbolic, forget plot] table[meta = label]{" << std::endl;
	latex_out_file << "x     y      label" << std::endl;

	/*
	Here the nodes
	*/

	/*First start and end*/
	latex_out_file << std::to_string(this->node_x_coord[start_node]) << "	" << std::to_string(this->node_y_coord[start_node]) << "	" << "start" << std::endl;
	latex_out_file << std::to_string(this->node_x_coord[end_node]) << "	" << std::to_string(this->node_y_coord[end_node]) << "	" << "end" << std::endl;
	/*Profitable nodes*/
	for (int n = first_pnode; n <= last_pnode; ++n)
	{
		if (bm_curr_sol[this->bm_v_w_i[n]] > EPS)
		{
			latex_out_file << std::to_string(this->node_x_coord[n]) << "	" << std::to_string(this->node_y_coord[n]) << "	" << "discovery_on" << std::endl;

		}
		else
		{
			latex_out_file << std::to_string(this->node_x_coord[n]) << "	" << std::to_string(this->node_y_coord[n]) << "	" << "discovery_off" << std::endl;

		}

	}


	latex_out_file << "};" << std::endl;
	latex_out_file << "\\end{axis}" << std::endl;
	latex_out_file << "\\end{tikzpicture}" << std::endl;

	latex_out_file.close();

	return true;
}

inline bool EXAsolver::mlp_update_discovery_parameters()
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
		std::cout << "Error in EXAsolver::mlp_update_discovery_parameters()" << std::endl;
		goto TERMINATE;
	}

	check = add_rows();
	if (!check)
	{
		std::cout << "Error in EXAsolver::mlp_update_discovery_parameters()" << std::endl;
		goto TERMINATE;
	}


	
	
TERMINATE:
	if (status || !check)
	{
		std::cout << "Error in EXAsolver::mlp_update_discovery_parameters()" << std::endl;
		return false;
	}
	else
		return true;
}

bool EXAsolver::rop_update_discovery_parameters()
{
	char lu[2];
	double bd[2];
	int cnt = 0;
	int status = 0;
	for (int i = this->first_pnode; i <= this->last_pnode; ++i)
	{
		if (this->w_val[i])
		{
			cnt = 0;

			this->matind[cnt] = this->v_gamma[i];
			lu[cnt] = 'L';
			bd[cnt] = -CPX_INFBOUND;
			++cnt;
			
			this->matind[cnt] = this->v_gamma[i];
			lu[cnt] = 'U';
			bd[cnt] = CPX_INFBOUND;
			++cnt;
			
			status = CPXchgbds(this->cpx_env_op, this->cpx_lp_op, cnt, this->matind.data(), lu, bd);
			if (checkCPXstatus(status, &cpx_env_op, &cpx_lp_op)) goto TERMINATE;

		}
		else
		{
			cnt = 0;
			this->matind[cnt] = this->v_gamma[i];
			lu[cnt] = 'B';
			bd[cnt] = 0;
			++cnt;

			status = CPXchgbds(this->cpx_env_op, this->cpx_lp_op, cnt, this->matind.data(), lu, bd);
			if (checkCPXstatus(status, &cpx_env_op, &cpx_lp_op)) goto TERMINATE;

		}
	}


TERMINATE:
	if (status)
		return false;
	else
		return true;
}

bool EXAsolver::bm_build_model()
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

bool EXAsolver::bm_set_cplex()
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

bool EXAsolver::bm_add_vars()
{

	unsigned int num_cols = CPXgetnumcols(this->cpx_env_bm, this->cpx_lp_bm);

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
	lb[0] = -CPX_INFBOUND;
	ub[0] = CPX_INFBOUND; // 0 profit
	status = CPXnewcols(this->cpx_env_bm, this->cpx_lp_bm, 1, obj, lb, ub, vartype, varsname);
	if (checkCPXstatus(status, &cpx_env_bm, &cpx_lp_bm)) goto TERMINATE;
	++num_cols;

	/*
	First stage variables W: discovery variables
	*/
	/*Only profitable nodes.*/
	for (int i = first_pnode; i <= last_pnode; ++i)
	{
		this->bm_v_w_i[i] = num_cols;
		obj[0] = 0;
			//this->disc_cost_n[i]; /*This can stay into objective.*/

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


		status = CPXnewcols(this->cpx_env_bm, this->cpx_lp_bm, 1, obj, lb, ub, vartype, varsname);
		if (checkCPXstatus(status, &cpx_env_bm, &cpx_lp_bm)) goto TERMINATE;

		++num_cols;
	}

TERMINATE:
	if (status)
		return false;
	else
		return true;
}

bool EXAsolver::bm_add_cnstr()
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


	if (matind.size() < CPXgetnumcols(cpx_env_bm, cpx_lp_bm))
	{
		matind.resize(CPXgetnumcols(cpx_env_bm, cpx_lp_bm));
		matval.resize(CPXgetnumcols(cpx_env_bm, cpx_lp_bm));

	}


	// Max number of discoveries
	rhs[0] = this->max_num_disc;
	sense[0] = 'L';
	nzc = 0;
	sprintf(cnstrname[0], "Max_num_discoveries");

	for (int i = this->first_pnode; i <= this->last_pnode; ++i)
	{
		matind[nzc] = this->bm_v_w_i[i];
		matval[nzc] = 1;

		++nzc;
	}
	status = CPXaddrows(cpx_env_bm, cpx_lp_bm, 0, 1, nzc, rhs, sense, matbeg, matind.data(), matval.data(), NULL, cnstrname);
	if (checkCPXstatus(status, &cpx_env_bm, &cpx_lp_bm)) goto TERMINATE;

TERMINATE:
	if (status)
		return false;
	else
		return true;
}

bool EXAsolver::set_branching_score()
{

	//CPXPARAM_MIP_Strategy_Order

	bm_branch_priority_score.assign(CPXgetnumcols(cpx_env_bm, cpx_lp_bm), 0);
	int status = 0;
	int	cnt = 0;
			for (int i = first_pnode; i <= last_pnode; ++i)
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

bool EXAsolver::bm_solve_model()
{
	bool is_ok = true;
	int status;
	int mip_status;
	CPXLONG      contextmask = 0;
	char errbuf[CPXMESSAGEBUFSIZE];
	int check_mip_stat;

	double start = 0;

	double end = 0;

	/*status = CPXwriteprob(this->cpx_env_bm, this->cpx_lp_bm, "model_bm.lp", "lp");
	if (checkCPXstatus(status, &cpx_env_bm, &cpx_lp_bm)) { is_ok = false; goto TERMINATE; };*/


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

bool EXAsolver::bm_add_starting_cuts()
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
	uint64_t curr_cut = 0;

	double start = 0;
	double end = 0;

	status = false;

	start = clock();

		/*Full discovery cut*/
		curr_cut = 0;
		for (int i = first_pnode; i <= last_pnode; ++i)
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


		this->bm_cut_lb = phi_val;
		this->bm_best_lb = phi_val;

		this->bm_add_logic_cut_to_pool(&curr_cut, &phi_val);
		curr_cut = 0;
		
		/*Now we save xi*/
		if(bm_use_scenarios_ccg_xib)
		//if (ccg_set_xib.getNumRows() > 0)
		{
			store_xbar();
		}
		else {


			nzc = 0;
			sense[0] = 'G';
			rhs[0] = phi_val;

			this->matind[nzc] = this->bm_v_phi;
			this->matval[nzc] = 1;
			++nzc;
			sprintf(cnstrname[0], "Full_discovery_cut");

			for (int i = first_pnode; i <= last_pnode; ++i)
				if (!w_val[i])
				{
					this->matind[nzc] = this->bm_v_w_i[i];
					this->matval[nzc] = (phi_val - this->bm_cut_lb);
					++nzc;

				}

			if (!this->bm_use_info_cuts)
				for (int i = first_pnode; i <= last_pnode; ++i)
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




		if (this->decision_depende_discovery) {
			curr_cut = 0;

			/*First 0 sensor*/
			for (int i = first_pnode; i <= last_pnode; ++i)
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

			/*Now we save xi*/
			//if (ccg_set_xib.getNumRows() > 0)
			if (bm_use_scenarios_ccg_xib)
			{
				store_xbar();
			}
			else{


			this->bm_add_logic_cut_to_pool(&curr_cut, &phi_val);
			nzc = 0;
			sense[0] = 'G';
			rhs[0] = phi_val;

			this->matind[nzc] = this->bm_v_phi;
			this->matval[nzc] = 1;
			++nzc;
			sprintf(cnstrname[0], "No_discovery_cut");

			for (int i = first_pnode; i <= last_pnode; ++i)
				if (!w_val[i])
				{

					this->matind[nzc] = this->bm_v_w_i[i];
					//this->matval[nzc] = -phi_val;
					this->matval[nzc] = (phi_val - this->bm_cut_lb);
					++nzc;


					//this->matind[nzc] = this->bm_v_w_i[i];
					//this->matval[nzc] = -phi_val;
					//++nzc;
				}

			if (!this->bm_use_info_cuts)
				for (int i = first_pnode; i <= last_pnode; ++i)
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


			if (false) {
				/*Cardinality Max -1*/
				for (int istar = first_pnode; istar <= last_pnode; ++istar)
				{
					curr_cut = 0;

					for (int i = first_pnode; i <= last_pnode; ++i)
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
					//if (ccg_set_xib.getNumRows() > 0)
					if (bm_use_scenarios_ccg_xib)
					{
						store_xbar();
					}

					this->bm_add_logic_cut_to_pool(&curr_cut, &phi_val);
					curr_cut = 0;


					nzc = 0;
					sense[0] = 'G';
					rhs[0] = phi_val;

					this->matind[nzc] = this->bm_v_phi;
					this->matval[nzc] = 1;
					++nzc;
					sprintf(cnstrname[0], "CrdMaxminusOne_discovery_cut_i%d", istar);

					for (int i = first_pnode; i <= last_pnode; ++i)
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
				for (int istar = first_pnode; istar <= last_pnode; ++istar)
				{
					curr_cut = 0;
					for (int i = first_pnode; i <= last_pnode; ++i)
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
					//if (ccg_set_xib.getNumRows() > 0)
					if (bm_use_scenarios_ccg_xib)
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

					for (int i = first_pnode; i <= last_pnode; ++i)
						if (!w_val[i])
						{

							this->matind[nzc] = this->bm_v_w_i[i];
							//this->matval[nzc] = -phi_val;
							this->matval[nzc] = (phi_val - this->bm_cut_lb);
							++nzc;

							/*this->matind[nzc] = this->bm_v_w_i[i];
							this->matval[nzc] = -phi_val;
							++nzc;*/
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

inline int CPXPUBLIC EXAsolver::bm_general_callback(CPXCALLBACKCONTEXTptr context, CPXLONG contextid, void* userdata)
{
	bool cut_added = false;
	bool w_integer = false;
	int status = 0;
	bool its_ok = true;
	EXAsolver* solver = (EXAsolver*)userdata;
	double val;
	double treshold;



	status = CPXcallbackgetinfodbl(context, CPXCALLBACKINFO_BEST_BND, &solver->bm_best_lb);
	if (status != 0)
	{
		solver->checkCPXstatus(status, &solver->cpx_env_bm, &solver->cpx_lp_bm);

		return status;
	}
	/*Update the bound for the cuts*/
	if (solver->bm_best_lb > solver->bm_cut_lb)
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
			its_ok = solver->bm_separate_logic_benders_cut(context, contextid, NULL, &cut_added);
		}

		/*if (!cut_added)
		{
			for (int xi = 0; xi < solver->num_ccg_xib; ++xi)
			{
				is_ok = solver->bm_separate_heuristic_subtours(context, contextid, 0.0, &cut_added, xi);
			}

		}*/
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
			w_integer = true;

			for (int i = solver->first_pnode; i <= solver->last_pnode; ++i)
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

			}


			/*cut per node if fractional.*/
		/*if (!w_integer)
		{*/
			if (node_id == solver->bm_curr_node)
			{

				solver->bm_curr_node = node_id;

				if (
					(node_depth <= 1
						&&
						/*fabs(val - solver->bm_lp_prev_iter) < (EXA_Benders_MIN_COEFF_LP_IMPR * fabs(solver->bm_lp_prev_iter))
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

			/*FRACTIONAL*/
			if (!w_integer)
			{
				if (!cut_added
					/*&&
					node_depth <= EXA_Benders_MAX_NODE_DEPTH_CUTS*/
					)
				{


						treshold = 0.25;
						for (int xi = 0; xi < solver->num_ccg_xib; ++xi)
						{
							while (!cut_added
								)
							{
								its_ok = solver->bm_separate_heuristic_subtours(context, contextid, treshold, &cut_added, xi);
								treshold += TRESHOLD_STEP;

								if (treshold > TRESHOLD_MAX)
									break;
							}

							if (!cut_added && node_depth <= 1)
							{
								its_ok = solver->bm_separate_exact_subtuors(context, contextid, &cut_added, xi);
							}

							if (cut_added)
							{
								//std::cout << "Cut Heuristics" << std::endl;
								break;
							}
						}
					
				
				}
				/*fractional ws at the end.*/
				if (!cut_added
					//&&
					//node_depth <= EXA_Benders_MAX_NODE_DEPTH_CUTS
					)
				{
					its_ok = solver->bm_separate_logic_benders_cut(context, contextid, NULL, &cut_added);
				}
		}
			else {
				its_ok = solver->bm_separate_logic_benders_cut(context, contextid, NULL, &cut_added);
			}
		}



	return status;

}

inline bool EXAsolver::bm_separate_logic_benders_cut(CPXCALLBACKCONTEXTptr context, CPXLONG contextid, double var_treshold, bool* cut_added)
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
	uint64_t curr_cut = 0;

	*cut_added = false;

		for (int i = first_pnode; i <= last_pnode; ++i)
		{
			//this->w_val[i] = (bool)std::ceil(this->bm_curr_sol[this->bm_v_w_i[i]]);

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

	if (obj > this->bm_curr_sol[this->bm_v_phi] + EXA_EPS)
	{

		nzc = 0;
		sense[0] = 'G';
		rhs[0] = obj;

		this->matind[nzc] = this->bm_v_phi;
		this->matval[nzc] = 1;
		++nzc;

		for (int i = first_pnode; i <= last_pnode; ++i)
			if(!w_val[i])
		{

				this->matind[nzc] = this->bm_v_w_i[i];
				this->matval[nzc] = (obj - this->bm_cut_lb);
				++nzc;

		
		}

		if(!this->bm_use_info_cuts)
			for (int i = first_pnode; i <= last_pnode; ++i)
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

inline void EXAsolver::bm_add_elem_to_logic_cut(uint64_t* set, unsigned int elem)
{
	*set |= (unsigned int)pow(2.00, (double)(elem - 1));
}

inline bool EXAsolver::bm_is_logic_cut_already_added(uint64_t* cut, int* index)
{
	for (int i = 0; i < bm_num_stored_logic_cuts; ++i)
		if (*cut == bm_logic_cuts_pool_set[i])
		{
			*index = i;

			return true;
		}

	return false;
}

inline void EXAsolver::bm_add_logic_cut_to_pool(uint64_t* cut, const double* phi)
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


