#pragma once
#include "K_Adapt_Problem.h"
#include "RobustKSolver.h"
#include "instance_spp.h"


class K_Adapt_Problem_SPP : public K_Adapt_Problem 
{
private:
	/** Specific instance data */
	SPP data;
	/*Maps arc with  index of y and w variables. .*/
	Matrix2D<double> map_arc_vind;
	/*Maps index-tple of y and w variables*/
	std::map<unsigned short, std::tuple<unsigned short, unsigned short>> map_vind_arc;
	/*Map nodes -> constraints*/
	std::vector<unsigned short> map_node_ycon;
	/*Map con -> node*/
	std::vector<unsigned short> map_ycon_node;

	unsigned short unc_budget_B;



	
	
	bool define_maps() override;

	bool define_sizes() override;

	/* Virtual Methods for  Matrices */
	bool define_C_matrix() override;
	bool define_D_matrix() override;
	bool define_Q_matrix() override;

	bool define_c_vector() override; //Will be 0
	bool define_d_vector() override; 
	bool define_q_vector() override; 

	bool define_lbq_vector() override;
	bool define_lbc_vector() override;


	bool define_X_matrix() override;
	bool define_W_matrix() override; 
	bool define_Y_matrix() override; 

	bool define_bx_vector() override; // Do not need this
	bool define_by_vector()	override; 
	bool define_bw_vector() override;


	bool define_A_matrix() override;
	bool define_b_vector() override;
	bool define_G_matrix() override; // All 0


	bool define_Ex_matrix() override;
	bool define_Ew_matrix() override; // Not necessary for this problem.
	bool define_Ey_matrix() override;
	bool define_g_vector() override;

	bool define_T_matrix() override;
	bool define_P_matrix() override;
	bool define_V_matrix() override;
	bool define_Txi_matrix() override;
	bool define_Vxi_matrix() override;
	bool define_Pxi_matrix()  override;
	bool define_H_matrix() override;
	bool define_h_vector() override;


	bool define_lb_xi() override;
	bool define_ub_xi() override;

	bool define_hmax_vector() override;
	bool define_v_vector() override;
	bool define_p_vector() override;
	bool define_t_vector() override;





public: 

/**
 * Set data of the instance that this object refers to
 * 
 @param data problem instance
 */
	
	bool setInstanceData(const SPP& data);

	SPP getInstanceData();

	bool print_out_file(RobustKSolver* solver);
};

