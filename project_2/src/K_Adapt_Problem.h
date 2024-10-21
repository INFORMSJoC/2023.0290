#pragma once
#include"Constant.h"
#include "Utilities.h"
#include "instance_knp.h"



#define CNSTR_UNC_BIG_M 1000
extern bool ddid_version; // if True, problem with Decision Dependent Information Discovery, otherwise, classical K-Adaptability problem.


/*
* Abstract class to represent the instance informartion 
  for a genereric TWO-stage RO problem.
*/
class K_Adapt_Problem
{
protected:

	std::string problem_name;
	unsigned short K;
	unsigned short Nx, Nw, Ny, Nxi, Nzeta; // Variables.
	unsigned short Lx, Ly, Lw, Lxi, R; // Number of cosntraints for X, Y, W, Uncertainty set and deterministic constraints linking X, Y, and W.
	unsigned short L; //Number of uncertanty constraints.
	
					  /*
	Deterministic Objective Matrices
	c, d, q
	*/
	 std::vector<double> c; // deterministic vector cost for x.
	 std::vector<double> d; // deterministic vector cost for w.
	 std::vector<double> q; // deterministic vector cost for y.

	/*
	Uncertainty Objective Matrices
	*/
	Matrix2D<double> C; 
	Matrix2D<double> D; // deterministic vector cost for w.
	Matrix2D<double> Q; // deterministic vector cost for y.
	
	/*
	
	*/



	/*
	* Min-Cost uncertainty vector: for Valid inequalities.
	*/
	std::vector<double> lbc; // deterministic vector cost for x.
	std::vector<double> lbd; // deterministic vector cost for w.
	std::vector<double> lbq; // deterministic vector cost for y.

	
	

	/*
	* Deterministic constraints matrices
	*/
	Matrix2D<double> X; // First stage decision variables (x) constraints matrix.
	std::vector<double> b_x; // Rhs for X.
		
	Matrix2D<double> W; // First stage decision .
	std::vector<double> b_w; // Rhs for W.
		
	Matrix2D<double> Y; // First discovery decision matrix.
	std::vector<double> b_y; // Rhs for W.

	/*
	Deterministic constraints linking first and second stage variables
	*/
	Matrix2D<double> Ex; // For x variables
	
	Matrix2D<double> Ew; // For w variables.
	
	Matrix2D<double> Ey; // For y variables
	std::vector<double> g; // Rhs general linking constrints.

	/*
	Uncertain constraints matrices
	*/
	Matrix2D<double> T; // For x variables

	Matrix3D<double> Txi; // Constraints matrix for x and xi variables.

	Matrix2D<double> V; // For w variables.

	Matrix3D<double> Vxi; // Constraints matrix for w and xi variables.

	Matrix2D<double> P; // For y variables

	Matrix3D<double> Pxi; // Constraints matrix for y and xi variables.

	Matrix2D<double> H; // For \xi  variables: All 0 for this problem.

	std::vector<double> h; // Rhs for uncertainty set.

	/*Relaxed uncertainty constraints.*/
	std::vector<double> t; // min coeff X variables.
	std::vector<double> v; // min coeff discovery variables.
	std::vector<double> p; // min coeff second stage variables. 
	std::vector<double> hmax; // max rhs 



		/*
	Uncertainty set
		*/
	Matrix2D<double> A; // Uncertainty set for xi vars.
	Matrix2D<double> G; // Uncertainty set for zeta vars size L_xi times Nzeta.
	std::vector<double> b; // Rhs for uncertainty set.


	/*
	Variables bounds
	*/
	std::vector<double> ub_x; // UB x
	std::vector<double> lb_x; // LB x

	std::vector<double> ub_w; // UB w
	std::vector<double> lb_w; // UB w

	std::vector<double> ub_y; // UB y
	std::vector<double> lb_y; // UB y

	std::vector<double> ub_xi; // UB xi
	std::vector<double> lb_xi; // UB xi

	std::vector<double> ub_zeta; // UB zeta
	std::vector<double> lb_zeta; // UB zeta

	
	/* Define indices lenght */
	virtual bool define_sizes() = 0;
	virtual bool define_maps() = 0;


	/* Virtual Methods for  Matrices */
	virtual bool define_C_matrix() = 0;
	virtual bool define_D_matrix() = 0;
	virtual bool define_Q_matrix() = 0;

	virtual bool define_c_vector() = 0;
	virtual bool define_d_vector() = 0;
	virtual bool define_q_vector() = 0;

	virtual bool define_lbc_vector() = 0;

	virtual bool define_lbq_vector() = 0;

	virtual bool define_X_matrix() = 0;
	virtual bool define_W_matrix() = 0;
	virtual bool define_Y_matrix() = 0;

	virtual bool define_bx_vector() = 0;
	virtual bool define_by_vector() = 0;
	virtual bool define_bw_vector() = 0;

	bool define_uncertainty_set();
	
	virtual bool define_A_matrix() = 0;
	virtual bool define_G_matrix() = 0;
	virtual bool define_b_vector() = 0;

	virtual bool define_Ex_matrix() = 0;
	virtual bool define_Ey_matrix() = 0;
	virtual bool define_Ew_matrix() = 0;
	virtual bool define_g_vector()  = 0;

	virtual bool define_T_matrix() = 0;
	virtual bool define_V_matrix() = 0;
	virtual bool define_P_matrix() = 0;
	virtual bool define_Txi_matrix() = 0;
	virtual bool define_Vxi_matrix() = 0;
	virtual bool define_Pxi_matrix() = 0;
	virtual bool define_H_matrix() = 0;
	virtual bool define_h_vector() = 0;

	virtual bool define_t_vector() = 0; // Vector that approximate the unc. constraint for First Stage vars.
	virtual bool define_v_vector() = 0; // Vector that approximate the unc. constraint for Discovery Stage vars.
	virtual bool define_p_vector() = 0; // Vector that approximate the unc. constraint for Second stage Stage vars.
	virtual bool define_hmax_vector() = 0; // RHS to approximate uncertainty constraints.


	virtual bool define_lb_xi() = 0;
	virtual bool define_ub_xi() = 0;




public:


	K_Adapt_Problem() = default; // Default constructor.

	~K_Adapt_Problem() = default; // Default distructor.

	bool setInstance();
	

	/*
	* get the number of first stage variables (Nx)
	*/
	unsigned short getNx();
	unsigned short getNy();
	unsigned short getNw();
	unsigned short getNxi();
	unsigned short getNzeta();
	unsigned short getK();
	unsigned short getL();
	unsigned short getLxi();
	unsigned short getLx();
	unsigned short getLy();
	unsigned short getLw();
	unsigned short getR(); 

	//***** Obejective matrices for primal variables
	// First stage uncertainty objective Matrix
	Matrix2D<double> get1S_uobj_matr_C();

	// Discovery variables uncertainty objective matrix
	Matrix2D<double> getDV_uobj_matr_D();

	// Second stage variables uncertainty objective matrix
	Matrix2D<double> get2S_uobj_matr_Q();

	// First stage variables deterministic objective vector.
	std::vector<double> get1S_dobj_vect_c();

	// First stage variables lb objective vector.
	std::vector<double> get1S_dobj_vect_lbc();

	// Discovery variables deterministic objective vector. 
	std::vector<double> getDV_dobj_vect_d();

	// Second stage variables deterministic objective vector.
	std::vector<double> get2S_dobj_vect_q();

	// Second stage variables lowerbound on the coefficient for each xi.
	std::vector<double> get2S_dobj_vect_lbq();
 
	//****** Constraints matrices for primal variables
	// First stage deterministic constraints Matrix
	Matrix2D<double> get1S_dcon_matr_X();

	// RHS of first stage decision constraints.
	std::vector<double> get1S_drhs_vect_bx();

	// Discovery variables deterministic constraints matrix
	Matrix2D<double> getDV_dcon_matr_W();
	
	// RHS of first stage decision constraints.
	std::vector<double> getDV_drhs_vect_bw();

	// Second stage matrix.
	Matrix2D<double> get2S_dcon_matr_Y();

	// RHS of second stage decision constraints.
	std::vector<double> get2S_drhs_vect_by();


	// Uncertainty set matrices

	// xi Matrix A
	Matrix2D<double> getUC_xicon_matr_A();

	// zeta Matrix G
	Matrix2D<double> getUC_zetacon_matr_G();

	// rhs vector for uncertainty set 
	std::vector<double> getUC_rhs_vect_b();


	/*Ex, Ey, Ew, g*/
	Matrix2D<double> getLink_dcon_matr_Ex();
	Matrix2D<double> getLink_dcon_matr_Ew();
	Matrix2D<double> getLink_dcon_matr_Ey();
	std::vector<double> getLink_dcon_vect_g();

	/*T, V, P, H, h		*/	
	Matrix2D<double> get1S_ucon_matr_T();
	Matrix2D<double> getDV_ucon_matr_V();
	Matrix2D<double> get2S_ucon_matr_P();
	Matrix2D<double> getUC_ucon_matr_H();
	std::vector<double> get_rhs_ucon_vect_h();
	
	/*Txi, Vxi, Pxi, 	*/
	Matrix3D<double> get1S_ucon_matr_Txi();
	Matrix3D<double> getDV_ucon_matr_Vxi();
	Matrix3D<double> get2S_ucon_matr_Pxi();

	/*t, v, p, hmax*/
	std::vector<double> get1S_ducon_vect_t();
	std::vector<double> getDV_ducon_vect_v();
	std::vector<double> get2S_ducon_vect_p();
	std::vector<double> get_rhs_ducon_vect_hmax();


	/*Bounds variables*/
	std::vector<double> get_ub_x();
	std::vector<double> get_ub_w();
	std::vector<double> get_ub_y();
	std::vector<double> get_ub_xi();
	std::vector<double> get_ub_zeta();

	std::vector<double> get_lb_x();
	std::vector<double> get_lb_w();
	std::vector<double> get_lb_y();
	std::vector<double> get_lb_xi();
	std::vector<double> get_lb_zeta();
	
	
	std::string getProbName();

};




