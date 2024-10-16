#include "K_Adapt_Problem_SPP.h"


bool K_Adapt_Problem_SPP::define_maps()
{
    unsigned short y_i = 0;

    // Map sizes.
    this->map_arc_vind.resizeMatrix2D(1 + data.N, 1 + data.N);
   // map_arc_vind.assign(1 + data.N, std::vector<unsigned short>(1 + data.N, 0));
    this->map_node_ycon.resize(this->data.N + 1);
    this->map_ycon_node.resize(this->data.N + 1);


    for(int i = 1; i <= data.N; ++i)
        for(int j = 1; j <= data.N; ++j)
            if (data.arc_exist[i][j])
            {

                this->map_arc_vind(i,j) = y_i;
                
                this->map_vind_arc.insert({ y_i, std::make_tuple(i, j) });
                
                ++y_i;
            }
    if (y_i > this->Ny)
    {
        std::cout << "Number of Y vars greater than expected in define_map_arc_Y_ind()" << std::endl;
        // Some getchar() here and there for debug mode.
        getchar();
        return false;
    }

  

    unsigned short ly = 0;
    for (int n = 1; n <= this->data.N; ++n)
    {
        this->map_node_ycon[n] = ly;
        this->map_ycon_node[ly] = n;

        if (ly >= this->Ly)
        {
            std::cout << "Number of Nodes vars greater than expected in define_maps" << std::endl;
            // Some getchar() here and there for debug mode.
            getchar();
            return false;
        }
        ++ly;       
    }
    
        return true;

}

bool K_Adapt_Problem_SPP::define_sizes()
{
    /*Num K*/
    this->K = data.K;

    this->unc_budget_B = this->data.B;

    // Number of first stage equal 0 in shortest path.
    this->Nx = 0;

    // Number of second stage equal to number of archs.
    this->Ny = data.A;

    // Number of discovery variables equal to number of archs. 
    this->Nw = data.A;
    
    // Number of xi variables and discovery variables equal to number of archs.
    this->Nxi = data.A;

    // Number of extra uncertainty param equal to 0.
    this->Nzeta = 0;

    // Number of x det constraints equal to 0.
    this->Lx = 0;

    //Number of y deterministic constraints equal to number of nodes of the graph.
    this->Ly = data.N;

    //Number of w deterministic constraints equal to 1 (max number of discovery). Set to 0 to do not have a budget.
    this->Lw = 1;

    // Number of uncertainty set constraints: 2 for each \xi var (bound [0,1]) plus one for the budget.
    this->Lxi = (2 * this->Nxi) + 1;

    // Number of other not specific variables set of deterministic constraints equal to 0.
    this->R = 0;



    return true;
}

bool K_Adapt_Problem_SPP::define_C_matrix()
{


    return false;
}

bool K_Adapt_Problem_SPP::define_D_matrix()
{
 
    return false;
}

bool K_Adapt_Problem_SPP::define_Q_matrix()
{
    // Vector size;
    //this->Q.resize(this->Nxi);

    // Matrix seixe
    this->Q.resizeMatrix2D(this->Nxi, this->Ny);

    unsigned short i, j;
    for (int row = 0; row < this->Nxi; ++row)
    {

        if (row >= Q.getNumRows())
        {
            std::cout << "Row greater than expected size of matrix Q in  bool K_Adapt_Problem_SPP::define_Q_matrix()" << std::endl;
            // Get char for debugging
            getchar();
            return false;
        }

        // Vector Size: Deprectaed
      //  this->Q[row].resize(this->Ny);

        for (int col = 0; col < this->Ny; ++col)
        {
        

            if (col >= Q.getNumCols())
            {
                std::cout << "Col greater than expected size of matrix Q in  bool K_Adapt_Problem_SPP::define_Q_matrix()" << std::endl;
                // Get char for debugging
                getchar();
                return false;
            }

            if (row == col)
            {
                i = std::get<0>(this->map_vind_arc.at(col));
                j = std::get<1>(this->map_vind_arc.at(col));
                
                this->Q(row,col) = this->data.costs[i][j] / 2.00;

                //Q[row][col] = 
            }
            else
            {
                this->Q(row, col) = 0;
            }
        }
}
    return true;
}

bool K_Adapt_Problem_SPP::define_c_vector()
{
    return false;
}

bool K_Adapt_Problem_SPP::define_d_vector()
{

    unsigned short i, j;

    this->d.resize(this->Nw);

    for (int wi = 0; wi < this->Nw; ++wi)
    {
        if (wi >= d.size())
        {

            std::cout << "Number of Y variables exceeds the size of q in  bool K_Adapt_Problem_SPP::define_q_vector()" << std::endl;
            // Some get char here and there in debug mode.
            getchar();

            return false;
        }

        i = std::get<0>(this->map_vind_arc.at(wi));
        j = std::get<0>(this->map_vind_arc.at(wi));


        this->d[wi] = 0;
            //this->data.costs[i][j];
    }
    return true;
}

bool K_Adapt_Problem_SPP::define_q_vector()
{

    this->q.resize(this->Ny);

    unsigned short i, j;
    for (int yi = 0; yi < this->Ny; ++yi)
    {
        if (yi >= q.size())
        {
            
            std::cout << "Number of Y variables exceeds the size of q in  bool K_Adapt_Problem_SPP::define_q_vector()" << std::endl;
            // Some get char here and there in debug mode.
            getchar();
            
            return false;
        }
        
        i = std::get<0>(this->map_vind_arc.at(yi));
        j = std::get<1>(this->map_vind_arc.at(yi));


        this->q[yi] = this->data.costs[i][j];
    }
    return true;
}

bool K_Adapt_Problem_SPP::define_lbq_vector()
{

    this->lbq.resize(this->Ny);

    unsigned short i, j;
    for (int yi = 0; yi < this->Ny; ++yi)
    {
        if (yi >= lbq.size())
        {

            std::cout << "Number of Y variables exceeds the size of q in  bool K_Adapt_Problem_SPP::define_q_vector()" << std::endl;
            // Some get char here and there in debug mode.
            getchar();

            return false;
        }

        i = std::get<0>(this->map_vind_arc.at(yi));
        j = std::get<1>(this->map_vind_arc.at(yi));


        this->lbq[yi] = this->data.costs[i][j];
    }
    return true;
}

bool K_Adapt_Problem_SPP::define_lbc_vector()
{
    return false;
}

bool K_Adapt_Problem_SPP::define_X_matrix()
{
    return false;
}

bool K_Adapt_Problem_SPP::define_W_matrix()
{

    this->W.resizeMatrix2D(this->Lw, this->Nw);

    unsigned short i, j;

    for (int row = 0; row < this->Lw; ++row)
    {
        if (row >= W.getNumRows())
        {
            std::cout << "Row greater than expected size of matrix W in  bool K_Adapt_Problem_SPP::define_W_matrix()" << std::endl;
            // Get char for debugging
            getchar();
            return false;
        }

        // Size of vector.
       // this->W[row].resize(this->Nw);
      //  this->W.getNumCols()
        for (int col = 0; col < this->Nw; ++col)
        {
            if (col >= W.getNumCols())
            {
                std::cout << "Col greater than expected size of matrix W in  bool K_Adapt_Problem_SPP::define_W_matrix()" << std::endl;
                // Get char for debugging
                getchar();
                return false;
            }

            this->W(row,col) = 1; // Every discovery appears in the cosntraints.
        }
    }
    return true;
    
}

bool K_Adapt_Problem_SPP::define_Y_matrix()
{
    //this->Y.resize(this->Ly);
    this->Y.resizeMatrix2D(this->Ly, this->Ny);

    unsigned short i, j, n;

    for (int row = 0; row < this->Ly; ++row)
    {
        if (row >= Y.getNumRows())
        {
            std::cout << "Row greater than expected size of matrix Y in  bool K_Adapt_Problem_SPP::define_Y_matrix()" << std::endl;
            // Get char for debugging
            getchar();
            return false;
        }
        // Vector size
       // this->Y[row].resize(this->Ny);

        for (int col = 0; col < this->Ny; ++col)
        {

            if (col >=this->Y.getNumCols())
            {
                std::cout << "Col greater than expected size of matrix Y in  bool K_Adapt_Problem_SPP::define_Y_matrix()" << std::endl;
                // Get char for debugging
                getchar();
                return false;
            }

            i = std::get<0>(this->map_vind_arc.at(col));
            j = std::get<1>(this->map_vind_arc.at(col));
            n = this->map_ycon_node.at(row);

            if (j == n)
            {
                this->Y(row, col) = 1;
               // this->Y[row][col] = 1;
            }
            else
            {
                if (i == n)
                {
                    this->Y(row, col) = -1;
                 //   this->Y[row][col] = -1;
                }
                else
                {
                    this->Y(row, col) = 0;
                   // this->Y[row][col] = 0;
                }
            }
        }
}
    return true;


}

bool K_Adapt_Problem_SPP::define_bx_vector()
{
    return false;
}

bool K_Adapt_Problem_SPP::define_by_vector()
{
    // Vector size
    this->b_y.resize(this->Ly);
    
    unsigned short n;
    for (int ly = 0; ly < this->Ly; ++ly)
    {
        n = this->map_ycon_node[ly];

        if (n == this->data.src ||
            n == this->data.tgt)
        {
            if(n == this->data.src)
             this->b_y[ly] = -1;
            else
             this->b_y[ly] = 1;
        }
        else
        {
            this->b_y[ly] = 0;
        }
    }

    return false;
}

bool K_Adapt_Problem_SPP::define_bw_vector()
{

    // Vector size
    this->b_w.resize(this->Lw);

    for (int lw = 0; lw < this->Lw; ++lw)
    {
        if (lw >= this->b_w.size())
        {

            std::cout << "Number of Y variables exceeds the size of q in  bool K_Adapt_Problem_SPP::define_q_vector()" << std::endl;
            // Some get char here and there in debug mode.
            getchar();

            return false;
        }

        //this->b_w[lw] = this->data.N;  //Number of discovery equal to the number of nodes.

         //this->b_w[lw] = ?  //Number of discovery equal to ? For the DDID version.
        if (ddid_version)
        {
            this->b_w[lw] = ceil((double)this->data.A * this->data.Delta);
           // this->b_w[lw] = (int)(this->data.A * DDID_FACTOR);
            //this->b_w[lw] = this->data.N;  //Max Number of discovery equal to the number of nodes in the shortest path.
        }
        else {
            this->b_w[lw] = this->Nw;  //Max Number of discovery equal to the number of w in K-Adapt.

        }


    }
    return true;

}



bool K_Adapt_Problem_SPP::define_A_matrix()
{
    unsigned short row, col;

    this->A.resizeMatrix2D(this->Lxi, this->Nxi);
   // this->A.resize(this->Lxi);

    // First the upper bound xi <= 1: therefore we stop at Nxi;
    for (row = 0; row < this->Nxi; ++row)
    {
        // Vector size;
        //this->A[row].resize(this->Nxi);

        for (col = 0; col < this->Nxi; ++col)
        {
            if (row == col)
            {
                this->A(row, col) = 1;
                //this->A[row][col] = 1;
            }
        }
    }
    // Now the lb, therefore we start from the second part of the matrix.
    for (row = this->Nxi; row < this->Lxi - 1; ++row)
    {
        // Vector size;
       // this->A[row].resize(this->Nxi);

        for (col = 0; col < this->Nxi; ++col)
        {
            if ((row - this->Nxi) == col)
            {
                this->A(row, col) = -1;
        
            }
        }
    }
    // Now the last constraints: the budget
    row = this->Lxi - 1;

    // Vector size;
   // this->A[row].resize(this->Nxi);
    for (col = 0; col < this->Nxi; ++col)
    {
        this->A(row, col) = 1;
        // this->A[row][col] = 1;
    }
    return true;
}

bool K_Adapt_Problem_SPP::define_b_vector()
{
    // Vector size
    this->b.resize(this->Lxi);


    // Same spirity of the A matrix
    unsigned short row;
    // First upper bound. Therefore 1.
    for (row = 0; row < this->Nxi; ++row)
    {
        this->b[row] = 1;
    }
    // Second lower bound. Therefore 0.
    for (row = this->Nxi; row < this->Lxi - 1; ++row)
    {
        this->b[row] = 0;
    }
    // Finally uncertainty budget.
    row = this->Lxi - 1;
    this->b[row] = this->unc_budget_B;
    
    return true;
}

bool K_Adapt_Problem_SPP::define_H_matrix()
{
    H.assignMatrix2D(this->L, this->Nxi, 0);
    return true;
}

bool K_Adapt_Problem_SPP::define_G_matrix()
{

    return false;
}

bool K_Adapt_Problem_SPP::define_Ex_matrix()
{
    return true;
}

bool K_Adapt_Problem_SPP::define_Ew_matrix()
{
    return true;
}

bool K_Adapt_Problem_SPP::define_Ey_matrix()
{
    return true;
}

bool K_Adapt_Problem_SPP::define_g_vector()
{
    return true;
}

bool K_Adapt_Problem_SPP::define_T_matrix()
{
    return true;
}

bool K_Adapt_Problem_SPP::define_P_matrix()
{
    return true;
}

bool K_Adapt_Problem_SPP::define_V_matrix()
{
    return true;
}

bool K_Adapt_Problem_SPP::define_Txi_matrix()
{
    return true;
}

bool K_Adapt_Problem_SPP::define_Vxi_matrix()
{
    return true;
}

bool K_Adapt_Problem_SPP::define_Pxi_matrix()
{
    return true;
}

bool K_Adapt_Problem_SPP::define_h_vector()
{
    return true;
}

bool K_Adapt_Problem_SPP::define_lb_xi()
{
    return true;
}

bool K_Adapt_Problem_SPP::define_ub_xi()
{
    return true;
}

bool K_Adapt_Problem_SPP::define_hmax_vector()
{
    return false;
}

bool K_Adapt_Problem_SPP::define_v_vector()
{
    return false;
}

bool K_Adapt_Problem_SPP::define_p_vector()
{
    return false;
}

bool K_Adapt_Problem_SPP::define_t_vector()
{
    return false;
}



bool K_Adapt_Problem_SPP::setInstanceData(const SPP& data)
{
    this->data = data;
    this->unc_budget_B = data.B;
    this->K = data.K;
    return true;
}

SPP K_Adapt_Problem_SPP::getInstanceData()
{
    return this->data;
}

bool K_Adapt_Problem_SPP::print_out_file(RobustKSolver* solver)
{
    std::ofstream outfile;

    outfile.open(this->data.solfilename);

    outfile << "Instance_Name: " << data.solfilename << std::endl;
    outfile << "K: " << data.K << std::endl;
    outfile << "U_Budget: " << data.B << std::endl;
    outfile << "Num_Nodes: " << data.N<< std::endl;
    outfile << "Num_Archs: " << data.A << std::endl;
    outfile << "Max_Num_Discoveries: " << this->b_w[0] << std::endl;
    outfile << "DDDI_Version: " << ddid_version << std::endl;
    outfile << "BUB: " << solver->bub << std::endl;
    outfile << "BLB: " << solver->blb << std::endl;
    outfile << "Gap: " << solver->cplex_gap << std::endl;
    outfile << "Elapsed_Time " << solver->comp_time << std::endl;

    //fprintf(out_file, "\n\nSolution:\n\n");
    outfile << "\n\nSolution:\n\n";
    outfile << "w[] variables\n:";
    for (int i = 0; i < this->Nw; ++i)
       if(solver->getCurrSol()[solver->get_v_w_ind()[i]] > INT_EPS)
        // if (solver->master_curr_sol[solver->getIndW()[i]] > mEPS)
        {
           outfile << "w_" << i << ": arc (" << std::get<0>(this->map_vind_arc[i]) << ", " << std::get<1>(this->map_vind_arc[i]) << ") = " << solver->getCurrSol()[solver->get_v_w_ind()[i]]<<std::endl;
          //  fprintf(out_file, "w_%d: arc (%d, %d) = %lf\n", i, std::get<0>(this->map_vind_arc[i]), std::get<1>(this->map_vind_arc[i]), solver->getBM_curr_sol()[solver->getIndW()[i]]);

        }
    outfile << "\n**********************\n";
    
    outfile << "y_k[] variables:\n";
    for(int k = 0; k < this->data.K; ++k)
    for (int i = 0; i < this->Ny; ++i)
        if (solver->getCurrSol()[solver->get_v_y_k_ind()(k,i)] > INT_EPS)
            // if (solver->master_curr_sol[solver->getIndW()[i]] > mEPS)
        {
            outfile << "y_k" << k <<"_i"<<i<< ": arc (" << std::get<0>(this->map_vind_arc[i]) << ", " << std::get<1>(this->map_vind_arc[i]) << ") = " << solver->getCurrSol()[solver->get_v_y_k_ind()(k,i)]<<std::endl;
            //  fprintf(out_file, "w_%d: arc (%d, %d) = %lf\n", i, std::get<0>(this->map_vind_arc[i]), std::get<1>(this->map_vind_arc[i]), solver->getBM_curr_sol()[solver->getIndW()[i]]);

        }
    outfile << "\n\n------------------------\n\n";


    outfile << "Legenda:Instance_Name   K   Num_Nodes   Num_Archs Seed  U_budget    Max_Num_Discoveries    DDDI_Version Delta    BUB BLB Gap Elapsed_Time" << std::endl;
    outfile << "Table:"<<data.solfilename<<"   "<<K<<"  "<<this->data.N<<" "<<this->data.A<<" "<<this->data.seed <<" "<<this->data.B<<"    "<< this->b_w[0]<<"  "<<
        ddid_version<<" "<<this->data.Delta<<" "<< solver->bub << " " << solver->blb <<
        "   "<<solver->cplex_gap<<" " <<solver->comp_time<<" "<< std::endl;




    outfile.close();



    return false;
}



