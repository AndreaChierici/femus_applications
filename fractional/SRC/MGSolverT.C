#include "Equations_conf.h"

// ============================================
#ifdef T_EQUATIONS  // 3D-2D Energy equation
// ============================================

// class local configuration -------
#include "MGSolverT.h"

// configuration files -----------
#include "Printinfo_conf.h"

// standard  library
#include <sstream>
#include "MGFE.h"
#include "MGGeomEl.h"
#include "MGUtils.h"
#include "MeshExtended.h"
#include "dense_matrixM.h"
#include "dense_vectorM.h"
#include "linear_solverM.h"
#include "numeric_vectorM.h"
#include "parallelM.h"
#include "sparse_matrixM.h"

// ======================================================
/// This function constructs the 3d-2D MGSolT class
// ==========================================================================
/*! This constructor needs    MGEquationsSystem &mg_equations_map_in object to be constructed.
 * This equation has 1 quadratic variable (T) defined in nvars_in[]=(0,0,1),
 * equation name "T", basic variable name "T"
 */
MGSolT::MGSolT(
    MGEquationsSystem& mg_equations_map_in,  ///<  mg_equations_map_in pointer
    const int nvars_in[],                    ///< KLQ number of variables
    std::string eqname_in,                   ///< equation name
    std::string varname_in                   ///< basic variable name
    )
    : MGSolDA(mg_equations_map_in, nvars_in, eqname_in, varname_in),
      _offset(_mgmesh._NoNodes[_NoLevels - 1]),  // mesh nodes
      _dt(stod(_mgutils._sim_config["dt"])),     // parameter  dt
      _uref(_mgutils._mat_prop["Uref"]),         // parameter  u reference
      _lref(_mgutils._mat_prop["Lref"]),         // parameter  l reference
      _rhof(_mgutils._mat_prop["rho0"]),         // parameter density
      _muf(_mgutils._mat_prop["mu0"]),           // parameter viscosity
      _Tref(_mgutils._mat_prop["Tref"]),         // parameter  temperature reference
      _cp0(_mgutils._mat_prop["cp0"]),           // parameter  Cp reference
      _kappa0(_mgutils._mat_prop["kappa0"]) {    // parameter  conductivity reference
  //  =========================================================================

  // READ PARAMETERS FROM CLASS T_parameter
  _T_parameter.read_param(_mgutils);

  /// A) reading parameters  for field coupling (in _FF_idx[])
  _nTdim = DIMENSION;
  for (int k_index = 0; k_index < 30; k_index++) { _FF_idx[k_index] = -1; }
  /// B) setting class variable name T (in _var_names[0]) and ref value T_ref (in _refvalue[0])
  _var_names[0] = varname_in;
  _refvalue[0] = _Tref;

  /// C ) Setting the  solver type (with _solver[l]->set_solver_type(SOLVERT))
  for (int l = 0; l < _NoLevels; l++) { _solver[l]->set_solver_type(_T_parameter._SolverType); }

  /// D) Setting nondimensional parameters _alpha _IPrdl _IRe .....
  _alpha = _kappa0 / (_rhof * _cp0);
  _IPrdl = _rhof * _alpha / _muf;
  _IRe = _muf / (_rhof * _uref * _lref);

  _IPrdl_turb = 1. / _T_parameter._Prt;
  _alpha_turb = 0.;

  _qheat = _mgutils._mat_prop["qheat"] * _lref / (_rhof * _cp0 * _Tref * _uref);
  _qs = _mgutils._mat_prop["qs"] / (_rhof * _cp0 * _Tref * _uref);
  _Wall_dist = _mgutils._geometry["Wall_dist"];

  _SolveT = (_mgutils._sim_config["SolveTemperature"].compare("yes") == 0) ? true : false;

  _Axisym = _mgutils._geometry["Axisym"];

  _NumRestartSol = _T_parameter._NumRestartSol;

  return;
}

//  ===============================================================================================
/// This function assembles the matrix and the rhs:
//  ===============================================================================================
void MGSolT::GenMatRhs(
    const double /**< time (in) */, const int Level /**< discretization Level (in) */,
    const int mode /**< y/n assemble rhs  (1=rhs) (0=only matrix) (in)*/
) {                // ===============================================
                   //   double Crank_Nicolson =1.;
  /// a) Set up

  // geometry ---------------------------------------------------------------------------------------
  const int offset = _mgmesh._NoNodes[_NoLevels - 1];  // mesh nodes
  const int el_sides = _mgmesh._GeomEl._n_sides[0];    // element sides
  int el_conn[NDOF_FEM];                               // element connectivity
  int el_neigh[NDOF_FEM];                              // bd element connectivity
  int el_conn_j[NDOF_FEM];                             // element connectivity (loop on j)
  int el_neigh_j[NDOF_FEM];                            // bd element connectivity (loop on j)
  int sur_toply[NDOF_FEMB];                            // boundary topology

  // gauss integration  -----------------------------------------------------------------------------
  double x_m[DIMENSION];
  double normal[DIMENSION];
  const int el_ngauss = _fe[2]->_NoGauss1[_nTdim - 1];   // elem gauss points
  const int elb_ngauss = _fe[2]->_NoGauss1[_nTdim - 2];  // bd elem gauss points
  const int el_ngauss_j = _fe[1]->_NoGauss1[_nTdim - 1]; // put as a variable

  double WallDist[NDOF_FEM], AlphaTurb[NDOF_FEM];

  // element dofs (costant[0]-linear[1]-quadratic[2]) -----------------------------------------------
  int el_ndof[3];
  el_ndof[0] = 1;
  int elb_ndof[3];
  elb_ndof[0] = 1;       // number of el dofs
  int el_mat_nrows = 0;  // number of mat rows (dofs)
  for (int ideg = 1; ideg < 3; ideg++) {
    el_ndof[ideg] = _fe[ideg]->_NoShape[_nTdim - 1];
    elb_ndof[ideg] = _fe[ideg]->_NoShape[_nTdim - 2];
    el_mat_nrows += _nvars[ideg] * el_ndof[ideg];
  };
  const int el_ndof2 = _fe[2]->_NoShape[_nTdim - 1];

  int el_mat_ncols = el_mat_nrows;                // square matrix
  std::vector<int> el_dof_indices(el_mat_ncols);  // element dof vector
  std::vector<int> el_dof_indices_j(el_mat_ncols);  // element dof vector

  // coupling  fields -------------------------------------------------------------------------------
  for (int k = 0; k < 30; k++) {  // coupling  basic system fields
    const int idx = _data_eq[2].tab_eqs[k];
    _FF_idx[k] = (idx >= 0) ? _data_eq[2].indx_ub[idx] : -1;
  }
  double vel_g[DIMENSION];
  for (int idim = 0; idim < _nTdim; idim++) {
    vel_g[idim] = 0.;  // velocity not coupled
  }

  // element matrix and rhs  (mode 0= matrix only or mode 1=matrix +rhs) ---------------------------
  A[Level]->zero();
  if (mode == 1) {
    b[Level]->zero();  // global matrix+rhs
  }

  _KeM.resize(el_mat_nrows, el_mat_ncols);
  _FeM.resize(el_mat_nrows);  // resize  local  matrix+rhs
  
  _CClocalII.resize(el_mat_nrows, el_mat_ncols);
  _CClocalIJ.resize(el_mat_nrows, el_mat_ncols);
  _CClocalJI.resize(el_mat_nrows, el_mat_ncols);
  _CClocalJJ.resize(el_mat_nrows, el_mat_ncols);
  _Res_nonlocalI.resize(el_mat_nrows);
  _Res_nonlocalJ.resize(el_mat_nrows);
  
  _CClocal_refined.resize(el_mat_nrows, el_mat_ncols);
  _Res_local_refined.resize(el_mat_nrows);

  int ndof_lev = 0;
  for (int pr = 0; pr < _mgmesh._iproc; pr++) {
    int delta = _mgmesh._off_el[0][pr * _NoLevels + Level + 1] - _mgmesh._off_el[0][pr * _NoLevels + Level];
    ndof_lev += delta;
  }

  
// //   Loop over jel, for non-local models
  for(unsigned kproc = 0; kproc < _mgmesh._n_subdom; kproc++) {
  const int nel_e_j = _mgmesh._off_el[0][Level + _NoLevels * kproc + 1];  // start element
  const int nel_b_j = _mgmesh._off_el[0][Level + _NoLevels * kproc];      // stop element
    for (unsigned jel = 0; jel < (nel_e_j - nel_b_j); jel++) {
        
      if(_iproc == kproc){
        /// 1. geometry and element  fields ------------------------------------
        // Element Connectivity (el_conn_j)  and coordinates (_xx_qnds_j)
        _mgmesh.get_el_nod_conn(0, Level, jel, el_conn_j, _xx_qnds_j);
        _mgmesh.get_el_neighbor(el_sides, 0, Level, jel, el_neigh_j);
        // set element-nodes variables  bc (bc_q_dofs)
        get_el_dof_bc(Level, jel + ndof_lev, el_ndof, el_conn_j, offset, el_dof_indices_j, _bc_vol_j, _bc_bd_j);
        // fill the node data vectors
        for (int deg = 0; deg < 3; deg++) {
            for (int eq = 0; eq < _data_eq[deg].n_eqs; eq++) {
                _data_eq[deg].mg_eqs[eq]->get_el_sol(
                0, 0, _data_eq[deg].indx_ub[eq + 1] - _data_eq[deg].indx_ub[eq], el_ndof[deg], el_conn_j, offset,
                _data_eq[deg].indx_ub[eq], _data_eq[deg].ub);
            }
        }
     _data_eq[2].mg_eqs[_data_eq[2].tab_eqs[T_F]]->get_el_sol(
        0, 0, 1, el_ndof[2], el_conn_j, offset, 0, _T_jel);  // time step -1

  
      }
      
      MPI_Bcast(& el_conn_j[0], NDOF_FEM, MPI_INT, kproc, MPI_COMM_WORLD);
      MPI_Bcast(& _xx_qnds_j[0], NDOF_FEM * DIMENSION, MPI_DOUBLE, kproc, MPI_COMM_WORLD);
      MPI_Bcast(& _bc_vol_j[0], NDOF_FEM, MPI_INT, kproc, MPI_COMM_WORLD);
      MPI_Bcast(& _bc_bd_j[0], NDOF_FEM, MPI_INT, kproc, MPI_COMM_WORLD);
      MPI_Bcast(& _T_jel[0], NDOF_FEM, MPI_DOUBLE, kproc, MPI_COMM_WORLD);
      
      
      vector < vector < double > > xg2(el_ngauss_j);
      vector <double> weight2(el_ngauss_j);
      vector < vector <double> > phi2(el_ngauss_j);  // local test function
      std::vector< double > solY(el_ngauss_j, 0.);

      for(unsigned jg = 0; jg < el_ngauss_j; jg++) {
        double det_j = _fe[2]->Jac(jg, _xx_qnds_j, _InvJacJ);          // Jacobian
        double JxW_gJ = det_j * _fe[2]->_weight1[_nTdim - 1][jg];    // weight
        _fe[2]->get_phi_gl_g(_nTdim, jg, _phi_g[2]);                // shape funct

        xg2[jg].assign(DIMENSION, 0.);
        solY[jg] = 0.;

        for(unsigned j = 0; j < el_ndof2; j++) {
          solY[jg] += _T_jel[j] * _phi_g[2][j];
          for(unsigned k = 0; k < DIMENSION; k++) {
            xg2[jg][k] += _xx_qnds_j[k * NDOF_FEM + j] * _phi_g[2][j];
          }
        }
      }
      

      
  
  /// b) Element  Loop over the volume (n_elem)
  const int nel_e = _mgmesh._off_el[0][Level + _NoLevels * _iproc + 1];  // start element
  const int nel_b = _mgmesh._off_el[0][Level + _NoLevels * _iproc];      // stop element
  for (unsigned iel = 0; iel < (nel_e - nel_b); iel++) {
    // set to zero matrix and rhs and center
    if(iel == jel){
      _KeM.zero();
      _FeM.zero();
      _CClocal_refined.zero();
      _Res_local_refined.zero();    //resize
    }
    _CClocalII.zero();
    _CClocalIJ.zero();
    _CClocalJI.zero();
    _CClocalJJ.zero();
    _Res_nonlocalI.zero();
    _Res_nonlocalJ.zero();

    /// 1. geometry and element  fields ------------------------------------
    // Element Connectivity (el_conn)  and coordinates (_xx_qnds)
    _mgmesh.get_el_nod_conn(0, Level, iel, el_conn, _xx_qnds);
    _mgmesh.get_el_neighbor(el_sides, 0, Level, iel, el_neigh);

    // set element-nodes variables  bc (bc_q_dofs)
    get_el_dof_bc(Level, iel + ndof_lev, el_ndof, el_conn, offset, el_dof_indices, _bc_vol, _bc_bd);
    // fill the node data vectors
    for (int deg = 0; deg < 3; deg++) {
      for (int eq = 0; eq < _data_eq[deg].n_eqs; eq++) {
        _data_eq[deg].mg_eqs[eq]->get_el_sol(
            0, 0, _data_eq[deg].indx_ub[eq + 1] - _data_eq[deg].indx_ub[eq], el_ndof[deg], el_conn, offset,
            _data_eq[deg].indx_ub[eq], _data_eq[deg].ub);
      }
    }

    _data_eq[2].mg_eqs[_data_eq[2].tab_eqs[T_F]]->get_el_sol(
        0, 0, 1, el_ndof[2], el_conn, offset, 0, _T_1ts);  // time step -1
    _data_eq[2].mg_eqs[_data_eq[2].tab_eqs[T_F]]->get_el_sol(
        1, 0, 1, el_ndof[2], el_conn, offset, 0, _T_2ts);  // time step -2

    // ----------------------------------------------------------------------------------
    /// 2. Boundary integration  (bc)
    // ----------------------------------------------------------------------------------

    if(iel == jel){
    for (int k = 0; k < el_ndof[2]; k++) {
      _bc_el[k] = (_bc_vol[k] / 10 == 0) ? 0 : 1;  // boundary condition
      _bc_bd[k] = 1;
    }
    for (int idim = 0; idim < _nTdim; idim++) {
      x_m[idim] = 0.;
      for (int idofk = 0; idofk < el_ndof[2]; idofk++)
        x_m[idim] += _xx_qnds[idim * NDOF_FEM + idofk] / NDOF_FEM;
    }

    for (int iside = 0; iside < el_sides; iside++) {
      if (el_neigh[iside] == -1) {
        for (int idof = 0; idof < elb_ndof[2]; idof++) {
          sur_toply[idof] = _mgmesh._GeomEl._surf_top[idof + NDOF_FEMB * iside];  // local nodes
          int idofb = sur_toply[idof];                                            // connectivity vector
          _bc_bd[idofb] = 0;
          for (int idim = 0; idim < _nTdim; idim++) {
            _xxb_qnds[idim * NDOF_FEMB + idof] = _xx_qnds[idim * NDOF_FEM + idofb];  // coordinates
          }
        }
        int sign_normal = 1;
        _fe[2]->normal_g(_xxb_qnds, x_m, normal, sign_normal);
        bc_set(sur_toply, el_ndof[2], elb_ndof[2], elb_ngauss, sign_normal);
      }
    }
    }
    // ----------------------------------------------------------------------------------
    //   3. Volume integration
    // ----------------------------------------------------------------------------------

//     if( iel == jel ) assembly_laplacian ( el_ndof2, el_ngauss, mode );
//     if( iel != jel ) assembly_frac_lap ( el_ndof2, el_ngauss, el_ndof2, el_ngauss, 0.5); // manca l'add matrix!!
    if( iel == jel ) adaptive_ref_frac(el_ndof2, el_ngauss, el_ndof2, el_ngauss, 0.5, 1, iel );
    if( iel == jel ) rhs_assembly( el_ndof2, el_ngauss );
    
//     	if (iel == jel ){
//         std::cout <<"iel"<<iel<<  "  KeM \n"  <<_CClocal_refined << " rhs \n" <<  _FeM <<"\n\n" <<std::endl;
//     }
        
    
    
//     A[Level]->add_matrix(_KeM, el_dof_indices);  // global matrix
    if(iel == jel){
        A[Level]->add_matrix( _KeM, el_dof_indices, el_dof_indices );
        A[Level]->add_matrix( _CClocal_refined, el_dof_indices, el_dof_indices );
        if (mode == 1) {
            b[Level]->add_vector( _FeM, el_dof_indices ); 
            b[Level]->add_vector( _Res_local_refined, el_dof_indices );
        }
    }
//     else{
//         A[Level]->add_matrix( _CClocalII, el_dof_indices, el_dof_indices );
//         A[Level]->add_matrix( _CClocalIJ, el_dof_indices, el_dof_indices_j );
//         A[Level]->add_matrix( _CClocalJI, el_dof_indices_j, el_dof_indices );
//         A[Level]->add_matrix( _CClocalJJ, el_dof_indices_j, el_dof_indices_j );
//         if (mode == 1) {
//             b[Level]->add_vector( _Res_nonlocalI, el_dof_indices );
//             b[Level]->add_vector( _Res_nonlocalJ, el_dof_indices_j );
//         }
//     }
  }  // end of iel element loop

  
} // end jel loop
 } // end kproc loop
  
  /// 5. clean
  el_dof_indices.clear();
  A[Level]->close();
  if (mode == 1) { b[Level]->close(); }
  //   A[Level]->print(); b[Level]->print();
#ifdef PRINT_INFO
  std::cout << " Matrix Assembled(T)  for  Level " << Level << " dofs " << A[Level]->n() << "\n";
#endif
  return;
}

// =========================================================================================
/// This function controls the assembly and the solution of the T_equation system:
void MGSolT::MGTimeStep(
    const double time,  ///< time
    const int /*iter*/  ///< Number of max inter
) {
  // =========================================================================================

//   if (_SolveT) {
    std::cout << std::endl
              << "\033[038;5;" << 196 << ";1m "
              << "--------------------------------------------------- \n\t" << _eqname.c_str()
              << " solution of problem " << _mgutils.get_name()
              << "\n ---------------------------------------------------\n  \033[0m";

    std::clock_t start_time = std::clock();
    GenMatRhs(time, _NoLevels - 1, 1);  // matrix and rhs
    for (int Level = 0; Level < _NoLevels - 1; Level++) {
      GenMatRhs(time, Level, 0);  // matrix
    }
    std::clock_t end_time = std::clock();
    MGSolve(1.e-6, 40);
    std::clock_t end_time2 = std::clock();

#if PRINT_TIME == 1
    std::cout << " Ass. time -----> =" << double(end_time - start_time) / CLOCKS_PER_SEC << "s "
              << " Ass. and sol. time: =" << double(end_time2 - start_time) / CLOCKS_PER_SEC << "s "
              << std::endl;
#endif

    x_old[0][_NoLevels - 1]->localize(*x_old[1][_NoLevels - 1]);
    x[_NoLevels - 1]->localize(*x_old[0][_NoLevels - 1]);
//   }
  return;
}  // =======================================================================================

#endif
// #endif // personal application

// kate: indent-mode cstyle; indent-width 4; replace-tabs on;
