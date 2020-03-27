#include "Equations_conf.h"

// ============================================
#ifdef T_EQUATIONS  // 3D-2D Energy equation
// ============================================

#include "MGFE.h"       // Mesh class,double vel_g[]
#include "MGSolverT.h"  // Navier-Stokes class header file

// ===============================================================================================
void MGSolT::assembly_laplacian(
    const int el_ndof2, const int el_ngauss,
    const int
        mode) {  // ==============================================================================================

  double xyz_g[DIMENSION];
  double T_1ts_g[1], T_2ts_g[1];
  // --------------------------------------------------------------------------------------------------------------------
  /// c) gaussian integration loop (n_gauss)
  // ------------------------------------------------------------------------------------------------------------------

  for (int qp = 0; qp < el_ngauss; qp++) {
    // shape functions at gaussian points -----------------------------------
    double det2 = _fe[2]->Jac(qp, _xx_qnds, _InvJac2);          // Jacobian
    double JxW_g2 = det2 * _fe[2]->_weight1[_nTdim - 1][qp];    // weight
    _fe[2]->get_phi_gl_g(_nTdim, qp, _phi_g[2]);                // shape funct
    _fe[2]->get_dphi_gl_g(_nTdim, qp, _InvJac2, _dphi_g[2]);    // global coord deriv
    _fe[2]->get_ddphi_gl_g(_nTdim, qp, _InvJac2, _ddphi_g[2]);  // local second deriv
    
    //  fields
    //  --------------------------------------------------------------------------------------------------------
    interp_el_sol(_xx_qnds, 0, _nTdim, _phi_g[2], el_ndof2, xyz_g);
    interp_el_sol(_T_1ts, 0, 1, _phi_g[2], el_ndof2, T_1ts_g);
    interp_el_sol(_T_2ts, 0, 1, _phi_g[2], el_ndof2, T_2ts_g);
    interp_el_sol(
        _data_eq[2].ub, 0, _data_eq[2].indx_ub[_data_eq[2].n_eqs], _phi_g[2], el_ndof2,
        _ub_g[2]);  // quadratic

    if (_Axisym == 1) JxW_g2 *= xyz_g[0];

    /// d) Local (element) assemblying energy equation
    // =====================================================================================================================
    for (int i = 0; i < el_ndof2; i++) {
      const double phii_g = _phi_g[2][i];
      double dtxJxW_g = JxW_g2;  // area with bc and weight
      if (_bc_el[i] == 1) {

        // Rhs Assemblying
        // ---------------------------------------------------------------------------------------------------
        if (mode == 1) {  // rhs
          double TimeDerivative = 0.;
          double SourceTerms = dtxJxW_g * (phii_g ) * ( _qheat );  // source term
          double AddedTerm = 1. * dtxJxW_g * phii_g;
          _FeM(i) += SourceTerms + AddedTerm;
        }

        // Matrix Assemblying
        // ------------------------------------------------------------------------------------------------
        for (int j = 0; j < el_ndof2; j++) {
          const double phij_g = _phi_g[2][j];
          double Lap = 0.;
          ;
          for (int idim = 0; idim < _nTdim; idim++) {
            const double dphiidxg = _dphi_g[2][i + idim * el_ndof2];
            const double dphijdxg = _dphi_g[2][j + idim * el_ndof2];
            Lap += dphijdxg * dphiidxg;  // diffusion
          }

          _KeM(i, j) += dtxJxW_g * Lap;
        }
      }
    }  // ----------------------------------------
  }    // end of the quadrature point qp-loop ***********************
  return;
}

// ===============================================================================================
void MGSolT::rhs_assembly( const int el_ndof2, const int el_ngauss ){  
// ==============================================================================================

  double xyz_g[DIMENSION];
  // --------------------------------------------------------------------------------------------------------------------
  /// c) gaussian integration loop (n_gauss)
  // ------------------------------------------------------------------------------------------------------------------

  for (int qp = 0; qp < el_ngauss; qp++) {
    // shape functions at gaussian points -----------------------------------
    double det2 = _fe[2]->Jac(qp, _xx_qnds, _InvJac2);          // Jacobian
    double JxW_g2 = det2 * _fe[2]->_weight1[_nTdim - 1][qp];    // weight
    _fe[2]->get_phi_gl_g(_nTdim, qp, _phi_g[2]);                // shape funct
    
    //  fields
    //  --------------------------------------------------------------------------------------------------------
    interp_el_sol(_xx_qnds, 0, _nTdim, _phi_g[2], el_ndof2, xyz_g);

    if (_Axisym == 1) JxW_g2 *= xyz_g[0];

    /// d) Local (element) assemblying energy equation
    // =====================================================================================================================
    for (int i = 0; i < el_ndof2; i++) {
      const double phii_g = _phi_g[2][i];
      double dtxJxW_g = JxW_g2;  // area with bc and weight
      if (_bc_el[i] == 1) {

        // Rhs Assemblying
        // ---------------------------------------------------------------------------------------------------
          double AddedTerm = 1. * dtxJxW_g * phii_g;
          _FeM(i) += AddedTerm;

      }
    }  // ----------------------------------------
  }    // end of the quadrature point qp-loop ***********************
  return;
}



void MGSolT::assembly_frac_lap( const int el_ndof_i, const int el_ngauss_i, 
                            const int el_ndof_j, const int el_ngauss_j, double s_frac ) {
  double phi1[el_ndof_i];
  double phi2[el_ndof_j];
  double xyz_g_i[DIMENSION];
  double xyz_g_j[DIMENSION];
  double solX[1];
  double solY[1];
  
  double C_ns = s_frac * pow(2, (2. * s_frac)) * tgamma((_nTdim + 2. * s_frac) / 2.) /
                (pow(M_PI, _nTdim / 2.) * tgamma(1 -  s_frac)) ;
  
  for (unsigned ig = 0; ig < el_ngauss_i; ig++) {
    double det_i = _fe[2]->Jac(ig, _xx_qnds, _InvJac2);          // Jacobian
    double JxW_g_i = det_i * _fe[2]->_weight1[_nTdim - 1][ig];    // weight
    
    _fe[2]->get_phi_gl_g(_nTdim, ig, phi1);                // shape funct

    //  fields
    //  --------------------------------------------------------------------------------------------------------
    interp_el_sol(_xx_qnds, 0, _nTdim, _phi_g[2], el_ndof_i, xyz_g_i);
    interp_el_sol(_T_1ts, 0, 1, _phi_g[2], el_ndof_i, solX);
    interp_el_sol(
        _data_eq[2].ub, 0, _data_eq[2].indx_ub[_data_eq[2].n_eqs], _phi_g[2], el_ndof_i,_ub_g[2]);  // quadratic
              
    for(unsigned jg = 0; jg < el_ngauss_j; jg++) {
      double det_j = _fe[2]->Jac(jg, _xx_qnds_j, _InvJacJ);          // Jacobian
      double JxW_g_j = det_j * _fe[2]->_weight1[_nTdim - 1][jg];    // weight
    
      _fe[2]->get_phi_gl_g(_nTdim, jg, phi2);                // shape funct
      
      interp_el_sol(_xx_qnds_j, 0, _nTdim, phi2, el_ndof_j, xyz_g_j);
      interp_el_sol(_T_jel, 0, 1, phi2, el_ndof_j, solY);
      
      double dist_xyz = 0;
      for(unsigned k = 0; k < _nTdim; k++) {
        dist_xyz += (xyz_g_i[k] - xyz_g_j[k]) * (xyz_g_i[k] - xyz_g_j[k]);
      }

      const double denom = pow(dist_xyz, (double)((_nTdim / 2.) + s_frac));
       

      for(unsigned i = 0; i < el_ndof_i; i++) {
        
        if (_bc_el[i] == 1) {

          _Res_nonlocalI( i )         +=  - (C_ns / 2.) * (solX[0] - solY[0]) * (phi1[i]) * JxW_g_i *JxW_g_j  / denom;

          _Res_nonlocalJ( i )         +=  - (C_ns / 2.)  * (solX[0] - solY[0]) * (- phi2[i]) * JxW_g_i *JxW_g_j  / denom;

          for(unsigned j = 0; j < el_ndof_j; j++) {

            _CClocalII(i, j) += (C_ns / 2.)  * phi1[j]  * phi1[i] * JxW_g_i *JxW_g_j / denom;
                              
            _CClocalIJ(i, j) += (C_ns / 2.)  * (- phi2[j]) * phi1[i] * JxW_g_i *JxW_g_j / denom;
                              
            _CClocalJI(i, j) += (C_ns / 2.)  * phi1[j] * (- phi2[i]) * JxW_g_i *JxW_g_j / denom;
                              
            _CClocalJJ(i, j) += (C_ns / 2.)  * (- phi2[j]) * (- phi2[i]) * JxW_g_i *JxW_g_j / denom;


          } // j-loop
        }
      }   // i-loop
    }     // jg-loop
  }       // ig-loop
  return;
}

// // Assembly for adaptive refinement when iel == jel
void MGSolT::adaptive_ref_frac( const int el_ndof_i, const int el_ngauss_i, 
                            const int el_ndof_j, const int el_ngauss_j, double s_frac, unsigned const Nsplit, int iel ){
  
  double phi3[el_ndof_j];
  double phi1[el_ndof_i];
  double xyz_g_i[DIMENSION];
  double solX[1];
  double solY3[1];
  
  double C_ns = s_frac * pow(2, (2. * s_frac)) * tgamma((_nTdim + 2. * s_frac) / 2.) /
                (pow(M_PI, _nTdim / 2.) * tgamma(1 -  s_frac)) ;
  
  std::cout.precision(14);
  std::vector< std::vector<std::vector<double>>> x3;
  std::vector <double > xg1;
  std::vector < std::vector <double > > xNodes;
  xg1.resize(_nTdim);
  xNodes.resize(_nTdim);
  for(unsigned k = 0; k < _nTdim; k++ ) xNodes[k].resize(el_ndof_j);
    
  for(unsigned ig = 0; ig < el_ngauss_i; ig++) {
    double det_i = _fe[2]->Jac(ig, _xx_qnds, _InvJac2);          // Jacobian
    double JxW_g_i = det_i * _fe[2]->_weight1[_nTdim - 1][ig];    // weight
    
    _fe[2]->get_phi_gl_g(_nTdim, ig, phi1);                // shape funct

    //  fields
    //  --------------------------------------------------------------------------------------------------------
    interp_el_sol(_xx_qnds, 0, _nTdim, _phi_g[2], el_ndof_i, xyz_g_i);
    interp_el_sol(_T_1ts, 0, 1, _phi_g[2], el_ndof_i, solX);
    interp_el_sol(
        _data_eq[2].ub, 0, _data_eq[2].indx_ub[_data_eq[2].n_eqs], _phi_g[2], el_ndof_i,_ub_g[2]);  // quadratic
    
    for(unsigned split = 0; split <= Nsplit; split++) {

      for(unsigned k = 0; k < _nTdim; k++){
        xg1[k] = xyz_g_i[k];
        for(unsigned dd = 0; dd < el_ndof_j; dd++){
          xNodes[k][dd] = _xx_qnds[ dd + k * el_ndof_j ];
        }
      }

      if(_nTdim == 1) GetElementPartition1D(xg1, xNodes, split, Nsplit, x3);
      else if(_nTdim == 2) GetElementPartitionQuad(xg1, xNodes, split, Nsplit, x3);
      else {
        std::cerr<< "Adaptive Gauss points: 3D NOT SUPPORTED!\n";
        break;
      }
      
      for(unsigned r = 0; r < x3.size(); r++) {
        
        for(unsigned jg = 0; jg < el_ngauss_j; jg++) {
          
          double x_loc[el_ndof_i * _nTdim]; /**< elem coords */
          for(unsigned d = 0; d < el_ndof_i; d++){
            for(unsigned k = 0; k < _nTdim; k++){
              x_loc[d + k * el_ndof_i] = x3[r][k][d];
            }
          }
          
          double det_ref = _fe[2]->Jac(jg, x_loc, _InvJacJ);          // Jacobian
          double JxW_g_ref = det_ref * _fe[2]->_weight1[_nTdim - 1][jg];    // weight
          _fe[2]->get_phi_gl_g(_nTdim, jg, phi3);                // shape funct
          
          vector < double > xg3(_nTdim, 0.);
          
          for(unsigned i = 0; i < el_ndof_i; i++) {
            for(unsigned k = 0; k < _nTdim; k++) {
              xg3[k] += x3[r][k][i] * phi3[i];
            }
          }
          
          std::vector <double> lin_nodes;
          std::vector <double> can_pos;
          lin_nodes.resize(4);
          for(int i = 0; i < _fe[1]->_NoShape[_nTdim - 1]; i++) lin_nodes[i] = _xx_qnds[i];
          for(int i = el_ndof_i; i < el_ndof_i + _fe[1]->_NoShape[_nTdim - 1]; i++) lin_nodes.push_back( _xx_qnds[i] );
          can_pos = _fe[2]->XiEtaChiCalc(_nTdim, xg3, lin_nodes);
          for(int i = 0; i < el_ndof_i; i++){
            phi3[i] = _fe[2]->QuadPhi(i, can_pos, 2);
          }
          
          interp_el_sol(_T_1ts, 0, 1, phi3, el_ndof_i, solY3);
          
          double dist_xyz3 = 0;
          for(unsigned k = 0; k < _nTdim; k++) {
            dist_xyz3 += (xyz_g_i[k] - xg3[k]) * (xyz_g_i[k] - xg3[k]);
          }
                     
          const double denom3 = pow(dist_xyz3, (double)((_nTdim / 2.) + s_frac));

          for(unsigned i = 0; i < el_ndof_i; i++) {
            
            if (_bc_el[i] == 1) {

              _Res_local_refined( i )  +=  - (C_ns / 2.) * ((solX[0] - solY3[0]) * (phi1[i] - phi3[i]) * JxW_g_ref / denom3 ) * JxW_g_i ;

              for(unsigned j = 0; j < el_ndof_j; j++) {
                _CClocal_refined( i , j ) += (C_ns / 2.) * ( (phi1[j] /*- phi3[j]*/) * (phi1[i] - phi3[i]) * 
                                                    JxW_g_ref / denom3 ) * JxW_g_i ;
                                                    
                                                    
              }
            }
          }
        }
      }
    }
    
  }
  return;
}

// // Function to get the coordiates of the points in adaptive split for fractional calculus in 1D
void MGSolT::GetElementPartition1D(const std::vector <double >  & xg1, const std::vector < std::vector <double > > & x1, 
                           const unsigned &split , const unsigned & Nsplit ,
                           std::vector < std::vector < std::vector<double>>> &x)
{
  unsigned dim = 1;
  unsigned left = 0;
  unsigned right = 1;

  if(split == 0) { //init
    x.resize(2);
    x[left].resize(dim);
    x[right].resize(dim);
    for(unsigned k = 0; k < dim; k++) {
      x[left][k].resize(x1[0].size());
      x[right][k].resize(x1[0].size());
//       for(unsigned k = 0; k < dim; k++) {
      x[left][k][0] = x1[k][0];
      x[left][k][1] = 0.5 * (x[left][k][0] + xg1[k]);
      x[left][k][2] = 0.5 * (x[left][k][0] + x[left][k][1]);
      x[right][k][1] = x1[k][1];
      x[right][k][0] = 0.5 * (x[right][k][1] + xg1[k]);
      x[right][k][2] = 0.5 * (x[right][k][0] + x[right][k][1]);
//       }
    }
  }
  else if(split == Nsplit) {
    for(unsigned k = 0; k < dim; k++) {
      x[left][k][0] = x[left][k][1];
      x[left][k][1] = xg1[k];
      x[left][k][2] = 0.5 * (x[left][k][0] + x[left][k][1]);

      x[right][k][1] = x[right][k][0];
      x[right][k][0] = xg1[k];
      x[right][k][2] = 0.5 * (x[right][k][0] + x[right][k][1]);
    }
  }
  else {
    for(unsigned k = 0; k < dim; k++) {
      x[left][k][0] = x[left][k][1];
      x[left][k][1] = 0.5 * (x[left][k][0] + xg1[k]);
      x[left][k][2] = 0.5 * (x[left][k][0] + x[left][k][1]);

      x[right][k][1] = x[right][k][0];
      x[right][k][0] = 0.5 * (x[right][k][1] + xg1[k]);
      x[right][k][2] = 0.5 * (x[right][k][0] + x[right][k][1]);
    }
  }
}
const unsigned ijndex[2][12][2] = {
  { {0, 0}, {3, 0}, {3, 3}, {0, 3},
    {1, 0}, {0, 1},
    {2, 0}, {3, 1},
    {2, 3}, {3, 2},
    {1, 3}, {0, 2}
  },
  {{0, 0}, {1, 0}, {1, 1}, {0, 1}}
};
// // Function to get the coordiates of the points in adaptive split for fractional calculus
void MGSolT::GetElementPartitionQuad(const std::vector <double >  & xg1, const std::vector < std::vector <double > > & xNodes, const unsigned & split, const unsigned & totalNumberofSplits,  std::vector < std::vector < std::vector<double>>> &x)
{
  unsigned dim = 2;

  unsigned solType;
  unsigned size = xNodes[0].size();

  if(size == 4) {
    solType = 0; //lagrange linear
  }
  else if(size == 8) {
    solType = 1; //lagrange serendipity
  }
  else if(size == 9) {
    solType = 2; //lagrange quadratic
  }
  else {
    std::cout << "abort in GetElementPartitionQuad" << std::endl;
    abort();
  }


  unsigned bl = 0; // bottom left
  unsigned br = 1; // bottom right
  unsigned tr = 2; // top left
  unsigned tl = 3; // top right

  std::vector < double > XX;
  std::vector < double > YY;

  unsigned size_part = 12;
  unsigned splitType = 0;

  if(split < totalNumberofSplits) { //init && update

    XX.resize(5);
    YY.resize(5);

    if(split == 0) { //init

      x.resize(size_part);
      for(unsigned j = 0; j < size_part; j++) {
        x[j].resize(dim);
        for(unsigned k = 0; k < dim; k++) {
          x[j][k].resize(size);
        }
      }

      XX[0] = xNodes[0][0];
      XX[4] = xNodes[0][1];
      YY[0] = xNodes[1][0];
      YY[4] = xNodes[1][3];

    }
    else { //update
      XX[0] = x[bl][0][1];
      XX[4] = x[br][0][0];
      YY[0] = x[bl][1][2];
      YY[4] = x[tl][1][1];
    }
    XX[2] = xg1[0];
    XX[1] = 0.5 * (XX[0] + XX[2]);
    XX[3] = 0.5 * (XX[2] + XX[4]);

    YY[2] = xg1[1];
    YY[1] = 0.5 * (YY[0] + YY[2]);
    YY[3] = 0.5 * (YY[2] + YY[4]);
  }
  else { //close

    XX.resize(3);
    YY.resize(3);

    XX[0] = x[bl][0][1];
    XX[1] = xg1[0];
    XX[2] = x[br][0][0];
    YY[0] = x[bl][1][2];
    YY[1] = xg1[1];
    YY[2] = x[tl][1][1];

    size_part = 4;
    splitType = 1;
    x.resize(size_part);
    for(unsigned j = 0; j < size_part; j++) {
      x[j].resize(dim);
      for(unsigned k = 0; k < dim; k++) {
        x[j][k].resize(size);
      }
    }
  }

  for(unsigned qq = 0; qq < size_part; qq++) {
    unsigned i = ijndex[splitType][qq][0];
    x[qq][0][0] = x[qq][0][3] = XX[i];
    x[qq][0][1] = x[qq][0][2] = XX[i + 1];

    unsigned j = ijndex[splitType][qq][1];
    x[qq][1][0] = x[qq][1][1] = YY[j];
    x[qq][1][2] = x[qq][1][3] = YY[j + 1];
  }
  if(solType > 0) {
    for(unsigned qq = 0; qq < size_part; qq++) {
      for(unsigned k = 0; k < dim; k++) { //middle point formula
        x[qq][k][4] = 0.5 * (x[qq][k][0] + x[qq][k][1]);
        x[qq][k][5] = 0.5 * (x[qq][k][1] + x[qq][k][2]);
        x[qq][k][6] = 0.5 * (x[qq][k][2] + x[qq][k][3]);
        x[qq][k][7] = 0.5 * (x[qq][k][3] + x[qq][k][0]);
      }
    }
  }

  if(solType > 1) {
    for(unsigned qq = 0; qq < size_part; qq++) {
      for(unsigned k = 0; k < dim; k++) { //middle point formula
        x[qq][k][8] = 0.5 * (x[qq][k][0] + x[qq][k][2]);
      }
    }
  }

}

#endif
