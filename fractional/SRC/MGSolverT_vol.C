#include "Equations_conf.h"

// ============================================
#ifdef T_EQUATIONS  // 3D-2D Energy equation
// ============================================

#include "MGFE.h"       // Mesh class,double vel_g[]
#include "MGSolverT.h"  // Navier-Stokes class header file

// ===============================================================================================
void MGSolT::vol_integral(
    const int el_ndof2, const int el_ngauss,
    const int
        mode) {  // ==============================================================================================

  double vel_g[DIMENSION];
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

    // Velocity field -> [NS_F] -> (quad, _indx_eqs[NS_F]) ----------------------------------------->
    for (int idim = 0; idim < _nTdim; idim++) vel_g[idim] = 0.;
    if (_FF_idx[NS_F] > -1)
      for (int idim = 0; idim < _nTdim; idim++)
        vel_g[idim] = _ub_g[2][_FF_idx[NS_F] + idim];  // velocity field

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

          _KeM(i, j) += dtxJxW_g * (                                       // energy-equation matrix
                                       + Lap                               // diffusion term
                                   );
        }
      }
    }  // ----------------------------------------
  }    // end of the quadrature point qp-loop ***********************
  return;
}

#endif
