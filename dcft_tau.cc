/*
 *@BEGIN LICENSE
 *
 * PSI4: an ab initio quantum chemistry software package
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 *@END LICENSE
 */

#include "dcft.h"
#include <libdpd/dpd.h>
#include <libqt/qt.h>
#include <libiwl/iwl.hpp>
#include <libmints/molecule.h>
#include <psifiles.h>
#include <libtrans/integraltransform.h>
#include "defines.h"

using namespace std;

namespace psi{ namespace dfdcft{

/**
 * Forms Tau in the MO basis from the Lambda tensors and transforms it back
 * to the SO basis.
 */
void
DCFTSolver::build_tau()
{
    dcft_timer_on("DCFTSolver::build_tau()");
    dpdbuf4 L1, L2;
    dpdfile2 T_OO, T_oo, T_VV, T_vv;

    // If the system is closed-shell, then...
    if(closedshell_){
        /*
         * The following 4 lines are not able to dump Tau <O|O>'s onto PSIF_DCFT_DPD
         * T_OO's are assigned zero at the first time but it has correct symmetry blocks (correct size)
         * Tau <O|O> are save on PSIF_DCFT_DPD after contract442
         */
        global_dpd_->file2_init(&T_OO, PSIF_DCFT_DPD, 0, ID('O'), ID('O'), "Tau <O|O>");
        global_dpd_->file2_init(&T_VV, PSIF_DCFT_DPD, 0, ID('V'), ID('V'), "Tau <V|V>");

        global_dpd_->buf4_init(&L1, PSIF_DCFT_DPD, 0,
                      _ints->DPD_ID("[O,O]"), _ints->DPD_ID("[V,V]"),
                      ID("[O,O]"), ID("[V,V]"), 0, "Lambda <OO|VV>");
        global_dpd_->buf4_init(&L2, PSIF_DCFT_DPD, 0,
                      _ints->DPD_ID("[O,O]"), _ints->DPD_ID("[V,V]"),
                      ID("[O,O]"), ID("[V,V]"), 0, "Lambda <OO|VV>");

        /*
         * Tau_IJ = -1/2 Lambda_IKAB Lambda_JKAB
         */
        global_dpd_->contract442(&L1, &L2, &T_OO, 0, 0, -0.5, 0.0);
        /*
         * Tau_AB = +1/2 Lambda_IJAC Lambda_IJBC
         */
        global_dpd_->contract442(&L1, &L2, &T_VV, 2, 2, 0.5, 0.0);
        global_dpd_->buf4_close(&L1);
        global_dpd_->buf4_close(&L2);

//        global_dpd_->buf4_init(&L1, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[V,v]"),
//                               ID("[O,o]"), ID("[V,v]"), 0, "Lambda <Oo|Vv>");
//        global_dpd_->buf4_init(&L2, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[V,v]"),
//                               ID("[O,o]"), ID("[V,v]"), 0, "Lambda <Oo|Vv>");

        global_dpd_->buf4_init(&L1, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                               ID("[O,O]"), ID("[V,V]"), 0, "Lambda SF <OO|VV>");
        global_dpd_->buf4_init(&L2, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                               ID("[O,O]"), ID("[V,V]"), 0, "Lambda SF <OO|VV>");

        /*
         * Tau_IJ -= 1/2 Lambda_IkAb Lambda_JkAb - 1/2 Lambda_IkaB Lambda_JkaB
         */
        global_dpd_->contract442(&L1, &L2, &T_OO, 0, 0, -1.0, 1.0);
        /*
         * Tau_AB += 1/2 Lambda_IjAc Lambda_IjBc + 1/2 Lambda_iJAc Lambda_iJBc
         */
        global_dpd_->contract442(&L1, &L2, &T_VV, 2, 2, 1.0, 1.0);

//        /*
//         * Tau_ij = Tau_IJ
//         */
//        global_dpd_->file2_copy(&T_OO, PSIF_DCFT_DPD, "Tau <o|o>");
//        global_dpd_->file2_init(&T_oo, PSIF_DCFT_DPD, 0, ID('o'), ID('o'), "Tau <o|o>");
//        /*
//         * Tau_ab = T_AB
//         */
//        global_dpd_->file2_copy(&T_VV, PSIF_DCFT_DPD, "Tau <v|v>");
//        global_dpd_->file2_init(&T_vv, PSIF_DCFT_DPD, 0, ID('v'), ID('v'), "Tau <v|v>");

//        global_dpd_->buf4_close(&L1);
//        global_dpd_->buf4_close(&L2);
        // Check:
//        global_dpd_->file2_print(&T_OO, "outfile");

        global_dpd_->file2_close(&T_OO);
//        global_dpd_->file2_close(&T_oo);
        global_dpd_->file2_close(&T_VV);
//        global_dpd_->file2_close(&T_vv);

        // Read MO-basis Tau from disk into the memory
        global_dpd_->file2_init(&T_OO, PSIF_DCFT_DPD, 0, ID('O'), ID('O'), "Tau <O|O>");
//        global_dpd_->file2_init(&T_oo, PSIF_DCFT_DPD, 0, ID('o'), ID('o'), "Tau <o|o>");
        global_dpd_->file2_init(&T_VV, PSIF_DCFT_DPD, 0, ID('V'), ID('V'), "Tau <V|V>");
//        global_dpd_->file2_init(&T_vv, PSIF_DCFT_DPD, 0, ID('v'), ID('v'), "Tau <v|v>");

        global_dpd_->file2_mat_init(&T_OO);
//        global_dpd_->file2_mat_init(&T_oo);
        global_dpd_->file2_mat_init(&T_VV);
//        global_dpd_->file2_mat_init(&T_vv);

        global_dpd_->file2_mat_rd(&T_OO);
//        global_dpd_->file2_mat_rd(&T_oo);
        global_dpd_->file2_mat_rd(&T_VV);
//        global_dpd_->file2_mat_rd(&T_vv);

        for(int h = 0; h < nirrep_; ++h){
            for(int i = 0; i < naoccpi_[h]; ++i){
                for(int j = 0; j < naoccpi_[h]; ++j){
                    aocc_tau_->set(h, i, j, T_OO.matrix[h][i][j]);
                }
            }
            for(int a = 0; a < navirpi_[h]; ++a){
                for(int b = 0; b < navirpi_[h]; ++b){
                    avir_tau_->set(h, a, b, T_VV.matrix[h][a][b]);
                }
            }
//            for(int i = 0; i < nboccpi_[h]; ++i){
//                for(int j = 0; j < nboccpi_[h]; ++j){
//                    bocc_tau_->set(h, i, j, T_oo.matrix[h][i][j]);
//                }
//            }
//            for(int a = 0; a < nbvirpi_[h]; ++a){
//                for(int b = 0; b < nbvirpi_[h]; ++b){
//                    bvir_tau_->set(h, a, b, T_vv.matrix[h][a][b]);
//                }
//            }
        }

        bocc_tau_->copy(aocc_tau_);
        bvir_tau_->copy(avir_tau_);

        global_dpd_->file2_close(&T_OO);
//        global_dpd_->file2_close(&T_oo);
        global_dpd_->file2_close(&T_VV);
//        global_dpd_->file2_close(&T_vv);

        dcft_timer_off("DCFTSolver::build_tau()");
    }

    // If the system is open-shell, then...
    else{
        /*
         * The following 4 lines are not able to dump Tau <O|O>'s onto PSIF_DCFT_DPD
         * T_OO's are assigned zero at the first time but it has correct symmetry blocks (correct size)
         * Tau <O|O> are save on PSIF_DCFT_DPD after contract442
         */
        global_dpd_->file2_init(&T_OO, PSIF_DCFT_DPD, 0, ID('O'), ID('O'), "Tau <O|O>");
        global_dpd_->file2_init(&T_oo, PSIF_DCFT_DPD, 0, ID('o'), ID('o'), "Tau <o|o>");
        global_dpd_->file2_init(&T_VV, PSIF_DCFT_DPD, 0, ID('V'), ID('V'), "Tau <V|V>");
        global_dpd_->file2_init(&T_vv, PSIF_DCFT_DPD, 0, ID('v'), ID('v'), "Tau <v|v>");

        global_dpd_->buf4_init(&L1, PSIF_DCFT_DPD, 0,
                      _ints->DPD_ID("[O,O]"), _ints->DPD_ID("[V,V]"),
                      ID("[O>O]-"), ID("[V>V]-"), 0, "Lambda <OO|VV>");
        global_dpd_->buf4_init(&L2, PSIF_DCFT_DPD, 0,
                      _ints->DPD_ID("[O,O]"), _ints->DPD_ID("[V,V]"),
                      ID("[O>O]-"), ID("[V>V]-"), 0, "Lambda <OO|VV>");
    //    global_dpd_->buf4_init(&L1, PSIF_DCFT_DPD, 0,
    //                  ID("[O,O]"), ID("[V,V]"),
    //                  ID("[O>O]-"), ID("[V>V]-"), 0, "Lambda <OO|VV>");
    //    global_dpd_->buf4_init(&L2, PSIF_DCFT_DPD, 0,
    //                  ID("[O,O]"), ID("[V,V]"),
    //                  ID("[O>O]-"), ID("[V>V]-"), 0, "Lambda <OO|VV>");

        /*
         * Tau_IJ = -1/2 Lambda_IKAB Lambda_JKAB
         */
        global_dpd_->contract442(&L1, &L2, &T_OO, 0, 0, -0.5, 0.0);
        /*
         * Tau_AB = +1/2 Lambda_IJAC Lambda_IJBC
         */
        global_dpd_->contract442(&L1, &L2, &T_VV, 2, 2, 0.5, 0.0);
        global_dpd_->buf4_close(&L1);
        global_dpd_->buf4_close(&L2);

        global_dpd_->buf4_init(&L1, PSIF_DCFT_DPD, 0,
                      _ints->DPD_ID("[o,o]"), _ints->DPD_ID("[v,v]"),
                      ID("[o>o]-"), ID("[v>v]-"), 0, "Lambda <oo|vv>");
        global_dpd_->buf4_init(&L2, PSIF_DCFT_DPD, 0,
                      _ints->DPD_ID("[o,o]"), _ints->DPD_ID("[v,v]"),
                      ID("[o>o]-"), ID("[v>v]-"), 0, "Lambda <oo|vv>");
        /*
         * Tau_ij = -1/2 Lambda_ikab Lambda_jkab
         */
        global_dpd_->contract442(&L1, &L2, &T_oo, 0, 0, -0.5, 0.0);
        /*
         * Tau_ab = +1/2 Lambda_ijac Lambda_ijbc
         */
        global_dpd_->contract442(&L1, &L2, &T_vv, 2, 2, 0.5, 0.0);
        global_dpd_->buf4_close(&L1);
        global_dpd_->buf4_close(&L2);

        global_dpd_->buf4_init(&L1, PSIF_DCFT_DPD, 0,
                      _ints->DPD_ID("[O,o]"), _ints->DPD_ID("[V,v]"),
                      _ints->DPD_ID("[O,o]"), _ints->DPD_ID("[V,v]"),
                      0, "Lambda <Oo|Vv>");
        global_dpd_->buf4_init(&L2, PSIF_DCFT_DPD, 0,
                      _ints->DPD_ID("[O,o]"), _ints->DPD_ID("[V,v]"),
                      _ints->DPD_ID("[O,o]"), _ints->DPD_ID("[V,v]"),
                      0, "Lambda <Oo|Vv>");
        /*
         * Tau_IJ -= 1/2 Lambda_IkAb Lambda_JkAb - 1/2 Lambda_IkaB Lambda_JkaB
         */
        global_dpd_->contract442(&L1, &L2, &T_OO, 0, 0, -1.0, 1.0);
        /*
         * Tau_ij -= 1/2 Lambda_KiAb Lambda_KjAb - 1/2 Lambda_KiaB Lambda_KjaB
         */
        global_dpd_->contract442(&L1, &L2, &T_oo, 1, 1, -1.0, 1.0);
        /*
         * Tau_AB += 1/2 Lambda_IjAc Lambda_IjBc + 1/2 Lambda_iJAc Lambda_iJBc
         */
        global_dpd_->contract442(&L1, &L2, &T_VV, 2, 2, 1.0, 1.0);
        /*
         * Tau_ab += 1/2 Lambda_IjCa Lambda_IjCb + 1/2 Lambda_iJCa Lambda_iJCb
         */
        global_dpd_->contract442(&L1, &L2, &T_vv, 3, 3, 1.0, 1.0);
        global_dpd_->buf4_close(&L1);
        global_dpd_->buf4_close(&L2);
        global_dpd_->file2_close(&T_OO);
        global_dpd_->file2_close(&T_oo);
        global_dpd_->file2_close(&T_VV);
        global_dpd_->file2_close(&T_vv);

        // Read MO-basis Tau from disk into the memory
        global_dpd_->file2_init(&T_OO, PSIF_DCFT_DPD, 0, ID('O'), ID('O'), "Tau <O|O>");
        global_dpd_->file2_init(&T_oo, PSIF_DCFT_DPD, 0, ID('o'), ID('o'), "Tau <o|o>");
        global_dpd_->file2_init(&T_VV, PSIF_DCFT_DPD, 0, ID('V'), ID('V'), "Tau <V|V>");
        global_dpd_->file2_init(&T_vv, PSIF_DCFT_DPD, 0, ID('v'), ID('v'), "Tau <v|v>");

        global_dpd_->file2_mat_init(&T_OO);
        global_dpd_->file2_mat_init(&T_oo);
        global_dpd_->file2_mat_init(&T_VV);
        global_dpd_->file2_mat_init(&T_vv);

        global_dpd_->file2_mat_rd(&T_OO);
        global_dpd_->file2_mat_rd(&T_oo);
        global_dpd_->file2_mat_rd(&T_VV);
        global_dpd_->file2_mat_rd(&T_vv);

        for(int h = 0; h < nirrep_; ++h){
            for(int i = 0; i < naoccpi_[h]; ++i){
                for(int j = 0; j < naoccpi_[h]; ++j){
                    aocc_tau_->set(h, i, j, T_OO.matrix[h][i][j]);
                }
            }
            for(int a = 0; a < navirpi_[h]; ++a){
                for(int b = 0; b < navirpi_[h]; ++b){
                    avir_tau_->set(h, a, b, T_VV.matrix[h][a][b]);
                }
            }
            for(int i = 0; i < nboccpi_[h]; ++i){
                for(int j = 0; j < nboccpi_[h]; ++j){
                    bocc_tau_->set(h, i, j, T_oo.matrix[h][i][j]);
                }
            }
            for(int a = 0; a < nbvirpi_[h]; ++a){
                for(int b = 0; b < nbvirpi_[h]; ++b){
                    bvir_tau_->set(h, a, b, T_vv.matrix[h][a][b]);
                }
            }
        }

        global_dpd_->file2_close(&T_OO);
        global_dpd_->file2_close(&T_oo);
        global_dpd_->file2_close(&T_VV);
        global_dpd_->file2_close(&T_vv);

        dcft_timer_off("DCFTSolver::build_tau()");
    }

}

void
DCFTSolver::transform_tau()
{
    // If the system is closed-shell, then ...
    if (closedshell_){
        dcft_timer_on("DCFTSolver::transform_tau()");

        dpdfile2 T_OO, T_VV;

        global_dpd_->file2_init(&T_OO, PSIF_DCFT_DPD, 0, ID('O'), ID('O'), "Tau <O|O>");
        global_dpd_->file2_init(&T_VV, PSIF_DCFT_DPD, 0, ID('V'), ID('V'), "Tau <V|V>");

        global_dpd_->file2_mat_init(&T_OO);
        global_dpd_->file2_mat_init(&T_VV);
        global_dpd_->file2_mat_rd(&T_OO);
        global_dpd_->file2_mat_rd(&T_VV);

        // Zero SO tau arrays before computing it in the MO basis
        tau_so_a_->zero();

        for(int h = 0; h < nirrep_; ++h){
            if(nsopi_[h] == 0) continue;

            double **temp = block_matrix(nsopi_[h], nsopi_[h]);
            /*
             * Backtransform the Tau matrices to the AO basis: soTau = C moTau Ct
             * Attn: the forward MO->AO transformation would be: moTau = Ct S soTau S C
             */
            double **paOccC = aocc_c_->pointer(h);
            double **paVirC = avir_c_->pointer(h);
            double **pa_tau_ = tau_so_a_->pointer(h);

            // Alpha occupied
            if(naoccpi_[h] && nsopi_[h]){
                C_DGEMM('n', 'n', nsopi_[h], naoccpi_[h], naoccpi_[h], 1.0, paOccC[0], naoccpi_[h],
                        T_OO.matrix[h][0], naoccpi_[h], 0.0, temp[0], nsopi_[h]);
                C_DGEMM('n', 't', nsopi_[h], nsopi_[h], naoccpi_[h], 1.0, temp[0], nsopi_[h],
                        paOccC[0], naoccpi_[h], 1.0, pa_tau_[0], nsopi_[h]);
            }

            // Alpha virtual
            if(navirpi_[h] && nsopi_[h]){
                C_DGEMM('n', 'n', nsopi_[h], navirpi_[h], navirpi_[h], 1.0, paVirC[0], navirpi_[h],
                        T_VV.matrix[h][0], navirpi_[h], 0.0, temp[0], nsopi_[h]);
                C_DGEMM('n', 't', nsopi_[h], nsopi_[h], navirpi_[h], 1.0, temp[0], nsopi_[h],
                        paVirC[0], navirpi_[h], 1.0, pa_tau_[0], nsopi_[h]);
            }

            free_block(temp);
        }

        global_dpd_->file2_close(&T_OO);
        global_dpd_->file2_close(&T_VV);

        // Copy tau_so_alpha to tau_so_beta
        tau_so_b_->copy(tau_so_a_);

        dcft_timer_off("DCFTSolver::transform_tau()");
    }
    // If the system is open-shell, then ...
    else{
        dcft_timer_on("DCFTSolver::transform_tau()");

        dpdfile2 T_OO, T_oo, T_VV, T_vv;

        global_dpd_->file2_init(&T_OO, PSIF_DCFT_DPD, 0, ID('O'), ID('O'), "Tau <O|O>");
        global_dpd_->file2_init(&T_oo, PSIF_DCFT_DPD, 0, ID('o'), ID('o'), "Tau <o|o>");
        global_dpd_->file2_init(&T_VV, PSIF_DCFT_DPD, 0, ID('V'), ID('V'), "Tau <V|V>");
        global_dpd_->file2_init(&T_vv, PSIF_DCFT_DPD, 0, ID('v'), ID('v'), "Tau <v|v>");

        global_dpd_->file2_mat_init(&T_OO);
        global_dpd_->file2_mat_init(&T_oo);
        global_dpd_->file2_mat_init(&T_VV);
        global_dpd_->file2_mat_init(&T_vv);
        global_dpd_->file2_mat_rd(&T_OO);
        global_dpd_->file2_mat_rd(&T_oo);
        global_dpd_->file2_mat_rd(&T_VV);
        global_dpd_->file2_mat_rd(&T_vv);

        // Zero SO tau arrays before computing it in the MO basis
        tau_so_a_->zero();
        tau_so_b_->zero();

        for(int h = 0; h < nirrep_; ++h){
            if(nsopi_[h] == 0) continue;

            double **temp = block_matrix(nsopi_[h], nsopi_[h]);
            /*
             * Backtransform the Tau matrices to the AO basis: soTau = C moTau Ct
             * Attn: the forward MO->AO transformation would be: moTau = Ct S soTau S C
             */
            double **paOccC = aocc_c_->pointer(h);
            double **pbOccC = bocc_c_->pointer(h);
            double **paVirC = avir_c_->pointer(h);
            double **pbVirC = bvir_c_->pointer(h);
            double **pa_tau_ = tau_so_a_->pointer(h);
            double **pb_tau_ = tau_so_b_->pointer(h);

            // Alpha occupied
            if(naoccpi_[h] && nsopi_[h]){
                C_DGEMM('n', 'n', nsopi_[h], naoccpi_[h], naoccpi_[h], 1.0, paOccC[0], naoccpi_[h],
                        T_OO.matrix[h][0], naoccpi_[h], 0.0, temp[0], nsopi_[h]);
                C_DGEMM('n', 't', nsopi_[h], nsopi_[h], naoccpi_[h], 1.0, temp[0], nsopi_[h],
                        paOccC[0], naoccpi_[h], 1.0, pa_tau_[0], nsopi_[h]);
            }
            // Check:
    //        tau_so_a_->print();

            // Beta occupied
            if(nboccpi_[h] && nsopi_[h]){
                C_DGEMM('n', 'n', nsopi_[h], nboccpi_[h], nboccpi_[h], 1.0, pbOccC[0], nboccpi_[h],
                        T_oo.matrix[h][0], nboccpi_[h], 0.0, temp[0], nsopi_[h]);
                C_DGEMM('n', 't', nsopi_[h], nsopi_[h], nboccpi_[h], 1.0, temp[0], nsopi_[h],
                        pbOccC[0], nboccpi_[h], 1.0, pb_tau_[0], nsopi_[h]);
            }
            // Alpha virtual
            if(navirpi_[h] && nsopi_[h]){
                C_DGEMM('n', 'n', nsopi_[h], navirpi_[h], navirpi_[h], 1.0, paVirC[0], navirpi_[h],
                        T_VV.matrix[h][0], navirpi_[h], 0.0, temp[0], nsopi_[h]);
                C_DGEMM('n', 't', nsopi_[h], nsopi_[h], navirpi_[h], 1.0, temp[0], nsopi_[h],
                        paVirC[0], navirpi_[h], 1.0, pa_tau_[0], nsopi_[h]);
            }
            // Check:
    //        tau_so_a_->print();

            // Beta virtual
            if(nbvirpi_[h] && nsopi_[h]){
                C_DGEMM('n', 'n', nsopi_[h], nbvirpi_[h], nbvirpi_[h], 1.0, pbVirC[0], nbvirpi_[h],
                        T_vv.matrix[h][0], nbvirpi_[h], 0.0, temp[0], nsopi_[h]);
                C_DGEMM('n', 't', nsopi_[h], nsopi_[h], nbvirpi_[h], 1.0, temp[0], nsopi_[h],
                        pbVirC[0], nbvirpi_[h], 1.0, pb_tau_[0], nsopi_[h]);
            }

            free_block(temp);
        }
        // Check:
    //    tau_so_a_->print();
    //    tau_so_b_->print();

        global_dpd_->file2_close(&T_OO);
        global_dpd_->file2_close(&T_oo);
        global_dpd_->file2_close(&T_VV);
        global_dpd_->file2_close(&T_vv);

        dcft_timer_off("DCFTSolver::transform_tau()");
    }

}

/**
 * Prints the occupation numbers from the OPDM
 */
void
DCFTSolver::print_opdm()
{
    // If the system is closed-shell, ...
    if(closedshell_){
        dpdbuf4 L1, L2;
        dpdfile2 T_OO, T_oo, T_VV, T_vv;
        global_dpd_->file2_init(&T_OO, PSIF_DCFT_DPD, 0, ID('O'), ID('O'), "Tau <O|O>");
//        global_dpd_->file2_init(&T_oo, PSIF_DCFT_DPD, 0, ID('o'), ID('o'), "Tau <o|o>");
        global_dpd_->file2_init(&T_VV, PSIF_DCFT_DPD, 0, ID('V'), ID('V'), "Tau <V|V>");
//        global_dpd_->file2_init(&T_vv, PSIF_DCFT_DPD, 0, ID('v'), ID('v'), "Tau <v|v>");

        global_dpd_->file2_mat_init(&T_OO);
//        global_dpd_->file2_mat_init(&T_oo);
        global_dpd_->file2_mat_init(&T_VV);
//        global_dpd_->file2_mat_init(&T_vv);

        global_dpd_->file2_mat_rd(&T_OO);
//        global_dpd_->file2_mat_rd(&T_oo);
        global_dpd_->file2_mat_rd(&T_VV);
//        global_dpd_->file2_mat_rd(&T_vv);

        std::vector<std::pair<double, int> > aPairs;

        for(int h = 0; h < nirrep_; ++h){
            for(int row = 0; row < T_OO.params->coltot[h]; ++row)
                aPairs.push_back(std::make_pair(1.0 + T_OO.matrix[h][row][row], h));
            for(int row = 0; row < T_VV.params->coltot[h]; ++row)
                aPairs.push_back(std::make_pair(T_VV.matrix[h][row][row], h));
//            for(int row = 0; row < T_oo.params->coltot[h]; ++row)
//                bPairs.push_back(std::make_pair(1.0 + T_oo.matrix[h][row][row], h));
//            for(int row = 0; row < T_vv.params->coltot[h]; ++row)
//                bPairs.push_back(std::make_pair(T_vv.matrix[h][row][row], h));
        }

        std::vector<std::pair<double, int> > bPairs(aPairs);

        global_dpd_->file2_close(&T_OO);
//        global_dpd_->file2_close(&T_oo);
        global_dpd_->file2_close(&T_VV);
//        global_dpd_->file2_close(&T_vv);

        sort(aPairs.begin(), aPairs.end(), greater<std::pair<double, int> >());
        sort(bPairs.begin(), bPairs.end(), greater<std::pair<double, int> >());

        int *aIrrepCount = init_int_array(nirrep_);
        int *bIrrepCount = init_int_array(nirrep_);
        char **irrepLabels = molecule_->irrep_labels();

        outfile->Printf( "\n\tOrbital occupations:\n\t\tAlpha occupied orbitals\n\t\t");
        for (int i = 0, count = 0; i < nalpha_; ++i, ++count) {
            int irrep = aPairs[i].second;
            outfile->Printf( "%4d%-4s%11.4f  ", ++aIrrepCount[irrep], irrepLabels[irrep], aPairs[i].first);
            if (count % 4 == 3 && i != nalpha_)
                outfile->Printf( "\n\t\t");
        }
        outfile->Printf( "\n\n\t\tBeta occupied orbitals\n\t\t");
        for (int i = 0, count = 0; i < nbeta_; ++i, ++count) {
            int irrep = bPairs[i].second;
            outfile->Printf( "%4d%-4s%11.4f  ", ++bIrrepCount[irrep], irrepLabels[irrep], bPairs[i].first);
            if (count % 4 == 3 && i != nbeta_)
                outfile->Printf( "\n\t\t");
        }
        outfile->Printf( "\n\n\t\tAlpha virtual orbitals\n\t\t");
        for (int i = nalpha_, count = 0; i < nmo_; ++i, ++count) {
            int irrep = aPairs[i].second;
            outfile->Printf( "%4d%-4s%11.4f  ", ++aIrrepCount[irrep], irrepLabels[irrep], aPairs[i].first);
            if (count % 4 == 3 && i != nmo_)
                outfile->Printf( "\n\t\t");
        }
        outfile->Printf( "\n\n\t\tBeta virtual orbitals\n\t\t");
        for (int i = nbeta_, count = 0; i < nmo_; ++i, ++count) {
            int irrep = bPairs[i].second;
            outfile->Printf( "%4d%-4s%11.4f  ", ++bIrrepCount[irrep], irrepLabels[irrep], bPairs[i].first);
            if (count % 4 == 3 && i != nmo_)
                outfile->Printf( "\n\t\t");
        }
        outfile->Printf( "\n\n");
        for (int h = 0; h < nirrep_; ++h)
            free(irrepLabels[h]);
        free(irrepLabels);
        free(aIrrepCount);
        free(bIrrepCount);

    }
    // If the system is open-shell, ...
    else{
        dpdbuf4 L1, L2;
        dpdfile2 T_OO, T_oo, T_VV, T_vv;
        global_dpd_->file2_init(&T_OO, PSIF_DCFT_DPD, 0, ID('O'), ID('O'), "Tau <O|O>");
        global_dpd_->file2_init(&T_oo, PSIF_DCFT_DPD, 0, ID('o'), ID('o'), "Tau <o|o>");
        global_dpd_->file2_init(&T_VV, PSIF_DCFT_DPD, 0, ID('V'), ID('V'), "Tau <V|V>");
        global_dpd_->file2_init(&T_vv, PSIF_DCFT_DPD, 0, ID('v'), ID('v'), "Tau <v|v>");

        global_dpd_->file2_mat_init(&T_OO);
        global_dpd_->file2_mat_init(&T_oo);
        global_dpd_->file2_mat_init(&T_VV);
        global_dpd_->file2_mat_init(&T_vv);

        global_dpd_->file2_mat_rd(&T_OO);
        global_dpd_->file2_mat_rd(&T_oo);
        global_dpd_->file2_mat_rd(&T_VV);
        global_dpd_->file2_mat_rd(&T_vv);

        std::vector<std::pair<double, int> > aPairs;
        std::vector<std::pair<double, int> > bPairs;

        for(int h = 0; h < nirrep_; ++h){
            for(int row = 0; row < T_OO.params->coltot[h]; ++row)
                aPairs.push_back(std::make_pair(1.0 + T_OO.matrix[h][row][row], h));
            for(int row = 0; row < T_VV.params->coltot[h]; ++row)
                aPairs.push_back(std::make_pair(T_VV.matrix[h][row][row], h));
            for(int row = 0; row < T_oo.params->coltot[h]; ++row)
                bPairs.push_back(std::make_pair(1.0 + T_oo.matrix[h][row][row], h));
            for(int row = 0; row < T_vv.params->coltot[h]; ++row)
                bPairs.push_back(std::make_pair(T_vv.matrix[h][row][row], h));
        }
        global_dpd_->file2_close(&T_OO);
        global_dpd_->file2_close(&T_oo);
        global_dpd_->file2_close(&T_VV);
        global_dpd_->file2_close(&T_vv);

        sort(aPairs.begin(), aPairs.end(), greater<std::pair<double, int> >());
        sort(bPairs.begin(), bPairs.end(), greater<std::pair<double, int> >());

        int *aIrrepCount = init_int_array(nirrep_);
        int *bIrrepCount = init_int_array(nirrep_);
        char **irrepLabels = molecule_->irrep_labels();

        outfile->Printf( "\n\tOrbital occupations:\n\t\tAlpha occupied orbitals\n\t\t");
        for (int i = 0, count = 0; i < nalpha_; ++i, ++count) {
            int irrep = aPairs[i].second;
            outfile->Printf( "%4d%-4s%11.4f  ", ++aIrrepCount[irrep], irrepLabels[irrep], aPairs[i].first);
            if (count % 4 == 3 && i != nalpha_)
                outfile->Printf( "\n\t\t");
        }
        outfile->Printf( "\n\n\t\tBeta occupied orbitals\n\t\t");
        for (int i = 0, count = 0; i < nbeta_; ++i, ++count) {
            int irrep = bPairs[i].second;
            outfile->Printf( "%4d%-4s%11.4f  ", ++bIrrepCount[irrep], irrepLabels[irrep], bPairs[i].first);
            if (count % 4 == 3 && i != nbeta_)
                outfile->Printf( "\n\t\t");
        }
        outfile->Printf( "\n\n\t\tAlpha virtual orbitals\n\t\t");
        for (int i = nalpha_, count = 0; i < nmo_; ++i, ++count) {
            int irrep = aPairs[i].second;
            outfile->Printf( "%4d%-4s%11.4f  ", ++aIrrepCount[irrep], irrepLabels[irrep], aPairs[i].first);
            if (count % 4 == 3 && i != nmo_)
                outfile->Printf( "\n\t\t");
        }
        outfile->Printf( "\n\n\t\tBeta virtual orbitals\n\t\t");
        for (int i = nbeta_, count = 0; i < nmo_; ++i, ++count) {
            int irrep = bPairs[i].second;
            outfile->Printf( "%4d%-4s%11.4f  ", ++bIrrepCount[irrep], irrepLabels[irrep], bPairs[i].first);
            if (count % 4 == 3 && i != nmo_)
                outfile->Printf( "\n\t\t");
        }
        outfile->Printf( "\n\n");
        for (int h = 0; h < nirrep_; ++h)
            free(irrepLabels[h]);
        free(irrepLabels);
        free(aIrrepCount);
        free(bIrrepCount);

    }
}

void
DCFTSolver::refine_tau() {

    // If the system is closed-shell, then ...
    if(closedshell_){
        dcft_timer_on("DCFTSolver::refine_tau()");

        // Read MO-basis Tau from disk into the memory
        dpdfile2 T_OO, T_VV;

        // Iteratively compute the exact Tau

        SharedMatrix aocc_tau_old(new Matrix("MO basis Tau (Alpha Occupied, old)", nirrep_, naoccpi_, naoccpi_));
        SharedMatrix avir_tau_old(new Matrix("MO basis Tau (Alpha Virtual, old)", nirrep_, navirpi_, navirpi_));
        SharedMatrix aocc_d(new Matrix("Non-idempotency of OPDM (Alpha Occupied, old)", nirrep_, naoccpi_, naoccpi_));
        SharedMatrix avir_d(new Matrix("Non-idempotency of OPDM (Alpha Virtual, old)", nirrep_, navirpi_, navirpi_));

        bool converged = false;
        bool failed = false;
        int cycle = 0;

        // Copy approximate Tau as the non-idempotency of OPDM
        aocc_d->copy(aocc_tau_);
        avir_d->copy(avir_tau_);

        while(!converged && !failed){

            // Save old tau from previous iteration
            aocc_tau_old->copy(aocc_tau_);
            avir_tau_old->copy(avir_tau_);

            // Tau_ij = d_ij
            // Tau_ab = -d_ab
            aocc_tau_->copy(aocc_d);
            avir_tau_->copy(avir_d);

            // Tau_ij -= Tau_ik * Tau_kj
            // Tau_ab += Tau_ac * Tau_cb
            aocc_tau_->gemm(false, false, -1.0, aocc_tau_old, aocc_tau_old, 1.0);
            avir_tau_->gemm(false, false, 1.0, avir_tau_old, avir_tau_old, 1.0);

            // Compute RMS
            aocc_tau_old->subtract(aocc_tau_);
            avir_tau_old->subtract(avir_tau_);

            double rms = aocc_tau_old->rms();
            rms += avir_tau_old->rms();
            rms *= 2.0;

            converged = (rms < cumulant_threshold_);
            failed    = (++cycle == maxiter_);

            if (print_ > 2) outfile->Printf( "\t Exact Tau Iterations: %-3d %20.12f\n", cycle, rms);
//            if (print_ > 0) outfile->Printf( "\t Exact Tau Iterations: %-3d %20.12f\n", cycle, rms);

        } // end of macroiterations

        // Test the trace of Tau
        // double trace = aocc_tau_->trace() + avir_tau_->trace() + bocc_tau_->trace() + bvir_tau_->trace();
        // outfile->Printf( "\t Trace of Tau: %8.7e\n", trace);

        // If exact tau iterations failed, throw a message about it and compute it non-iteratively
        if (failed) {
            outfile->Printf( "\t Exact Tau didn't converge. Evaluating it non-iteratively\n");
            // Set old tau matrices to identity
            aocc_tau_old->identity();
            avir_tau_old->identity();
            // Scale the non-idempotency elements
            aocc_d->scale(4.0);
            avir_d->scale(-4.0);
            // Add them to the old tau
            aocc_tau_old->add(aocc_d);
            avir_tau_old->add(avir_d);
            // Zero out new tau
            aocc_tau_->zero();
            avir_tau_->zero();
            // Diagonalize and take a square root
            SharedMatrix aocc_evecs(new Matrix("Eigenvectors (Alpha Occupied)", nirrep_, naoccpi_, naoccpi_));
            SharedMatrix avir_evecs(new Matrix("Eigenvectors (Alpha Virtual)", nirrep_, navirpi_, navirpi_));
            SharedVector aocc_evals(new Vector("Eigenvalues (Alpha Occupied)", nirrep_, naoccpi_));
            SharedVector avir_evals(new Vector("Eigenvalues (Alpha Virtual)", nirrep_, navirpi_));
            aocc_tau_old->diagonalize(aocc_evecs, aocc_evals);
            avir_tau_old->diagonalize(avir_evecs, avir_evals);

            for(int h = 0; h < nirrep_; ++h){
                if(nsopi_[h] == 0) continue;

                // Alpha occupied
                for(int p = 0 ; p < naoccpi_[h]; ++p) aocc_tau_->set(h, p, p, (-1.0 + sqrt(aocc_evals->get(h, p))) / 2.0);

                // Alpha virtual
                for(int p = 0 ; p < navirpi_[h]; ++p) avir_tau_->set(h, p, p, (1.0 - sqrt(avir_evals->get(h, p))) / 2.0);
            }

            // Back-transform the diagonal Tau to the original basis
            aocc_tau_->back_transform(aocc_evecs);
            avir_tau_->back_transform(avir_evecs);
        }

        // Copy Tau_alpha to Tau_beta
        bocc_tau_->copy(aocc_tau_);
        bvir_tau_->copy(avir_tau_);

        // Write the exact tau back to disk

        global_dpd_->file2_init(&T_OO, PSIF_DCFT_DPD, 0, ID('O'), ID('O'), "Tau <O|O>");
        global_dpd_->file2_init(&T_VV, PSIF_DCFT_DPD, 0, ID('V'), ID('V'), "Tau <V|V>");

        global_dpd_->file2_mat_init(&T_OO);
        global_dpd_->file2_mat_init(&T_VV);

        for(int h = 0; h < nirrep_; ++h){
            for(int i = 0; i < naoccpi_[h]; ++i){
                for(int j = 0; j < naoccpi_[h]; ++j){
                    T_OO.matrix[h][i][j] = aocc_tau_->get(h, i, j);
                }
            }
            for(int a = 0; a < navirpi_[h]; ++a){
                for(int b = 0; b < navirpi_[h]; ++b){
                    T_VV.matrix[h][a][b] = avir_tau_->get(h, a, b);
                }
            }
        }

        global_dpd_->file2_mat_wrt(&T_OO);
        global_dpd_->file2_mat_wrt(&T_VV);

        // Copy exact tau_alpha to tau_beta in disk
//        global_dpd_->file2_copy(&T_OO, PSIF_DCFT_DPD, "Tau <o|o>");
//        global_dpd_->file2_copy(&T_VV, PSIF_DCFT_DPD, "Tau <v|v>");

        global_dpd_->file2_close(&T_OO);
        global_dpd_->file2_close(&T_VV);

        dcft_timer_off("DCFTSolver::refine_tau()");

    }
    // If the system is open-shell, then ...
    else{
        // Read MO-basis Tau from disk into the memory
        dpdfile2 T_OO, T_oo, T_VV, T_vv;

        // Iteratively compute the exact Tau

        SharedMatrix aocc_tau_old(new Matrix("MO basis Tau (Alpha Occupied, old)", nirrep_, naoccpi_, naoccpi_));
        SharedMatrix bocc_tau_old(new Matrix("MO basis Tau (Beta Occupied, old)", nirrep_, nboccpi_, nboccpi_));
        SharedMatrix avir_tau_old(new Matrix("MO basis Tau (Alpha Virtual, old)", nirrep_, navirpi_, navirpi_));
        SharedMatrix bvir_tau_old(new Matrix("MO basis Tau (Beta Virtual, old)", nirrep_, nbvirpi_, nbvirpi_));
        SharedMatrix aocc_d(new Matrix("Non-idempotency of OPDM (Alpha Occupied, old)", nirrep_, naoccpi_, naoccpi_));
        SharedMatrix bocc_d(new Matrix("Non-idempotency of OPDM (Beta Occupied, old)", nirrep_, nboccpi_, nboccpi_));
        SharedMatrix avir_d(new Matrix("Non-idempotency of OPDM (Alpha Virtual, old)", nirrep_, navirpi_, navirpi_));
        SharedMatrix bvir_d(new Matrix("Non-idempotency of OPDM (Beta Virtual, old)", nirrep_, nbvirpi_, nbvirpi_));

        bool converged = false;
        bool failed = false;
        int cycle = 0;

        // Copy approximate Tau as the non-idempotency of OPDM
        aocc_d->copy(aocc_tau_);
        avir_d->copy(avir_tau_);
        bocc_d->copy(bocc_tau_);
        bvir_d->copy(bvir_tau_);

        while(!converged && !failed){

            // Save old tau from previous iteration
            aocc_tau_old->copy(aocc_tau_);
            avir_tau_old->copy(avir_tau_);
            bocc_tau_old->copy(bocc_tau_);
            bvir_tau_old->copy(bvir_tau_);

            // Tau_ij = d_ij
            // Tau_ab = -d_ab
            aocc_tau_->copy(aocc_d);
            avir_tau_->copy(avir_d);
            bocc_tau_->copy(bocc_d);
            bvir_tau_->copy(bvir_d);

            // Tau_ij -= Tau_ik * Tau_kj
            // Tau_ab += Tau_ac * Tau_cb
            aocc_tau_->gemm(false, false, -1.0, aocc_tau_old, aocc_tau_old, 1.0);
            avir_tau_->gemm(false, false, 1.0, avir_tau_old, avir_tau_old, 1.0);
            bocc_tau_->gemm(false, false, -1.0, bocc_tau_old, bocc_tau_old, 1.0);
            bvir_tau_->gemm(false, false, 1.0, bvir_tau_old, bvir_tau_old, 1.0);

            // Compute RMS
            aocc_tau_old->subtract(aocc_tau_);
            avir_tau_old->subtract(avir_tau_);
            bocc_tau_old->subtract(bocc_tau_);
            bvir_tau_old->subtract(bvir_tau_);

            double rms = aocc_tau_old->rms();
            rms += avir_tau_old->rms();
            rms += bocc_tau_old->rms();
            rms += bvir_tau_old->rms();

            converged = (rms < cumulant_threshold_);
            failed    = (++cycle == maxiter_);

            if (print_ > 2) outfile->Printf( "\t Exact Tau Iterations: %-3d %20.12f\n", cycle, rms);

        } // end of macroiterations

        // Test the trace of Tau
        // double trace = aocc_tau_->trace() + avir_tau_->trace() + bocc_tau_->trace() + bvir_tau_->trace();
        // outfile->Printf( "\t Trace of Tau: %8.7e\n", trace);

        // If exact tau iterations failed, throw a message about it and compute it non-iteratively
        if (failed) {
            outfile->Printf( "\t Exact Tau didn't converge. Evaluating it non-iteratively\n");
            // Set old tau matrices to identity
            aocc_tau_old->identity();
            bocc_tau_old->identity();
            avir_tau_old->identity();
            bvir_tau_old->identity();
            // Scale the non-idempotency elements
            aocc_d->scale(4.0);
            bocc_d->scale(4.0);
            avir_d->scale(-4.0);
            bvir_d->scale(-4.0);
            // Add them to the old tau
            aocc_tau_old->add(aocc_d);
            bocc_tau_old->add(bocc_d);
            avir_tau_old->add(avir_d);
            bvir_tau_old->add(bvir_d);
            // Zero out new tau
            aocc_tau_->zero();
            avir_tau_->zero();
            bocc_tau_->zero();
            bvir_tau_->zero();
            // Diagonalize and take a square root
            SharedMatrix aocc_evecs(new Matrix("Eigenvectors (Alpha Occupied)", nirrep_, naoccpi_, naoccpi_));
            SharedMatrix bocc_evecs(new Matrix("Eigenvectors (Beta Occupied)", nirrep_, nboccpi_, nboccpi_));
            SharedMatrix avir_evecs(new Matrix("Eigenvectors (Alpha Virtual)", nirrep_, navirpi_, navirpi_));
            SharedMatrix bvir_evecs(new Matrix("Eigenvectors (Beta Virtual)", nirrep_, nbvirpi_, nbvirpi_));
            SharedVector aocc_evals(new Vector("Eigenvalues (Alpha Occupied)", nirrep_, naoccpi_));
            SharedVector bocc_evals(new Vector("Eigenvalues (Beta Occupied)", nirrep_, nboccpi_));
            SharedVector avir_evals(new Vector("Eigenvalues (Alpha Virtual)", nirrep_, navirpi_));
            SharedVector bvir_evals(new Vector("Eigenvalues (Beta Virtual)", nirrep_, nbvirpi_));
            aocc_tau_old->diagonalize(aocc_evecs, aocc_evals);
            bocc_tau_old->diagonalize(bocc_evecs, bocc_evals);
            avir_tau_old->diagonalize(avir_evecs, avir_evals);
            bvir_tau_old->diagonalize(bvir_evecs, bvir_evals);

            for(int h = 0; h < nirrep_; ++h){
                if(nsopi_[h] == 0) continue;

                // Alpha occupied
                for(int p = 0 ; p < naoccpi_[h]; ++p) aocc_tau_->set(h, p, p, (-1.0 + sqrt(aocc_evals->get(h, p))) / 2.0);

                // Beta occupied
                for(int p = 0 ; p < nboccpi_[h]; ++p) bocc_tau_->set(h, p, p, (-1.0 + sqrt(bocc_evals->get(h, p))) / 2.0);

                // Alpha virtual
                for(int p = 0 ; p < navirpi_[h]; ++p) avir_tau_->set(h, p, p, (1.0 - sqrt(avir_evals->get(h, p))) / 2.0);

                // Beta virtual
                for(int p = 0 ; p < nbvirpi_[h]; ++p) bvir_tau_->set(h, p, p, (1.0 - sqrt(bvir_evals->get(h, p))) / 2.0);
            }

            // Back-transform the diagonal Tau to the original basis
            aocc_tau_->back_transform(aocc_evecs);
            bocc_tau_->back_transform(bocc_evecs);
            avir_tau_->back_transform(avir_evecs);
            bvir_tau_->back_transform(bvir_evecs);
        }

        // Write the exact tau back to disk

        global_dpd_->file2_init(&T_OO, PSIF_DCFT_DPD, 0, ID('O'), ID('O'), "Tau <O|O>");
        global_dpd_->file2_init(&T_oo, PSIF_DCFT_DPD, 0, ID('o'), ID('o'), "Tau <o|o>");
        global_dpd_->file2_init(&T_VV, PSIF_DCFT_DPD, 0, ID('V'), ID('V'), "Tau <V|V>");
        global_dpd_->file2_init(&T_vv, PSIF_DCFT_DPD, 0, ID('v'), ID('v'), "Tau <v|v>");

        global_dpd_->file2_mat_init(&T_OO);
        global_dpd_->file2_mat_init(&T_oo);
        global_dpd_->file2_mat_init(&T_VV);
        global_dpd_->file2_mat_init(&T_vv);

        for(int h = 0; h < nirrep_; ++h){
            for(int i = 0; i < naoccpi_[h]; ++i){
                for(int j = 0; j < naoccpi_[h]; ++j){
                    T_OO.matrix[h][i][j] = aocc_tau_->get(h, i, j);
                }
            }
            for(int a = 0; a < navirpi_[h]; ++a){
                for(int b = 0; b < navirpi_[h]; ++b){
                    T_VV.matrix[h][a][b] = avir_tau_->get(h, a, b);
                }
            }
            for(int i = 0; i < nboccpi_[h]; ++i){
                for(int j = 0; j < nboccpi_[h]; ++j){
                    T_oo.matrix[h][i][j] = bocc_tau_->get(h, i, j);
                }
            }
            for(int a = 0; a < nbvirpi_[h]; ++a){
                for(int b = 0; b < nbvirpi_[h]; ++b){
                    T_vv.matrix[h][a][b] = bvir_tau_->get(h, a, b);
                }
            }
        }

        global_dpd_->file2_mat_wrt(&T_OO);
        global_dpd_->file2_mat_wrt(&T_oo);
        global_dpd_->file2_mat_wrt(&T_VV);
        global_dpd_->file2_mat_wrt(&T_vv);

        global_dpd_->file2_close(&T_OO);
        global_dpd_->file2_close(&T_oo);
        global_dpd_->file2_close(&T_VV);
        global_dpd_->file2_close(&T_vv);

    }
}

}} // Namespaces


