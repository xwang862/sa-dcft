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
#include <psifiles.h>
#include <libpsio/psio.hpp>
#include <libtrans/integraltransform.h>
#include "defines.h"

namespace psi{ namespace dfdcft{

/**
 * Uses the intermediates to compute the energy
 */
void
DCFTSolver::compute_dcft_energy()
{

    // If the system is closed-shell, then...
    if(closedshell_){
        /*
         * E = lambda_IjAb * ( 2 M_IjAb - M_JiAb )
         * where M_IjAb = G_IjAb + gbar_IjAb
         */
        dcft_timer_on("DCFTSolver::compute_dcft_energy()");

        dpdbuf4 L, G, M, temp;
        double E_test;
        psio_->open(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);

//        global_dpd_->buf4_init(&L, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[V,v]"),
//                               ID("[O,o]"), ID("[V,v]"), 0, "Lambda <Oo|Vv>");
        global_dpd_->buf4_init(&L, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                               ID("[O,O]"), ID("[V,V]"), 0, "Lambda SF <OO|VV>");

        // M_IjAb = G_IjAb
//        global_dpd_->buf4_init(&G, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[V,v]"),
//                               ID("[O,o]"), ID("[V,v]"), 0, "G <Oo|Vv>");
        global_dpd_->buf4_init(&G, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                               ID("[O,O]"), ID("[V,V]"), 0, "G <OO|VV>");
//        global_dpd_->buf4_copy(&G, PSIF_DCFT_DPD, "M <Oo|Vv>");
        global_dpd_->buf4_copy(&G, PSIF_DCFT_DPD, "M <OO|VV>");
        global_dpd_->buf4_close(&G);
        // M_IjAb += gbar_IjAb
//        global_dpd_->buf4_init(&M, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[V,v]"),
//                               ID("[O,o]"), ID("[V,v]"), 0, "M <Oo|Vv>");
        global_dpd_->buf4_init(&M, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                               ID("[O,O]"), ID("[V,V]"), 0, "M <OO|VV>");
//        global_dpd_->buf4_init(&temp, PSIF_LIBTRANS_DPD, 0, ID("[O,o]"), ID("[V,v]"),
//                               ID("[O,o]"), ID("[V,v]"), 0, "MO Ints <Oo|Vv>");
        global_dpd_->buf4_init(&temp, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                               ID("[O,O]"), ID("[V,V]"), 0, "MO Ints <OO|VV>");
        dpd_buf4_add(&M, &temp, 1.0);
        global_dpd_->buf4_close(&M);
        global_dpd_->buf4_close(&temp);
//        // M(temp) <Oo|Vv> = M_JiAb
//        global_dpd_->buf4_init(&M, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[V,v]"),
//                               ID("[O,o]"), ID("[V,v]"), 0, "M <Oo|Vv>");
//        global_dpd_->buf4_sort(&M, PSIF_DCFT_DPD, qprs, ID("[O,o]"), ID("[V,v]"), "M(temp) <Oo|Vv>");
//        global_dpd_->buf4_close(&M);
//        // Energy = lambda_IjAb * ( 2 M_IjAb - M_JiAb )
//        global_dpd_->buf4_init(&M, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[V,v]"),
//                               ID("[O,o]"), ID("[V,v]"), 0, "M <Oo|Vv>");
//        global_dpd_->buf4_init(&temp, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[V,v]"),
//                               ID("[O,o]"), ID("[V,v]"), 0, "M(temp) <Oo|Vv>");
//        global_dpd_->buf4_scm(&M, 2.0);
//        dpd_buf4_add(&M, &temp, -1.0);

        global_dpd_->buf4_init(&M, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                               ID("[O,O]"), ID("[V,V]"), 1, "M <OO|VV>");
        global_dpd_->buf4_copy(&M, PSIF_DCFT_DPD, "M(temp) <OO|VV>");
        global_dpd_->buf4_close(&M);
        global_dpd_->buf4_init(&M, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                               ID("[O,O]"), ID("[V,V]"), 0, "M(temp) <OO|VV>");
        global_dpd_->buf4_init(&temp, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                               ID("[O,O]"), ID("[V,V]"), 0, "M <OO|VV>");
        dpd_buf4_add(&M, &temp, 1.0);

        E_test = global_dpd_->buf4_dot(&L, &M);
        global_dpd_->buf4_close(&M);
        global_dpd_->buf4_close(&temp);
        global_dpd_->buf4_close(&L);

        psio_->close(PSIF_LIBTRANS_DPD, 1);

        lambda_energy_ = E_test;

        // Check:
//        outfile->Printf( "\t* %d Lambda Energy = %20.16f\n", macro_cycle, lambda_energy_);

        dcft_timer_off("DCFTSolver::compute_dcft_energy()");
    }
    // If the system is open-shell, then...
    else{
        dcft_timer_on("DCFTSolver::compute_dcft_energy()");

        dpdbuf4 L, G, I;
        double eGaa, eGab, eGbb, eIaa, eIab, eIbb;

        psio_->open(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);

        // E += 1/4 L_IJAB G_IJAB
        global_dpd_->buf4_init(&L, PSIF_DCFT_DPD, 0, ID("[O>O]-"), ID("[V>V]-"),
                      ID("[O>O]-"), ID("[V>V]-"), 0, "Lambda <OO|VV>");
        global_dpd_->buf4_init(&G, PSIF_DCFT_DPD, 0, ID("[O>O]-"), ID("[V>V]-"),
                      ID("[O>O]-"), ID("[V>V]-"), 0, "G <OO|VV>");
        eGaa = global_dpd_->buf4_dot(&G, &L);
        global_dpd_->buf4_close(&G);

        // E += 1/4 gbar_IJAB L_IJAB
        global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O>O]-"), ID("[V>V]-"),
                      ID("[O,O]"), ID("[V,V]"), 1, "MO Ints <OO|VV>");
        eIaa = global_dpd_->buf4_dot(&I, &L);
        global_dpd_->buf4_close(&I);
        global_dpd_->buf4_close(&L);

        // E += L_IjAb G_IjAb
        global_dpd_->buf4_init(&L, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                      ID("[O,o]"), ID("[V,v]"), 0, "Lambda <Oo|Vv>");
        global_dpd_->buf4_init(&G, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                      ID("[O,o]"), ID("[V,v]"), 0, "G <Oo|Vv>");
        eGab =  global_dpd_->buf4_dot(&G, &L);
        global_dpd_->buf4_close(&G);

        // E += gbar_IjAb L_IjAb
        global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                      ID("[O,o]"), ID("[V,v]"), 0, "MO Ints <Oo|Vv>");
        eIab = global_dpd_->buf4_dot(&I, &L);
        global_dpd_->buf4_close(&I);
        global_dpd_->buf4_close(&L);

        // E += 1/4 L_ijab G_ijab
        global_dpd_->buf4_init(&L, PSIF_DCFT_DPD, 0, ID("[o>o]-"), ID("[v>v]-"),
                      ID("[o>o]-"), ID("[v>v]-"), 0, "Lambda <oo|vv>");
        global_dpd_->buf4_init(&G, PSIF_DCFT_DPD, 0, ID("[o>o]-"), ID("[v>v]-"),
                      ID("[o>o]-"), ID("[v>v]-"), 0, "G <oo|vv>");
        eGbb = global_dpd_->buf4_dot(&G, &L);
        global_dpd_->buf4_close(&G);


        // E += 1/4 gbar_ijab L_ijab
        global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[o>o]-"), ID("[v>v]-"),
                      ID("[o,o]"), ID("[v,v]"), 1, "MO Ints <oo|vv>");
        eIbb = global_dpd_->buf4_dot(&I, &L);
        global_dpd_->buf4_close(&I);
        global_dpd_->buf4_close(&L);
        psio_->close(PSIF_LIBTRANS_DPD, 1);

    #if PRINT_ENERGY_COMPONENTS
        outfile->Printf( "\tAA G Energy = %20.12f\n", eGaa);
        outfile->Printf( "\tAB G Energy = %20.12f\n", eGab);
        outfile->Printf( "\tBB G Energy = %20.12f\n", eGbb);
        outfile->Printf( "\tAA I Energy = %20.12f\n", eIaa);
        outfile->Printf( "\tAB I Energy = %20.12f\n", eIab);
        outfile->Printf( "\tBB I Energy = %20.12f\n", eIbb);
        outfile->Printf( "\tTotal G Energy = %20.12f\n", eGaa + eGab + eGbb);
        outfile->Printf( "\tTotal I Energy = %20.12f\n", eIaa + eIab + eIbb);
    #endif

        lambda_energy_ = eGaa + eGab + eGbb + eIaa + eIab + eIbb;

        // Check:
//        outfile->Printf( "\t* %d Lambda Energy = %20.16f\n", macro_cycle, lambda_energy_);

        dcft_timer_off("DCFTSolver::compute_dcft_energy()");

    }
}

void
DCFTSolver::compute_cepa0_energy()
{
    dcft_timer_on("DCFTSolver::compute_dcft_energy()");

    psio_->open(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);

    /*
     * Compute the CEPA-0 correlation energy
     * E = 1/4 L_IJAB <IJ||AB>
     *        +L_IjAb <Ij|Ab>
     *    +1/4 L_ijab <ij||ab>
     */
    dpdbuf4 I, L;
    // Alpha - Alpha
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 1, "MO Ints <OO|VV>");
    global_dpd_->buf4_init(&L, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O>O]-"), ID("[V>V]-"), 0, "Lambda <OO|VV>");
    double eAA = 0.25 * global_dpd_->buf4_dot(&L, &I);
    global_dpd_->buf4_close(&I);
    global_dpd_->buf4_close(&L);

    // Alpha - Beta
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "MO Ints <Oo|Vv>");
    global_dpd_->buf4_init(&L, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "Lambda <Oo|Vv>");
    double eAB = global_dpd_->buf4_dot(&L, &I);
    global_dpd_->buf4_close(&I);
    global_dpd_->buf4_close(&L);

    // Beta - Beta
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 1, "MO Ints <oo|vv>");
    global_dpd_->buf4_init(&L, PSIF_DCFT_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o>o]-"), ID("[v>v]-"), 0, "Lambda <oo|vv>");
    double eBB = 0.25 * global_dpd_->buf4_dot(&L, &I);
    global_dpd_->buf4_close(&I);
    global_dpd_->buf4_close(&L);

    psio_->close(PSIF_LIBTRANS_DPD, 1);

    lambda_energy_ = eAA + eAB + eBB;

    dcft_timer_off("DCFTSolver::compute_dcft_energy()");
}

}} // Namespaces



