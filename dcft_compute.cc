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
#include <cmath>
#include <libdpd/dpd.h>
#include <libtrans/integraltransform.h>
#include <libdiis/diismanager.h>
#include <libpsio/psio.hpp>
#include <libpsio/psio.h>

#include "defines.h"

using namespace boost;

namespace psi{ namespace dfdcft{

/**
 * Computes the DCFT density matrix and energy
 */
double
DCFTSolver::compute_energy()
{
    if (Process::environment.wavefunction()->same_a_b_orbs()){
        outfile->Printf("\tTransform Type: Restricted.\n\n");
    }
    else outfile->Printf("\tTransform Type: Unrestricted.\n\n");

    if (same_a_b_orbs()){
        outfile->Printf("\tTransform Type: Restricted.\n\n");
    }
    else outfile->Printf("\tTransform Type: Unrestricted.\n\n");

    orbitalsDone_    = false;
    cumulantDone_ = false;
    densityConverged_ = false;
    energyConverged_ = false;

    // If closedshell_ = true, spin adaptation will be used
    if (options_.get_str("REFERENCE") == "UHF")
        closedshell_ = false;
    else if((same_a_b_dens() && same_a_b_orbs()) || options_.get_str("REFERENCE") == "RHF")
        closedshell_ = true;
    else
        closedshell_ = false;

    // Not implemented
    if (options_.get_str("DCFT_GUESS") == "DCFT" && closedshell_)
        throw FeatureNotImplemented("Spin-adapted RHF-reference ODC-12", "DCFT_GUESS = DCFT", __FILE__, __LINE__);

    // Perform SCF guess for the orbitals
    scf_guess();

    // Perform DF-MP2 guess for the cumulant
    if(0){
        outfile->Printf( "\n\tDF-MP2 guess is running... \n\n");
        dfmp2_guess();
    }

    // Perform MP2 guess for the cumulant
    mp2_guess();


    // Print out information about the job
    outfile->Printf( "\n\tDCFT Functional:    \t\t %s", options_.get_str("DCFT_FUNCTIONAL").c_str());
    outfile->Printf( "\n\tAlgorithm:          \t\t %s", options_.get_str("ALGORITHM").c_str());
    outfile->Printf( "\n\tAO-Basis Integrals: \t\t %s", options_.get_str("AO_BASIS").c_str());
    if (options_.get_str("ALGORITHM") == "QC") {
        outfile->Printf( "\n\tQC type:            \t\t %s", options_.get_str("QC_TYPE").c_str());
        outfile->Printf( "\n\tQC coupling:        \t\t %s", options_.get_bool("QC_COUPLING") ? "TRUE" : "FALSE");
    }

    // Things that are not implemented yet...
    if (options_.get_str("DERTYPE") == "FIRST" && (options_.get_str("DCFT_FUNCTIONAL") == "DC-12"))
        throw FeatureNotImplemented("DC-12 functional", "Analytic gradients", __FILE__, __LINE__);
    if (options_.get_str("AO_BASIS") == "DISK" && options_.get_str("DCFT_FUNCTIONAL") == "CEPA0")
        throw FeatureNotImplemented("CEPA0", "AO_BASIS = DISK", __FILE__, __LINE__);
    if (options_.get_str("AO_BASIS") == "DISK" && options_.get_str("ALGORITHM") == "QC" && options_.get_str("QC_TYPE") == "SIMULTANEOUS")
        throw FeatureNotImplemented("Simultaneous QC", "AO_BASIS = DISK", __FILE__, __LINE__);
    if (!(options_.get_str("ALGORITHM") == "TWOSTEP") && options_.get_str("DCFT_FUNCTIONAL") == "CEPA0")
        throw FeatureNotImplemented("CEPA0", "Requested DCFT algorithm", __FILE__, __LINE__);

    // Orbital-optimized stuff
    if (options_.get_str("ALGORITHM") == "TWOSTEP" && orbital_optimized_)
        throw PSIEXCEPTION("Two-step algorithm cannot be run for the orbital-optimized DCFT methods");

    // Choose a paricular algorithm and solve the equations
    if(options_.get_str("ALGORITHM") == "TWOSTEP") {
        run_twostep_dcft();
    }
    else if (options_.get_str("ALGORITHM") == "SIMULTANEOUS") {
        if (!orbital_optimized_) {
            run_simult_dcft();
        }
        else {
            run_simult_dcft_oo();
        }
    }
    else if (options_.get_str("ALGORITHM") == "QC") {
        run_qc_dcft();
    }
    else {
        throw PSIEXCEPTION("Unknown DCFT algoritm");
    }

    // If not converged -> Break
    if(!orbitalsDone_ || !cumulantDone_ || !densityConverged_)
        throw ConvergenceError<int>("DCFT", maxiter_, cumulant_threshold_, cumulant_convergence_, __FILE__, __LINE__);

    outfile->Printf( "\n\t*DCFT SCF Energy                                 = %20.15f\n", scf_energy_);
    outfile->Printf(   "\t*DCFT Lambda Energy                              = %20.15f\n", lambda_energy_);
    outfile->Printf(   "\t*DCFT Total Energy                               = %20.15f\n", new_total_energy_);

    Process::environment.globals["CURRENT ENERGY"] = new_total_energy_;
    Process::environment.globals["DCFT TOTAL ENERGY"] = new_total_energy_;
    Process::environment.globals["DCFT SCF ENERGY"] = scf_energy_;
    Process::environment.globals["DCFT LAMBDA ENERGY"] = lambda_energy_;

    if(!options_.get_bool("MO_RELAX")){
        outfile->Printf( "Warning!  The orbitals were not relaxed\n");
    }

    print_opdm();

    // Check
//    if(options_.get_bool("TPDM")) outfile->Printf( "\n\t TPDM = True \n", scf_energy_);
//    else outfile->Printf( "\n\t TPDM = False \n", scf_energy_);

    if(options_.get_bool("TPDM")) dump_density(); // By default, TPDM = FALSE
//    check_n_representability();

    // Compute the analytic gradients, if requested
    if(options_.get_str("DERTYPE") == "FIRST") {
        // Shut down the timers
        tstop();
        // Start the timers
        tstart();
        // Solve the response equations, compute relaxed OPDM and TPDM and dump them to disk
        compute_gradient();
    }

    if (Process::environment.wavefunction()->same_a_b_orbs()){
        outfile->Printf("\tTransform Type: Restricted.\n\n");
    }
    else outfile->Printf("\tTransform Type: Unrestricted.\n\n");

    if (Process::environment.wavefunction()->density_fitted()){
        outfile->Printf("\tDensity Fitting: Yes!\n\n");
    }
    else outfile->Printf("\tDensity Fitting: No!\n\n");

    // Free up memory and close files
    finalize();

    return(new_total_energy_);
}

void
DCFTSolver::run_twostep_dcft()
{
    // This is the two-step update - in each macro iteration, update the orbitals first, then update lambda
    // to self-consistency, until converged.  When lambda is converged and only one scf cycle is needed to reach
    // the desired cutoff, we're done

    // Check:
    macro_cycle = 0;

    int cycle = 0;
    outfile->Printf( "\n\n\t*=================================================================================*\n"
                         "\t* Cycle  RMS [F, Kappa]   RMS Lambda Error   delta E        Total Energy     DIIS *\n"
                         "\t*---------------------------------------------------------------------------------*\n");

    // Set up the DIIS manager for the density cumulant and SCF iterations
    old_ca_->copy(Ca_);
    old_cb_->copy(Cb_);
    // Save F0 = H + G * Kappa for the Fock intermediate update in lambda iterations
    moF0a_->copy(Fa_);
    moF0b_->copy(Fb_);
    moF0a_->transform(Ca_);
    moF0b_->transform(Cb_);
    // Check:
//    psio_->open(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);
//    dpdfile2 F_VV;
//    global_dpd_->file2_init(&F_VV, PSIF_LIBTRANS_DPD, 0, ID('V'), ID('V'), "F <V|V>");
//    global_dpd_->file2_print(&F_VV, "outfile");
//    global_dpd_->file2_close(&F_VV);
//    psio_->close(PSIF_LIBTRANS_DPD, 1);

    // Just so the correct value is printed in the first macro iteration
    orbitals_convergence_ = compute_scf_error_vector();
    // Start macro-iterations
    while((!orbitalsDone_ || !cumulantDone_) && cycle++ < maxiter_){
        outfile->Printf( "\t                          *** Macro Iteration %d ***\n"
                         "\tCumulant Iterations\n",cycle);
        // If it's the first iteration and the user requested to relax guess orbitals, then skip the density cumulant update
        if ((cycle != 1) || !options_.get_bool("RELAX_GUESS_ORBITALS")) { // The default RELAX_GUESS_ORBITALS is TRUE !
            run_twostep_dcft_cumulant_updates();
            // Check:
            macro_cycle++;
        }
        else outfile->Printf( "\tSkipping the cumulant update to relax guess orbitals\n");

        // Break if it's a CEPA0 computation
        if (options_.get_str("DCFT_FUNCTIONAL") == "CEPA0") {
            orbitalsDone_ = true;
            cumulantDone_ = true;
            densityConverged_ = true;
            break;
        }
        // Build new Tau from the density cumulant in the MO basis and transform it the SO basis
        build_tau();
        // Compute tau exactly if requested
        if (exact_tau_) {
            refine_tau();
        }
        transform_tau();
        run_twostep_dcft_orbital_updates();
    }

    // Check: PASSED
    if(0){
        /*
         * Test: if lambda_IJAB = lambda_IjAb - lambda_JiAb
         */
        dpdbuf4 L_ref, L_test, temp;
//        psio_->open(PSIF_DCFT_DPD, PSIO_OPEN_OLD);
        psio_->tocprint(PSIF_DCFT_DPD);
        // L_ref = lambda_IJAB
        global_dpd_->buf4_init(&L_ref, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                               ID("[O>O]-"), ID("[V>V]-"), 0, "Lambda <OO|VV>");
        // L_test = lambda_IjAb
        global_dpd_->buf4_init(&L_test, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                               ID("[O,o]"), ID("[V,v]"), 0, "Lambda <Oo|Vv>");
        global_dpd_->buf4_init(&temp, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                               ID("[O,o]"), ID("[V,v]"), 0, "Lambda <Oo|Vv>");
//        global_dpd_->buf4_print(&temp, "outfile", 1);
        global_dpd_->buf4_sort(&temp, PSIF_DCFT_DPD, qprs, ID("[O,o]"), ID("[V,v]"), "Lambda temp <Oo|Vv>");
        global_dpd_->buf4_close(&temp);
        global_dpd_->buf4_init(&temp, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                               ID("[O,o]"), ID("[V,v]"), 0, "Lambda temp <Oo|Vv>");
//        global_dpd_->buf4_print(&temp, "outfile", 1);
        // L_test -= lambda_JiAb
        dpd_buf4_add(&L_test, &temp, -1.0);

        // 0 ?= L_ref - L_test
        dpd_buf4_add(&L_ref, &L_test, -1.0);
        global_dpd_->buf4_print(&L_ref, "outfile", 1);

        global_dpd_->buf4_close(&temp);
        global_dpd_->buf4_close(&L_ref);
        global_dpd_->buf4_close(&L_test);


//        psio_->close(PSIF_DCFT_DPD, 1);
    }

    // Check: PASSED
    if(0){
        /*
         * Test: if G_IJAB = G_IjAb - G_JiAb
         */
        dpdbuf4 G_ref, G_test, temp;
        psio_->tocprint(PSIF_DCFT_DPD);
        // G_ref = G_IJAB
        global_dpd_->buf4_init(&G_ref, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                               ID("[O>O]-"), ID("[V>V]-"), 0, "G <OO|VV>");
        // G_test = G_IjAb
        global_dpd_->buf4_init(&G_test, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                               ID("[O,o]"), ID("[V,v]"), 0, "G <Oo|Vv>");
        global_dpd_->buf4_init(&temp, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                               ID("[O,o]"), ID("[V,v]"), 0, "G <Oo|Vv>");
        global_dpd_->buf4_sort(&temp, PSIF_DCFT_DPD, qprs, ID("[O,o]"), ID("[V,v]"), "G temp <Oo|Vv>");
        global_dpd_->buf4_close(&temp);
        global_dpd_->buf4_init(&temp, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                               ID("[O,o]"), ID("[V,v]"), 0, "G temp <Oo|Vv>");
        // G_test -= G_JiAb
        dpd_buf4_add(&G_test, &temp, -1.0);

        // 0 ?= G_ref - G_test
        dpd_buf4_add(&G_ref, &G_test, -1.0);
        global_dpd_->buf4_print(&G_ref, "outfile", 1);

        global_dpd_->buf4_close(&temp);
        global_dpd_->buf4_close(&G_ref);
        global_dpd_->buf4_close(&G_test);

    }

    outfile->Printf( "\t*=================================================================================*\n");

}

int
DCFTSolver::run_twostep_dcft_cumulant_updates() {

    // If the system is closed-shell, then...
    if(closedshell_){
        // Set up DIIS
        dpdbuf4 Lab;
        global_dpd_->buf4_init(&Lab, PSIF_LIBTRANS_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                      ID("[O,o]"), ID("[V,v]"), 0, "Lambda <Oo|Vv>"); //Lab are all assigned zero (but with correct dimensions)
        DIISManager lambdaDiisManager(maxdiis_, "DCFT DIIS Lambdas",DIISManager::LargestError,DIISManager::InCore);
        if ((nalpha_ + nbeta_) > 1) {
            lambdaDiisManager.set_error_vector_size(1, DIISEntry::DPDBuf4, &Lab);
            lambdaDiisManager.set_vector_size(1, DIISEntry::DPDBuf4, &Lab);
        }
        global_dpd_->buf4_close(&Lab);

        cumulantDone_ = false;
        int nLambdaIterations = 0;
        // Start density cumulant (lambda) iterations
        while((!cumulantDone_ || !energyConverged_) && nLambdaIterations++ < maxiter_){
            std::string diisString;
            // Build new Tau from current Lambda
            if (options_.get_str("DCFT_FUNCTIONAL") != "CEPA0") {
                // If not CEPA0
                if (options_.get_bool("RELAX_TAU")) { // The default RELAX_TAU is true
                    build_tau();
                    // Compute tau exactly if requested
                    if (exact_tau_) { // The default exact_tau is false
                        refine_tau();
                    }
                    if (options_.get_str("AO_BASIS") == "DISK") { // The default AO_BASIS is DISK !
                        // Transform new Tau to the SO basis
                        transform_tau();

                        // Build Density-Fitted <VV||VV>, <vv||vv>, and <Vv|Vv> terms in the G intermediate
                        if (nLambdaIterations == -1 && macro_cycle == 0){
    //                    if(1){
                            build_DF_tensors();
                        }

                        // Build SO basis tensors for the <VV||VV>, <vv||vv>, and <Vv|Vv> terms in the G intermediate
                        build_AO_tensors();
                        // Check:
//                        outfile->Printf( "\n\tDoing well...\n\n");

                    }
                    else {
                        // Compute GTau contribution for the Fock operator
                        build_gtau();
                    }
                    // Update Fock operator for the F intermediate
                    update_fock();
                }
                else {
                    if (options_.get_str("AO_BASIS") == "DISK") {
                        // Build SO basis tensors for the <VV||VV>, <vv||vv>, and <Vv|Vv> terms in the G intermediate
                        build_AO_tensors();
                    }
                }
            }
            // Build G and F intermediates needed for the density cumulant residual equations and DCFT energy computation
            build_cumulant_intermediates();
            // Compute the residuals for density cumulant equations
            cumulant_convergence_ = compute_cumulant_residual();
            // Check:
//            outfile->Printf( "\n\tcumulant_convergence_ = %f\n\n", cumulant_convergence_);

            // Update density cumulant tensor
            update_cumulant_jacobi();
            if(cumulant_convergence_ < diis_start_thresh_ && (nalpha_ + nbeta_) > 1){
                //Store the DIIS vectors
                dpdbuf4 Lab, Rab;
                global_dpd_->buf4_init(&Rab, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                              ID("[O,o]"), ID("[V,v]"), 0, "R <Oo|Vv>");
                global_dpd_->buf4_init(&Lab, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                              ID("[O,o]"), ID("[V,v]"), 0, "Lambda <Oo|Vv>");

                if(lambdaDiisManager.add_entry(2, &Rab, &Lab)){
                    diisString += "S";
                }
                if(lambdaDiisManager.subspace_size() >= mindiisvecs_ && maxdiis_ > 0){
                    diisString += "/E";
                    lambdaDiisManager.extrapolate(1, &Lab);
                }
                global_dpd_->buf4_close(&Rab);
                global_dpd_->buf4_close(&Lab);
            }
            // Save old DCFT energy
            old_total_energy_ = new_total_energy_;
            // Compute new DCFT energy (lambda contribution)
            if (options_.get_str("DCFT_FUNCTIONAL") == "CEPA0") {
                compute_cepa0_energy();
            } else {
                compute_dcft_energy();
            }
            new_total_energy_ = scf_energy_ + lambda_energy_;
            // Check convergence for density cumulant iterations
            cumulantDone_ = cumulant_convergence_ < cumulant_threshold_;
            energyConverged_ = fabs(new_total_energy_ - old_total_energy_) < cumulant_threshold_;
            if (options_.get_str("ALGORITHM") == "TWOSTEP") {
                outfile->Printf( "\t* %-3d   %12.3e      %12.3e   %12.3e  %21.15f  %-3s *\n",
                        nLambdaIterations, orbitals_convergence_, cumulant_convergence_, new_total_energy_ - old_total_energy_,
                        new_total_energy_, diisString.c_str());
            }
            if (fabs(cumulant_convergence_) > 100.0) throw PSIEXCEPTION("DCFT density cumulant equations diverged");

        }

        return nLambdaIterations;

    }

    // If the system is open-shell, then...
    else{
        // Set up DIIS
        dpdbuf4 Laa, Lab, Lbb;
        global_dpd_->buf4_init(&Laa, PSIF_LIBTRANS_DPD, 0, ID("[O>O]-"), ID("[V>V]-"),
                      ID("[O>O]-"), ID("[V>V]-"), 0, "Lambda <OO|VV>");
        global_dpd_->buf4_init(&Lab, PSIF_LIBTRANS_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                      ID("[O,o]"), ID("[V,v]"), 0, "Lambda <Oo|Vv>");
        global_dpd_->buf4_init(&Lbb, PSIF_LIBTRANS_DPD, 0, ID("[o>o]-"), ID("[v>v]-"),
                      ID("[o>o]-"), ID("[v>v]-"), 0, "Lambda <oo|vv>"); //Laa, Lab, and Lbb are all assigned zero (but with correct dimensions)

        DIISManager lambdaDiisManager(maxdiis_, "DCFT DIIS Lambdas",DIISManager::LargestError,DIISManager::InCore);
        if ((nalpha_ + nbeta_) > 1) {
            lambdaDiisManager.set_error_vector_size(3, DIISEntry::DPDBuf4, &Laa,
                                                    DIISEntry::DPDBuf4, &Lab,
                                                    DIISEntry::DPDBuf4, &Lbb);
            lambdaDiisManager.set_vector_size(3, DIISEntry::DPDBuf4, &Laa,
                                              DIISEntry::DPDBuf4, &Lab,
                                              DIISEntry::DPDBuf4, &Lbb);
        }
        global_dpd_->buf4_close(&Laa);
        global_dpd_->buf4_close(&Lab);
        global_dpd_->buf4_close(&Lbb);
        cumulantDone_ = false;
        int nLambdaIterations = 0;
        // Start density cumulant (lambda) iterations
        while((!cumulantDone_ || !energyConverged_) && nLambdaIterations++ < maxiter_){
            std::string diisString;
            // Build new Tau from current Lambda
            if (options_.get_str("DCFT_FUNCTIONAL") != "CEPA0") {
                // If not CEPA0
                if (options_.get_bool("RELAX_TAU")) { // The default RELAX_TAU is true
                    build_tau();
                    // Compute tau exactly if requested
                    if (exact_tau_) {
                        refine_tau();
                    }
                    if (options_.get_str("AO_BASIS") == "DISK") { // The default AO_BASIS is DISK !
                        // Transform new Tau to the SO basis
                        transform_tau();

                        // Build Density-Fitted <VV||VV>, <vv||vv>, and <Vv|Vv> terms in the G intermediate
                        if (nLambdaIterations == -1 && macro_cycle == 0){
    //                    if(1){
                            build_DF_tensors();
                        }

                        // Build SO basis tensors for the <VV||VV>, <vv||vv>, and <Vv|Vv> terms in the G intermediate
                        build_AO_tensors();
                    }
                    else {
                        // Compute GTau contribution for the Fock operator
                        build_gtau();
                    }
                    // Update Fock operator for the F intermediate
                    update_fock();
                }
                else {
                    if (options_.get_str("AO_BASIS") == "DISK") {
                        // Build SO basis tensors for the <VV||VV>, <vv||vv>, and <Vv|Vv> terms in the G intermediate
                        build_AO_tensors();
                    }
                }
            }
            // Build G and F intermediates needed for the density cumulant residual equations and DCFT energy computation
            build_cumulant_intermediates();
            // Compute the residuals for density cumulant equations
            cumulant_convergence_ = compute_cumulant_residual();
            // Update density cumulant tensor
            update_cumulant_jacobi();
            if(cumulant_convergence_ < diis_start_thresh_ && (nalpha_ + nbeta_) > 1){
                //Store the DIIS vectors
                dpdbuf4 Laa, Lab, Lbb, Raa, Rab, Rbb;
                global_dpd_->buf4_init(&Raa, PSIF_DCFT_DPD, 0, ID("[O>O]-"), ID("[V>V]-"),
                              ID("[O>O]-"), ID("[V>V]-"), 0, "R <OO|VV>");
                global_dpd_->buf4_init(&Rab, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                              ID("[O,o]"), ID("[V,v]"), 0, "R <Oo|Vv>");
                global_dpd_->buf4_init(&Rbb, PSIF_DCFT_DPD, 0, ID("[o>o]-"), ID("[v>v]-"),
                              ID("[o>o]-"), ID("[v>v]-"), 0, "R <oo|vv>");
                global_dpd_->buf4_init(&Laa, PSIF_DCFT_DPD, 0, ID("[O>O]-"), ID("[V>V]-"),
                              ID("[O>O]-"), ID("[V>V]-"), 0, "Lambda <OO|VV>");
                global_dpd_->buf4_init(&Lab, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                              ID("[O,o]"), ID("[V,v]"), 0, "Lambda <Oo|Vv>");
                global_dpd_->buf4_init(&Lbb, PSIF_DCFT_DPD, 0, ID("[o>o]-"), ID("[v>v]-"),
                              ID("[o>o]-"), ID("[v>v]-"), 0, "Lambda <oo|vv>");

                if(lambdaDiisManager.add_entry(6, &Raa, &Rab, &Rbb, &Laa, &Lab, &Lbb)){
                    diisString += "S";
                }
                if(lambdaDiisManager.subspace_size() >= mindiisvecs_ && maxdiis_ > 0){
                    diisString += "/E";
                    lambdaDiisManager.extrapolate(3, &Laa, &Lab, &Lbb);
                }
                global_dpd_->buf4_close(&Raa);
                global_dpd_->buf4_close(&Rab);
                global_dpd_->buf4_close(&Rbb);
                global_dpd_->buf4_close(&Laa);
                global_dpd_->buf4_close(&Lab);
                global_dpd_->buf4_close(&Lbb);
            }
            // Save old DCFT energy
            old_total_energy_ = new_total_energy_;
            // Compute new DCFT energy (lambda contribution)
            if (options_.get_str("DCFT_FUNCTIONAL") == "CEPA0") {
                compute_cepa0_energy();
            } else {
                compute_dcft_energy();
            }
            new_total_energy_ = scf_energy_ + lambda_energy_;
            // Check convergence for density cumulant iterations
            cumulantDone_ = cumulant_convergence_ < cumulant_threshold_;
            energyConverged_ = fabs(new_total_energy_ - old_total_energy_) < cumulant_threshold_;
            if (options_.get_str("ALGORITHM") == "TWOSTEP") {
                outfile->Printf( "\t* %-3d   %12.3e      %12.3e   %12.3e  %21.15f  %-3s *\n",
                        nLambdaIterations, orbitals_convergence_, cumulant_convergence_, new_total_energy_ - old_total_energy_,
                        new_total_energy_, diisString.c_str());
            }
            if (fabs(cumulant_convergence_) > 100.0) throw PSIEXCEPTION("DCFT density cumulant equations diverged");

        }

        return nLambdaIterations;

    }
}


void
DCFTSolver::run_twostep_dcft_orbital_updates() {

    SharedMatrix tmp = SharedMatrix(new Matrix("temp", nirrep_, nsopi_, nsopi_));

    // Set up DIIS
    DIISManager scfDiisManager(maxdiis_, "DCFT DIIS Orbitals",DIISManager::LargestError,DIISManager::InCore);
    if ((nalpha_ + nbeta_) > 1) {
        scfDiisManager.set_error_vector_size(2, DIISEntry::Matrix, scf_error_a_.get(),
                                             DIISEntry::Matrix, scf_error_b_.get());
        scfDiisManager.set_vector_size(2, DIISEntry::Matrix, Fa_.get(),
                                       DIISEntry::Matrix, Fb_.get());
    }
    // Update the orbitals
    int nSCFCycles = 0;
    // Reset the booleans that control the convergence
    densityConverged_ = false;
    energyConverged_ = false;
    outfile->Printf( "\tOrbital Updates\n");
    while((!densityConverged_ || !orbitalsDone_ || !energyConverged_) && (nSCFCycles++ < maxiter_)){
        std::string diisString;
        // Copy core hamiltonian into the Fock matrix array: F = H
        Fa_->copy(so_h_);
        Fb_->copy(so_h_);
        // Build the new Fock matrix from the SO integrals: F += Gbar * Kappa
        process_so_ints();
        // Save F0 = H + G * Kappa for the Fock intermediate update in lambda iterations
        moF0a_->copy(Fa_);
        moF0b_->copy(Fb_);
        moF0a_->transform(Ca_);
        moF0b_->transform(Cb_);
        // Save old SCF energy
        old_total_energy_ = new_total_energy_;
        // Add non-idempotent density contribution (Tau) to the Fock matrix: F += Gbar * Tau
        Fa_->add(g_tau_a_);
        Fb_->add(g_tau_b_);
        // Back up the SO basis Fock before it is symmetrically orthogonalized to transform it to the MO basis
        moFa_->copy(Fa_);
        moFb_->copy(Fb_);
        // Compute new SCF energy
        compute_scf_energy();
        // Check SCF convergence
        orbitals_convergence_ = compute_scf_error_vector();
        orbitalsDone_ = orbitals_convergence_ < orbitals_threshold_;
        if(orbitals_convergence_ < diis_start_thresh_ && (nalpha_ + nbeta_) > 1){
            if(scfDiisManager.add_entry(4, scf_error_a_.get(), scf_error_b_.get(), Fa_.get(), Fb_.get()))
                diisString += "S";
            if(scfDiisManager.subspace_size() > mindiisvecs_ && (nalpha_ + nbeta_) > 1){
                diisString += "/E";
                scfDiisManager.extrapolate(2, Fa_.get(), Fb_.get());
            }
        }
        // Transform the Fock matrix to the symmetrically orhogonalized basis set and digonalize it
        // Obtain new orbitals
        Fa_->transform(s_half_inv_);
        Fa_->diagonalize(tmp, epsilon_a_);
        old_ca_->copy(Ca_);
        Ca_->gemm(false, false, 1.0, s_half_inv_, tmp, 0.0);
        Fb_->transform(s_half_inv_);
        Fb_->diagonalize(tmp, epsilon_b_);
        old_cb_->copy(Cb_);
        Cb_->gemm(false, false, 1.0, s_half_inv_, tmp, 0.0);
        // Make sure that the orbital phase is retained
        correct_mo_phases(false);
        // Update SCF density (Kappa) and check its RMS
        densityConverged_ = update_scf_density() < orbitals_threshold_;
        // Compute the DCFT energy
        new_total_energy_ = scf_energy_ + lambda_energy_;
        // Check convergence of the total DCFT energy
        energyConverged_ = fabs(new_total_energy_ - old_total_energy_) < cumulant_threshold_;
        outfile->Printf( "\t* %-3d   %12.3e      %12.3e   %12.3e  %21.15f  %-3s *\n",
                nSCFCycles, orbitals_convergence_, cumulant_convergence_, new_total_energy_ - old_total_energy_,
                new_total_energy_, diisString.c_str());
        if (fabs(orbitals_convergence_) > 100.0) throw PSIEXCEPTION("DCFT orbital updates diverged");

    }
    // Write orbitals to the checkpoint file
    write_orbitals_to_checkpoint();
    orbitalsDone_ = nSCFCycles == 1;
    energyConverged_ = false;
    // Transform the Fock matrix to the MO basis
    moFa_->transform(Ca_);
    moFb_->transform(Cb_);
    // Transform two-electron integrals to the MO basis using new orbitals, build denominators
    transform_integrals();
}

void
DCFTSolver::run_simult_dcft()
{
    // If the system is closed-shell, then ...
    if(closedshell_){
        // This is the simultaneous orbital/lambda update algorithm
        int cycle = 0;
        outfile->Printf( "\n\n\t*=================================================================================*\n"
                             "\t* Cycle  RMS [F, Kappa]   RMS Lambda Error   delta E        Total Energy     DIIS *\n"
                             "\t*---------------------------------------------------------------------------------*\n");

        SharedMatrix tmp = SharedMatrix(new Matrix("temp", nirrep_, nsopi_, nsopi_));
        // Set up the DIIS manager
        DIISManager diisManager(maxdiis_, "DCFT DIIS vectors");
        dpdbuf4 Lab;
        global_dpd_->buf4_init(&Lab, PSIF_LIBTRANS_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                      ID("[O,o]"), ID("[V,v]"), 0, "Lambda <Oo|Vv>");
        diisManager.set_error_vector_size(3, DIISEntry::Matrix, scf_error_a_.get(),
                                             DIISEntry::Matrix, scf_error_b_.get(),
                                             DIISEntry::DPDBuf4, &Lab);
        diisManager.set_vector_size(3, DIISEntry::Matrix, Fa_.get(),
                                       DIISEntry::Matrix, Fb_.get(),
                                       DIISEntry::DPDBuf4, &Lab);
        global_dpd_->buf4_close(&Lab);
        while((!orbitalsDone_ || !cumulantDone_ || !densityConverged_ || !energyConverged_)
                && cycle++ < maxiter_){
            std::string diisString;
            // Save the old energy
            old_total_energy_ = new_total_energy_;
            // Build new Tau from the density cumulant in the MO basis and transform it the SO basis
            build_tau();
            if (exact_tau_) {
                refine_tau();
            }
            transform_tau();
            // Copy core hamiltonian into the Fock matrix array: F = H
            Fa_->copy(so_h_);
            Fb_->copy(so_h_);
            // Build the new Fock matrix from the SO integrals: F += Gbar * Kappa
            process_so_ints();
            // Add non-idempotent density contribution (Tau) to the Fock matrix: F += Gbar * Tau
            Fa_->add(g_tau_a_);
            Fb_->add(g_tau_b_);
            // Back up the SO basis Fock before it is symmetrically orthogonalized to transform it to the MO basis
            moFa_->copy(Fa_);
            moFb_->copy(Fb_);
            // Transform the Fock matrix to the MO basis
            moFa_->transform(Ca_);
            moFb_->transform(Cb_);
            // Compute new SCF energy
            compute_scf_energy();
            // Add SCF energy contribution to the total DCFT energy
            new_total_energy_ = scf_energy_;
            // Check SCF convergence
            orbitals_convergence_ = compute_scf_error_vector();
            orbitalsDone_ = orbitals_convergence_ < orbitals_threshold_;
            // Build G and F intermediates needed for the density cumulant residual equations and DCFT energy computation
            build_cumulant_intermediates();
            // Compute the residuals for density cumulant equations
            cumulant_convergence_ = compute_cumulant_residual();
            if (fabs(cumulant_convergence_) > 100.0) throw PSIEXCEPTION("DCFT density cumulant equations diverged");
            // Check convergence for density cumulant iterations
            cumulantDone_ = cumulant_convergence_ < cumulant_threshold_;
            // Update density cumulant tensor
            update_cumulant_jacobi();
            // Compute new DCFT energy (lambda contribution)
            compute_dcft_energy();
            // Add lambda energy to the DCFT total energy
            new_total_energy_ += lambda_energy_;
            // Check convergence of the total DCFT energy
            energyConverged_ = fabs(old_total_energy_ - new_total_energy_) < cumulant_threshold_;
            if(orbitals_convergence_ < diis_start_thresh_ && cumulant_convergence_ < diis_start_thresh_){
                //Store the DIIS vectors
                dpdbuf4 Lab, Rab;
                global_dpd_->buf4_init(&Rab, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                              ID("[O,o]"), ID("[V,v]"), 0, "R <Oo|Vv>");
                global_dpd_->buf4_init(&Lab, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                              ID("[O,o]"), ID("[V,v]"), 0, "Lambda <Oo|Vv>");
                if(diisManager.add_entry(6, scf_error_a_.get(), scf_error_b_.get(), &Rab,
                                           Fa_.get(), Fb_.get(), &Lab)){
                    diisString += "S";
                }
                if(diisManager.subspace_size() > mindiisvecs_){
                    diisString += "/E";
                    diisManager.extrapolate(3, Fa_.get(), Fb_.get(), &Lab);
                }
                global_dpd_->buf4_close(&Rab);
                global_dpd_->buf4_close(&Lab);
            }

            // Transform the Fock matrix to the symmetrically orhogonalized basis set and digonalize it
            // Obtain new orbitals
            Fa_->transform(s_half_inv_);
            Fa_->diagonalize(tmp, epsilon_a_);
            old_ca_->copy(Ca_);
            Ca_->gemm(false, false, 1.0, s_half_inv_, tmp, 0.0);
            Fb_->transform(s_half_inv_);
            Fb_->diagonalize(tmp, epsilon_b_);
            old_cb_->copy(Cb_);
            Cb_->gemm(false, false, 1.0, s_half_inv_, tmp, 0.0);
            // Make sure that the orbital phase is retained
            if(!correct_mo_phases(false)){
                outfile->Printf("\t\tThere was a problem correcting the MO phases.\n"
                                "\t\tIf this does not converge, try ALGORITHM=TWOSTEP\n");
            }
            // Write orbitals to the checkpoint file
            write_orbitals_to_checkpoint();
            // Transform two-electron integrals to the MO basis using new orbitals, build denominators
            transform_integrals();
            // Update SCF density (Kappa) and check its RMS
            densityConverged_ = update_scf_density() < orbitals_threshold_;
            // If we've performed enough lambda updates since the last orbitals
            // update, reset the counter so another SCF update is performed
            outfile->Printf( "\t* %-3d   %12.3e      %12.3e   %12.3e  %21.15f  %-3s *\n",
                    cycle, orbitals_convergence_, cumulant_convergence_, new_total_energy_ - old_total_energy_,
                    new_total_energy_, diisString.c_str());

        }

        outfile->Printf( "\t*=================================================================================*\n");

    }
    // If the system is open-shell, then ...
    else{
        // This is the simultaneous orbital/lambda update algorithm
        int cycle = 0;
        outfile->Printf( "\n\n\t*=================================================================================*\n"
                             "\t* Cycle  RMS [F, Kappa]   RMS Lambda Error   delta E        Total Energy     DIIS *\n"
                             "\t*---------------------------------------------------------------------------------*\n");

        SharedMatrix tmp = SharedMatrix(new Matrix("temp", nirrep_, nsopi_, nsopi_));
        // Set up the DIIS manager
        DIISManager diisManager(maxdiis_, "DCFT DIIS vectors");
        dpdbuf4 Laa, Lab, Lbb;
        global_dpd_->buf4_init(&Laa, PSIF_LIBTRANS_DPD, 0, ID("[O>O]-"), ID("[V>V]-"),
                      ID("[O>O]-"), ID("[V>V]-"), 0, "Lambda <OO|VV>");
        global_dpd_->buf4_init(&Lab, PSIF_LIBTRANS_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                      ID("[O,o]"), ID("[V,v]"), 0, "Lambda <Oo|Vv>");
        global_dpd_->buf4_init(&Lbb, PSIF_LIBTRANS_DPD, 0, ID("[o>o]-"), ID("[v>v]-"),
                      ID("[o>o]-"), ID("[v>v]-"), 0, "Lambda <oo|vv>");
        diisManager.set_error_vector_size(5, DIISEntry::Matrix, scf_error_a_.get(),
                                             DIISEntry::Matrix, scf_error_b_.get(),
                                             DIISEntry::DPDBuf4, &Laa,
                                             DIISEntry::DPDBuf4, &Lab,
                                             DIISEntry::DPDBuf4, &Lbb);
        diisManager.set_vector_size(5, DIISEntry::Matrix, Fa_.get(),
                                       DIISEntry::Matrix, Fb_.get(),
                                       DIISEntry::DPDBuf4, &Laa,
                                       DIISEntry::DPDBuf4, &Lab,
                                       DIISEntry::DPDBuf4, &Lbb);
        global_dpd_->buf4_close(&Laa);
        global_dpd_->buf4_close(&Lab);
        global_dpd_->buf4_close(&Lbb);
        while((!orbitalsDone_ || !cumulantDone_ || !densityConverged_ || !energyConverged_)
                && cycle++ < maxiter_){
            std::string diisString;
            // Save the old energy
            old_total_energy_ = new_total_energy_;
            // Build new Tau from the density cumulant in the MO basis and transform it the SO basis
            build_tau();
            if (exact_tau_) {
                refine_tau();
            }
            transform_tau();
            // Copy core hamiltonian into the Fock matrix array: F = H
            Fa_->copy(so_h_);
            Fb_->copy(so_h_);
            // Build the new Fock matrix from the SO integrals: F += Gbar * Kappa
            process_so_ints();
            // Add non-idempotent density contribution (Tau) to the Fock matrix: F += Gbar * Tau
            Fa_->add(g_tau_a_);
            Fb_->add(g_tau_b_);
            // Back up the SO basis Fock before it is symmetrically orthogonalized to transform it to the MO basis
            moFa_->copy(Fa_);
            moFb_->copy(Fb_);
            // Transform the Fock matrix to the MO basis
            moFa_->transform(Ca_);
            moFb_->transform(Cb_);
            // Compute new SCF energy
            compute_scf_energy();
            // Add SCF energy contribution to the total DCFT energy
            new_total_energy_ = scf_energy_;
            // Check SCF convergence
            orbitals_convergence_ = compute_scf_error_vector();
            orbitalsDone_ = orbitals_convergence_ < orbitals_threshold_;
            // Build G and F intermediates needed for the density cumulant residual equations and DCFT energy computation
            build_cumulant_intermediates();
            // Compute the residuals for density cumulant equations
            cumulant_convergence_ = compute_cumulant_residual();
            if (fabs(cumulant_convergence_) > 100.0) throw PSIEXCEPTION("DCFT density cumulant equations diverged");
            // Check convergence for density cumulant iterations
            cumulantDone_ = cumulant_convergence_ < cumulant_threshold_;
            // Update density cumulant tensor
            update_cumulant_jacobi();
            // Compute new DCFT energy (lambda contribution)
            compute_dcft_energy();
            // Add lambda energy to the DCFT total energy
            new_total_energy_ += lambda_energy_;
            // Check convergence of the total DCFT energy
            energyConverged_ = fabs(old_total_energy_ - new_total_energy_) < cumulant_threshold_;
            if(orbitals_convergence_ < diis_start_thresh_ && cumulant_convergence_ < diis_start_thresh_){
                //Store the DIIS vectors
                dpdbuf4 Laa, Lab, Lbb, Raa, Rab, Rbb;
                global_dpd_->buf4_init(&Raa, PSIF_DCFT_DPD, 0, ID("[O>O]-"), ID("[V>V]-"),
                              ID("[O>O]-"), ID("[V>V]-"), 0, "R <OO|VV>");
                global_dpd_->buf4_init(&Rab, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                              ID("[O,o]"), ID("[V,v]"), 0, "R <Oo|Vv>");
                global_dpd_->buf4_init(&Rbb, PSIF_DCFT_DPD, 0, ID("[o>o]-"), ID("[v>v]-"),
                              ID("[o>o]-"), ID("[v>v]-"), 0, "R <oo|vv>");
                global_dpd_->buf4_init(&Laa, PSIF_DCFT_DPD, 0, ID("[O>O]-"), ID("[V>V]-"),
                              ID("[O>O]-"), ID("[V>V]-"), 0, "Lambda <OO|VV>");
                global_dpd_->buf4_init(&Lab, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                              ID("[O,o]"), ID("[V,v]"), 0, "Lambda <Oo|Vv>");
                global_dpd_->buf4_init(&Lbb, PSIF_DCFT_DPD, 0, ID("[o>o]-"), ID("[v>v]-"),
                              ID("[o>o]-"), ID("[v>v]-"), 0, "Lambda <oo|vv>");
                if(diisManager.add_entry(10, scf_error_a_.get(), scf_error_b_.get(), &Raa, &Rab, &Rbb,
                                           Fa_.get(), Fb_.get(), &Laa, &Lab, &Lbb)){
                    diisString += "S";
                }
                if(diisManager.subspace_size() > mindiisvecs_){
                    diisString += "/E";
                    diisManager.extrapolate(5, Fa_.get(), Fb_.get(), &Laa, &Lab, &Lbb);
                }
                global_dpd_->buf4_close(&Raa);
                global_dpd_->buf4_close(&Rab);
                global_dpd_->buf4_close(&Rbb);
                global_dpd_->buf4_close(&Laa);
                global_dpd_->buf4_close(&Lab);
                global_dpd_->buf4_close(&Lbb);
            }

            // Transform the Fock matrix to the symmetrically orhogonalized basis set and digonalize it
            // Obtain new orbitals
            Fa_->transform(s_half_inv_);
            Fa_->diagonalize(tmp, epsilon_a_);
            old_ca_->copy(Ca_);
            Ca_->gemm(false, false, 1.0, s_half_inv_, tmp, 0.0);
            Fb_->transform(s_half_inv_);
            Fb_->diagonalize(tmp, epsilon_b_);
            old_cb_->copy(Cb_);
            Cb_->gemm(false, false, 1.0, s_half_inv_, tmp, 0.0);
            // Make sure that the orbital phase is retained
            if(!correct_mo_phases(false)){
                outfile->Printf("\t\tThere was a problem correcting the MO phases.\n"
                                "\t\tIf this does not converge, try ALGORITHM=TWOSTEP\n");
            }
            // Write orbitals to the checkpoint file
            write_orbitals_to_checkpoint();
            // Transform two-electron integrals to the MO basis using new orbitals, build denominators
            transform_integrals();
            // Update SCF density (Kappa) and check its RMS
            densityConverged_ = update_scf_density() < orbitals_threshold_;
            // If we've performed enough lambda updates since the last orbitals
            // update, reset the counter so another SCF update is performed
            outfile->Printf( "\t* %-3d   %12.3e      %12.3e   %12.3e  %21.15f  %-3s *\n",
                    cycle, orbitals_convergence_, cumulant_convergence_, new_total_energy_ - old_total_energy_,
                    new_total_energy_, diisString.c_str());

        }

        outfile->Printf( "\t*=================================================================================*\n");

    }
}

}} // Namespaces

