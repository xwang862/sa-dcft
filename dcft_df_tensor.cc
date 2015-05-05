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

#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <sstream>

#include <psifiles.h>
#include <libciomr/libciomr.h>
#include <libpsio/psio.h>
#include <libchkpt/chkpt.h>
#include <libpsio/psio.hpp>
#include <libchkpt/chkpt.hpp>
#include <libiwl/iwl.h>
#include <libqt/qt.h>
#include <libmints/mints.h>
#include <psi4-dec.h>

#include "dcft.h"
#include "defines.h"
#include <vector>
#include <liboptions/liboptions.h>

#include <libtrans/integraltransform.h>
#include <libtrans/mospace.h>
#include <libdpd/dpd.h>
#include <libdiis/diismanager.h>

#include <lib3index/3index.h>
#include <libmints/mints.h>
#include <libmints/sieve.h>
#include <libfock/jk.h>
#include <libfock/apps.h>
#include <physconst.h>

using namespace boost;

namespace psi{

namespace dfdcft{

/**
 * Compute the density-fitted <VV||VV>, <vv||vv>, and <Vv|Vv> tensors in G intermediates
 * and contract with lambda_ijcd
 */

void DCFTSolver::build_DF_tensors()
{
    dcft_timer_on("DCFTSolver::build_df_tensors()");
    outfile->Printf( "\n\tBuilding density-fitted tensors <VV||VV> ... \n\n");

    Caocc_a_ = Ca_subset("AO","ACTIVE_OCC");
    Cavir_a_ = Ca_subset("AO","ACTIVE_VIR");
    Caocc_b_ = Cb_subset("AO","ACTIVE_OCC");
    Cavir_b_ = Cb_subset("AO","ACTIVE_VIR");

    eps_aocc_a_ = epsilon_a_subset("AO","ACTIVE_OCC");
    eps_avir_a_ = epsilon_a_subset("AO","ACTIVE_VIR");
    eps_aocc_b_ = epsilon_b_subset("AO","ACTIVE_OCC");
    eps_avir_b_ = epsilon_b_subset("AO","ACTIVE_VIR");

    // (ab|cd) = (ab|Q)(Q|cd)
    // Form (A|ac) = (A|mn) C_ma C_nc
    dcft_timer_on("DCFTSolver::df_form_Aac");
    df_tensor_form_Aac();
    dcft_timer_off("DCFTSolver::df_form_Aac");

    // Apply the fitting (Q|ac) = J_QA^-1/2 (A|ac)
    dcft_timer_on("DCFTSolver::df_form_Qac");
    df_tensor_form_Qac();
    dcft_timer_off("DCFTSolver::df_form_Qac");

    // Form df_tau(temp)_ijab = Sum_cd ( gbar_cdab * lambda_ijcd ) . Choose a better name later...
    dcft_timer_on("DCFTSolver::df_form_tau(temp)_OOVV");
    df_tensor_form_tau_temp_OOVV();
    dcft_timer_off("DCFTSolver::df_form_tau(temp)_OOVV");

    dcft_timer_off("DCFTSolver::build_df_tensors()");
}

void DCFTSolver::df_tensor_form_Aac()
{
    // Schwarz Sieve
    boost::shared_ptr<ERISieve> sieve(new ERISieve(basisset_, options_.get_double("INTS_TOLERANCE")));
    double test_tol = options_.get_double("INTS_TOLERANCE");
    const std::vector<std::pair<int,int> >& shell_pairs = sieve->shell_pairs();
    const size_t npairs = shell_pairs.size();

    // ERI objects
    // Parallel version of DF-tensor building can be added later
    int nthread = 1;

    boost::shared_ptr<IntegralFactory> factory(new IntegralFactory(ribasis_, BasisSet::zero_ao_basis_set(), basisset_, basisset_));
    std::vector<boost::shared_ptr<TwoBodyAOInt> > eri;
    std::vector<const double*> buffer;
    for (int thread = 0; thread < nthread; ++thread){
        eri.push_back(boost::shared_ptr<TwoBodyAOInt>(factory->eri()));
        buffer.push_back(eri[thread]->buffer());
    }

    // Sizing
    int nso = basisset_->nbf();
    int naux = ribasis_->nbf();
    int navir_a = Cavir_a_->colspi()[0];
    int navir_b = Cavir_b_->colspi()[0];
    int navir = (navir_a > navir_b ? navir_a : navir_b);
    int maxQ = ribasis_->max_function_per_shell();

    // Max block size in naux
    ULI Amn_cost_per_row = nso * (ULI) nso;
    ULI Ama_cost_per_row = nso * (ULI) navir;
    ULI Aac_cost_per_row = navir * (ULI) navir;
    ULI total_cost_per_row = Amn_cost_per_row + Ama_cost_per_row + Aac_cost_per_row;
    /* A new name for this memory factor needs to be created ... */
    ULI doubles = ((ULI) (options_.get_double("DFMP2_MEM_FACTOR") * memory_ / 8L));
    ULI max_temp = doubles / (total_cost_per_row);
    int max_naux = (max_temp > (ULI) naux ? naux : max_temp);
    max_naux = (max_naux < maxQ ? maxQ : max_naux);

    // Block extens
    std::vector<int> block_Q_starts;
    int counter = 0;
    block_Q_starts.push_back(0);
    for (int Q = 0; Q < ribasis_->nshell(); ++Q){
        int nQ = ribasis_->shell(Q).nfunction();
        if (counter + nQ > max_naux){
            counter = 0;
            block_Q_starts.push_back(Q);
        }
        counter += Q;
    }
    block_Q_starts.push_back(ribasis_->nshell());

    // Tensor blocks
    SharedMatrix Amn(new Matrix("(A|mn) Block", max_naux, nso * (ULI) nso));
    SharedMatrix Ama(new Matrix("(A|ma) Block", max_naux, nso * (ULI) navir));
    SharedMatrix Aac(new Matrix("(A|ac) Block", max_naux, navir * (ULI) navir));
    double** Amnp = Amn->pointer();
    double** Amap = Ama->pointer();
    double** Aacp = Aac->pointer();

    // C Matrices
    double** Cavirap = Cavir_a_->pointer();
    double** Cavirbp = Cavir_b_->pointer();

    psio_->open(PSIF_DFMP2_AIA, PSIO_OPEN_OLD);
    psio_address next_AAC = PSIO_ZERO;
    psio_->open(PSIF_DFMP2_QIA, PSIO_OPEN_OLD);
    psio_address next_QAC = PSIO_ZERO;

    // Loop over blocks of Qshell
    for (int block = 0; block < block_Q_starts.size() - 1; block++){
        // Block sizing/offests
        int Qstart = block_Q_starts[block];
        int Qstop = block_Q_starts[block+1];
        int qoff = ribasis_->shell(Qstart).function_index();
        int nrows = (Qstop == ribasis_->nshell() ?
                         ribasis_->nbf() - ribasis_->shell(Qstart).function_index() :
                         ribasis_->shell(Qstop).function_index() - ribasis_->shell(Qstart).function_index());

        // Clear Amn for Schwarz sieve
        ::memset((void*) Amnp[0], '\0', sizeof(double) * nrows * nso * nso);

        // Compute TEI tensor block (A|mn)
        dcft_timer_on("DCFTSolver::DF tensor (A|mn)");
        for (long int QMN = 0L; QMN < (Qstop - Qstart) * (ULI) npairs; QMN++) {

            int thread = 0;

            int Q =  QMN / npairs + Qstart;
            int MN = QMN % npairs;

            std::pair<int,int> pair = shell_pairs[MN];
            int M = pair.first;
            int N = pair.second;

            int nq = ribasis_->shell(Q).nfunction();
            int nm = basisset_->shell(M).nfunction();
            int nn = basisset_->shell(N).nfunction();

            int sq =  ribasis_->shell(Q).function_index();
            int sm =  basisset_->shell(M).function_index();
            int sn =  basisset_->shell(N).function_index();

            eri[thread]->compute_shell(Q,0,M,N);

            for (int oq = 0; oq < nq; oq++) {
                for (int om = 0; om < nm; om++) {
                    for (int on = 0; on < nn; on++) {
                        Amnp[sq + oq - qoff][(om + sm) * nso + (on + sn)] =
                        Amnp[sq + oq - qoff][(on + sn) * nso + (om + sm)] =
                        buffer[thread][oq * nm * nn + om * nn + on];
                    }
                }
            }
        }
        dcft_timer_off("DCFTSolver::DF tensor (A|mn)");

        // => Alpha Case <= //

        // Compute (A|ma) tensor block (A|mn) C_na
        dcft_timer_on("DCFTSolver::DF tensor (A|ma)");
        C_DGEMM('N', 'N', nrows*(ULI)nso, navir_a, nso, 1.0, Amnp[0], nso, Cavirap[0], navir_a, 0.0, Amap[0], navir);
        dcft_timer_off("DCFTSolver::DF tensor (A|ma)");

        // Compute (A|ac) tensor block (A|ma) C_mc
        dcft_timer_on("DCFTSolver::DF tensor (A|ac)");
        for(int row = 0; row < nrows; ++row){
            C_DGEMM('T', 'N', navir_a, navir_a, nso, 1.0, Amap[row], navir, Cavirap[0], navir_a, 0.0, &Aacp[0][row*(ULI)navir_a*navir_a], navir_a);
        }
        dcft_timer_off("DCFTSolver::DF tensor (A|ac)");

        // Stripe (A|ac) out to disk
        dcft_timer_on("DCFTSolver::DF tensor Aac Write");
        psio_->write(PSIF_DFMP2_AIA, "(A|ac)", (char*)Aacp[0], sizeof(double)*nrows*navir_a*navir_a, next_AAC, &next_AAC);
        dcft_timer_off("DCFTSolver::DF tensor Aac Write");

        // => Beta Case => //

        // Compute (A|ma) tensor block (A|mn) C_na
        dcft_timer_on("DCFTSolver::DF tensor (A|ma)");
        C_DGEMM('N', 'N', nrows*(ULI)nso, navir_b, nso, 1.0, Amnp[0], nso, Cavirbp[0], navir_b, 0.0, Amap[0], navir);
        dcft_timer_off("DCFTSolver::DF tensor (A|ma)");

        // Compute (A|ac) tensor block (A|ma) C_mc
        dcft_timer_on("DCFTSolver::DF tensor (A|ac)");
        for (int row = 0; row < nrows; ++row){
            C_DGEMM('T', 'N', navir_b, navir_b, nso, 1.0, Amap[row], navir, Cavirbp[0], navir_b, 0.0, &Aacp[0][row*(ULI)navir_b*navir_b], navir_b);
        }
        dcft_timer_off("DCFTSolver::DF tensor (A|ac)");

        // Stripe (A|ac) out to disk
        dcft_timer_on("DCFTSolver::DF tensor Aac Write");
        psio_->write(PSIF_DFMP2_QIA, "(A|ac)", (char*)Aacp[0], sizeof(double)*nrows*navir_b*navir_b, next_QAC, &next_QAC);
        dcft_timer_off("DCFTSolver::DF tensor Aac Write");
    }

    psio_->close(PSIF_DFMP2_AIA, 1);
    psio_->close(PSIF_DFMP2_QIA, 1);
}

void DCFTSolver::df_tensor_form_Qac()
{
    SharedMatrix Jm12 = form_inverse_metric();
    apply_fitting(Jm12, PSIF_DFMP2_AIA, "(A|ac)", "(Q|ac)", ribasis_->nbf(), Cavir_a_->colspi()[0] * (ULI) Cavir_a_->colspi()[0]);
    apply_fitting(Jm12, PSIF_DFMP2_QIA, "(A|ac)", "(Q|ac)", ribasis_->nbf(), Cavir_b_->colspi()[0] * (ULI) Cavir_b_->colspi()[0]);
    // Check:
//    psio_->open(PSIF_DFMP2_AIA, PSIO_OPEN_OLD);
//    psio_->open(PSIF_DFMP2_QIA, PSIO_OPEN_OLD);
//    psio_->tocprint(PSIF_DFMP2_AIA);
//    psio_->tocprint(PSIF_DFMP2_QIA);
//    psio_->close(PSIF_DFMP2_AIA, 1);
//    psio_->close(PSIF_DFMP2_QIA, 1);

}

void DCFTSolver::apply_fitting(SharedMatrix Jm12, unsigned int file, const char *key1, const char *key2, ULI naux, ULI npq)
{
    // Memory constraints
    ULI Jmem = naux * naux;
    ULI doubles = (ULI) (options_.get_double("DFMP2_MEM_FACTOR") * (memory_ / 8L));
    if (doubles < 2L * Jmem) {
        throw PSIEXCEPTION("DFMP2: More memory required for tractable disk transpose");
    }
    ULI rem = (doubles - Jmem) / 2L;
    ULI max_npq = (rem / naux);
    max_npq = (max_npq > npq ? npq : max_npq);
    max_npq = (max_npq < 1L ? 1L : max_npq);

    // Block sizing
    std::vector<ULI> pq_starts;
    pq_starts.push_back(0);
    for (ULI pq = 0L; pq < npq; pq+=max_npq) {
        if (pq + max_npq >= npq) {
            pq_starts.push_back(npq);
        } else {
            pq_starts.push_back(pq + max_npq);
        }
    }
    //block_status(pq_starts, __FILE__,__LINE__);

    // Tensor blocks
    SharedMatrix Apq(new Matrix("Apq", naux, max_npq));
    SharedMatrix Qpq(new Matrix("Qpq", max_npq, naux));
    double** Apqp = Apq->pointer();
    double** Qpqp = Qpq->pointer();
    double** Jp   = Jm12->pointer();

    // Loop through blocks
    psio_->open(file, PSIO_OPEN_OLD);
    psio_address next_APQ = PSIO_ZERO;
    psio_address next_QPQ = PSIO_ZERO;
    for (int block = 0; block < pq_starts.size() - 1; block++) {

        // Sizing
        ULI pq_start = pq_starts[block];
        ULI pq_stop  = pq_starts[block+1];
        ULI ncols = pq_stop - pq_start;

        // Read Apq
        dcft_timer_on("DCFTSolver::DFMP2 Apq Read");
        for (ULI Q = 0; Q < naux; Q++) {
            next_APQ = psio_get_address(PSIO_ZERO,sizeof(double)*(Q*npq+pq_start));
            psio_->read(file,key1,(char*)Apqp[Q],sizeof(double)*ncols,next_APQ,&next_APQ);
        }
        dcft_timer_off("DCFTSolver::DFMP2 Apq Read");

        if (debug_){
            psio_->tocprint(file);
        }

        // Apply Fitting
        dcft_timer_on("DCFTSolver::DFMP2 (Q|A)(A|pq)");
        C_DGEMM('T','N',ncols,naux,naux,1.0,Apqp[0],max_npq,Jp[0],naux,0.0,Qpqp[0],naux);
        dcft_timer_off("DCFTSolver::DFMP2 (Q|A)(A|pq)");

        // Write Qpq
        dcft_timer_on("DCFTSolver::DFMP2 Qpq Write");
        psio_->write(file,key2,(char*)Qpqp[0],sizeof(double)*ncols*naux,next_QPQ,&next_QPQ);
        dcft_timer_off("DCFTSolver::DFMP2 Qpq Write");

        if (debug_){
            psio_->tocprint(file);
        }

    }
    psio_->close(file, 1);
}

void DCFTSolver::df_tensor_form_tau_temp_OOVV()
{
    dpdfile4 tau_temp;
    dpdbuf4 lambda;

    /* Make DPD - Non-DPD orbital pairs */
    map_orbitals();

    /*
     * tau(temp)_ijab = 1/2 Sum_cd (gbar_cdab lambda_ijcd)
     * Using density fitting, gbar_cdab = (ca|Q)(Q|db) - (cb|Q)(Q|da)
     *                                  = (ac|Q)(Q|bd) - (ad|Q)(Q|bc)
     */

    /* => tau(temp)_IJAB <= */ {
    // tau(temp)_IJAB = 1/2 Sum_CD ( gbar_CDAB lambda_IJCD )
    global_dpd_->buf4_init(&lambda, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                           ID("[O>O]-"), ID("[V>V]-"), 0, "Lambda <OO|VV>");
    global_dpd_->file4_init(&tau_temp, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"), "DF tau(temp) <OO|VV>");
    // Check:
    //    global_dpd_->file4_print(&tau_temp, "outfile");

    // Sizing
    int navir = Cavir_a_->colspi()[0];
    int naux = ribasis_->nbf();

    // Read Qac in memory. Will treat memory limit later ...
    SharedMatrix Qac (new Matrix("Qac", navir * (ULI) navir, naux));
    SharedMatrix Qbd (new Matrix("Qbd", navir * (ULI) navir, naux));
    double** Qacp = Qac->pointer();
    double** Qbdp = Qbd->pointer();

    // Read Qac all
    psio_->open(PSIF_DFMP2_AIA, PSIO_OPEN_OLD);
    psio_->read_entry(PSIF_DFMP2_AIA, "(Q|ac)", (char*) Qacp[0], sizeof(double) * (ULI) navir * navir * naux);
    psio_->read_entry(PSIF_DFMP2_AIA, "(Q|ac)", (char*) Qbdp[0], sizeof(double) * (ULI) navir * navir * naux);
    psio_->close(PSIF_DFMP2_AIA, 1);

    SharedMatrix Icd = SharedMatrix(new Matrix("Icd", navir, navir));
    for(int h = 0; h < nirrep_; ++h){
        global_dpd_->file4_mat_irrep_init(&tau_temp, h);
        global_dpd_->buf4_mat_irrep_init(&lambda, h);
        global_dpd_->buf4_mat_irrep_rd(&lambda, h);
        for (int IJ = 0; IJ < tau_temp.params->rowtot[h]; ++IJ){
            for (int AB = 0; AB < tau_temp.params->coltot[h]; ++AB){
                int A = tau_temp.params->colorb[h][AB][0];
                int B = tau_temp.params->colorb[h][AB][1];
                int DFA = aVirMap[A];
                int DFB = aVirMap[B];
                double** Icdp = Icd->pointer();
                C_DGEMM('N', 'T', navir, navir, naux, 1.0, Qacp[DFA*navir], naux, Qbdp[DFB*navir], naux, 0.0, Icdp[0], navir);
                double IJAB = 0.0;
                for (int CD = 0; CD < lambda.params->coltot[h]; ++CD){
                    int C = lambda.params->colorb[h][CD][0];
                    int D = lambda.params->colorb[h][CD][1];
                    int DFC = aVirMap[C];
                    int DFD = aVirMap[D];
                    double ACBD = Icdp[DFC][DFD];
                    double ADBC = Icdp[DFD][DFC];
                    IJAB += (ACBD - ADBC) * lambda.matrix[h][IJ][CD];
                }
                IJAB *= 0.5;
                tau_temp.matrix[h][IJ][AB] = IJAB;
            }
        }
        global_dpd_->buf4_mat_irrep_close(&lambda,h);
        global_dpd_->file4_mat_irrep_wrt(&tau_temp, h);
        global_dpd_->file4_mat_irrep_close(&tau_temp, h);
    }
    //Check:
//  global_dpd_->file4_print(&tau_temp, "outfile");
    global_dpd_->file4_close(&tau_temp);
    global_dpd_->buf4_close(&lambda);

    /* End of tau(temp)_IJAB */}

    //Check:
//    dpdbuf4 test;
//    global_dpd_->buf4_init(&test, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
//                           ID("[O,O]"), ID("[V,V]"), 0, "DF tau(temp) <OO|VV>");
//    global_dpd_->buf4_print(&test, "outfile", 1);
//    global_dpd_->buf4_close(&test);

    /* => tau(temp)_ijab <= */ {
    // tau(temp)_ijab = 1/2 Sum_cd ( gbar_cdab lambda_ijcd )
    global_dpd_->buf4_init(&lambda, PSIF_DCFT_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                           ID("[o>o]-"), ID("[v>v]-"), 0, "Lambda <oo|vv>");
    global_dpd_->file4_init(&tau_temp, PSIF_DCFT_DPD, 0, ID("[o,o]"), ID("[v,v]"), "DF tau(temp) <oo|vv>");
    // Check:
    //    global_dpd_->file4_print(&tau_temp, "outfile");

    // Sizing
    int navir = Cavir_b_->colspi()[0];
    int naux = ribasis_->nbf();

    // Read Qac in memory. Will treat memory limit later ...
    SharedMatrix Qac (new Matrix("Qac", navir * (ULI) navir, naux));
    SharedMatrix Qbd (new Matrix("Qbd", navir * (ULI) navir, naux));
    double** Qacp = Qac->pointer();
    double** Qbdp = Qbd->pointer();

    // Read Qac all
    psio_->open(PSIF_DFMP2_QIA, PSIO_OPEN_OLD);
    psio_->read_entry(PSIF_DFMP2_QIA, "(Q|ac)", (char*) Qacp[0], sizeof(double) * (ULI) navir * navir * naux);
    psio_->read_entry(PSIF_DFMP2_QIA, "(Q|ac)", (char*) Qbdp[0], sizeof(double) * (ULI) navir * navir * naux);
    psio_->close(PSIF_DFMP2_QIA, 1);

    SharedMatrix Icd = SharedMatrix(new Matrix("Icd", navir, navir));
    for(int h = 0; h < nirrep_; ++h){
        global_dpd_->file4_mat_irrep_init(&tau_temp, h);
        global_dpd_->buf4_mat_irrep_init(&lambda, h);
        global_dpd_->buf4_mat_irrep_rd(&lambda, h);
        for (int ij = 0; ij < tau_temp.params->rowtot[h]; ++ij){
            for (int ab = 0; ab < tau_temp.params->coltot[h]; ++ab){
                int a = tau_temp.params->colorb[h][ab][0];
                int b = tau_temp.params->colorb[h][ab][1];
                int DFa = bVirMap[a];
                int DFb = bVirMap[b];
                double** Icdp = Icd->pointer();
                C_DGEMM('N', 'T', navir, navir, naux, 1.0, Qacp[DFa*navir], naux, Qbdp[DFb*navir], naux, 0.0, Icdp[0], navir);
                double ijab = 0.0;
                for (int cd = 0; cd < lambda.params->coltot[h]; ++cd){
                    int c = lambda.params->colorb[h][cd][0];
                    int d = lambda.params->colorb[h][cd][1];
                    int DFc = bVirMap[c];
                    int DFd = bVirMap[d];
                    double acbd = Icdp[DFc][DFd];
                    double adbc = Icdp[DFd][DFc];
                    ijab += (acbd - adbc) * lambda.matrix[h][ij][cd];
                }
                ijab *= 0.5;
                tau_temp.matrix[h][ij][ab] = ijab;
            }
        }
        global_dpd_->buf4_mat_irrep_close(&lambda,h);
        global_dpd_->file4_mat_irrep_wrt(&tau_temp, h);
        global_dpd_->file4_mat_irrep_close(&tau_temp, h);
    }
    //Check:
//  global_dpd_->file4_print(&tau_temp, "outfile");
    global_dpd_->file4_close(&tau_temp);
    global_dpd_->buf4_close(&lambda);

    /* End of tau(temp)_ijab */}

    //Check:
//    dpdbuf4 test;
//    global_dpd_->buf4_init(&test, PSIF_DCFT_DPD, 0, ID("[o,o]"), ID("[v,v]"),
//                           ID("[o,o]"), ID("[v,v]"), 0, "DF tau(temp) <oo|vv>");
//    global_dpd_->buf4_print(&test, "outfile", 1);
//    global_dpd_->buf4_close(&test);

    /* => tau(temp)_IjAb <= */ {
    // tau(temp)_IjAb = Sum_Cd ( g_CdAb lambda_IjCd )
    global_dpd_->buf4_init(&lambda, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                           ID("[O,o]"), ID("[V,v]"), 0, "Lambda <Oo|Vv>");
    global_dpd_->file4_init(&tau_temp, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[V,v]"), "DF tau(temp) <Oo|Vv>");
    // Check:
    //    global_dpd_->file4_print(&tau_temp, "outfile");

    // Sizing
    int navir_a = Cavir_a_->colspi()[0];
    int navir_b = Cavir_b_->colspi()[0];
    int naux = ribasis_->nbf();

    // Read Qac in memory. Will treat memory limit later ...
    SharedMatrix Qac (new Matrix("Qac", navir_a * (ULI) navir_a, naux));
    SharedMatrix Qbd (new Matrix("Qbd", navir_b * (ULI) navir_b, naux));
    double** Qacp = Qac->pointer();
    double** Qbdp = Qbd->pointer();

    // Read Qac all
    psio_->open(PSIF_DFMP2_AIA, PSIO_OPEN_OLD);
    psio_->open(PSIF_DFMP2_QIA, PSIO_OPEN_OLD);
    psio_->read_entry(PSIF_DFMP2_AIA, "(Q|ac)", (char*) Qacp[0], sizeof(double) * (ULI) navir_a * navir_a * naux);
    psio_->read_entry(PSIF_DFMP2_QIA, "(Q|ac)", (char*) Qbdp[0], sizeof(double) * (ULI) navir_b * navir_b * naux);
    psio_->close(PSIF_DFMP2_AIA, 1);
    psio_->close(PSIF_DFMP2_QIA, 1);

    SharedMatrix Icd = SharedMatrix(new Matrix("Icd", navir_a, navir_b));
    for(int h = 0; h < nirrep_; ++h){
        global_dpd_->file4_mat_irrep_init(&tau_temp, h);
        global_dpd_->buf4_mat_irrep_init(&lambda, h);
        global_dpd_->buf4_mat_irrep_rd(&lambda, h);
        for (int Ij = 0; Ij < tau_temp.params->rowtot[h]; ++Ij){
            for (int Ab = 0; Ab < tau_temp.params->coltot[h]; ++Ab){
                int A = tau_temp.params->colorb[h][Ab][0];
                int b = tau_temp.params->colorb[h][Ab][1];
                int DFA = aVirMap[A];
                int DFb = bVirMap[b];
                double** Icdp = Icd->pointer();
                C_DGEMM('N', 'T', navir_a, navir_b, naux, 1.0, Qacp[DFA*navir_a], naux, Qbdp[DFb*navir_b], naux, 0.0, Icdp[0], navir_b);
                double IjAb = 0.0;
                for (int Cd = 0; Cd < lambda.params->coltot[h]; ++Cd){
                    int C = lambda.params->colorb[h][Cd][0];
                    int d = lambda.params->colorb[h][Cd][1];
                    int DFC = aVirMap[C];
                    int DFd = bVirMap[d];
                    double ACbd = Icdp[DFC][DFd];
                    IjAb += ACbd * lambda.matrix[h][Ij][Cd];
                }
                tau_temp.matrix[h][Ij][Ab] = IjAb;
            }
        }
        global_dpd_->buf4_mat_irrep_close(&lambda,h);
        global_dpd_->file4_mat_irrep_wrt(&tau_temp, h);
        global_dpd_->file4_mat_irrep_close(&tau_temp, h);
    }
    //Check:
//  global_dpd_->file4_print(&tau_temp, "outfile");
    global_dpd_->file4_close(&tau_temp);
    global_dpd_->buf4_close(&lambda);

    /* End of tau(temp)_IjAb */}

    //Check:
//    dpdbuf4 test;
//    global_dpd_->buf4_init(&test, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[V,v]"),
//                           ID("[O,o]"), ID("[V,v]"), 0, "DF tau(temp) <Oo|Vv>");
//    global_dpd_->buf4_print(&test, "outfile", 1);
//    global_dpd_->buf4_close(&test);

}

}}// End of namespace
