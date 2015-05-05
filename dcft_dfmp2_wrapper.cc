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

//#include "dcft_dfmp2.h"
//#include "dcft_DCFTSolver_dfmp2.h"

#include "dcft.h"
#include "defines.h"
#include <vector>
#include <liboptions/liboptions.h>
#include <libpsio/psio.hpp>
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

namespace psi {

namespace dfdcft{

/**
 * Compute the density-fitted MP2 energy as an initial guess
 * This code is responsible for initializing the integral transformation, too.
 */

void DCFTSolver::dfmp2_guess()
{
    dcft_timer_on("DCFTSolver::dfmp2_guess()");

    if (options_.get_str("REFERENCE") == "UHF" || options_.get_str("REFERENCE") == "UKS"){
        dfmp2_common_init();
        dfmp2_print_header();
        dcft_timer_on("DCFTSolver::DFMP2 Singles");
        dfmp2_form_singles();
        dcft_timer_off("DCFTSolver::DFMP2 Singles");

        dcft_timer_on("DCFTSolver::DFMP2 Aia");
        dfmp2_form_Aia();
        dcft_timer_off("DCFTSolver::DFMP2 Aia");

        dcft_timer_on("DCFTSolver::DFMP2 Qia");
        dfmp2_form_Qia();
        dcft_timer_off("DCFTSolver::DFMP2 Qia");

        dcft_timer_on("DCFTSolver::DFMP2 Energy");
        dfmp2_form_energy();
        dcft_timer_off("DCFTSolver::DFMP2 Energy");

        dfmp2_print_energies();
    }
    else {
        throw PSIEXCEPTION("DFMP2 guess: Unrecognized reference. DF-DCFT only works for UHF reference.");
    }

    dcft_timer_off("DCFTSolver::dfmp2_guess()");
}

void DCFTSolver::dfmp2_common_init()
{
    //outfile->Printf( "\n\n\t This is a test for DFMP2 function in DFDCFT plugin ! \n\n"); 

    //Wavefunction::print_ = options_.get_int("PRINT");
    //Wavefunction::debug_ = options_.get_int("DEBUG");

    set_debug(0);

    // Check. Error: Pointer being freed was not allocated. Neither works.
    //set_name("DF-MP2");
    //name_ = "DF-MP2";
    //outfile->Printf( "\tThis is the wavefunction's name: %s", name_.c_str()); 

    energies_["Singles Energy"] = 0.0;
    energies_["Opposite-Spin Energy"] = 0.0;
    energies_["Same-Spin Energy"] = 0.0;
    energies_["Reference Energy"] = reference_wavefunction_->reference_energy();

    //double eSCF = Process::environment.wavefunction()->reference_energy();
    //outfile->Printf( "\n\n\t This is SCF energy: %20.15f \n\n", energy_); 

    sss_ = options_.get_double("MP2_SS_SCALE");
    oss_ = options_.get_double("MP2_OS_SCALE");

    boost::shared_ptr<BasisSetParser> parser(new Gaussian94BasisSetParser());
    ribasis_ = BasisSet::construct(parser, molecule_, "DF_BASIS_MP2");

    Caocc_a_ = Ca_subset("AO","ACTIVE_OCC");
    Cavir_a_ = Ca_subset("AO","ACTIVE_VIR");
    Caocc_b_ = Cb_subset("AO","ACTIVE_OCC");
    Cavir_b_ = Cb_subset("AO","ACTIVE_VIR");

    eps_aocc_a_ = epsilon_a_subset("AO","ACTIVE_OCC");
    eps_avir_a_ = epsilon_a_subset("AO","ACTIVE_VIR");
    eps_aocc_b_ = epsilon_b_subset("AO","ACTIVE_OCC");
    eps_avir_b_ = epsilon_b_subset("AO","ACTIVE_VIR");

    aOccOrbsPI = new int [nirrep_];
    aVirOrbsPI = new int [nirrep_];
    bOccOrbsPI = new int [nirrep_];
    bVirOrbsPI = new int [nirrep_];
    for (int h = 0; h < nirrep_; ++h){
        aOccOrbsPI[h] = doccpi_[h] + soccpi_[h] - frzcpi_[h];
        bOccOrbsPI[h] = doccpi_[h] - frzcpi_[h];
        aVirOrbsPI[h] = nmopi_[h] - doccpi_[h] - soccpi_[h] - frzvpi_[h];
        bVirOrbsPI[h] = nmopi_[h] - doccpi_[h] -frzvpi_[h];
    }

}

void DCFTSolver::dfmp2_print_header()
{
    outfile->Printf( "\t --------------------------------------------------------\n");
    outfile->Printf( "\t                DF-MP2 Module in DF-DCFT                 \n");
    outfile->Printf( "\t      2nd-Order Density-Fitted Moller-Plesset Theory     \n");
    outfile->Printf( "\t                    UMP2 Wavefunction                    \n");
    outfile->Printf( "\t                                                         \n");
    outfile->Printf( "\t        Rob Parrish, Justin Turney, Andy Simmonett,      \n");
    outfile->Printf( "\t           Ed Hohenstein, and C. David Sherrill          \n");
    outfile->Printf( "\t   and Xiao (floccinaucinihipilificatious contribution   \n");
    outfile->Printf( "\t --------------------------------------------------------\n");
    outfile->Printf( "\n");

    int focc_a = frzcpi_.sum();
    int fvir_a = frzvpi_.sum();
    int aocc_a = Caocc_a_->colspi()[0];
    int avir_a = Cavir_a_->colspi()[0];
    int occ_a = focc_a + aocc_a;
    int vir_a = fvir_a + avir_a;

    int focc_b = frzcpi_.sum();
    int fvir_b = frzvpi_.sum();
    int aocc_b = Caocc_b_->colspi()[0];
    int avir_b = Cavir_b_->colspi()[0];
    int occ_b = focc_b + aocc_b;
    int vir_b = fvir_b + avir_b;

    outfile->Printf( "\t --------------------------------------------------------\n");
    outfile->Printf( "\t                 NBF = %5d, NAUX = %5d\n", basisset_->nbf(), ribasis_->nbf());
    outfile->Printf( "\t --------------------------------------------------------\n");
    outfile->Printf( "\t %7s %7s %7s %7s %7s %7s %7s\n", "CLASS", "FOCC", "OCC", "AOCC", "AVIR", "VIR", "FVIR");
    outfile->Printf( "\t %7s %7d %7d %7d %7d %7d %7d\n", "ALPHA", focc_a, occ_a, aocc_a, avir_a, vir_a, fvir_a);
    outfile->Printf( "\t %7s %7d %7d %7d %7d %7d %7d\n", "BETA", focc_b, occ_b, aocc_b, avir_b, vir_b, fvir_b);
    outfile->Printf( "\t --------------------------------------------------------\n\n");

}

void DCFTSolver::dfmp2_form_singles()
{
    double E_singles_a = 0.0;
    double E_singles_b = 0.0;

    SharedMatrix Caocc_a = Ca_subset("SO","ACTIVE_OCC");
    SharedMatrix Cavir_a = Ca_subset("SO","ACTIVE_VIR");
    SharedMatrix Caocc_b = Cb_subset("SO","ACTIVE_OCC");
    SharedMatrix Cavir_b = Cb_subset("SO","ACTIVE_VIR");

    SharedVector eps_aocc_a = epsilon_a_subset("SO","ACTIVE_OCC");
    SharedVector eps_avir_a = epsilon_a_subset("SO","ACTIVE_VIR");
    SharedVector eps_aocc_b = epsilon_b_subset("SO","ACTIVE_OCC");
    SharedVector eps_avir_b = epsilon_b_subset("SO","ACTIVE_VIR");

    SharedMatrix Fia_a(new Matrix("Fia a", Caocc_a->colspi(), Cavir_a->colspi()));
    SharedMatrix Fia_b(new Matrix("Fia b", Caocc_b->colspi(), Cavir_b->colspi()));


    double* temp = new double[Fa_->max_nrow() * (ULI) (Cavir_a->max_ncol() > Cavir_b->max_ncol() ? Cavir_a->max_ncol() : Cavir_b->max_ncol())];

    // Fia a
    for (int h = 0; h < Caocc_a->nirrep(); h++) {

        int nso = Fa_->rowspi()[h];
        int naocc = Caocc_a->colspi()[h];
        int navir = Cavir_a->colspi()[h];

        if (!nso || !naocc || !navir) continue;

        double** Fsop = Fa_->pointer(h);
        double** Fmop = Fia_a->pointer(h);
        double** Cip = Caocc_a->pointer(h);
        double** Cap = Cavir_a->pointer(h);

        C_DGEMM('N','N',nso,navir,nso,1.0,Fsop[0],nso,Cap[0],navir,0.0,temp,navir);
        C_DGEMM('T','N',naocc,navir,nso,1.0,Cip[0],naocc,temp,navir,0.0,Fmop[0],navir);

        double* eps_i = eps_aocc_a->pointer(h);
        double* eps_a = eps_avir_a->pointer(h);

        for (int i = 0; i < naocc; i++) {
            for (int a = 0; a < navir; a++) {
                E_singles_a -= Fmop[i][a] * Fmop[i][a] / (eps_a[a] - eps_i[i]);
            }
        }
    }

    // Fia b
    for (int h = 0; h < Caocc_b->nirrep(); h++) {

        int nso = Fb_->rowspi()[h];
        int naocc = Caocc_b->colspi()[h];
        int navir = Cavir_b->colspi()[h];

        if (!nso || !naocc || !navir) continue;

        double** Fsop = Fb_->pointer(h);
        double** Fmop = Fia_b->pointer(h);
        double** Cip = Caocc_b->pointer(h);
        double** Cap = Cavir_b->pointer(h);

        double* eps_i = eps_aocc_b->pointer(h);
        double* eps_a = eps_avir_b->pointer(h);

        C_DGEMM('N','N',nso,navir,nso,1.0,Fsop[0],nso,Cap[0],navir,0.0,temp,navir);
        C_DGEMM('T','N',naocc,navir,nso,1.0,Cip[0],naocc,temp,navir,0.0,Fmop[0],navir);

        for (int i = 0; i < naocc; i++) {
            for (int a = 0; a < navir; a++) {
                E_singles_b -= Fmop[i][a] * Fmop[i][a] / (eps_a[a] - eps_i[i]);
            }
        }
    }

    delete[] temp;

    energies_["Singles Energy"] = E_singles_a + E_singles_b;

    if (debug_) {
        outfile->Printf( "\n\t==================== dfmp2_form_singles part ====================\n");
        Caocc_a->print();
        Cavir_a->print();
        eps_aocc_a->print();
        eps_avir_a->print();
        Caocc_b->print();
        Cavir_b->print();
        eps_aocc_b->print();
        eps_avir_b->print();

        Fia_a->print();
        Fia_b->print();
        outfile->Printf( "  Alpha singles energy = %24.16E\n", E_singles_a);
        outfile->Printf( "  Beta  singles energy = %24.16E\n\n", E_singles_b);
    }

}

void DCFTSolver::dfmp2_form_Aia()
{
    // Schwarz Sieve
    boost::shared_ptr<ERISieve> sieve(new ERISieve(basisset_, options_.get_double("INTS_TOLERANCE")));
    double test_tol = options_.get_double("INTS_TOLERANCE");
    const std::vector<std::pair<int,int> >& shell_pairs = sieve->shell_pairs();
    const size_t npairs = shell_pairs.size();

    // ERI objects
    // Parallel version of DF-MP2 can be added later
    int nthread = 1;

    boost::shared_ptr<IntegralFactory> factory(new IntegralFactory(ribasis_,BasisSet::zero_ao_basis_set(), basisset_, basisset_));
    std::vector<boost::shared_ptr<TwoBodyAOInt> > eri;
    std::vector<const double*> buffer;
    for (int thread = 0; thread < nthread; thread++) {
        eri.push_back(boost::shared_ptr<TwoBodyAOInt>(factory->eri()));
        buffer.push_back(eri[thread]->buffer());
    }

    // Sizing
    int nso = basisset_->nbf();
    int naux = ribasis_->nbf();
    int naocc_a = Caocc_a_->colspi()[0];
    int navir_a = Cavir_a_->colspi()[0];
    int naocc_b = Caocc_b_->colspi()[0];
    int navir_b = Cavir_b_->colspi()[0];
    int naocc = (naocc_a > naocc_b ? naocc_a : naocc_b);
    int navir = (navir_a > navir_b ? navir_a : navir_b);
    int maxQ = ribasis_->max_function_per_shell();

    // Max block size in naux
    ULI Amn_cost_per_row = nso * (ULI) nso;
    ULI Ami_cost_per_row = nso * (ULI) naocc;
    ULI Aia_cost_per_row = naocc * (ULI) navir;
    ULI total_cost_per_row = Amn_cost_per_row + Ami_cost_per_row + Aia_cost_per_row;
    ULI doubles = ((ULI) (options_.get_double("DFMP2_MEM_FACTOR") * memory_ / 8L));
    ULI max_temp = doubles / (total_cost_per_row);
    int max_naux = (max_temp > (ULI) naux ? naux : max_temp);
    max_naux = (max_naux < maxQ ? maxQ : max_naux);

    // Block extents
    std::vector<int> block_Q_starts;
    int counter = 0;
    block_Q_starts.push_back(0);
    for (int Q = 0; Q < ribasis_->nshell(); Q++) {
        int nQ = ribasis_->shell(Q).nfunction();
        if (counter + nQ > max_naux) {
            counter = 0;
            block_Q_starts.push_back(Q);
        }
        counter += nQ;
    }
    block_Q_starts.push_back(ribasis_->nshell());
    //block_status(block_Q_starts, __FILE__,__LINE__);

    // Tensor blocks
    SharedMatrix Amn(new Matrix("(A|mn) Block", max_naux, nso * (ULI) nso));
    SharedMatrix Ami(new Matrix("(A|mi) Block", max_naux, nso * (ULI) naocc));
    SharedMatrix Aia(new Matrix("(A|ia) Block", max_naux, naocc * (ULI) navir));
    double** Amnp = Amn->pointer();
    double** Amip = Ami->pointer();
    double** Aiap = Aia->pointer();

    // C Matrices
    double** Caoccap = Caocc_a_->pointer();
    double** Cavirap = Cavir_a_->pointer();
    double** Caoccbp = Caocc_b_->pointer();
    double** Cavirbp = Cavir_b_->pointer();

    psio_->open(PSIF_DFMP2_AIA,PSIO_OPEN_NEW);
    psio_address next_AIA = PSIO_ZERO;
    psio_->open(PSIF_DFMP2_QIA,PSIO_OPEN_NEW);
    psio_address next_QIA = PSIO_ZERO;

    // Loop over blocks of Qshell
    for (int block = 0; block < block_Q_starts.size() - 1; block++) {

        // Block sizing/offsets
        int Qstart = block_Q_starts[block];
        int Qstop  = block_Q_starts[block+1];
        int qoff   = ribasis_->shell(Qstart).function_index();
        int nrows  = (Qstop == ribasis_->nshell() ?
                     ribasis_->nbf() - ribasis_->shell(Qstart).function_index() :
                     ribasis_->shell(Qstop).function_index() - ribasis_->shell(Qstart).function_index());

        // Clear Amn for Schwarz sieve
        ::memset((void*) Amnp[0], '\0', sizeof(double) * nrows * nso * nso);

        // Compute TEI tensor block (A|mn)
        dcft_timer_on("DCFTSolver::DFMP2 (A|mn)");
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
        dcft_timer_off("DCFTSolver::DFMP2 (A|mn)");

        // => Alpha Case <= //

        // Compute (A|mi) tensor block (A|mn) C_ni
        dcft_timer_on("DCFTSolver::DFMP2 (A|mn)C_mi");
        C_DGEMM('N','N',nrows*(ULI)nso,naocc_a,nso,1.0,Amnp[0],nso,Caoccap[0],naocc_a,0.0,Amip[0],naocc);
        dcft_timer_off("DCFTSolver::DFMP2 (A|mn)C_mi");

        // Compute (A|ia) tensor block (A|ia) = (A|mi) C_ma
        dcft_timer_on("DCFTSolver::DFMP2 (A|mi)C_na");
        for (int row = 0; row < nrows; row++) {
            C_DGEMM('T','N',naocc_a,navir_a,nso,1.0,Amip[row],naocc,Cavirap[0],navir_a,0.0,&Aiap[0][row*(ULI)naocc_a*navir_a],navir_a);
        }
        dcft_timer_off("DCFTSolver::DFMP2 (A|mi)C_na");

        // Stripe (A|ia) out to disk
        dcft_timer_on("DCFTSolver::DFMP2 Aia Write");
        psio_->write(PSIF_DFMP2_AIA, "(A|ia)", (char*)Aiap[0], sizeof(double)*nrows*naocc_a*navir_a, next_AIA, &next_AIA);
        dcft_timer_off("DCFTSolver::DFMP2 Aia Write");

        // => Beta Case <= //

        // Compute (A|mi) tensor block (A|mn) C_ni
        dcft_timer_on("DCFTSolver::DFMP2 (A|mn)C_mi");
        C_DGEMM('N','N',nrows*(ULI)nso,naocc_b,nso,1.0,Amnp[0],nso,Caoccbp[0],naocc_b,0.0,Amip[0],naocc);
        dcft_timer_off("DCFTSolver::DFMP2 (A|mn)C_mi");

        // Compute (A|ia) tensor block (A|ia) = (A|mi) C_ma
        dcft_timer_on("DCFTSolver::DFMP2 (A|mi)C_na");
        for (int row = 0; row < nrows; row++) {
            C_DGEMM('T','N',naocc_b,navir_b,nso,1.0,Amip[row],naocc,Cavirbp[0],navir_b,0.0,&Aiap[0][row*(ULI)naocc_b*navir_b],navir_b);
        }
        dcft_timer_off("DCFTSolver::DFMP2 (A|mi)C_na");

        // Stripe (A|ia) out to disk
        dcft_timer_on("DCFTSolver::DFMP2 Aia Write");
        psio_->write(PSIF_DFMP2_QIA,"(A|ia)",(char*)Aiap[0],sizeof(double)*nrows*naocc_b*navir_b,next_QIA,&next_QIA); // Save (A|ia) beta-beta on PSIF_DFMP2_QIA
        dcft_timer_off("DCFTSolver::DFMP2 Aia Write");
    }

    psio_->close(PSIF_DFMP2_AIA,1);
    psio_->close(PSIF_DFMP2_QIA,1);

    if (debug_) {
         outfile->Printf( "\n\t==================== dfmp2_form_Aia part ====================\n");
         outfile->Printf( "\n\tINTS_TOLERENCE  = %20.15f\n", test_tol);
         outfile->Printf( "\n\t%8s %8s %8s %8s %8s %8s %8s %8s %8s %12s \n", "nso", "naux", "naocc_a", "navir_a", "naocc_b", "navir_b", "naocc", "navir", "maxQ", "nshell_aux");
         outfile->Printf(   "\t%8d %8d %8d %8d %8d %8d %8d %8d %8d %12d \n", nso, naux, naocc_a, navir_a, naocc_b, navir_b, naocc, navir, maxQ, ribasis_->nshell());
    }
}


void DCFTSolver::dfmp2_form_Qia()
{
    SharedMatrix Jm12 = form_inverse_metric();
    apply_fitting(Jm12, PSIF_DFMP2_AIA, "(A|ia)", "(Q|ia)", ribasis_->nbf(), Caocc_a_->colspi()[0] * (ULI) Cavir_a_->colspi()[0]);
    apply_fitting(Jm12, PSIF_DFMP2_QIA, "(A|ia)", "(Q|ia)", ribasis_->nbf(), Caocc_b_->colspi()[0] * (ULI) Cavir_b_->colspi()[0]);

    // Check:
//    psio_->open(PSIF_DFMP2_QIA, PSIO_OPEN_OLD);
//    psio_->open(PSIF_DFMP2_AIA, PSIO_OPEN_OLD);
//    psio_->tocprint(PSIF_DFMP2_QIA);
//    psio_->tocprint(PSIF_DFMP2_AIA);
//    psio_->close(PSIF_DFMP2_QIA, 1);
//    psio_->close(PSIF_DFMP2_AIA, 1);
}

void DCFTSolver::dfmp2_form_energy()
{
    // Energy registers
    double e_ss = 0.0;
    double e_os = 0.0;

    // Open PSIF_MO_xx_TEI files to save transformed tei
    psio_->open(PSIF_MO_AA_TEI, PSIO_OPEN_OLD);
    psio_->open(PSIF_MO_AB_TEI, PSIO_OPEN_OLD);
    psio_->open(PSIF_MO_BB_TEI, PSIO_OPEN_OLD);
    //psio_->tocprint(PSIF_MO_AA_TEI);
    //psio_->tocprint(PSIF_MO_AB_TEI);
    //psio_->tocprint(PSIF_MO_BB_TEI);
    psio_address next_MO_AA = PSIO_ZERO;
    psio_address next_MO_AB = PSIO_ZERO;
    psio_address next_MO_BB = PSIO_ZERO;

    // Check: Open PSIF_LIBTRANS_DPD
    if(0){
        psio_->open(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);
        psio_->tocprint(PSIF_LIBTRANS_DPD);
        psio_->close(PSIF_LIBTRANS_DPD, 1);
    }

    // Check: Try to find some information about symmetry of AOs/MOs
    if(0){
        nsopi_.print();
        nmopi_.print();
        doccpi_.print();
        soccpi_.print();

        std::vector<std::pair<double, int> > aPairs;
        std::vector<std::pair<double, int> > bPairs;

        for(int h = 0; h < nirrep_; ++h){
            for (int i = 0; i < nsopi_[h]; ++i){
                aPairs.push_back(std::make_pair(epsilon_a_->get(h, i), h));
                bPairs.push_back(std::make_pair(epsilon_b_->get(h, i), h));
            }
        }
        sort(aPairs.begin(), aPairs.end());
        sort(bPairs.begin(), bPairs.end());

        char **irrepLabels = chkpt_->rd_irr_labs();
        outfile->Printf( "\n\n\t%10s %10s %10s\n", "aOrbital", "Irrep", "Energy");
        for(int i = 0, count = 0; i < nmo_; ++i, ++count){
            int irrep = aPairs[i].second;
            outfile->Printf( "\t%10d %10s %10.5f\n", i+1, irrepLabels[irrep], aPairs[i].first);
            if (count % 4 == 3 && i != nmo_)
                outfile->Printf( "\n");
        }
    }

    /* Check: Initialize the integral transformation object, _ints,
     * and assign dpd ID
     */
    if(1){
        // Set up the occupied and virtual space
        std::vector<boost::shared_ptr<MOSpace> > spaces;
        spaces.push_back(MOSpace::occ);
        spaces.push_back(MOSpace::vir);
        // This wavefunction is really the global reference wavefunction
        boost::shared_ptr<Wavefunction> wfn = Process::environment.wavefunction();
        _ints = new IntegralTransform(wfn, spaces, IntegralTransform::Unrestricted);
        _ints->set_keep_iwl_so_ints(true);
        _ints->set_keep_dpd_so_ints(true);
        dpd_set_default(_ints->get_dpd_id());
    }

    // Check: Open a new dpdfile2 object with correct dimensions and zero values
    if(0){
        dpdfile2 DF_T;
        global_dpd_->file2_init(&DF_T, PSIF_DCFT_DPD, 0, ID('O'), ID('O'), "DF Tau <O|O>");
        global_dpd_->file2_print(&DF_T, "outfile");
        global_dpd_->file2_close(&DF_T);
    }

    /* Check: Open a new dpdfile4 object with correct dimensions and zero values
     *  then, assign values
     */
    if(0){
        dpdfile4 M;

        psio_->open(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);
        psio_->tocprint(PSIF_LIBTRANS_DPD);

        global_dpd_->file4_init(&M, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[O,V]"), "DF MO Ints (OV|OV)");

        global_dpd_->file4_print(&M, "outfile");

        /* (OV|OV) block */
        outfile->Printf( "\n\t%2s %2s %2s %2s %20s %20s\n\n", "I", "A", "J", "B", "(IA|JB)", "Denom_IJ^AB");
        for(int h = 0; h < nirrep_; ++h){
            global_dpd_->file4_mat_irrep_init(&M, h);
            outfile->Printf( "\n\tHere comes the file4_mat_irrep_rd() function.\n");
//            global_dpd_->file4_mat_irrep_rd(&M, h);
            for(int IA = 0; IA < M.params->rowtot[h]; ++IA){
                int I = M.params->roworb[h][IA][0];
                int A = M.params->roworb[h][IA][1];
                for (int JB = 0; JB < M.params->coltot[h]; ++JB){
                    int J = M.params->colorb[h][JB][0];
                    int B = M.params->colorb[h][JB][1];
                    M.matrix[h][IA][JB] = 1.0 ;
                    outfile->Printf( "\t%2d %2d %2d %2d %20.10f %20.10f\n", I, A, J, B, M.matrix[h][IA][JB], 1.0/(epsilon_a_->get(h, I) + epsilon_a_->get(h, J) - epsilon_a_->get(h, A) - epsilon_a_->get(h, B)));
                }
            }
            global_dpd_->file4_mat_irrep_wrt(&M, h);
//            global_dpd_->file4_print(&M, "outfile");
            global_dpd_->file4_mat_irrep_close(&M, h);
        }

        global_dpd_->file4_print(&M, "outfile");

        global_dpd_->file4_close(&M);
        psio_->tocprint(PSIF_LIBTRANS_DPD);
        psio_->close(PSIF_LIBTRANS_DPD, 1);
    }

    /*
     * Compute lambda_ijab: mp2 guess
     */
    if (1) {
        outfile->Printf( "\tComputing MP2 amplitude guess...\n\n");  

        /* Make DPD - Non-DPD orbital pairs */
        map_orbitals();

        /*
        * L_ijab = <ij||ab> / D_ijab
        */

        /* => L_IJAB <= */{
        // L_IJAB = <IJ||AB> / D_IJAB = ((IA|Q)(Q|JB) - (IB|Q)(Q|JA)) / D_IJAB

        dpdfile4 L;
        global_dpd_->file4_init(&L, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"), "DF Lambda <OO|VV>");
//        global_dpd_->file4_print(&L, "outfile");

        // Sizing
        int naocc = Caocc_a_->colspi()[0];
        int navir = Cavir_a_->colspi()[0];
        int naux  = ribasis_->nbf();

        // Read Qia into memory. Will treat memory limit later
        SharedMatrix Qia (new Matrix("Qia0", naocc * (ULI) navir, naux));
        SharedMatrix Qjb (new Matrix("Qjb0", naocc * (ULI) navir, naux));
        double ** Qiap = Qia->pointer();
        double ** Qjbp = Qjb->pointer();

        // Read Qia all
        psio_->open(PSIF_DFMP2_AIA, PSIO_OPEN_OLD);
        psio_->read_entry(PSIF_DFMP2_AIA, "(Q|ia)", (char*) Qiap[0], sizeof(double) * (ULI) naocc * navir * naux);
        psio_->read_entry(PSIF_DFMP2_AIA, "(Q|ia)", (char*) Qjbp[0], sizeof(double) * (ULI) naocc * navir * naux);
        psio_->close(PSIF_DFMP2_AIA, 1);

        // Orbital energies
        double* eps_aoccp = eps_aocc_a_->pointer();
        double* eps_avirp = eps_avir_a_->pointer();

        // Check: orbital orders
        if(0){
            // DPD form
            nsopi_.print();
            std::vector<std::pair<double, int> > aPairs;
            for(int h = 0; h < nirrep_; ++h){
                for (int i = 0; i < nsopi_[h]; ++i){
                    aPairs.push_back(std::make_pair(epsilon_a_->get(h, i), h));
                }
            }
            //sort(aPairs.begin(), aPairs.end());
            char **irrepLabels = chkpt_->rd_irr_labs();
            outfile->Printf( "\n\n\t%10s %10s %10s\n", "aOrbital", "Irrep", "Energy");
            for(int i = 0, count = 0; i < nmo_; ++i, ++count){
                int irrep = aPairs[i].second;
                outfile->Printf( "\t%10d %10s %10.5f\n", i+1, irrepLabels[irrep], aPairs[i].first);
                if (count % 4 == 3 && i != nmo_)
                    outfile->Printf( "\n");
            }
            // DF form
            outfile->Printf( "\n\n\t%10s %10s %10s\n", "aOcc", " ", "Energy");
            for(int i = 0; i < naocc; ++i){
                outfile->Printf( "\t%10d %10s %10.5f\n", i+1, " ", eps_aoccp[i]);
            }
            outfile->Printf( "\n\n\t%10s %10s %10s\n", "aVir", " ", "Energy");
            for(int a = 0; a < navir; ++a){
                outfile->Printf( "\t%10d %10s %10.5f\n", a+1, " ", eps_avirp[a]);
            }
        }

        // Assign MP2-guess values for lambda matrix
        SharedMatrix Iab = SharedMatrix(new Matrix("Iab0", navir, navir));
//        outfile->Printf( "\n\t%5s %5s %5s %5s %5s %5s %5s %5s \n", "I", "J", "A", "B", "DFI", "DFJ", "DFA", "DFB");
        for(int h = 0; h < nirrep_; ++h){
            global_dpd_->file4_mat_irrep_init(&L, h);
            for ( int IJ = 0; IJ < L.params->rowtot[h]; ++IJ){
                int I = L.params->roworb[h][IJ][0];
                int J = L.params->roworb[h][IJ][1];
                int DFI = aOccMap[I];
                int DFJ = aOccMap[J];
                double** Iabp = Iab->pointer();
                C_DGEMM('N', 'T', navir, navir, naux, 1.0, Qiap[DFI*navir], naux, Qjbp[DFJ*navir], naux, 0.0, Iabp[0], navir);
                //outfile->Printf( "\tNumber of SOs in Irrep_[%2d]: %d\n", h, DFi);
                for (int AB = 0; AB < L.params->coltot[h]; ++AB){
                    int A = L.params->colorb[h][AB][0];
                    int B = L.params->colorb[h][AB][1];
                    int DFA = aVirMap[A];
                    int DFB = aVirMap[B];
                    double IAJB = Iabp[DFA][DFB];
                    double IBJA = Iabp[DFB][DFA];
                    double denom = 1.0 / (eps_aoccp[DFI] + eps_aoccp[DFJ] - eps_avirp[DFA] - eps_avirp[DFB]);
                    L.matrix[h][IJ][AB] = (IAJB - IBJA) * denom;
//                    outfile->Printf( "\t%5d %5d %5d %5d %5d %5d %5d %5d \n", I, J, A, B, DFI, DFJ, DFA, DFB);
                }
            }
            global_dpd_->file4_mat_irrep_wrt(&L, h);
            global_dpd_->file4_mat_irrep_close(&L, h);
        }
//        global_dpd_->file4_print(&L, "outfile");
        global_dpd_->file4_close(&L);

        /* End of L_IJAB */}



        /* => L_IjAb <= */ {
        // L_IjAb = <Ij|Ab> / D_IjAb = (IA|Q)(Q|jb) / D_IjAb

        dpdfile4 L;
        global_dpd_->file4_init(&L, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[V,v]"), "DF Lambda <Oo|Vv>");
//        global_dpd_->file4_print(&L, "outfile");

        // Sizing
        int naocc_a = Caocc_a_->colspi()[0];
        int naocc_b = Caocc_b_->colspi()[0];
        int navir_a = Cavir_a_->colspi()[0];
        int navir_b = Cavir_b_->colspi()[0];
        int naux  = ribasis_->nbf();

        // Read Qia into memory. Will treat memory limit later
        SharedMatrix Qia (new Matrix("Qia", naocc_a * (ULI) navir_a, naux));
        SharedMatrix Qjb (new Matrix("Qjb", naocc_b * (ULI) navir_b, naux));
        double ** Qiap = Qia->pointer();
        double ** Qjbp = Qjb->pointer();

        // Read Qia all
        psio_->open(PSIF_DFMP2_AIA, PSIO_OPEN_OLD);
        psio_->open(PSIF_DFMP2_QIA, PSIO_OPEN_OLD);
        psio_->read_entry(PSIF_DFMP2_AIA, "(Q|ia)", (char*) Qiap[0], sizeof(double) * (ULI) naocc_a * navir_a * naux);
        psio_->read_entry(PSIF_DFMP2_QIA, "(Q|ia)", (char*) Qjbp[0], sizeof(double) * (ULI) naocc_b * navir_b * naux);
        psio_->close(PSIF_DFMP2_AIA, 1);
        psio_->close(PSIF_DFMP2_QIA, 1);

        // Orbital energies
        double* eps_aoccap = eps_aocc_a_->pointer();
        double* eps_aoccbp = eps_aocc_b_->pointer();
        double* eps_avirap = eps_avir_a_->pointer();
        double* eps_avirbp = eps_avir_b_->pointer();

        // Assign MP2-guess values for lambda matrix
        SharedMatrix Iab = SharedMatrix(new Matrix("Iab", navir_a, navir_b));
//        outfile->Printf( "\n\t%5s %5s %5s %5s %5s %5s %5s %5s \n", "I", "J", "A", "B", "DFI", "DFJ", "DFA", "DFB");
        for(int h = 0; h < nirrep_; ++h){
            global_dpd_->file4_mat_irrep_init(&L, h);
            for ( int Ij = 0; Ij < L.params->rowtot[h]; ++Ij){
                int I = L.params->roworb[h][Ij][0];
                int j = L.params->roworb[h][Ij][1];
                int DFI = aOccMap[I];
                int DFj = bOccMap[j];
                double** Iabp = Iab->pointer();
                C_DGEMM('N', 'T', navir_a, navir_b, naux, 1.0, Qiap[DFI*navir_a], naux, Qjbp[DFj*navir_b], naux, 0.0, Iabp[0], navir_b);
                //outfile->Printf( "\tNumber of SOs in Irrep_[%2d]: %d\n", h, DFi);
                for (int Ab = 0; Ab < L.params->coltot[h]; ++Ab){
                    int A = L.params->colorb[h][Ab][0];
                    int b = L.params->colorb[h][Ab][1];
                    int DFA = aVirMap[A];
                    int DFb = bVirMap[b];
                    double IAjb = Iabp[DFA][DFb];
                    double IbjA = Iabp[DFb][DFA];
                    double denom = 1.0 / (eps_aoccap[DFI] + eps_aoccbp[DFj] - eps_avirap[DFA] - eps_avirbp[DFb]);
//                    L.matrix[h][Ij][Ab] = (IAjb - IbjA) * denom;
                    L.matrix[h][Ij][Ab] = IAjb * denom;
//                    outfile->Printf( "\t%5d %5d %5d %5d %5d %5d %5d %5d \n", I, J, A, B, DFI, DFJ, DFA, DFB);
                }
            }
            global_dpd_->file4_mat_irrep_wrt(&L, h);
            global_dpd_->file4_mat_irrep_close(&L, h);
        }
//        global_dpd_->file4_print(&L, "outfile");
        global_dpd_->file4_close(&L);

        /* End of L_IjAb */ }


        /* => L_ijab <= */ {
        // L_ijab = <ij||ab> / D_ijab = ((ia|Q)(Q|jb) - (ib|Q)(Q|ja)) / D_ijab

        dpdfile4 L;
        global_dpd_->file4_init(&L, PSIF_DCFT_DPD, 0, ID("[o,o]"), ID("[v,v]"), "DF Lambda <oo|vv>");
//        global_dpd_->file4_print(&L, "outfile");

        // Sizing
        int naocc = Caocc_b_->colspi()[0];
        int navir = Cavir_b_->colspi()[0];
        int naux  = ribasis_->nbf();

        // Read Qia into memory. Will treat memory limit later
        SharedMatrix Qia (new Matrix("Qia", naocc * (ULI) navir, naux));
        SharedMatrix Qjb (new Matrix("Qjb", naocc * (ULI) navir, naux));
        double ** Qiap = Qia->pointer();
        double ** Qjbp = Qjb->pointer();

        // Read Qia all
        psio_->open(PSIF_DFMP2_QIA, PSIO_OPEN_OLD);
        psio_->read_entry(PSIF_DFMP2_QIA, "(Q|ia)", (char*) Qiap[0], sizeof(double) * (ULI) naocc * navir * naux);
        psio_->read_entry(PSIF_DFMP2_QIA, "(Q|ia)", (char*) Qjbp[0], sizeof(double) * (ULI) naocc * navir * naux);
        psio_->close(PSIF_DFMP2_QIA, 1);

        // Orbital energies
        double* eps_aoccp = eps_aocc_b_->pointer();
        double* eps_avirp = eps_avir_b_->pointer();

        // Assign MP2-guess values for lambda matrix
        SharedMatrix Iab = SharedMatrix(new Matrix("Iab", navir, navir));
        for(int h = 0; h < nirrep_; ++h){
            global_dpd_->file4_mat_irrep_init(&L, h);
            for ( int ij = 0; ij < L.params->rowtot[h]; ++ij){
                int i = L.params->roworb[h][ij][0];
                int j = L.params->roworb[h][ij][1];
                int DFi = bOccMap[i];
                int DFj = bOccMap[j];
                double** Iabp = Iab->pointer();
                C_DGEMM('N', 'T', navir, navir, naux, 1.0, Qiap[DFi*navir], naux, Qjbp[DFj*navir], naux, 0.0, Iabp[0], navir);
                //outfile->Printf( "\tNumber of SOs in Irrep_[%2d]: %d\n", h, DFi);
                for (int ab = 0; ab < L.params->coltot[h]; ++ab){
                    int a = L.params->colorb[h][ab][0];
                    int b = L.params->colorb[h][ab][1];
                    int DFa = bVirMap[a];
                    int DFb = bVirMap[b];
                    double iajb = Iabp[DFa][DFb];
                    double ibja = Iabp[DFb][DFa];
                    double denom = 1.0 / (eps_aoccp[DFi] + eps_aoccp[DFj] - eps_avirp[DFa] - eps_avirp[DFb]);
                    L.matrix[h][ij][ab] = (iajb - ibja) * denom;
//                    outfile->Printf( "\t%5d %5d %5d %5d %5d %5d %5d %5d \n", I, J, A, B, DFI, DFJ, DFA, DFB);
                }
            }
            global_dpd_->file4_mat_irrep_wrt(&L, h);
            global_dpd_->file4_mat_irrep_close(&L, h);
        }
//        global_dpd_->file4_print(&L, "outfile");
        global_dpd_->file4_close(&L);

        /* End of L_ijab */}


    }

    /*
     * E = 1/4 L_IJAB * L_IJAB * D_IJAB
     *        +L_IjAb * L_IjAb * D_IjAb
     *    +1/4 L_ijab * L_ijab * D_ijab
     */
//    if(0){
//        dpdbuf4 L, I;
//        // Alpha - Alpha
//        global_dpd_->buf4_init(&L, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
//                               ID("[O,O]"), ID("[V,V]"), 0, "DF Lambda <OO|VV>");
//        global_dpd_->buf4_print(&L, "outfile", 1);
//        global_dpd_->buf4_close(&L);
//    }


    /*
     * Compute MP2 energy
     **/
    if(1){

        /* => AA Terms <= */ {

        // Sizing
        int naocc = Caocc_a_->colspi()[0];
        int navir = Cavir_a_->colspi()[0];
        int naux  = ribasis_->nbf();

//        outfile->Printf( "\n\tnaocc = %d, navir = %d, naux = %d\n", naocc, navir, naux);

        // Read Qia into memory
        // will treat memory limit later
        SharedMatrix Qia (new Matrix("Qia0", naocc * (ULI) navir, naux));
        SharedMatrix Qjb (new Matrix("Qjb0", naocc * (ULI) navir, naux));
        double ** Qiap = Qia->pointer();
        double ** Qjbp = Qjb->pointer();

        // Read Qia all
        psio_->open(PSIF_DFMP2_AIA, PSIO_OPEN_OLD);
        psio_->read_entry(PSIF_DFMP2_AIA, "(Q|ia)", (char*) Qiap[0], sizeof(double) * (ULI) naocc * navir * naux);
        psio_->read_entry(PSIF_DFMP2_AIA, "(Q|ia)", (char*) Qjbp[0], sizeof(double) * (ULI) naocc * navir * naux);
        psio_->close(PSIF_DFMP2_AIA, 1);

        //Check: Qia and Qjb
//        Qia->print();
//        Qjb->print();

        double* eps_aoccp = eps_aocc_a_->pointer();
        double* eps_avirp = eps_avir_a_->pointer();

        // Alpha-Alpha MP2 energy contribution
        double Eaa = 0.0;

        SharedMatrix Iab = SharedMatrix(new Matrix("Iab0", navir, navir));

        for (long int ij = 0L; ij < naocc * naocc; ++ij){
            // Sizing
            ULI i = ij/naocc;
            ULI j = ij%naocc;
            if (j > i) continue;

            double** Iabp = Iab->pointer();

            // Form the integral block (ia|jb) = (ia|Q)(Q|jb)
            C_DGEMM('N', 'T', navir, navir, naux, 1.0, Qiap[i*navir], naux, Qjbp[j*navir], naux, 0.0, Iabp[0], navir);

            double perm_factor = ( i == j ? 1.0 : 2.0);

            for (int a = 0; a < navir; ++a){
                for (int b = 0; b < navir; ++b){
//                    double iajb = 0.0, ibja = 0.0;
//                    for(int Q = 0; Q < naux; ++Q){
//                        iajb += Qiap[i*a][Q] * Qiap[j*b][Q];
//                        ibja += Qiap[i*b][Q] * Qiap[j*a][Q];
//                    }
                    double iajb = Iabp[a][b];
                    double ibja = Iabp[b][a];

                    double denom = - perm_factor / (eps_avirp[a] + eps_avirp[b] - eps_aoccp[i] - eps_aoccp[j]);

                    Eaa += 0.5*(iajb*iajb - iajb*ibja) * denom;
                }
            }
        }

        outfile->Printf( "\n\tThe test MP2 alpha-alpha energy is %20.15f\n", Eaa);

        /* End AA Terms */ }







    }



    /* => AA Terms <= */ {

    // Sizing
    int naux  = ribasis_->nbf();
    int naocc = Caocc_a_->colspi()[0];
    int navir = Cavir_a_->colspi()[0];

    // Thread considerations
    int nthread = 1;

    // Memory
    ULI Iab_memory = navir * (ULI) navir;
    ULI Qa_memory  = naux  * (ULI) navir; // ??
    ULI doubles = ((ULI) (options_.get_double("DFMP2_MEM_FACTOR") * memory_ / 8L));
    if (doubles < nthread * Iab_memory) {
        throw PSIEXCEPTION("DFMP2: Insufficient memory for Iab buffers. Reduce OMP Threads or increase memory.");
    }
    ULI remainder = doubles - nthread * Iab_memory;
    ULI max_i = remainder / (2L * Qa_memory);
    max_i = (max_i > naocc? naocc : max_i);
    max_i = (max_i < 1L ? 1L : max_i);

    // Blocks
    std::vector<ULI> i_starts;
    i_starts.push_back(0L);
    for (ULI i = 0; i < naocc; i += max_i) {
        if (i + max_i >= naocc) {
            i_starts.push_back(naocc);
        } else {
            i_starts.push_back(i + max_i);
        }
    }
    //block_status(i_starts, __FILE__,__LINE__);

    // Tensor blocks
    SharedMatrix Qia (new Matrix("Qia", max_i * (ULI) navir, naux));
    SharedMatrix Qjb (new Matrix("Qjb", max_i * (ULI) navir, naux));
    double** Qiap = Qia->pointer();
    double** Qjbp = Qjb->pointer();

    std::vector<SharedMatrix> Iab;
    for (int i = 0; i < nthread; i++) {
        Iab.push_back(SharedMatrix(new Matrix("Iab",navir,navir)));
    }

    double* eps_aoccp = eps_aocc_a_->pointer();
    double* eps_avirp = eps_avir_a_->pointer();

    // Loop through pairs of blocks
    psio_->open(PSIF_DFMP2_AIA,PSIO_OPEN_OLD);
    psio_address next_AIA = PSIO_ZERO;

    for (int block_i = 0; block_i < i_starts.size() - 1; block_i++) {

        // Sizing
        ULI istart = i_starts[block_i];
        ULI istop  = i_starts[block_i+1];
        ULI ni     = istop - istart;

        // Read iaQ chunk
        dcft_timer_on("DCFTSolver::DFMP2 iaQ alpha-alpha Read");
        next_AIA = psio_get_address(PSIO_ZERO,sizeof(double)*(istart * navir * naux));
        psio_->read(PSIF_DFMP2_AIA,"(Q|ia)",(char*)Qiap[0],sizeof(double)*(ni * navir * naux),next_AIA,&next_AIA);
        dcft_timer_off("DCFTSolver::DFMP2 iaQ alpha-alpha Read");

        for (int block_j = 0; block_j <= block_i; block_j++) {

            // Sizing
            ULI jstart = i_starts[block_j];
            ULI jstop  = i_starts[block_j+1];
            ULI nj     = jstop - jstart;

            // Read iaQ chunk (if unique)
            dcft_timer_on("DCFTSolver::DFMP2 Qia alpha-alpha Read");
            if (block_i == block_j) {
                ::memcpy((void*) Qjbp[0], (void*) Qiap[0], sizeof(double)*(ni * navir * naux));
            } else {
                next_AIA = psio_get_address(PSIO_ZERO,sizeof(double)*(jstart * navir * naux));
                psio_->read(PSIF_DFMP2_AIA,"(Q|ia)",(char*)Qjbp[0],sizeof(double)*(nj * navir * naux),next_AIA,&next_AIA);
            }
            dcft_timer_off("DCFTSolver::DFMP2 Qia alpha-alpha Read");

            // Check Qia and Qjb
//            Qia->print();
//            Qjb->print();

            for (long int ij = 0L; ij < ni * nj; ij++) {

                // Sizing
                ULI i = ij / nj + istart;
                ULI j = ij % nj + jstart;
                if (j > i) continue;

                double perm_factor = (i == j ? 1.0 : 2.0);

                // Which thread is this?
                int thread = 0;

                double** Iabp = Iab[thread]->pointer();

                // Form the integral block (ia|jb) = (ia|Q)(Q|jb)
                C_DGEMM('N','T',navir,navir,naux,1.0,Qiap[(i-istart)*navir],naux,Qjbp[(j-jstart)*navir],naux,0.0,Iabp[0],navir);

                // Add the MP2 energy contributions
                for (int a = 0; a < navir; a++) {
                    for (int b = 0; b < navir; b++) {
                        double iajb = Iabp[a][b];
                        double ibja = Iabp[b][a];

                        double denom = - perm_factor / (eps_avirp[a] + eps_avirp[b] - eps_aoccp[i] - eps_aoccp[j]);

                        double buf = 0.5*(iajb*iajb - iajb*ibja) * denom;
                        e_ss += buf;

                        //Check: print all iterations
//                        outfile->Printf( "\t%3lu %3lu %3d %3d %20.10f %20.10f\n", i, j, a, b, denom, buf);
                    }
                }

                // Save transformed 2-e integrals (ia|jb) and (ib|ja)
                psio_->write(PSIF_MO_AA_TEI, "<IJ|AB>", (char*)Iabp[0], navir * navir, next_MO_AA, &next_MO_AA);
                //psio_->tocprint(PSIF_MO_AA_TEI);

            }

            // Check: print Alpha-Alpha energy
            outfile->Printf( "\n\tThe reference MP2 alpha-alpha energy is %20.15f\n", e_ss);
        }
    }


    psio_->close(PSIF_DFMP2_AIA,1);

    psio_->close(PSIF_MO_AA_TEI, 1);

    /* End AA Terms */ }

    /* => BB Terms <= */ {

    // Sizing
    int naux  = ribasis_->nbf();
    int naocc = Caocc_b_->colspi()[0];
    int navir = Cavir_b_->colspi()[0];

    // Thread considerations
    int nthread = 1;

    // Memory
    ULI Iab_memory = navir * (ULI) navir;
    ULI Qa_memory  = naux  * (ULI) navir;
    ULI doubles = ((ULI) (options_.get_double("DFMP2_MEM_FACTOR") * memory_ / 8L));
    if (doubles < nthread * Iab_memory) {
        throw PSIEXCEPTION("DFMP2: Insufficient memory for Iab buffers. Reduce OMP Threads or increase memory.");
    }
    ULI remainder = doubles - nthread * Iab_memory;
    ULI max_i = remainder / (2L * Qa_memory);
    max_i = (max_i > naocc? naocc : max_i);
    max_i = (max_i < 1L ? 1L : max_i);

    // Blocks
    std::vector<ULI> i_starts;
    i_starts.push_back(0L);
    for (ULI i = 0; i < naocc; i += max_i) {
        if (i + max_i >= naocc) {
            i_starts.push_back(naocc);
        } else {
            i_starts.push_back(i + max_i);
        }
    }
    //block_status(i_starts, __FILE__,__LINE__);

    // Tensor blocks
    SharedMatrix Qia (new Matrix("Qia", max_i * (ULI) navir, naux));
    SharedMatrix Qjb (new Matrix("Qjb", max_i * (ULI) navir, naux));
    double** Qiap = Qia->pointer();
    double** Qjbp = Qjb->pointer();

    std::vector<SharedMatrix> Iab;
    for (int i = 0; i < nthread; i++) {
        Iab.push_back(SharedMatrix(new Matrix("Iab",navir,navir)));
    }

    double* eps_aoccp = eps_aocc_b_->pointer();
    double* eps_avirp = eps_avir_b_->pointer();

    // Loop through pairs of blocks
    psio_->open(PSIF_DFMP2_QIA,PSIO_OPEN_OLD);
    psio_address next_AIA = PSIO_ZERO;
    for (int block_i = 0; block_i < i_starts.size() - 1; block_i++) {

        // Sizing
        ULI istart = i_starts[block_i];
        ULI istop  = i_starts[block_i+1];
        ULI ni     = istop - istart;

        // Read iaQ chunk
        dcft_timer_on("DCFTSolver::DFMP2 iaQ beta-beta Read");
        next_AIA = psio_get_address(PSIO_ZERO,sizeof(double)*(istart * navir * naux));
        psio_->read(PSIF_DFMP2_QIA,"(Q|ia)",(char*)Qiap[0],sizeof(double)*(ni * navir * naux),next_AIA,&next_AIA);
        dcft_timer_off("DCFTSolver::DFMP2 iaQ beta-beta Read");

        for (int block_j = 0; block_j <= block_i; block_j++) {

            // Sizing
            ULI jstart = i_starts[block_j];
            ULI jstop  = i_starts[block_j+1];
            ULI nj     = jstop - jstart;

            // Read iaQ chunk (if unique)
            dcft_timer_on("DCFTSolver::DFMP2 Qia beta-beta Read");
            if (block_i == block_j) {
                ::memcpy((void*) Qjbp[0], (void*) Qiap[0], sizeof(double)*(ni * navir * naux));
            } else {
                next_AIA = psio_get_address(PSIO_ZERO,sizeof(double)*(jstart * navir * naux));
                psio_->read(PSIF_DFMP2_QIA,"(Q|ia)",(char*)Qjbp[0],sizeof(double)*(nj * navir * naux),next_AIA,&next_AIA);
            }
            dcft_timer_off("DCFTSolver::DFMP2 Qia beta-beta Read");

            for (long int ij = 0L; ij < ni * nj; ij++) {

                // Sizing
                ULI i = ij / nj + istart;
                ULI j = ij % nj + jstart;
                if (j > i) continue;

                double perm_factor = (i == j ? 1.0 : 2.0);

                // Which thread is this?
                int thread = 0;

                double** Iabp = Iab[thread]->pointer();

                // Form the integral block (ia|jb) = (ia|Q)(Q|jb)
                C_DGEMM('N','T',navir,navir,naux,1.0,Qiap[(i-istart)*navir],naux,Qjbp[(j-jstart)*navir],naux,0.0,Iabp[0],navir);

                // Add the MP2 energy contributions
                for (int a = 0; a < navir; a++) {
                    for (int b = 0; b < navir; b++) {
                        double iajb = Iabp[a][b];
                        double ibja = Iabp[b][a];
                        double denom = - perm_factor / (eps_avirp[a] + eps_avirp[b] - eps_aoccp[i] - eps_aoccp[j]);
                        // Add psio_write to save transformed 2-e integrals (ia|jb) and (ib|ja)
                        //
                        e_ss += 0.5*(iajb*iajb - iajb*ibja) * denom;
                    }
                }
            }
        }
    }
    psio_->close(PSIF_DFMP2_QIA,1);

    /* End BB Terms */ }

    /* => AB Terms <= */ {

    // Sizing
    int naux  = ribasis_->nbf();
    int naocc_a = Caocc_a_->colspi()[0];
    int navir_a = Cavir_a_->colspi()[0];
    int naocc_b = Caocc_b_->colspi()[0];
    int navir_b = Cavir_b_->colspi()[0];
    int naocc = (naocc_a > naocc_b ? naocc_a : naocc_b);
    int navir = (navir_a > navir_b ? navir_a : navir_b);

    // Thread considerations
    int nthread = 1;

    // Memory
    ULI Iab_memory = navir_a * (ULI) navir_b;
    ULI Qa_memory  = naux  * (ULI) navir_a;
    ULI Qb_memory  = naux  * (ULI) navir_b;
    ULI doubles = ((ULI) (options_.get_double("DFMP2_MEM_FACTOR") * memory_ / 8L));
    if (doubles < nthread * Iab_memory) {
        throw PSIEXCEPTION("DFMP2: Insufficient memory for Iab buffers. Reduce OMP Threads or increase memory.");
    }
    ULI remainder = doubles - nthread * Iab_memory;
    ULI max_i = remainder / (Qa_memory + Qb_memory);
    max_i = (max_i > naocc? naocc : max_i);
    max_i = (max_i < 1L ? 1L : max_i);

    // Blocks
    std::vector<ULI> i_starts_a;
    i_starts_a.push_back(0L);
    for (ULI i = 0; i < naocc_a; i += max_i) {
        if (i + max_i >= naocc_a) {
            i_starts_a.push_back(naocc_a);
        } else {
            i_starts_a.push_back(i + max_i);
        }
    }
    std::vector<ULI> i_starts_b;
    i_starts_b.push_back(0L);
    for (ULI i = 0; i < naocc_b; i += max_i) {
        if (i + max_i >= naocc_b) {
            i_starts_b.push_back(naocc_b);
        } else {
            i_starts_b.push_back(i + max_i);
        }
    }

    // Tensor blocks
    SharedMatrix Qia (new Matrix("Qia", max_i * (ULI) navir_a, naux));
    SharedMatrix Qjb (new Matrix("Qjb", max_i * (ULI) navir_b, naux));
    double** Qiap = Qia->pointer();
    double** Qjbp = Qjb->pointer();

    std::vector<SharedMatrix> Iab;
    for (int i = 0; i < nthread; i++) {
        Iab.push_back(SharedMatrix(new Matrix("Iab",navir_a,navir_b)));
    }

    double* eps_aoccap = eps_aocc_a_->pointer();
    double* eps_avirap = eps_avir_a_->pointer();
    double* eps_aoccbp = eps_aocc_b_->pointer();
    double* eps_avirbp = eps_avir_b_->pointer();

    // Loop through pairs of blocks
    psio_->open(PSIF_DFMP2_AIA,PSIO_OPEN_OLD);
    psio_->open(PSIF_DFMP2_QIA,PSIO_OPEN_OLD);
    psio_address next_AIA = PSIO_ZERO;
    psio_address next_QIA = PSIO_ZERO;
    for (int block_i = 0; block_i < i_starts_a.size() - 1; block_i++) {

        // Sizing
        ULI istart = i_starts_a[block_i];
        ULI istop  = i_starts_a[block_i+1];
        ULI ni     = istop - istart;

        // Read iaQ chunk
        dcft_timer_on("DCFTSolver::DFMP2 iaQ alpha-beta Read");
        next_AIA = psio_get_address(PSIO_ZERO,sizeof(double)*(istart * navir_a * naux));
        psio_->read(PSIF_DFMP2_AIA,"(Q|ia)",(char*)Qiap[0],sizeof(double)*(ni * navir_a * naux),next_AIA,&next_AIA);
        dcft_timer_off("DCFTSolver::DFMP2 iaQ alpha-beta Read");

        for (int block_j = 0; block_j < i_starts_b.size() - 1; block_j++) {

            // Sizing
            ULI jstart = i_starts_b[block_j];
            ULI jstop  = i_starts_b[block_j+1];
            ULI nj     = jstop - jstart;

            // Read iaQ chunk
            dcft_timer_on("DCFTSolver::DFMP2 Qia alpha-beta Read");
            next_QIA = psio_get_address(PSIO_ZERO,sizeof(double)*(jstart * navir_b * naux));
            psio_->read(PSIF_DFMP2_QIA,"(Q|ia)",(char*)Qjbp[0],sizeof(double)*(nj * navir_b * naux),next_QIA,&next_QIA);
            dcft_timer_off("DCFTSolver::DFMP2 Qia alpha-beta Read");

            for (long int ij = 0L; ij < ni * nj; ij++) {

                // Sizing
                ULI i = ij / nj + istart;
                ULI j = ij % nj + jstart;

                // Which thread is this?
                int thread = 0;

                double** Iabp = Iab[thread]->pointer();

                // Form the integral block (ia|jb) = (ia|Q)(Q|jb)
                C_DGEMM('N','T',navir_a,navir_b,naux,1.0,Qiap[(i-istart)*navir_a],naux,Qjbp[(j-jstart)*navir_b],naux,0.0,Iabp[0],navir_b);

                // Add the MP2 energy contributions
                for (int a = 0; a < navir_a; a++) {
                    for (int b = 0; b < navir_b; b++) {
                        double iajb = Iabp[a][b];
                        double denom = - 1.0 / (eps_avirap[a] + eps_avirbp[b] - eps_aoccap[i] - eps_aoccbp[j]);
                        // Add psio_write to save transformed 2-e integrals (ia|jb) and (ib|ja)
                        //
                        e_os += (iajb*iajb) * denom;
                    }
                }
            }
        }
    }
    psio_->close(PSIF_DFMP2_AIA,1);
    psio_->close(PSIF_DFMP2_QIA,1);

    /* End BB Terms */ }

    energies_["Same-Spin Energy"] = e_ss;
    energies_["Opposite-Spin Energy"] = e_os;

}

void DCFTSolver::dfmp2_print_energies()
{
    energies_["Correlation Energy"] = energies_["Opposite-Spin Energy"] + energies_["Same-Spin Energy"] + energies_["Singles Energy"];
    energies_["Total Energy"] = energies_["Reference Energy"] + energies_["Correlation Energy"];

    outfile->Printf( "\t----------------------------------------------------------\n");
    outfile->Printf( "\t ==================> DF-MP2 Energies <=================== \n");
    outfile->Printf( "\t----------------------------------------------------------\n");
    outfile->Printf( "\t %-25s = %24.16f [H]\n", "Reference Energy",         energies_["Reference Energy"]);
    outfile->Printf( "\t %-25s = %24.16f [H]\n", "Singles Energy",           energies_["Singles Energy"]);
    outfile->Printf( "\t %-25s = %24.16f [H]\n", "Same-Spin Energy",         energies_["Same-Spin Energy"]);
    outfile->Printf( "\t %-25s = %24.16f [H]\n", "Opposite-Spin Energy",     energies_["Opposite-Spin Energy"]);
    outfile->Printf( "\t %-25s = %24.16f [H]\n", "Correlation Energy",       energies_["Correlation Energy"]);
    outfile->Printf( "\t %-25s = %24.16f [H]\n", "Total Energy",             energies_["Total Energy"]);
    outfile->Printf( "\t----------------------------------------------------------\n");

    outfile->Printf( "\n");
     

}

SharedMatrix DCFTSolver::form_inverse_metric()
{
    dcft_timer_on("DCFTSolver::DFMP2 Inverse Metric");

    int naux = ribasis_->nbf();

    // Load inverse metric from the SCF three-index integral file if it exists
    if (options_.get_str("DF_INTS_IO") == "LOAD") {

        SharedMatrix Jm12(new Matrix("SO Basis Fitting Inverse (Eig)", naux, naux));
        outfile->Printf("\t Will attempt to load fitting metric from file %d.\n\n", PSIF_DFSCF_BJ);  
        psio_->open(PSIF_DFSCF_BJ, PSIO_OPEN_OLD);
        psio_->read_entry(PSIF_DFSCF_BJ, "DFMP2 Jm12", (char*) Jm12->pointer()[0], sizeof(double) * naux * naux);
        psio_->close(PSIF_DFSCF_BJ, 1);

        dcft_timer_off("DCFTSolver::DFMP2 Inverse Metric");

        return Jm12;

    } else {

        // Form the inverse metric manually
        boost::shared_ptr<FittingMetric> metric(new FittingMetric(ribasis_, true));
        metric->form_eig_inverse(1.0E-10);
        SharedMatrix Jm12 = metric->get_metric();

        // Save inverse metric to the SCF three-index integral file if it exists
        if (options_.get_str("DF_INTS_IO") == "SAVE") {
            outfile->Printf("\t Will save fitting metric to file %d.\n\n", PSIF_DFSCF_BJ);  
            psio_->open(PSIF_DFSCF_BJ, PSIO_OPEN_OLD);
            psio_->write_entry(PSIF_DFSCF_BJ, "DFMP2 Jm12", (char*) Jm12->pointer()[0], sizeof(double) * naux * naux);

            //psio_address next_BJ = PSIO_ZERO;
            //psio_->write(PSIF_DFSCF_BJ, "DFMP2 Jm12", (char*) Jm12->pointer()[0], sizeof(double) * naux * naux, next_BJ, &next_BJ);
            psio_->close(PSIF_DFSCF_BJ, 1);
        }
        //psio_->tocprint(PSIF_DFSCF_BJ);
        dcft_timer_off("DCFTSolver::DFMP2 Inverse Metric");

        return Jm12;
    }
}

//void DCFTSolver::apply_fitting(SharedMatrix Jm12, unsigned int file, ULI naux, ULI nia)
//{
//    // Memory constraints
//    ULI Jmem = naux * naux;
//    ULI doubles = (ULI) (options_.get_double("DFMP2_MEM_FACTOR") * (memory_ / 8L));
//    if (doubles < 2L * Jmem) {
//        throw PSIEXCEPTION("DFMP2: More memory required for tractable disk transpose");
//    }
//    ULI rem = (doubles - Jmem) / 2L;
//    ULI max_nia = (rem / naux);
//    max_nia = (max_nia > nia ? nia : max_nia);
//    max_nia = (max_nia < 1L ? 1L : max_nia);

//    // Block sizing
//    std::vector<ULI> ia_starts;
//    ia_starts.push_back(0);
//    for (ULI ia = 0L; ia < nia; ia+=max_nia) {
//        if (ia + max_nia >= nia) {
//            ia_starts.push_back(nia);
//        } else {
//            ia_starts.push_back(ia + max_nia);
//        }
//    }
//    //block_status(ia_starts, __FILE__,__LINE__);

//    // Tensor blocks
//    SharedMatrix Aia(new Matrix("Aia", naux, max_nia));
//    SharedMatrix Qia(new Matrix("Qia", max_nia, naux));
//    double** Aiap = Aia->pointer();
//    double** Qiap = Qia->pointer();
//    double** Jp   = Jm12->pointer();

//    // Loop through blocks
//    psio_->open(file, PSIO_OPEN_OLD);
//    psio_address next_AIA = PSIO_ZERO;
//    psio_address next_QIA = PSIO_ZERO;
//    for (int block = 0; block < ia_starts.size() - 1; block++) {

//        // Sizing
//        ULI ia_start = ia_starts[block];
//        ULI ia_stop  = ia_starts[block+1];
//        ULI ncols = ia_stop - ia_start;

//        // Read Aia
//        dcft_timer_on("DCFTSolver::DFMP2 Aia Read");
//        for (ULI Q = 0; Q < naux; Q++) {
//            next_AIA = psio_get_address(PSIO_ZERO,sizeof(double)*(Q*nia+ia_start));
//            psio_->read(file,"(A|ia)",(char*)Aiap[Q],sizeof(double)*ncols,next_AIA,&next_AIA);
//        }
//        dcft_timer_off("DCFTSolver::DFMP2 Aia Read");

//        if (debug_){
//            psio_->tocprint(file);
//        }

//        // Apply Fitting
//        dcft_timer_on("DCFTSolver::DFMP2 (Q|A)(A|ia)");
//        C_DGEMM('T','N',ncols,naux,naux,1.0,Aiap[0],max_nia,Jp[0],naux,0.0,Qiap[0],naux);
//        dcft_timer_off("DCFTSolver::DFMP2 (Q|A)(A|ia)");

//        // Write Qia
//        dcft_timer_on("DCFTSolver::DFMP2 Qia Write");
//        psio_->write(file,"(Q|ia)",(char*)Qiap[0],sizeof(double)*ncols*naux,next_QIA,&next_QIA);
//        dcft_timer_off("DCFTSolver::DFMP2 Qia Write");

//        if (debug_){
//            psio_->tocprint(file);
//        }

//    }
//    psio_->close(file, 1);
//}

void DCFTSolver::block_status(std::vector<int> inds, const char* file, int line)
{
    bool gimp = false;
    if (inds.size() > 2) {
        gimp = ((inds[inds.size() - 1] - inds[inds.size() - 2]) != (inds[1] - inds[0]));
    }
    printf("%s:%d %zu %s %d %d\n", file, line, inds.size(), (gimp ? "GIMP" : "NOT GIMP"), inds[1] - inds[0], inds[inds.size() - 1] - inds[inds.size()-2]);
}
void DCFTSolver::block_status(std::vector<unsigned long int> inds, const char* file, int line)
{
    bool gimp = false;
    if (inds.size() > 2) {
        gimp = ((inds[inds.size() - 1] - inds[inds.size() - 2]) != (inds[1] - inds[0]));
    }
    printf("%s:%d %zu %s %zu %zu\n", file, line, inds.size(), (gimp ? "GIMP" : "NOT GIMP"), inds[1] - inds[0], inds[inds.size() - 1] - inds[inds.size()-2]);
}

void DCFTSolver::map_orbitals()
{
    int naocc = Caocc_a_->colspi()[0];
    int navir = Cavir_a_->colspi()[0];
    int nbocc = Caocc_b_->colspi()[0];
    int nbvir = Cavir_b_->colspi()[0];

    std::vector<std::pair<double, int> > temp;
    // alpha occupied
    for (int h = 0, offset = 0; h < nirrep_; ++h){
        if ( aOccOrbsPI[h] == 0) continue;
        for (int i = 0; i < aOccOrbsPI[h]; ++i){
            int DPDi = i + offset + frzcpi_[h];
            temp.push_back(std::make_pair(epsilon_a_->get(h, i), DPDi));
        }
        offset += aOccOrbsPI[h];
    }
    std::sort(temp.begin(), temp.end());
//    outfile->Printf( "\n\t%15s %15s\n", "DPD-orb", "NonDPD-orb");
    for (int i = 0; i < naocc; ++i){
        for (int j = 0; j < naocc; ++j){
            if(i == temp[j].second){
                aOccMap[i] = j;
            }
        }
//        outfile->Printf( "\t%15d %15d\n", i, aOccMap[i]);
    }
    temp.clear();

    // alpha virtual
    for (int h = 0, offset = 0; h < nirrep_; ++h){
        if (aVirOrbsPI[h] == 0) continue;
        for (int a = 0; a < aVirOrbsPI[h]; ++a){
            int DPDa = a + offset;
            temp.push_back(std::make_pair(epsilon_a_->get(h, aOccOrbsPI[h] + a), DPDa));
        }
        offset += aVirOrbsPI[h];
    }
    std::sort(temp.begin(), temp.end());
//    outfile->Printf( "\n\t%15s %15s\n", "DPD-orb", "NonDPD-orb");
    for (int a = 0; a < navir; ++a){
        for ( int b = 0; b < navir; ++b){
            if (a == temp[b].second){
                aVirMap[a] = b;
            }
        }
//        outfile->Printf( "\t%15d %15d\n", a, aVirMap[a]);
    }
    temp.clear();

    // beta occupied
    for(int h = 0, offset = 0; h < nirrep_; ++h){
        if ( bOccOrbsPI[h] == 0) continue;
        for (int i = 0; i < bOccOrbsPI[h]; ++i){
            int DPDi = i + offset + frzcpi_[h];
            temp.push_back(std::make_pair(epsilon_b_->get(h, i), DPDi));
        }
        offset += bOccOrbsPI[h];
    }
    std::sort(temp.begin(), temp.end());
    for (int i = 0; i < nbocc; ++i){
        for (int j = 0; j < nbocc; ++j){
            if(i == temp[j].second){
                bOccMap[i] = j;
            }
        }
    }
    temp.clear();

    // beta virtual
    for (int h = 0, offset = 0; h < nirrep_; ++h){
        if (bVirOrbsPI[h] == 0) continue;
        for (int a = 0; a < bVirOrbsPI[h]; ++a){
            int DPDa = a + offset;
            temp.push_back(std::make_pair(epsilon_b_->get(h, bOccOrbsPI[h] + a), DPDa));
        }
        offset += bVirOrbsPI[h];
    }
    std::sort(temp.begin(), temp.end());
    for(int a = 0; a < nbvir; ++a){
        for (int b = 0; b < nbvir; ++b){
            if (a == temp[b].second){
                bVirMap[a] = b;
            }
        }
    }
    temp.clear();
}



}}// End of namespace

