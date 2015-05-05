#include <libplugin/plugin.h>
#include <psi4-dec.h>
#include <libparallel/parallel.h>
#include <liboptions/liboptions.h>
#include <libmints/mints.h>
#include <libpsio/psio.hpp>

INIT_PLUGIN

using namespace boost;

namespace psi{ namespace dfdcft {

extern "C"
int read_options(std::string name, Options& options)
{
    if (name == "DFDCFT"|| options.read_globals()) {
        /*- Reference wavefunction type -*/
        options.add_str("REFERENCE", "UHF", "UHF");
        /*- The algorithm to use for the density cumulant and orbital updates in the DCFT energy computation.
          Two-step algorithm (default) is usually more efficient for small
          systems, but for large systems the simultaneous algorithm is recommended.
          In the cases where the convergence problems are encountered (especially
          for highly symmetric systems) QC algorithm can be used. -*/
        options.add_str("ALGORITHM", "SIMULTANEOUS", "TWOSTEP SIMULTANEOUS QC");
        /*- The algorithm to use for the solution of the response equations for the analytic gradients and properties-*/
        options.add_str("RESPONSE_ALGORITHM", "TWOSTEP", "TWOSTEP SIMULTANEOUS");
        /*- Chooses the type of the quadratically-convergent algorithm (effective for ALGORITHM = QC).
          If set to TWOSTEP the Newton-Raphson equations are only solved for the orbital updates,
          the cumulant is updated using the standard Jacobi algorithm. If set to SIMULTANEOUS both cumulant
          and orbitals are updated in a single Newton-Raphson step. -*/
        options.add_str("QC_TYPE", "SIMULTANEOUS", "TWOSTEP SIMULTANEOUS");
        /*- Convergence criterion for the RMS of the residual vector in the density cumulant updates, as well as
          the solution of the density cumulant and orbital response equations. In the orbital updates controls
          the RMS of the SCF error vector -*/
        options.add_double("R_CONVERGENCE", 1e-10);
        /*- Convergence criterion for the density cumulant and orbital guess for the
          variationally orbital-optimized DCFT methods. Currently only available for ALGORITHM = SIMULTANEOUS. -*/
        options.add_double("GUESS_R_CONVERGENCE", 1e-3);
        /*- Maximum number of the macro- or micro-iterations for both the energy and the solution of the response equations -*/
        options.add_int("MAXITER", 40);
        /*- Value of RMS of the density cumulant residual and SCF error vector below which DIIS extrapolation starts.
          Same keyword controls the DIIS extrapolation for the solution of the response equations. -*/
        options.add_double("DIIS_START_CONVERGENCE", 1e-3);
        /*- Maximum number of error vectors stored for DIIS extrapolation !expert-*/
        options.add_int("DIIS_MAX_VECS", 6);
        /*- Minimum number of error vectors stored for DIIS extrapolation !expert-*/
        options.add_int("DIIS_MIN_VECS", 3);
        /*- Controls whether to avoid the AO->MO transformation of the two-electron integrals for the four-virtual case
          (<VV||VV>) by computing the corresponding terms in the AO basis. AO_BASIS = DISK algorithm reduces the memory
          requirements and can significantly reduce the cost of the energy computation if SIMULTANEOUS
          algorithm is used. For the TWOSTEP algorithm, however, AO_BASIS = DISK
          option is not recommended due to the extra I/O. -*/
        options.add_str("AO_BASIS", "DISK", "NONE DISK");
        /*- The amount (percentage) of damping to apply to the orbital update procedure:
          0 will result in a full update, 100 will completely stall the
          update. A value around 20 (which corresponds to 20\% of the previous
          iteration's density being mixed into the current iteration)
          can help in cases where oscillatory convergence is observed. !expert-*/
        options.add_double("DAMPING_PERCENTAGE",0.0);
        /*- The shift applied to the denominator in the density cumulant update iterations !expert-*/
        options.add_double("TIKHONOW_OMEGA", 0.0);
        /*- The shift applied to the denominator in the orbital update iterations !expert-*/
        options.add_double("ORBITAL_LEVEL_SHIFT", 0.0);
        /*- Controls whether to relax the orbitals during the energy computation or not (for debug puproses only).
          For practical applications only the default must be used !expert-*/
        options.add_bool("MO_RELAX", true);
        /*- Controls whether to ignore terms containing non-idempotent contribution to OPDM or not (for debug puproses only).
          For practical applications only the default must be used !expert-*/
        options.add_bool("IGNORE_TAU", false);
        /*- Controls how to cache quantities within the DPD library !expert-*/
        options.add_int("CACHELEVEL", 2);
        /*- Minimum absolute value below which integrals are neglected !expert-*/
        options.add_double("INTS_TOLERANCE", 1e-14);
        /*- Whether to read the orbitals from a previous computation, or to compute
          an MP2 guess !expert -*/
        options.add_str("DCFT_GUESS", "MP2", "CC BCC MP2 DCFT");
        /*- Whether to perform a guess DC-06 or DC-12 computation for ODC-06 or ODC-12 methods, respectively.
          Currently only available for ALGORITHM = SIMULTANEOUS. -*/
        options.add_bool("ODC_GUESS", false);
        /*- Controls whether to relax the guess orbitals by taking the guess density cumulant
          and performing orbital update on the first macroiteration (for ALOGRITHM = TWOSTEP only) !expert-*/
        options.add_bool("RELAX_GUESS_ORBITALS", false);
        /*- Controls whether to include the coupling terms in the DCFT electronic Hessian (for ALOGRITHM = QC
          with QC_TYPE = SIMULTANEOUS only) -*/
        options.add_bool("QC_COUPLING", false);
        /*- Performs stability analysis of the DCFT energy !expert-*/
        options.add_bool("STABILITY_CHECK", false);
        /*- The value of the rms of the residual in Schmidt orthogonalization which is used as a threshold
          for augmenting the vector subspace in stability check !expert-*/
        options.add_double("STABILITY_AUGMENT_SPACE_TOL", 0.1);
        /*- Controls the convergence of the Davidson's diagonalization in stability check !expert-*/
        options.add_double("STABILITY_CONVERGENCE", 1e-4);
        /*- The number of vectors that can be added simultaneously into the subspace for Davidson's diagonalization in stability check !expert-*/
        options.add_int("STABILITY_ADD_VECTORS", 20);
        /*- The number of guess vectors used for Davidson's diagonalization in stability check !expert-*/
        options.add_int("STABILITY_N_GUESS_VECTORS", 20);
        /*- The number of Hessian eigenvalues computed during the stability check !expert-*/
        options.add_int("STABILITY_N_EIGENVALUES", 3);
        /*- The maximum size of the subspace for the stability check. The program will terminate if this parameter is exceeded
          and the convergence (STABILITY_CONVERGENCE) is not satisfied !expert-*/
        options.add_int("STABILITY_MAX_SPACE_SIZE", 200);
        /*- Controls whether to relax tau during the cumulant updates or not !expert-*/
        options.add_bool("RELAX_TAU", true);
        /*- Chooses appropriate DCFT method -*/
        options.add_str("DCFT_FUNCTIONAL", "DC-06", "DC-06 DC-12 ODC-06 ODC-12 CEPA0");
    }

    return true;
}

PsiReturnType dcft(Options&);

extern "C"
PsiReturnType dfdcft(Options& options)
{
    return dcft(options);
}

}} // End namespaces

