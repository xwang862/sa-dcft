import psi4
import re
import os
import inputparser
import math
import warnings
from driver import *
from wrappers import *
from molutil import *
import p4util
from p4xcpt import *
#from psiexceptions import *


def run_dfdcft(name, **kwargs):
    r"""Function encoding sequence of PSI module and plugin calls so that
    dfdcft can be called via :py:func:`~driver.energy`. For post-scf plugins.

    >>> energy('dfdcft')

    """
    lowername = name.lower()
    kwargs = p4util.kwargs_lower(kwargs)

#    psi4.set_local_option('SCF', 'REFERENCE', 'UHF')
#    psi4.set_local_option('DFDCFT', 'REFERENCE', 'UHF')

    # Your plugin's psi4 run sequence goes here
    scf_helper(name, **kwargs)
    returnvalue = psi4.plugin('dfdcft.so')
    psi4.set_variable('CURRENT ENERGY', returnvalue)

def run_dfdcft_gradient(name, **kwargs):
    """Function encoding sequence of PSI module calls for
    DCFT gradient calculation.

    """
    optstash = p4util.OptionsState(
        ['GLOBALS', 'DERTYPE'])

    psi4.set_global_option('DERTYPE', 'FIRST')
    run_dfdcft(name, **kwargs)
    psi4.deriv()

    optstash.restore()


# Integration with driver routines
procedures['energy']['dfdcft'] = run_dfdcft
procedures['gradient']['dfdcft'] = run_dfdcft_gradient


def exampleFN():
    # Your Python code goes here
    pass
