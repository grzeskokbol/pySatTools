import importlib.machinery
import os, sys
from astropy.utils.iers import iers

def _load_module(path):
    """
     Helper to load a Python file at path and return as a module
    :param path: path to a file
    :return: module
    """

    module_name = os.path.splitext(os.path.basename(path))[0]

    module = None
    if sys.version_info.minor < 5:
        loader = importlib.machinery.SourceFileLoader(module_name, path)
        module = loader.load_module()
    else:
        spec = importlib.util.spec_from_file_location(module_name, path)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)

    return module

def GetIersTable():
    """

    :return: iers table
    """

    # we need the IERS values to correct for the precession/nutation of the Earth
    try:
        iers_b_file = iers.download_file(iers.IERS_B_URL, cache=True)
        iers_b = iers.IERS_B.open(iers_b_file)
        iers.IERS.iers_table = iers_b
    except:
        try:
            iers_a_file = iers.download_file(iers.IERS_A_URL, cache=True)
            iers_a = iers.IERS_A.open(iers_a_file)
            iers.IERS.iers_table=iers_a
        except:
            warnings.warn('Could not load IERS_B nor IERS_A !')

    return iers.IERS.iers_table