import os
from pytest import fixture
from nplinker.genomics import GCF
from nplinker.metabolomics import MolecularFamily
from nplinker.metabolomics import Spectrum
from nplinker.nplinker import NPLinker
from nplinker.scoring import MetcalfScoring
from nplinker.strain import Strain
from nplinker.strain import StrainCollection
from .. import CONFIG_FILE_LOCAL_MODE


@fixture(scope="session")
def strains_list() -> tuple[Strain, Strain, Strain]:
    return Strain("strain1"), Strain("strain2"), Strain("strain3")


@fixture(scope="session")
def strains(strains_list) -> StrainCollection:
    strains = StrainCollection()
    for strain in strains_list:
        strains.add(strain)
    return strains


@fixture(scope="session")
def gcfs(strains_list) -> tuple[GCF, GCF, GCF]:
    gcf1 = GCF("gcf1")
    gcf1.strains.add(strains_list[0])
    gcf2 = GCF("gcf2")
    gcf2.strains.add(strains_list[1])
    gcf3 = GCF("gcf3")
    gcf3.strains.add(strains_list[0])
    gcf3.strains.add(strains_list[1])
    return gcf1, gcf2, gcf3


@fixture(scope="session")
def spectra(strains_list) -> tuple[Spectrum, Spectrum, Spectrum]:
    spectrum1 = Spectrum("spectrum1", [1], [1], 10.0)
    spectrum1.strains.add(strains_list[0])
    spectrum2 = Spectrum("spectrum2", [1], [1], 10.0)
    spectrum2.strains.add(strains_list[1])
    spectrum3 = Spectrum("spectrum3", [1], [1], 10.0)
    spectrum3.strains.add(strains_list[0])
    spectrum3.strains.add(strains_list[1])
    return spectrum1, spectrum2, spectrum3


@fixture(scope="session")
def mfs(spectra) -> tuple[MolecularFamily, MolecularFamily, MolecularFamily]:
    """For simplicity, we just use one Spectrum object for each MolecularFamily
    object, and notice that they are not SingletonFamily object.
    """
    mf1 = MolecularFamily("mf1")
    mf1.add_spectrum(spectra[0])
    mf2 = MolecularFamily("mf2")
    mf2.add_spectrum(spectra[1])
    mf3 = MolecularFamily("mf3")
    mf3.add_spectrum(spectra[2])
    return mf1, mf2, mf3


@fixture(scope="function")
def npl(gcfs, spectra, mfs, strains, tmp_path) -> NPLinker:
    """Constructed NPLinker object.

    This NPLinker object does not do loading `npl.load_data()`, instead we
    manually set its attributes to the values we want to test.

    The config file `nplinker_demo1.toml` does not affect the tests, just
    making sure the NPLinker object can be created successfully.
    """
    os.environ["NPLINKER_ROOT_DIR"] = str(tmp_path)  # Create a temporary root dir for NPLinker
    npl = NPLinker(CONFIG_FILE_LOCAL_MODE)
    npl._strains = strains
    npl._gcf_dict = {gcf.id: gcf for gcf in gcfs}
    npl._mf_dict = {mf.id: mf for mf in mfs}
    npl._spec_dict = {spec.id: spec for spec in spectra}
    return npl


@fixture(scope="function")
def mc(npl) -> MetcalfScoring:
    """MetcalfScoring object."""
    mc = MetcalfScoring()
    mc.setup(npl)
    return mc
