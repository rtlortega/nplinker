# Copyright 2021 The NPLinker Authors
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

import os
import glob

from ..logconfig import LogConfig

logger = LogConfig.getLogger(__file__)


# load Chem_class_predictions (canopus, molnetenhancer are loaded)
# for canopus, check if results can be converted with canopus_treemap
# otherwise use the pre-existing output of canopus
class ChemClassPredictions:
    """Class for storing results for chemical class predictions of spectra

    Currently, CANOPUS and MolNetEnhancer results are loaded
    """
    def __init__(self, canopus_dir, mne_dir):
        """Load classes with CanopusResults, MolNetEnhancerResults

        Args:
            canopus_dir: str, canopus_dir found in root_dir of nplinker project
            mne_dir: str, mne_dir found in root_dir of nplinker project
        """
        # todo: use canopus_treemap to convert canopus result
        self._canopus = CanopusResults(canopus_dir)
        self._molnetenhancer = MolNetEnhancerResults(mne_dir)

        class_predict_options = []
        if self._canopus.spectra_classes:
            class_predict_options.append('canopus')
        if self._molnetenhancer.spectra2molfam:
            class_predict_options.append('molnetenhancer')
        if class_predict_options:
            class_predict_options = ['mix', 'main'] + class_predict_options
        self._class_predict_options = class_predict_options

    @property
    def canopus(self):
        return self._canopus

    @property
    def molnetenhancer(self):
        return self._molnetenhancer

    @property
    def class_predict_options(self):
        """The available class predictions"""
        return self._class_predict_options


class CanopusResults:
    """Class for storing canopus results

    The two input files from input_dir are read for the spectra and molfams, respectively:
        -cluster_index_classifications.txt
        -component_index_classifications.txt
    """

    def __init__(self, canopus_dir):
        """Read the class info from root_dir/canopus

        Args:
            canopus_dir: str, canopus_dir found in root_dir of nplinker project
        """
        spectra_classes_names, spectra_classes = self._read_spectra_classes(
            canopus_dir)
        self._spectra_classes = spectra_classes
        self._spectra_classes_names = spectra_classes_names
        self._spectra_classes_names_inds = {elem: i for i, elem in
                                            enumerate(spectra_classes_names)}

        molfam_classes_names, molfam_classes = self._read_molfam_classes(
            canopus_dir)
        self._molfam_classes = molfam_classes
        self._molfam_classes_names = molfam_classes_names
        self._molfam_classes_names_inds = {elem: i for i, elem in
                                           enumerate(molfam_classes_names)}

    def _read_spectra_classes(self, canopus_dir):
        """Read canopus classes for spectra, return classes_names, classes

        Args:
            root_dir: str, root_dir of nplinker project
        Returns:
            Tuple of:
            - ci_classes_names: list of str - the names of each different level
            - ci_classes: dict of {str, lists of tuple(str, float)} - per spectrum (key) the classes for each level
                where each level is a list of (class_name, score) sorted on best choice so index 0 is the best
                class prediction for a level
        """
        input_file = glob.glob(os.path.join(
            canopus_dir, '*cluster_index_classifications.txt'))[0]

        ci_classes = {}  # make a dict {ci: [[(class,score)]]}
        ci_classes_header = None
        ci_classes_names = []

        if os.path.isfile(input_file):
            logger.info(f"reading canopus results from {canopus_dir}")
            with open(input_file) as inf:
                ci_classes_header = inf.readline().strip().split("\t")
                for line in inf:
                    line = line.strip('\n').split("\t")
                    classes_list = []
                    for lvl in line[3:]:
                        lvl_list = []
                        for l_class in lvl.split("; "):
                            if l_class:
                                l_class = l_class.split(":")
                                c_tup = tuple([l_class[0], float(l_class[1])])
                            else:
                                c_tup = None  # default value for class value
                            lvl_list.append(c_tup)
                        classes_list.append(lvl_list)
                    ci_classes[line[1]] = classes_list
        else:
            logger.warn(
                'could not load cluster_index_classifications.txt; missing from canopus_dir')

        if ci_classes_header:
            #  todo: rename the output from the canopus script directly
            ci_classes_names = [f"cf_{elem}" for elem in
                                ci_classes_header[3:-3]] + \
                               [f"npc_{elem}" for elem in
                                ci_classes_header[-3:]]
        return ci_classes_names, ci_classes

    def _read_molfam_classes(self, canopus_dir):
        """Read canopus classes for molfams, return classes_names, classes

        Args:
            root_dir: str, root_dir of nplinker project
        Returns:
            Tuple of:
            - compi_classes_names: list of str - the names of each different level
            - compi_classes: dict of {str: lists of tuple(str, float)} - per molfam (key) the classes for each level
                where each level is a list of (class_name, fraction) sorted on best choice so index 0 is the best
                class prediction for a level
        """
        input_file = glob.glob(os.path.join(
            canopus_dir, '*component_index_classifications.txt'))[0]

        compi_classes = {}  # make a dict {compi: [[(class,score)]]}
        compi_classes_header = None
        compi_classes_names = []

        if os.path.isfile(input_file):
            with open(input_file) as inf:
                compi_classes_header = inf.readline().strip().split("\t")
                for line in inf:
                    line = line.strip('\n').split("\t")
                    classes_list = []
                    for lvl in line[2:]:
                        lvl_list = []
                        for l_class in lvl.split("; "):
                            if l_class:
                                l_class = l_class.split(":")
                                c_tup = tuple([l_class[0], float(l_class[1])])
                            else:
                                c_tup = None  # default value for class value
                            lvl_list.append(c_tup)
                        classes_list.append(lvl_list)
                    compi_classes[line[0]] = classes_list
        else:
            logger.warn(
                'could not load component_index_classifications.txt; missing from canopus_dir')

        if compi_classes_header:
            #  todo: rename the output from the canopus script directly
            compi_classes_names = [f"cf_{elem}" for elem in
                                   compi_classes_header[2:-3]] + \
                                  [f"npc_{elem}" for elem in
                                   compi_classes_header[-3:]]
        return compi_classes_names, compi_classes

    @property
    def spectra_classes(self):
        return self._spectra_classes

    @property
    def spectra_classes_names(self):
        return self._spectra_classes_names

    @property
    def spectra_classes_names_inds(self):
        return self._spectra_classes_names_inds

    @property
    def molfam_classes(self):
        return self._molfam_classes

    @property
    def molfam_classes_names(self):
        return self._molfam_classes_names

    @property
    def molfam_classes_names_inds(self):
        return self._molfam_classes_names_inds


class MolNetEnhancerResults:
    """Class for storing MolNetEnhancer results

    The input file for ClassyFire results is read from the molnetenhancer directory:
        - ClassyFireResults_Network.txt
    """

    def __init__(self, mne_dir):
        """Read the class info from file in root_dir/molnetenhancer/

        Args:
            mne_dir: str, mne_dir found in root_dir of nplinker project
        """
        cf_classes_names, molfam_classes, spectra2molfam = self._read_cf_classes(
            mne_dir)
        self._spectra2molfam = spectra2molfam
        self._molfam_classes = molfam_classes
        self._spectra_classes_names = cf_classes_names  # if NPC gets implemented, add here
        self._spectra_classes_names_inds = {elem: i for i, elem in
                                            enumerate(cf_classes_names)}

    def _read_cf_classes(self, mne_dir):
        """Read ClassyFireResults_Network.txt in molnetenhancer dir

        Args:
            mne_dir: str, mne_dir found in root_dir of nplinker project
        Returns:
            tuple of:
            -list of str - names of the classes in order
            -dict of {str: [(str, float)]} - linking molfams to (classes, scores) in order of names,
                singleton families are denoted with S[\d]+
            -dict of {str:str} - linking spectra to molfams
        """
        columns = []
        mne_component_dict = {}
        mne_cluster2component = {}
        # flexible finding of CF results
        input_file = 'not_found'
        wanted_file = "ClassyFireResults_Network.txt"
        possible_files = glob.glob(os.path.join(mne_dir, "*")) + \
                         glob.glob(os.path.join(mne_dir, "*", "*"))
        try:
            input_file = [pos_file for pos_file in possible_files
                          if pos_file.endswith(wanted_file)][0]
        except IndexError:
            pass

        if not os.path.isfile(input_file):
            logger.warn(f"no molnetenhancer input found in {mne_dir}")
            return columns, mne_component_dict, mne_cluster2component

        with open(input_file) as inf:
            logger.info(f"reading molnetenhancer results from {mne_dir}")
            header = inf.readline().strip().split("\t")
            # get the columns that are of interest to us
            columns = [
                'cf_direct_parent' if col == 'CF_Dparent' else col.lower()
                for i, col in enumerate(header[3:]) if i % 2 == 0]
            for line in inf:
                line = line.strip('\n').split("\t")
                cluster = line.pop(0)
                component = line.pop(0)
                nr_nodes = line.pop(0)
                # todo: make it easier to query classes of singleton families
                # if singleton family, format like '-1_spectrum-id' like canopus results
                if nr_nodes == '1':
                    component = f'-1_{cluster}'
                class_info = []
                # get class names and scores in (class, score)
                for i in range(0, len(line), 2):
                    class_tup = (line[i], float(line[i + 1]))
                    class_info.append(class_tup)
                if component not in mne_component_dict:
                    mne_component_dict[component] = class_info
                mne_cluster2component[cluster] = component

        return columns, mne_component_dict, mne_cluster2component

    def spectra_classes(self, spectrum_id):
        """Return classes by relating spectrum_id in the molfam_classes

        Args:
            spectrum_id: int/str, spectrum_id - ints will be converted to str
        """
        classes = []
        if isinstance(spectrum_id, int):
            spectrum_id = str(spectrum_id)
        molfam_id = self.spectra2molfam.get(spectrum_id)
        if molfam_id:
            classes = self.molfam_classes.get(molfam_id)
        return classes

    @property
    def spectra2molfam(self):
        return self._spectra2molfam

    @property
    def spectra_classes_names(self):
        return self._spectra_classes_names

    @property
    def spectra_classes_names_inds(self):
        return self._spectra_classes_names_inds

    @property
    def molfam_classes(self):
        return self._molfam_classes
