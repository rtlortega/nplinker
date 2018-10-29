import csv,glob
import numpy as np

import aa_pred


class Strain(object):
    def __init__(self,name):
        self.name = name
    def __str__(self):
        return self.name

class BGC(object):
    def __init__(self,strain,name,bigscape_class,product_prediction):
        self.strain = strain
        self.name = name
        self.bigscape_class = bigscape_class
        self.product_prediction = product_prediction

        self.antismash_file = None
        self._aa_predictions = None

    def __str__(self):
        return self.name + "(" + str(self.strain) + ")"

    @property
    def aa_predictions(self):
        # Load aa predictions and cache them
        if self._aa_predictions is None:
            self._aa_predictions = []
            if self.antismash_file is not None:
                for p in aa_pred.predict_aa(self.antismash_file):
                    self._aa_predictions.append(p)
        return self._aa_predictions


class GCF(object):
    def __init__(self,gcf_id):
        self.gcf_id = gcf_id
        self.bgc_list = []
        self.random_gcf = None

        self._aa_predictions = None

    def add_bgc(self, bgc):
        self.bgc_list.append(bgc)

    def has_strain(self,strain):
        for bgc in self.bgc_list:
            if bgc.strain == strain:
                return True
        return False

    def add_random(self,strain_list):
        self.random_gcf = RandomGCF(self,strain_list)

    def get_mibig_bgcs(self):
        mibig = []
        for bgc in self.bgc_list:
            if type(bgc) == MiBIGBGC:
                mibig.append(bgc)
        return mibig

    @property
    def aa_predictions(self):
        """
        Return the predicted AAs for the GCF
        """
        if self._aa_predictions is None:
            # Make sure that we record a 0 probability if an AA is predicted
            # for _some_ but not _all_ BGCs in a GCF
            bgc_aa_prob = {}

            for bgc_count, bgc in enumerate(self.bgc_list):
                for aa, p_aa in bgc.aa_predictions:
                    # If we come across a new AA, set it to zero for all previous BGCs
                    if aa not in bgc_aa_prob:
                        aa_prob = [0.0]
                        bgc_aa_prob[aa] = aa_prob
                    bgc_aa_prob[aa].append(p_aa)

            for aa in list(bgc_aa_prob.keys()):
                # Replace the prob list with the mean
                bgc_aa_prob[aa] = np.mean(bgc_aa_prob[aa])

            self._aa_predictions = list(bgc_aa_prob.items())

        return self._aa_predictions


class RandomGCF(object):
    def __init__(self,real_gcf,strain_list):
        n_strains = 0
        self.real_gcf = real_gcf
        for s in strain_list:
            if self.real_gcf.has_strain(s):
                n_strains += 1
        # select n_strains from strain_list
        self.strain_set = set(np.random.choice(strain_list,n_strains,replace = False))
    def has_strain(self,strain):
        if strain in self.strain_set:
            return True
        else:
            return False


class MiBIGBGC(BGC):
    def __init__(self,name,product_prediction):
        super(MiBIGBGC,self).__init__(None,name,None,product_prediction)


def loadBGC_from_cluster_files(network_file_list,ann_file_list,antismash_dir = None,antismash_format = 'flat'):
    strain_id_dict = {}
    strain_dict = {}
    gcf_dict = {}
    gcf_list = []
    strain_list = []
    bgc_list = []
    with open('strain_ids.csv','r') as f:
        reader = csv.reader(f)
        for line in reader:
            strain_id_dict[line[0]] = line[1]
    metadata = {}
    for a in ann_file_list:
        with open(a,'rU') as f:
            reader =  csv.reader(f,delimiter = '\t')
            heads = next(reader)
            for line in reader:
                metadata[line[0]] = line

    for filename in network_file_list:
        with open(filename,'rU') as f:
            reader = csv.reader(f,delimiter = '\t')
            heads = next(reader)
            for line in reader:
                name = line[0]
                family = filename + ":" + line[1]
                if name.startswith("BGC"):
                    strain_name = 'MiBIG'
                else:
                    try:
                        try:
                            strain_name = strain_id_dict[name.split('_')[0]]
                        except:
                            strain_name = strain_id_dict[name.split('.')[0]]
                    except:
                        print("NO STRAIN")
                if not strain_name in strain_dict:
                    new_strain = Strain(strain_name)
                    strain_dict[strain_name] = new_strain
                    strain_list.append(new_strain)
                strain = strain_dict[strain_name]
                tokens = name.split('.')
                clusterid = tokens[-1]
                rest = '.'.join(tokens[:-1])

                metadata_line = metadata[name]
                description = metadata_line[2]
                bigscape_class = metadata_line[4]
                product_prediction = metadata_line[3]

                # make a BGC object
                # the same BGC objects might be made more than once
                # because they appear in multiple clusterings
                if not strain_name == 'MiBIG':
                    new_bgc = BGC(strain,name,bigscape_class,product_prediction)
                    if antismash_dir:
                        if antismash_format == 'flat':
                            new_bgc.antismash_file = find_antismash_file_flat(antismash_dir,new_bgc.name)
                        else:
                            new_bgc.antismash_file = find_antismash_file(antismash_dir,new_bgc.name)
                else:
                    new_bgc = MiBIGBGC(name,product_prediction)
                bgc_list.append(new_bgc)

                if not family in gcf_dict:
                    new_gcf = GCF(family)
                    gcf_dict[family] = new_gcf
                    gcf_list.append(new_gcf)
                gcf_dict[family].add_bgc(new_bgc)

    return gcf_list,bgc_list,strain_list

# this is really slow. But hey, it works. Finally.
def find_antismash_file_flat(antismash_dir,bgc_name):
    import glob,os
    all_gbk_files = glob.glob(antismash_dir + os.sep + '*.gbk')
    last_bit = [o.split(os.sep)[-1] for o in all_gbk_files]
    bgc_file_name = bgc_name + '.gbk'
    if bgc_file_name in last_bit:
        idx = last_bit.index(bgc_file_name)
        return all_gbk_files[idx]
    else:
        print("NOOO",bgc_name)
        return None


def find_antismash_file(antismash_dir,bgc_name):
    import glob,os
    subdirs = [s.split(os.sep)[-1] for s in glob.glob(antismash_dir + os.sep+'*')]
    if bgc_name.startswith('BGC'):
        print("No file for MiBIG BGC")
        return None # MiBIG BGC
    # this code is nasty... :-)
    name_tokens = bgc_name.split('_')
    found = False
    for i in range(len(name_tokens)):
        sub_name = '_'.join(name_tokens[:i])
        if sub_name in subdirs:
            found = True
            found_name = sub_name
    if not found:
        name_tokens = bgc_name.split('.')[0]
        for i in range(len(name_tokens)):
            sub_name = '.'.join(name_tokens[:i])
            if sub_name in subdirs:
                found = True
                found_name = sub_name
    if not found:
        print("Can't find antiSMASH info for ",bgc_name)
        return None
#     print found_name
    dir_contents = glob.glob(antismash_dir + os.sep + found_name + os.sep + '*.gbk')
    cluster_names = [d.split('.')[-2] for d in dir_contents]
    this_name = bgc_name.split('.')[-1]
    try:
        antismash_name = dir_contents[cluster_names.index(bgc_name.split('.')[-1])]
    except:
        print(bgc_name)
        print(cluster_names)
        print()
        print()
        return None
    return antismash_name

def loadBGC_from_node_files(file_list):
    strain_id_dict = {}
    with open('strain_ids.csv','r') as f:
        reader = csv.reader(f)
        for line in reader:
            strain_id_dict[line[0]] = line[1]
    
    strain_dict = {}
    gcf_dict = {}
    bgc_list = []
    gcf_list = []
    strain_list = []
    for filename in file_list:
        with open(filename,'rU') as f:
            reader = csv.reader(f)
            heads = next(reader)
            name_pos = heads.index("shared name")
            description_pos = heads.index("Description")
            bigscape_class_pos = heads.index("BiG-SCAPE class")
            product_prediction_pos = heads.index("Product Prediction")
            family_pos = heads.index("Family Number")
            for line in reader:
                name = line[name_pos]
                try:
                    try:
                        strain_name = strain_id_dict[name.split('_')[0]]
                    except:
                        strain_name = strain_id_dict[name.split('.')[0]]
                except:
                    # it's a MiBIG one
                    strain_name = 'MiBIG'
                if not strain_name in strain_dict:
                    new_strain = Strain(strain_name)
                    strain_dict[strain_name] = new_strain
                    strain_list.append(new_strain)
                strain = strain_dict[strain_name]


                tokens = name.split('.')
                clusterid = tokens[-1]
                rest = '.'.join(tokens[:-1])
                # print name,rest,clusterid
                description = line[description_pos]
                bigscape_class = line[bigscape_class_pos]
                product_prediction = line[product_prediction_pos]
                family = filename + " " + line[family_pos]

                # make a BGC object
                if not strain_name == 'MiBIG':
                    new_bgc = BGC(strain,name,bigscape_class,product_prediction)
                else:
                    new_bgc = MiBIGBGC(name,product_prediction)
                bgc_list.append(new_bgc)

                if not family in gcf_dict:
                    new_gcf = GCF(family)
                    gcf_dict[family] = new_gcf
                    gcf_list.append(new_gcf)
                gcf_dict[family].add_bgc(new_bgc)


                
    return gcf_list,bgc_list,strain_list
    

def load_mibig_map(filename = 'mibig_gnps_links_q3_loose.csv'):
    mibig_map = {}
    with open(filename,'rU') as f:
        reader = csv.reader(f)
        heads = next(reader)
        
        for line in reader:
            bgc = line[0]
            if bgc in mibig_map:
                mibig_map[bgc].append(line[3])
            else:
                mibig_map[bgc] = [line[3]]
    return mibig_map


def load_mibig_library_json(mibig_json_directory):
    import glob,os,json
    mibig = {}
    files = glob.glob(mibig_json_directory + os.sep + '*.json')
    print("Found {} files".format(len(files)))
    for file in files:
        with open(file,'r') as f:
            bgc_id = file.split(os.sep)[-1].split('.')[0]
            mibig[bgc_id] = json.load(f)
    return mibig


