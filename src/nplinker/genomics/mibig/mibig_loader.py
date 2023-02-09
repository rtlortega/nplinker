import os.path
from nplinker.logconfig import LogConfig
from nplinker.utils import list_files
from nplinker.strains import Strain
from .mibig_bgc import MibigBGC
from .mibig_metadata import MibigMetadata
from ..abc import BGCLoaderBase

logger = LogConfig.getLogger(__name__)


class MibigBGCLoader:
    def __init__(self, data_dir: str):
        """Parse MIBiG metadata files and return MibigBGC objects

        MIBiG metadata file (json) contains annotations/metadata information
        for each BGC. See https://mibig.secondarymetabolites.org/download.

        Args:
            data_dir(str): Path to the directory of MIBiG metadata json files
        """
        self.data_dir = data_dir
        self._file_dict = self.parse_data_dir(self.data_dir)
        self._metadata_dict = self._parse_metadatas()
        self._bgc_dict = self._parse_bgcs()

    def get_files(self) -> dict[str, str]:
        """Get the path of all MIBiG metadata json files.

        Returns:
            dict[str, str]: key is metadata file name (BGC accession), value is
                path to the metadata json file
        """
        return self._file_dict

    @staticmethod
    def parse_data_dir(data_dir: str) -> dict[str, str]:
        """Parse metadata directory and return pathes to all metadata json
            files.

        Args:
            data_dir(str): path to the directory of MIBiG metadata json files

        Returns:
            dict[str, str]: key is metadata file name (BGC accession), value is
                 path to the metadata json file
        """
        file_dict = {}
        json_files = list_files(data_dir, prefix="BGC", suffix=".json")
        for file in json_files:
            fname = os.path.splitext(os.path.basename(file))[0]
            file_dict[fname] = file
        return file_dict

    def get_metadatas(self) -> dict[str, MibigMetadata]:
        """Get MibigMetadata objects.

        Returns:
            dict[str, MibigMetadata]: key is BGC accession (file name) and
                value is :class:`nplinker.genomics.mibig.MibigMetadata` object
        """
        return self._metadata_dict

    def _parse_metadatas(self) -> dict[str, MibigMetadata]:
        """Parse all metadata files and return MibigMetadata objects.

        Returns:
            dict[str, MibigMetadata]: key is BGC accession (file name) and
                value is :class:`nplinker.genomics.mibig.MibigMetadata` object
        """
        metadata_dict = {}
        for name, file in self._file_dict.items():
            metadata = MibigMetadata(file)
            metadata_dict[name] = metadata
        return metadata_dict

    def get_bgcs(self) -> dict[str, MibigBGC]:
        """Get MibigBGC objects.

        Returns:
            dict[str, MibigBGC]: key is BGC name and value is
                :class:`nplinker.genomics.mibig.MibigBGC` object
        """
        return self._bgc_dict

    def _parse_bgcs(self) -> dict[str, MibigBGC]:
        """Parse all metadata files as MibigBGC objects

        Returns:
            dict[str, MibigBGC]: key is BGC accession (file name) and value is
                MibigBGC object
        """
        bgc_dict = {}
        i = 0
        for name, file in self._file_dict.items():
            bgc = parse_bgc_metadata_json(file)
            bgc.id = i
            bgc_dict[name] = bgc
            i += 1
        return bgc_dict


def parse_bgc_metadata_json(file: str) -> MibigBGC:
    """Parse MIBiG metadata file and return MibigBGC object

    Note:
        Index of MibigBGC object (`MibigBGC.id`) is set to `-1`.

    Args:
        file(str): Path to the MIBiG metadata json file

    Returns:
        MibigBGC: :class:`nplinker.genomics.mibig.MibigBGC` object
    """
    metadata = MibigMetadata(file)
    strain = Strain(metadata.mibig_accession)
    mibig_bgc = MibigBGC(-1, strain, metadata.mibig_accession,
                            metadata.biosyn_class)
    return mibig_bgc


# register as virtual class to prevent metaclass conflicts
BGCLoaderBase.register(MibigBGCLoader)