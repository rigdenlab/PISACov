"""
This is PISACov, a PISA extension to infer quaternary structure
of proteins from evolutionary covariance.
"""

from pisacov import __prog__, __description__, __version__
from pisacov import __author__, __date__, __copyright__

import copy
import logging

import xml.etree.ElementTree as ET

def parse_interface_xml(interface_xml_path, assembly_xml_path = None):
    """Read interface.xml pisa file, parse it, and return number of interfaces.

    :param xml_path: Path to xml file.
    :type xml_path: str
    :return: List of interface info.
    :rtype: list [:class:`~pisacov.core.interfaces.interface`]

    """
    try:
        xmlparse = ET.parse(interface_xml_path)
    except Exception:
        logging.critical("ERROR: Unable to open XML file at " +
                         interface_xml_path)
        raise OSError

    root = xmlparse.getroot()
    ifinfolist = []
    for iface in root.iter('interface'):
        ifid = iface.find('id').text
        ifinfolist.append(interface(name = ifid))
        for mol in iface.iter('molecule'):
            cid = mol.find('chain_id').text
            ctype = mol.find('class').text
            ifinfolist[-1].chains[cid] = ctype

    if assembly_xml_path is not None:
        try:
            xmlparse = ET.parse(assembly_xml_path)
        except Exception:
            logging.critical("ERROR: Unable to open XML file at " +
                             assembly_xml_path)
            raise OSError
        root = xmlparse.getroot()
        for asm_set in root.iter('asm_set'):
            for assembly in asm_set.iter('assembly'):
                for ifaces in assembly.iter('interfaces'):
                    for iface in ifaces.iter('interface'):
                        iid = iface.find('id')
                        diss = iface.find('dissociates')
                        if ifinfolist[int(iid)-1].name == ifinfolist[iid]:
                            if diss == 'No':
                                ifinfolist[int(iid)-1].stable = True
                            elif diss == 'Yes':
                                ifinfolist[int(iid)-1].stable = False
                        else:
                            for iff in ifinfolist:
                                if iff.name == iid:
                                    if diss == 'No':
                                        iff.stable = True
                                    elif diss == 'Yes':
                                        iff.stable = False

    return ifinfolist

class interface:
    """A :class:`~pisacov.core.interfaces.interface` object containing information
    like composition and stability.

    :param name: Name of the interface.
    :type seqid: str
    :param structure: Save parsed structure here, defaults to None.
    :type oligomer: Any structure format.
    :param chains: A dictionary containing the chain ids as keys and the type as values.
    :type seq: dict [str: str]

    :ivar stable: Stability of the interface
    :vartype stable: bool, str
    """
    _kind = 'Interface'
    __slots__ = ['name', 'chains', 'stable', 'structure']

    def __init__(self, name, structure=None, chains=None):
        self.name = name
        self.chains = {}
        self.structure = None
        self.stable = False

    def __repr__(self):
        chainstring = ''
        typestring = ''
        for key, value in self.chains.items():
            chainstring += key + ','
            typestring += value +'-'
        chainstring = chainstring[:-1]
        typestring = typestring[:-1]
        string = (self._kind + " object " + self.name +
                  " (chains="+chainstring + ", type=" + typestring +
                  ", stable=" + str(self.stable) + ")")
        return string

    def __iter__(self):
        return iter(self.seqs['mainseq'].values())

    def copy(self):
        return copy.copy(self)

    def deepcopy(self):
        return copy.deepcopy(self)
