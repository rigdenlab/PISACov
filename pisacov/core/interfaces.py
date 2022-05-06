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
            pid = mol.find('id').text
            if pid == '1':
                did = 'A'
            elif pid == '2':
                did = 'B'
            newchain=chain_info(pisa_id=pid,
                                monomer_id=cid,
                                biotype=ctype,
                                dimer_id=did)
            ifinfolist[-1].chains.append(newchain)

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
                        iid = iface.find('id').text
                        diss = iface.find('dissociates').text
                        for ifobject in ifinfolist:
                            if ifobject.name == iid:
                                if diss == 'No':
                                    if ifinfolist[int(iid)-1].stable is True:
                                        pass
                                    else:
                                        ifinfolist[int(iid)-1].stable = True
                                elif diss == 'Yes':
                                    pass

    return ifinfolist

class chain_info:
    _kind = 'Chain Info'
    __slots__ = ['pisa_id', 'dimer_id', 'monomer_id', 'seq_id', 'type']

    def __init__(self, dimer_id=None, pisa_id=None, monomer_id=None,
                 seqid=None, biotype=None):
        self.dimer_id = dimer_id # ID IN INPUT DIMER (not equal to monomer_id)
        self.pisa_id = pisa_id # ID ASSIGNED BY PISA (1, 2)
        self.monomer_id = monomer_id # ID IN ORIGINAL STRUCTURE AND SEQUENCE
        self.seq_id = seqid
        self.type = biotype

    def __repr__(self):

        string = (self._kind + " object: (Dimer ID = " + str(self.dimer_id))
        if self.pisa_id is not None:
            string += ", " + "PISA xml ID = " + str(self.pisa_id)
        if self.monomer_id is not None:
            ", Monomer ID in assymetric unit = " + str(self.monomer_id)
        if self.seq_id is not None:
            ", Sequence ID = " + str(self.seq_id)
        if self.type is not None:
            ", Biotype = " + str(self.type)
        string += ")"

        return string

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
    __slots__ = ['name', 'chains', 'stable', 'structure', 'contactmap']

    def __init__(self, name, structure=None, chains=None):
        self.name = name
        self.chains = []
        self.structure = None
        self.contactmap = None
        self.stable = False

    def __repr__(self):
        chainstring = ''
        typestring = ''
        for chain in self.chains:
            chainstring += chain.monomer_id + ','
            typestring += chain.type +'-'
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
