"""
This is PISACov, a program designed to infer quaternary structure
of proteins from evolutionary covariance.
"""

from pisacov import __prog__, __description__, __version__
from pisacov import __author__, __date__, __copyright__

import copy
import logging

from pisacov import iomod as pio
import xml.etree.ElementTree as ET

def parse_interface_xml(interface_xml_path, assembly_xml_path = None):
    """Read interface.xml pisa file, parse it, and return interface list.

    :param xml_path: Path to xml file.
    :type xml_path: str
    :return: List of interface info.
    :rtype: list [:class:`~conkit.core.contactmap.ContactMap`]

    """
    try:
        # xmlparse = ET.parse(interface_xml_path)
        xmlparse = pio.read(interface_xml_path, 'xml')
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
                                crystal_id=cid,
                                biotype=ctype,
                                dimer_id=did)
            ifinfolist[-1].chains.append(newchain)

    if assembly_xml_path is not None:
        try:
            # xmlparse = ET.parse(assembly_xml_path)
            xmlparse = pio.read(assembly_xml_path, 'xml')
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
                                    if ifinfolist[int(iid)-1].stable is True:
                                        pass
                                    else:
                                        ifinfolist[int(iid)-1].stable = False

    return ifinfolist

class chain_info:
    """A :class:`~pisacov.core.interfaces.chain_info` object contains several identifiers for a chain.

    :param dimer_id: ID for chain in PISA's interface PDB structure, defaults to None.
    :type dimer_id: str, optional
    :param pisa_id: ID assigned by PISA in xml file (1, 2), defaults to None.
    :type pisa_id: str, optional
    :param crystal_id: ID for chain in original PDB structure, defaults to None.
    :type crystal_id: str, optional
    :param seqid: Sequence ID assigned by CROPS (e.g. Seq ID = 1 for chains A, B, C), defaults to None.
    :type seqid: str, optional
    :param biotype: Biological classification of chain, defaults to None.
    :type biotype: str, optional

    :ivar dimer_id: ID for chain in PISA's interface PDB structure.
    :vartype dimer_id: str
    :ivar pisa_id: ID assigned by PISA in xml file (1, 2).
    :vartype pisa_id: str
    :ivar crystal_id: ID for chain in original PDB structure.
    :vartype crystal_id: str
    :ivar seq_id: Sequence ID assigned by CROPS (e.g. Seq ID = 1 for chains A, B, C).
    :vartype seqid: str
    :ivar type: Biological classification of chain.
    :vartype type: str

    :example:

    >>> from pisacov.core import interfaces as pci
    >>> mychain = pci.chain_info(dimer_id='B', pisa_id='1', crystal_id='A', seqid='1', biotype='Protein')
    >>> mychain
    Chain Info object: (Dimer ID = B, PISA xml ID = 1, Monomer ID in assymmetric unit = A, Sequence ID = 1, Biotype = Protein)
    """
    _kind = 'Chain Info'
    __slots__ = ['pisa_id', 'dimer_id', 'crystal_id', 'seq_id', 'type']

    def __init__(self, dimer_id=None, pisa_id=None, crystal_id=None,
                 seqid=None, biotype=None):
        self.dimer_id = dimer_id
        self.pisa_id = pisa_id
        self.crystal_id = crystal_id
        self.seq_id = seqid
        self.type = biotype

    def __repr__(self):

        string = (self._kind + " object: (Dimer ID = " + str(self.dimer_id))
        if self.pisa_id is not None:
            string += ", " + "PISA xml ID = " + str(self.pisa_id)
        if self.crystal_id is not None:
            string += ", Monomer ID in assymmetric unit = " + str(self.crystal_id)
        if self.seq_id is not None:
            string += ", Sequence ID = " + str(self.seq_id)
        if self.type is not None:
            string += ", Biotype = " + str(self.type)
        string += ")"

        return string

class interface:
    """A :class:`~pisacov.core.interfaces.interface` object containing information
    like composition and stability.

    :param name: Name of the interface.
    :type name: str
    :param structure: Save parsed Structure Contact maps here, defaults to None.
    :type structure: list [:class:`~conkit.core.contactmap.ContactMap`] or :class:`~conkit.core.contactfile.ContactFile`, optional
    :param chains: A list containing the two chain objects making up the interface, defaults to [].
    :type seq: list [:class:`~pisacov.core.interfaces.chain_info`], optional

    :ivar name: Name of the interface.
    :vartype name: str
    :ivar structure: Save parsed Structure Contact maps here.
    :vartype structure: list [:class:`~conkit.core.contactmap.ContactMap`] or :class:`~conkit.core.contactfile.ContactFile`
    :ivar chains: A list containing the two chain objects making up the interface.
    :vartype seq: list [:class:`~pisacov.core.interfaces.chain_info`]
    :ivar contactmap: Numpy-read contact list (generated from pdb structure).
    :vartype contactmap: :class:`~numpy.ndarray`
    :ivar stable: Stability of the interface
    :vartype stable: bool, str

    :example:

    >>> from pisacov.core import interfaces as pci
    >>> from conkit import io as ckio
    >>> mycontacts = ckio.read('mypath/mystructure.pdb', 'pdb')
    >>> mymaplist = []
    >>> for element in mycontacts
            mymaplist.append(element.as_contactmap())
    >>> mychains = []
    >>> for n in len(mycontacts[1].id):
            did = mycontacts[1].id[n]
            mychains.append(pci.chain_info(dimer_id=did, seqid='1', biotype='Protein')
    >>> myIF = pci.interface('My interface',
                             structure=mymaplist[1],
                             chains=mychains)
    >>> myIF
    Interface object My interface (chains=AB(in dimer), type=Protein-Protein, stable=Unknown)

    """
    _kind = 'Interface'
    __slots__ = ['name', 'chains', 'stable', 'structure', 'contactmap']

    def __init__(self, name, structure=None, chains=None):
        self.name = name
        if chains is None:
            self.chains = []
        else:
            self.chains = chains
        self.structure = None
        self.contactmap = None
        self.stable = None

    def __repr__(self):
        chainstring = ''
        chuse = 'crystal'
        for chain in self.chains:
            if chain.crystal_id is None:
                chuse = 'dimer'
                for chain2 in self.chains:
                    if chain2.dimer_id is None:
                        chuse = 'Unknown'
                        break
                break

        if chuse == 'Unknown':
            chainstring = chuse
        else:
            for chain in self.chains:
                if chuse == 'crystal':
                    chainstring += chain.crystal_id
                elif chuse == 'dimer':
                    chainstring += chain.dimer_id
            if chuse == 'crystal':
                chainstring += '(in crystal)'
            elif chuse == 'dimer':
                chainstring += '(in dimer)'

        typestring = ''
        if len(self.chains) < 2:
            typestring = 'None'
        elif len(self.chains) == 2:
            typestring = ''
            for chain in self.chains:
                typestring += chain.type + '-'
            typestring = typestring[:-1]
        else:
            logging.error('Interface contains more than 2 chains.')
            raise ValueError

        if self.stable is None:
            stability = 'Unknown'
        else:
            stability = str(self.stable)

        string = (self._kind + " object " + self.name +
                  " (chains="+chainstring + ", type=" + typestring +
                  ", stable=" + stability + ")")
        return string

    def __iter__(self):
        return iter(self.seqs['mainseq'].values())

    def copy(self):
        return copy.copy(self)

    def deepcopy(self):
        return copy.deepcopy(self)

    def itype(self):
        """Produces interface type string from information in self.chains.

        :raises ValueError: Interface contains more than 2 chains.
        :return: Interface type (only if 2 chains are defined.)
        :rtype: str

        """
        if len(self.chains) < 2:
            typestring = 'None'
        elif len(self.chains) == 2:
            typestring = ''
            for chain in self.chains:
                typestring += chain.type + '-'
            typestring = typestring[:-1]
        else:
            logging.error('Interface contains more than 2 chains.')
            raise ValueError
