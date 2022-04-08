"""
This is PISACov, a PISA extension to infer quaternary structure
of proteins from evolutionary covariance.
"""

from pisacov import __prog__, __description__, __version__
from pisacov import __author__, __date__, __copyright__

import copy
import logging

import xml.etree.ElementTree as ET

def n_int_xml(xml_path):
    """Read interface.xml pisa file, parse it, and return number of interfaces.

    :param xml_path: Path to xml file.
    :type xml_path: str
    :return: Number of interfaces.
    :rtype: int

    """
    try:
        xmlparse = ET.parse(xml_path)
    except Exception:
        logging.critical("ERROR: Unable to open XML file at " + xml_path)
        raise OSError

    root = xmlparse.getroot()
    n = int(root.find('n_interfaces').text)

    return n

def parse_interface_xml(xml_path):
    """Read interface.xml pisa file, parse it, and return number of interfaces.

    :param xml_path: Path to xml file.
    :type xml_path: str
    :return: Number of interfaces.
    :rtype: int

    """
    try:
        xmlparse = ET.parse(xml_path)
    except Exception:
        logging.critical("ERROR: Unable to open XML file at " + xml_path)
        raise OSError

    root = xmlparse.getroot()
    n = int(root.find('n_interfaces').text)

    return n

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

    def get_stability(self, xmlpath):
        # PARSE assembly.xml and set the stability values
