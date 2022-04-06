"""
This is PISACov, a PISA extension to infer quaternary structure
of proteins from evolutionary covariance.
"""

from pisacov import __prog__, __description__, __version__
from pisacov import __author__, __date__, __copyright__

import copy
import logging

from crops.elements import sequences as pes
from crops.io import taggers as ctg
from conkit.core import contactmap as ckc

def backmapping(conpredmap, sequence):
    """
    Return the contact prediction map with the original residue numbers.

    :param conpredmap: Contact prediction map
    :type conpredmap: :class:`~conkit.core.contactmap.ContactMap`
    :param sequence: Sequence.
    :type sequence: :class:`~crops.elements.sequences.sequence`
    :return: Contact prediction map with backmapped ids.
    :rtype: :class:`~conkit.core.contactmap.ContactMap`

    """
    if isinstance(conpredmap, ckc.ContactMap) is False:
        logging.critical('First argument must be a Conkit ContactMap object.')
        raise TypeError
    if isinstance(sequence, pes.sequence) is False:
        logging.critical('Second argument must be a CROPS sequence object.')
        raise TypeError

    conpredout = conpredmap.deepcopy()
    for contact in conpredmap:
        c1 = sequence.cropbackmap(contact.res1_seq)
        c2 = sequence.cropbackmap(contact.res2_seq)
        nid = "(" + str(c1) + ", " + str(c2) + ")"
        contact.res1_seq = c1
        contact.res2_seq = c2
        contact.id = nid

    return conpredout


def remove_neighbours(conpredmap, mindist=2):
    """
    Return :class:`~conkit.core.contactmap.ContactMap` without neighbouring pairs.

    :param conpredmap: Contact map.
    :type conpredmap: :class:`~conkit.core.contactmap.ContactMap`
    :param mindist: Minimum allowed distance, defaults to 2.
    :type mindist: int, optional
    :return: Contact map.
    :rtype: :class:`~conkit.core.contactmap.ContactMap`

    """
    if isinstance(conpredmap, ckc.ContactMap) is False:
        logging.critical('First argument must be a Conkit ContactMap object.')
        raise TypeError

    return conpredmap.remove_neighbors(min_distance=mindist, inplace=True)

class contact_atlas:
    """A :class:`~pisacov.core.contacts.contact_atlas` object containing information from
    sequences, contact preduction maps and structure contacts.
    The :class:`~pisacov.core.contacts.contact_atlas` class represents a data structure to hold
    contact maps, matched and unmatched with structure contacts and sequence.

    :param seqid: Sequence identifier. Can be used alone or together with oligomer ID.
    :type seqid: str
    :param oligomer: Oligomer identifier. Sometimes as important as seqid.
    :type oligomer: str, optional
    :param seq: Sequence string.
    :type seq: str, optional

    :param header: Standard .fasta header, starting with ">".
    :type header: str, optional
    :ivar info: Useful information of the :class:`~crops.elements.sequence.monomer_sequence`.
    :vartype info: dict [str, any]
    :ivar seqs: The set of sequences, including default "mainseq", in :class:`~crops.elements.sequence.monomer_sequence`.
    :vartype seqs: dict [str, str]

    :example:

    >>> from crops.elements import sequences as csq
    >>> myseq = csq.sequence('exampleID')
    >>> mysq.mainseq('GATTACA')
    >>> myseq.mainseq()
    'GATTACA'
    >>> myseq.addseq('gapseq','GAT--C-')
    >>> myseq.addseq('cobra','TACATACA')
    >>> myseq.length()
    7
    >>> myseq.ngaps('gapseq')
    3
    >>> myseq.guess_biotype()
    'DNA'
    >>> print(myseq)
    Sequence object: ('>exampleID|Chain exampleID', seq='GATTACA', type='DNA', length=7)

    :example:

    >>> from crops.elements import sequences as csq
    >>> from crops.io import parsers as csp
    >>> myseq = csp.parseseqfile('7M6C.fasta')
    >>> myseq
    Sequence object: (>7M6C_1|Chain A, seq=MRTLWIMAVL[...]KPLCKKADPC, type=Undefined, length=138)
    >>> myseq.guess_biotype()
    'Protein'
    >>> myseq
    Sequence object: (>7M6C_1|Chain A, seq=MRTLWIMAVL[...]KPLCKKADPC, type=Protein, length=138)
    """

    _kind = 'Contact Atlas'
    __slots__ = ['name', 'chain1', 'chain2', 'sequence',
                 'conpred_raw', 'conpred', 'strmap_raw', 'strmap',
                 'conkitmatch_raw', 'conkitmatch', 'intramap',
                 'tp_raw', 'tn_raw', 'fp_raw', 'fn_raw',
                 'tp', 'tn', 'fp', 'fn', 'npotential']
    def __init__(self, name, conpredmap, strmap, sequence,
                 minneigh=2, removeintra=True, backmapping=False):
        self.name = name
        self.chain1 = None
        self.chain2 = None
        self.sequence = sequence
        self.conpred_raw = conpredmap
        self.conpred = None
        self.strmap_raw = strmap
        self.strmap = None
        self.conkitmatch = None
        self.tp = 0
        self.fp = 0
        self.tn = 0
        self.fn = 0
        self.npotential = 0
        # Set sequence values and Number of total potential contacts
        lseq = sequence.length()
        self.npotential = lseq**2 - lseq
        for n in range(1, minneigh):
            self.npotential -= 2*(lseq - n)
        self.npotential = int(self.potential / 2)
        # Set Contact prediction list
        self.conpred = remove_neighbours(conpredmap, mindist=minneigh)
        self.conpred.sort('raw_score', reverse=True, inplace=True)
        if backmapping is True:
            self.conpred = backmapping(self.conpred, sequence)
        self.conpred.seq = self.sequence.seqs['conkit']
        self.conpred.set_sequence_register()
        # Set Structure map
        self.strmap = strmap[1]
        self.chain1 = tuple(self.strmap.id)[0]
        self.chain2 = tuple(self.strmap.id)[1]
        # Remove intramolecular chains
        if removeintra is True:
            for m in [0,3]:
                intra = strmap[m]
                for contact1 in intra:
                    c1 = str(contact1.id)[1:-1].split(', ')
                    for contact2 in reversed(self.conpred):
                        c2 = str(contact2.id)[1:-1].split(', ')
                        if ((c1[0] == c2[0] and c1[1] == c2[1]) or
                                (c1[1] == c2[0] and c1[0] == c2[1])):
                            self.conpred.remove(c2)  # CHECK THAT REMOVAL INSIDE LOOP IS OK
        # MATCH
        self.conkitmatch = self.conpred.deepcopy().match(self.strmap, add_false_negatives=True, inplace=False)
        for contact in self.conkitmatch:
            if contact.true_positive:
                self.tp += 1
            elif contact.false_positive:
                self.fp += 1
            elif contact.false_negative:
                self.fn += 1
            elif contact.true_negative:
                logging.warning('True negatives appearing in conkit match.')
            else:
                logging.warning('Contact not evaluated.')

        self.tn = self.npotential - self.tp - self.fp - self.fn

    def __repr__(self):
        chtype = self.biotype if self.biotype is not None else 'Undefined'
        if 'mainseq' not in self.seqs:
            raise ValueError("'mainseq' sequence not found.")
        if len(self.seqs['mainseq']) <= 20:
            showseq = self.seqs['mainseq']
        else:
            showseq = (self.seqs['mainseq'][:10]+'[...]' +
                       self.seqs['mainseq'][len(self.seqs['mainseq'])-10:])
        tempolig = self.oligomer_id if self.oligomer_id is not None else 'NOID'
        shortid = ctg.makeheader(mainid=tempolig, seqid=self.name,
                                 chains=self.chains, short=True)
        string = (self._kind+" object "+shortid+" (seq="+str(showseq) +
                  ", type=" + chtype + ", length=" +
                  str(len(self.seqs['mainseq']))+")")
        return string

    def __iter__(self):
        return iter(self.seqs['mainseq'].values())

    def copy(self):
        return copy.copy(self)

    def deepcopy(self):
        return copy.deepcopy(self)
