"""
This is PISACov, a PISA extension to infer quaternary structure
of proteins from evolutionary covariance.
"""

from pisacov import __prog__, __description__, __version__
from pisacov import __author__, __date__, __copyright__

import copy

def backmapping(conpredmap, sequence):
    conpredout = conpredmap.deepcopy()
    for contact in conpredmap:
        c1 = sequence.cropbackmap(contact.res1_seq)
        c2 = sequence.cropbackmap(contact.res2_seq)
        nid = "(" + str(c1) + ", " + str(c2) + ")"
        contact.res1_seq = c1
        contact.res2_seq = c2
        contact.id = nid

    return conpredout


def match_maps(conpredmap, strmap,  )

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
                 'tp', 'tn', 'fp', 'fn']
    def __init__(self, name=None, conpredmap=None, strmap=None, sequence=None,
                 minneigh=2, removeintra=True):
        self.name = None
        self.chain1 = None
        self.chain2 = None
        self.sequence = None
        self.conpred_raw = None
        self.conpred = None
        self.strmap_raw = {}
        self.strmap = None
        self.intramap = None
        self.conkitmatch = None
        self.tp_raw = None
        self.tp = None
        self.fp_raw = None
        self.fp = None
        self.tn_raw = None
        self.tn = None
        self.fn_raw = None
        self.fn = None

        if name is not None:
            self.name = name

        if conpredmap is not None:
            self.conpred_raw = conpredmap[1]

        if strmap is not None:
            self.strmap_raw = strmap

        if sequence is not None:
            self.sequence = sequence

        if conpredmap is not None and strmap is not None:
            matchedmap = conpredInt_pdb.match(intmap[nif][1], add_false_negatives=False, inplace=False)
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
        shortid = makeheader(mainid=tempolig, seqid=self.name,
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
