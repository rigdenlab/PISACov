"""
This is PISACov, a PISA extension to infer quaternary structure
of proteins from evolutionary covariance.
"""

from pisacov import __prog__, __description__, __version__
from pisacov import __author__, __date__, __copyright__

from pisacov.io import _conf_ops as pco
from pisacov.core import _psicov_modes as PSICOV_modes

import copy
import logging
import os

from crops.elements import sequences as pes
from crops.io import taggers as ctg
from conkit.core import contactmap as ckc
from conkit.core.contactmap import ContactMap
from conkit import plot as ckplot

def backmapping(cmap, sequence):
    """
    Return the contact prediction map with the original residue numbers.

    :param cmap: Contact prediction map
    :type cmap: :class:`~conkit.core.contactmap.ContactMap`
    :param sequence: Sequence.
    :type sequence: :class:`~crops.elements.sequences.sequence`
    :return: Contact prediction map with backmapped ids.
    :rtype: :class:`~conkit.core.contactmap.ContactMap`

    """
    if ((isinstance(cmap, ckc.ContactMap) is False) and
            (isinstance(cmap, ContactMap) is False)):
        logging.critical('First argument must be a Conkit ContactMap object.')
        raise TypeError
    if isinstance(sequence, pes.sequence) is False:
        logging.critical('Second argument must be a CROPS sequence object.')
        raise TypeError

    conpredout = cmap.deepcopy()
    for n in range(len(cmap)):
        c1 = sequence.cropbackmap[cmap[n].res1_seq]
        c2 = sequence.cropbackmap[cmap[n].res2_seq]
        if c1 < c2:
            conpredout[n].res1_seq = c1
            conpredout[n].res2_seq = c2
            nid = c1, c2
        else:
            conpredout[n].res2_seq = c1
            conpredout[n].res1_seq = c2
            nid = c2, c1
        conpredout[n].id = nid

    return conpredout

def filter_contacts(cmap, threshold=0.2):
    """
    Remove low score contacts from contact prediction.

    :param cmap: Contact prediction map
    :type cmap: :class:`~conkit.core.contactmap.ContactMap`
    :param threshold: Threshold, defaults to 0.2
    :type threshold: float, optional

    """
    cmap.sort('raw_score', reverse=True, inplace=True)
    cnt = 0
    for contact in cmap:
        if contact.raw_score < threshold:
            break
        else:
            cnt=cnt+1

    return cmap[:cnt-1]


class contact_atlas:
    """A :class:`~pisacov.core.contacts.contact_atlas` object containing information from
    sequences, contact prediction maps and structure contacts.
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
    __slots__ = ['name', 'interface', 'sequence',
                 'conpred_raw', 'conpred', 'conpred_source',
                 'conkitmatch', 'ckplotmatch',
                 'tp', 'tn', 'fp', 'fn', 'npotential']
    def __init__(self, name, sequence, dimer_interface, conpredmap, conpredtype):
        self.name = name
        self.interface = dimer_interface
        self.conpred_raw = conpredmap
        self.conpred = self.conpred_raw.deepcopy()
        self.conpred_source = conpredtype
        self.sequence = sequence
        self.conkitmatch = {}
        self.ckplotmatch = {}
        self.tp = {}
        self.fp = {}
        self.tn = {}
        self.fn = {}
        self.npotential = None
        # Set sequence values and Number of total potential contacts
        lseq = sequence.length()
        self.npotential = (lseq**2 - lseq) / 2
        self.npotential = int(self.npotential / 2)

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

    def remove_neighbours(self, mindist=2):
        """
        Return :class:`~conkit.core.contactmap.ContactMap` without neighbouring pairs.

        :param mindist: Minimum allowed distance, defaults to 2.
        :type mindist: int, optional

        """
        lseq = self.sequence.length()
        self.npotential = lseq**2 - lseq
        for n in range(1, mindist):
            self.npotential -= 2*(lseq - n)
        self.npotential = int(self.npotential / 2)

        self.conpred.remove_neighbors(min_distance=mindist, inplace=True)

    def remove_intra(self):
        """
        Remove intramolecular contacts from intermolecular :class:`~conkit.core.contactmap.ContactMap`.

        """
        for m in [0,3]:
            intra = self.interface.structure[m]
            for contact1 in intra:
                c1 = str(contact1.id)[1:-1].split(', ')
                for contact2 in reversed(self.conpred):
                    c2 = str(contact2.id)[1:-1].split(', ')
                    if ((c1[0] == c2[0] and c1[1] == c2[1]) or
                            (c1[1] == c2[0] and c1[0] == c2[1])):
                        self.conpred.remove(contact2.id)  # CHECK THAT REMOVAL INSIDE LOOP IS OK

    def set_cropmap(self):
        """
        Renumber :class:`~conkit.core.contactmap.ContactMap` according to cropped sequence.

        """
        self.conpred = backmapping(self.conpred, self.sequence)


    def set_conpred_seq(self, sequence=None):
        """
        Renumber :class:`~conkit.core.contactmap.ContactMap` according to cropped sequence.

        :param sequence: Conkit type sequence, defaults to self.sequence.seqs['conkit'].
        :type mindist: :class:`~conkit.core.sequence.Sequence`

        """
        seq_in = self.sequence.seqs['conkit'] if sequence is None else sequence
        self.conpred.sequence = seq_in
        self.conpred.set_sequence_register()


    def make_match(self, filterout=None):
        """
        Match Structure and contact prediction :class:`~conkit.core.contactmap.ContactMap`.

        :param filterout: Threshold, defaults to None
        :type filterout: float, optional

        """
        self.conkitmatch['raw'] = self.conpred.deepcopy()
        if self.conpred_source == 'psicov':
            rscmin = 0.0
            rscmax = 0.0
            for contact in self.conkitmatch['raw']:
                if contact.raw_score < rscmin:
                    rscmin = contact.raw_score
                if contact.raw_score > rscmax:
                    rscmax = contact.raw_score
        if filterout is not None:
            self.conkitmatch['raw'] = filter_contacts(self.conpred.deepcopy(),
                                                      threshold=filterout)
        else:
            self.conkitmatch['raw'] = self.conpred.deepcopy()

        if self.conpred_source == 'psicov':
            psicovmodes = PSICOV_modes()
            for pm in psicovmodes:
                self.conkitmatch[pm] = self.conpred.deepcopy()
            for contact in self.conkitmatch['shifted']:
                contact.raw_score -= rscmin
            for contact in self.conkitmatch['norm']:
                contact.raw_score -= rscmin
                contact.raw_score /= rscmax
            for contact in self.conkitmatch['abs']:
                contact.raw_score = abs(contact.raw_score)

            if filterout is not None:
                for pm in psicovmodes:
                    self.conkitmatch[pm] = filter_contacts(self.conkitmatch[pm],
                                                           threshold=filterout)

        self.tp['raw'] = 0
        self.fp['raw'] = 0
        self.tn['raw'] = 0
        self.fn['raw'] = 0
        if self.conpred_source == 'psicov':
            for pm in psicovmodes:
                self.tp[pm] = 0
                self.fp[pm] = 0
                self.tn[pm] = 0
                self.fn[pm] = 0
        structuremap = self.interface.structure[1].deepcopy()
        for altsc, cmap in self.conkitmatch.items():
            logging.info('Structure: ' + str(structuremap) +
                         ', Conpred: ' + str(cmap) +
                         ', Source: ' + self.conpred_source +
                         ', Mode: ' + altsc)
            if len(cmap) > 0 and len(structuremap) > 0:
                self.ckplotmatch[altsc] = cmap.deepcopy()
                cmap.match(structuremap, add_false_negatives=True, inplace=True)
                self.ckplotmatch[altsc].match(structuremap,
                                              match_other=True,
                                              remove_unmatched=True,
                                              renumber=True)
                for contact in cmap:
                    if contact.true_positive:
                        self.tp[altsc] += 1
                    elif contact.false_positive:
                        self.fp[altsc] += 1
                    elif contact.false_negative:
                        self.fn[altsc] += 1
                    elif contact.true_negative:
                        logging.warning('True negatives appearing in conkit match.')
                    else:
                        logging.warning('Contact not evaluated.')
            else:
                logging.info('Contact map contains no contacts.')

        self.tn[altsc] = (self.npotential -
                          self.tp[altsc] - self.fp[altsc] - self.fn[altsc])


    def plot_map(self, outpath, mode='raw'):
        """
        Plot matched contact map.

        :param outpath: Outpath.
        :type outpath: str
        :param mode: Mode, if any, defaults to 'raw'.
        :type mode: str, optional

        """
        try:
            fig = ckplot.ContactMapFigure(self.ckplotmatch[mode],
                                          reference=self.interface.structure[1],
                                          legend=True)
            fig.savefig(outpath, overwrite=True)
            del fig
        except Exception:
            logging.warning('Something went wrong with ConKit ' +
                            'and Contact Plot was not produced.')
