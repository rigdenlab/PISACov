"""
This is PISACov, a PISA extension to infer quaternary structure
of proteins from evolutionary covariance.
"""

from pisacov import __prog__, __description__, __version__
from pisacov import __author__, __date__, __copyright__

from pisacov.io import _conf_ops as pco

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
    for n in range(len(conpredmap)):
        c1 = sequence.cropbackmap[conpredmap[n].res1_seq]
        c2 = sequence.cropbackmap[conpredmap[n].res2_seq]
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


def merge_interfaces(contactmap1, contactmap2):
    """
    Return a contact map that is the result of merging two of them.

    :param contactmap1: Base contact map.
    :type contactmap1: :class:`~conkit.core.contactmap.ContactMap`
    :param contactmap2: Secondary contact map.
    :type contactmap2: :class:`~conkit.core.contactmap.ContactMap`
    :return: Base contact map including also contacts from secondary contact map.
    :rtype: :class:`~conkit.core.contactmap.ContactMap`

    """
    cnt0 = 0
    cnt1 = 0
    cnt2 = 0
    outmap = contactmap1.deepcopy()
    print(outmap)
    print('===')
    for n2 in range(len(contactmap2)):
        repeated = False
        k1 = contactmap2[n2].res1_seq
        k2 = contactmap2[n2].res2_seq
        for n1 in range(len(contactmap1)):
            c1 = contactmap1[n1].res1_seq
            c2 = contactmap1[n1].res2_seq
            if (k1 == c1 and k2 == c2) or (k2 == c1 and k1 == c2):
                repeated = True
                break
        if repeated is False:
            cnt2 += 1
            outmap.add(contactmap2[n2])
    print(outmap)
    return outmap

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
                 'conpred_raw', 'conpred', 'conpred_source',
                 'strmap_raw', 'strmap', 'conkitmatch', 'intramap',
                 'tp_raw', 'tn_raw', 'fp_raw', 'fn_raw',
                 'tp', 'tn', 'fp', 'fn', 'npotential']
    def __init__(self, name, conpredmap, structure, strmap, conpredtype,
                 sequence, minneigh=2, removeintra=True, backmap=False):
        self.name = name
        self.chain1 = None
        self.chain2 = None
        self.sequence = sequence
        self.conpred_raw = conpredmap
        self.conpred = None
        self.conpred_source = conpredtype
        self.strmap_raw = structure
        self.strmap = {}
        self.conkitmatch = {}
        self.tp = {}
        self.fp = {}
        self.tn = {}
        self.fn = {}
        self.npotential = None
        # Set sequence values and Number of total potential contacts
        lseq = sequence.length()
        self.npotential = lseq**2 - lseq
        for n in range(1, minneigh):
            self.npotential -= 2*(lseq - n)
        self.npotential = int(self.npotential / 2)
        # Set Contact prediction list
        self.conpred = remove_neighbours(conpredmap, mindist=minneigh)
        self.conpred.sort('raw_score', reverse=True, inplace=True)
        if backmap is True:
            self.conpred = backmapping(self.conpred, sequence)
        self.conpred.sequence = self.sequence.seqs['conkit']
        self.conpred.set_sequence_register()
        # Set Structure map
        self.chain1 = tuple(self.strmap_raw[1].id)[0]
        self.chain2 = tuple(self.strmap_raw[1].id)[1]
        # Remove intramolecular contacts
        if removeintra is True:
            for m in [0,3]:
                intra = self.strmap_raw[m]
                for contact1 in intra:
                    c1 = str(contact1.id)[1:-1].split(', ')
                    for contact2 in reversed(self.conpred):
                        c2 = str(contact2.id)[1:-1].split(', ')
                        if ((c1[0] == c2[0] and c1[1] == c2[1]) or
                                (c1[1] == c2[0] and c1[0] == c2[1])):
                            self.conpred.remove(contact2.id)  # CHECK THAT REMOVAL INSIDE LOOP IS OK
        # Contact prediction maps
        self.conkitmatch['raw'] = self.conpred.deepcopy()
        if self.conpred_source == 'psicov':
            self.conkitmatch['abs'] = self.conpred.deepcopy()
            self.conkitmatch['shifted'] = self.conpred.deepcopy()
            self.conkitmatch['norm'] = self.conpred.deepcopy()
            rscmin = 0.0
            rscmax = 0.0
            for contact in self.conkitmatch['raw']:
                if contact.raw_score < rscmin:
                    rscmin = contact.raw_score
                if contact.raw_score > rscmax:
                    rscmax = contact.raw_score
            for contact in self.conkitmatch['shifted']:
                contact.raw_score -= rscmin
            for contact in self.conkitmatch['norm']:
                contact.raw_score -= rscmin
                contact.raw_score /= rscmax
            for contact in self.conkitmatch['abs']:
                contact.raw_score = abs(contact.raw_score)
#        threshold = pco._sources()[self.conpred_source][3]
#        for contact in reversed(self.conkitmatch):
#            if contact.raw_score < threshold:
#                self.conkitmatch.remove(contact.id)
#        if self.conpred_source == 'psicov':
#            for cmap in self.conkitmatch_alt.values():
#                for contact in reversed(cmap):
#                    if contact.raw_score < threshold:
#                        self.conkitmatch.remove(contact.id)
        # Interface map
        self.strmap['raw'] = ckc.ContactMap(self.strmap_raw[1].id)
        if self.conpred_source == 'psicov':
            self.strmap['abs'] = ckc.ContactMap(self.strmap_raw[1].id)
            self.strmap['shifted'] = ckc.ContactMap(self.strmap_raw[1].id)
            self.strmap['raw'] = ckc.ContactMap(self.strmap_raw[1].id)
        npdbcontacts = 1 if len(strmap) == 1 else strmap.shape[0]
        for altsc, cmap in self.conkitmatch.items():
            for iconpred in cmap:
                if npdbcontacts == 1:
                    if ((int(iconpred.res1_seq) == int(strmap[0]) and
                            int(iconpred.res2_seq) == int(strmap[1])) or
                            (int(iconpred.res1_seq) == int(strmap[1]) and
                            int(iconpred.res2_seq) == int(strmap[0]))):
                        try:
                            self.strmap[altsc].add(iconpred)
                        except Exception:
                            pass
                elif npdbcontacts > 1:
                    for ipdb in range(npdbcontacts):
                        if ((int(iconpred.res1_seq) == int(strmap[ipdb][0]) and
                                int(iconpred.res2_seq) == int(strmap[ipdb][1])) or
                                (int(iconpred.res1_seq) == int(strmap[ipdb][1]) and
                                int(iconpred.res2_seq) == int(strmap[ipdb][0]))):
                            try:
                                self.strmap[altsc].add(iconpred)
                            except Exception:
                                pass
            #cmap.sequence = self.sequence.seqs['conkit']
            #cmap.set_sequence_register()
        # MATCH
        for altsc, cmap in self.conkitmatch.items():
            cmap.match(self.strmap[altsc], add_false_negatives=True, inplace=True)
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

        self.tn[altsc] = (self.npotential -
                          self.tp[altsc] - self.fp[altsc] - self.fn[altsc])

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
