"""
This is PISACov, a program designed to infer quaternary structure
of proteins from evolutionary covariance.
"""

from pisacov import __prog__, __description__, __version__
from pisacov import __author__, __date__, __copyright__

from pisacov.iomod import _conf_ops as pco
from pisacov.core import _psicov_modes as PSICOV_modes

import copy
import logging
import os

from crops.elements import sequences as pes
from crops.iomod import taggers as ctg
from conkit.core import contactmap as ckc
from conkit.core.contactmap import ContactMap
from conkit.core.contact import Contact as ckContact
from conkit.core.mappings import ContactMatchState as ckCMS
from conkit import plot as ckplot

from matplotlib import pyplot as plt


def backmapping(cmap, sequence):
    """Return the contact prediction map with the original residue numbers.

    :param cmap: Contact prediction map.
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
    """Remove low score contacts from contact prediction list.

    :param cmap: Contact prediction map.
    :type cmap: :class:`~conkit.core.contactmap.ContactMap`
    :param threshold: Threshold, defaults to 0.2.
    :type threshold: float, str, optional

    """
    cmap.sort('raw_score', reverse=True, inplace=True)
    cnt = 0
    if threshold == "inf":
        threshold = 1000000000000
    elif threshold == "-inf":
        threshold = -1000000000000
    for contact in cmap:
        if contact.raw_score < threshold:
            break
        else:
            cnt = cnt+1

    return cmap[:cnt-1]

def map_intersection(conpredmap, strcon):
    """Generate a Contact map with the common contacts of conpredmap and strcon.

    :param conpredmap: Conkit Contact Map (contact prediction).
    :type conpredmap: :class:`~conkit.core.contactmap.ContactMap`
    :param strcon: Structural contact map (generated from pdb structure).
    :type strcon: :class:`~conkit.core.contactmap.ContactMap`
    :return: Intersection map.
    :rtype: :class:`~conkit.core.contactmap.ContactMap`

    """
    newmap = ContactMap(id=strcon.id)

    for ic2 in conpredmap:
        for ic1 in strcon:
            if (int(ic2.res1_seq) == int(ic1.res1_seq) and
                    int(ic2.res2_seq) == int(ic1.res2_seq)):
                try:
                    newmap.add(ic2)
                except Exception:
                    pass
            elif (int(ic2.res2_seq) == int(ic1.res1_seq) and
                  int(ic2.res1_seq) == int(ic1.res2_seq)):
                try:
                    newmap.add(ic2)
                except Exception:
                    pass

    newmap.sort("raw_score", reverse=True, inplace=True)

    return newmap

def map_intersection_array(conpredmap, strconarray):
    """Generate a Contact map with the common contacts of conpredmap and strconarray.

    :param conpredmap: Conkit Contact Map (contact prediction).
    :type conpredmap: :class:`~conkit.core.contactmap.ContactMap`
    :param strconarray: Numpy-read contact list (generated from pdb structure).
    :type strconarray: :class:`~numpy.ndarray`
    :return: Intersection map.
    :rtype: :class:`~conkit.core.contactmap.ContactMap`

    """
    newmap = ContactMap(id=strconarray.id)

    ncst = 1 if len(strconarray.shape) == 1 else strconarray.shape[0]
    for ic2 in conpredmap:
        if ncst == 1:
            if (int(ic2.res1_seq) == int(strconarray[0]) and
                    int(ic2.res2_seq) == int(strconarray[1])):
                try:
                    newmap.add(ic2)
                except Exception:
                    pass
            elif (int(ic2.res2_seq) == int(strconarray[0]) and
                  int(ic2.res1_seq) == int(strconarray[1])):
                try:
                    newmap.add(ic2)
                except Exception:
                    pass
        elif ncst > 1:
            for ic1 in range(ncst):
                if (int(ic2.res1_seq) == int(strconarray[ic1][0]) and
                        int(ic2.res2_seq) == int(strconarray[ic1][1])):
                    try:
                        newmap.add(ic2)
                    except Exception:
                        pass
                elif (int(ic2.res2_seq) == int(strconarray[ic1][0]) and
                      int(ic2.res1_seq) == int(strconarray[ic1][1])):
                    try:
                        newmap.add(ic2)
                    except Exception:
                        pass

    newmap.sort("raw_score", reverse=True, inplace=True)

    return newmap


class contact_atlas:
    """
    A :class:`~pisacov.core.contacts.contact_atlas` object containing information from
    sequences, contact prediction maps and structure contacts.
    The :class:`~pisacov.core.contacts.contact_atlas` class represents a data structure to hold
    contact maps, matched and unmatched with structure contacts and sequence.

    :param name: Atlas identifier, defaults to None.
    :type name: str, optional
    :param sequence: A sequence object, defaults to None.
    :type sequence: :class:`~crops.elements.sequences.sequence`, optional
    :param dimer_interface: An interface object, defaults to None.
    :type dimer_interface: :class:`~pisacov.core.interfaces.interface`, optional
    :param conpredmap: A Contact prediction map, defaults to None.
    :type conpredmap: :class:`~conkit.core.contactmap.ContactMap`, optional
    :param conpredtype: Source of contact prediction list (one of 'psicov', 'deepmetapsicov', 'ccmpred'), defaults to None.
    :type conpredtype: str, optional

    :ivar name: Atlas identifier.
    :vartype name: str
    :ivar sequence: A sequence object.
    :vartype sequence: :class:`~crops.elements.sequences.sequence`
    :ivar interface: An interface object.
    :vartype interface: :class:`~pisacov.core.interfaces.interface`
    :ivar conpred_raw: A Contact prediction map, as originally parsed.
    :vartype conpred_raw: :class:`~conkit.core.contactmap.ContactMap`
    :ivar conpred: A Contact prediction map, after processing.
    :vartype conpred: :class:`~conkit.core.contactmap.ContactMap`
    :ivar conpred_source: Source of contact prediction list (one of 'psicov', 'deepmetapsicov', 'ccmpred').
    :vartype conpred_source: str
    :ivar conkitmatch: Dictionary of prediction-structure matched contact maps.
    :vartype conkitmatch: dict [str : :class:`~conkit.core.contactmap.ContactMap`]
    :ivar ckplotmatch: Dictionary of prediction-structure matched contact maps (for plotting).
    :vartype ckplotmatch: dict [str : :class:`~conkit.core.contactmap.ContactMap`]
    :ivar tp: Dictionary of number of true positives.
    :vartype tp: dict [str : int]
    :ivar tn: Dictionary of number of false positives.
    :vartype tn: dict [str : int]
    :ivar fp: Dictionary of number of true negatives.
    :vartype fp: dict [str : int]
    :ivar fn: Dictionary of number of false negatives.
    :vartype fn: dict [str : int]
    :ivar npotential: Number of potential contacts.
    :vartype npotential: int

    :example:

    >>> from pisacov.core import contacts as pcc
    >>> from pisacov.core import interfaces as pci
    >>> from crops.iomod import parsers as cps
    >>> from conkit import io as ckio
    >>> myseq = cps.parseseqfile('7M6C.fasta')  # Parse sequence with CROPS
    >>> myseq['7m6c']
    Multiple sequence object: (id=7m6c, sequences = {'1': Sequence object >7M6C_1|Chain A (seq=MRTLWIMAVL[...]KPLCKKADPC, type=Undefined, length=138)})
    >>> myseq['7m6c'].imer['1'].seqs['conkit'] = ckio.read('7M6C.fasta', 'fasta')  # Add to sequence the Conkit-parsed version for later use
    >>> myseq['7m6c'].imer['1'].seqs
    {'mainseq': 'MRTLWIMAVLLVGVEGSLVELGKMILQETGKNPVTSYGAYGCNCGVLGRGKPKDATDRCCSVHKCCYKKMTGCNPKKDRYSYSWKDKTIVCDENNPCLKELCECDKAVAICLRENLDTYNEKYKKYYKKPLCKKADPC',
     'conkit': Sequence(id="7M6C_1|Chain A|Basic phospholipase A2|Bothrops atrox (8725)" seq="MRTLW...KADPC" seq_len=138)}
    >>> myinterfacelist = []
    >>> myinterfacelist.append(pci.interface(name='1'))  # Create new interface and add chains
    >>> myinterfacelist[0].chains.append(pci.chain_info(dimer_id='A',
                                                        pisa_id='1',
                                                        crystal_id='A',
                                                        seqid='1',
                                                        biotype='Protein'))
    >>> myinterfacelist[0].chains.append(pci.chain_info(dimer_id='B',
                                                        pisa_id='1',
                                                        crystal_id='A',
                                                        seqid='1',
                                                        biotype='Protein'))
    >>> myinterfacelist[0].stable = False  # Define the stability of interface
    >>> myinterfacelist[0]
    Interface object 1 (chains=A,A, type=Protein-Protein, stable=False)
    >>> myinterfacelist[0].chains[0]
    Chain Info object: (Dimer ID = A, PISA xml ID = 1, Monomer ID in assymmetric unit = A, Sequence ID = 1, Biotype = Protein)
    >>> inputmap = ckio.read('7m6c.interface.1.pdb', 'pdb')  # Parse interface dimer structure
    >>> myinterfacelist[0].structure = []
    >>> for m in range(len(inputmap)):  # Add interface structure to interface object
            myinterfacelist[0].structure.append(inputmap[m].as_contactmap())
            myinterfacelist[0].structure[m].id = inputmap[m].id
    >>> myinterfacelist[0].structure
    [ContactMap(id="A", ncontacts=572),
     ContactMap(id="AB", ncontacts=41),
     ContactMap(id="BA", ncontacts=41),
     ContactMap(id="B", ncontacts=572)]
    >>> myconpredmap = ckio.read('7m6c.psicov', 'psicov')  # Parse contact prediction list
    >>> myatlas = pcc.contact_atlas(name='7M6C_1',  # Create atlas
                                    sequence=myseq['7m6c'].imer['1'],
                                    dimer_interface=myinterfacelist[0],
                                    conpredmap=myconpredmap,
                                    conpredtype='psicov')
    >>> myatlas
    Contact Atlas object 7M6C_1 (Interface=Interface object 1 (chains=A,A, type=Protein-Protein, stable=False), Contact Prediction=ContactMap(id="map_1", ncontacts=8298))
    >>> myatlas.set_conpred_seq() # Assign sequence to contact prediction list
    >>> myatlas.remove_neighbours()
    >>> myatlas.remove_intra()
    >>> myatlas.make_match(filterout=0.2, tpmap=False)  # Match Structure and Prediction
    >>> myatlas
    Contact Atlas object 7M6C_1 (Interface=Interface object 1 (chains=A,A, type=Protein-Protein, stable=False), Contact Prediction=ContactMap(id="map_1", ncontacts=8092), True Positives=0)
    >>> myatlas.tp
    >>> myatlas.fp
    >>> myatlas.fn
    >>> myatlas.plotmap('myoutpath/plotfile.png')

    """

    _kind = 'Contact Atlas'
    __slots__ = ['name', 'interface', 'sequence',
                 'conpred_raw', 'conpred', 'conpred_source',
                 'conkitmatch', 'ckplotmatch',
                 'tp', 'tn', 'fp', 'fn', 'npotential']

    def __init__(self, name=None, sequence=None, dimer_interface=None,
                 conpredmap=None, conpredtype=None):

        self.name = name
        self.interface = dimer_interface
        self.conpred_raw = conpredmap
        if conpredmap is not None:
            self.conpred = self.conpred_raw.deepcopy()
        else:
            self.conpred = None
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
        if sequence is not None:
            lseq = sequence.length()
            self.npotential = (lseq**2 - lseq) / 2
            self.npotential = int(self.npotential / 2)

    def __repr__(self):
        string = (self._kind+" object "+self.name +
                  " (Interface=" + str(self.interface) + ", " +
                  "Contact Prediction=" + str(self.conpred))
        if 'raw' in self.tp:
            string += ", True Positives=" + str(self.tp['raw'])
        string += ")"
        return string

    def copy(self):
        return copy.copy(self)

    def deepcopy(self):
        return copy.deepcopy(self)

    def remove_neighbors(self, mindist=2):
        md = mindist
        self.remove_neighbours(mindist=md)

    def remove_neighbours(self, mindist=2):
        """Return :class:`~conkit.core.contactmap.ContactMap` without neighbouring pairs.

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
        """Remove intramolecular contacts from intermolecular :class:`~conkit.core.contactmap.ContactMap`."""
        for m in [0, 3]:
            intra = self.interface.structure[m]
            for contact1 in intra:
                if isinstance(contact1.id, tuple):
                    c1 = contact1.id
                elif isinstance(contact1.id, str):
                    c1 = contact1.id[1:-1].split(', ')
                for contact2 in reversed(self.conpred):
                    if isinstance(contact2.id, tuple):
                        c2 = contact2.id
                    elif isinstance(contact2.id, str):
                        c2 = contact2.id[1:-1].split(', ')
                    if ((c1[0] == c2[0] and c1[1] == c2[1]) or
                            (c1[1] == c2[0] and c1[0] == c2[1])):
                        self.conpred.remove(contact2.id)
                for contact2 in reversed(self.interface.contactmap):
                    if isinstance(contact2.id, tuple):
                        c2 = contact2.id
                    elif isinstance(contact2.id, str):
                        c2 = contact2.id[1:-1].split(', ')
                    if ((c1[0] == c2[0] and c1[1] == c2[1]) or
                            (c1[1] == c2[0] and c1[0] == c2[1])):
                        self.interface.contactmap.remove(contact2.id)

    def set_sequence(self, sequence):
        """
        Set Atlas sequence.

        This is the correct way to update the sequence.
        Storing directly in :attr:`~pisacov.core.contacts.contact_atlas.sequence` will have sequence-dependent values not updated.
        :param sequence: A sequence object.
        :type sequence: :class:`~crops.elements.sequences.sequence`

        """
        lseq = sequence.length()
        self.npotential = (lseq**2 - lseq) / 2
        self.npotential = int(self.npotential / 2)

    def set_cropmap(self):
        """Renumber :class:`~conkit.core.contactmap.ContactMap` according to :class:`~crops.elements.sequences.sequence`."""
        self.conpred = backmapping(self.conpred, self.sequence)

    def set_conpred_seq(self, sequence=None):
        """Set contact prediction sequence.

        :param sequence: Conkit type sequence, defaults to self.sequence.seqs['conkit'].
        :type mindist: :class:`~conkit.core.sequence.Sequence`, optional

        """
        seq_in = self.sequence.seqs['conkit'] if sequence is None else sequence
        self.conpred.sequence = seq_in
        self.conpred.set_sequence_register()

    def make_match(self, filterout=None, tpmap=False):
        """Match Structure and contact prediction maps.

        :param filterout: Threshold score below which contacts are filtered out, defaults to None.
        :type filterout: float, str, optional
        :param tpmap: If True, only consider conpred's TPs, defaults to False.
        :type tpmap: bool, optional

        """
        if tpmap is False:
            self.conkitmatch['raw'] = self.conpred.deepcopy()
        else:
            self.conkitmatch['raw'] = map_intersection(self.conpred,
                                                       self.interface.contactmap)
        if self.conpred_source == 'psicov':
            rscmin = 0.0
            rscmax = 0.0
            for contact in self.conkitmatch['raw']:
                if contact.raw_score < rscmin:
                    rscmin = contact.raw_score
                if contact.raw_score > rscmax:
                    rscmax = contact.raw_score
        if filterout is not None:
            self.conkitmatch['raw'] = filter_contacts(self.conkitmatch['raw'].deepcopy(),
                                                      threshold=filterout)
        else:
            pass

        if self.conpred_source == 'psicov':
            psicovmodes = PSICOV_modes()
            for pm in psicovmodes:
                self.conkitmatch[pm] = self.conkitmatch['raw'].deepcopy()
            for contact in self.conkitmatch['shifted']:
                contact.raw_score -= rscmin
            for contact in self.conkitmatch['norm']:
                contact.raw_score -= rscmin
                contact.raw_score /= rscmax
            for contact in self.conkitmatch['abs']:
                contact.raw_score = abs(contact.raw_score)

            if filterout is not None:
                for pm in psicovmodes:
                    self.conkitmatch[pm] = filter_contacts(self.conkitmatch[pm].deepcopy(),
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
        structuremap = self.interface.contactmap
        for altsc, cmap in self.conkitmatch.items():
            logging.info('Structure: ' + str(structuremap) +
                         ', Conpred: ' + str(cmap) +
                         ', Source: ' + self.conpred_source +
                         ', Mode: ' + altsc)
            if len(cmap) > 0 and len(structuremap) > 0:
                ##self.conkitmatch[altsc] = ContactMap()
                self.conkitmatch[altsc] = cmap.deepcopy()
                #self.conkitmatch[altsc] = cmap.match_naive(structuremap,
                #                                          add_false_negatives=True,
                #                                         inplace=False,
                #                                         match_other=True)
                for c1 in cmap:
                    if c1.id[0] != c1.id[1]:
                        self.conkitmatch[altsc].add(ckContact(c1.id[1], c1.id[0], c1.raw_score))
                self.conkitmatch[altsc].sort('raw_score', reverse=True, inplace=True)
                # Conkit ckCMS.status: {0: Unk, 1: TP, 2: TN, 3: FP, 4: FN}
                for c2 in self.interface.contactmap:
                    c2.status = ckCMS.unknown
                for c1 in self.conkitmatch[altsc]:
                    c1.status = ckCMS.unknown
                    for c2 in self.interface.contactmap:
                        if c1.id[0] == c2.id[0] and c1.id[1] == c2.id[1]:
                            c1.status = ckCMS.true_positive
                            c2.status = ckCMS.true_positive
                            self.tp[altsc] += 1
                for c2 in self.interface.contactmap:
                    if c2.true_positive is False:
                        c2.status = ckCMS.false_negative
                        self.conkitmatch[altsc].add(c2)
                        self.fn[altsc] += 1
                for c1 in self.conkitmatch[altsc]:
                    if c1.true_positive is False and c1.false_negative is False:
                        c1.status = ckCMS.false_positive
                        self.fp[altsc] += 1

                    #if contact.true_positive:
                    #    self.tp[altsc] += 1
                    #elif contact.false_positive:
                    #    self.fp[altsc] += 1
                    #elif contact.false_negative:
                    #    self.fn[altsc] += 1
                    #elif contact.true_negative:
                    #    logging.warning('True negatives appearing in conkit match.')
                    #else:
                    #    logging.warning('Contact ' + str(contact.id) + ' not evaluated.')
            else:
                logging.info('Contact map contains no contacts.')

            self.tn[altsc] = (self.npotential -
                              self.tp[altsc] - self.fp[altsc] - self.fn[altsc])

    def plot_map(self, outpath, mode='raw'):
        """Plot matched contact map.

        :param outpath: Path to output file.
        :type outpath: str
        :param mode: Mode, if any, defaults to 'raw'.
        :type mode: str, optional

        """
        try:
            fig = ckplot.ContactMapFigure(self.conkitmatch[mode],
            # fig = ckplot.ContactMapFigure(self.ckplotmatch[mode],
                                          reference=self.interface.structure[1],
                                          legend=True,
                                          lim=(1, self.sequence.full_length()))
            fig.savefig(outpath)
            plt.close(fig.fig)
        except Exception:
            logging.warning('Something went wrong with ConKit ' +
                            'and Contact Plot was not produced.')


    def plot_map_alt(self, outpath, mode='raw', plot_type='png', xL=None):
        """Plot matched contact map.

        :param outpath: Path to output file.
        :type outpath: str
        :param mode: Mode, if any, defaults to 'raw'.
        :type mode: str, optional
        :param plot_type: Plot either as a 'png' image, 'eps' or 'svg' vector image, or 'agr' in raw grace format, defaults to 'png'.
        :type plot_type: str, optional
        :param xL: Number of contacts plotted as a function of L, defaults to None (all contacts).
        :type xL: int or float, optional

        """
        import matplotlib.pyplot as plt

        if xL is not None:
            nc = round(self.sequence.length()*xL)
        else:
            nc = len(self.conkitmatch[mode])

        fpx = []
        fpy = []
        tpx = []
        tpy = []
        fnx = []
        fny = []
        # TO DO: Adapt title to non-crystal dimer inputs
        title = (self.name + ', ' + 'Interface ' + self.interface.name +
                 ', Chains ' + self.interface.chains[0].crystal_id
                 + self.interface.chains[1].crystal_id +
                 ', ' + self.conpred_source)
        if self.conpred_source == 'psicov':
            title += ' (' + mode + ')'

        n = 0
        for contact in self.conkitmatch[mode]:
            c1 = contact.id[0]
            c2 = contact.id[1]
            if contact.true_positive is True and n < nc:
                n += 1
                tpx.append(c1)
                tpy.append(c2)
            elif contact.false_positive is True and n < nc:
                n += 1
                fpx.append(c1)
                fpy.append(c2)
            else:
                fnx.append(c1)
                fny.append(c2)

        if plot_type == 'png':
            fig, ax = plt.subplots(dpi=141)
        elif plot_type == 'eps' or plot_type == 'svg':
            fig, ax = plt.subplots(dpi=1200)
        elif plot_type == 'agr':
            tpx, tpy = (list(e) for e in zip(*sorted(zip(tpx, tpy))))
            fpx, fpy = (list(e) for e in zip(*sorted(zip(fpx, fpy))))
            fnx, fny = (list(e) for e in zip(*sorted(zip(fnx, fny))))
            xdat = [tpx, fpx, fnx]
            ydat = [tpy, fpy, fny]
            with open(outpath, 'w') as fout:
                for p in range(3):
                    fout.write('@target G0.S' + str(p) + '\n@type xy\n')
                    for n in range(len(xdat[p])):
                        fout.write(str(xdat[p][n]) + '  ' + str(ydat[p][n]) + '\n')
                    fout.write('&\n')
            return

        ax.set_title(title, y=1.08)

        vmin = 1
        vmax = self.sequence.length()
        ax.axis([vmin, vmax, vmin, vmax])
        ax.set_xlim(vmin - 0.5, vmax + 0.5)
        ax.set_ylim(vmin - 0.5, vmax + 0.5)

        ax.set_xlabel('Residues from Chain 1')
        ax.set_ylabel('Residues from Chain 2')

        s = ((ax.get_window_extent().width  / (vmax-vmin+1.) * 50./fig.dpi) ** 2)
        ax.scatter(tpx, tpy, s=s, marker='o', linewidth=0, c='k', label='Matched (TP)')
        ax.scatter(fpx, fpy, s=s, marker='o', linewidth=0, c='r', label='Unmatched (FP)')
        ax.scatter(fnx, fny, s=s, marker='o', linewidth=0, c='lightgrey', label='Structure (FN)')

        ax.legend(numpoints=1,
                  fontsize=10,
                  bbox_to_anchor=(0.0, 1.02, 1.0, 0.102),
                  loc=3,
                  ncol=3,
                  mode="expand",
                  borderaxespad=0.0)
        for n in range(len(ax.legend_.legendHandles)):
            ax.legend_.legendHandles[n]._sizes = [30]
        if plot_type == 'png' or plot_type == 'eps':
            fig.savefig(outpath, format=plot_type)
            plt.close(fig)
