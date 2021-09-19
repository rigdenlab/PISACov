# -*- coding: utf-8 -*-

import conkit
import os
import numpy as np

from .sequence import biofile

from .init2 import pdbid
from .output import output_tmpdir
from .init1 import minneigh
from .init1 import scorethreshold

#########################################################
##### CONTACT-RELATED FUNCTIONS #########################

# CONTACT LIST FILTERING

def filter_contacts(contactlist, source):
    """
    Remove neighbours and low-score contacts from contact list

    Parameters
    ----------
    contactlist : Conkit Contact map
        Original contact map.
    source : str
        Source of the contact list ("deepmetapsicov" or "psicov")

    Returns
    -------
    contactlist : Conkit Contact Map
        Filtered contact map.

    """

    if minneigh()>1:
        contactlist.remove_neighbors(min_distance=minneigh(), inplace=True)

    contactlist.sort("raw_score", reverse=True, inplace=True)

    if scorethreshold(source) > 0.0:
        cnt=0
        cntprev=1
        for contact in contactlist:
            if contact.raw_score < scorethreshold(source) and cntprev >= scorethreshold(source):
                cntprev=contact.raw_score
            elif contact.raw_score>=scorethreshold(source):
                cnt=cnt+1
            else:
                cntprev=0
        contactlist=contactlist[:cnt-1]

    return contactlist

def import_interfacepdb(niface,isitbio,outdir=output_tmpdir("pisacov")):
    """
    Import Interface PDB files created by PISA and renumbered afterwards.

    Parameters
    ----------
    niface : int
        Loop index in range(n_interfaces)
    isitbio : bool
        Specify whether using the ".bio" label in files or not.
    outdir : TYPE, optional
        Path to directory containing PDB files. The default is output_tmpdir("pisacov").

    Returns
    -------
    mapping : ConKit contact map.
        Contact map generated from the PDB file by ConKit.

    """
    pdbfile=pdbid()+".interface."+str(niface+1)+biofile(isitbio)+".fasta.pdb"

    tmppath = os.path.join(outdir, pdbfile)
    mapping=conkit.io.read(tmppath, 'pdb')

    return mapping

def remove_intra_contacts(conpred_list,intrapdb):
    """
    Remove contacts from the prediction list that are also identified as intramolecular contacts

    Parameters
    ----------
    conpred_list : Conkit Contact Map (prediction list)
        Contact prediction list.
    intrapdb : Conkit Contact Map (pdb map)
        Intramolecular contact map.

    Returns
    -------
    outmap : Conkit Contact Map (prediction list)
        Contact prediction list without unwanted contacts.

    """
    outmap=conpred_list.deepcopy()

    for contact in intrapdb:
        contactid = contact.id
        for contact2 in conpred_list:
            contact2id = contact2.id
            if contactid[0] == contact2id[0] and contactid[1] == contact2id[1]:
                outmap.remove(tuple(contact2id))
            elif contactid[0] == contact2id[1] and contactid[1] == contact2id[0]:
                outmap.remove(tuple(contact2id))

    return outmap

def pred_pdb_matchlist(conpred_in,pdbcontactfilepath):
    """
    Takes in a contact prediction list and a pdb contact prediction list and returns a Conkit Contact Map containing matching contacts.

    Parameters
    ----------
    conpred_in : Conkit Contact Map
        Contact prediction list
    pdbcontactfilepath : str
        File path to contact list file containing interface pdb contacts

    Returns
    -------
    conpred_out : Conkit Contact Map
        Matching contact list.

    """
    pdbcontacts=np.loadtxt(pdbcontactfilepath)
    conpred_out=conkit.core.ContactMap("tmp")
    npdbcontacts = 1 if len(pdbcontacts.shape) == 1 else pdbcontacts.shape[0]

    for iconpred in conpred_in:
        if npdbcontacts == 1:
            if int(iconpred.res1_seq)==int(pdbcontacts[0]) and int(iconpred.res2_seq)==int(pdbcontacts[1]):
                try:
                    conpred_out.add(iconpred)
                except:
                    pass
            elif int(iconpred.res2_seq)==int(pdbcontacts[0]) and int(iconpred.res1_seq)==int(pdbcontacts[1]):
                try:
                    conpred_out.add(iconpred)
                except:
                    pass
        elif npdbcontacts > 1:
            for ipdb in range(npdbcontacts):
                if int(iconpred.res1_seq)==int(pdbcontacts[ipdb][0]) and int(iconpred.res2_seq)==int(pdbcontacts[ipdb][1]):
                    try:
                        conpred_out.add(iconpred)
                    except:
                        pass
                elif int(iconpred.res2_seq)==int(pdbcontacts[ipdb][0]) and int(iconpred.res1_seq)==int(pdbcontacts[ipdb][1]):
                    try:
                        conpred_out.add(iconpred)
                    except:
                        pass

    conpred_out.sort("raw_score", reverse=True, inplace=True)

    return conpred_out

def match_maps(conpredlist,pdbmap):
    """
    Obtain several scores for an interface matched contact prediction map

    Parameters
    ----------
    conpredlist : Conkit Contact Map
        Contact prediction list
    pdbmap : str
        File path to contact list file containing interface pdb contacts

    Returns
    -------
    matchedmap : Conkit Contact Map
        The matched map (pdb interface + contact prediction list).
    cnt : int
        The number of True Positives.
    jaccardindex : float
        The Jaccard Index
    AvScore : float
        The average score of all True Positives.
    AccScore : float
        The accumulated score of all True Positives.
    probability : float
        The interface score.
    """
    matchedmap = conpredlist.match(pdbmap, add_false_negatives=False,inplace=False)
    cnt=0
    AccScore=0.0
    probability = 1.0
    for contact in matchedmap:
        if contact.true_positive:
            cnt += 1
            AccScore += contact.raw_score
            probability *= (1.0 - contact.raw_score)
    jaccardindex=conpredlist.get_jaccard_index(pdbmap)
    probability = 1.0 - probability
    if cnt > 0:
        AvScore = AccScore / cnt
    else:
        AvScore=0

    return matchedmap, cnt, jaccardindex, AvScore, AccScore, probability
