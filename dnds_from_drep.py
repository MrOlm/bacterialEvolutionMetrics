#!/usr/bin/env python

# Version 0.2.1 - 7.14.20
# Updated the substitution_matrix that seems to be different now?
# Added some additional debugging

# Version 0.1.3
# Matt Olm
# mattolm@berkeley.edu
# 02.20.19
# https://biotite.berkeley.edu/j/user/mattolm/notebooks/OtherProjects/RefSeq/goANI_2_development_3_dRep.ipynb

# CHANGED IN VERSION 0.1.3:
# Fixed a too many alignments bug with a workaround https://biotite.berkeley.edu/j/user/mattolm/notebooks/OtherProjects/RefSeq/goANI_2_development_4_debugging.ipynb

# CHANGED IN VERSION 0.2.0:
# updated to match https://biotite.berkeley.edu/j/user/mattolm/notebooks/OtherProjects/RefSeq/_SupplementalNotebooks_1_dnds_HGT_calculations.ipynb

__author__ = "Matt Olm"
__version__ = "0.2.1"
__license__ = "MIT"

import os
import sys
import glob
import scipy
import pickle
import shutil
import argparse
import textwrap
import datetime
import traceback
import numpy as np
import pandas as pd

import warnings
warnings.simplefilter("ignore")

import numpy as np
from math import log
import concurrent.futures
from itertools import permutations
from collections import defaultdict

import Bio
from Bio import Align
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.SubsMat import MatrixInfo
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from Bio.codonalign.codonalphabet import default_codon_table

from Bio.Align import substitution_matrices

def main(args):
    '''
    Main entry point
    '''
    wd = args.drep_dir
    assert os.path.isdir(wd)
    if wd[-1] != '/':
        wd = wd + '/'

    dnds_from_dRep(wd, p=args.processes, debug=args.debug)

def dnds_from_dRep(wd, p=6, debug=False):
    '''
    wrapper to calculate a bunch of dnds things from a dRep folder run with goANI

    arguments:
        wd = location of dRep work directory (must end in /)
        p = number of threads to use
        debug = a number of things that you normally don't want
    '''
    # Gather the .nsimscan files from processing (you should actually do this based on Ndb!)
    Ndb = pd.read_csv(wd + 'data_tables/Ndb.csv')
    if debug:
        Ndb = Ndb.head(10)
    ndb = Ndb[[True if (q < r) else False for q, r in zip(\
                Ndb['querry'], Ndb['reference'])]]

    # Make an output folder
    of = wd + 'data/dNdS_data/'
    if os.path.isdir(of):
        pass
    else:
        os.mkdir(of)

    # Thread the commands
    Rdb = pd.DataFrame()

    j = 0
    for Ndb in iterate_Ndbs(ndb, chunkSize=p*10):
        RNdb = pd.DataFrame()
        print("Pooling {0} dN/dS commands (which take ~5 min each)...".format(int(len(Ndb))))
        try:
            with concurrent.futures.ProcessPoolExecutor(max_workers=p) as executor:
                for cmd, result in executor.map(goANI_wrapper, iterate_commands(wd, Ndb)):
                    rdb = cmd.geneDb
                    rdb['querry'] = cmd.row['querry']
                    rdb['reference'] = cmd.row['reference']
                    rdb = pd.merge(rdb, result, on=['qry_id', 'sbj_id'], how='outer')
                    #rdb['querry'] = rdb['querry'].astype('category')
                    #rdb['reference'] = rdb['reference'].astype('category')

                    print("{2}: finished {0} vs {1}".format(cmd.row['querry'], cmd.row['reference'],
                        str(datetime.datetime.now())))

                    RNdb = RNdb.append(rdb)

            print('finished chunk {0}'.format(j))
            Rdb = Rdb.append(RNdb)
            RNdb.to_csv(of + 'detailed_dNdS_info_chunk_{0}.csv'.format(j), index=False)
        except Exception as e:
            print('Chunk {0} failed!!'.format(j))
            traceback.print_exc()

        j += 1

    # Save output
    Rdb.to_csv(of + 'detailed_dNdS_info.csv', index=False)
    print('Base calculations done- output saved to {0}'.format(of + 'detailed_dNdS_info.csv'))


    GDdb = parse_to_genomeLevel(Rdb)
    GDdb.to_csv(of + 'genomeWide_dNdS_info.csv', index=False)
    print('Overall calculations done- output saved to {0}'.format(of + 'genomeWide_dNdS_info.csv'))
    return Rdb, GDdb

def parse_to_genomeLevel(Ddb):
    Tdbs = []
    for query, qdb in Ddb.groupby('querry'):
        for ref, db in qdb.groupby('reference'):
            tdb = summarize_dnds_HGT(db, noFilt=True)
            Tdbs.append(tdb)
    return pd.concat(Tdbs)

from tqdm import tqdm

def summarize_dnds_HGT(db, noFilt=False, fastANI=False):
    table = defaultdict(list)

    # Basic info
    table['reference'].append(db['reference'].tolist()[0])
    table['querry'].append(db['querry'].tolist()[0])

    # Number of comparisons
    table['total_comps'].append(len(db))
    table['failed_comps'].append(len(db[db['N_sites'] == 0]))
    table['successful_comps'].append(len(db[db['N_sites'] > 0]))
    table['considered_bases'].append(db[(db['N_sites'] > 0)]['al_len'].sum())
    table['percet_successful'].append((len(db[db['N_sites'] > 0]) / len(db)) * 100)

    # Filter to successful comparisons
    db = db[db['N_sites'] > 0]

    # dN/dS stuff
    N = db['N_sites'].sum()
    table['N_sites'].append(N)
    dN = db['N_changed'].sum()
    table['N_changed'].append(dN)
    S = db['S_sites'].sum()
    table['S_sites'].append(S)
    dS = db['S_changed'].sum()
    table['S_changed'].append(dS)
    if (N > 0) & (S > 0) & (dS > 0):
        table['dN/dS'].append(((dN/N) / (dS/S)))
    else:
        table['dN/dS'].append(np.nan)

    # Calculate ANI if you must
    if fastANI == False:
        fastANI = sum([p * l for p, l in zip(db['p_inden'], db['al_len'])]) \
                    / db['al_len'].sum()
    table['fast_ani'].append(fastANI)

    # Filter to only ones with a full alignment
    if noFilt:
        hdb = db[(db['al_len'] > 500)]
    else:
        hdb = db[(db['al_len'] > 500) & (db['al_len'] < 1000)]

    ex_iden = _calc_expected(hdb, fastANI)
    table['counted_comps'].append(len(hdb))
    table['identical_comps'].append(len(hdb[hdb['p_inden'] > 99.99]))
    table['expected_identicle'].append(ex_iden)

    try:
        table['percent_enriched'].append(((len(hdb[hdb['p_inden'] > 99.99]) - int(ex_iden)) / len(hdb))*100)
        hani = sum([p * l for p, l in zip(hdb['p_inden'], hdb['al_len'])]) \
                / hdb['al_len'].sum()
        table['filtered_ani'].append(hani)
        table['filtered_af'].append(db['al_len'].sum() / db['qry_len'].sum())
    except:
        table['percent_enriched'].append(0)
        if not noFilt:
            table['filtered_ani'].append(0)
            table['filtered_af'].append(0)

    Ndb = pd.DataFrame(table)
    return Ndb

def _calc_expected(db, ani):
    expected = 0
    prob_diff = (100-ani) / 100
    prob_same = ani / 100
    for i, row in db.iterrows():
        # probability of this row being 0
        #p = prob_diff * int(row['al_len'])
        p = prob_same ** int(row['al_len'])
        if p > 1:
            expected += 1
        else:
            expected += p
    return expected

def parse_to_genomeLevel_HGT(Ddb, noFilt=False):
    Tdbs = []
    for query, qdb in tqdm(Ddb.groupby('querry')):
        for ref, db in qdb.groupby('reference'):
            # get the FastANI
            try:
                fastANI = _get_ani(FAdb, db['reference'].tolist()[0], db['querry'].tolist()[0])
            except:
                fastANI = 0

            tdb = summarize_dnds_HGT(db, noFilt=noFilt, fastANI=fastANI)
            Tdbs.append(tdb)
    return pd.concat(Tdbs)

def summarize_dnds(db):
    table = defaultdict(list)

    table['reference'].append(db['reference'].tolist()[0])
    table['querry'].append(db['querry'].tolist()[0])

    table['total_comps'].append(len(db))
    table['failed_comps'].append(len(db[db['N_sites'] == 0]))
    table['successful_comps'].append(len(db[db['N_sites'] > 0]))
    table['considered_bases'].append(db[(db['N_sites'] > 0)]['al_len'].sum())
    table['percet_successful'].append((len(db[db['N_sites'] > 0]) / len(db)) * 100)

    N = db['N_sites'].sum()
    table['N_sites'].append(N)
    dN = db['N_changed'].sum()
    table['N_changed'].append(dN)
    S = db['S_sites'].sum()
    table['S_sites'].append(S)
    dS = db['S_changed'].sum()
    table['S_changed'].append(dS)

    table['homolog_ani'].append(sum([p * l for p, l in zip(db['p_inden'], db['al_len'])]) \
                    / db['al_len'].sum())
    table['homolog_af'].append(db['al_len'].sum() / db['qry_len'].sum())

    if (N > 0) & (S > 0) & (dS > 0):
        table['dN/dS'].append(((dN/N) / (dS/S)))
    else:
        table['dN/dS'].append(np.nan)

    Ndb = pd.DataFrame(table)
    return Ndb

def iterate_Ndbs(df, chunkSize=100):
    '''
    Break up Ndbs into chunks
    '''
    numberChunks = len(df) // chunkSize + 1
    for i in range(numberChunks):
        yield (df[i*chunkSize:(i+1)*chunkSize])

def iterate_commands(wd, Ndb):
    '''
    Make and iterate commands

    This makes it require WAY less RAM
    '''
    # Make this into a list of commands
    for i, row in Ndb.iterrows():
        # You only want each combo once; you don't want self
        if row['querry'] >= row['reference']:
            continue

        # get the locations of the .nimsan files
        nsimscan_locs = _get_nsim_files(wd, row)

        # get a list of gense to align
        gedb = parse_goANI(nsimscan_locs[0], nsimscan_locs[1])
        #print(len(gedb))

        # get the fna and faa files
        Qfna, Qfaa, Rfna, Rfaa = _get_protein_files(wd, row)

        # make a command
        cmd = dnDs_command()
        cmd.geneDb = gedb
        cmd.row = row
        cmd.Qfna = Qfna
        cmd.Qfaa = Qfaa
        cmd.Rfna = Rfna
        cmd.Rfaa = Rfaa

        # verify command
        cmd.verify()

        # append it to the list
        print('generated command {0} vs {1}'.format(row['querry'], row['reference']))
        yield cmd

def _get_nsim_files(wd, row):
    '''
    From a wd and a row of Ndb, get the simscan files
    '''
    # Get the first direction
    cur_folder = wd + 'data/goANI_files/{0}'.format(row['querry'])
    fn1 = "{0}/{1}_vs_{2}.sim".format(cur_folder, \
                            row['querry'], row['reference'])

    # Get the reverse direction
    cur_folder = wd + 'data/goANI_files/{0}'.format(row['reference'])
    fn2 = "{0}/{1}_vs_{2}.sim".format(cur_folder, \
                            row['reference'], row['querry'])

    return [fn1, fn2]

def _get_protein_files(wd, row):
    '''
    From a wd and a row of Ndb, get prodigal files
    '''
    # Get the first files
    Qfna = "{0}data/prodigal/{1}.fna".format(wd, row['querry'])
    Qfaa = "{0}data/prodigal/{1}.faa".format(wd, row['querry'])

    # Get the second files
    Rfna = "{0}data/prodigal/{1}.fna".format(wd, row['reference'])
    Rfaa = "{0}data/prodigal/{1}.faa".format(wd, row['reference'])

    # Return
    return Qfna, Qfaa, Rfna, Rfaa

def parse_goANI(loc1, loc2):
    '''
    RIGHT NOW THIS IGNORES LOC2!!!

    Fromt the locations of two locations, return a list of genes to calculate dn/ds for
    '''
    # Parse first data-frame
    db1 = _quick_parse_filter(loc1)
    #db2 = _quick_parse_filter(loc2)

    # Get a list of best hits

    return db1

def _quick_parse_filter(loc):
    db1 = pd.read_csv(loc, sep='\t')
    db1 = db1.rename(columns={'#qry_id':'qry_id'})

    # Filter like gANI
    db1 = db1.sort_values(['al_len','qry_id', 'sbj_id'], ascending=False)
    db1['af'] = [a/min(o,t) for a,o,t in zip(db1['al_len'],
                                             db1['qry_len'],
                                             db1['sbj_len'])]
    db1 = db1[(db1['af'] >= 0.7) & (db1['p_inden'] >= 70)]

    # Only keep reciprical best hits
    q2b = db1.sort_values('sw_score', ascending=False).drop_duplicates(subset=['qry_id'])\
            .set_index('qry_id')['sbj_id'].to_dict()
    s2b = db1.sort_values('sw_score', ascending=False).drop_duplicates(subset=['sbj_id'])\
            .set_index('sbj_id')['qry_id'].to_dict()

    db1['best_sbj'] = db1['qry_id'].map(q2b)
    db1['best_qry'] = db1['sbj_id'].map(s2b)

    db1 = db1[(db1['sbj_id'] == db1['best_sbj']) & (db1['qry_id'] == db1['best_qry'])]
    del db1['best_sbj']
    del db1['best_qry']

    # Return your list
    return db1

class dnDs_command():
    '''
    This class just holds all of the mumbo-jumbo needed to calculate dn/ds
    '''

    def __init__(self):
        pass

    def verify(self):
        for att in ['Qfna', 'Qfaa', 'Rfna', 'Rfaa']:
            assert os.path.isfile(getattr(self, att))

def goANI_wrapper(cmd):
    '''
    Fromt that command object, call the real method
    '''
    print('running {0} vs {1}'.format(cmd.row['querry'], cmd.row['reference']))
    try:
        return cmd, goANI_dnds_calculation(cmd.Qfna, cmd.Qfaa, cmd.Rfna, cmd.Rfaa, cmd.geneDb)
    except Exception as e:
        print("whole file exception- {0}".format(str(e)))
        traceback.print_exc()
        return cmd, pd.DataFrame({'Failed':[True]})

def goANI_dnds_calculation(fna1, faa1, fna2, faa2, gedb, debug=False):
    '''
    This is a threadable command to determine the dn/ds of two genomes
    based on a list of genes

    Arguments:
        fna1 : .fna file of genome1
        faa1 : .faa file of genome1
        fna2 : .fna file of genome2
        faa2 : .faa file of genome2
        gedb : datatable listing the genes to align and calculate dn/ds for

    Returns:
        dndb : data-table containing raw dn/ds information
    '''
    # load .fasta files
    g1n = SeqIO.to_dict(SeqIO.parse(fna1, 'fasta', alphabet=IUPAC.IUPACUnambiguousDNA()))
    g1a = SeqIO.to_dict(SeqIO.parse(faa1, 'fasta', alphabet=IUPAC.protein))
    g2n = SeqIO.to_dict(SeqIO.parse(fna2, 'fasta', alphabet=IUPAC.IUPACUnambiguousDNA()))
    g2a = SeqIO.to_dict(SeqIO.parse(faa2, 'fasta', alphabet=IUPAC.protein))

    # set up aligner
    aligner = Align.PairwiseAligner()
    aligner.mode = 'global'
    #print(MatrixInfo.blosum62)
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    aligner.open_gap_score = -12
    aligner.extend_gap_score = -3

    # set up table
    table = defaultdict(list)

    # for every gene-pair to align
    j = 0
    for i, row in gedb.iterrows():
        try:
            # get the sequences
            a1 = g1a[row['qry_id']]
            if a1[-1] == '*':
                a1 = a1[:-1]
            a2 = g2a[row['sbj_id']]
            if a2[-1] == '*':
                a2 = a2[:-1]

            s1 = g1n[row['qry_id']]
            s2 = g2n[row['sbj_id']]

            # alingn them
            alignments = aligner.align(a1.seq, a2.seq)

            # Arbitrary cutoff to make sure this doesn't bug out
            if len(alignments) > 10000:
                print("Ahhh! {0} vs {1} has {2} alignments".format(row['qry_id'], row['sbj_id'],len(alignments)))
                raise Exception('Too many alignments exception')

            # convert to multi-sequence alignment
            alignment = min(alignments)
            ass = str(alignment).split('\n')
            msa = MultipleSeqAlignment([SeqRecord(Seq(ass[0], alphabet=IUPAC.protein)),
                                       SeqRecord(Seq(ass[-2], alphabet=IUPAC.protein))])

            # convert to codon alignment
            codon_aln = Bio.codonalign.build(msa, [s1, s2])

            # calculate dn/ds on the codon alignment
            dS, S, dN, N = custom_dn_ds(codon_aln._records[0], codon_aln._records[1])

            # save
            table['qry_id'].append(row['qry_id'])
            table['sbj_id'].append(row['sbj_id'])
            table['S_changed'].append(dS)
            table['S_sites'].append(S)
            table['N_changed'].append(dN)
            table['N_sites'].append(N)

            j += 1
            if debug:
                if j >= 10:
                    break

        except Exception as e:
            print("Alignment exception- {0}".format(e))
            table['qry_id'].append(row['qry_id'])
            table['sbj_id'].append(row['sbj_id'])
            table['S_changed'].append(0)
            table['S_sites'].append(0)
            table['N_changed'].append(0)
            table['N_sites'].append(0)

    dnDb = pd.DataFrame(table)

    return dnDb

def custom_dn_ds(codon_seq1, codon_seq2, method="NG86",
                codon_table=default_codon_table, k=1, cfreq=None):
    """
    http://biopython.org/DIST/docs/api/Bio.codonalign.codonseq-pysrc.html#cal_dn_ds

    Calculate dN and dS of the given two sequences.

    Available methods:
      - NG86  - `Nei and Gojobori (1986)`_ (PMID 3444411).
      - LWL85 - `Li et al. (1985)`_ (PMID 3916709).
      - ML    - `Goldman and Yang (1994)`_ (PMID 7968486).
      - YN00  - `Yang and Nielsen (2000)`_ (PMID 10666704).

    .. _`Nei and Gojobori (1986)`: http://www.ncbi.nlm.nih.gov/pubmed/3444411
    .. _`Li et al. (1985)`: http://www.ncbi.nlm.nih.gov/pubmed/3916709
    .. _`Goldman and Yang (1994)`: http://mbe.oxfordjournals.org/content/11/5/725
    .. _`Yang and Nielsen (2000)`: https://doi.org/10.1093/oxfordjournals.molbev.a026236

    Arguments:
    - codon_seq1 - CodonSeq or or SeqRecord that contains a CodonSeq
    - codon_seq2 - CodonSeq or or SeqRecord that contains a CodonSeq
    - w  - transition/transversion ratio
    - cfreq - Current codon frequency vector can only be specified
     when you are using ML method. Possible ways of
     getting cfreq are: F1x4, F3x4 and F61.

    """
    # Convert to codonseq
    if isinstance(codon_seq1, Bio.codonalign.codonseq.CodonSeq) \
            and isinstance(codon_seq2, Bio.codonalign.codonseq.CodonSeq):
        pass
    elif isinstance(codon_seq1, SeqRecord) and isinstance(codon_seq2, SeqRecord):
        codon_seq1 = codon_seq1.seq
        codon_seq2 = codon_seq2.seq
    else:
        raise TypeError("cal_dn_ds accepts two CodonSeq objects or SeqRecord "
                      "that contains CodonSeq as its seq!")
    if len(codon_seq1.get_full_rf_table()) != len(codon_seq2.get_full_rf_table()):
        raise RuntimeError("full_rf_table length of seq1 ({0}) and seq2 ({1}) "
                         "are not the same".format(
                             len(codon_seq1.get_full_rf_table()),
                             len(codon_seq2.get_full_rf_table()))
                         )
    if cfreq is None:
        cfreq = 'F3x4'
    elif cfreq is not None and method != 'ML':
        raise RuntimeError("cfreq can only be specified when you "
                         "are using ML method")
    if cfreq not in ('F1x4', 'F3x4', 'F61'):
        import warnings
        warnings.warn("Unknown cfreq ({0}). Only F1x4, F3x4 and F61 are "
                    "acceptable. Use F3x4 in the following.".format(cfreq))
        cfreq = 'F3x4'

    # get codon lists
    seq1_codon_lst = get_codon_list(codon_seq1)
    seq2_codon_lst = get_codon_list(codon_seq2)

    # filter for stop codons
    if ('TGA' in seq1_codon_lst) or ('TGA' in seq2_codon_lst):
        return 0,0,0,0

    # remove gaps in seq_codon_lst
    seq1 = []
    seq2 = []
    for i, j in zip(seq1_codon_lst, seq2_codon_lst):
        if ('-' not in i) and ('-' not in j):
            seq1.append(i)
            seq2.append(j)

    dS, S, dN, N = _ng86_custom(seq1, seq2, k, codon_table)
    return dS, S, dN, N

def _ng86_custom(seq1, seq2, k, codon_table):
    """
    NG86 method main function (PRIVATE).
    http://biopython.org/DIST/docs/api/Bio.codonalign.codonseq-pysrc.html#cal_dn_ds

    MATT - You changed this to return details instead of dN, dS.

    It now returns changed S, total S, changed N, total N

    """
    S_sites1, N_sites1 = count_sites(seq1,
                                            codon_table=codon_table, k=k)
    S_sites2, N_sites2 = count_sites(seq2,
                                            codon_table=codon_table, k=k)
    S_sites = (S_sites1 + S_sites2) / 2.0
    N_sites = (N_sites1 + N_sites2) / 2.0
    SN = [0, 0]
    for i, j in zip(seq1, seq2):
        SN = [m + n for m, n in zip(SN, _count_diff_NG86(i, j,
                                                           codon_table=codon_table))]
    #print(SN)
    ps = SN[0] / S_sites
    pn = SN[1] / N_sites
    if ps < 3 / 4:
        dS = abs(-3.0 / 4 * log(1 - 4.0 / 3 * ps))
    else:
        dS = -1
    if pn < 3 / 4:
        dN = abs(-3.0 / 4 * log(1 - 4.0 / 3 * pn))
    else:
        dN = -1
    return SN[0], S_sites, SN[1], N_sites


def get_codon_list(codonseq):
    """List of codons according to full_rf_table for counting (PRIVATE)."""
    full_rf_table = codonseq.get_full_rf_table()
    codon_lst = []
    for i, k in enumerate(full_rf_table):
        if isinstance(k, int):
            start = k
            try:
                end = int(full_rf_table[i+1])
            except IndexError:
                end = start+3
            this_codon = str(codonseq[start:end])
            if len(this_codon) == 3:
                codon_lst.append(this_codon)
            else:
                codon_lst.append(str(this_codon.ungap()))
        elif str(codonseq[int(k):int(k)+3]) == "---":
            codon_lst.append("---")
        else:
          # this may be problematic, as normally no codon shoud
          # fall into this condition
            codon_lst.append(codonseq[int(k):int(k)+3])
    return codon_lst


def count_sites(codon_lst, k=1, codon_table=default_codon_table):
    S_site = 0.0  # synonymous sites
    N_site = 0.0  # non-synonymous sites
    purine = ('A', 'G')
    pyrimidine = ('T', 'C')
    base_tuple = ('A', 'T', 'C', 'G')
    for codon in codon_lst:
        neighbor_codon = {'transition': [], 'transversion': []}
        # classify neighbor codons
        codon = codon.replace('U', 'T')
        if codon == '---':
            continue
        for n, i in enumerate(codon):
            for j in base_tuple:
                if i == j:
                    pass
                elif i in purine and j in purine:
                    codon_chars = [c for c in codon]
                    codon_chars[n] = j
                    this_codon = ''.join(codon_chars)
                    neighbor_codon['transition'].append(this_codon)
                elif i in pyrimidine and j in pyrimidine:
                    codon_chars = [c for c in codon]
                    codon_chars[n] = j
                    this_codon = ''.join(codon_chars)
                    neighbor_codon['transition'].append(this_codon)
                else:
                    codon_chars = [c for c in codon]
                    codon_chars[n] = j
                    this_codon = ''.join(codon_chars)
                    neighbor_codon['transversion'].append(this_codon)
        # count synonymous and non-synonymous sites
        #codon = codon.replace('T', 'U')
        if (codon == 'TAG'):
            print("STOP DETECTED")
            continue
        aa = codon_table.forward_table[codon]
        this_codon_N_site = this_codon_S_site = 0
        for neighbor in neighbor_codon['transition']:
            if neighbor in codon_table.stop_codons:
                this_codon_N_site += 1
            elif codon_table.forward_table[neighbor] == aa:
                this_codon_S_site += 1
            else:
                this_codon_N_site += 1
        for neighbor in neighbor_codon['transversion']:
            if neighbor in codon_table.stop_codons:
                this_codon_N_site += k
            elif codon_table.forward_table[neighbor] == aa:
                this_codon_S_site += k
            else:
                this_codon_N_site += k
        norm_const = (this_codon_N_site + this_codon_S_site)/3
        S_site += float(this_codon_S_site) / float(norm_const)
        N_site += float(this_codon_N_site) / float(norm_const)
    return (S_site, N_site)

def _count_diff_NG86(codon1, codon2, codon_table=default_codon_table):
    """Count differences between two codons (three-letter string; PRIVATE).

      The function will take multiple pathways from codon1 to codon2
      into account.
      """

    if not isinstance(codon1, str) or not isinstance(codon2, str):
        raise TypeError('_count_diff_NG86 accepts string object to represent codon ({0}, {1} detected)'.format(type(codon1),
                        type(codon2)))
    if len(codon1) != 3 or len(codon2) != 3:
        raise RuntimeError('codon should be three letter string ({0}, {1} detected)'.format(len(codon1),
                           len(codon2)))
    SN = [0, 0]  # synonymous and nonsynonymous counts
    if codon1 == '---' or codon2 == '---':
        return SN
    base_tuple = ('A', 'C', 'G', 'T')
    if not all(i in base_tuple for i in codon1):
        raise RuntimeError('Unrecognized character detected in codon1 {0} (Codons consist of A, T, C or G)'.format(codon1))
    if not all(i in base_tuple for i in codon2):
        raise RuntimeError('Unrecognized character detected in codon2 {0} (Codons consist of A, T, C or G)'.format(codon2))
    if codon1 == codon2:
        return SN
    else:
        diff_pos = []
        for (i, k) in enumerate(zip(codon1, codon2)):
            if k[0] != k[1]:
                diff_pos.append(i)

        def compare_codon(
            codon1,
            codon2,
            codon_table=default_codon_table,
            weight=1,
            ):
            """Compare two codon accounting for different pathways."""

            sd = nd = 0
            if len(set(map(codon_table.forward_table.get, [codon1,
                   codon2]))) == 1:
                # If both codons end up with the same aa
#                 for x in map(codon_table.forward_table.get, [codon1,
#                    codon2]):
#                     print(x)
#                 print(codon1, codon2)
                sd += weight
            else:
                nd += weight

            return (sd, nd)

        if len(diff_pos) == 1:
            SN = [i + j for (i, j) in zip(SN, compare_codon(codon1,
                  codon2, codon_table=codon_table))]
        elif len(diff_pos) == 2:
            codon2_aa = codon_table.forward_table[codon2]
            for i in diff_pos:
                temp_codon = codon1[:i] + codon2[i] + codon1[i + 1:]
                SN = [i + j for (i, j) in zip(SN, compare_codon(codon1,
                      temp_codon, codon_table=codon_table, weight=0.5))]
                SN = [i + j for (i, j) in zip(SN,
                      compare_codon(temp_codon, codon2,
                      codon_table=codon_table, weight=0.5))]
        elif len(diff_pos) == 3:
            codon2_aa = codon_table.forward_table[codon2]
            paths = list(permutations([0, 1, 2], 3))
            tmp_codon = []
            for p in paths:
                tmp1 = codon1[:p[0]] + codon2[p[0]] + codon1[p[0] + 1:]
                tmp2 = tmp1[:p[1]] + codon2[p[1]] + tmp1[p[1] + 1:]
                tmp_codon.append((tmp1, tmp2))
                SN = [i + j for (i, j) in zip(SN, compare_codon(codon1,
                      tmp1, codon_table, weight=0.5 / 3))]
                SN = [i + j for (i, j) in zip(SN, compare_codon(tmp1,
                      tmp2, codon_table, weight=0.5 / 3))]
                SN = [i + j for (i, j) in zip(SN, compare_codon(tmp2,
                      codon2, codon_table, weight=0.5 / 3))]
    return SN

class test_dnds_from_drep():
    '''
    Tests
    '''
    def setUp(self):
        self.test_dir = '/Users/mattolm/BioScripts/test_suite/test_backend/testdir/'

        if os.path.isdir(self.test_dir):
            shutil.rmtree(self.test_dir)
        os.mkdir(self.test_dir)

        shutil.copytree('/Users/mattolm/BioScripts/test_suite/test_files/goANI_test/',
                    '/Users/mattolm/BioScripts/test_suite/test_backend/testdir/goANI_test')
        self.wd = '/Users/mattolm/BioScripts/test_suite/test_backend/testdir/goANI_test/'

    def tearDown(self):
        if os.path.isdir(self.test_dir):
            shutil.rmtree(self.test_dir)

    def run(self):
        self.setUp()
        self.test0()
        self.tearDown()

        self.setUp()
        self.test1()
        self.tearDown()

        # self.setUp()
        # self.test3()
        # self.tearDown()

        print('all good')

    def test0(self):
        '''
        Test the split function
        '''
        Ndb = pd.read_csv(self.wd + 'data_tables/Ndb.csv')

        splits = []
        for ndb in iterate_Ndbs(Ndb, chunkSize=3):
            splits.append(ndb)
        nndb = pd.concat(splits)
        assert len(nndb) == len(nndb.drop_duplicates(subset=['reference', 'querry']))
        assert len(Ndb) == len(nndb), (len(Ndb), len(nndb))

        splits = []
        for ndb in iterate_Ndbs(Ndb, chunkSize=100):
            splits.append(ndb)
        nndb = pd.concat(splits)
        assert len(nndb) == len(nndb.drop_duplicates(subset=['reference', 'querry']))
        assert len(Ndb) == len(nndb), (len(Ndb), len(nndb))

    def test1(self):
        dnds_from_dRep(self.wd)

        # Load results
        Rdb = pd.read_csv(self.wd + 'data/dNdS_data/detailed_dNdS_info.csv')
        GDdb = pd.read_csv(self.wd + 'data/dNdS_data/genomeWide_dNdS_info.csv')

        failed = len(Rdb[Rdb['S_sites'] == 0])
        print("{0:.1f}% failed ({1} of {2})".format((failed / len(Rdb))*100,
                failed, len(Rdb)))

        assert len(Rdb) > 7000
        assert len(GDdb) == 3, len(GDdb)

    # def test3(self):
    #     dnds_from_dRep(self.wd)
    #
    #     # Load results
    #     Rdb = pd.read_csv(self.wd + 'data/dNdS_data/detailed_dNdS_info.csv')
    #
    #     failed = len(Rdb[Rdb['S_sites'] == 0])
    #     print("{0:.1f}% failed ({1} of {2})".format((failed / len(Rdb))*100,
    #             failed, len(Rdb)))
    #
    #     assert len(Rdb) > 7000


if __name__ == "__main__":
    """ This is executed when run from the command line """
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        usage=textwrap.dedent('''\
         \n
         dnds_from_drep.py combines many lists of scaffolds into a single,
         dereplicated list. \n\n
         -wd Takes a dRep work directory. THIS MUST HAVE BEEN GENERATED USING
         --S_algorithm goANI (available in v2.3.0). The command compare works just fine.

         -p is the number of processors

         Output will be made inside the input folder
         '''))

    parser.add_argument("-wd", "--drep_dir", help='fasta files to dereplicate on a per-scaffold level')
    parser.add_argument("-p", "--processes", help='number of processes to use', default=6, type=int)
    parser.add_argument('--debug', help='you probably dont want this- will only run on a subset of the data',
                        action='store_true')
    parser.add_argument('--test', help='run tests', action='store_true')

    # Specify output of "--version"
    parser.add_argument(
        "--version",
        action="version",
        version="%(prog)s (version {version})".format(version=__version__))

    args = parser.parse_args()

    if args.test:
        test_dnds_from_drep().run()
    else:
        main(args)
