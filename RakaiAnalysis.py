__author__ = 'jayna'

import numpy as np
import os

from matplotlib import pyplot as plt
from Bio import AlignIO
from Bio import SeqIO
from Bio.Align import AlignInfo
from Bio.Alphabet import generic_dna

import pandas as pd
import statsmodels.formula.api as smf
from Bio.Seq import Seq

import matplotlib as mpl
import matplotlib.cm as cm

import random
import array
import re


def get_consensus_from_file(fastaout):
    aln = AlignIO.read('%s' % fastaout, 'fasta')
    return aln[0]


def atoi(text):
    return int(text) if text.isdigit() else text

def natural_keys(text):
    '''
    alist.sort(key=natural_keys) sorts in human order
    http://nedbatchelder.com/blog/200712/human_sorting.html
    (See Toothy's implementation in the comments)
    '''
    return [atoi(c) for c in re.split('(\d+)', text)]


def number_of_N_and_S_sites(fastain, sites):
    fastaseq = AlignIO.read('%s' % fastain, 'fasta')

    sequence = AlignInfo.SummaryInfo(fastaseq).dumb_consensus().upper()

    print sequence
    # print len(sequence)


    site1 = 0
    site2 = 0

    bases = ['A', 'T', 'G', 'C']
    if sites == None:

        site1 = 0
        site2 = len(sequence)

    else:

        site1 = sites[0] - 1
        site2 = sites[1] - 1

    total_nonsyn = 0
    total_syn = 0

    for i in xrange(site1, site2):

        non_syn = 0
        syn = 0

        codon = sequence[0:3]
        codon_pos = 0

        if i % 3 == 0:

            codon = sequence[i:i + 3]
            codon_pos = 0

        elif i % 3 == 1:

            codon = sequence[i - 1:(i - 1) + 3]
            codon_pos = 1

        elif i % 3 == 2:

            codon = sequence[(i - 2):(i - 2) + 3]
            codon_pos == 3

        #print i, codon

        for b in bases:

            if("X" in str(codon)):
                break
            codon_string = list(codon)
            codon_string[codon_pos] = b

            new_codon_seq = Seq("".join(codon_string))

            new_aa = new_codon_seq.translate()
            aa = codon.translate()

            if str(new_aa) == str(aa):

                syn += 1
            else:
                non_syn += 1

        total_nonsyn += float(non_syn / 4.0)
        total_syn += float(syn / 4.0)
        # print float(non_syn/4.0), float(syn/4.0)

    return total_nonsyn, total_syn


def split(fastain, n):

    file = SeqIO.parse('%s'%fastain, 'fasta')
    dates = []
    for each in file:

        seq = each.seq
        name = each.name
        parts = name.split('_')
        date = parts[len(parts)-n] # when date is at the end of the sequence name
        #print int(re.match(r'\d+', parts[0]))
        #print parts[0]
        #date_numbers = map(int, re.findall(r'\d+', parts[0]))
        #print date_numbers[1]

        dates.append(date)
    unique_dates = set(dates)
    #print unique_dates

    # date_list = []
    # for d in unique_dates:
    #     print d
    #
    #     if d.isdigit() == True:
    #         date_list.append(float(d))
    #
    # date_list.sort()
    # print str(date_list)+","+str(len(date_list))

    all_timepoints = {}

    for d in unique_dates:

        all_timepoints[str(d)] = []

    file = SeqIO.parse('%s'%fastain, 'fasta')

    # dividing the alignment by different timepojnts
    for each in file:
        seq = each.seq
        name = each.name
        parts = name.split('_')
        # date_numbers = map(int, re.findall(r'\d+', parts[0]))
        # date = date_numbers[1]
        #if parts[len(parts)-1].isdigit():
        date = str(parts[len(parts)-n])

        for d in unique_dates:
            if str(d) in str(date):
                 all_timepoints[str(d)].append(each)


    sorted_keys = all_timepoints.keys()
    sorted_keys.sort(key=natural_keys)
    return all_timepoints


#For HIV Rakai cohort
def divergence(fastain, translate, date_part, patient_id, sites):
    #fasta = open('%s' % filename, 'r')

    seqs_by_timepoint = split(fastain, date_part)

    # conseq = consensus.seq[(sites_pos[0]-1):(sites_pos[1]-1)]
    # conseq = Seq(str(consensus).replace('-','N'))
    # consensus = Seq(conseq.seq.tostring().replace('-','N'))

    # seq_length = len(consensus)
    mean_divergence = []
    median_divergence = []

    lower_divergence_25 = []
    upper_divergence_75 = []
    lower_divergence_5 = []
    upper_divergence_95 = []
    divergence_std = []
    mean_N_divergence = []
    median_N_divergence = []

    lower_N_divergence_25 = []
    upper_N_divergence_75 = []
    lower_N_divergence_5 = []
    upper_N_divergence_95 = []
    N_divergence_std = []
    mean_S_divergence = []
    median_S_divergence = []
    lower_S_divergence_25 = []
    upper_S_divergence_75 = []
    lower_S_divergence_5 = []
    upper_S_divergence_95 = []
    S_divergence_std = []
    dN = []
    dN_med = []
    dN_lower_25 = []
    dN_upper_75 = []
    dN_lower_5 = []
    dN_upper_95 = []
    dN_std = []
    dS = []
    dS_med = []
    dS_lower_25 = []
    dS_upper_75 = []
    dS_lower_5 = []
    dS_upper_95 = []
    dS_std = []
    patient = []

    # parts = str.split(fastain, "/")
    # parts2 = str.split(parts[len(parts)-1], "_")


    patient.append(patient_id)

    nonsyn_sites, syn_sites = number_of_N_and_S_sites(fastain, None)
    print nonsyn_sites, syn_sites

    sorted_timepoints = seqs_by_timepoint.keys()
    sorted_timepoints.sort(key=natural_keys)

    print sorted_timepoints
    first_timepoint = AlignIO.MultipleSeqAlignment(seqs_by_timepoint[sorted_timepoints[0]])

    consensus = AlignInfo.SummaryInfo(first_timepoint).dumb_consensus(threshold=0.01).upper()

    sampleTimes = []
    for t in sorted_timepoints:
        sampleTimes.append(float(t))


    #for f in filelist:
    for t in range(0,len(sorted_timepoints)):

        divergence = []
        divergence_N = []
        divergence_S = []
        divergence_dN = []
        divergence_dS = []
        # diff = 0


        seqs_at_t = seqs_by_timepoint[sorted_timepoints[t]]

        for each in seqs_at_t:

            parts = str.split(each.name, "_")
            freq = 1#int(parts[2].strip())
            diff = 0
            diff_N = 0
            diff_S = 0

            seq = Seq(str(each.seq).upper().replace('-', 'N'))[sites[0]:sites[1]]

            codon_pos_start = 0
            codon_pos_end = 2

            A_i = str(seq).find('A')
            T_i = str(seq).find('T')
            G_i = str(seq).find('G')
            C_i = str(seq).find('C')

            start = [A_i, T_i, G_i, C_i]

            A_ii = str(seq).rfind('A')
            T_ii = str(seq).rfind('T')
            G_ii = str(seq).rfind('G')
            C_ii = str(seq).rfind('C')

            end = [A_ii, T_ii, G_ii, C_ii]

            start_i = min(start)
            end_i = max(end)

            if start_i > -1 and end_i > -1:
                # print start_i, end_i

                remainder_1 = start_i % 3
                remainder_2 = end_i % 3

                if remainder_1 != 0:
                    b = remainder_1 != codon_pos_start
                    # print start_i, start_i + (3-remainder_1)
                    start_i = start_i + (3 - remainder_1)

                if remainder_2 != 2:
                    # tprint end_i, end_i + (3-remainder_2)
                    end_i = end_i + (2 - remainder_2)

                seq = seq[start_i: end_i + 1]
                gaps = str(seq).count('N')

                seq_length = len(seq)
                aa_length = seq_length / 3

                #if ((len(seq) % 3) > 0):
                    #print 'x', start_i, end_i, (len(seq) % 3), remainder_2, len(seq)

                # print start_i, end_i

                # print consensus.seq[start_i: end_i].translate()

                #conseq = consensus.seq[(start_i) + sites_pos[0] - 1: (start_i + sites_pos[0] - 1 + length)]

                conseq = Seq(str(consensus).replace('X', 'N'))[sites[0]:sites[1]]



                translated_seq = seq.translate()

                gaps_con = str(conseq).count('N')

                if gaps_con == seq_length:

                    print ("all gaps in conseq")
                    break

                else:
                    # if (seq_length >= length and (float(gaps) / float(seq_length)) < 0.05 and
                    #             (float(gaps_con) / float(seq_length)) < 0.05):

                        # print translated_seq, conseq.translate(),
                        count = 0
                        if (translate):
                            seq = each.seq.translate()

                        # count +=1

                        for a in range(seq_length):

                            i = a

                            if (str(conseq[i]) != "N"):

                                if (str(seq[i]) != "N"):

                                    count = count + 1

                                    if (conseq[i] != seq[i]):

                                        codon = []

                                        if (i % 3 == 0):
                                            cp = i
                                            cp_a = i + 1
                                            cp_b = i + 2

                                            codon = [cp, cp_a, cp_b]


                                        elif (i % 3 == 1):
                                            cp_a = i - 1
                                            cp = i
                                            cp_b = i + 1

                                            codon = [cp_a, cp, cp_b]

                                        else:

                                            cp_a = i - 2
                                            cp_b = i - 1
                                            cp = i

                                            codon = [cp_a, cp_b, cp]



                                        consensus_aa = conseq[codon[0]:(codon[2] + 1)].translate()
                                        current_aa = seq[codon[0]:(codon[2] + 1)].translate()

                                        #print(str(consensus_aa), str(current_aa))
                                        if 'X' in conseq[codon[0]:(codon[2] + 1)]:
                                            break

                                        if (str(consensus_aa) != str(current_aa)):

                                            diff_N += 1
                                        else:
                                            diff_S += 1

                                        diff += 1

                        # print diff/count, diff/seq_length


                        divergence.extend([float(diff) / float(count)]*freq)
                        divergence_N.extend([float(diff_N) / float(count)]*freq)
                        divergence_S.extend([float(diff_S) / float(count)]*freq)
                        divergence_dN.extend([float(diff_N) / float(nonsyn_sites)]*freq)
                        divergence_dS.extend([float(diff_S) / float(syn_sites)]*freq)

        #print(t, len(divergence))

        if len(divergence) < 100:

            mean_divergence.append(float('nan'))
            median_divergence.append(float('nan'))
            lower_divergence_25.append(float('nan'))
            upper_divergence_75.append(float('nan'))
            lower_divergence_5.append(float('nan'))
            upper_divergence_95.append(float('nan'))
            divergence_std.append(float(1000))

            mean_N_divergence.append(float('nan'))
            median_N_divergence.append(float('nan'))
            lower_N_divergence_25.append(float('nan'))
            upper_N_divergence_75.append(float('nan'))
            lower_N_divergence_5.append(float('nan'))
            upper_N_divergence_95.append(float('nan'))
            N_divergence_std.append(float(1000))

            mean_S_divergence.append(float('nan'))
            median_S_divergence.append(float('nan'))
            lower_S_divergence_25.append(float('nan'))
            upper_S_divergence_75.append(float('nan'))
            lower_S_divergence_5.append(float('nan'))
            upper_S_divergence_95.append(float('nan'))
            S_divergence_std.append(float(1000))

            dN.append(float('nan'))
            dN_med.append(float('nan'))
            dN_lower_25.append(float('nan'))
            dN_upper_75.append(float('nan'))
            dN_lower_5.append(float('nan'))
            dN_upper_95.append(float('nan'))
            dN_std.append(float('nan'))

            dS.append(float('nan'))
            dS_med.append(float('nan'))
            dS_lower_25.append(float('nan'))
            dS_upper_75.append(float('nan'))
            dS_lower_5.append(float('nan'))
            dS_upper_95.append(float('nan'))
            dS_std.append(float(1000))


        else:

            #print divergence
            mean_divergence.append(np.mean(divergence))
            median_divergence.append(np.percentile(divergence, 50))
            lower_divergence_25.append(np.percentile(divergence, 25))
            upper_divergence_75.append(np.percentile(divergence, 75))
            lower_divergence_5.append(np.percentile(divergence, 5))
            upper_divergence_95.append(np.percentile(divergence, 95))
            divergence_std.append(np.std(divergence))

            mean_N_divergence.append(np.mean(divergence_N))
            median_N_divergence.append(np.percentile(divergence_N, 50))
            lower_N_divergence_25.append(np.percentile(divergence_N, 25))
            upper_N_divergence_75.append(np.percentile(divergence_N, 75))
            lower_N_divergence_5.append(np.percentile(divergence_N, 5))
            upper_N_divergence_95.append(np.percentile(divergence_N, 95))
            N_divergence_std.append(np.std(divergence_N))

            mean_S_divergence.append(np.mean(divergence_S))
            median_S_divergence.append(np.percentile(divergence_S, 50))
            lower_S_divergence_25.append(np.percentile(divergence_S, 25))
            upper_S_divergence_75.append(np.percentile(divergence_S, 75))
            lower_S_divergence_5.append(np.percentile(divergence_S, 5))
            upper_S_divergence_95.append(np.percentile(divergence_S, 95))
            S_divergence_std.append(np.std(divergence_S))

            dN.append(np.mean(divergence_dN))
            dN_med.append(np.percentile(divergence_dN, 50))
            dN_lower_25.append(np.percentile(divergence_dN, 25))
            dN_upper_75.append(np.percentile(divergence_dN, 75))
            dN_lower_5.append(np.percentile(divergence_dN, 5))
            dN_upper_95.append(np.percentile(divergence_dN, 95))
            dN_std.append(np.std(divergence_dN))

            dS.append(np.mean(divergence_dS))
            dS_med.append(np.percentile(divergence_dS, 50))
            dS_lower_25.append(np.percentile(divergence_dS, 25))
            dS_upper_75.append(np.percentile(divergence_dS, 75))
            dS_lower_5.append(np.percentile(divergence_dS, 5))
            dS_upper_95.append(np.percentile(divergence_dS, 95))
            dS_std.append(np.std(divergence_dS))

    # if (type(median_divergence[0]) == type(float('nan'))) == True:
    #     mean_divergence[0] = 0.0
    #     lower_divergence_5[0] = 0.0
    #     lower_divergence_25[0] = 0.0
    #     upper_divergence_95[0] = 0.0
    #     upper_divergence_75[0] = 0.0
    #     mean_S_divergence[0] = 0.0
    #     lower_S_divergence_5[0] = 0.0
    #     lower_S_divergence_25[0] = 0.0
    #     upper_S_divergence_95[0] = 0.0
    #     upper_S_divergence_75[0] = 0.0
    #     mean_N_divergence[0] = 0.0
    #     lower_N_divergence_5[0] = 0.0
    #     lower_N_divergence_25[0] = 0.0
    #     upper_N_divergence_95[0] = 0.0
    #     upper_N_divergence_75[0] = 0.0
    #     divergence_std[0] = 0.0000000001
    #     N_divergence_std[0] = 0.0000000001
    #     S_divergence_std[0] = 0.0000000001
    #
    # if (type(median_divergence[1]) == type(float('nan'))) == True:
    #     mean_divergence[1] = 0.0
    #     lower_divergence_5[1] = 0.0
    #     lower_divergence_25[1] = 0.0
    #     upper_divergence_95[1] = 0.0
    #     upper_divergence_75[1] = 0.0
    #     mean_S_divergence[1] = 0.0
    #     lower_S_divergence_5[1] = 0.0
    #     lower_S_divergence_25[1] = 0.0
    #     upper_S_divergence_95[1] = 0.0
    #     upper_S_divergence_75[1] = 0.0
    #     mean_N_divergence[1] = 0.0
    #     lower_N_divergence_5[1] = 0.0
    #     lower_N_divergence_25[1] = 0.0
    #     upper_N_divergence_95[1] = 0.0
    #     upper_N_divergence_75[1] = 0.0
    #     divergence_std[1] = 1000
    #     N_divergence_std[1] = 1000
    #     S_divergence_std[1] = 1000


    df = pd.DataFrame.from_items(
        [('Times', sampleTimes), ('Divergence_median', median_divergence), ('Divergence_mean', mean_divergence),
         ('Divergence_N_med', median_N_divergence), ('Divergence_N_mean', mean_N_divergence),
         ('Divergence_S_med', median_S_divergence), ('Divergence_S_mean', mean_S_divergence),
         ('Divergence_lower_25', lower_divergence_25), ('Divergence_upper_75', upper_divergence_75),
         ('Divergence_lower_5', lower_divergence_5), ('Divergence_upper_95', upper_divergence_95),
         ("Divergence_std", divergence_std),
         ('Divergence_N_lower_25', lower_N_divergence_25), ('Divergence_N_upper_75', upper_N_divergence_75),
         ('Divergence_N_lower_5', lower_N_divergence_5), ('Divergence_N_upper_95', upper_N_divergence_95),
         ("Divergence_N_std", N_divergence_std),
         ('Divergence_S_lower_25', lower_S_divergence_25), ('Divergence_S_upper_75', upper_S_divergence_75),
         ('Divergence_S_lower_5', lower_S_divergence_5), ('Divergence_S_upper_95', upper_S_divergence_95),
         ("Divergence_S_std", S_divergence_std), ('Patients', patient*len(sampleTimes))])

    # csvfilename = filename.replace("filelist_", "divergence_results_by_birth_patristic_sites_"+str(sites_pos[0])+"_to_"+str(sites_pos[1])+"_")
    # csvfilename = csvfilename.replace(".txt",".csv")
    # df.to_csv(csvfilename)

    df = df.dropna()

    # print sampleTimes - (np.ones(len(sampleTimes))*sampleTimes[0])


    # print mean_divergence
    # print divergence_std

    if len(df) > 1:

        lm1 = smf.wls(formula='Divergence_mean ~ Times-1', data=df,
                      missing='drop').fit()
        lm2 = smf.wls(formula='Divergence_N_mean ~ Times-1', data=df,
                      missing='drop').fit()
        lm3 = smf.wls(formula='Divergence_S_mean ~ Times-1', data=df,
                      missing='drop').fit()

        lm4 = smf.wls(formula='Divergence_lower_5 ~ Times-1', data=df,
                      missing='drop').fit()
        lm5 = smf.wls(formula='Divergence_N_lower_5 ~ Times-1', data=df,
                      missing='drop').fit()
        lm6 = smf.wls(formula='Divergence_S_lower_5 ~ Times-1', data=df,
                      missing='drop').fit()

        lm7 = smf.wls(formula='Divergence_upper_95 ~ Times-1', data=df,
                      missing='drop').fit()
        lm8 = smf.wls(formula='Divergence_N_upper_95 ~ Times-1', data=df,
                      missing='drop').fit()
        lm9 = smf.wls(formula='Divergence_S_upper_95 ~ Times-1', data=df,
                      missing='drop').fit()

        lm10 = smf.wls(formula='Divergence_lower_25 ~ Times-1', data=df,
                       missing='drop').fit()
        lm11 = smf.wls(formula='Divergence_N_lower_25 ~ Times-1', data=df,
                       missing='drop').fit()
        lm12 = smf.wls(formula='Divergence_S_lower_25 ~ Times-1', data=df,
                       missing='drop').fit()

        lm13 = smf.wls(formula='Divergence_upper_75 ~ Times-1', data=df,
                       missing='drop').fit()
        lm14 = smf.wls(formula='Divergence_N_upper_75 ~ Times-1', data=df,
                       missing='drop').fit()
        lm15 = smf.wls(formula='Divergence_S_upper_75 ~ Times-1', data=df,
                       missing='drop').fit()

        total_div = ['total', lm1.params[0], lm4.params[0], lm7.params[0], lm10.params[0], lm13.params[0]]
        N_div = ['N', lm2.params[0], lm5.params[0], lm8.params[0], lm11.params[0], lm14.params[0]]
        S_div = ['S', lm3.params[0], lm6.params[0], lm9.params[0], lm12.params[0], lm15.params[0]]
    else:

        total_div = ['total', float('nan'), float('nan'), float('nan'), float('nan'), float('nan')]
        N_div = ['N', float('nan'), float('nan'), float('nan'), float('nan'), float('nan')]
        S_div = ['S', float('nan'), float('nan'), float('nan'), float('nan'), float('nan')]
    # dN_div = ['dN', dN_mean[0], dN_5[0], dN_95[0], dN_25[0], dN_75[0]]
    # dS_div = ['dS', dS_mean[0], dS_5[0], dS_95[0], dS_25[0], dS_75[0]]  # print('lower_5')

    return df, total_div, N_div, S_div, nonsyn_sites, syn_sites

def pwd_per_timepoint(fastain):

    seqs_by_timepoint = split(fastain)

    mean_diversity = []


    patient=[]

    parts = str.split(fastain, "/")
    parts2 = str.split(parts[len(parts) - 1], "_")

    patient.append(parts2[0])

    nonsyn_sites, syn_sites = number_of_N_and_S_sites(fastain, None)
    print nonsyn_sites, syn_sites

    sorted_timepoints = seqs_by_timepoint.keys()
    sorted_timepoints.sort(key=natural_keys)

    sampleTimes = []
    for t in sorted_timepoints:
        sampleTimes.append(int(t))

    for t in range(0,len(sorted_timepoints)):

        diversity = []

        seqs_at_t = seqs_by_timepoint[sorted_timepoints[t]]


        for i in range(0,len(seqs_at_t)-1):

            for j in range(1, len(seqs_at_t)):
                diff = 0
                seq1 = seqs_at_t[i]
                seq2 = seqs_at_t[j]

                parts1 = str.split(seq1.name, "_")
                freq1 = int(parts1[2].strip())

                parts2 = str.split(seq2.name, "_")
                freq2 = int(parts2[2].strip())



                for a in range(len(seq1.seq)):

                    if (str(seq1.seq)[a] != '-'):

                        if (str(seq2.seq)[a] != '-'):

                            if (str(seq1.seq[a]).upper() != str(seq2.seq[a]).upper()):
                                #print str(seq1.seq[a]).upper(), str(seq2.seq[a]).upper()
                                diff += 1

                pwd_i = (float(diff) / len(seq1.seq))

                diversity.extend([pwd_i]*freq1*freq2)
                #print t, i, j, pwd_i, freq1, freq2, pwd_i*freq1*freq2/(freq1+freq2), len(seq1.seq)
                    # pwd.append(pwd_i)
                    # raw_diff.append(float(diff) / float(len(sub_aln)))


        print sorted_timepoints[t], np.mean(diversity)
        mean_diversity.append(np.mean(diversity))

    df = pd.DataFrame.from_items(
        [('Times', sampleTimes), ('Diversity_mean', mean_diversity),('Patients', patient*len(sampleTimes))])

    return df
