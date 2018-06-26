__author__ = 'jayna'

import numpy as np
from Bio import AlignIO
from Bio import SeqIO
from Bio.Align import AlignInfo

import pandas as pd
from Bio.Seq import Seq
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

        # print i, codon

        for b in bases:

            if ("X" in str(codon)):
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
    file = SeqIO.parse('%s' % fastain, 'fasta')
    dates = []
    total_seq = 0
    for each in file:
        seq = each.seq
        name = each.name
        parts = name.split('_')
        date = parts[len(parts) - n]  # when date is at the end of the sequence name
        freq = int(parts[2].strip())
        total_seq += freq


        dates.append(date)

    unique_dates = set(dates)

    all_timepoints = {}

    for d in unique_dates:
        all_timepoints[str(d)] = []

    file = SeqIO.parse('%s' % fastain, 'fasta')

    # dividing the alignment by different timepojnts
    for each in file:
        seq = each.seq
        name = each.name
        parts = name.split('_')
        # date_numbers = map(int, re.findall(r'\d+', parts[0]))
        # date = date_numbers[1]
        # if parts[len(parts)-1].isdigit():
        date = str(parts[len(parts) - n])

        for d in unique_dates:
            if str(d) in str(date):
                all_timepoints[str(d)].append(each)

    sorted_keys = all_timepoints.keys()
    sorted_keys.sort(key=natural_keys)
    return all_timepoints, total_seq


# This method estimates the mean divergence - subs per nucleotide.
# Key difference to "divergenceFromFounders" method is that it
# assumes a single founder sequence (consensus sequence of the first
# time point). As a consequence, divergence estimates from this
# method are expected to be higher when compared to estimates
# from "divergenceFromFounders" method

def divergence(fastain, patient_id, cutoff):
    # fasta = open('%s' % filename, 'r')

    split_fasta = split(fastain, 1)
    seqs_by_timepoint = split_fasta[0]
    total_seq = split_fasta[1]

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

    sorted_timepoints = seqs_by_timepoint.keys()
    sorted_timepoints.sort(key=natural_keys)

    print sorted_timepoints
    first_timepoint = AlignIO.MultipleSeqAlignment(seqs_by_timepoint[sorted_timepoints[0]])

    consensus = AlignInfo.SummaryInfo(first_timepoint).dumb_consensus(threshold=0.01).upper()
    conseq = Seq(str(consensus).replace('X','N'))

    prot = ""
    if "gag" in fastain:
        prot = "gag"
    else:
        prot = "gp41"


    sampleTimes = []
    for t in sorted_timepoints:
        sampleTimes.append(float(t))

    # for f in filelist:
    for t in range(0, len(sorted_timepoints)):

        divergence = []
        divergence_N = []
        divergence_S = []
        divergence_dN = []
        divergence_dS = []
        # diff = 0

        seqs_at_t = seqs_by_timepoint[sorted_timepoints[t]]

        seq_length = len(seqs_at_t[0].seq)

        seq_freq = get_seq_freq(seqs_at_t)


        seqs_at_t_array = np.asarray(seqs_at_t)



        # i want to calculate derived freq wrt to consequence not minor freq per site
        #for c in xrange(0,len(consensus_seqs)):

        full_der_freq = []

        total_site_freq = []

        for i in range(seq_length):

            site_a = seqs_at_t_array[:, i]

            anc_freq = 0
            der_freq = 0

            #gap_count = "".join(site_a).count('-')


            for j in range(0, len(seq_freq)):

                if site_a[j] != '-':
                    if conseq[i].lower()==site_a[j]:
                        anc_freq += seq_freq[j]
                    else:
                        der_freq += seq_freq[j]


                # if (site_a[j] == 'a'):
                #     A += seq_freq[j]
                # elif (site_a[j] == 'c'):
                #     C += seq_freq[j]
                # elif (site_a[j] == 't'):
                #     T += seq_freq[j]
                # elif (site_a[j] == 'g'):
                #     G += seq_freq[j]

            total_seq = sum([der_freq, anc_freq])

            full_der_freq.append(der_freq)

            total_site_freq.append(total_seq)

            #print [der_freq, anc_freq], total_seq
            #total_site_freq_per_consensus.append(total_site_freq)
            #full_der_freq_per_consensus.append(full_der_freq)



        #for c in xrange(0, len(consensus_seqs)):
        for i in range(seq_length):

            # print i, full_der_freq[i], patient_id, sorted_timepoints[t], total_seq, float(
            #     full_der_freq[i]) / float(total_seq)
            diff = 0
            diff_N = 0
            diff_S = 0
            count = total_site_freq[i]
            count1 = 0
            if full_der_freq[i] > cutoff * total_seq:


                for each in seqs_at_t:

                    parts = str.split(each.name, "_")
                    freq = int(parts[2].strip())


                    seq = Seq(str(each.seq).upper().replace('-', 'N'))

                    if (str(conseq[i]) != "N"):

                        if (str(seq[i]) != "N"):


                            count1+=freq

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

                                # print(str(consensus_aa), str(current_aa))
                                if 'X' in conseq[codon[0]:(codon[2] + 1)]:
                                    break

                                if (str(consensus_aa) != str(current_aa)):

                                    diff_N += freq
                                else:
                                    diff_S += freq


                                #print i, current_aa, consensus_aa, diff_N, diff_S, each.name, freq
                                diff += freq

                        #print each.name, sorted_timepoints[t], "d", float(diff), i, seq_length, count

            print(count, count1, i, diff, diff_N, diff_S)
            #
            # if((count-count1) != 0):
            #     print(count, count1, i, diff, diff_N, diff_S)

            if count > 0:

                #print i, patient_id, diff, count
                divergence.extend([float(diff)/float(count)])
                divergence_N.extend([float(diff_N)/float(count)])
                divergence_S.extend([float(diff_S)/float(count)])
                divergence_dN.extend([float(diff_N)/float(nonsyn_sites)/float(count)])
                divergence_dS.extend([float(diff_S)/float(syn_sites)/float(count)])


        if len(divergence)>1:
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

            if ("gag" in fastain):
                csvfile_gag_b.write(patient_id + "," + str(sorted_timepoints[t]) + "," +
                                        str(np.mean(divergence)) + "," + str(np.percentile(divergence, 50)) + "," +
                                        str(np.percentile(divergence, 5)) + "," + str(np.percentile(divergence, 95)) + "," +
                                        str(np.mean(divergence_N)) + "," + str(np.percentile(divergence_N, 50)) + "," +
                                        str(np.percentile(divergence_N, 5)) + "," + str(np.percentile(divergence_N, 95)) + "," +
                                        str(np.mean(divergence_S)) + "," + str(np.percentile(divergence_S, 50)) + "," +
                                        str(np.percentile(divergence_S, 5)) + "," + str(np.percentile(divergence_S, 95)) + "\n")

                csvfile_gag_b.flush()

            elif ("gp41" in fastain):
                csvfile_gp41_b.write(patient_id + "," + str(sorted_timepoints[t]) + "," +
                                        str(np.mean(divergence)) + "," + str(np.percentile(divergence, 50)) + "," +
                                        str(np.percentile(divergence, 5)) + "," + str(np.percentile(divergence, 95)) + "," +
                                        str(np.mean(divergence_N)) + "," + str(np.percentile(divergence_N, 50)) + "," +
                                        str(np.percentile(divergence_N, 5)) + "," + str(np.percentile(divergence_N, 95)) + "," +
                                        str(np.mean(divergence_S)) + "," + str(np.percentile(divergence_S, 50)) + "," +
                                        str(np.percentile(divergence_S, 5)) + "," + str(np.percentile(divergence_S, 95)) + "\n")


        else:
            print "xxx", patient_id, sorted_timepoints[t]

        print patient_id, sorted_timepoints[t], len(divergence)


# Method to estimate divergence - used in the main paper
# Multiple founders are possible in this analysis, which have been inferred using BAPS.
# The method calculates the minimum distance to the founder strain for each sequence,
# and returns the average divergence in units of substitutions (from the founder) per nucleotide site
def divergenceFromFounders(fastain, patient_id, cutoff):
    # fasta = open('%s' % filename, 'r')

    split_fasta = split(fastain, 1)
    seqs_by_timepoint = split_fasta[0]
    total_seq = split_fasta[1]

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

    patient = []

    # parts = str.split(fastain, "/")
    # parts2 = str.split(parts[len(parts)-1], "_")


    patient.append(patient_id)

    #nonsyn_sites, syn_sites = number_of_N_and_S_sites(fastain, None)

    sorted_timepoints = seqs_by_timepoint.keys()
    sorted_timepoints.sort(key=natural_keys)

    print sorted_timepoints
    #first_timepoint = AlignIO.MultipleSeqAlignment(seqs_by_timepoint[sorted_timepoints[0]])
    # consensus = AlignInfo.SummaryInfo(first_timepoint).dumb_consensus(threshold=0.01).upper()
    # conseq = Seq(str(consensus).replace('X','N'))

    prot = ""
    if "gag" in fastain:
        prot = "gag"
    else:
        prot = "gp41"

    consensus_seq_file = "~/consensus_sequences_1st_timepoint/"+p+"_"+prot+"_firsttimepoint_consensus.fasta"
    consensus_seqs = AlignIO.read(consensus_seq_file, 'fasta')

    sampleTimes = []
    for t in sorted_timepoints:
        sampleTimes.append(float(t))

    # for f in filelist:
    for t in range(0, len(sorted_timepoints)):

        divergence = []
        divergence_N = []
        divergence_S = []


        seqs_at_t = seqs_by_timepoint[sorted_timepoints[t]]

        seq_length = len(seqs_at_t[0].seq)

        seq_freq = get_seq_freq(seqs_at_t)


        seqs_at_t_array = np.asarray(seqs_at_t)


        full_der_freq_per_consensus = []

        total_site_freq_per_consensus = []

        #conseq = consensus_seqs[c]


        for c in range(len(consensus_seqs)):

            full_der_freq = []

            total_site_freq = []

            conseq = consensus_seqs[c]

            for i in range(seq_length):

                site_a = seqs_at_t_array[:, i]

                anc_freq = 0
                der_freq = 0

                #gap_count = "".join(site_a).count('-')


                for j in range(0, len(seq_freq)): #comparing each site from each sequence to consensus site (nucleotide comparison)


                    if site_a[j] != '-':

                        if conseq[i]==site_a[j]:
                            anc_freq += seq_freq[j]
                        else:
                            der_freq += seq_freq[j]


                total_seq = sum([der_freq, anc_freq])

                full_der_freq.append(der_freq)

                total_site_freq.append(total_seq)

            #print [der_freq, anc_freq], total_seq
            total_site_freq_per_consensus.append(total_site_freq)
            full_der_freq_per_consensus.append(full_der_freq)


        for each in seqs_at_t:

            list_diffN = []
            list_diffS = []
            list_diff = []

            parts = str.split(each.name, "_")
            freq = int(parts[2].strip())

            seq = Seq(str(each.seq).upper().replace('-', 'N'))

            for c in range(len(consensus_seqs)):


                conseq = consensus_seqs[c].seq

                conseq = Seq(str(conseq).upper().replace('-', 'N'))

                full_der_freq = full_der_freq_per_consensus[c]
                total_site_freq = total_site_freq_per_consensus[c]

                diff_N = 0
                diff_S = 0
                diff = 0

                for i in range(seq_length):


                    if full_der_freq[i] > cutoff*total_site_freq[i]:

                        if (str(conseq[i]) != "N"):

                            if (str(seq[i]) != "N"):


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

                                    # print(str(consensus_aa), str(current_aa))
                                    if 'X' in conseq[codon[0]:(codon[2] + 1)]:
                                        break

                                    if (str(consensus_aa) != str(current_aa)):

                                        diff_N += 1
                                    else:
                                        diff_S += 1
                                        #print(str(c)+" site "+str(i)+" "+str(conseq[codon[0]:(codon[2] + 1)])+" "+str(seq[codon[0]:(codon[2] + 1)]))

                                    # print i, current_aa, consensus_aa, diff_N, diff_S, each.name, freq
                                    diff += 1

                    print c, each.name, sorted_timepoints[t], "d", float(diff), i, seq_length, total_site_freq[i]

                list_diff.append(diff)
                list_diffN.append(diff_N)
                list_diffS.append(diff_S)


            # minimum differences to founder strain per site
            m_diff = min(list_diff)
            m_diffN = min(list_diffN)
            m_diffS = min(list_diffS)


            print(each.name, m_diff, m_diffN, m_diffS, freq, m_diff/seq_length)

            divergence.extend([float(m_diff)/float(seq_length)])
            divergence_N.extend([float(m_diffN)/float(seq_length)])
            divergence_S.extend([float(m_diffS)/float(seq_length)])


        #print(t, p)
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



        # if ("gag" in fastain):
        #     csvfile_gag_b.write(patient_id + "," + str(sorted_timepoints[t]) + "," +
        #                             str(np.mean(divergence)) + "," + str(np.percentile(divergence, 50)) + "," +
        #                             str(np.percentile(divergence, 5)) + "," + str(np.percentile(divergence, 95)) + "," +
        #                             str(np.mean(divergence_N)) + "," + str(np.percentile(divergence_N, 50)) + "," +
        #                             str(np.percentile(divergence_N, 5)) + "," + str(np.percentile(divergence_N, 95)) + "," +
        #                             str(np.mean(divergence_S)) + "," + str(np.percentile(divergence_S, 50)) + "," +
        #                             str(np.percentile(divergence_S, 5)) + "," + str(np.percentile(divergence_S, 95)) + "\n")
        #
        #     csvfile_gag_b.flush()
        #
        # elif ("gp41" in fastain):
        #     csvfile_gp41_b.write(patient_id + "," + str(sorted_timepoints[t]) + "," +
        #                             str(np.mean(divergence)) + "," + str(np.percentile(divergence, 50)) + "," +
        #                             str(np.percentile(divergence, 5)) + "," + str(np.percentile(divergence, 95)) + "," +
        #                             str(np.mean(divergence_N)) + "," + str(np.percentile(divergence_N, 50)) + "," +
        #                             str(np.percentile(divergence_N, 5)) + "," + str(np.percentile(divergence_N, 95)) + "," +
        #                             str(np.mean(divergence_S)) + "," + str(np.percentile(divergence_S, 50)) + "," +
        #                             str(np.percentile(divergence_S, 5)) + "," + str(np.percentile(divergence_S, 95)) + "\n")

        # else:
        #     print "xxx", patient_id, sorted_timepoints[t]
        #
        # print patient_id, sorted_timepoints[t], len(divergence)


def pwd_per_timepoint(fastain):
    # split fasta by time
    splitFasta = split(fastain, 1)

    seqs_by_timepoint = splitFasta[0]
    total_seq = splitFasta[1]

    mean_diversity = []

    patient = []

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

    for t in range(0, len(sorted_timepoints)):

        diversity = []

        seqs_at_t = seqs_by_timepoint[sorted_timepoints[t]]

        for i in range(0, len(seqs_at_t) - 1):

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
                                # print str(seq1.seq[a]).upper(), str(seq2.seq[a]).upper()
                                diff += 1

                pwd_i = (float(diff) / len(seq1.seq))

                diversity.extend([pwd_i] * freq1 * freq2)
                # print t, i, j, pwd_i, freq1, freq2, pwd_i*freq1*freq2/(freq1+freq2), len(seq1.seq)
                # pwd.append(pwd_i)
                # raw_diff.append(float(diff) / float(len(sub_aln)))

        print sorted_timepoints[t], np.mean(diversity)
        mean_diversity.append(np.mean(diversity))

    df = pd.DataFrame.from_items(
        [('Times', sampleTimes), ('Diversity_mean', mean_diversity), ('Patients', patient * len(sampleTimes))])

    return df


def get_seq_freq(seqs_at_t):
    seq_freq = []

    for x in range(0, len(seqs_at_t)):
        seq_i = seqs_at_t[x]
        parts_i = str.split(seq_i.name, "_")
        freq_i = int(parts_i[2].strip())
        seq_freq.append(freq_i)

    return seq_freq


def get_minor_freq(seqs_at_t, site, seq_freq):
    seqs_at_t_array = np.asarray(seqs_at_t)

    site_a = seqs_at_t_array[:, site]

    full_site_a = []

    for i in range(0, len(site_a)):
        full_site_a.extend([site_a[i]] * seq_freq[i])

    site_a_str = ''.join(full_site_a)

    A = site_a_str.count('a')
    C = site_a_str.count('c')
    T = site_a_str.count('t')
    G = site_a_str.count('g')
    N = site_a_str.count('-')

    site_freqs = [A, C, T, G]

    major_freq = max(site_freqs)

    total_seq = sum(site_freqs)

    minor_freq = total_seq - major_freq

    # print total_seq, minor_freq, site, major_freq+minor_freq, len(site_a_str)

    return minor_freq


# diversity over time by each codon position
def pwd_per_timepoint_cp(fastain, cutoff):
    split_fasta = split(fastain, 1)
    seqs_by_timepoint = split_fasta[0]
    total_seq = split_fasta[1]

    mean_diversity_1cp = []
    mean_diversity_2cp = []
    mean_diversity_3cp = []
    mean_diversity = []

    patient = []

    parts = str.split(fastain, "/")
    parts2 = str.split(parts[len(parts) - 1], "_")

    patient.append(parts2[0])

    # nonsyn_sites, syn_sites = number_of_N_and_S_sites(fastain, None)
    # print nonsyn_sites, syn_sites

    sorted_timepoints = seqs_by_timepoint.keys()
    sorted_timepoints.sort(key=natural_keys)

    sampleTimes = []
    for t in sorted_timepoints:
        sampleTimes.append(int(t))

    for t in range(0, len(sorted_timepoints)):

        diversity_1cp = []
        diversity_2cp = []
        diversity_3cp = []
        diversity = []

        seqs_at_t = seqs_by_timepoint[sorted_timepoints[t]]

        seq_length = len(seqs_at_t[0].seq)

        seq_freq = get_seq_freq(seqs_at_t)

        full_minor_freq = []

        full_site_freq = []

        for i in range(seq_length):

            seqs_at_t_array = np.asarray(seqs_at_t)

            site_a = seqs_at_t_array[:, i]

            A = 0;
            C = 0;
            T = 0;
            G = 0;

            for j in range(0, len(seq_freq)):

                if (site_a[j] == 'a'):
                    A += seq_freq[j]
                elif (site_a[j] == 'c'):
                    C += seq_freq[j]
                elif (site_a[j] == 't'):
                    T += seq_freq[j]
                elif (site_a[j] == 'g'):
                    G += seq_freq[j]

            total_seq = sum([A, C, T, G])
            minor_freq = total_seq - max([A, C, T, G])

            full_minor_freq.append(minor_freq)

            full_site_freq.append([A, C, T, G])

        # diversity at 1cp
        for a in range(0, seq_length, 3):

            minor_freq = full_minor_freq[a]

            if minor_freq > total_seq * cutoff:

                site_a = full_site_freq[a]
                total_seq = sum(site_a)

                diff_1cp = 1

                for i in site_a:

                    if i > 0:
                        diff_1cp *= float(i) / float(total_seq)

                diversity_1cp.append(diff_1cp)

            else:
                diversity_1cp.append(0)

            mean_diversity_1cp.append(np.mean(diversity_1cp))

        # diversity at 2cp
        for b in range(1, seq_length, 3):

            minor_freq = full_minor_freq[b]

            if minor_freq > total_seq * cutoff:

                site_a = full_site_freq[b]
                total_seq = sum(site_a)

                diff_2cp = 1

                for i in site_a:

                    if i > 0:
                        diff_2cp *= float(i) / float(total_seq)

                diversity_2cp.append(diff_2cp)

            else:
                diversity_2cp.append(0)

        mean_diversity_2cp.append(np.mean(diversity_2cp))

        # diversity at 3cp
        for c in range(2, seq_length, 3):

            minor_freq = full_minor_freq[c]

            if minor_freq > total_seq * cutoff:

                site_a = full_site_freq[c]
                total_seq = sum(site_a)

                diff_3cp = 1

                for i in site_a:

                    if i > 0:
                        diff_3cp *= float(i) / float(total_seq)

                diversity_3cp.append(diff_3cp)

            else:

                diversity_3cp.append(0)

        mean_diversity_3cp.append(np.mean(diversity_3cp))

        # diversity at all sites
        for d in range(0, seq_length):

            minor_freq = full_minor_freq[d]

            if minor_freq > total_seq * cutoff:

                site_a = full_site_freq[d]
                total_seq = sum(site_a)

                diff = 1

                for i in site_a:

                    if i > 0:
                        diff *= float(i) / float(total_seq)

                diversity.append(diff)

            else:
                diversity.append(0)

        mean_diversity.append(np.mean(diversity))

        print parts2[0], str(sorted_timepoints[t])

        # if("gag" in fastain):
        #     csvfile_gag.write(parts2[0] +","+ str(sorted_timepoints[t])+ "," + str(np.mean(diversity_1cp)) +","+
        #                  str(np.mean(diversity_2cp))+","+str(np.mean(diversity_3cp))+","+str(np.mean(diversity))+"\n")
        #
        # elif("gp41" in fastain):
        #     csvfile_gp41.write(parts2[0] +","+ str(sorted_timepoints[t])+ "," + str(np.mean(diversity_1cp)) +","+
        #                  str(np.mean(diversity_2cp))+","+str(np.mean(diversity_3cp))+","+str(np.mean(diversity))+"\n")



        # print sorted_timepoints[t], parts2[0], np.mean(diversity_1cp), np.mean(diversity_2cp), np.mean(diversity_3cp)
        # df = pd.DataFrame.from_items(
        #     [('Times', sampleTimes), ('Diversity_mean_1cp', mean_diversity_1cp),('Diversity_mean_2cp', mean_diversity_2cp), ('Diversity_mean_1cp', mean_diversity_3cp), ('Patients', patient*len(sampleTimes))])
        # return df



patients = ["p1","p2","p3","p4","p5","p6","p7","p8","p9","p10",
            "p11","p12","p13","p14","p15","p16","p17","p18","p19",
            "p20","p21","p22","p23","p24","p25","p26","p27","p28",
            "p29","p30","p31","p32","p33","p34"]



#New diversity analysis for revisions

# fileout_gag = "/Users/jayna/Dropbox/Latency_in_Rakai_ms/Revisions/diversity_at_cp_p24_through_time.csv"
# csvfile_gag = open('%s' % fileout_gag, 'w')
# csvfile_gag.write("Patient,Timepoint(days),mean_diversity_1cp,mean_diversity_2cp,mean_diversity_3cp,mean_diversity\n")
#
# fileout_gp41 = "/Users/jayna/Dropbox/Latency_in_Rakai_ms/Revisions/diversity_at_cp_gp41_through_time.csv"
# csvfile_gp41 = open('%s' % fileout_gp41, 'w')
# csvfile_gp41.write("Patient,Timepoint(days),mean_diversity_1cp,mean_diversity_2cp,mean_diversity_3cp,mean_diversity\n")

#New divergence analysis for revisions
fileout_gag_b = "~/results/divergence_at_cp_p24_through_time_noBAPS.csv"
csvfile_gag_b = open('%s' % fileout_gag_b, 'w')
csvfile_gag_b.write("Patient,Timepoint(days),mean_divergence,median_divergence,LQ_divergence,UQ_divergence,mean_nonsyn_divergence,"+
                    "median_nonsyn_divergence,LQ_nonsyn_divergence,UQ_nonsyn_divergence,mean_syn_divergence,median_syn_divergence,"+
                     "LQ_syn_divergence,UQ_syn_divergence\n")

fileout_gp41_b = "~/results/divergence_at_cp_gp41_through_time_noBAPS.csv"
csvfile_gp41_b = open('%s' % fileout_gp41_b, 'w')
csvfile_gp41_b.write("Patient,Timepoint(days),mean_divergence,median_divergence,LQ_divergence,UQ_divergence,mean_nonsyn_divergence,"+
                    "median_nonsyn_divergence,LQ_nonsyn_divergence,UQ_nonsyn_divergence,mean_syn_divergence,median_syn_divergence,"+
                     "LQ_syn_divergence,UQ_syn_divergence\n")



for p in patients:

    filein = "~/results/"+p+"_gag.fasta"

    xx = divergence(filein, p, 0.56/100)

csvfile_gag_b.close()

for p in patients:

    filein = "~/results/"+p+"_gp41.fasta"


    xx = divergence(filein, p, 0.56/100)


csvfile_gp41_b.close()
