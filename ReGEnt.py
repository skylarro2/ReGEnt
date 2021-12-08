#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct  9 15:44:24 2018
@author: skylar
"""

import argparse
import math
import time
import sys
import numpy as np
from Bio import AlignIO               # Read fasta files
from sklearn.cluster import SpectralClustering

def ewhole(grp, r_prob, aln):
    """Calculates family entropy."""
    col = len(aln[0].seq)
    row = [len(aln)]*col
    for i in range(col):
        row[i] = row[i] * (1-grp["gap"][i])
    aa_p = [[0] * 20 for i in range(col)]
    aa_ajus = [[0] * 20 for i in range(col)]
    fam_ent = [0]*col

    for i in range(col):
        for j in range(20):
            if grp["gap"][i] > GAPS:
                break
            aa_p[i][j] = grp["ps"][i]*float(r_prob[j])
            aa_ajus[i][j] = (row[i] / (row[i] + grp["ps"][i]) * grp["aa_ct"][i][j] / row[i] +
                             grp["ps"][i] / (row[i] + grp["ps"][i]) * aa_p[i][j] / grp["ps"][i])
    for i in range(col):
        for j in range(20):
            if grp["gap"][i] > GAPS:
                break
            fam_ent[i] += aa_ajus[i][j]*math.log(aa_ajus[i][j]/float(r_prob[j]), 2)
    grp["aa_p_fam"] = fam_ent


def grpent(in_grp, out_grp):
    """Calculates Group Entropy"""
    global REFRENCE
    col = in_grp["col"]
    in_ent = [0]*col
    out_ent = [0]*col
    tot_ent = [0]*col

    for i in range(col):               # Finds Family Entropy
        if in_grp["gap"][i] < GAPS and out_grp["gap"][i] < GAPS:
            tot_ent[i] = np.corrcoef(in_grp["aa_ajus"][i],out_grp["aa_ajus"][i])[0,1]

    in_grp["in_ent"] = in_ent
    in_grp["out_ent"] = out_ent
    in_grp["tot_ent"] = tot_ent

def cert(in_grp, out_grp, all_grouped, certn, num_in):
    col = in_grp["col"]
    grps = all_grouped
    imp_pos = [None]*len(grps)
    for j in range(len(grps)):
        imp_pos[j] = [0]*5
    count = 0
    cert = [0]*col
    for i in range(col):               # Finds Family Entropy
        if in_grp["gap"][i] < GAPS and out_grp["gap"][i] < GAPS and in_grp["tot_ent"][i] < certn:
            for j in range(len(grps)):
                imp_pos[j][count] = ret_num(grps[j].seq[i])
            count+=1

    clstr = SpectralClustering(n_clusters=2,
                               assign_labels='discretize',
                               random_state=0).fit(imp_pos).labels_
    sure = 0
    for i in range(len(grps)):
        if i < num_in and clstr[i] == 1 or i>= num_in and clstr[i] == 0:
            sure +=1
    
    sure = sure / len(grps)
    for i in range(col):               # Finds Family Entropy
        if in_grp["gap"][i] < GAPS and out_grp["gap"][i] < GAPS and in_grp["tot_ent"][i] < certn:
            cert[i] = sure
    in_grp["cert"] = cert
def ret_num(char_aa):
    """Converts char to a numerical representation of the Amino acid"""
    if char_aa == 'X':
        char_aa = 'A'
    temp = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L',
            'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V', '-', '.']
    str(char_aa).capitalize()
    for i in range(22):
        if char_aa == temp[i]:
            aa_val = i
    return aa_val


def ret_aa(int_ascii):
    """Converts numerical representation back to a char."""
    temp = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L',
            'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V', '-', '.']
    return temp[int_ascii]


def get_raw(aln, ret):
    """Passes the gap percentages and Amino Acid counts back by reference."""
    global GAPS
    rows = len(aln)
    col = ret["col"]
    gap = [0]*col
    temp = ""
    num_gaps = 0
    aa_ct = [[0]*22 for i in range(col)]
    for i in range(col):
        for j in range(rows):
            char = (aln[j].seq)[i]
            if char in ('.', '-'):
                num_gaps += 1
        gap[i] = num_gaps/rows
        num_gaps = 0

    for i in range(col):
        temp = ""
        for j in range(rows):
            char = (aln[j].seq)[i]

            if gap[i] <= GAPS:
                temp += char
                aa_ct[i][ret_num(char)] += 1
    ret["aa_ct"] = aa_ct
    ret["gap"] = gap


def get_consensus(ret):
    """Passes back consensus Amino Acid by reference."""
    global GAPS
    col = ret["col"]
    aa_ct = ret["aa_ct"]
    gap = ret["gap"]
    con = ['']*col
    for i in range(col):
        max_v = 0
        temp = ''
        for j in range(20):       # Finds Consensus AA
            if gap[i] > GAPS:
                con[i] = '.'
                break
            elif max_v <= aa_ct[i][j]:
                temp = ret_aa(j)
                max_v = aa_ct[i][j]
                con[i] = temp
    ret["con"] = con


def get_counts(ret):
    """Passes Pseudocount and amino acid totals by reference."""
    global PC_MUL
    global GAPS
    col = ret["col"]
    aa_ct = ret["aa_ct"]
    gap = ret["gap"]
    pseudocount = [0]*col
    aa_tot = [0]*col

    for i in range(col):
        acids = 0
        for j in range(20):
            if gap[i] > GAPS:
                break
            if aa_ct[i][j] > 0:
                acids += 1
        aa_tot[i] = acids
        pseudocount[i] = acids * PC_MUL
    ret["aa_tot"] = aa_tot
    ret["ps"] = pseudocount


def ajust_counts(sequences, odds, ret):
    """Passes back both Pseudocount totals and the Ajusted total"""
    global GAPS
    col = ret["col"]
    aa_pc = ret["ps"]
    aa_ct = ret["aa_ct"]
    gap = ret["gap"]
    row = [len(sequences)]*col
    for i in range(col):
        row[i] = row[i] * (1-gap[i])
    aa_ps_tot = [[0] * 20 for i in range(col)]
    aa_ajus = [[0] * 20 for i in range(col)]

    for i in range(col):
        for j in range(20):  # j is the positon
            if gap[i] > GAPS:
                break
            aa_ps_tot[i][j] = aa_pc[i] * float(odds[j])  # there is a much better forumla one might be able to use here.
            aa_ajus[i][j] = (row[i] / (row[i] + aa_pc[i]) * aa_ct[i][j] / row[i] +
                             aa_pc[i] / (row[i] + aa_pc[i]) * aa_ps_tot[i][j] / aa_pc[i])
    ret["aa_ps_tot"] = aa_ps_tot
    ret["aa_ajus"] = aa_ajus


def ref_position():
    """Finds position within reference sequence."""
    global REFRENCE
    col = len(REFRENCE.seq)
    out = [0]*col
    count = 0
    for i in range(col):
        if REFRENCE.seq[i] != '.' and REFRENCE.seq[i] != '-':
            count += 1
            out[i] = count
    return out


class Groups():
    """Handles data storage."""
    grp_dict = dict()
    dat_dict = dict()
    grp_data_list = dict()
    aln = None
    odds = None
    certianty = None

    def __init__(self, alignments, delimiter, groups, odds,cert):
        temp_dict = dict()
        for i in alignments:
            temp_dict[i.name] = i
        temp_int = 0
        global REFRENCE
        for line in groups:
            if line == delimiter + "\n":
                temp = "Group" + str(temp_int)
                temp_int += 1
                self.grp_dict[temp] = list()
            elif line[0] == ">" and line.find(' ') != -1:
                self.grp_dict[temp].append(temp_dict[line[1:line.find(' ')]])
            elif line[0] == ">":
                self.grp_dict[temp].append(temp_dict[line[1:len(line)-1]])
        REFRENCE = temp_dict[REFRENCE]
        self.aln = alignments
        self.odds = odds
        self.certianty = cert

    def proc_data(self):
        """Calls all calculation functions."""
        self.dat_dict["col"] = len(self.aln[0])
        get_raw(self.aln, self.dat_dict)
        get_consensus(self.dat_dict)
        get_counts(self.dat_dict)
        ewhole(self.dat_dict, self.odds, self.aln)

        for grp in list(self.grp_dict.keys()):
            in_grp = self.get_group(grp)[0]
            out_grp = self.get_group(grp)[1]
            self.grp_data_list[grp + "_in"] = dict()
            self.grp_data_list[grp + "_out"] = dict()
            self.grp_data_list[grp + "_in"]["col"] = len(in_grp[0].seq)
            self.grp_data_list[grp + "_out"]["col"] = len(out_grp[0].seq)
            get_raw(in_grp, self.grp_data_list[grp + "_in"])
            get_raw(out_grp, self.grp_data_list[grp + "_out"])
            get_counts(self.grp_data_list[grp + "_in"])
            get_counts(self.grp_data_list[grp + "_out"])
            get_consensus(self.grp_data_list[grp + "_in"])
            ajust_counts(in_grp, self.odds, self.grp_data_list[grp + "_in"])
            ajust_counts(out_grp, self.odds, self.grp_data_list[grp + "_out"])
            grpent(self.grp_data_list[grp + "_in"], 
                   self.grp_data_list[grp + "_out"])
            cert(self.grp_data_list[grp + "_in"], 
                   self.grp_data_list[grp + "_out"],self.get_all_grouped(),self.certianty, len(in_grp))

    def get_group(self, group_name):
        """Returns in_group and out_group."""
        out_group = []
        for k in self.grp_dict:
            if k != group_name:
                out_group += self.grp_dict[k]
        return [self.grp_dict[group_name], out_group]

    def get_all_grouped(self):
        """Returns all group sequences."""
        out_group = []
        for k in list(self.grp_dict.keys()):
            out_group += self.grp_dict[k]
        return out_group


def get_top(grp_ent, ref):
    """Returns top residues."""
    i = sorted(range(len(grp_ent)), key=lambda i: grp_ent[i])[0:]
    ret = []
    i = i[::-1]
    for j in range(10):
        ret.append(ref[i[j]])
    return ret


def main():
    """handles IO and calls necessary functions from the group class"""
    global REFRENCE
    global GAPS
    start_time = time.time()
    parser = argparse.ArgumentParser(
        description='Calculate group entropy')
    parser.add_argument('Alignment', type=str,
                        help='Alignment file in fasta format')
    parser.add_argument('GrDel', type=str,
                        help='What character(s) delimit different groups')
    parser.add_argument('PsMultiplier', type=float,
                        help='Pseudocount Multiplier')
    parser.add_argument('RSequence', type=str,
                        help='Sequence to base positions on')
    parser.add_argument('gapScore', type=float,
                        help='Allowed percentage of gaps')
    parser.add_argument('grpEntCutoff', type=float,
                        help='below this point, a certianty score will be added')
    global PC_MUL
    REFRENCE = parser.parse_args().RSequence
    PC_MUL = float(parser.parse_args().PsMultiplier)
    GAPS = parser.parse_args().gapScore
    delimiter = parser.parse_args().GrDel
    aln = open("aln.temp", "w")
    fasta = open(parser.parse_args().Alignment, 'r')

    # Taken from LG matrix
    matrix = list({0.079066, 0.055941, 0.041977, 0.053052, 0.012937, 0.040767,
                   0.071586, 0.057337, 0.022355, 0.062157, 0.099081, 0.064600,
                   0.022951, 0.042302, 0.044040, 0.061197, 0.053287, 0.012066,
                   0.034155, 0.069147})

    # removes delimiter from file for alignment processing
    file = ""
    for line in fasta:
        if delimiter not in line:
            file += line
    aln.write(file)

    aln = AlignIO.read("aln.temp", "fasta")
    fasta = open(parser.parse_args().Alignment, 'r')

    # Does all the calculations
    grp = Groups(aln, delimiter, fasta, matrix,parser.parse_args().grpEntCutoff)
    grp.proc_data()

    # Prints out all the groups to csv files
    for group in list(grp.grp_dict.keys()):
        out = open(group + ".csv", 'w')
        out.write("Position,Family Entropy,Group Entropy," +
                  "Certianty,Highest Group AA," +
                  "Highest Family AA,Ref. Position\n")
        for i in range(grp.dat_dict["col"]):
            if not grp.dat_dict["gap"][i] > GAPS:
                out.write(str(i + 1) + ',' +
                          str(grp.dat_dict["aa_p_fam"][i]) + ',' +
                          str(grp.grp_data_list[group + "_in"]["tot_ent"][i]) + ',' +
                          str(grp.grp_data_list[group + "_in"]["cert"][i]) + ',' +
                          str(grp.grp_data_list[group + "_in"]["con"][i]) + ',' +
                          str(grp.dat_dict["con"][i])+',' +
                          str(ref_position()[i]) + "\n")

    sys.stdout.write("GEnt Complete in " + str(time.time()-start_time) + " seconds\n")
    out = open("GEnt.out", 'w')
    out.write("GEnt.py\n")
    out.write("Written by Skylar Olson\n")
    out.write("Based on An algorithm for identification and ranking of " +
              "family-specific residues, applied to the ALDH3 family\n")
    out.write("Written by John Hempel, John Perozich, Troy Wymore, and Hugh B. Nicholas\n")


if __name__ == '__main__':
    main()
