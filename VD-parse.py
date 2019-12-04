from collections import defaultdict
from scipy.stats import binom
from random import choice
from random import choices
from random import random
from Bio.Seq import reverse_complement
from Bio.Seq import complement
import seaborn
import argparse

#create a parser
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='VD_model')
#add arguments
    parser.add_argument('--p_art', help='Probability of binomial distribution of Artemis stop', default=0.5, type = float)
    #p = 0.5  # вероятность для биномиального распределения остановки Артемиды
    parser.add_argument('--p_compl', help='Probability of two nucleotides joining complementary strands', default=0.8, type = float)
    #p_compl = 0.8  # вероятность того, что хватит двух комплементарных нуклеотидов подряд, чтобы соединить цепочки
    parser.add_argument('--p_number_seq', help='Probability of nucleotides joining to the first strand', default=0.5, type = float)
    #p_number_seq = 0.5  # вероятность того, что нуклеотиды пристраиваются к первой цепочке
    parser.add_argument('--A_weig', help='Weight of adenine', default=15, type = int)
    parser.add_argument('--T_weig', help='Weight of thymine', default=15, type = int)
    parser.add_argument('--G_weig', help='Weight of guanine', default=35, type = int)
    parser.add_argument('--C_weig', help='Weight of cytosine', default=35, type = int)
    #A_weig = 15  # вес нуклеотида А
    #C_weig = 35  # вес нуклеотида C
    #T_weig = 15  # вес нуклеотида T
    #G_weig = 35  # вес нуклеотида G
    parser.add_argument('--exo_start_delete', help='Probability of exonuclease start', default=0.4, type = float)
    #exo_start_delete = 0.4  # вероятность того, что экзонуклеаза станет удалять нуклеотиды
    parser.add_argument('--exo_stop_delete', help='Probability of exonuclease stop', default=0.7, type = float)
    #exo_stop_delete = 0.7  # вероятность того, что экзонуклеаза престанет удалять нуклеотиды
    parser.add_argument('--seq_num', help = 'Desired number of sequences to be displayed', default=100, type = int)
    parser.add_argument('--visual_flag', help='Flag to build the plot of sequences distribution', action='store_true')
    args = parser.parse_args()
    
#connect arguments with variables in functions
    p = args.p_art
    p_compl = args.p_compl
    p_number_seq = args.p_number_seq
    A_weig = args.A_weig
    T_weig = args.T_weig
    G_weig = args.G_weig
    C_weig = args.C_weig
    exo_start_delete = args.exo_start_delete
    exo_stop_delete = args.exo_stop_delete
    seq_num = args.seq_num
    visual_flag = args.visual_flag

#insert the program
    sequences = []  # variable for a list of processed nucleotide sequences. program will fill it with # seq_num of sequences.

    def add_nucl(S_1, S_2, p_number_seq, count_nucl):
        version_for_1 = {}
        version_for_2 = {}
        for i in range(len(S_1) - count_nucl + 1):
            version_for_1.setdefault(S_1[i:(i + count_nucl)], []).append(i)
            version_for_2.setdefault(S_2[i:(i + count_nucl)], []).append(len(S_1) - i - count_nucl)

        while (complement(S_1[-count_nucl:]) not in version_for_2) and (
                complement(S_2[:count_nucl]) not in version_for_1):
            if random() < p_number_seq:
                S_1 = S_1 + choices(['A', 'C', 'T', 'G'], weights=[A_weig, C_weig, T_weig, G_weig])[0]
                version_for_1.setdefault(S_1[-count_nucl:], []).append(len(S_1) - count_nucl)
            else:
                S_2 = choices(['A', 'C', 'T', 'G'], weights=[A_weig, C_weig, T_weig, G_weig])[0] + S_2
                version_for_2.setdefault(S_2[:count_nucl], []).append(len(S_2) - count_nucl)

        if complement(S_1[-count_nucl:]) in version_for_2:
            index = choice(version_for_2.get(complement(S_1[-count_nucl:])))
            S = S_1 + complement(S_2[(-index - 1):])

        else:
            index = choice(version_for_1.get(complement(S_2[:count_nucl])))
            S = S_1[:index] + complement(S_2)
        return S


    def exo_work(S, exo_start_delete, exo_stop_delete):
        n = len(S)
        flag_arr = []
        status = 1
        if (exo_start_delete + exo_stop_delete) == 0:
            print("Impossible probabilities")
            return "Error"
        if random() < (exo_start_delete) / (exo_start_delete + exo_stop_delete):
            status = 2
        for i in range(n):
            if status == 1:
                flag_arr.append(1)
                if random() < exo_start_delete:
                    status = 2
            else:
                flag_arr.append(0)
                if random() < exo_stop_delete:
                    status = 1
        S_new = ""
        for i in range(n):
            if flag_arr[i]:
                S_new = S_new + S[i]
        return S_new


    def main_fun(p, p_compl, exo_start_delete, exo_stop_delete):
        n = 3
        DNA_pk = binom.rvs(n, p, 1)

        S_1 = ''.join(choice(['A', 'C', 'T', 'G']) for _ in range(DNA_pk))
        S_2 = ''.join(choice(['A', 'C', 'T', 'G']) for _ in range(DNA_pk))

        S_1 = S_1 + reverse_complement(S_1)
        S_2 = S_2 + reverse_complement(S_2)

        base = (DNA_pk, S_1, S_2)

        if random() < p_compl:
            count_nucl = 2
        else:
            count_nucl = 3

        S_add = add_nucl(S_1, S_2, p_number_seq, count_nucl)
        S_del = exo_work(S_add, exo_start_delete, exo_stop_delete)
        sequences.append(S_del)

    #quantify results to seq_num with for loop
    for el in range(1,seq_num):
        main_fun(p, p_compl, exo_start_delete, exo_stop_delete)

    # visualize
    def visual(visual_flag):
        seq_len = []
        for el in sequences:
            seq_len.append(len(el))
        if visual_flag:
            print(seaborn.distplot(seq_len, axlabel='Length of sequences', bins=20))
        else:
            seq_len_c = [[el, seq_len.count(el)] for el in set(seq_len)]
            for el in seq_len_c:
                el[1] = float(el[1]) / seq_num
            print(seq_len_c)
