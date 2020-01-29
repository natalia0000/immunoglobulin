from scipy.stats import binom
from random import choice
from random import choices
from random import random
from Bio.Seq import reverse_complement
from Bio.Seq import complement


p_number_seq = 0.5  # the probability that nucleotides attach to the first strand
A_weig = 15  # nucleotide A weight
C_weig = 35  # nucleotide C weight
T_weig = 15  # nucleotide T weight
G_weig = 35  # nucleotide G weight


def add_nucl(S_1, S_2, p_number_seq, count_nucl):
    """
    Adds nucleotides until the linking happens

    :param S_1: the first sequence
    :param S_2: the second sequence
    :param p_number_seq: the probability that nucleotides attach to the first strand
    :param count_nucl: the number of nucleotides that enough for linkage
    :return: linked sequence
    """
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
    """
    Produces sequence as it was reduced by exonuclease.
    Function bases on markov chain with two states: delete nucleotides and not delete.

    :param S: the sequence to be reduced
    :param exo_start_delete: the probability that exonuclease will begin to remove nucleotides
    :param exo_stop_delete: the probability that exonuclease stops removing nucleotides
    :return: sequence as it was reduced by exonuclease
    """
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
    """
    Gives sequence of nucleotides that can be possible junction between V and D regions

    :param p: the probability for the binomial distribution of Artemis stop
    :param p_compl: the probability that two complementary nucleotides in a row is enough to connect the chains
    :param exo_start_delete: the probability that exonuclease will begin to remove nucleotides
    :param exo_stop_delete: the probability that exonuclease stops removing nucleotides
    :return: sequence of nucleotides that can be possible junction between V and D regions
    """
    n = 3
    DNA_pk = binom.rvs(n, p, 1)

    S_1 = ''.join(choice(['A', 'C', 'T', 'G']) for _ in range(DNA_pk))
    S_2 = ''.join(choice(['A', 'C', 'T', 'G']) for _ in range(DNA_pk))

    S_1 = S_1 + reverse_complement(S_1)
    S_2 = S_2 + reverse_complement(S_2)

    if random() < p_compl:
        count_nucl = 2
    else:
        count_nucl = 3

    S = add_nucl(S_1, S_2, p_number_seq, count_nucl)
    S = exo_work(S, exo_start_delete, exo_stop_delete)
    return S
