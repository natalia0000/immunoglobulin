from scipy.stats import binom
from random import choice
from random import choices
from random import random
from Bio.Seq import reverse_complement
from Bio.Seq import complement


def add_nucl_2(S_1, S_2, p_number_seq):
    version_for_1 = {}
    version_for_2 = {}
    version_for_1[S_1[-2:]] = len(S_1) - 2
    version_for_2[S_2[:2]] = len(S_2) - 2
    while (version_for_2.get(complement(S_1[-2:])) is None) and (version_for_1.get(complement(S_2[:2])) is None):
        if random() < p_number_seq:
            S_1 = S_1 + choices(['A', 'C', 'T', 'G'], weights=[15, 35, 15, 35])[0]
            version_for_1[S_1[-2:]] = len(S_1) - 2
        else:
            S_2 = choices(['A', 'C', 'T', 'G'], weights=[15, 35, 15, 35])[0] + S_2
            version_for_2[S_2[:2]] = len(S_2) - 2

    if version_for_2.get(complement(S_1[-2:])) is not None:
        index = version_for_2.get(complement(S_1[-2:]))
        S = S_1 + complement(S_2[(-index - 1):])

    else:
        index = version_for_1.get(complement(S_2[:2]))
        S = S_1[:index] + complement(S_2)

    return S


def add_nucl_3(S_1, S_2, p_number_seq):
    version_for_1 = {}
    version_for_2 = {}
    version_for_1[S_1[-3:]] = len(S_1) - 3
    version_for_2[S_2[:3]] = len(S_2) - 3
    while (version_for_2.get(complement(S_1[-3:])) is None) and (version_for_1.get(complement(S_2[:3])) is None):
        if random() < p_number_seq:
            S_1 = S_1 + choice(['A', 'C', 'T', 'G'])
            version_for_1[S_1[-3:]] = len(S_1) - 3
        else:
            S_2 = choice(['A', 'C', 'T', 'G']) + S_2
            version_for_2[S_2[:3]] = len(S_2) - 3

    if version_for_2.get(complement(S_1[-3:])) is not None:
        index = version_for_2.get(complement(S_1[-3:]))
        S = S_1 + complement(S_2[(-index - 1):])

    else:
        index = version_for_1.get(complement(S_2[:3]))
        S = S_1[:index] + complement(S_2)

    return S


def exo_work(S, exo_start_delete, exo_stop_delete):
    status = 1
    if random() < exo_start_delete:
        status = 2
    for i in range(len(S)):
        if status == 1:
            if random() < exo_start_delete:
                status = 2
        else:
            S = S[:i] + S[(i + 1):]
            i -= 1
            if random() < exo_stop_delete:
                status = 1

    return S


n = 3
print("Введите вероятность для биномиального распределения остановки Артемиды:")
p = float(input())
# p = 0.5

DNA_pk = binom.rvs(n, p, 1)

S_1 = []
S_2 = []
for i in range(DNA_pk + 1):
    S_1.append(choice(['A', 'C', 'T', 'G']))
    S_2.append(choice(['A', 'C', 'T', 'G']))
S_1 = ''.join(S_1)
S_2 = ''.join(S_2)

S_1 = S_1 + reverse_complement(S_1[1:(DNA_pk + 1)])
S_2 = S_2 + reverse_complement(S_2[1:(DNA_pk + 1)])

print("Введите вероятность того, что хватит двух комплементарных нуклеотидов подряд, чтобы соединить цепочки:")
p_compl = float(input())
# p_compl = 0.8
print("Введите вероятность того, что нуклеотиды пристраиваются к первой цепочке:")
p_number_seq = float(input())
# p_number_seq = 0.5

if random() < p_compl:
    count_nukl = 2
else:
    count_nukl = 3

if count_nukl == 2:
    S = add_nucl_2(S_1, S_2, p_number_seq)

else:
    S = add_nucl_3(S_1, S_2, p_number_seq)
print(S)

print("Введите вероятность того, что экзонуклеаза станет удалять нуклеотиды:")
exo_start_delete = float(input())
print("Введите вероятность того, что экзонуклеаза престанет удалять нуклеотиды:")
exo_stop_delete = float(input())
# exo_start_delete = 0.4
# exo_stop_delete = 0.55

S = exo_work(S, exo_start_delete, exo_stop_delete)
print(S)
