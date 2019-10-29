from scipy.stats import binom
from random import choice
from random import choices
from random import random
from Bio.Seq import reverse_complement
from Bio.Seq import complement

n = 2
print("Введите вероятность для биномиального распределения остановки Артемиды:")
p = float(input())

DNA_pk = binom.rvs(n, p, 1)

S_1 = []
S_2 = []
for i in range(DNA_pk+1):
    S_1.append(choice(['A', 'C', 'T', 'G']))
    S_2.append(choice(['A', 'C', 'T', 'G']))
S_1 = ''.join(S_1)
S_2 = ''.join(S_2)

S_1 = S_1 + reverse_complement(S_1[1:(DNA_pk+1)])
S_2 = S_2 + reverse_complement(S_2[1:(DNA_pk+1)])


print("Введите вероятность того, что хватит двух комплементарных нуклеотидов подряд, чтобы соединить цепочки:")
p_compl = float(input())
print("Введите вероятность того, что нуклеотиды пристраиваются к первой цепочке:")
p_number_seq = float(input())

if random() < p_compl:
    count_nukl = 2
else:
    count_nukl = 3

version_for_1 = {}
version_for_2 = {}
if count_nukl == 2:
    version_for_1[S_1[-2:]] = len(S_1)-2
    version_for_2[S_2[:2]] = len(S_2)-2
else:
    version_for_1[S_1[-3:]] = len(S_1)-3
    version_for_2[S_2[:3]] = len(S_2)-3


if count_nukl == 2:
    while (version_for_2.get(complement(S_1[-2:])) == None) and (version_for_1.get(complement(S_2[:2])) == None):
        if random() < p_number_seq:
            S_1 = S_1 + choice(['A', 'C', 'T', 'G'])
            version_for_1[S_1[-2:]] = len(S_1)-2
        else:
            S_2 = choice(['A', 'C', 'T', 'G']) + S_2
            version_for_2[S_2[:2]] = len(S_2)-2
        
    if version_for_2.get(complement(S_1[-2:])) != None:
        index = version_for_2.get(complement(S_1[-2:]))
        S = S_1 + complement(S_2[(-index-1):])
        
    else:
        index = version_for_1.get(complement(S_2[:2]))
        S = S_1[:index] + complement(S_2)
        
else:
    while (version_for_2.get(complement(S_1[-3:])) == None) and (version_for_1.get(complement(S_2[:3])) == None):
        if random() < p_number_seq:
            S_1 = S_1 + choice(['A', 'C', 'T', 'G'])
            version_for_1[S_1[-3:]] = len(S_1)-3
        else:
            S_2 = choice(['A', 'C', 'T', 'G']) + S_2
            version_for_2[S_2[:3]] = len(S_2)-3
        
    if version_for_2.get(complement(S_1[-3:])) != None:
        index = version_for_2.get(complement(S_1[-3:]))
        S = S_1 + complement(S_2[(-index-1):])
        
    else:
        index = version_for_1.get(complement(S_2[:3]))
        S = S_1[:index] + complement(S_2)
        
print(S)



