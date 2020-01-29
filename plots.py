from math import log
from VD import main_fun
import matplotlib.pyplot as plt
import numpy as np
from numpy import genfromtxt


my_data = genfromtxt('model_data.csv', delimiter=',')
my_data = np.vstack(my_data).T
index = 4

p = 0.21477354
p_compl = 0.10367869
exo_start_delete = 0.21482267
exo_stop_delete = 0.20595069
freq_norm = np.concatenate((my_data[index], [0]*100), axis=None)
number_of_loop = 100000
arr_with_len = [0]*100
number_of_more_1 = 0
max_len = 0
for i in range(number_of_loop):
    n = len(main_fun(p, p_compl, exo_start_delete, exo_stop_delete))
    if n > 0:
        number_of_more_1 += 1
        arr_with_len[n-1] += 1
    if n > max_len:
        max_len = n
for i in range(max_len+1):
    arr_with_len[i-1] = arr_with_len[i-1] / number_of_more_1
    # if arr_with_len[i] != 0:
    #     arr_with_len[i] = log(arr_with_len[i])
    # else:
    #     arr_with_len[i] = min(arr_with_len)
    # if freq_norm[i] != 0:
    #     freq_norm[i] = log(freq_norm[i])
    # else:
    #     freq_norm[i] = min(freq_norm) #for log plot


Summ = 0
for i in range(max_len):
    Summ = Summ + (arr_with_len[i] - freq_norm[i]) ** 2
print(Summ ** 0.5)

x = list(np.arange(1, max_len+1))

fig = plt.subplots()
plt.plot(x, arr_with_len[:max_len], color="red")
plt.plot(x, freq_norm[:max_len], color="blue")
plt.xlabel("Length of N1")
plt.ylabel("Prob")
plt.savefig('test.png')
plt.show()
