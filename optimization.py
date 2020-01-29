"""
Optimization script can be used for finding optimal VD modal parameters
for given distribution of junction lengths
"""
from VD import main_fun
from scipy.optimize import minimize
from scipy.optimize import basinhopping
import numpy as np
from numpy import genfromtxt

# p = 0.5  # the probability for the binomial distribution of Artemis stop
# p_compl = 0.8  # the probability that two complementary nucleotides in a row is enough to connect the chains
# exo_start_delete = 0.4  # the probability that exonuclease will begin to remove nucleotides
# exo_stop_delete = 0.7  # the probability that exonuclease stops removing nucleotides
index = 4  # data set number from 0 to 12(model_data.csv)


def mod_for_opt(p_arr):
    global freq_norm
    p = p_arr[0] % 1
    p_compl = p_arr[1] % 1
    exo_start_delete = p_arr[2] % 1
    exo_stop_delete = p_arr[3] % 1
    number_of_loop = 1000
    number_of_more_1 = 0
    dict_with_len = {}
    Summ = 0
    for i in range(number_of_loop):
        n = len(main_fun(p, p_compl, exo_start_delete, exo_stop_delete))
        if n > 0:
            number_of_more_1 += 1
            dict_with_len[n - 1] = dict_with_len.get(n - 1, 0) + 1
    length = max(dict_with_len.keys())
    for i in range(length):
        Summ = Summ + (dict_with_len.get(i, 0) / number_of_more_1 - freq_norm[i]) ** 2
    return Summ ** 0.5


def make_points(step):
    matrix = np.meshgrid(np.arange(0.01, 1, step), np.arange(0.01, 1, step), np.arange(0.01, 1, step),
                         np.arange(0.01, 1, step))
    coord = [each_one.ravel() for each_one in matrix]
    points = np.vstack(coord).T
    return points


points = make_points(0.1)
smaller = 0.03
x_arr = []
fun_arr = []

my_data = genfromtxt('model_data.csv', delimiter=',')
my_data = np.vstack(my_data).T
freq_norm = my_data[index]
freq_norm = np.concatenate((freq_norm, [0]*100), axis=None)

for point in points:
    result = minimize(mod_for_opt, point, method='SLSQP',
                   bounds=((0.01, 0.999), (0.01, 0.999), (0.01, 0.999), (0.01, 0.999)))
    print(result.fun, result.x)
    if result.fun < smaller:
        x_arr.append(result.x)
        fun_arr.append(result.fun)
print(x_arr)
print(fun_arr)


# uncomment following code to run additional optimizer for refinement
# result2 = basinhopping(mod_for_opt, [0.21, 0.11, 0.21, 0.21], stepsize=0.01)
# print(result2.fun, result2.x)
