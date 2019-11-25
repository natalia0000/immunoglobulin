from VD import main_fun
from scipy.optimize import minimize
import numpy as np

p = 0.5  # вероятность для биномиального распределения остановки Артемиды
p_compl = 0.8  # вероятность того, что хватит двух комплементарных нуклеотидов подряд, чтобы соединить цепочки
exo_start_delete = 0.4  # вероятность того, что экзонуклеаза станет удалять нуклеотиды
exo_stop_delete = 0.7  # вероятность того, что экзонуклеаза престанет удалять нуклеотиды


def mod_for_opt(p_arr):
    p = p_arr[0]
    p_compl = p_arr[1]
    exo_start_delete = p_arr[2]
    exo_stop_delete = p_arr[3]
    freq_norm = [0.0632506, 0.088070456, 0.094875901, 0.103682946,
                 0.085668535, 0.085668535, 0.076461169, 0.065652522,
                 0.058046437, 0.055244195, 0.040032026, 0.040432346,
                 0.032826261, 0.026821457, 0.017614091, 0.012810248,
                 0.008807046, 0.008406725, 0.008406725, 0.004803843,
                 0.006004804, 0.002401922, 0.002802242, 0.001200961,
                 0.000800641, 0.002802242, 0.001200961, 0.000800641,
                 0.00040032, 0.000800641, 0.001200961, 0.00040032,
                 0.00040032, 0, 0.00040032, 0, 0.00040032, 0.00040032,
                 0, 0, 0, 0, 0, 0] + [0] * 100
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
    matrix = np.meshgrid(np.arange(0.01, 1, step), np.arange(0.01, 1, step), np.arange(0.01, 1, step), np.arange(0.01, 1, step))
    coord = [each_one.ravel() for each_one in matrix]
    points = np.vstack(coord).T
    return points

points = make_points(0.1)
for point in points:
    result = minimize(mod_for_opt, point, method='SLSQP', bounds=((0.01, 0.999), (0.01, 0.999), (0.01, 0.999), (0.01, 0.999)))
    print(result.x, result.fun)
