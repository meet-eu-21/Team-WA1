import math
import warnings

import matplotlib.pyplot as plt
from scipy.stats import wilcoxon

import topdom.import_matrix as import_matrix

warnings.filterwarnings("ignore")


def l2d_to_1d(list2d):
    return [j for i in list2d for j in i]


def generate_binsignal(matrix, k):
    n = len(matrix)
    binsignal = []
    for i in range(n):
        sig = matrix[i - k:i, i:i + k].mean()
        binsignal.append(0 if math.isinf(sig) else sig)
    return binsignal[k: n - k]


def flatten_binsignal(binsignal, sensitivity, bins):
    ################################ POCZĄTEK WYPŁASZCZANIA
    # diff to różnica, na podstawie której znajdowane będą obszary zawierające ekstrema, których
    # wartości wahają się o małe odchylenia
    # zmniejszając go z np 0.05 do 0.025 wzrasta czułość wykrywania tadów dla danych testowych,
    # polecam sprawdzić dla różnych parametrów i się przekonać, dla 0.025 były najlepsze rezulatay,
    # poniżej zaczynało już tworzyć dziwne arefakty
    min_bin = min(binsignal)
    max_bin = max(binsignal)
    diff = (max_bin - min_bin) * sensitivity
    # wydobycie WSZYSTKICH ekstremów lokalnych z binsignalu
    extremes_positions = []
    extremes_values = []
    for i in range(1, len(binsignal) - 1):
        if binsignal[i - 1] > binsignal[i] < binsignal[i + 1]:
            extremes_positions.append(bins[i])
            extremes_values.append(binsignal[i])
        elif binsignal[i - 1] < binsignal[i] > binsignal[i + 1]:
            extremes_positions.append(bins[i])
            extremes_values.append(binsignal[i])
    # sztuczne dodanie na końcu minimum lokalnego o takiej wartości, że nie zostanie
    # poźniej przybliżone poziomą prostą
    extremes_values[-1] = extremes_values[-2] - 2 * diff
    # wydobycie z ekstremów informacji o tym, które rejony można przybliżyć
    # poziomą prostą
    flat_areas = []
    flat_area_beginning = 0
    for i in range(1, len(extremes_positions)):
        if abs(extremes_values[i] - extremes_values[i - 1]) > diff:
            if extremes_positions[i - 1] - extremes_positions[flat_area_beginning] > 1:
                flat_areas.append((flat_area_beginning, i - 1))
            flat_area_beginning = i
    # stworzenie nowej listy ekstremów, w której obszary płaskie będą reprezentowane jako
    # dwa punktu z początku i końca przedziału, mające taką samą wartość
    # (równą średniej z ekstremów znajdujących się na początku i końcu przedziału)
    new_extremes_positions = []
    new_extremes_values = []
    last = 0
    for i in flat_areas:
        new_extremes_positions += extremes_positions[last + 1:i[0]]
        new_extremes_values += extremes_values[last + 1:i[0]]

        new_extremes_positions.append(extremes_positions[i[0]])
        new_extremes_positions.append(extremes_positions[i[1]])
        new_extremes_values.append((extremes_values[i[0]] + extremes_values[i[1]]) / 2)
        new_extremes_values.append((extremes_values[i[0]] + extremes_values[i[1]]) / 2)
        last = i[1]
    ################################ KONIEC WYPŁASZCZANIA
    return new_extremes_positions, new_extremes_values


def find_minimums(posits, vals):
    min_posits = []
    min_vals = []
    for i in range(1, len(vals) - 1):
        if vals[i - 1] > vals[i] <= vals[i + 1]:
            min_posits.append(posits[i])
            min_vals.append(vals[i])
    return min_posits, min_vals


def statistical_filtering(matrix, min_coords, wsize, msize):
    p_values = []
    for i in min_coords:
        lower = max(1, i - wsize)
        up_mtx = matrix[lower:i, lower:i]
        upper = min(i + wsize, msize)
        down_mtx = matrix[i:upper, i:upper]
        middle_mtx = matrix[i:i + wsize, lower:upper]
        p_value = wilcoxon(l2d_to_1d(middle_mtx), l2d_to_1d(down_mtx) + l2d_to_1d(up_mtx)).pvalue
        p_values.append(p_value)
    return p_values


def filter_coords(coords, p_values, p_limit):
    delete = []
    for i in range(len(coords)):
        if max(p_values[i], p_values[i + 1]) > p_limit:
            delete.append(i)
    delete.reverse()
    for i in delete:
        coords.pop(i)


def vizualize_tads(matrix, min_pos):
    plt.imshow(matrix, cmap='Wistia')
    for i in range(len(min_pos) - 1):
        plt.hlines(min_pos[i], min_pos[i], min_pos[i + 1])
        plt.vlines(min_pos[i], min_pos[i], min_pos[i + 1])
        plt.hlines(min_pos[i + 1], min_pos[i], min_pos[i + 1])
        plt.vlines(min_pos[i + 1], min_pos[i], min_pos[i + 1])
    plt.show()


def vizualize_signal(exr_pos, exr_vals, min_pos, min_vals):
    plt.plot(exr_pos, exr_vals, color='orange')
    plt.grid()
    plt.scatter(min_pos, min_vals, color='red', s=8)
    # plt.savefig('name_here.png', dpi=300)
    plt.show()


def visualize_pvalue(min_pos, pvalues, pvalue_cut):
    if len(min_pos) > 0:
        plt.hlines(pvalue_cut, min_pos[0], min_pos[len(min_pos) - 1], color='red')
        plt.scatter(min_pos, pvalues, color='green', s=8)
        plt.show()


def vizualize_coords(coords, color):
    for coo in coords:
        plt.hlines(coo[0], coo[0], coo[1], colors=color)
        plt.vlines(coo[0], coo[0], coo[1], colors=color)
        plt.hlines(coo[1], coo[0], coo[1], colors=color)
        plt.vlines(coo[1], coo[0], coo[1], colors=color)


def topdom_visual(np_matrix, window_size, sensitivity):
    binsignal = generate_binsignal(np_matrix, window_size)
    print(binsignal)
    plt.plot(binsignal)
    plt.show()
    bins = list(range(window_size, len(np_matrix) - window_size))
    new_extremes_positions, new_extremes_values = flatten_binsignal(binsignal, sensitivity, bins)
    # TODO implement better min search
    minima_positions, minima_values = find_minimums(new_extremes_positions, new_extremes_values)
    min_coords = []
    for i in range(len(minima_positions) - 1):
        min_coords.append((minima_positions[i], minima_positions[i + 1]))
    print("minima: " + str(len(min_coords)))
    p_values = statistical_filtering(np_matrix, minima_positions, window_size, len(np_matrix) - window_size)
    vizualize_signal(new_extremes_positions, new_extremes_values, minima_positions, minima_values)
    visualize_pvalue(minima_positions, p_values, 0.05)
    vizualize_tads(np_matrix, minima_positions)


def topdom(np_matrix, window_size, sensitivity, pval_limit):
    binsignal = generate_binsignal(np_matrix, window_size)
    bins = list(range(window_size, len(np_matrix) - window_size))
    new_extremes_positions, new_extremes_values = flatten_binsignal(binsignal, sensitivity, bins)
    # TODO implement better min search
    minima_positions, minima_values = find_minimums(new_extremes_positions, new_extremes_values)
    min_coords = []
    for i in range(len(minima_positions) - 1):
        min_coords.append((minima_positions[i], minima_positions[i + 1]))
    p_values = statistical_filtering(np_matrix, minima_positions, window_size, len(np_matrix) - window_size)
    filter_coords(min_coords, p_values, pval_limit)
    return min_coords


def run(matrix_filepath, R, alpha):
    mtx = import_matrix.import_matrix(matrix_filepath, R)
    # mtx = import_matrix.import_normalized_matrix(matrix_filepath, R, alpha)
    topdom_coords = topdom(mtx, 5, 0.04, 0.05)
    return mtx, topdom_coords


if __name__ == "__main__":
    a = import_matrix.import_matrix('./data/chr21_25kb.RAWobserved.txt', 25000)
    for win in [5]:
        for sns in [0.04]:
            print("Window size: " + str(win) + ", sensitivity: " + str(sns))
            topdom_visual(a, win, sns)
