import glob
import os

import arrowhead.core as awhd
import topdom.core as tpd
from assess.evaluator import evaluate_results


def final_algorithm_with_evaluation(matrix_filepath, expected_results, results_filepath):
    # topdom_images_filepath = f"{results_filepath}/topdom/{chromosome}.results.images.png"
    topdom_images_filepath = None

    topdom_metrics_filepath = f"{results_filepath}/topdom/{chromosome}.results.metrics.txt"
    # topdom_metrics_filepath = None

    topdom_found_filepath = f"{results_filepath}/topdom/{chromosome}.results.found.txt"
    # topdom_found_filepath = None

    topdom_expected_filepath = f"{results_filepath}/topdom/{chromosome}.results.expected.txt"
    # topdom_expected_filepath = None

    imported_and_adjusted_matrix_topdom, results_topdom = tpd.run(matrix_filepath, R, 0.227)
    evaluate_results(algorithm_results=results_topdom,
                     show=False,
                     mtx=imported_and_adjusted_matrix_topdom,
                     expected_results=expected_results,
                     metrics_filepath=topdom_metrics_filepath,
                     images_filepath=topdom_images_filepath,
                     found_filepath=topdom_found_filepath,
                     expected_filepath=topdom_expected_filepath)

    # imported_and_adjusted_matrix_arrowhead, results_arrowhead = awhd.run(matrix_filepath)
    # evaluate_results(algorithm_results=results_arrowhead,
    #                  show=False,
    #                  mtx=imported_and_adjusted_matrix_arrowhead,
    #                  expected_results=expected_results,
    #                  metrics_filepath=f"{results_filepath}/arrowhead/{chromosome}.results.txt",
    #                  images_filepath=f"{results_filepath}/arrowhead/{chromosome}.results.png")


if __name__ == "__main__":
    R = 100000
    data_path = '../data/www.lcqb.upmc.fr/meetu/dataforstudent/HiC/GM12878/25kb_resolution_intrachromosomal'
    chromosomes = {}
    with open(
            '../data/www.lcqb.upmc.fr/meetu/dataforstudent/TAD/GSE63525_GM12878_primary'
            '+replicate_Arrowhead_domainlist.txt',
            'r') as results_file:
        next(results_file)
        for line in results_file:
            str = line.split('\t', 3)
            if not str[0] in chromosomes:
                chromosomes[str[0]] = []
            chromosomes[str[0]].append((int(str[1]) / R, int(str[2]) / R))
    for chromosome in chromosomes:
        for filename in glob.glob(os.path.join(data_path, 'chr' + chromosome + '_*.RAWobserved')):
            final_algorithm_with_evaluation(filename, chromosomes[chromosome], "../results")
