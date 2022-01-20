import argparse
import glob
import os

import arrowhead.core as awhd
import topdom.core as tpd
from assess.evaluator import evaluate_results


def final_algorithm_with_evaluation(matrix_filepath, expected_results, results_filepath,
                                    should_run_arrowhead, should_run_topdom,
                                    should_dump_metrics, should_dump_coordinates, should_dump_expected_coordinates,
                                    should_dump_images, topdom_normalization_alpha, topdom_sensitivity,
                                    topdom_window_size, topdom_pval_limit, should_show):
    if should_run_topdom:
        topdom_images_filepath =\
            f"{results_filepath}/topdom/{chromosome}.results.images.png" if should_dump_images else None

        topdom_metrics_filepath =\
            f"{results_filepath}/topdom/{chromosome}.results.metrics.txt" if should_dump_metrics else None

        topdom_found_filepath =\
            f"{results_filepath}/topdom/{chromosome}.results.found.txt" if should_dump_coordinates else None

        topdom_expected_filepath =\
            f"{results_filepath}/topdom/{chromosome}.results.expected.txt" if should_dump_expected_coordinates else None

        imported_and_adjusted_matrix_topdom, results_topdom = tpd.run(
            matrix_filepath=matrix_filepath,
            R=R,
            alpha=topdom_normalization_alpha,
            window_size=topdom_window_size,
            sensitivity=topdom_sensitivity,
            pval_limit=topdom_pval_limit)
        evaluate_results(algorithm_results=results_topdom,
                         show=should_show,
                         mtx=imported_and_adjusted_matrix_topdom,
                         expected_results=expected_results,
                         metrics_filepath=topdom_metrics_filepath,
                         images_filepath=topdom_images_filepath,
                         found_filepath=topdom_found_filepath,
                         expected_filepath=topdom_expected_filepath)
    if should_run_arrowhead:
        arrowhead_images_filepath =\
            f"{results_filepath}/arrowhead/{chromosome}.results.images.png" if should_dump_images else None

        arrowhead_metrics_filepath = \
            f"{results_filepath}/arrowhead/{chromosome}.results.metrics.txt" if should_dump_metrics else None

        arrowhead_found_filepath =\
            f"{results_filepath}/arrowhead/{chromosome}.results.found.txt" if should_dump_coordinates else None

        arrowhead_expected_filepath =\
            f"{results_filepath}/arrowhead/{chromosome}.results.expected.txt" if should_dump_expected_coordinates else None

        imported_and_adjusted_matrix_arrowhead, results_arrowhead = awhd.run(matrix_filepath)
        evaluate_results(algorithm_results=results_arrowhead,
                         show=should_show,
                         mtx=imported_and_adjusted_matrix_arrowhead,
                         expected_results=expected_results,
                         metrics_filepath=arrowhead_metrics_filepath,
                         images_filepath=arrowhead_images_filepath,
                         found_filepath=arrowhead_found_filepath,
                         expected_filepath=arrowhead_expected_filepath)


def init_parser(parser):
    parser.add_argument('--run-arrowhead', type=bool, default=False,
                        help="Indicate whether the Arrowhead should run, omitting means False")

    parser.add_argument('--run-topdom', type=bool, default=False,
                        help="Indicate whether the default algorithm TopDom should run, omitting means False")

    parser.add_argument('--results-path', type=str, default="../results",
                        help="Relative or absolute path to directory where results should be stored")

    parser.add_argument('--chromosomes', type=str, default="22",
                        help="Comma separated list of chromosomes to be processes, eg. 1,22,3,X")

    parser.add_argument('--with-metrics-results', type=bool, default=False,
                        help="Should dump metrics of found TADs, omitting means False")

    parser.add_argument('--with-results-coordinates', type=bool, default=False,
                        help="Should dump found TADs placements, omitting means False")

    parser.add_argument('--with-expected-results-coordinates', type=bool, default=False,
                        help="Should dump expected TADs placements that we compare to, omitting means False")

    parser.add_argument('--with-images-results', type=bool, default=False,
                        help="Should dump images of found TADs, omitting means False")

    parser.add_argument('--topdom-normalization-alpha', type=float, default=None,
                        help="Alpha used for normalization, no normalization if not provided")

    parser.add_argument('--topdom-sensitivity', type=float, default=0.04,
                        help="Sensitivity used in TopDom")

    parser.add_argument('--topdom-window-size', type=int, default=5,
                        help="Window size used in TopDom")

    parser.add_argument('--topdom-pval-limit', type=float, default=0.05,
                        help="Limit for Pval in TopDom")

    parser.add_argument('--should-show', type=bool, default=False,
                        help="Should show each results, omitting means False")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Team WA1 - TAD Sniffer')
    init_parser(parser)
    args = parser.parse_args()
    chosen_chromosomes_to_processing = args.chromosomes.split(sep=",")
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
        if chromosome in chosen_chromosomes_to_processing:
            for filename in glob.glob(os.path.join(data_path, 'chr' + chromosome + '_*.RAWobserved')):
                final_algorithm_with_evaluation(matrix_filepath=filename,
                                                expected_results=chromosomes[chromosome],
                                                results_filepath=args.results_path,
                                                should_run_arrowhead=args.run_arrowhead,
                                                should_run_topdom=args.run_topdom,
                                                should_dump_images=args.with_images_results,
                                                should_dump_metrics=args.with_metrics_results,
                                                should_dump_coordinates=args.with_results_coordinates,
                                                should_dump_expected_coordinates=args.with_expected_results_coordinates,
                                                topdom_normalization_alpha=args.topdom_normalization_alpha,
                                                topdom_sensitivity=args.topdom_sensitivity,
                                                topdom_window_size=args.topdom_window_size,
                                                topdom_pval_limit=args.topdom_pval_limit,
                                                should_show=args.should_show)
