import sys
import os
import timeit
import numpy

sys.path.append(os.curdir)

try:
    from analyzer.analyzer import calculate_selected_reaction_propensities, calculate_sum_of_reaction_propensities
except ImportError:
    print("analyzer not found in current working dir", file=sys.stderr)
    raise


def main():
    length = 1000000
    iterations = 100
    components = numpy.random.random_sample((2, length))
    reaction_k = numpy.array([0.5, 0.8])
    reaction_reactants = numpy.array([[0], [-1]])
    reaction_events = numpy.random.randint(0, 2, size=length)

    def bench(): return calculate_selected_reaction_propensities(
        components, reaction_k, reaction_reactants, reaction_events)

    timing = timeit.timeit(bench, setup=bench, number=iterations)

    print('calculate_selected_reaction_propensities took {} s per iteration'.format(
        timing/iterations))

    def bench_sum(): return calculate_sum_of_reaction_propensities(
        components, reaction_k, reaction_reactants)

    timing = timeit.timeit(bench_sum, setup=bench_sum, number=iterations)

    print('calculate_sum_of_reaction_propensities took {} s per iteration'.format(
        timing/iterations))

    # propensities = numpy.empty(reaction_events.shape)
    # def bench_gpu(): return gpu_selected_reaction_propensities(
    #     components, reaction_k, reaction_reactants, reaction_events, propensities)

    # timing = timeit.timeit(bench_gpu, setup=bench_gpu, number=iterations)

    # print('gpu_selected_reaction_propensities took {} s per iteration'.format(
    #     timing/iterations))


if __name__ == '__main__':
    main()
