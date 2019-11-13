import sys
import os
import timeit
import numpy

sys.path.append(os.curdir)

try:
    from analyzer.analyzer import calculate_selected_reaction_propensities, calculate_sum_of_reaction_propensities, gpu_selected_reaction_propensities, time_average, evaluate_trajectory_at
    from analyzer.stochastic_sim import ReactionNetwork
except ImportError:
    print("analyzer not found in current working dir", file=sys.stderr)
    raise


def main():
    length = 100000
    iterations = 500
    components = numpy.random.random_sample((2, length))

    reactions = ReactionNetwork(2)
    reactions.k = numpy.array([0.5, 0.8])
    reactions.reactants = numpy.array([[0], [-1]], dtype=numpy.int32)
    reaction_events = numpy.random.randint(0, 2, size=length)

    total = 0.0

    def bench(): return calculate_selected_reaction_propensities(
        components, reaction_events, reactions)

    timing = timeit.timeit(bench, setup=bench, number=iterations)
    total += timing

    print('calculate_selected_reaction_propensities took {} s per iteration'.format(
        timing/iterations))

    def bench_sum(): return calculate_sum_of_reaction_propensities(components, reactions)

    timing = timeit.timeit(bench_sum, setup=bench_sum, number=iterations)
    total += timing

    print('calculate_sum_of_reaction_propensities took {} s per iteration'.format(
        timing/iterations))

    signal_timestamps = -numpy.log(numpy.random.random_sample(size=length))
    response_timestamps = -numpy.log(numpy.random.random_sample(size=length))

    output = numpy.zeros(length-1)
    evaluated = numpy.zeros(length-1)
    def bench_average(): return time_average(
        components[0], signal_timestamps, response_timestamps, out=output, evaluated=evaluated)

    timing = timeit.timeit(
        bench_average, setup=bench_average, number=iterations)
    total += timing

    print('time_average took {} s per iteration'.format(
        timing/iterations))

    def bench_evaluate(): return evaluate_trajectory_at(
        components[0], signal_timestamps, response_timestamps[1:], out=output)

    timing = timeit.timeit(
        bench_evaluate, setup=bench_evaluate, number=iterations)
    total += timing

    print('evaluate_trajectory_at took {} s per iteration'.format(
        timing/iterations))

    print('total: {} s per iteration'.format(total/iterations))

    # propensities = numpy.empty(reaction_events.shape)
    # def bench_gpu(): return gpu_selected_reaction_propensities(
    #     components, reaction_events, reactions, propensities)

    # timing = timeit.timeit(bench_gpu, setup=bench_gpu, number=iterations)

    # print('gpu_selected_reaction_propensities took {} s per iteration'.format(
    #     timing/iterations))


if __name__ == '__main__':
    main()
