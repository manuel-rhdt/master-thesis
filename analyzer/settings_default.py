configuration = {}

configuration['mutual_information'] = {
    'num_signals': 100,
    'num_responses': 100,

    # there will be (signal_duration * signal_resolution) timestamps in the signal
    'signal_duration': 2500,  # in seconds
    'signal_resolution': 100,  # in 1/second

    'signal_component_name': 'S',
    'signal_mean': 4000,
    'signal_correlation_time': 200,
    'signal_diffusion': 20,
    'output': '$SIM_OUTPUT/run10'
}
