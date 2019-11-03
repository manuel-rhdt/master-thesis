configuration = {
    # Path to the Gillespie executable
    'executable': '../CPP/build/Gillespie',
    'gillespie_cwd': '../runs/response'
}

configuration['mutual_information'] = {
    'input_name': 'response.inp',
    'signals': '../data/signal/sig{sig}.traj',
    'responses': '../data/res{sig}-{res}.traj',
    'num_signals': 30,
    'num_responses': 30,

    # there will be (signal_duration * signal_resolution) timestamps in the signal
    'signal_duration': 3000,  # in seconds
    'signal_resolution': 100,  # in 1/second

    'signal_component_name': 'S',
    'signal_mean': 4000,
    'signal_correlation_time': 100,
    'signal_diffusion': 40,
    'output': '../data/run0/mi'
}
