configuration = {
    # Path to the Gillespie executable
    'executable': '../CPP/build/Gillespie',
    'gillespie_cwd': '../runs/response'
}

configuration['mutual_information'] = {
    'input_name': 'response.inp',
    'signals': '../data/signal/sig{sig}.traj',
    'responses': '../data/res{sig}-{res}.traj',
    'num_signals': 100,
    'num_responses': 30,
    'signal_duration': 1000, # in seconds
    'signal_component_name': 'S',
    'signal_mean': 1990,
    'signal_correlation_time': 10,
    'signal_diffusion': 10,
    'output': '../data/mutual_information.npz'
}
