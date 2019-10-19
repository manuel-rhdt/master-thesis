configuration = {
    # Path to the Gillespie executable
    'executable': '/dev/shm/Git/Gillespie/CPP/cmake-build-release/Gillespie'
}

configuration['mutual_information'] = {
    'signals': '/data/signal/sig{sig}.traj',
    'responses': '/dev/shm/res{sig}-{res}.traj',
    'num_signals': 50,
    'num_responses': 50,
    'signal_component_name': 'S',
    'signal_mean': 199,
    'signal_correlation_time': 10,
    'signal_diffusion': 10,
    'output': '/data/mutual_information.npz'
}
