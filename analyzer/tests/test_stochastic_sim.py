import unittest

import numpy as np

from analyzer import stochastic_sim


class TestSim(unittest.TestCase):
    def test_select_reaction(self):
        i = stochastic_sim.select_reaction([0.1, 0.2, 0.3])
        self.assertIn(i, [0, 1, 2])

    def test_simulate(self):
        length = 100

        timestamps = np.zeros(length)
        trajectory = np.zeros((1, length))
        reaction_events = np.zeros(length, dtype='i4')
        trajectory[:, 0] = np.array([5.0])
        ext_components = np.full((2, 5), 10.0)
        ext_timestamps = np.array([0.0, 0.5, 1.0, 1.5, 2.0])

        reaction_k = np.array([0.5, 0.8])
        reaction_reactants = np.array([[-1], [0]])
        reaction_products = np.array([[2, -1], [0, 2]])

        stochastic_sim.simulate_one(timestamps, trajectory, ext_components, ext_timestamps,
                                    reaction_k, reaction_reactants, reaction_products, reaction_events)

        self.assertEqual(timestamps[0], 0.0)
        self.assertListEqual(trajectory[:, 0].tolist(), [5.0])

        for i in range(1, length):
            self.assertGreater(timestamps[i], timestamps[i - 1])


if __name__ == '__main__':
    unittest.main()
