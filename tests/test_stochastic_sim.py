import unittest

import numpy as np

from analyzer import stochastic_sim


class TestSim(unittest.TestCase):
    def test_select_reaction(self):
        i = stochastic_sim.select_reaction([0.1, 0.2, 0.3])
        self.assertIn(i, [0, 1, 2])

    def test_update_concentration(self):
        reactions = stochastic_sim.ReactionNetwork(2)
        reactions.reactants = np.array([[-1], [0]], dtype=np.int32)
        reactions.products = np.array([[0], [-1]], dtype=np.int32)

        components = np.array([10.0])

        stochastic_sim.update_components(0, components, reactions)
        self.assertAlmostEqual(components[0], 11.0)
        stochastic_sim.update_components(1, components, reactions)
        self.assertAlmostEqual(components[0], 10.0)

    def test_simulate(self):
        length = 100

        timestamps = np.zeros(length)
        trajectory = np.zeros((1, length), dtype=np.uint16)
        reaction_events = np.zeros(length - 1, dtype="i4")
        trajectory[:, 0] = np.array([5])
        ext_components = np.full((2, 5), 10, dtype=np.uint16)
        ext_timestamps = np.array([0.0, 0.5, 1.0, 1.5, 2.0])

        reactions = stochastic_sim.ReactionNetwork(2)
        reactions.k = np.array([0.5, 0.8], dtype=np.single)
        reactions.reactants = np.array([[-1], [0]], dtype=np.int32)
        reactions.products = np.array([[2, -1], [0, 2]], dtype=np.int32)

        stochastic_sim.simulate_one(
            timestamps,
            trajectory,
            reaction_events,
            reactions,
            ext_timestamps,
            ext_components,
        )

        self.assertEqual(timestamps[0], 0.0)
        self.assertListEqual(trajectory[:, 0].tolist(), [5])

        for i in range(1, length):
            self.assertGreater(timestamps[i], timestamps[i - 1])


if __name__ == "__main__":
    unittest.main()
