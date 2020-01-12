import unittest

import numpy as np

from gillespie import simulate

reactions1 = {"k": [0.5, 500.0], "reactants": [[0], []], "products": [[], [0]]}
reactions2 = {"k": [1.0, 4.0], "reactants": [[0], [1]], "products": [[0, 1], []]}


class TestSim(unittest.TestCase):
    def test_simulate_simple(self):
        traj = simulate.simulate(10, 100, reactions1, None, initial_value=1000)
        self.assertEqual(traj.components.shape[0], 10)
        self.assertEqual(traj.components.shape[1], 1)
        self.assertEqual(traj.components.shape[2], 100)
        self.assertEqual(traj.components[0, 0, 0], 1000)

    def test_simulate_statistics_correct(self):
        traj = simulate.simulate(1, 1_000_000, reactions1, None, initial_value=1000)
        self.assertAlmostEqual(traj.components[0, 0].mean() / 1000.0, 1.0, delta=1.0)
        self.assertAlmostEqual(traj.components[0, 0].var() / 1000.0, 1.0, delta=1.0)

        responses = simulate.simulate(
            1, 1_000_000, reactions2, traj, initial_value=1000
        )
        self.assertAlmostEqual(responses.components[0, 0].mean() / 250.0, 1.0, delta=1.0)

    def test_simulate_with_ext_trajectory(self):
        signals = simulate.simulate(10, 100, reactions1, None, initial_value=1000)
        responses = simulate.simulate(10, 100, reactions2, signals, initial_value=1000)
        self.assertEqual(responses.components.shape[0], 10)
        self.assertEqual(responses.components.shape[1], 1)
        self.assertEqual(responses.components.shape[2], 100)
        self.assertEqual(responses.components[0, 0, 0], 1000)


if __name__ == "__main__":
    unittest.main()
