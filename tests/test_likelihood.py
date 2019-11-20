import unittest

import numpy as np

from gillespie import likelihood, stochastic_sim


class TestAnalyzer(unittest.TestCase):
    def test_reaction_propensities1(self):
        components = np.array([[1, 2]])
        reactions = stochastic_sim.ReactionNetwork(1)
        reactions.k = np.array([1.0], dtype=np.single)
        reactions.reactants = np.array([[0]], dtype=np.int32)
        reaction_events = np.array([0, 0])
        result = likelihood.calculate_selected_reaction_propensities(
            components, reaction_events, reactions
        )
        self.assertListEqual(result.tolist(), [1, 2])

    def test_reaction_propensities2(self):
        components = np.array([[1, 2]])
        reactions = stochastic_sim.ReactionNetwork(2)
        reactions.k = np.array([2.0, 5.0], dtype=np.single)
        reactions.reactants = np.array([[0], [-1]], dtype=np.int32)

        result = likelihood.calculate_sum_of_reaction_propensities(
            components, reactions
        )

        self.assertListEqual(result.tolist(), [5 + 1 * 2, 5 + 2 * 2])

    def test_evaluate_trajectory(self):
        old_ts = np.array([1.0, 2.0, 3.0])
        traj = np.array([10.0, 15.0, 10.0])
        new_ts = np.array([0.0, 1.5, 2.5, 4.0])

        result = likelihood.evaluate_trajectory_at(traj, old_ts, new_ts)
        self.assertListEqual(result.tolist(), [10.0, 10.0, 15.0, 10.0])

    def test_resample_trajectory2(self):
        old_ts = np.array([1.0, 2.0, 3.0])
        traj = np.array([10.0, 15.0, 10.0])
        new_ts = np.array([0.0, 1.5, 2.5, 4.0])

        result = likelihood.evaluate_trajectory_at(traj, old_ts, new_ts)
        self.assertListEqual(result.tolist(), [10.0, 10.0, 15.0, 10.0])

    def test_time_average_trajectory_trivial(self):
        old_ts = np.array([1.0, 2.0, 3.0])
        traj = np.array([10.0, 15.0, 10.0])
        result = likelihood.time_average(traj, old_ts, old_ts)
        self.assertListEqual(result.tolist(), [10.0, 15.0])

    def test_time_average_trajectory(self):
        old_ts = np.array([1.0, 2.0, 3.0])
        traj = np.array([10.0, 15.0, 10.0])
        new_ts = np.array([0.0, 1.5, 2.5, 3.5])

        result = likelihood.time_average(traj, old_ts, new_ts)
        self.assertListEqual(result.tolist(), [10.0, 12.5, 12.5])

        new_ts = np.array([2.0, 2.5, 3.5, 4.0])
        result = likelihood.time_average(traj, old_ts, new_ts)
        self.assertListEqual(result.tolist(), [15.0, 12.5, 10.0])


if __name__ == "__main__":
    unittest.main()
