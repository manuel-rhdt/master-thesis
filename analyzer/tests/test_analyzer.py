import unittest

import numpy as np

from analyzer import analyzer


class TestAnalyzer(unittest.TestCase):

    def test_reaction_propensities1(self):
        components = np.array([[1, 2]])
        reaction_k = np.array([1])
        reaction_reactants = np.array([[0]])
        reaction_events = np.array([0, 0])
        result = analyzer.calculate_selected_reaction_propensities(components, reaction_k, reaction_reactants, reaction_events)
        self.assertListEqual(result.tolist(), [1, 2])

    def test_reaction_propensities2(self):
        reactions = [{'k': 2, 'reactants': ['C1']}, {'k': 5, 'reactants': []}]

        components = np.array([[1, 2]])
        reaction_k = np.array([2, 5])
        reaction_reactants = np.array([[0], [-1]])

        result = analyzer.calculate_sum_of_reaction_propensities(
            components, reaction_k, reaction_reactants)
        
        self.assertListEqual(result.tolist(), [5 + 1 * 2, 5 + 2 * 2])

    def test_evaluate_trajectory(self):
        old_ts = np.array([1.0, 2.0, 3.0])
        traj = np.array([10.0, 15.0, 10.0])
        new_ts = np.array([0.0, 1.5, 2.5, 4.0])

        result = analyzer.evaluate_trajectory_at(traj, old_ts, new_ts)
        self.assertListEqual(result.tolist(), [10.0, 10.0, 15.0, 10.0])

    def test_resample_trajectory2(self):
        old_ts = np.array([1.0, 2.0, 3.0])
        traj = np.array([10.0, 15.0, 10.0])
        new_ts = np.array([0.0, 1.5, 2.5, 4.0])

        result = analyzer.evaluate_trajectory_at(traj, old_ts, new_ts)
        self.assertListEqual(
            result.tolist(), [10.0, 10.0, 15.0, 10.0])

    def test_time_average_trajectory_trivial(self):
        old_ts = np.array([1.0, 2.0, 3.0])
        traj = np.array([10.0, 15.0, 10.0])
        result = analyzer.time_average(traj, old_ts, old_ts)
        self.assertListEqual(result.tolist(), [10.0, 15.0])

    def test_time_average_trajectory(self):
        old_ts = np.array([1.0, 2.0, 3.0])
        traj = np.array([10.0, 15.0, 10.0])
        new_ts = np.array([0.0, 1.5, 2.5, 3.5])

        result = analyzer.time_average(traj, old_ts, new_ts)
        self.assertListEqual(result.tolist(), [10.0, 12.5, 12.5])

        new_ts = np.array([2.0, 2.5, 3.5, 4.0])
        result = analyzer.time_average(traj, old_ts, new_ts)
        self.assertListEqual(result.tolist(), [15.0, 12.5, 10.0])


if __name__ == '__main__':
    unittest.main()
