import unittest

import numpy as np

from analyzer import analyzer


class TestAnalyzer(unittest.TestCase):

    def test_reaction_propensities1(self):
        reactions = [{'k': 1, 'reactants': ['C1']}]
        result = analyzer.calculate_reaction_propensities(
            reactions, {'C1': [1, 2]})
        self.assertListEqual(result[0].tolist(), [1, 2])

    def test_reaction_propensities2(self):
        reactions = [{'k': 2, 'reactants': ['C1']}, {'k': 5, 'reactants': []}]
        result = analyzer.calculate_reaction_propensities(
            reactions, {'C1': [1, 2]})
        # broadcasting is necessary since the second array in `result` will actually
        # be just a number (namely 5.0)
        result = np.array(np.broadcast_arrays(*result))
        self.assertListEqual(result.tolist(), [[2, 4], [5, 5]])

    def test_resample_trajectory(self):
        old_ts = np.array([1.0, 2.0, 3.0])
        traj = np.array([10.0, 15.0, 10.0])
        new_ts = np.array([0.0, 1.5, 2.5, 4.0])

        result = analyzer.resample_trajectory(traj, old_ts, new_ts)
        self.assertListEqual(result.tolist(), [10.0, 10.0, 15.0, 10.0])

    def test_resample_trajectory2(self):
        old_ts = np.array([1.0, 2.0, 3.0])
        traj = np.array([10.0, 15.0, 10.0])
        new_ts = np.array([0.0, 1.5, 2.5, 4.0])

        result = analyzer.resample_trajectory(traj, old_ts, new_ts)
        self.assertListEqual(
            result.tolist(), [10.0, 10.0, 15.0, 10.0])

    def test_resample_averaged_trajectory_trivial(self):
        old_ts = np.array([1.0, 2.0, 3.0])
        traj = np.array([10.0, 15.0, 10.0])
        result = analyzer.resample_averaged_trajectory(traj, old_ts, old_ts)
        self.assertListEqual(result.tolist(), [10.0, 15.0, 10.0])

    def test_resample_averaged_trajectory(self):
        old_ts = np.array([1.0, 2.0, 3.0])
        traj = np.array([10.0, 15.0, 10.0])
        new_ts = np.array([0.0, 1.5, 2.5, 3.5])

        result = analyzer.resample_averaged_trajectory(traj, old_ts, new_ts)
        self.assertListEqual(result.tolist(), [10.0, 12.5, 12.5, 10.0])

        new_ts = np.array([2.0, 2.5, 3.5, 4.0])
        result = analyzer.resample_averaged_trajectory(traj, old_ts, new_ts)
        self.assertListEqual(result.tolist(), [15.0, 12.5, 10.0, 10.0])


if __name__ == '__main__':
    unittest.main()
