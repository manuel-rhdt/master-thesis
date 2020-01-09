import dask
import dask.threaded
import dask.array
import numpy
from dask.base import DaskMethodsMixin
from dask.highlevelgraph import HighLevelGraph
from dask.base import tokenize
from dask.utils import funcname
from dask.delayed import Delayed


def _expand_dims(array, dims=2):
    if array is None:
        return None

    array_dim = len(array.shape)
    if array_dim > dims:
        raise Exception("Array has too many dimensions!")

    shape = array.shape
    return array.reshape(((1,) * (dims - array_dim) + shape))


def concatenate(trajectories):
    if all(dask.is_dask_collection(t) for t in trajectories):
        return _concat_distributed(trajectories)
    t = numpy.concatenate([t.timestamps for t in trajectories])
    c = numpy.concatenate([t.components for t in trajectories])
    if not any(t.reaction_events is None for t in trajectories):
        re = numpy.concatenate([t.reaction_events for t in trajectories])
    else:
        re = None
    return TrajectoryArray(t, c, re)


def _concat_distributed(trajectories):
    chunks = (sum((t.chunks[0] for t in trajectories), ()), trajectories[0].chunks[1])
    names = [t.name for t in trajectories]
    name = "concatenate-" + tokenize(names)
    allpartitions = sum((t.__dask_keys__() for t in trajectories), [])
    dsk = {(name, i): n for i, n in enumerate(allpartitions)}
    graph = HighLevelGraph.from_collections(name, dsk, dependencies=trajectories)
    return DistributedTrajectoryArray(graph, name, chunks, meta=None)


def from_dict(dict):
    return TrajectoryArray(
        dict["timestamps"], dict["components"], dict["reaction_events"]
    )


def from_delayed(value, num_traj, length, num_components=1, meta=None, name=None):
    name = name or "from-value-" + tokenize(value, num_traj, length, meta)
    dsk = {(name, 0): value.key}
    graph = HighLevelGraph.from_collections(name, dsk, dependencies=[value])
    chunks = ((num_traj,), (length,))
    return DistributedTrajectoryArray(graph, name, chunks, meta=None)

class TrajectoryArray:
    def __init__(self, timestamps, components, reaction_events=None):
        self.timestamps = _expand_dims(timestamps, 2)
        self.components = _expand_dims(components, 3)
        self.reaction_events = _expand_dims(reaction_events, 2)

        assert self.timestamps.shape[0] == self.components.shape[0]

    def __len__(self):
        return self.timestamps.shape[0]

    def __getitem__(self, key):
        if not (key >= 0 and key < len(self)):
            raise IndexError()
        t = self.timestamps[key]
        c = self.components[key]
        if self.reaction_events is not None:
            re = self.reaction_events[key]
        else:
            re = None
        return TrajectoryArray(t, c, re)


class DistributedTrajectoryArray(DaskMethodsMixin):
    __dask_scheduler__ = staticmethod(dask.threaded.get)

    def __init__(self, dsk, name, chunks, meta):
        if not isinstance(dsk, HighLevelGraph):
            dsk = HighLevelGraph.from_collections(name, dsk, dependencies=[])
        self.dask = dsk
        self.meta = meta
        self._name = name
        self.chunks = chunks

    def __dask_graph__(self):
        return self.dask

    def __dask_keys__(self):
        return [(self._name, i) for i in range(self.npartitions)]

    def __len__(self):
        return sum(self.chunks[0])

    def __getitem__(self, key):
        if not (key >= 0 and key < len(self)):
            raise IndexError()

        # find the chunk to index
        for chunk, chunk_key in zip(self.chunks[0], self.__dask_keys__()):
            if chunk <= key:
                key -= chunk
                assert key >= 0
            else:
                # found the correct chunk. Now create result.
                name = self.name + "-" + str(key)
                dsk = {(name, 0): (TrajectoryArray.__getitem__, chunk_key, key)}
                graph = HighLevelGraph.from_collections(name, dsk, dependencies=[self])
                return DistributedTrajectoryArray(
                    graph, name, ((1,), self.chunks[1]), self.meta
                )

    @property
    def npartitions(self):
        """Return number of partitions"""
        return len(self.chunks[0])

    def _get_array(self, attr, chunks, meta):
        name = self._name + "-" + attr
        dsk = {}
        for i in range(self.npartitions):
            dsk[(name, i) + (0,) * (len(chunks) - 1)] = (getattr, (self._name, i), attr)
        graph = HighLevelGraph.from_collections(name, dsk, dependencies=[self])
        return dask.array.Array(graph, name, chunks, meta=meta)

    @property
    def timestamps(self):
        """Return timestamps as a dask array"""
        meta = numpy.empty((0, 0), dtype=numpy.double)
        return self._get_array("timestamps", self.chunks, meta)

    @property
    def components(self):
        """Return components as a dask array"""
        meta = numpy.empty((0, 0, 0), dtype=numpy.uint16)
        chunks = (self.chunks[0], (1,), self.chunks[1])
        return self._get_array("components", chunks, meta)

    @property
    def reaction_events(self):
        """Return reaction_events as a dask array"""
        meta = numpy.empty((0, 0), dtype=numpy.uint8)
        return self._get_array("reaction_events", self.chunks, meta)

    @property
    def name(self):
        return self._name

    def __dask_layers__(self):
        return (self.name,)

    def __dask_postcompute__(self):
        return (concatenate, ())

    def __dask_postpersist__(self):
        return (DistributedTrajectoryArray, (self.name, self.chunks, self.meta))

    def __dask_tokenize__(self):
        return self.name

    def to_delayed(self):
        keys = self.__dask_keys__()
        graph = self.__dask_graph__()
        return [Delayed(k, graph) for k in keys]
