"""Lazy Nikon ND2 input for single-channel, single-plane time series.

Readers are scoped to each iterator/call, so previews and the fitting read-ahead
thread never share a file handle. Pixels are neither scaled nor reoriented.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from operator import index as integer
from pathlib import Path
import warnings

import numpy as np

from .tiff import ImageSource


def is_nd2(path):
    """Recognise a standalone ND2 path without importing the optional reader."""
    path = Path(path)
    return path.suffix.lower() == ".nd2" and not path.is_dir()


def reject_live_nd2(path):
    if is_nd2(path):
        raise ValueError("Live ND2 reading is not supported; use a completed ND2 file")


def _reader(path):
    try:
        from nd2 import ND2File
    except ImportError as exc:
        raise ImportError("ND2 input requires nd2>=0.10; install it with "
                          "your SMAPpy Python: python -m pip install 'nd2>=0.10'") from exc
    return ND2File(path)


@dataclass
class ND2Source(ImageSource):
    # Already normalized camera fields; do not interpret Nikon gain fields as
    # Micro-Manager camera settings. Unknown camera parameters stay manual.
    camera_values: dict = field(default_factory=dict)

    def frames(self, chunk=100, start=0, stop=None):
        chunk, start = integer(chunk), integer(start)
        stop = self.n_frames if stop is None else integer(stop)
        if chunk < 1 or start < 0 or stop < 0:
            raise ValueError("chunk must be positive; start and stop must be non-negative")
        stop = min(stop, self.n_frames)
        if start >= stop:
            return
        with _reader(self.files[0]) as reader:
            for first in range(start, stop, chunk):
                block = np.empty((min(chunk, stop - first), *self.shape), self.dtype)
                for j in range(len(block)):
                    block[j] = self._read(reader, first + j)
                yield first, block

    def _read(self, reader, index):
        frame = reader.read_frame(index)
        # nd2 squeezes singleton spatial axes too; restore Y,X explicitly.
        if frame.size != np.prod(self.shape) or frame.dtype != self.dtype:
            raise ValueError(f"ND2 frame {index} has inconsistent shape or dtype")
        return frame.reshape(self.shape)

    def frame(self, index):
        index = integer(index)
        if not 0 <= index < self.n_frames:
            raise IndexError(f"frame {index} outside stack ({self.n_frames})")
        with _reader(self.files[0]) as reader:
            # Own the memory even if the library uses a view of the file.
            return self._read(reader, index).copy()

    def watch(self, chunk=100, **kwargs):
        raise ValueError("Live ND2 reading is not supported; use a completed ND2 file")


def open_nd2(path):
    path = Path(path)
    with _reader(path) as reader:
        sizes = dict(reader.sizes)
        if reader.is_rgb or any(n != 1 for axis, n in sizes.items()
                                if axis not in ("T", "Y", "X")):
            raise ValueError(f"ND2 supports only single-channel YX or TYX data; got {sizes}")
        attrs = reader.attributes
        shape = (int(attrs.heightPx), int(attrs.widthPx))
        count = int(attrs.sequenceCount)
        if min(shape) < 1 or count < 1 or count != sizes.get("T", 1):
            raise ValueError(f"ND2 has empty or inconsistent time dimensions: {sizes}")
        dtype = np.dtype(reader.dtype)
        summary = {"Format": "Nikon ND2", "ND2-Sizes": sizes,
                   "ND2-BitsPerComponentSignificant": attrs.bitsPerComponentSignificant}
        values = {}
        try:
            channels = reader.metadata.channels
            if channels:
                channel = channels[0]
                for key, value in (
                    ("ND2-Channel", getattr(getattr(channel, "channel", None), "name", None)),
                    ("ND2-Objective", getattr(getattr(channel, "microscope", None), "objectiveName", None)),
                ):
                    if value is not None:
                        summary[key] = value
                # voxel_size() may return 1 for uncalibrated axes. Only use
                # values explicitly marked calibrated, never that fallback.
                calibrated = channel.volume.axesCalibrated[:2]
                if len(calibrated) == 2 and all(calibrated):
                    voxel = reader.voxel_size()
                    xy = (float(voxel.x), float(voxel.y))
                    if all(np.isfinite(v) and v > 0 for v in xy):
                        values["pixelsize_um"] = xy[0] if xy[0] == xy[1] else xy
                        summary["ND2-PixelSizeUmXY"] = xy
        except (AttributeError, TypeError, ValueError, IndexError) as exc:
            warnings.warn(f"ND2 optional metadata unavailable: {exc}", RuntimeWarning)
    return ND2Source(files=[path], shape=shape, dtype=dtype, n_frames=count,
                     n_frames_declared=count, mm_metadata={}, summary=summary,
                     camera_values=values)
