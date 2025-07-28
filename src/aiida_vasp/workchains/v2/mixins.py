"""
Some convenience mixins
"""

# ruff: noqa: PLC0415
from __future__ import annotations

from typing import Any


class WithBuilderUpdater:
    @classmethod
    def get_builder_updater(cls, *args: Any, **kwargs: Any) -> Any:
        """
        Return the corresponding builder updater class for the workchain.

        The arguments are passed directly to the underling `BuilderUpdater` constructor.
        """
        from .common.builder_updater import (
            VaspBandUpdater,
            VaspBuilderUpdater,
            VaspConvUpdater,
            VaspHybridBandUpdater,
            VaspMultiStageRelaxUpdater,
            VaspNEBUpdater,
            VaspRelaxUpdater,
        )

        if cls.__name__ == 'VaspWorkChain':
            return VaspBuilderUpdater(*args, **kwargs)
        elif cls.__name__ == 'VaspRelaxWorkChain':
            return VaspRelaxUpdater(*args, **kwargs)
        elif cls.__name__ == 'VaspMultiStageRelaxWorkChain':
            return VaspMultiStageRelaxUpdater(*args, **kwargs)
        elif cls.__name__ == 'VaspBandsWorkChain':
            return VaspBandUpdater(*args, **kwargs)
        elif cls.__name__ == 'VaspHybridBandsWorkChain':
            return VaspHybridBandUpdater(*args, **kwargs)
        elif cls.__name__ == 'VaspConvergenceWorkChain':
            return VaspConvUpdater(*args, **kwargs)
        elif cls.__name__ == 'VaspNEBWorkChain':
            return VaspNEBUpdater(*args, **kwargs)
        raise NotImplementedError('No builder updater found for workchain {}'.format(cls.__name__))
