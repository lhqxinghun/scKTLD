"""Factory function for Cost classes."""

from .base import BaseCost
from .costkernel import CostRbf


def cost_factory(model, *args, **kwargs):
    for cls in BaseCost.__subclasses__():
        if cls.model == model:
            return cls(*args, **kwargs)
    raise ValueError("Not such model: {}".format(model))
