"""Repository-tracked entry, discovered by SMAPpy's extra plugin roots."""
from smappy_centerline.plugin import CenterlineRadial as CenterlineRadialPlugin
from smappy.plugins import register


@register("Analysis/Measure/Radial distance and density")
class CenterlineRadial(CenterlineRadialPlugin):
    name = "Radial distance and density"
