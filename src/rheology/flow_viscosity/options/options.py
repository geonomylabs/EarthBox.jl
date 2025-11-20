#!/usr/bin/env python3
# -*- coding: utf-8 -*-
""" Classes used to manage marker flow law types.
"""
from typing import Callable
from dataclasses import dataclass
from earthbox.model import ModelData
from earthbox.utilities.option_tools import OptionCollection
from earthbox.utilities.option_tools import Option
from ..utils.update import update_marker_flow_viscosity_blankenbach89
from ..utils.update import update_marker_flow_viscosity_couette_flow_benchmark
from ..utils.update import update_marker_flow_viscosity_isoviscous
from ..utils.update import update_marker_flow_viscosity_x_dependent
from ..utils.composite_viscosity import update_marker_composite_viscosity



UpdateViscosityFuncType = Callable[[ModelData], None]


class UpdateViscosityFuncs:
    """ Class used to define flow viscosity update functions.
    """
    x_coordinate_dependent: UpdateViscosityFuncType
    isoviscous: UpdateViscosityFuncType
    composite: UpdateViscosityFuncType
    temperature_dependent_couette_benchmark: UpdateViscosityFuncType
    temperature_dependent_convection_benchmark: UpdateViscosityFuncType

    def __init__(self):
        self.x_coordinate_dependent = (
            update_marker_flow_viscosity_x_dependent
            )
        self.isoviscous = (
            update_marker_flow_viscosity_isoviscous
            )
        self.composite = (
            update_marker_composite_viscosity
            )
        self.temperature_dependent_couette_benchmark = (
            update_marker_flow_viscosity_couette_flow_benchmark
            )
        self.temperature_dependent_convection_benchmark = (
            update_marker_flow_viscosity_blankenbach89
            )


class FlowLawOption(Option):
    """ Class used to define a flow law option.
    """
    update_viscosity_func: UpdateViscosityFuncType | None

    def __init__(
            self,
            option_name: str | None = None,
            description: str | None = None,
            update_viscosity_func: UpdateViscosityFuncType | None = None,
            flow_type_ids_used_in_function: tuple[int,...] | None = None
    ):
        super().__init__(option_name, description)
        self.update_viscosity_func = update_viscosity_func
        self.flow_type_ids_used_in_function = flow_type_ids_used_in_function

@dataclass
class FlowLawOptionNames:
    """ Class used to define a collection of option names.
    """
    x_coordinate_dependent: str = 'XCoordinateDependent'
    isoviscous: str = 'Isoviscous'
    composite: str = 'Composite'
    temperature_dependent_couette_benchmark: str = (
        'TemperatureDependentCouetteBenchmark')
    temperature_dependent_convection_benchmark: str = (
        'TemperatureDependentConvectionBenchmark')


class FlowLawOptionCollection(OptionCollection):
    """ Class used to define a collection of options.
    """
    options: dict[int, FlowLawOption] | None
    option_names: FlowLawOptionNames
    update_viscosity_funcs: UpdateViscosityFuncs

    def __init__(
        self,
        name: str,
        options: dict[int, FlowLawOption] | None = None
    ):
        super().__init__(name)
        self.options = options
        self.option_names = FlowLawOptionNames()
        self.update_viscosity_funcs = UpdateViscosityFuncs()
