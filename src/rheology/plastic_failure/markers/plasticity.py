#!/usr/bin/env python3
# -*- coding: utf-8 -*-
""" Plastic failure model.
"""
from numba import njit
from earthbox.utilities.printfuncs import timeit_memit
from earthbox.utilities.initialization_tools import InitializationTools
from earthbox.model import PyModelData
from earthbox.model import ModelData
from .options.options import MarkerPlasticityOptionCollection
from .options.options import MarkerPlasticityOption
# pylint:disable=R0913


class MarkerPlasticity(InitializationTools):
    """ Marker plasticity model.
    """
    options: MarkerPlasticityOptionCollection
    option_active: MarkerPlasticityOption | None

    def __init__(self):
        super().__init__()
        self.options = MarkerPlasticityOptionCollection(
            'MarkerPlasticityTypes')
        self.options.set_options_dict(self.get_options())
        self.option_active = None

    def get_options(self) -> dict[int, MarkerPlasticityOption]:
        """ Build dictionary of basic itype_plasticity id's and options.

        Returns
        -------
        options: dict[int, MarkerPlasticityOption]
            Integer keys equal to itype_plasticity.

        """
        options = {
            0: MarkerPlasticityOption(
                option_name=self.options.option_names.
                viscoelastic_stress_forecast,
                description=(
                    'Viscoelastic plasticity model where stress is forecasted '
                    'using a viscoelastic stress forecast.'
                    ),
                update_func=self.options.update_funcs.viscoelastic_forecast
                ),
            1: MarkerPlasticityOption(
                option_name=self.options.option_names.
                pure_elastic_stress_forecast,
                description=(
                    'Viscoelastic plasticity model where stress is forecasted '
                    'using a purely elastic stress forecast.'
                    ),
                update_func=self.options.update_funcs.elastic_forecast
                )
            }
        return options

    @timeit_memit('Finished initializing marker plasticity', level=1)
    def initialize(
            self,
            pymodel: PyModelData,
            marker_plasticity_model: int | str | None = None
    ) -> None:
        """ Initialize marker plasticity model.
        """
        self.update_option_id_using_stype(
            pymodel, get_stype_from_model, update_option_id)
        option_id = self.get_option_id(
            marker_plasticity_model, pymodel,
            get_option_id_from_model, update_option_id
            )
        self.option_active = self.options.options[option_id]

    def print_option(self, pymodel: PyModelData) -> None:
        """ Print active option.
        """
        option_id = get_option_id_from_model(pymodel)
        self.options.print_option(option_id)


def get_option_id_from_model(pymodel: PyModelData) -> int:
    """ Get itype_plasticity from model.
    """
    return pymodel.materials.parameters.material_description.\
            itype_plasticity.value


def get_stype_from_model(pymodel: PyModelData) -> int:
    """ Get stype_plasticity from model.
    """
    return pymodel.materials.parameters.material_description.\
            stype_plasticity.value


@njit
def update_option_id(
        model: ModelData, # type: ignore
        option_id: int
) -> None:
    """ Update itype_plasticity in model.
    """
    model.materials.parameters.material_description.itype_plasticity.value = (
        option_id)
