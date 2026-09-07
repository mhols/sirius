import inspect
import os
import pathlib
import nextrapol as nx
from dotenv import load_dotenv

import nextrapol.units as nu
import nextrapol.continuum as continuum


DIR = pathlib.Path(__file__).parent.absolute


class SettingsSiriusEspadons(nx.settings_espadons.SettingsReferenceEspadons):
    IS_REFERENCE = False
    VOIE_METHOD = "OPTIMAL_EXTRACT" #SUM_DIVIDE_CENTRALROW"

    BIG_PSEUDO_FLAT = False # do not use the BIG_PSEUDO_FLAT

    ORDERS = list(range(24, 57))

    CONTINUUM_METHOD_CLASS = continuum.SigmaClippingContinuum

    STORE_PATH = DIR


def get_kwargs():
    """
    convenient shortcut
    """
    return SettingsSiriusEspadons.get_kwargs()