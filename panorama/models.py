#!/usr/bin/env python3
# coding:utf-8

# default libraries
from __future__ import annotations
import logging
from pathlib import Path
from typing import Dict, List, Generator, Set, Tuple, Union
from tqdm import tqdm
import json

# TODO Try to add those variable in class
suprules_params = ['min_mandatory', 'max_forbidden', 'min_total', "max_mandatory", "max_total"]
keys_param = suprules_params.append('max_separation')
rule_keys = ['name', 'parameters', 'presence']
accept_type = ['mandatory', 'accessory', 'forbidden', 'neutral']


def check_key(parameter_data: Dict, need_key: List):
    """
    Global function to check if all the keys are present in the dict.
    This function is applied to check model, functional unit and family consistency

    :param parameter_data: Dictonnary which define model, functional unit or family
    :param need_key: List of all the key needed
    """
    if not all([key in parameter_data for key in need_key]):
        raise KeyError(f"All the following key are necessary : {need_key}")


def check_parameters(param_dict: Dict[str, int], mandatory_keys: List[str]):
    """Check if all parameters are inside and with the good presence

    :param param_dict: Dictionnary with all the parameters for the rule
    :param mandatory_keys: list of the mandatory keys

    :raise KeyError: One or more mandatory parameters are missing
    :raise TypeError: One or more parameters are with a non-acceptable presence
    :raise ValueError: One or more parameters are with a non-acceptable value
    :raise Exception: Manage an unexpected error
    """
    try:
        check_key(param_dict, mandatory_keys)
    except KeyError:
        raise KeyError("One or more attribute are missing in parameters")
    except Exception as error:
        raise Exception(f"Unexpected Error: {error} to get parameters")
    else:
        for key, value in param_dict.items():
            if key == "max_separation":
                if not isinstance(value, int):
                    raise TypeError("The max_separation value is not an integer")
                if value < -1:
                    raise ValueError("The max_separation value must be positive or "
                                     "egal to -1 to indicate an undefinied distance.")
            elif key in ["min_mandatory", "min_total", "max_mandatory", "max_total"]:
                if not isinstance(value, int):
                    raise TypeError(f"The {key} value is not an integer")
                if value < -1:
                    raise ValueError(f"The {key} value is not positive. "
                                     "You can also use -1 value to skip the check on this poarameter.")
            elif key == "max_forbidden":
                if not isinstance(value, int):
                    raise TypeError("The max_forbidden value is not an integer")
                if value < 0:
                    raise ValueError("The max_forbidden value is not positive")
            elif key == 'duplicate':
                if not isinstance(value, int):
                    raise TypeError("duplicate value from family in json must be an int")
                if value < 0:
                    raise ValueError("Dupplication must be positive")
            elif key == 'overlap':
                if not isinstance(value, bool):
                    raise TypeError("Overlap value from family in json must be a boolean")
            else:
                raise KeyError(f"{key} is not an acceptable attribute in parameters")


def check_dict(data_dict: Dict[str, Union[str, int, list, Dict[str, int]]], mandatory_keys: List[str],
               param_keys: List[str] = None):
    """Check if all keys and values are present and with good presence before to add in model

    :param data_dict: Dictionnary with all model information
    :param mandatory_keys: list of the mandatory keys
    :param param_keys: list of the mandatory keys for parameters

    :raise KeyError: One or more keys are missing or non-acceptable
    :raise TypeError: One or more value are not with good presence
    :raise ValueError: One or more value are not non-acceptable
    :raise Exception: Manage unexpected error
    """
    param_keys = [] if param_keys is None else param_keys

    try:
        check_key(data_dict, mandatory_keys)
    except KeyError:
        raise KeyError("One or more keys are missing")
    except Exception as error:
        raise Exception(f"Unexpected Error: {error}")
    else:
        for key, value in data_dict.items():
            if key == 'name':
                if not isinstance(value, str):
                    raise TypeError("The name value must be str presence")
            elif key == 'presence':
                if not isinstance(value, str):
                    raise TypeError("The presence attribute must be a string")
                if value not in accept_type:
                    raise ValueError(f"Accepted presence must be in {accept_type}")
            elif key == 'parameters':
                if value is not None:
                    check_parameters(value, param_keys)
            elif key == 'func_units':
                if not isinstance(value, list):
                    raise TypeError("func_unit value in json must be an array")
                if len(value) < 1:
                    raise ValueError("Model need at least one functional unit")
            elif key == 'families':
                if not isinstance(value, list):
                    raise TypeError("families value in json must be an array")
                if len(value) < 1:
                    raise ValueError("Functional unit need at least one families")
            elif key == 'canonical':
                if not isinstance(value, list):
                    raise TypeError("canonical value in json must be an array")
            elif key == 'exchangeables':
                if not isinstance(value, list):
                    raise TypeError("exchangeables value from family in json must be a list")
                if not all(isinstance(elem, str) for elem in value):
                    raise ValueError("Exchangeables families must be a string")
            else:
                raise KeyError(f"{key} is not an acceptable attribute")


class Models:
    """
    :param models: A set of model defining system
    """

    def __init__(self, models: Set[Model] = None):
        """Constructor Method
        """
        self._model_getter = models if models is not None else {}

    def __iter__(self) -> Generator[Model, None, None]:
        for _model in self._model_getter.values():
            yield _model

    @property
    def value(self) -> List[Model]:
        """Return all models added. Useful if you need a list ang not a generator"""
        return list(self)

    @property
    def size(self) -> int:
        """Get the number of model added

        :return: number of model inside
        """
        return len(self.value)

    @property
    def func_units(self) -> Dict[str, Model]:
        """Get all functional units in models link to referent model

        :return: Dictionnary of function unit link to referent model
        """
        func_unit_dict = {}
        for model in self:
            for fu in model.func_units:
                if fu.name not in func_unit_dict:
                    func_unit_dict[fu.name] = {model}
                else:
                    func_unit_dict[fu.name].add(model)
        return func_unit_dict

    @property
    def families(self) -> Dict[str, Model]:
        """Get all families in models link to referent model

        :return: Dictionnary of families link to referent model
        """
        families_dict = {}
        for model in self.__iter__():
            for family in model.families:
                if family.name not in families_dict:
                    families_dict[family.name] = [model]
                else:
                    families_dict[family.name].append(model)
        return families_dict

    def read(self, models_path: Path, disable_bar: bool = False) -> Models:
        """Read all json files models in the directory

        :param models_path: path of models directory
        :param disable_bar: Disable progress bar

        :raise KeyError: One or more keys are missing or non-acceptable
        :raise TypeError: One or more value are not with good presence
        :raise ValueError: One or more value are not non-acceptable
        :raise Exception: Manage unexpected error
        """
        for file in tqdm(list(models_path.glob("*.json")), unit='model', desc="Read model", disable=disable_bar):
            with open(file.resolve().as_posix()) as json_file:
                data = json.load(json_file)
                try:
                    model = Model.read_model(data)
                except KeyError:
                    raise KeyError(f"Problem with one or more key in {file} are missing.")
                except TypeError:
                    raise TypeError(f"One or more attribute are not with the good presence in {file}.")
                except ValueError:
                    raise ValueError(f"One or more attribute are not with an acceptable value in {file}.")
                except Exception:
                    raise Exception(f"Unexpected problem to read json {file}")
                else:
                    self.add_model(model)

    def get_model(self, name: str) -> Model:
        """
        Get a model by his name

        :param name: name to find

        :raise KeyError: Model not present
        """
        try:
            model = self._model_getter[name]
        except KeyError:
            raise KeyError("Model not present in set of value")
        else:
            return model

    def add_model(self, model: Model):
        """
        Add model

        :param model: Complete model object

        :raise Exception: A model with the same name is already present in system
        """
        try:
            self.get_model(model.name)
        except KeyError:
            model.check_model()
            self._model_getter[model.name] = model
        else:
            raise Exception(f"Model {model.name} already in set of value")


class _BasicFeatures:

    def __init__(self, name: str = "", max_separation: int = 0):
        """Constructor Method
        """
        self.name = name
        self.max_separation = max_separation

    def __repr__(self):
        return f"{self.__class__} name : {self.name}"

    def __str__(self):
        return f"{self.__class__} name : {self.name}"

    def read_parameters(self, parameters: Dict[str, Union[str, int, bool]], mandatory: List[str] = None,
                        default: Dict[str, Union[str, int, bool]] = None):
        """Check parameters consistency

        :raise Exception: Model is not consistent
        """
        if mandatory is not None:
            for param in mandatory:
                if param in parameters:
                    self.__setattr__(param, parameters[param])
                else:
                    try:
                        if hasattr(self, '_parent'):
                            parent = self.__getattribute__('_parent')
                            if hasattr(parent, param):
                                self.__setattr__(param, parent.__getattribute__(param))
                            else:
                                raise AttributeError
                        else:
                            raise AttributeError
                    except AttributeError:
                        raise AttributeError(f"There is no value set for {param}")
        if default is not None:
            for param in default.keys():
                if param in parameters:
                    self.__setattr__(param, parameters[param])
                else:
                    if hasattr(self, '_parent'):
                        parent = self.__getattribute__('_parent')
                        if hasattr(parent, param):
                            self.__setattr__(param, parent.__getattribute__(param))
                        else:
                            self.__setattr__(param, default[param])
                    else:
                        self.__setattr__(param, default[param])


class _FuFamFeatures:

    def __init__(self, presence: str = "", parent: Union[FuncUnit, Model] = None, duplicate: int = 0):
        """Constructor Method
        """
        self.presence = presence
        self.duplicate = duplicate
        self._parent = parent


class _ModFuFeatures:
    def __init__(self, mandatory: Set[FuncUnit, Family] = None, min_mandatory: int = 1, max_mandatory: int = None,
                 accessory: Set[FuncUnit, Family] = None, min_total: int = 1, max_total: int = None,
                 neutral: Set[FuncUnit, Family] = None, forbidden: Set[FuncUnit, Family] = None, max_forbidden: int = 0):
        """Constructor Method
        """
        self.mandatory = mandatory if mandatory is not None else set()
        self.min_mandatory = min_mandatory
        self.max_mandatory = max_mandatory if max_mandatory is not None else len(self.mandatory)
        self.accessory = accessory if accessory is not None else set()
        self.min_total = min_total
        self.max_total = max_total if max_total is not None else len(self.mandatory) + len(self.accessory)
        self.max_forbidden = max_forbidden
        self.forbidden = forbidden if forbidden is not None else set()
        self.neutral = neutral if neutral is not None else set()
        self._child_type = "Functional unit" if type(self) == Model else "Family"

    @property
    def _children(self):
        for child in self.mandatory.union(self.accessory, self.forbidden, self.neutral):
            yield child

    def _duplicate(self, filter_type: str = None):
        """Access to all families that are duplicated in functional unit

        :return: A generator with all families
        """
        assert filter_type in [None, 'mandatory', 'accessory', 'forbidden', 'neutral']
        if filter_type is None:
            select_children = self._children
        elif filter_type == "mandatory":
            select_children = self.mandatory
        elif filter_type == "forbidden":
            select_children = self.forbidden
        elif filter_type == "accessory":
            select_children = self.accessory
        elif filter_type == "neutral":
            select_children = self.neutral
        else:
            raise Exception("Unexpected error")
        for child in select_children:
            if child.duplicate >= 1:
                yield child

    def _check(self):
        """Check model consistency

        :raise Exception: Model is not consistent
        """
        if self.min_mandatory > len(self.mandatory) + sum([child.duplicate for child in self._duplicate("mandatory")]):
            raise Exception(f"There is less mandatory {self._child_type} than the minimum mandatory")
        if self.max_forbidden > len(self.forbidden) + sum([child.duplicate for child in self._duplicate("forbidden")]):
            raise Exception(f"There is less forbidden {self._child_type} than the maximum forbidden accepted")
        if self.min_total > len(list(self._children)) + sum([child.duplicate for child in self._duplicate()]):
            raise Exception(f"There is less {self._child_type} than the minimum total")
        if self.min_mandatory > self.min_total:
            raise Exception(f"Minimun mandatory {self._child_type} value is greater than minimum total.")
        if len(self.mandatory) == 0:
            raise Exception(f"There is not mandatory {self._child_type}."
                            f"You should have at least one mandatory {self._child_type} mandatory presence.")

    def add(self, child: Union[FuncUnit, Family]):
        """Add a function unit in model

        :param child: function unit
        """
        if isinstance(child, FuncUnit):
            child.check_func_unit()
        if child.presence == "mandatory":
            self.mandatory.add(child)
        elif child.presence == "accessory":
            self.accessory.add(child)
        elif child.presence == "forbidden":
            self.forbidden.add(child)
        else:
            self.neutral.add(child)


class Model(_BasicFeatures, _ModFuFeatures):
    """Represent Model rules which describe biological system

     :param name: Neme of the element
     :param mandatory: Set of mandatory sub element
     :param accessory: Set of accessory sub element
     :param forbidden: Set of forbidden sub element
     :param neutral: Set of neutral sub element
     :param canonical: List of canonical models
     """

    def __init__(self, name: str = "", mandatory: set = None, min_mandatory: int = 1, accessory: set = None,
                 neutral: set = None, min_total: int = 1, forbidden: set = None, max_forbidden: int = 0,
                 max_separation: int = 0, canonical: list = None):
        """Constructor Method
        """

        super().__init__(name=name, max_separation=max_separation)
        super(_BasicFeatures, self).__init__(mandatory=mandatory, min_mandatory=min_mandatory,
                                             accessory=accessory, neutral=neutral, min_total=min_total,
                                             forbidden=forbidden, max_forbidden=max_forbidden)
        self.canonical = canonical if canonical is not None else []

    @property
    def func_units(self) -> Generator[FuncUnit, None, None]:
        """ Access to all functional units in models

        :return: A generator with all functional units
        """
        yield from self._children

    @property
    def families(self) -> Generator[Family, None, None]:
        """Access to all families in models

        :return: A generator with all families
        """
        for func_unit in self.func_units:
            yield from func_unit.families

    @property
    def size(self) -> Tuple[int, int]:
        """Get number of elements in model

        :return: number of functional unit and number of families
        """
        return len(list(self.func_units)), len(list(self.families))

    def duplicate_fu(self, filter_type: str = None):
        """Access to all families that are duplicated in functional unit

        :return: A generator with all families
        """
        yield from self._duplicate(filter_type)

    def check_model(self):
        """Check model consistency

        :raise Exception: Model is not consistent
        """
        try:
            self._check()
        except Exception:
            raise Exception(f"Consistency not respected  in {self.name}")

    def _read(self, data_model: dict):
        mandatory_key = ['name', 'parameters', 'func_units']
        param_mandatory = ['max_separation', 'min_mandatory', 'max_forbidden', 'min_total', "max_mandatory", "max_total"]

        check_dict(data_model, mandatory_keys=mandatory_key, param_keys=param_mandatory)

        self.name = data_model["name"]
        self.read_parameters(data_model["parameters"], mandatory=param_mandatory)
        for dict_fu in data_model["func_units"]:
            f_unit = FuncUnit.read_func_unit(dict_fu)
            f_unit.model = self
            self.add(f_unit)
        if 'canonical' in data_model and data_model["canonical"] is not None:
            self.canonical = data_model["canonical"]

    @staticmethod
    def read_model(data_model: dict) -> Model:
        """Read model to parse in self attributes

        :param data_model: json data dictionary
        """
        model = Model()
        model._read(data_model)
        return model


class FuncUnit(_BasicFeatures, _FuFamFeatures, _ModFuFeatures):
    """Represent functional unit definition rule

    :param name: Neme of the element
    :param presence: Type of the rule (mandatory, accessory, forbidden or neutral)
    :param mandatory: Set of mandatory sub element
    :param accessory: Set of accessory sub element
    :param forbidden: Set of forbidden sub element
    :param neutral: Set of neutral sub element
    """

    def __init__(self, name: str = "", presence: str = "", mandatory: set = None, min_mandatory: int = 1,
                 accessory: set = None, neutral: set = None, min_total: int = 1, max_separation: int = 0,
                 forbidden: set = None, max_forbidden: int = 0, duplicate: int = 0, model: Model = None):
        """Constructor Method
        """
        super().__init__(name=name, max_separation=max_separation)
        super(_BasicFeatures, self).__init__(presence=presence, duplicate=duplicate, parent=model)
        super(_FuFamFeatures, self).__init__(mandatory=mandatory, min_mandatory=min_mandatory,
                                             accessory=accessory, neutral=neutral, min_total=min_total,
                                             forbidden=forbidden, max_forbidden=max_forbidden)

    @property
    def model(self) -> Model:
        return self._parent

    @model.setter
    def model(self, model: Model):
        self._parent = model

    @model.deleter
    def model(self):
        del self._parent

    @property
    def families(self) -> Generator[Family, None, None]:
        """Access to all families in functional unit

        :return: A generator with all families
        """
        yield from self._children

    @property
    def size(self) -> int:
        """Get number of families in model

        :return: Number of families
        """
        return len(list(self.families))

    def duplicate_fam(self, filter_type: str = None):
        """Access to all families that are duplicated in functional unit

        :return: A generator with all families
        """
        yield from self._duplicate(filter_type)

    def check_func_unit(self):
        """Check functional unit consistency

        :raise Exception: Model is not consistent
        """
        try:
            self._check()
        except Exception:
            raise Exception(f"Consistency not respected  in model {self.model.name} at functional unit {self.name}")

    def _read(self, data_fu: dict):
        """Read functional unit

        :param data_fu: data json file of all function units
        """
        mandatory_key = ['name', 'families', 'presence']
        check_dict(data_fu, mandatory_keys=mandatory_key)

        self.name = data_fu["name"]
        self.presence = data_fu["presence"]
        for fam_dict in data_fu["families"]:
            family = Family.read_family(fam_dict)
            family.func_unit = self
            self.add(family)
        if "parameters" in data_fu:
            param_default = {'duplicate': 0, 'min_total': -1, 'min_mandatory': -1, "max_forbidden": 0,
                             'max_separation': -1, 'max_mandatory': -1, 'max_total': -1}
            self.read_parameters(data_fu["parameters"], default=param_default)

    @staticmethod
    def read_func_unit(data_fu: dict) -> FuncUnit:
        func_unit = FuncUnit()
        func_unit._read(data_fu)
        return func_unit


class Family(_BasicFeatures, _FuFamFeatures):
    """Represent family model definition rule

    :param name: Neme of the element
    :param max_separation: maximum intergene distance
    :param presence: Type of the rule (mandatory, accessory, forbidden or neutral)
    :param func_unit: Functional unit in which is the family
    :param exchangeables: List of exhangeable families
    """

    def __init__(self, name: str = "", max_separation: int = 0, presence: str = "", func_unit: FuncUnit = None,
                 duplicate: int = 0, exchangeables: List[str] = None, overlap: bool = False):
        """Constructor Method
        """
        super().__init__(name=name, max_separation=max_separation)
        super(_BasicFeatures, self).__init__(presence=presence, duplicate=duplicate, parent=func_unit)
        self.exchangeables = exchangeables if exchangeables is not None else []
        self.overlap = overlap

    @property
    def func_unit(self) -> FuncUnit:
        return self._parent

    @func_unit.setter
    def func_unit(self, model: FuncUnit):
        self._parent = model

    @func_unit.deleter
    def func_unit(self):
        del self._parent

    @property
    def model(self) -> Model:
        """Get the model in which is family"""
        return self.func_unit.model

    def _read(self, data_fam: dict):
        param_mandatory = ['max_separation']
        param_default = {"overlap": False, "duplicate": 0}

        check_dict(data_fam, mandatory_keys=['name', 'presence'],
                   param_keys=param_mandatory + list(param_default.keys()))
        self.name = data_fam['name']
        self.presence = data_fam['presence']
        if "parameters" in data_fam:
            self.read_parameters(data_fam["parameters"], mandatory=param_mandatory, default=param_default)
        if 'exchangeables' in data_fam:
            self.exchangeables = data_fam["exchangeables"]

    @staticmethod
    def read_family(data_fam: dict) -> Family:
        """Read family

        :param data_fam: data json file with families
        """
        fam = Family()
        fam._read(data_fam)
        return fam
