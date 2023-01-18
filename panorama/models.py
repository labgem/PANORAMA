#!/usr/bin/env python3
# coding:utf-8

# default libraries
from __future__ import annotations
import logging
from pathlib import Path
from typing import Dict, List, Generator, Set, Tuple, Union
from tqdm import tqdm
import json

suprules_params = ['min_mandatory', 'max_forbidden', 'min_total']
keys_param = suprules_params.append('max_separation')
rule_keys = ['name', 'parameters', 'type']
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
    """Check if all parameters are inside and with the good type

    :param param_dict: Dictionnary with all the parameters for the rule
    :param mandatory_keys: list of the mandatory keys

    :raise KeyError: One or more mandatory parameters are missing
    :raise TypeError: One or more parameters are with a non-acceptable type
    :raise ValueError: One or more parameters are with a non-acceptable value
    :raise Exception: Manage an unexpected error
    """
    try:
        check_key(param_dict, mandatory_keys)
    except KeyError:
        raise KeyError(f"One or more attribute are missing in parameters")
    except Exception as error:
        raise Exception(f"Unexpected Error: {error} to get parameters")
    else:
        for key, value in param_dict.items():
            if key == "max_separation":
                if not isinstance(value, int):
                    raise TypeError("The max_separation value is not an integer")
                if value < 0:
                    raise ValueError("The max_separation value is not positive")
            if key == "min_mandatory":
                if not isinstance(value, int):
                    raise TypeError("The min_mandatory value is not an integer")
                if value < 0:
                    raise ValueError("The min_mandatory value is not positive")
            if key == "max_forbidden":
                if not isinstance(value, int):
                    raise TypeError("The max_forbidden value is not an integer")
                if value < 0:
                    raise ValueError("The max_forbidden value is not positive")
            if key == "min_total":
                if not isinstance(value, int):
                    raise TypeError("The min_total value is not an integer")
                if value < 0:
                    raise ValueError("The min_total value is not positive")


def check_dict(data_dict: Dict[str, Union[str, int, list, Dict[str, int]]], mandatory_keys: List[str],
               param_keys: List[str] = None):
    """Check if all keys and values are present and with good type before to add in model

    :param data_dict: Dictionnary with all model information
    :param mandatory_keys: list of the mandatory keys
    :param param_keys: list of the mandatory keys for parameters

    :raise KeyError: One or more keys are missing or non-acceptable
    :raise TypeError: One or more value are not with good type
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
                    raise TypeError(f"The name value must be str type")
            elif key == 'type':
                if not isinstance(value, str):
                    raise TypeError(f"The type attribute must be a string")
                if value not in accept_type:
                    raise ValueError(f"Accepted type must be in {accept_type}")
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
            elif key == 'relation':
                if value is not None:
                    if not isinstance(value, str):
                        raise TypeError("relation value from family in json must be an str or null")
                    if value not in ['homologs', 'analogs']:
                        raise ValueError("Family relation must be null homologs or analogs")
            elif key == 'duplicate':
                if not isinstance(value, int):
                    raise TypeError("duplicate value from family in json must be an int or null")
                if value < 0:
                    raise ValueError("Dupplication must be positive")
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
        :raise TypeError: One or more value are not with good type
        :raise ValueError: One or more value are not non-acceptable
        :raise Exception: Manage unexpected error
        """
        for file in tqdm(list(models_path.glob("*.json")), unit='model', desc="Read model", disable=disable_bar):
            with open(file.resolve().as_posix()) as json_file:
                data = json.load(json_file)
                model = Model()
                try:
                    model.read_model(data)
                except KeyError:
                    raise KeyError(f"One or more key in {file} are missing.")
                except TypeError:
                    raise TypeError(f"One or more attribute are not with the good type in {file}.")
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


class Model:
    """Represent Model rules which describe biological system

     :param name: Neme of the element
     :param parameters: Dictionary with the parameters name and value
     :param mandatory: Set of mandatory sub element
     :param accessory: Set of accessory sub element
     :param forbidden: Set of forbidden sub element
     :param neutral: Set of neutral sub element
     :param canonical: List of canonical models
     """

    def __init__(self, name: str = "", parameters: dict = None, mandatory: set = None, accessory: set = None,
                 forbidden: set = None, neutral: set = None, canonical: list = None):
        """Constructor Method
        """

        self.name = name
        self.parameters = parameters if parameters is not None else dict()
        self.mandatory = mandatory if mandatory is not None else set()
        self.accessory = accessory if accessory is not None else set()
        self.forbidden = forbidden if forbidden is not None else set()
        self.neutral = neutral if neutral is not None else set()
        self.canonical = canonical if canonical is not None else []

    def __repr__(self):
        return f"name: {self.name}, number of FU: {self.size[0]}, number of families: {self.size[1]}"

    def __str__(self):
        out_dic = {"name": self.name}
        out_dic.update(
            {k: v for k, v in self.parameters.items()} if self.parameters is not None else {"parameters": None})
        return ", ".join([f"{k}: {v}" for k, v in out_dic.items()]) + "\n\t" + \
            "\n\t".join(fu.__str__(ident="\n\t\t") for fu in self.func_units)

    @property
    def size(self) -> Tuple[int, int]:
        """Get number of elements in model

        :return: number of functional unit and number of families
        """
        return len(list(self.func_units)), len(list(self.families))

    @property
    def func_units(self) -> Generator[FuncUnit, None, None]:
        """ Access to all functional units in models

        :return: A generator with all functional units
        """
        for func_unit in self.mandatory.union(self.accessory, self.forbidden, self.neutral):
            yield func_unit

    @property
    def families(self) -> Generator[Family, None, None]:
        """Access to all families in models

        :return: A generator with all families
        """
        for func_unit in self.func_units:
            yield from func_unit.families

    def duplicate_fu(self, filter_type: str = None):
        """Access to all families that are duplicated in functional unit

        :return: A generator with all families
        """
        assert filter_type in [None, 'mandatory', 'forbidden']
        if filter_type is None:
            select_fu = self.func_units
        elif filter_type == "mandatory":
            select_fu = self.mandatory
        elif filter_type == "forbidden":
            select_fu = self.forbidden
        else:
            raise Exception("Unexpected error")
        for fam in select_fu:
            if fam.duplicate >= 1:
                yield fam

    def check_model(self):
        """Check model consistency

        :raise Exception: Model is not consistent
        """
        if self.parameters["min_mandatory"] > len(self.mandatory) + len(list(self.duplicate_fu("mandatory"))):
            raise Exception(f"There is less mandatory functional units than the minimum mandatory in {self.name}")
        if self.parameters["max_forbidden"] > len(self.forbidden) + len(list(self.duplicate_fu("forbidden"))):
            raise Exception(f"There is less forbidden functional units than the maximum forbidden accepted in {self.name}")
        if self.parameters["min_total"] > self.size[0] + len(list(self.duplicate_fu())):
            raise Exception(f"There is less functional units than the minimum total"
                            f" of functional unit needed in {self.name}")
        if len(self.neutral) == 1 and self.size[0] == 1:
            logging.getLogger().warning(f"There is only one family in {self.name}. "
                                        f"It should be a mandatory type.")

    def read_model(self, data_model: dict):
        """Read model to parse in self attributes

        :param data_model: json data dictionary
        """
        mandatory_key = ['name', 'parameters', 'func_units']
        param_keys = ['max_separation', 'min_mandatory', 'max_forbidden', 'min_total']

        check_dict(data_model, mandatory_keys=mandatory_key, param_keys=param_keys)

        self.name = data_model["name"]
        self.parameters = data_model["parameters"]
        for dict_fu in data_model["func_units"]:
            f_unit = FuncUnit()
            f_unit.model = self
            f_unit.read_func_unit(dict_fu)
            self.add_func_unit(f_unit)
        if 'canonical' in data_model and data_model["canonical"] is not None:
            self.canonical = data_model["canonical"]

    def add_func_unit(self, func_unit: FuncUnit):
        """Add a function unit in model

        :param func_unit: function unit
        """
        func_unit.check_func_unit()
        if func_unit.type == "mandatory":
            self.mandatory.add(func_unit)
        elif func_unit.type == "accessory":
            self.accessory.add(func_unit)
        elif func_unit.type == "forbidden":
            self.forbidden.add(func_unit)
        else:
            self.neutral.add(func_unit)


class FuncUnit:
    """Represent functional unit definition rule

    :param name: Neme of the element
    :param parameters: Dictionary with the parameters name and value
    :param mtype: Type of the rule (mandatory, accessory, forbidden or neutral)
    :param mandatory: Set of mandatory sub element
    :param accessory: Set of accessory sub element
    :param forbidden: Set of forbidden sub element
    :param neutral: Set of neutral sub element
    """

    def __init__(self, name: str = "", parameters: dict = None, mtype: str = "", duplicate: int = 0,
                 mandatory: set = None, accessory: set = None, forbidden: set = None, neutral: set = None,
                 model: Model = None):
        """Constructor Method
        """
        self.name = name
        self.parameters = parameters if parameters is not None else dict()
        self.type = mtype
        self.duplicate = duplicate
        self.mandatory = mandatory if mandatory is not None else set()
        self.accessory = accessory if accessory is not None else set()
        self.forbidden = forbidden if forbidden is not None else set()
        self.neutral = neutral if neutral is not None else set()
        self.model = model

    def __repr__(self):
        return f"name: {self.name}, from: {self.model.name}, number of families: {self.size}"

    def __str__(self, ident: str = "\n\t"):
        out_dic = {"name": self.name, "type": self.type,
                   "model": self.model.name if self.model is not None else None}
        out_dic.update(
            {k: v for k, v in self.parameters.items()} if self.parameters is not None else {"parameters": None})
        return ", ".join([f"{k}: {v}" for k, v in out_dic.items()]) + ident + \
            ident.join(str(fam) for fam in self.families)

    @property
    def size(self) -> int:
        """Get number of families in model

        :return: Number of families
        """
        return len(list(self.families))

    @property
    def families(self) -> Generator[Family, None, None]:
        """Access to all families in functional unit

        :return: A generator with all families
        """
        for fam in self.mandatory.union(self.accessory, self.forbidden, self.neutral):
            yield fam

    def duplicate_fam(self, filter_type: str = None):
        """Access to all families that are duplicated in functional unit

        :return: A generator with all families
        """
        assert filter_type in [None, 'mandatory', 'forbidden']
        if filter_type is None:
            select_fam = self.families
        elif filter_type == "mandatory":
            select_fam = self.mandatory
        elif filter_type == "forbidden":
            select_fam = self.forbidden
        else:
            raise Exception("Unexpected error")
        for fam in select_fam:
            if fam.duplicate >= 1:
                yield fam

    def check_func_unit(self):
        """Check functional unit consistency

        :raise Exception: Model is not consistent
        """
        if self.parameters["min_mandatory"] > len(self.mandatory) + len(list(
                self.duplicate_fam(filter_type="mandatory"))):
            raise Exception(f"There is less mandatory families than the minimum mandatory in {self.name}")
        if self.parameters["max_forbidden"] > len(self.forbidden):
            raise Exception(f"There is less forbidden families than the maximum forbidden accepted in {self.name}")
        if self.parameters["min_total"] > self.size + len(list(self.duplicate_fam())):
            raise Exception(f"There is less families than the minimum total"
                            f" of functional unit needed in {self.name}")
        if len(self.neutral) == 1 and self.size == 1:
            logging.getLogger().warning(f"There is only one family in {self.name}. "
                                        f"It should be a mandatory type.")

    def check_parameters(self):
        """Check parameters consistency

        :raise Exception: Model is not consistent
        """
        if "min_mandatory" not in self.parameters:
            self.parameters["min_mandatory"] = self.model.parameters["min_mandatory"]
            logging.getLogger().warning(f"No min_mandatory parameter given in functional unit: {self.name}. "
                                        f"min_mandatory from system {self.model.name} will be used")
        if "max_forbidden" not in self.parameters:
            self.parameters["max_forbidden"] = self.model.parameters["max_forbidden"]
            logging.getLogger().warning(f"No max_forbidden parameter given in functional unit: {self.name}. "
                                        f"max_forbidden from system {self.model.name} will be used")
        if "min_total" not in self.parameters:
            self.parameters["min_total"] = self.model.parameters["min_total"]
            logging.getLogger().warning(f"No min_total parameter given in functional unit: {self.name}. "
                                        f"min_total from system {self.model.name} will be used")
        if "max_separation" not in self.parameters:
            self.parameters["max_separation"] = self.model.parameters["max_separation"]
            logging.getLogger().warning(f"No max_separation parameter given in functional unit: {self.name}. "
                                        f"max_separation from system {self.model.name} will be used")

    def read_func_unit(self, data_fu: dict):
        """Read functional unit

        :param data_fu: data json file of all function units
        """
        mandatory_key = ['name', 'families', 'type']

        check_dict(data_fu, mandatory_keys=mandatory_key)

        self.name = data_fu["name"]
        self.type = data_fu["type"]
        if "parameters" in data_fu and data_fu["parameters"] is not None:
            self.parameters = data_fu["parameters"]
        self.check_parameters()
        for fam_dict in data_fu["families"]:
            family = Family(func_unit=self)
            family.read_family(fam_dict)
            self.add_family(family)

    def add_family(self, family: Family):
        """Add family in families functional unit

        :param family: a family
        """
        if family.type == 'mandatory':
            self.mandatory.add(family)
        elif family.type == "accessory":
            self.accessory.add(family)
        elif family.type == "forbidden":
            self.forbidden.add(family)
        else:
            self.neutral.add(family)


class Family:
    """Represent family model definition rule

    :param name: Neme of the element
    :param parameters: Dictionary with the parameters name and value
    :param mtype: Type of the rule (mandatory, accessory, forbidden or neutral)
    :param relation:
    :param func_unit: Functional unit in which is the family
    """

    def __init__(self, name: str = "", parameters: dict = None, mtype: str = "", relation: str = "",
                 duplicate: int = 0, func_unit: FuncUnit = None):
        """Constructor Method
        """
        self.name = name
        self.parameters = parameters if parameters is not None else dict()
        self.type = mtype
        self.duplicate = duplicate
        self.relation = relation
        self.func_unit = func_unit

    def __repr__(self):
        return f"Family name : {self.name}, type : {self.type}, " \
               f"Functional Unit : {self.func_unit.name if self.func_unit is not None else 'Unknow'}, " \
               f"Model : {self.model.name if self.model is not None else 'Unknow'}"

    def __str__(self):
        out_dic = {"name": self.name, "type": self.type, "relation": self.relation,
                   "functional unit": self.func_unit.name if self.func_unit is not None else None,
                   "model": self.model.name if self.model is not None else None}
        out_dic.update(
            {k: v for k, v in self.parameters.items()} if self.parameters is not None else {"parameters": None})
        return ", ".join([f"{k}: {v}" for k, v in out_dic.items()])

    @property
    def model(self) -> Model:
        """Get the model in which is family"""
        return self.func_unit.model

    def check_parameters(self):
        """Check parameters consistency

        :raise Exception: Model is not consistent
        """
        if "max_separation" not in self.parameters:
            self.parameters["max_separation"] = self.model.parameters["max_separation"]
            logging.getLogger().warning(f"No max_separation parameter given in family: {self.name}."
                                        f"max_separation from functional unit: {self.func_unit.name} will be used")

    def read_family(self, data_fam: dict):
        """Read family

        :param data_fam: data json file with families
        """
        check_dict(data_fam, mandatory_keys=['name', 'type'])
        self.name = data_fam['name']
        self.type = data_fam['type']
        if "parameters" in data_fam and data_fam["parameters"] is not None:
            self.parameters = data_fam["parameters"]
        self.check_parameters()
        if 'relation' in data_fam and data_fam["relation"] is not None:
            self.relation = data_fam['relation']
        if 'duplicate' in data_fam and data_fam["duplicate"] is not None:
            self.duplicate = data_fam['duplicate']
