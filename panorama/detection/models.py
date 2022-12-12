#!/usr/bin/env python3
# coding:utf-8

# default libraries
from __future__ import annotations
import logging
from typing import Dict, List, Generator, Set, Union


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
            self._model_getter[model.name] = model
        else:
            raise Exception(f"Model {model.name} already in set of value")


class Rule:
    """Super class for Model and FuncUnit and Family

    :param name: Name of the element
    :param parameters: Dictionary with the parameters name and value
    :param mtype: Type of the rule (mandatory, accessory, forbidden or neutral)
    :param duplicate: If item is duplicate, number of duplication
    """

    def __init__(self, name: str = "", parameters: dict = None, mtype: str = "", duplicate: int = 0):
        """Constructor Method
        """
        self.name = name
        self.parameters = parameters if parameters is not None else dict()
        self.type = mtype
        self.duplicate = duplicate if duplicate is not None else 0

    def check_parameters(self, param_dict):
        """Check if all parameters are inside and with the good type

        :param param_dict: Dictionnary with all the parameters for the rule

        :raise KeyError: One or more mandatory parameters are missing
        :raise TypeError: One or more parameters are with a non-acceptable type
        :raise ValueError: One or more parameters are with a non-acceptable value
        :raise Exception: Manage an unexpected error
        """
        try:
            max_sep = param_dict['max_separation']
        except KeyError:
            raise KeyError("No max_separation in data rules")
        except Exception:
            raise Exception("Unexpected Error")
        else:
            if not isinstance(max_sep, int):
                raise TypeError("The max_separation value is not an integer")
            if max_sep < 0:
                raise ValueError("The max_separation value is not positive int type")

    def check_rule_key(self, rule_dict, key: str):
        """Check if key is inside and with the good type

        :param rule_dict: Dictionnary with all the parameters for the rule
        :param key: key which must be checked

        :raise KeyError: The given key is non-acceptable
        :raise TypeError: One or more parameters are with a non-acceptable type
        :raise ValueError: One or more parameters are with a non-acceptable value
        """
        if key not in rule_keys:
            raise KeyError(f"Unreadable key : {key}. Accepeted key are {rule_keys}")
        else:
            if key == 'name':
                if not isinstance(rule_dict[key], str):
                    raise TypeError(f"The name value must be str type")
            elif key == 'type':
                if rule_dict[key] not in accept_type:
                    raise ValueError(f"The type value {rule_dict[key]} is incorrect in {type(self)} {self.name}")
            elif key == 'parameters' and rule_dict[key] is not None:
                self.check_parameters(rule_dict[key])

    def read_rule(self, key: str, value: Union[str, dict]):
        """Change value of key attributes

        :param key: Attribute to change
        :param value: New value
        """
        if key == 'name':
            self.name = value
        elif key == 'type':
            self.type = value
        elif key == 'parameters':
            self.parameters = value


class SupRule(Rule):
    """Super class for Model and FuncUnit


    :param name: Name of the element
    :param parameters: Dictionary with the parameters name and value
    :param mtype: Type of the rule (mandatory, accessory, forbidden or neutral)
    :param mandatory: Set of mandatory sub element
    :param accessory: Set of accessory sub element
    :param forbidden: Set of forbidden sub element
    :param neutral: Set of neutral sub element
    """

    def __init__(self, name: str = "", parameters: dict = None, mtype: str = "", mandatory: set = None,
                 accessory: set = None, forbidden: set = None, neutral: set = None):
        """Constructor Method
        """
        super().__init__(name, parameters, mtype)
        self.mandatory = mandatory if mandatory is not None else set()
        self.accessory = accessory if accessory is not None else set()
        self.forbidden = forbidden if forbidden is not None else set()
        self.neutral = neutral if neutral is not None else set()

    def __repr__(self):
        out_dic = {"name": self.name, "type": self.type, "mandatory": self.mandatory, "accessory": self.accessory,
                   "forbidden": self.forbidden, "neutral": self.neutral}
        return out_dic.update(
            {k: v for k, v in self.parameters.items()} if self.parameters is not None else {"parameters": None})

    @property
    def size(self) -> int:
        """Get the number of items added

        :return: number of model inside
        """
        return len(self.mandatory.union(self.accessory, self.forbidden, self.neutral))

    def check_parameters(self, param_dict: dict):
        """Check if all parameters are inside and with the good type

        :param param_dict: Dictionnary with all the parameters for the rule

        :raise KeyError: One or more mandatory parameters are missing
        :raise TypeError: One or more parameters are with a non-acceptable type
        :raise ValueError: One or more parameters are with a non-acceptable value
        :raise Exception: Manage an unexpected error
        """
        for param in suprules_params:
            try:
                param_val = param_dict[param]
            except KeyError:
                raise KeyError(f"No {param} in data rules")
            except Exception:
                raise Exception("Unexpected Error")
            else:
                if not isinstance(param_val, int) or param_val < 0:
                    raise ValueError(f"The {param} value is not positive int type")

    def check_suprule_key(self, rule_dict, key: str):
        """Check if key is inside and with the good type

        :param rule_dict: Dictionnary with all the parameters for the rule
        :param key: key which must be checked

        :raise KeyError: The given key is non-acceptable
        :raise TypeError: One or more parameters are with a non-acceptable type
        :raise ValueError: One or more parameters are with a non-acceptable value
        """
        super().check_rule_key(rule_dict, key)
        if key == 'parameters' and rule_dict[key] is not None:
            self.check_parameters(rule_dict[key])


class Model(SupRule):
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
        super().__init__(name=name, parameters=parameters, mandatory=mandatory, accessory=accessory,
                         forbidden=forbidden, neutral=neutral)
        self.bool_param = False
        self.canonical = canonical

    def __str__(self):
        out_dic = {"name": self.name, "type": self.type}
        out_dic.update(
            {k: v for k, v in self.parameters.items()} if self.parameters is not None else {"parameters": None})
        return ", ".join([f"{k}: {v}" for k, v in out_dic.items()]) + "\n\t" + \
               "\n\t".join(fu.__str__(ident="\n\t\t") for fu in self.func_units)

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

    def check_model_dict(self, model_dict):
        """Check if all keys and values are present and with good type before to add in model

        :param model_dict: Dictionnary with all model information

        :raise KeyError: One or more keys are missing
        :raise TypeError: One or more value are not with good type
        :raise ValueError: One or more value are not non-acceptable
        """
        try:
            check_key(model_dict, rule_keys[:-1])
        except KeyError:
            raise KeyError("One or more keys are missing in your file at model level")
        else:
            for key in model_dict.keys():
                if key == 'func_units':
                    if not isinstance(model_dict[key], list):
                        raise TypeError("func_unit value in json must be an array")
                    if len(model_dict[key]) < 1:
                        raise ValueError("Model need at least one functional unit")
                elif key == 'canonical':
                    if not isinstance(model_dict[key], list):
                        raise TypeError("canonical value in json must be an array")
                else:
                    self.check_rule_key(model_dict, key)

    def check_model(self):
        """Check model consistency

        :raise Exception: Model is not consistent
        """
        if len(self.mandatory) < self.parameters["min_mandatory"]:
            raise Exception(f"There is less mandatory than the minimum mandatory in {self.name}")
        if len(self.forbidden) < self.parameters["max_forbidden"]:
            raise Exception(f"There is less forbidden than the maximum forbidden accepted in {self.name}")
        if self.size < self.parameters["min_total"]:
            raise Exception(f"There is less functional units than the minimum total"
                            f" of functional unit needed in {self.name}")
        if len(self.neutral) == 1 and self.size == 1:
            logging.getLogger().warning(f"There is only one function unit in {self.name}. "
                                        f"It should be a mandatory type.")

    def read_model(self, data_model: dict):
        """Read model to parse in self attributes

        :param data_model: json data dictionary
        """
        self.check_model_dict(data_model)
        fu_set = set()
        for key, value in data_model.items():
            if key == 'func_units':
                for dict_fu in value:
                    f_unit = FuncUnit()
                    f_unit.model = self
                    f_unit.read_func_unit(dict_fu)
                    fu_set.add(f_unit)
            elif key == 'canonical':
                self.canonical = value
            else:
                self.read_rule(key, value)
        self.add_func_units(fu_set)
        self.check_model()

    def add_func_units(self, func_units: Set[FuncUnit]):
        """Add a set of functional unit in the system

        :param func_units: Functional unit which must be added
        """
        for fu in func_units:
            if fu.parameters is None:
                fu.parameters = self.parameters
            fu.check_func_unit()
            self.add_func_unit(fu)

    def add_func_unit(self, func_unit: FuncUnit):
        """Add a function unit in model

        :param func_unit: function unit
        """
        if func_unit.type == "mandatory":
            self.mandatory.add(func_unit)
        elif func_unit.type == "accessory":
            self.accessory.add(func_unit)
        elif func_unit.type == "forbidden":
            self.forbidden.add(func_unit)
        elif func_unit.type == "neutral":
            self.neutral.add(func_unit)


class FuncUnit(SupRule):
    """Represent functional unit definition rule

    :param name: Neme of the element
    :param parameters: Dictionary with the parameters name and value
    :param mtype: Type of the rule (mandatory, accessory, forbidden or neutral)
    :param mandatory: Set of mandatory sub element
    :param accessory: Set of accessory sub element
    :param forbidden: Set of forbidden sub element
    :param neutral: Set of neutral sub element
    """

    def __init__(self, name: str = "", parameters: dict = None, mtype: str = "", mandatory: set = None,
                 accessory: set = None, forbidden: set = None, neutral: set = None, model: Model = None):
        """Constructor Method
        """
        super().__init__(name, parameters, mtype, mandatory, accessory, forbidden, neutral)
        self.model = model

    def __str__(self, ident: str = "\n\t"):
        out_dic = {"name": self.name, "type": self.type,
                   "model": self.model.name if self.model is not None else None}
        out_dic.update(
            {k: v for k, v in self.parameters.items()} if self.parameters is not None else {"parameters": None})
        return ", ".join([f"{k}: {v}" for k, v in out_dic.items()]) + ident + \
               ident.join(str(fam) for fam in self.families)

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

    def check_fu_dict(self, fu_dict):
        """Check if all keys and values are present and with good type before to add in functional unit

        :param fu_dict: Dictionnary with all functional unit information

        :raise KeyError: One or more keys are missing
        :raise TypeError: One or more value are not with good type
        :raise ValueError: One or more value are not non-acceptable
        """
        try:
            check_key(fu_dict, rule_keys + ['families'])
        except KeyError:
            raise KeyError("One or more keys are missing in your file at functional unit level")
        else:
            for key in fu_dict.keys():
                if key == 'families':
                    if not isinstance(fu_dict[key], list):
                        raise TypeError("func_unit value in json must be an array")
                    if len(fu_dict[key]) < 1:
                        raise ValueError("Model need at least one functional unit")
                else:
                    self.check_rule_key(fu_dict, key)

    def check_func_unit(self):
        """Check functional unit consistency

        :raise Exception: Model is not consistent
        """
        if self.parameters["min_mandatory"] > len(self.mandatory) + len(list(
                self.duplicate_fam(filter_type="mandatory"))):
            raise Exception(f"There is less mandatory than the minimum mandatory in {self.name}")
        if self.parameters["max_forbidden"] > len(self.forbidden):
            raise Exception(f"There is less forbidden than the maximum forbidden accepted in {self.name}")
        if self.parameters["min_total"] > self.size + len(list(self.duplicate_fam())):
            raise Exception(f"There is less families than the minimum total"
                            f" of functional unit needed in {self.name}")
        if len(self.neutral) == 1 and self.size == 1:
            logging.getLogger().warning(f"There is only one function unit in {self.name}. "
                                        f"It should be a mandatory type.")

    def read_func_unit(self, data_fu: dict):
        """Read functional unit

        :param data_fu: data json file of all function units
        """
        self.check_fu_dict(data_fu)
        for key, value in data_fu.items():
            if key == 'families':
                for fam_dict in value:
                    family = Family(func_unit=self)
                    family.read_family(fam_dict)
                    self.add_family(family)
            else:
                self.read_rule(key, value)

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


class Family(Rule):
    """Represent family model definition rule

    :param name: Neme of the element
    :param parameters: Dictionary with the parameters name and value
    :param mtype: Type of the rule (mandatory, accessory, forbidden or neutral)
    :param relation:
    :param func_unit: Functional unit in which is the family
    """
    def __init__(self, name: str = "", parameters: dict = None, mtype: str = "", relation: str = "",
                 func_unit: FuncUnit = None):
        """Constructor Method
        """
        super().__init__(name, parameters, mtype)
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

    def check_fam_dict(self, fam_dict):
        """Check if all keys and values are present and with good type before to add in family

        :param fam_dict: Dictionnary with all family information

        :raise KeyError: One or more keys are missing
        :raise TypeError: One or more value are not with good type
        :raise ValueError: One or more value are not non-acceptable
        """
        try:
            check_key(fam_dict, rule_keys + ['relation'])
        except KeyError:
            raise KeyError("One or more keys are missing in your file at family level")
        else:
            for key in fam_dict.keys():
                if key == 'relation':
                    if fam_dict[key] is not None:
                        if not isinstance(fam_dict[key], str):
                            raise TypeError("relation value from family in json must be an str or null")
                        if fam_dict[key] not in ['homologs', 'analogs']:
                            raise ValueError("Family relation must be null homologs or analogs")
                elif key == 'duplicate':
                    if not isinstance(fam_dict[key], int):
                        raise TypeError("duplicate value from family in json must be an int or null")
                else:
                    self.check_rule_key(fam_dict, key)

    def read_family(self, data_fam: dict):
        """Read family

        :param data_fam: data json file with families
        """
        self.check_fam_dict(data_fam)
        for key, value in data_fam.items():
            if key == 'relation':
                self.relation = value
            elif key == 'duplicate':
                self.duplicate = value
            else:
                self.read_rule(key, value)
