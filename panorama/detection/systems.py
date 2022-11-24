#!/usr/bin/env python3
# coding:utf-8

# default libraries
from __future__ import annotations
import logging
from typing import Set


suprules_params = ['min_mandatory', 'max_forbidden', 'min_total']
keys_param = suprules_params.append('max_separation')
rule_keys = ['name', 'parameters', 'type']
accept_type = ['mandatory', 'accessory', 'forbidden', 'neutral']


def check_key(dict, need_key):
    if not all([key in dict for key in need_key]):
        raise KeyError(f"All the following key are necessary : {need_key}")


class Systems:
    """
    :param systems: The internal identifier to give to the gene family
    """

    def __init__(self, systems: Set[System] = None):
        """Constructor Method
        """
        self._system_getter = systems if systems is not None else {}

    # def __str__(self):
    #     """
    #     Print all systems predicted
    #
    #     """
    #     for system in self.__iter__():
    #         return system

    def __iter__(self):
        for sys in self._system_getter.values():
            yield sys

    @property
    def systems(self):
        return list(self)

    @property
    def size(self):
        return len(self.systems)

    @property
    def func_units(self):
        func_unit_dict = {}
        for sys in self.__iter__():
            for fu in sys.func_units:
                if fu.name not in func_unit_dict:
                    func_unit_dict[fu.name] = [sys.name]
                else:
                    func_unit_dict[fu.name].append(sys.name)
        return func_unit_dict

    @property
    def families(self):
        families_dict = {}
        for sys in self.__iter__():
            for family in sys.families:
                if family.name not in families_dict:
                    families_dict[family.name] = [sys.name]
                else:
                    families_dict[family.name].append(sys.name)
        return families_dict

    def get_sys(self, name: str):
        """
        Get name system

        :param name: name to find
        """
        try:
            sys = self._system_getter[name]
        except KeyError:
            raise KeyError("System not present in set of systems")
        else:
            return sys

    def add_sys(self, sys: System):
        """
        Add system

        :param sys: system
        """
        try:
            self.get_sys(sys.name)
        except KeyError:
            self._system_getter[sys.name] = sys
        else:
            raise Exception(f"System {sys.name} already in set of systems")


class Rule:
    """Super class for System and FuncUnit and Family

    :param name: Neme of the element
    :param parameters: Dictionary with the parameters name and value
    :param type: Type of the rule (mandatory, accessory, forbidden or neutral)
    """

    def __init__(self, name: str = "", parameters: dict = None, type: str = ""):
        """Constructor Method
        """
        self.name = name
        self.parameters = parameters if parameters is not None else dict()
        self.type = type

    def check_param(self, param_dict):
        try:
            max_sep = param_dict['max_separation']
        except KeyError:
            raise KeyError("No max_separation in data rules")
        except Exception:
            raise Exception("Unexpected Error")
        else:
            if not isinstance(max_sep, int) or max_sep < 0:
                raise ValueError("The max_separation value is not positive int type")

    def check_rule_key(self, rule_dict, key: str):
        if key not in rule_keys:
            raise KeyError(f"Unreadable key : {key}. Accepeted key are {rule_keys}")
        else:
            if key == 'name':
                if not isinstance(rule_dict[key], str):
                    raise Exception(f"The name value must be str type")
            elif key == 'type':
                if rule_dict[key] not in accept_type:
                    raise KeyError(f"The type value {rule_dict[key]} is incorrect in {type(self)}")
            elif key == 'parameters' and rule_dict[key] is not None:
                self.check_param(rule_dict[key])

    def read_rule(self, key: str, value):
        if key == 'name':
            self.name = value
        elif key == 'type':
            self.type = value
        elif key == 'parameters':
            self.parameters = value


class SupRule(Rule):
    """Super class for System and FuncUnit


    :param name: Neme of the element
    :param parameters: Dictionary with the parameters name and value
    :param type: Type of the rule (mandatory, accessory, forbidden or neutral)
    :param mandatory: Set of mandatory sub element
    :param accessory: Set of accessory sub element
    :param forbidden: Set of forbidden sub element
    :param neutral: Set of neutral sub element
    """

    def __init__(self, name: str = "", parameters: dict = None, type: str = "", mandatory: set = None,
                 accessory: set = None, forbidden: set = None, neutral: set = None):
        """Constructor Method
        """
        super().__init__(name, parameters, type)
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
    def size(self):
        return len(self.mandatory.union(self.accessory, self.forbidden, self.neutral))

    def mandatory_name(self):
        return [elem.name for elem in self.mandatory]

    def accessory_name(self):
        return [elem.name for elem in self.accessory]

    def forbidden_name(self):
        return [elem.name for elem in self.forbidden]

    def neutral_name(self):
        return [elem.name for elem in self.neutral]

    def check_param(self, param_dict: dict):
        """
        Verify with boolean the condition parameters

        :param param_dict: Message to print if problem with parameters
        """
        super(SupRule, self).check_param(param_dict)
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

    def check_suprules(self):
        if len(self.mandatory) < self.parameters["min_mandatory"]:
            raise Exception(f"There is less mandatory than the minimum mandatory in {self.name}")
        if len(self.forbidden) < self.parameters["max_forbidden"]:
            raise Exception(f"There is less forbidden than the maximum forbidden accepted in {self.name}")


class System(SupRule):
    """Represent System rules which describe biological system

     :param name: Neme of the element
     :param parameters: Dictionary with the parameters name and value
     :param mandatory: Set of mandatory sub element
     :param accessory: Set of accessory sub element
     :param forbidden: Set of forbidden sub element
     :param neutral: Set of neutral sub element
     """

    def __init__(self, name: str = "", parameters: dict = None, mandatory: set = None, accessory: set = None,
                 forbidden: set = None, neutral: set = None):
        """Constructor Method
        """
        super().__init__(name, parameters, mandatory, accessory, forbidden, neutral)
        self.bool_param = False

    def __str__(self):
        out_dic = {"name": self.name, "type": self.type}
        out_dic.update(
            {k: v for k, v in self.parameters.items()} if self.parameters is not None else {"parameters": None})
        return ", ".join([f"{k}: {v}" for k, v in out_dic.items()]) + "\n\t" + \
               "\n\t".join(fu.__str__(ident="\n\t\t") for fu in self.func_units)

    @property
    def func_units(self):
        for func_unit in self.mandatory.union(self.accessory, self.forbidden, self.neutral):
            yield func_unit

    @property
    def families(self):
        for func_unit in self.func_units:
            yield from func_unit.families

    def func_unit_name(self):
        return [fu.name for fu in self.func_units]

    def families_name(self):
        return [fam.name for fam in self.families]

    def check_system_dict(self, system_dict):
        try:
            check_key(system_dict, rule_keys[:-1])
        except KeyError:
            raise KeyError("One or more keys are missing in your file at system level")
        else:
            for key in system_dict.keys():
                if key == 'func_units':
                    if not isinstance(system_dict[key], list):
                        raise TypeError("func_unit value in json must be an array")
                    if len(system_dict[key]) < 1:
                        raise ValueError("System need at least one functional unit")
                else:
                    self.check_rule_key(system_dict, key)

    def check_system(self):
        super(System, self).check_suprules()
        if self.size < self.parameters["min_total"]:
            raise Exception(f"There is less functional units than the minimum total"
                            f" of functional unit needed in {self.name}")
        if len(self.neutral) == 1 and self.size == 1:
            logging.getLogger().warning(f"There is only one function unit in {self.name}. "
                                        f"It should be a mandatory type.")

    def read_system(self, data_sys: dict):
        """
        Read system to parse in self attributes

        :param data_sys: json data dictionary
        """
        self.check_system_dict(data_sys)
        fu_set = set()
        for key, value in data_sys.items():
            if key == 'func_units':
                for dict_fu in value:
                    f_unit = FuncUnit()
                    f_unit.system = self
                    f_unit.read_func_unit(dict_fu)
                    fu_set.add(f_unit)
            else:
                self.read_rule(key, value)
        self.add_func_units(fu_set)
        self.check_system()

    def add_func_units(self, func_units: Set[FuncUnit]):
        for fu in func_units:
            if fu.parameters is None:
                fu.parameters = self.parameters
            fu.check_func_unit()
            self.add_func_unit(fu)

    def add_func_unit(self, func_unit: FuncUnit):
        """
        Add function unit in function units dictionary of system

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
    :param type: Type of the rule (mandatory, accessory, forbidden or neutral)
    :param mandatory: Set of mandatory sub element
    :param accessory: Set of accessory sub element
    :param forbidden: Set of forbidden sub element
    :param neutral: Set of neutral sub element
    """

    def __init__(self, name: str = "", parameters: dict = None, type: str = "", mandatory: set = None,
                 accessory: set = None, forbidden: set = None, neutral: set = None, system: System = None):
        """Constructor Method
        """
        super().__init__(name, parameters, type, mandatory, accessory, forbidden, neutral)
        self.system = system

    def __str__(self, ident: str = "\n\t"):
        out_dic = {"name": self.name, "type": self.type,
                   "system": self.system.name if self.system is not None else None}
        out_dic.update(
            {k: v for k, v in self.parameters.items()} if self.parameters is not None else {"parameters": None})
        return ", ".join([f"{k}: {v}" for k, v in out_dic.items()]) + ident + \
               ident.join(str(fam) for fam in self.families)

    @property
    def families(self):
        for fam in self.mandatory.union(self.accessory, self.forbidden, self.neutral):
            yield fam

    def families_name(self):
        return [fam.name for fam in self.families]

    def check_fu_dict(self, fu_dict):
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
                        raise ValueError("System need at least one functional unit")
                else:
                    self.check_rule_key(fu_dict, key)

    def check_func_unit(self):
        super(FuncUnit, self).check_suprules()
        if self.size < self.parameters["min_total"]:
            raise Exception(f"There is less families than the minimum total"
                            f" of functional unit needed in {self.name}")
        if len(self.neutral) == 1 and self.size == 1:
            logging.getLogger().warning(f"There is only one function unit in {self.name}. "
                                        f"It should be a mandatory type.")

    def read_func_unit(self, data_fu: dict):
        """
        Read function unit

        :param data_fu: data json file of all function units

        """
        self.check_fu_dict(data_fu)
        for key, value in data_fu.items():
            if key == 'families':
                for fam_dict in value:
                    family = Family(func_unit=self)
                    family.read_family(fam_dict)
                    self.add_fam(family)
            else:
                self.read_rule(key, value)

    def add_fam(self, fam: Family):
        """
        Add family in families dictionary of function unit

        :param fam: a family
        """
        if fam.type == 'mandatory':
            self.mandatory.add(fam)
        elif fam.type == "accessory":
            self.accessory.add(fam)
        elif fam.type == "forbidden":
            self.forbidden.add(fam)
        else:
            self.neutral.add(fam)


class Family(Rule):
    def __init__(self, name: str = "", parameters: dict = None, type: str = "", relation: str = "",
                 func_unit: FuncUnit = None):
        """Constructor Method
        """
        super().__init__(name, parameters, type)
        self.relation = relation
        self.func_unit = func_unit

    def __repr__(self):
        return f"Family name : {self.name}, type : {self.type}, " \
               f"Functional Unit : {self.func_unit.name if self.func_unit is not None else 'Unknow'}, " \
               f"System : {self.system.name if self.system is not None else 'Unknow'}"

    def __str__(self):
        out_dic = {"name": self.name, "type": self.type, "relation": self.relation,
                   "functional unit": self.func_unit.name if self.func_unit is not None else None,
                   "system": self.system.name if self.system is not None else None}
        out_dic.update(
            {k: v for k, v in self.parameters.items()} if self.parameters is not None else {"parameters": None})
        return ", ".join([f"{k}: {v}" for k, v in out_dic.items()])

    @property
    def system(self):
        return self.func_unit.system

    def check_fam_dict(self, fam_dict):
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
                else:
                    self.check_rule_key(fam_dict, key)

    def read_family(self, data_fam: dict):
        """
        Read family

        :param data_fam: data json file with families
        """
        self.check_fam_dict(data_fam)
        for key, value in data_fam.items():
            if key == 'relation':
                self.relation = value
            else:
                self.read_rule(key, value)


if __name__ == "__main__":
    import json

    with open("./rule_one.json", 'r') as json_file:
        data = json.load(json_file)
        system = System()
        system.read_system(data)
        print(system.func_units)