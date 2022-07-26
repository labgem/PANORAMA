#!/usr/bin/env python3
# coding:utf-8

# default libraries
from __future__ import annotations
import logging
from typing import Dict

suprules_params = ['min_mandatory', 'max_forbidden', 'min_total']
keys_param = suprules_params.append('max_separation')
rule_keys = ['type', 'parameters']


class Systems:
    """
    :param systems: The internal identifier to give to the gene family
    """
    def __init__(self, systems: Dict[System] = None):
        """Constructor Method
        """
        self.systems = systems

    def __str__(self):

        """
        Print all systems predicted

        """
        for system in self.systems.values():
            system.print_system()

    def in_sys(self, name: str):
        if name in self.systems:
            return True
        else:
            raise KeyError(f"The {name} system is not in the list of systems")

    def get_sys(self, name: str):
        """
        Get name system

        :param name: name to find
        """
        for system in self.systems:
            if system == name:
                return system

    def get_sys2fam(self, name_fam: str):
        """
        Get name system of family

        :param name_fam: name of family to find its system
        :return: system
        """
        for obj_fam, system in self.dict_families.items():
            if obj_fam.name == name_fam:
                return system

    def add_sys(self, sys: System):
        """
        Add system

        :param sys: system
        """
        if sys.check_param() is True:
            self.systems[sys.name] = sys
            for family in sys.families.values():
                self.dict_families[family] = sys.name


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

    def check_param(self, message):
        """
        Verify with boolean the condition parameters

        :param message: Message to print if problem with parameters
        """
        if self.parameters is not None:
            dict_bool = dict()

    def check_dict(self, rule_dict):
        for key in rule_keys:
            try:
                _ = rule_dict[key]
            except KeyError:
                raise KeyError(f"Requires {key} in {self.name}")
            else:
                if key == 'type':
                    if not rule_dict[key] in ['mandatory', 'accessory', 'forbidden']:
                        raise KeyError(f"The type value {rule_dict[key]} is incorrect in {self.name}")
                elif key == 'parameters' and rule_dict[key] is not None:
                    for parameter, value in rule_dict['parameters'].items():
                        if parameter != 'max_separation':
                            raise KeyError(f"The parameter {parameter} is incorrect in {self.name}")
                        elif not (isinstance(value, int)) or value < 0:
                            raise ValueError(f"The value {value} is not positive int type in {self.name}")

    def read_element(self, data: dict):
        for key, value in data.items():
            if isinstance(key, str):
                self.name = key
            else:
                raise Exception(f"The name value must be str type")
            if isinstance(value, dict):
                return value
            else:
                raise Exception(f"The element are not readable")


class SupRule(Rule):
    """Super class for System and FuncUnit

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

    def check_param(self, message):
        """
        Verify with boolean the condition parameters

        :param message: Message to print if problem with parameters
        """
        if self.parameters is not None:
            dict_bool = dict()
            for parameter, value in self.parameters.items():
                if parameter == "min_mandatory":
                    if len(self.mandatory) >= value:
                        dict_bool[parameter] = True
                    else:
                        dict_bool[parameter] = False
                elif parameter == "max_forbidden":
                    if len(self.forbidden) >= value:
                        dict_bool[parameter] = True
                    else:
                        dict_bool[parameter] = False
                elif parameter == "min_total":
                    if len(self.families) >= value:
                        dict_bool[parameter] = True
                    else:
                        dict_bool[parameter] = False
                elif parameter == "max_separation":
                    dict_bool[parameter] = True
            if all(x for x in dict_bool.values()):
                return True
            else:
                for param, value in dict_bool.items():
                    if value is False:
                        logging.getLogger().warning(message + f"{param}")
                return False
        else:
            return True

    def check_dict(self, data_dict: dict):
        super().check_dict(data_dict)
        for parameter, value in data_dict['parameters'].items():
            if parameter not in suprules_params:
                raise KeyError(f"The parameter {parameter} is incorrect in {self.name}")
            elif not (isinstance(value, int)) or value < 0:
                raise ValueError(f"The value {value} is not positive int type in {self.name}")

    def read_element(self, data: dict):
        sub_dict = super(SupRule, self).read_element(data)
        self.check_dict(sub_dict)
        self.type = sub_dict['type'] if sub_dict is not None else 'neutral'
        self.parameters = sub_dict['parameters']
        return sub_dict
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

    def __repr__(self):
        return [self.name, self.mandatory, self.accessory, self.forbidden, self.neutral, self.parameters]

    def __str__(self):
        "\t".join([self.name, ",".join(self.mandatory), ",".join(self.accessory), ",".join(self.forbidden),
                   ",".join(self.neutral), ",".join(self.parameters.values())])

    def func_units(self):
        return self.mandatory.union(self.accessory, self.forbidden, self.neutral)

    def families(self):
        families = dict()
        for func_unit in self.func_units():
            families.update(func_unit.families)
        return families

    def read_system(self, data: dict):
        """
        Read system to parse in self attributes

        :param data: json data dictionary
        """
        sys_dict = super().read_element(data)
        if not any([True for k in sys_dict.keys() if k in ['func_units', 'parameters']]):
            raise KeyError(f"Requires func_units and parameters keys in {self.name}")
        for key, value in sys_dict.items():
            if key == 'func_units':
                for dict_fu in value.values():
                    f_unit = FuncUnit()
                    f_unit.system = self.name
                    f_unit.read_func_unit(dict_fu)
                    self.add_func_unit(f_unit)
            elif key == 'parameters':
                if not all([True for k in value.keys() if k in keys_param]):
                    raise Exception(f"Not all the following key are in {keys_param} in {self.name}")
                self.parameters = value
            else:
                raise Exception(f"{key} in {self.name} are not readable in systems")
        self.check_system()
        if len(self.neutral) == 1 and len(self.func_units) == 1:
            logging.getLogger().warning(f"There is only one function unit in {self.name}. "
                                        f"It should be a mandatory type.")

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

    def check_system(self):
        """
        Verify system coherence

        :param dict_sys: data json file
        """
        pass
        # keys_fu = ['name', 'func_units', 'parameters']
        # for key in keys_fu:
        #     try:
        #         sub_dict = dict_sys[key]
        #     except KeyError:
        #         raise KeyError(f"Requires {key} in {self.name}")
        #     else:
        #         if key == 'name':
        #             if sub_dict is None:
        #                 raise Exception(f"The system must be to have a name")
        #         elif key == 'parameters':
        #             if sub_dict is None:
        #                 for func_unit in self.func_units.values():
        #                     if func_unit.parameters is None:
        #                         raise Exception("It must to have the four parameters (min_mandatory, max_forbidden, "
        #                                         "min_total and max_separation) in system's parameters or in all "
        #                                         "func units's parameters")
        #                     else:
        #                         k_param = ['min_mandatory', 'max_forbidden', 'min_total', 'max_separation']
        #                         for parameter in k_param:
        #                             try:
        #                                 _ = func_unit.parameters[parameter]
        #                             except KeyError:
        #                                 raise KeyError("It must to have the four parameters (min_mandatory, "
        #                                                "max_forbidden, min_total and max_separation) in system's "
        #                                                "parameters or in all func units's parameters")
        #             else:
        #                 keys_param = ['min_mandatory', 'max_forbidden', 'min_total', 'max_separation']
        #                 for parameter in keys_param:
        #                     try:
        #                         _ = sub_dict[parameter]
        #                     except:
        #                         for func_unit in self.func_units.values():
        #                             if func_unit.parameters is None:
        #                                 raise Exception(
        #                                     "It must to have the four parameters (min_mandatory, max_forbidden, "
        #                                     "min_total and max_separation) in system's parameters or in all "
        #                                     "func units's parameters")
        #                             else:
        #                                 k_param = ['min_mandatory', 'max_forbidden', 'min_total', 'max_separation']
        #                                 for parameter in k_param:
        #                                     try:
        #                                         _ = func_unit.parameters[parameter]
        #                                     except KeyError:
        #                                         raise KeyError(
        #                                             "It must to have the four parameters (min_mandatory, max_forbidden,"
        #                                             " min_total and max_separation) in system's parameters or in all "
        #                                             "func units's parameters")

    def check_param(self):
        """
        Verify with boolean the condition parameters

        """
        bool_param_fu = False
        for f_unit in self.func_units.values():
            bool_param_fu = f_unit.check_param()
        if bool_param_fu is not False:
            if self.parameters is not None:
                dict_bool = dict()
                for parameter, value in self.parameters.items():
                    if parameter == "min_mandatory":
                        if len(self.mandatory) >= value:
                            dict_bool[parameter] = True
                        else:
                            dict_bool[parameter] = False
                    elif parameter == "max_forbidden":
                        if len(self.forbidden) >= value:
                            dict_bool[parameter] = True
                        else:
                            dict_bool[parameter] = False
                    elif parameter == "min_total":
                        if len(self.families) >= value:
                            dict_bool[parameter] = True
                        else:
                            dict_bool[parameter] = False
                    elif parameter == "max_separation":
                        dict_bool[parameter] = True
                if all(x for x in dict_bool.values()):
                    return True
                else:
                    for param, value in dict_bool.items():
                        if value is False:
                            print(f"The system {self.name} doesn't respect the parameter {param}.")
                    return False
            else:
                return True
        else:
            return False

    def print_system(self):
        """
        Print system

        """
        self.pre_print()
        print(f"{self.name} :")
        print(f"\tparameters: {self.parameters}")
        if self.mandatory:
            print(f"\tmandatory :")
            for f_unit in self.mandatory:
                f_unit.print_func_unit()
        if self.accessory:
            print(f"\taccessory :")
            for f_unit in self.accessory:
                f_unit.print_func_unit()
        if self.forbidden:
            print(f"\tforbidden :")
            for f_unit in self.forbidden:
                f_unit.print_func_unit()
        if self.neutral:
            print(f"\tneutral :")
            for f_unit in self.neutral:
                f_unit.print_func_unit()

    def pre_print(self):
        """
        If type is empty, it's none

        """
        if len(self.mandatory) == 0:
            self.mandatory = None
        if len(self.forbidden) == 0:
            self.forbidden = None
        if len(self.accessory) == 0:
            self.accessory = None


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
        self.bool_param = False

    def __repr__(self):
        return f"Functional unit name : {self.name}, type : {self.type}"

    @property
    def families(self):
        for fam in self.mandatory.union(self.accessory, self.forbidden, self.neutral):
            yield fam

    def check_dict(self, fu_dict):
        try:
            super(FuncUnit, self).check_dict(fu_dict)
        except Exception:
            raise Exception(f"Functional unit dictionnary unreadable in "
                            f"{self.name} from {self.system.name if self.system is not None else 'Unknow system'}")
        else:
            try:
                _ = fu_dict['families']
            except KeyError:
                raise KeyError(f"Functional unit {self.name} from "
                               f"{self.system.name if self.system is not None else 'Unknow system'}"
                               f"can't be without families.")

    def read_func_unit(self, data_fu: dict):
        """
        Read function unit

        :param data_fu: data json file of all function units

        """
        fu_dict = super().read_element(data_fu)
        for family_dict in fu_dict['families'].values():
            curent_family = Family()
            curent_family.func_unit = self.name
            curent_family.read_family(family_dict)
            self.add_fam(curent_family)

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

    def check_func_unit(self, dict_fu: dict):
        """
        Verify function units of system

        :param dict_fu: data function unit of json file
        """
        keys_fu = ['type', 'families', 'parameters']
        for key in keys_fu:
            try:
                sub_dict = dict_fu[key]
            except KeyError:
                raise KeyError(f"Requires {key} in {self.name} of {self.system}")
            else:
                if key == 'type':
                    if sub_dict is not None:
                        if not (sub_dict == 'mandatory' or sub_dict == 'accessory'
                                or sub_dict == 'forbidden'):
                            raise KeyError(f"The type value is incorrect in {self.name} of {self.system}, "
                                           f"it must be mandatory, accessory or forbidden")
                elif key == 'parameters':
                    if sub_dict is not None:
                        for parameter, value in sub_dict.items():
                            if not (parameter == 'min_mandatory' or parameter == 'min_total' or parameter ==
                                    'max_forbidden' or parameter == 'max_separation'):
                                raise Exception(f"{parameter} doesn't exist in {self.name} of {self.system}")
                            elif not (isinstance(value, int)) or value < 0:
                                raise ValueError(f"The value {value} is not positive int type in "
                                                 f"{self.name} of {self.system}")

    def check_param(self):
        """
         Verify with boolean the condition parameters

        """
        bool_param_fam = False
        for fam in self.families.values():
            bool_param_fam = fam.check_param()
        if bool_param_fam is not False:
            if self.parameters is not None:
                dict_bool = dict()
                for parameter, value in self.parameters.items():
                    if parameter == "min_mandatory":
                        if len(self.mandatory) >= value:
                            dict_bool[parameter] = True
                        else:
                            dict_bool[parameter] = False
                    elif parameter == "max_forbidden":
                        if len(self.forbidden) >= value:
                            dict_bool[parameter] = True
                        else:
                            dict_bool[parameter] = False
                    elif parameter == "min_total":
                        if len(self.families) >= value:
                            dict_bool[parameter] = True
                        else:
                            dict_bool[parameter] = False
                    elif parameter == "max_separation":
                        dict_bool[parameter] = True
                if all(x for x in dict_bool.values()):
                    return True
                else:
                    for param, value in dict_bool.items():
                        if value is False:
                            print(f"The function unit {self.name} of system {self.system} doesn't respect the parameter"
                                  f" {param}.")
                    return False
            else:
                return None
        else:
            return False

    def print_func_unit(self):
        """
        Print function units

        """
        self.pre_print()
        print(f"\t\t{self.name} :")
        print(f"\t\t\tparameters : {self.parameters}")
        print(f"\t\t\tmandatory : {self.pre_print_str(self.mandatory)}")
        print(f"\t\t\taccessory : {self.pre_print_str(self.accessory)}")
        print(f"\t\t\tforbidden : {self.pre_print_str(self.forbidden)}")

    def pre_print_str(self, dict_type: dict):
        """
        If type is empty, it's none.
        Print families

        :param dict_type: type dictionary of families
        """
        if dict_type is not None:
            string = ""
            for name, parameters in dict_type.items():
                if parameters is not None:
                    string += f"{name}(max_separation : {parameters})  "
                else:
                    string += f"{name}  "
            return string

    def pre_print(self):
        """
        If type is empty, it's none

        """
        if len(self.mandatory) == 0:
            self.mandatory = None
        if len(self.forbidden) == 0:
            self.forbidden = None
        if len(self.accessory) == 0:
            self.accessory = None

    def get_sys(self, system_name: str, systems: Systems):
        """
        Get system of function unit

        :param system_name: name system to find
        :param systems: class Systems with all systems
        :return name system
        """
        return systems.systems[system_name]


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

    def check_dict(self, family_dict: dict):
        try:
            super(Family, self).check_dict(family_dict)
        except Exception:
            raise Exception(f"Family dictionnary unreadable in "
                            f"{self.func_unit.name if self.func_unit is not None else 'Unknow functional unit'} from "
                            f"{self.system.name if self.system is not None else 'Unknow system'}")
        else:
            try:
                _ = family_dict['relation']
            except KeyError:
                raise KeyError(f"The relation attribute does not exist in {self.name} from "
                               f"{self.func_unit.name if self.func_unit is not None else 'Unknow functional unit'} "
                               f"and {self.system.name if self.system is not None else 'Unknow system'}")
            else:
                if family_dict['relation'] not in [None, 'homologs', 'analogs']:
                    raise ValueError(f"The relation value is incorrect in {self.name} from "
                                     f"{self.func_unit.name if self.func_unit is not None else 'Unknow func unit'} "
                                     f"and {self.system.name if self.system is not None else 'Unknow system'}. "
                                     f"It must be None, homologs or analogs")

    def read_family(self, data_dict: dict):
        """
        Read family

        :param family_dict: data json file with families
        """
        family_dict = self.read_element(data_dict)
        self.check_dict(family_dict)
        for attrib, value in family_dict.items():
            if attrib == 'type':
                self.type = value if value is not None else 'neutral'
            elif attrib == 'relation':
                self.relation = value
            elif attrib == 'parameters':
                self.parameters = value

    @property
    def system(self):
        return self.func_unit.system


if __name__ == "__main__":
    """Test Family class"""
    fam_dict = {"type": "mandatory", "relation": 'a', "parameters": {"max_separation": 4}}
    family = Family(name="test_fam")
    family.read_family(fam_dict)
    print(family)
