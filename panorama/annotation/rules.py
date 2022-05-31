#!/usr/bin/env python3
# coding:utf-8
from __future__ import annotations


class Systems:

    def __init__(self, systems: dict = {}):
        self.systems = systems

    def print_systems(self):
        for system in self.systems.values():
            system.print_system()

    def get_sys(self, name: str):
        for rule in self.systems:
            if rule == name:
                return rule

    def add_sys(self, sys: System):
        if sys.check_param() is True:
            self.systems[sys.name] = sys


class Element:
    def __init__(self):
        self.name = ""
        self.parameters = dict()


class System(Element):
    def __init__(self):
        super().__init__()
        self.func_units = dict()
        self.bool_param = False
        self.families = set()
        self.mandatory = set()
        self.accessory = set()
        self.forbidden = set()
        self.neutral = set()

    def read_system(self, data):
        for key, value in data.items():
            if key == 'name':
                if isinstance(value, str):
                    self.name = value
                else:
                    raise Exception("The system name must to be str type.")
            elif key == 'func_units':
                for func_unit, dict_fu in data['func_units'].items():
                    f_unit = FuncUnit()
                    f_unit.name = func_unit
                    f_unit.system = self.name
                    f_unit.read_func_unit(dict_fu)
                    self.add_func_unit(f_unit)
                self.add_type_fu()
            elif key == 'parameters':
                self.parameters = value
            else:
                Exception(f"{key} doesn't exist in {self.name}")
        self.check_system(data)
        if len(self.neutral) == 1 and len(self.func_units) == 1:
            self.mandatory = self.neutral
            self.neutral = None

    def add_func_unit(self, func_unit: FuncUnit):
        self.func_units[func_unit.name] = func_unit

    def add_type_fu(self):
        for name, fu in self.func_units.items():
            if fu.type == "mandatory":
                self.mandatory.add(fu)
            elif fu.type == "accessory":
                self.accessory.add(fu)
            elif fu.type == "forbidden":
                self.forbidden.add(fu)
            elif fu.type == "neutral":
                self.neutral.add(fu)
        self.add_fam()

    def add_fam(self):
        for func_unit in self.func_units.values():
            self.families.update(func_unit.families)

    def check_system(self, dict_sys: dict):
        keys_fu = ['name', 'func_units', 'parameters']
        for key in keys_fu:
            try:
                sub_dict = dict_sys[key]
            except KeyError:
                raise KeyError(f"Requires {key} in {self.name}")
            else:
                if key == 'name':
                    if sub_dict is None:
                        raise Exception(f"The system must be to have a name")
                elif key == 'parameters':
                    if sub_dict is None:
                        for func_unit in self.func_units.values():
                            if func_unit.parameters is None:
                                raise Exception("It must to have the four parameters (min_mandatory, max_forbidden, "
                                                "min_total and max_separation) in system's parameters or in all "
                                                "func units's parameters")
                            else:
                                k_param = ['min_mandatory', 'max_forbidden', 'min_total', 'max_separation']
                                for parameter in k_param:
                                    try:
                                        _ = func_unit.parameters[parameter]
                                    except KeyError:
                                        raise KeyError("It must to have the four parameters (min_mandatory, "
                                                       "max_forbidden, min_total and max_separation) in system's "
                                                       "parameters or in all func units's parameters")
                    else:
                        keys_param = ['min_mandatory', 'max_forbidden', 'min_total', 'max_separation']
                        for parameter in keys_param:
                            try:
                                _ = sub_dict[parameter]
                            except:
                                for func_unit in self.func_units.values():
                                    if func_unit.parameters is None:
                                        raise Exception(
                                            "It must to have the four parameters (min_mandatory, max_forbidden, "
                                            "min_total and max_separation) in system's parameters or in all "
                                            "func units's parameters")
                                    else:
                                        k_param = ['min_mandatory', 'max_forbidden', 'min_total', 'max_separation']
                                        for parameter in k_param:
                                            try:
                                                _ = func_unit.parameters[parameter]
                                            except KeyError:
                                                raise KeyError(
                                                    "It must to have the four parameters (min_mandatory, max_forbidden,"
                                                    " min_total and max_separation) in system's parameters or in all "
                                                    "func units's parameters")

    def check_param(self):
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
                        if len(self.forbidden) <= value:
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
        if len(self.mandatory) == 0:
            self.mandatory = None
        if len(self.forbidden) == 0:
            self.forbidden = None
        if len(self.accessory) == 0:
            self.accessory = None


class FuncUnit(Element):
    def __init__(self):
        super().__init__()
        self.type = ""
        self.families = dict()
        self.system = ""
        self.bool_param = False
        self.mandatory = dict()
        self.accessory = dict()
        self.forbidden = dict()

    def read_func_unit(self, data_fu: dict):
        fam_fu = Family()
        fam_fu.name = self.name
        fam_fu.parameters = None
        if data_fu['families']:
            for family, dict_fam in data_fu['families'].items():
                fam = Family()
                fam.name = family
                fam.func_unit = self.name
                fam.system = self.system
                fam.read_family(dict_fam)
                self.add_fam(fam)
            self.add_type_fam()
        self.check_func_unit(data_fu)
        for key, value in data_fu.items():
            if key == "type":
                if value is None:
                    self.type = "neutral"
                else:
                    if value == 'mandatory' or value == 'accessory' or value == 'forbidden':
                        self.type = value
                    else:
                        raise Exception(f"{value} doesn't exist, in {self.name} of {self.system}")
                    self.type = value
                fam_fu.type = self.type
            elif key == "families":
                pass
            elif key == "parameters":
                self.parameters = value
            else:
                Exception(f"{key} doesn't exist in {self.name} of {self.system}")
        self.families[fam_fu.name] = fam_fu

    def add_type_fam(self):
        for name, fam in self.families.items():
            if fam.parameters is not None:
                if fam.type == "mandatory":
                    self.mandatory[name] = int(fam.parameters.get('max_separation'))
                elif fam.type == "accessory":
                    self.accessory[name] = int(fam.parameters.get('max_separation'))
                elif fam.type == "forbidden":
                    self.forbidden[name] = int(fam.parameters.get('max_separation'))
            else:
                if fam.type == "mandatory":
                    self.mandatory[name] = None
                elif fam.type == "accessory":
                    self.accessory[name] = None
                elif fam.type == "forbidden":
                    self.forbidden[name] = None

    def add_fam(self, fam: Family):
        self.families[fam.name] = fam

    def check_func_unit(self, dict_fu: dict):
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
                            elif not(isinstance(value, int)) or value < 0:
                                raise ValueError(f"The value {value} is not positive int type in "
                                                 f"{self.name} of {self.system}")

    def check_param(self):
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
                        if len(self.forbidden) <= value:
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
        self.pre_print()
        print(f"\t\t{self.name} :")
        print(f"\t\t\tparameters : {self.parameters}")
        print(f"\t\t\tmandatory : {self.pre_print_str(self.mandatory)}")
        print(f"\t\t\taccessory : {self.pre_print_str(self.accessory)}")
        print(f"\t\t\tforbidden : {self.pre_print_str(self.forbidden)}")

    def pre_print_str(self, dict_type: dict):
        if dict_type is not None:
            string = ""
            for name, parameters in dict_type.items():
                if parameters is not None:
                    string += f"{name}(max_separation : {parameters})  "
                else:
                    string += f"{name}  "
            return string

    def pre_print(self):
        if len(self.mandatory) == 0:
            self.mandatory = None
        if len(self.forbidden) == 0:
            self.forbidden = None
        if len(self.accessory) == 0:
            self.accessory = None

    def get_sys(self, system_name: str, systems: Systems):
        return systems.systems[system_name]


class Family(Element):
    def __init__(self):
        super().__init__()
        self.type = ""
        self.relation = ""
        self.func_unit = ""
        self.system = ""

    def read_family(self, family: dict):
        self.check_family(family)
        for attrib, value in family.items():
            if attrib == 'type':
                self.type = value
            elif attrib == 'relation':
                self.relation = value
            elif attrib == 'parameters':
                self.parameters = value

    def check_family(self, dict_fam: dict):
        keys_fam = ['type', 'relation', 'parameters']
        for key in keys_fam:
            try:
                _ = dict_fam[key]
            except KeyError:
                raise KeyError(f"Requires {key} in {self.func_unit} {self.name} of {self.system}")
            else:
                if key == 'type':
                    if not (dict_fam['type'] == 'mandatory' or dict_fam['type'] == 'accessory'
                            or dict_fam['type'] == 'forbidden'):
                        raise KeyError(f"The type value {dict_fam['type']} is incorrect in {self.func_unit} {self.name}"
                                       f" of {self.system}, it must be mandatory, accessory or forbidden")
                elif key == 'relation':
                    if not (dict_fam['relation'] is None or dict_fam['relation'] == 'homologs'
                            or dict_fam['relation'] == 'analogs'):
                        raise KeyError(f"The relation value {dict_fam['relation']} is incorrect in {self.func_unit} "
                                       f"{self.name} of {self.system}, it must be None, homologs or analogs")
                elif key == 'parameters':
                    if dict_fam['parameters'] is not None:
                        for parameter, value in dict_fam['parameters'].items():
                            if not(parameter == 'max_separation'):
                                raise KeyError(f"The parameter {parameter} is incorrect. In family parameters, there is"
                                               f" only max_separation in {self.func_unit} {self.name} of {self.system}")
                            elif not(isinstance(value, int)) or value < 0:
                                raise ValueError(f"The value {value} is not positive int type in {self.func_unit} "
                                                 f"{self.name} of {self.system}")
                else:
                    raise Exception(f"{key} doesn't exit in {self.func_unit} {self.name} of {self.system}, only type, "
                                    "relation and parameters exist")

    def check_param(self):
        return True

    def get_func_unit(self, func_unit: str, systems: Systems):
        return self.get_sys(self.system, systems).func_units[func_unit]

    def get_sys(self, system_name: str, systems: Systems):
        return systems.systems[system_name]
