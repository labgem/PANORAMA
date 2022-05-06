#!/usr/bin/env python3
# coding:utf-8


class Rule:
    """
        A class used to represent a gene context
    """
    def __init__(self, name_rule=str):
        self.mandatory = {}
        self.accessory = {}
        self.forbidden = {}
        self.free = {}
        self.name = name_rule
        self.legal = {}

    def print_attr(self, attr: dict = None):
        """
        If attribute given, print the attribute rule. Else, print all rules attributes.

        :param attr: attributes dictionary
        """
        if attr is None:
            attr = {'mandatory': True, 'accessory': True, 'forbidden': True, 'free': True}
        for key, value in attr.items():
            if value:
                if key == 'mandatory':
                    print(f"{key} : {self.mandatory}")
                elif key == 'accessory':
                    print(f"{key} : {self.accessory}")
                elif key == 'forbidden':
                    print(f"{key} : {self.forbidden}")
                elif key == 'free':
                    print(f"{key} : {self.free}")

    def fill_attr(self, data: dict):
        """
        Transmits json file rule values in the arguments attributes, calculate if parameters have been respected and
        print parameters not respected

        :param data: dictionary of content json file rule
        """

        self.name = data.get('name')    # transmit values
        for key, value in data['families'].items():
            if key == 'mandatory':
                self.mandatory = set(value)
            elif key == 'accessory':
                self.accessory = set(value)
            elif key == 'forbidden':
                self.forbidden = set(value)
            else:
                self.free = set(value)
        for key, value in data['parameters'].items():  # calculate parameters
            if key == 'min_mandatory':
                if len(self.mandatory) >= value:
                    self.legal['min_mandatory'] = True
                else:
                    self.legal['min_mandatory'] = False
            elif key == 'min_total':
                if sum(len(x) for x in [self.mandatory, self.accessory,
                                        self.forbidden, self.free]) >= value:
                    self.legal['min_total'] = True
                else:
                    self.legal['min_total'] = False
            elif key == 'max_forbidden':
                if len(self.forbidden) <= value:
                    self.legal['max_forbidden'] = True
                else:
                    self.legal['max_forbidden'] = False
            elif key == 'max_free':
                if len(self.free) <= value:
                    self.legal['max_free'] = True
                else:
                    self.legal['max_free'] = False
            else:
                Exception(f"The parameter {key} doesn't exist, use this parameters : min_mandatory, min_total,"
                          f" max_forbidden, max_free")
        if any(not x for x in self.legal.values()):     # print parameters not respected
            print(f"\n{self.name}")
            for parameter, boolean in self.legal.items():
                if boolean is False:
                    print(f"{parameter} isn't respected")


class Rules:

    def __init__(self):
        self.rules = {}

    def add_rule(self, rule: Rule):
        """
        Add one rule in rules dictionary

        :param rule: the json file rule in object Rule
        """
        if all(x for x in rule.legal.values()):
            self.rules[rule.name] = rule

    def show_rules(self):
        """
        Print all the rules and the number of rules in dictionary
        """
        for r_name, rule in self.rules.items():
            print(f"\n{r_name}")
            rule.print_attr()
        print(f"\nnumber of rules : {len(self.rules)}")
