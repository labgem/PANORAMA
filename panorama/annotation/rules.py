#!/usr/bin/env python3
# coding:utf-8
import json


# Define association functions of rules
class Rule:
    # def __init__(self):
    def __init__(self, mandatory: set = None, accessory: set = None, forbidden: set = None, other: set = None):
        self.mandatory = mandatory
        self.accessory = accessory
        self.forbidden = forbidden
        self.other = other

    # Print rule attribute
    def print_rule(self, attr: str = None):
        if attr is None:
            print(f"mandatory : {str(self.mandatory)}\n"
                  f"accessory : {str(self.accessory)}\n"
                  f"forbidden : {str(self.forbidden)}\n"
                  f"other : {str(self.other)}")
        elif attr == "mandatory":
            print(attr+" : "+str(self.mandatory))
        elif attr == "accessory":
            print(attr+" : "+str(self.accessory))
        elif attr == "forbidden":
            print(attr+" : "+str(self.forbidden))
        elif attr == "other":
            print(attr+" : "+str(self.other))
        else:
            raise Exception("This option doesn't exist")

    # read json file given
    def read_json(self, json_file):
        with open(json_file) as my_file:
            data = json.load(my_file)
        # distribute values to attributes
        for key, value in data['families'].items():
            if key == 'mandatory':
                self.mandatory = set(value)
            elif key == 'accessory':
                self.accessory = set(value)
            elif key == 'forbidden':
                self.forbidden = set(value)
            else:
                self.other = set(value)
        for key, value in data['Parameters'].items():
            if key == 'min_mandatory':
                if len(self.mandatory) <= value:
                    print("min_mandatory is respected")
                else:
                    print("Error : min_mandatory isn't respected")
            if key == 'min_total':
                with value in data['families'].items():
                    if len(value):
                        print("min_total is respected")
                    else:
                        print("Error : min_total isn't respected")
            if key == 'max_forbidden':
                if len(self.forbidden) > value:
                    print("max_forbidden is respected")
                else:
                    print("Error : max_forbidden isn't respected")
            if key == 'max_other':
                if len(self.other) <= value:
                    print("max_other is respected")
                else:
                    print("Error : max_other isn't respected")


