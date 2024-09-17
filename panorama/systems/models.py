#!/usr/bin/env python3
# coding:utf-8

"""
This module provides functions to describe rules used to detect biological systems
"""

# default libraries
from __future__ import annotations
from pathlib import Path
from typing import Dict, List, Generator, Set, Tuple, Union
import json

# TODO Try to add those variable in class
suprules_params = ['min_mandatory', 'min_total']
keys_param = suprules_params + ['transitivity', 'window']
rule_keys = ['name', 'parameters', 'presence']
accept_type = ['mandatory', 'accessory', 'forbidden', 'neutral']


def check_key(parameter_data: Dict, need_key: List):
    """
    Check if all the keys are present in the dictionary.

    This function is applied to check model, functional unit and family consistency.

    Args:
        parameter_data (Dict): Dictionary which defines model, functional unit or family.
        need_key (List): List of all the keys needed.

    Raises:
        KeyError: If all the following keys are not present.
    """
    if not all([key in parameter_data for key in need_key]):
        raise KeyError(f"All the following keys are necessary: {need_key}")


def check_parameters(param_dict: Dict[str, int], mandatory_keys: List[str]):
    """
    Check if all parameters are inside and with the correct presence.

    Args:
        param_dict (Dict[str, int]): Dictionary with all the parameters for the rule.
        mandatory_keys (List[str]): List of the mandatory keys.

    Raises:
        KeyError: If one or more mandatory parameters are missing.
        TypeError: If one or more parameters are with an incorrect presence.
        ValueError: If one or more parameters are with an incorrect value.
        Exception: If an unexpected error occurs.
    """
    try:
        check_key(param_dict, mandatory_keys)
    except KeyError:
        raise KeyError("One or more attributes are missing in parameters")
    except Exception as error:
        raise Exception(f"Unexpected Error: {error} to get parameters")
    else:
        for key, value in param_dict.items():
            if key in ["min_mandatory", "min_total", "transitivity", 'window']:
                if not isinstance(value, int):
                    raise TypeError(f"The {key} value is not an integer")
                if value < -1:
                    raise ValueError(f"The {key} value is not positive. "
                                     "You can also use -1 value to skip the check on this parameter.")
            elif key in ["duplicate"]:
                if not isinstance(value, int):
                    raise TypeError(f"The {key} value is not an integer")
                if value < 0:
                    raise ValueError(f"The {key} value is not positive")
            elif key in ['multi_system', "multi_model"]:
                if not isinstance(value, bool):
                    raise TypeError("Overlap value from family in json must be a boolean")
            else:
                raise KeyError(f"{key} is not an acceptable attribute in parameters")


def check_dict(data_dict: Dict[str, Union[str, int, list, Dict[str, int]]], mandatory_keys: List[str],
               param_keys: List[str] = None):
    """
    Check if all keys and values are present and with the correct presence before adding to the model.

    Args:
        data_dict (Dict[str, Union[str, int, list, Dict[str, int]]]): Dictionary with all model information.
        mandatory_keys (List[str]): List of the mandatory keys.
        param_keys (List[str], optional): List of the mandatory keys for parameters. Defaults to None.

    Raises:
        KeyError: If one or more keys are missing or not acceptable.
        TypeError: If one or more values are not with the correct presence.
        ValueError: If one or more values are not acceptable.
        Exception: If an unexpected error occurs.
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
                    raise TypeError("The name value must be a string")
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
                    raise ValueError("Model needs at least one functional unit")
            elif key == 'families':
                if not isinstance(value, list):
                    raise TypeError("families value in json must be an array")
                if len(value) < 1:
                    raise ValueError("Functional unit needs at least one family")
            elif key == 'canonical':
                if not isinstance(value, list):
                    raise TypeError("canonical value in json must be an array")
            elif key == 'exchangeable':
                if not isinstance(value, list):
                    raise TypeError("exchangeable value from family in json must be a list")
                if not all(isinstance(elem, str) for elem in value):
                    raise ValueError("Exchangeable families must be a string")
            elif key == 'same_strand':
                if not isinstance(value, bool):
                    raise TypeError("same_strand in json must be a boolean")
            else:
                raise KeyError(f"{key} is not an acceptable attribute")


class Models:
    """
    A set of models defining a system.

    Args:
        models (Set[Model], optional): A set of models defining a system. Defaults to None.
    """

    def __init__(self, models: Set[Model] = None):
        """
        Constructor method to create Models object.

        Args:
            models (Set[Model], optional): A set of models defining a system. Defaults to None.
        """
        self._model_getter = models if models is not None else {}

    def __iter__(self) -> Generator[Model, None, None]:
        """
        Iterate over the models.

        Yields:
            Model: A model in the set.
        """
        for _model in self._model_getter.values():
            yield _model

    @property
    def value(self) -> List[Model]:
        """
        Return all models added. Useful if you need a list and not a generator.

        Returns:
            List[Model]: A list of all models added.
        """
        return list(self)

    @property
    def size(self) -> int:
        """
        Get the number of models added.

        Returns:
            int: The number of models inside.
        """
        return len(self.value)

    @property
    def func_units(self) -> Generator[FuncUnit, None, None]:
        """
        Get all functional units in models.

        Yields:
            FuncUnit: All functional units.
        """
        func_units = set()
        for model in self:
            for fu in model.func_units:
                func_units.add(fu)
        yield from func_units

    def func_units_to_model(self) -> Dict[FuncUnit, Model]:
        """
        Get all functional units in models and link them to the corresponding model.

        Returns:
            Dict[FuncUnit, Model]: All functional units linked to their model.
        """
        fu2model = {}
        for fu in self.func_units:
            fu2model[fu] = fu.model
        return fu2model

    @property
    def families(self) -> Generator[Family, None, None]:
        """
        Get all families in models.

        Yields:
            Family: All families in models.
        """
        families = set()
        for model in self:
            for family in model.families:
                families.add(family)
        yield from families

    def families_to_model(self) -> Dict[Family, Model]:
        """
        Get all families in models and link them to the corresponding model.

        Returns:
            Dict[Family, Model]: All families in models linked to their corresponding model.
        """
        fam2model = {}
        for family in self.families:
            fam2model[family] = family.model
        return fam2model

    def read(self, model_path: Path):
        """
        Read all JSON files models in the directory.

        Args:
            model_path (Path): Path to model.

        Raises:
            KeyError: If one or more keys are missing or not acceptable.
            TypeError: If one or more attributes are not with the correct presence.
            ValueError: If one or more attributes are not with an acceptable value.
            Exception: If an unexpected problem occurs to read JSON.
        """
        with open(model_path.resolve().as_posix()) as json_file:
            data = json.load(json_file)
            try:
                model = Model.read_model(data)
            except KeyError:
                raise KeyError(f"Problem with one or more keys in {model_path} are missing.")
            except TypeError:
                raise TypeError(f"One or more attributes are not with the correct presence in {model_path}.")
            except ValueError:
                raise ValueError(f"One or more attributes are not with an acceptable value in {model_path}.")
            except Exception:
                raise Exception(f"Unexpected problem to read JSON {model_path}")
            else:
                self.add_model(model)

    def get_model(self, name: str) -> Model:
        """
        Get a model by its name.

        Args:
            name (str): Name to find.

        Raises:
            KeyError: If the model is not present.

        Returns:
            Model: The model with the given name.
        """
        try:
            model = self._model_getter[name]
        except KeyError:
            raise KeyError("Model not present in set of value")
        else:
            return model

    def add_model(self, model: Model):
        """
        Add a model.

        Args:
            model (Model): Complete model object.

        Raises:
            Exception: If a model with the same name is already present in the system.
        """
        try:
            self.get_model(model.name)
        except KeyError:
            model.check_model()
            self._model_getter[model.name] = model
        else:
            raise Exception(f"Model {model.name} already in set of value")


class _BasicFeatures:
    """
    Basic features for Model, FuncUnit, and Family classes.

    Args:
        name (str, optional): Name of the element. Defaults to "".
        transitivity (int, optional): Size of the transitive closure used to build the graph. Defaults to 0.
        window (int, optional): Number of neighboring genes that are considered on each side of a gene of interest when searching for conserved genomic contexts. Defaults to 1.
    """

    def __init__(self, name: str = "", transitivity: int = 0, window: int = 1):
        """
        Constructor method to create a Basic features object

        Args:
            name (str, optional): Name of the element. Defaults to "".
            transitivity (int, optional): Size of the transitive closure used to build the graph. Defaults to 0.
            window (int, optional): Number of neighboring genes that are considered on each side of a gene of interest when searching for conserved genomic contexts. Defaults to 1.
        """
        self.name = name
        self.transitivity = transitivity
        self.window = window

    def __repr__(self):
        return f"{self.__class__.__name__} name: {self.name}"

    def __str__(self):
        return f"{self.__class__.__name__} name: {self.name}"

    def read_parameters(self, parameters: Dict[str, Union[str, int, bool]], param_keys: List[str]):
        """
        Check parameters' consistency.

        Raises:
            Exception: If the model is not consistent.
        """
        for param in param_keys:
            if param in parameters:
                self.__setattr__(param, parameters[param])
            else:
                if hasattr(self, '_parent'):
                    parent = self.__getattribute__('_parent')
                    if hasattr(parent, param):
                        self.__setattr__(param, parent.__getattribute__(param))


class _FuFamFeatures:
    """
    Features for FuncUnit and Family classes.

    Args:
        presence (str, optional): Type of the rule (mandatory, accessory, forbidden or neutral). Defaults to "".
        parent (Union[FuncUnit, Model], optional): Functional unit or model in which is the family. Defaults to None.
        duplicate (int, optional): Number of duplicates. Defaults to 0.
        exchangeable (Set[str], optional): List of exchangeable families. Defaults to None.
        multi_system (bool, optional): If the family can be present in multiple systems. Defaults to False.
        multi_model (bool, optional): If the family can be present in multiple models. Defaults to True.
    """

    def __init__(self, presence: str = "", parent: Union[FuncUnit, Model] = None, duplicate: int = 0,
                 exchangeable: Set[str] = None, multi_system: bool = False, multi_model: bool = True):
        """
        Constructor method to create a FuFamfeature object

        Args:
            presence (str, optional): Type of the rule (mandatory, accessory, forbidden or neutral). Defaults to "".
            parent (Union[FuncUnit, Model], optional): Functional unit or model in which is the family. Defaults to None.
            duplicate (int, optional): Number of duplicates. Defaults to 0.
            exchangeable (Set[str], optional): List of exchangeable families. Defaults to None.
            multi_system (bool, optional): If the family can be present in multiple systems. Defaults to False.
            multi_model (bool, optional): If the family can be present in multiple models. Defaults to True.
        """

        self.presence = presence
        self.duplicate = duplicate
        self.exchangeable = exchangeable if exchangeable is not None else set()
        self.multi_system = multi_system
        self.multi_model = multi_model

        self._parent = parent


class _ModFuFeatures:
    """
    Features for Model and FuncUnit classes.

    Args:
        mandatory (Set[FuncUnit, Family], optional): Set of mandatory sub-elements. Defaults to None.
        min_mandatory (int, optional): Minimum number of mandatory sub-elements. Defaults to 1.
        accessory (Set[FuncUnit, Family], optional): Set of accessory sub-elements. Defaults to None.
        min_total (int, optional): Minimum number of total sub-elements. Defaults to 1.
        forbidden (Set[FuncUnit, Family], optional): Set of forbidden sub-elements. Defaults to None.
        neutral (Set[FuncUnit, Family], optional): Set of neutral sub-elements. Defaults to None.
        same_strand (bool, optional): If the sub-elements must be on the same strand. Defaults to False.
    """

    def __init__(self, mandatory: Set[FuncUnit, Family] = None, accessory: Set[FuncUnit, Family] = None,
                 forbidden: Set[FuncUnit, Family] = None, neutral: Set[FuncUnit, Family] = None,
                 min_mandatory: int = 1, min_total: int = 1, same_strand: bool = False):
        """
        Constructor method to create a ModFuFeatures object.

        Args:
            mandatory (Set[FuncUnit, Family], optional): Set of mandatory sub-elements. Defaults to None.
            min_mandatory (int, optional): Minimum number of mandatory sub-elements. Defaults to 1.
            accessory (Set[FuncUnit, Family], optional): Set of accessory sub-elements. Defaults to None.
            min_total (int, optional): Minimum number of total sub-elements. Defaults to 1.
            forbidden (Set[FuncUnit, Family], optional): Set of forbidden sub-elements. Defaults to None.
            neutral (Set[FuncUnit, Family], optional): Set of neutral sub-elements. Defaults to None.
            same_strand (bool, optional): If the sub-elements must be on the same strand. Defaults to False.
        """
        self.mandatory = mandatory if mandatory is not None else set()
        self.min_mandatory = min_mandatory
        self.accessory = accessory if accessory is not None else set()
        self.min_total = min_total
        self.forbidden = forbidden if forbidden is not None else set()
        self.neutral = neutral if neutral is not None else set()
        self.same_strand = same_strand
        self._child_type = _FuFamFeatures
        self._child_getter = None

    @property
    def _children(self) -> Generator[Union[FuncUnit, Family], None, None]:
        """
        Get all child elements.

        Yields:
            Union[FuncUnit, Family]: All child elements.
        """
        for child in self.mandatory.union(self.accessory, self.forbidden, self.neutral):
            yield child

    def _child_names(self, presence: str):
        """
        Get all child elements names.

        Args:
            presence (str): Type of the rule (mandatory, accessory, forbidden or neutral).

        Returns:
            Set[str]: All child elements names.
        """
        if presence is None:
            return {child.name for child in self._children}
        else:
            return {child.name for child in self._children if child.presence == presence}

    def _duplicate(self, filter_type: str = None):
        """
        Access to all families that are duplicated in functional unit.

        Args:
            filter_type (str, optional): Type of the rule (mandatory, accessory, forbidden or neutral). Defaults to None.

        Yields:
            Union[FuncUnit, Family]: All families that are duplicated in functional unit.
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
        """
        Check model consistency.

        Raises:
            Exception: If the model is not consistent.
        """
        if self.min_mandatory > len(self.mandatory) + sum([child.duplicate for child in self._duplicate("mandatory")]):
            raise Exception(f"There are less mandatory {self.child_type} than the minimum mandatory")
        if self.min_total > len(list(self._children)) + sum([child.duplicate for child in self._duplicate()]):
            raise Exception(f"There are less {self.child_type} than the minimum total")
        if self.min_mandatory > self.min_total:
            raise Exception(f"Minimum mandatory {self.child_type} value is greater than minimum total.")
        if len(self.mandatory) == 0:
            raise Exception(f"There are no mandatory {self.child_type}. "
                            f"You should have at least one mandatory {self.child_type} with mandatory presence.")

    def add(self, child: Union[FuncUnit, Family]):
        """
        Add a functional unit to the model.

        Args:
            child (Union[FuncUnit, Family]): Functional unit or family.
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
        child._parent = self

    def _mk_child_getter(self):
        self._child_getter = {}
        for child in self._children:
            self._child_getter[child.name] = child

    def get(self, name: str) -> Union[FuncUnit, Family]:
        """
        Get a child from his name

        Args:
            name: name of the child to get.

        Returns:
            Union[FuncUnit, Family]: The child element

        Raises:
            KeyError: If the child is not found.
        """
        if self._child_getter is None:
            self._mk_child_getter()
        try:
            child = self._child_getter[name]
        except KeyError:
            raise KeyError(f"No such {self.child_type} with {name} in {type(self)}")
        else:
            return child

    @property
    def child_type(self):
        return "_FuFamFeature"


class Model(_BasicFeatures, _ModFuFeatures):
    """
    Represents Model rules which describe a biological system.

    Args:
        name (str, optional): Name of the element. Defaults to "".
        mandatory (Set[FuncUnit, Family], optional): Set of mandatory sub-elements. Defaults to None.
        min_mandatory (int, optional): Minimum number of mandatory sub-elements. Defaults to 1.
        accessory (Set[FuncUnit, Family], optional): Set of accessory sub-elements. Defaults to None.
        min_total (int, optional): Minimum number of total sub-elements. Defaults to 1.
        forbidden (Set[FuncUnit, Family], optional): Set of forbidden sub-elements. Defaults to None.
        neutral (Set[FuncUnit, Family], optional): Set of neutral sub-elements. Defaults to None.
        same_strand (bool, optional): If the sub-elements must be on the same strand. Defaults to False.
        transitivity (int, optional): Size of the transitive closure used to build the graph. Defaults to 0.
        window (int, optional): Number of neighboring genes that are considered on each side of a gene of interest when searching for conserved genomic contexts. Defaults to 1.
        canonical (list, optional): List of canonical models. Defaults to None.
    """

    def __init__(self, name: str = "", mandatory: Set[Union[FuncUnit, Family]] = None,
                 accessory: Set[Union[FuncUnit, Family]] = None, forbidden: Set[Union[FuncUnit, Family]] = None,
                 neutral: Set[Union[FuncUnit, Family]] = None, min_mandatory: int = 1, min_total: int = 1,
                 transitivity: int = 0, window: int = 1, same_strand: bool = False, canonical: list = None):
        """
        Constructor method to create a Model rules which describe a biological system.

        Args:
            name (str, optional): Name of the element. Defaults to "".
            mandatory (Set[FuncUnit, Family], optional): Set of mandatory sub-elements. Defaults to None.
            min_mandatory (int, optional): Minimum number of mandatory sub-elements. Defaults to 1.
            accessory (Set[FuncUnit, Family], optional): Set of accessory sub-elements. Defaults to None.
            min_total (int, optional): Minimum number of total sub-elements. Defaults to 1.
            forbidden (Set[FuncUnit, Family], optional): Set of forbidden sub-elements. Defaults to None.
            neutral (Set[FuncUnit, Family], optional): Set of neutral sub-elements. Defaults to None.
            same_strand (bool, optional): If the sub-elements must be on the same strand. Defaults to False.
            transitivity (int, optional): Size of the transitive closure used to build the graph. Defaults to 0.
            window (int, optional): Number of neighboring genes that are considered on each side of a gene of interest when searching for conserved genomic contexts. Defaults to 1.
            canonical (list, optional): List of canonical models. Defaults to None.
        """

        super().__init__(name=name, transitivity=transitivity, window=window)
        super(_BasicFeatures, self).__init__(mandatory=mandatory, min_mandatory=min_mandatory,
                                             accessory=accessory, neutral=neutral, min_total=min_total,
                                             forbidden=forbidden, same_strand=same_strand)
        self._child_type = FuncUnit
        self.canonical = canonical if canonical is not None else []

    def __hash__(self) -> int:
        return hash((self.name, self.min_mandatory, self.min_total, self.same_strand))

    def __eq__(self, other: Model) -> bool:
        if not isinstance(other, Model):
            raise TypeError(f"A Model object is expected, you give a {type(other)}")
        if hash(self) == hash(other):
            return set(self.func_units) == set(other.func_units)
        else:
            return False

    def __getstate__(self):
        state = self.__dict__.copy()
        return state

    def __setstate__(self, state):
        self.__dict__.update(state)
        for func_unit in self.func_units:
            func_unit.model = self

    @property
    def func_units(self) -> Generator[FuncUnit, None, None]:
        """
        Access to all functional units in models.

        Yields:
            FuncUnit: All functional units.
        """
        yield from self._children

    def func_units_names(self, presence: str):
        """
        Get all functional units names.

        Args:
            presence (str): Type of the rule (mandatory, accessory, forbidden or neutral).

        Returns:
            Set[str]: All functional units names.
        """
        return self._child_names(presence)

    @property
    def families(self) -> Generator[Family, None, None]:
        """
        Access to all families in models.

        Yields:
            Family: All families.
        """
        for func_unit in self.func_units:
            yield from func_unit.families

    @property
    def size(self) -> Tuple[int, int]:
        """
        Get the number of elements in the model.

        Returns:
            Tuple[int, int]: Number of functional units and number of families.
        """
        return len(list(self.func_units)), len(list(self.families))

    def duplicate_fu(self, filter_type: str = None):
        """
        Access to all families that are duplicated in functional unit.

        Args:
            filter_type (str, optional): Type of the rule (mandatory, accessory, forbidden or neutral). Defaults to None.

        Yields:
            Family: All families that are duplicated in functional unit.
        """
        yield from self._duplicate(filter_type)

    def get(self, name: str) -> Union[FuncUnit]:
        """
        Get a functional unit in the model with his name

        Args:
            name: name of the functional unit to get.

        Returns:
            FuncUnit: The functional unit element

        Raises:
            KeyError: If the functional unit with the name is not found.
        """
        return super().get(name)

    def check_model(self):
        """
        Check model consistency.

        Raises:
            Exception: If the model is not consistent.
        """
        try:
            self._check()
        except Exception as err:
            raise Exception(f"Consistency not respected in {self.name}. {err}")

    def read(self, data_model: dict):
        """
        Read model to parse into self attributes.

        Args:
            data_model (dict): JSON data dictionary.
        """
        mandatory_key = ['name', 'parameters', 'func_units']
        param_mandatory = ['transitivity', 'min_mandatory', 'min_total']

        check_dict(data_model, mandatory_keys=mandatory_key, param_keys=param_mandatory)

        self.name = data_model["name"]
        self.read_parameters(data_model["parameters"], param_keys=param_mandatory)
        self.window = data_model["window"] if "window" in data_model else self.transitivity + 1
        for dict_fu in data_model["func_units"]:
            func_unit = FuncUnit()
            func_unit.model = self
            func_unit.read(dict_fu)
            self.add(func_unit)
        if 'canonical' in data_model and data_model["canonical"] is not None:
            self.canonical = data_model["canonical"]

    @staticmethod
    def read_model(data_model: dict) -> Model:
        """
        Read model to parse into self attributes.

        Args:
            data_model (dict): JSON data dictionary.

        Returns:
            Model: Model object.
        """
        model = Model()
        model.read(data_model)
        return model

    @property
    def child_type(self):
        return "FuncUnit"


class FuncUnit(_BasicFeatures, _FuFamFeatures, _ModFuFeatures):
    """
    Represents functional unit definition rule.

    Args:
        name (str, optional): Name of the element. Defaults to "".
        presence (str, optional): Type of the rule (mandatory, accessory, forbidden or neutral). Defaults to "".
        mandatory (Set[FuncUnit, Family], optional): Set of mandatory sub-elements. Defaults to None.
        min_mandatory (int, optional): Minimum number of mandatory sub-elements. Defaults to 1.
        accessory (Set[FuncUnit, Family], optional): Set of accessory sub-elements. Defaults to None.
        min_total (int, optional): Minimum number of total sub-elements. Defaults to 1.
        forbidden (Set[FuncUnit, Family], optional): Set of forbidden sub-elements. Defaults to None.
        neutral (Set[FuncUnit, Family], optional): Set of neutral sub-elements. Defaults to None.
        same_strand (bool, optional): If the sub-elements must be on the same strand. Defaults to False.
        min_total (int, optional): Minimum number of total sub-elements. Defaults to 1.
        transitivity (int, optional): Size of the transitive closure used to build the graph. Defaults to 0.
        window (int, optional): Number of neighboring genes that are considered on each side of a gene of interest when searching for conserved genomic contexts. Defaults to 1.
        duplicate (int, optional): Number of duplicates. Defaults to 0.
        model (Model, optional): Model in which is the functional unit. Defaults to None.
        exchangeable (Set[str], optional): List of exchangeable families. Defaults to None.
        multi_system (bool, optional): If the functional unit can be present in multiple systems. Defaults to False.
        multi_model (bool, optional): If the functional unit can be present in multiple models. Defaults to False.
    """

    def __init__(self, name: str = "", presence: str = "", mandatory: Set[FuncUnit, Family] = None,
                 accessory: Set[FuncUnit, Family] = None, forbidden: Set[FuncUnit, Family] = None,
                 neutral: Set[FuncUnit, Family] = None, min_mandatory: int = 1, min_total: int = 1,
                 same_strand: bool = False, transitivity: int = 0, window: int = 1,
                 duplicate: int = 0, model: Model = None, exchangeable: Set[str] = None,
                 multi_system: bool = False, multi_model: bool = False):
        """
        Constructor method to create a functional unit definition rule.

        Args:
            name (str, optional): Name of the element. Defaults to "".
            presence (str, optional): Type of the rule (mandatory, accessory, forbidden or neutral). Defaults to "".
            mandatory (Set[FuncUnit, Family], optional): Set of mandatory sub-elements. Defaults to None.
            min_mandatory (int, optional): Minimum number of mandatory sub-elements. Defaults to 1.
            accessory (Set[FuncUnit, Family], optional): Set of accessory sub-elements. Defaults to None.
            min_total (int, optional): Minimum number of total sub-elements. Defaults to 1.
            forbidden (Set[FuncUnit, Family], optional): Set of forbidden sub-elements. Defaults to None.
            neutral (Set[FuncUnit, Family], optional): Set of neutral sub-elements. Defaults to None.
            same_strand (bool, optional): If the sub-elements must be on the same strand. Defaults to False.
            min_total (int, optional): Minimum number of total sub-elements. Defaults to 1.
            transitivity (int, optional): Size of the transitive closure used to build the graph. Defaults to 0.
            window (int, optional): Number of neighboring genes that are considered on each side of a gene of interest when searching for conserved genomic contexts. Defaults to 1.
            duplicate (int, optional): Number of duplicates. Defaults to 0.
            model (Model, optional): Model in which is the functional unit. Defaults to None.
            exchangeable (Set[str], optional): List of exchangeable families. Defaults to None.
            multi_system (bool, optional): If the functional unit can be present in multiple systems. Defaults to False.
            multi_model (bool, optional): If the functional unit can be present in multiple models. Defaults to False.
        """
        super().__init__(name=name, transitivity=transitivity, window=window)
        super(_BasicFeatures, self).__init__(presence=presence, duplicate=duplicate, parent=model,
                                             exchangeable=exchangeable, multi_system=multi_system,
                                             multi_model=multi_model)
        super(_FuFamFeatures, self).__init__(mandatory=mandatory, min_mandatory=min_mandatory,
                                             accessory=accessory, neutral=neutral, min_total=min_total,
                                             forbidden=forbidden, same_strand=same_strand)
        self._child_type = Family

    def __hash__(self) -> int:
        return hash((self.name, self.presence, self.min_mandatory, self.min_total, self.duplicate,
                     self.multi_system, self.multi_model, self.same_strand))

    def __eq__(self, other: FuncUnit) -> bool:
        if not isinstance(other, FuncUnit):
            raise TypeError(f"A FuncUnit object is expected, you give a {type(other)}")
        if hash(self) == hash(other):
            return set(self.families) == set(other.families)
        else:
            return False

    def __getstate__(self):
        state = self.__dict__.copy()
        state['_parent'] = None
        return state

    def __setstate__(self, state):
        self.__dict__.update(state)
        for family in self.families:
            family.func_unit = self

    @property
    def model(self) -> Model:
        """
        Get the model in which is the functional unit.

        Returns:
            Model: Model in which is the functional unit.
        """
        return self._parent

    @model.setter
    def model(self, model: Model):
        """
        Set the model in which is the functional unit.

        Args:
            model (Model): Model in which is the functional unit.
        """
        if not isinstance(model, Model):
            raise TypeError(f'A Model object is expected to be added to FuncUnit, you give {type(model)}')
        self._parent = model

    @model.deleter
    def model(self):
        """
        Delete the model in which is the functional unit.
        """
        del self._parent

    @property
    def families(self) -> Generator[Family, None, None]:
        """
        Access to all families in functional unit.

        Yields:
            Family: All families.
        """
        yield from self._children

    def get(self, name: str) -> Union[Family]:
        """
        Get a family in the functional unit with his name

        Args:
            name: name of the family to get.

        Returns:
            FuncUnit: The family element

        Raises:
            KeyError: If the family with the name is not found.
        """
        return super().get(name)

    def families_names(self, presence: str = None):
        """
        Get all families names.

        Args:
            presence (str, optional): Type of the rule (mandatory, accessory, forbidden or neutral). Defaults to None.

        Returns:
            Set[str]: All families names.
        """
        return self._child_names(presence)

    @property
    def size(self) -> int:
        """
        Get the number of families in the model.

        Returns:
            int: Number of families.
        """
        return len(list(self.families))

    def duplicate_fam(self, filter_type: str = None):
        """
        Access to all families that are duplicated in functional unit.

        Args:
            filter_type (str, optional): Type of the rule (mandatory, accessory, forbidden or neutral). Defaults to None.

        Yields:
            Family: All families that are duplicated in functional unit.
        """
        yield from self._duplicate(filter_type)

    def check_func_unit(self):
        """
        Check functional unit consistency.

        Raises:
            Exception: If the functional unit is not consistent.
        """
        try:
            self._check()
        except Exception:
            raise Exception(f"Consistency not respected in model {self.model.name} at functional unit {self.name}")

    def read(self, data_fu: dict):
        """
        Read functional unit.

        Args:
            data_fu (dict): Data JSON file of all functional units.
        """
        mandatory_key = ['name', 'families', 'presence']
        fu_params = ['duplicate', 'min_total', 'min_mandatory', 'transitivity', "multi_system", "multi_model"]
        check_dict(data_fu, mandatory_keys=mandatory_key)

        self.name = data_fu["name"]
        self.presence = data_fu["presence"]
        self.read_parameters(data_fu["parameters"] if "parameters" in data_fu else {}, param_keys=fu_params)
        self.window = data_fu["window"] if "window" in data_fu else self.transitivity + 1
        for fam_dict in data_fu["families"]:
            family = Family()
            family.func_unit = self
            family.read(fam_dict)
            self.add(family)

    @staticmethod
    def read_func_unit(data_fu: dict) -> FuncUnit:
        """
        Read functional unit.

        Args:
            data_fu (dict): Data JSON file of all functional units.

        Returns:
            FuncUnit: Functional unit object.
        """
        func_unit = FuncUnit()
        func_unit.read(data_fu)
        return func_unit

    @property
    def child_type(self):
        return "Family"


class Family(_BasicFeatures, _FuFamFeatures):
    """
    Represents family model definition rule.

    Args:
        name (str, optional): Name of the element. Defaults to "".
        transitivity (int, optional): Size of the transitive closure used to build the graph. Defaults to 0.
        window (int, optional): Number of neighboring genes that are considered on each side of a gene of interest when searching for conserved genomic contexts. Defaults to 1.
        presence (str, optional): Type of the rule (mandatory, accessory, forbidden or neutral). Defaults to "".
        func_unit (FuncUnit, optional): Functional unit in which is the family. Defaults to None.
        duplicate (int, optional): Number of duplicates. Defaults to 0.
        exchangeable (Set[str], optional): List of exchangeable families. Defaults to None.
        multi_system (bool, optional): If the family can be present in multiple systems. Defaults to False.
        multi_model (bool, optional): If the family can be present in multiple models. Defaults to False.
    """

    def __init__(self, name: str = "", transitivity: int = 0, window: int = 1, presence: str = "",
                 func_unit: FuncUnit = None,
                 duplicate: int = 0, exchangeable: Set[str] = None, multi_system: bool = False,
                 multi_model: bool = False):
        """
        Constructor method to create a family model definition rule.

        Args:
            name (str, optional): Name of the element. Defaults to "".
            transitivity (int, optional): Size of the transitive closure used to build the graph. Defaults to 0.
            window (int, optional): Number of neighboring genes that are considered on each side of a gene of interest when searching for conserved genomic contexts. Defaults to 1.
            presence (str, optional): Type of the rule (mandatory, accessory, forbidden or neutral). Defaults to "".
            func_unit (FuncUnit, optional): Functional unit in which is the family. Defaults to None.
            duplicate (int, optional): Number of duplicates. Defaults to 0.
            exchangeable (Set[str], optional): List of exchangeable families. Defaults to None.
            multi_system (bool, optional): If the family can be present in multiple systems. Defaults to False.
            multi_model (bool, optional): If the family can be present in multiple models. Defaults to False.
        """
        super().__init__(name=name, transitivity=transitivity, window=window)
        super(_BasicFeatures, self).__init__(presence=presence, duplicate=duplicate, parent=func_unit,
                                             exchangeable=exchangeable, multi_system=multi_system,
                                             multi_model=multi_model)

    def __hash__(self) -> int:
        return hash((self.name, self.presence, self.transitivity, self.window, self.duplicate,
                     self.multi_model, self.multi_system))

    def __eq__(self, other: Family) -> bool:
        if not isinstance(other, Family):
            raise TypeError("A Family instance can only be compared to another Family instance.")
        return (self.name == other.name and
                self.presence == other.presence and
                self.transitivity == other.transitivity and
                self.window == other.window and
                self.duplicate == other.duplicate and
                self.multi_model == other.multi_model and
                self.multi_system == other.multi_system and
                self.exchangeable == other.exchangeable)

    def __getstate__(self):
        state = self.__dict__.copy()
        state['_parent'] = None
        return state

    def __setstate__(self, state):
        self.__dict__.update(state)

    @property
    def func_unit(self) -> FuncUnit:
        """
        Get the functional unit in which is the family.

        Returns:
            FuncUnit: Functional unit in which is the family.
        """
        return self._parent

    @func_unit.setter
    def func_unit(self, func_unit: FuncUnit):
        """
        Set the functional unit in which is the family.

        Args:
            func_unit (FuncUnit): Functional unit in which is the family.
        """
        if not isinstance(func_unit, FuncUnit):
            raise TypeError(f"An object of type FuncUnit is expected, you give a {type(func_unit)}")
        self._parent = func_unit

    @func_unit.deleter
    def func_unit(self):
        """
        Delete the functional unit in which is the family.
        """
        del self._parent

    @property
    def model(self) -> Model:
        """
        Get the model in which is the family.

        Returns:
            Model: Model in which is the family.
        """
        return self.func_unit.model

    def read(self, data_fam: dict):
        """
        Read family.

        Args:
            data_fam (dict): Data JSON file with families.
        """
        fam_param = ['transitivity', "duplicate", "multi_system", "multi_model"]

        check_dict(data_fam, mandatory_keys=['name', 'presence'])
        self.name = data_fam['name']
        self.presence = data_fam['presence']
        self.read_parameters(data_fam["parameters"] if "parameters" in data_fam else {}, param_keys=fam_param)
        self.window = data_fam["window"] if "window" in data_fam else self.transitivity + 1
        if 'exchangeable' in data_fam:
            self.exchangeable = set(data_fam["exchangeable"])

    @staticmethod
    def read_family(data_fam: dict) -> Family:
        """
        Read family.

        Args:
            data_fam (dict): Data JSON file with families.

        Returns:
            Family: Family object.
        """
        fam = Family()
        fam.read(data_fam)
        return fam
