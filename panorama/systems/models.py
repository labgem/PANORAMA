#!/usr/bin/env python3
# coding:utf-8

"""
This module provides tools to define and validate rules used to detect biological systems.
"""

from __future__ import annotations
from pathlib import Path
from typing import Dict, List, Generator, Set, Tuple, Union, Iterable
import json

# Constants (could be encapsulated later in a config class)
SUPRULES_PARAMS = ["min_mandatory", "min_total"]
PARAM_KEYS = SUPRULES_PARAMS + ["transitivity", "window"]
RULE_KEYS = ["name", "parameters", "presence"]
ACCEPTED_PRESENCE = ["mandatory", "accessory", "forbidden", "neutral"]


def check_key(data: Dict, required_keys: Set[str]) -> None:
    """
    Validates that all required keys are present in the given dictionary.

    Args:
        data (Dict): The dictionary to validate.
        required_keys (Set[str]): Set of keys that must be present in the dictionary.

    Raises:
        KeyError: If any required keys are missing from the dictionary.
    """

    missing_keys = required_keys - set(data.keys())

    if missing_keys:
        raise KeyError(
            f"the following keys are required: {', '.join(sorted(missing_keys))}"
        )


def check_parameters(param_dict: Dict[str, int], mandatory_keys: Set[str]) -> None:
    """
    Validates the parameters in the given dictionary against the required keys and their rules.

    This function ensures that all the required keys are present in the dictionary and that their
    values adhere to the specified type and value constraints. If any rule is violated, an exception
    is raised. The function is designed to handle specific parameter keys with associated rules,
    raising meaningful errors for incorrect usage.

    Args:
        param_dict (Dict[str, int]): The dictionary containing parameter keys and their values.
        mandatory_keys (List[str]): The list of keys that are required to be present in param_dict.

    Raises:
        KeyError: If one or more required keys are missing from param_dict or an unexpected
            parameter key is found.
        TypeError: If a parameter value in param_dict does not match the expected type.
        ValueError: If a parameter value in param_dict does not meet the defined constraints.
        Exception: If an unexpected error occurs during parameter validation.

    Todo:
        look at dataclasses to encapsulate validation rules (note ask Jerome for the example code)
    """
    try:
        check_key(param_dict, mandatory_keys)
    except KeyError:
        raise KeyError("One or more required parameters are missing.")
    except Exception as e:
        raise Exception(f"Unexpected error checking parameters: {e}")

    for key, value in param_dict.items():
        if key in SUPRULES_PARAMS + ["transitivity", "window"]:
            if not isinstance(value, int):
                raise TypeError(f"{key} must be an integer.")
            if value < -1:
                raise ValueError(f"{key} must be >= -1.")
        elif key == "duplicate":
            if not isinstance(value, int):
                raise TypeError("duplicate must be a non-negative integer.")
            if value < 0:
                raise ValueError("duplicate must be a non-negative integer.")
        elif key in ["multi_system", "multi_model"]:
            if not isinstance(value, bool):
                raise TypeError(f"{key} must be a boolean.")
        else:
            raise KeyError(f"Unexpected parameter key: {key}")


def check_dict(
    data_dict: Dict[str, Union[str, int, list, Dict[str, int]]],
    mandatory_keys: Set[str],
    param_keys: Set[str] = None,
) -> None:
    """
    Performs validation of the provided dictionary by checking for required keys,
    key types, value types, and specific constraints based on the content of the
    dictionary.

    This function takes a data dictionary and validates its structure and contents
    against the provided required keys and optional parameter keys. Validation checks
    include ensuring required keys exist, verifying data types for key values, and
    matching key values against accepted criteria. Raises specific exceptions if
    any checks fail.

    Args:
        data_dict (Dict[str, Union[str, int, list, Dict[str, int]]]): The dictionary
            that needs to be validated.
        mandatory_keys (List[str]): A list of keys that must be present in the
            dictionary.
        param_keys (List[str], optional): An optional list of keys required in the
            'parameters' field of the dictionary, if present. Defaults to an empty list.

    Raises:
        KeyError: If required top-level keys are missing, or unexpected keys are found
            in the dictionary.
        TypeError: If a field contains a value of an unexpected type.
        ValueError: If a field contains a value that does not meet the required
            constraints.
        Exception: If any unexpected issues occur during the validation process.
    Todo:
        split in multiple functions: one to validate type, one to validate value, and keep this function at the main one.
        (ask Jerome for the example code)
    """
    param_keys = set() if param_keys is None else param_keys

    try:
        check_key(data_dict, mandatory_keys)
    except KeyError:
        raise KeyError("One or more required top-level keys are missing.")
    except Exception as error:
        raise Exception(f"Unexpected error during validation: {error}")

    for key, value in data_dict.items():
        if key == "name":
            if not isinstance(value, str):
                raise TypeError("The 'name' field must be a string.")
        elif key == "presence":
            if not isinstance(value, str):
                raise TypeError("The 'presence' field must be a string.")
            if value not in ACCEPTED_PRESENCE:
                raise ValueError(f"'presence' must be one of {ACCEPTED_PRESENCE}.")
        elif key == "parameters" and value is not None:
            check_parameters(value, param_keys)
        elif key in ["func_units", "families"]:
            if not isinstance(value, list):
                raise TypeError(f"{key} must be a list.")
            if not value:
                raise ValueError(f"{key} must be a non-empty list.")
        elif key in ["canonical", "exchangeable"]:
            if not isinstance(value, list):
                raise TypeError(f"{key} must be a list.")
            if key == "exchangeable" and not all(isinstance(v, str) for v in value):
                raise ValueError("All 'exchangeable' values must be strings.")
        elif key == "same_strand":
            if not isinstance(value, bool):
                raise TypeError("'same_strand' must be a boolean.")
        elif key not in mandatory_keys:
            raise KeyError(f"Unexpected key found: {key}")


class Models:
    """
    Represents a collection of models and provides interfaces for interaction.

    This class is designed to manage a collection of `Model` objects. It
    provides functionality for iterating over the models, accessing model
    details, generating mappings between function units and their models,
    and handling families associated with models. Additionally, it allows
    reading model configurations from JSON files and adding new models to
    the collection.
    """

    def __init__(self, models: Set[Model] = None):
        """
        Initializes an instance of the class.

        Args:
            models (Set[Model], optional): A set of models to be used. Defaults to None.
        """
        self._model_getter = {model.name: model for model in models} if models is not None else {}

    def __iter__(self) -> Generator[Model, None, None]:
        """
        Yields Model instances from the internal model getter.

        This method provides an iterator over the `Model` objects
        retrieved from the `self._model_getter`. Each model instance
        is yielded one at a time.

        Yields:
            Model: The next model instance in the sequence.
        """
        for _model in self._model_getter.values():
            yield _model

    @property
    def value(self) -> List[Model]:
        """
        Gets the list of Model instances that are currently available.

        This property provides access to all the Model instances contained within the object.

        Returns:
            List[Model]: A list of Model instances.
        """
        return list(self)

    @property
    def size(self) -> int:
        """
        Gets the size of the `value` attribute.

        This property computes the size of the `value` attribute by returning its length.

        Returns:
            int: The length of the `value` attribute.
        """
        return len(self.value)

    @property
    def func_units(self) -> Generator[FuncUnit, None, None]:
        """
        Gets a generator of distinct functional units (FuncUnit) across all models.

        Yields:
            Generator[FuncUnit, None, None]: A generator that iterates over unique
            functional units found within all models.
        """
        func_units = set()
        for model in self:
            for fu in model.func_units:
                func_units.add(fu)
        yield from func_units

    def func_units_to_model(self) -> Dict[FuncUnit, Model]:
        """
        Creates a mapping between func units and their corresponding models.

        This method iterates through the collection of function units and creates a
        dictionary where each function unit is associated with its corresponding
        model. It provides a convenient structure to access models by their related
        function units.

        Returns:
            Dict[FuncUnit, Model]: A dictionary mapping function units to their
            respective models.
        """
        fu2model = {}
        for fu in self.func_units:
            fu2model[fu] = fu.model
        return fu2model

    @property
    def families(self) -> Generator[Family, None, None]:
        """
        Provides a generator that yields unique Family instances from models.

        This property iterates over all models, collects all unique Family
        instances associated with the models, and yields them one by one.

        Yields:
            Generator[Family, None, None]: A generator of Family instances.

        """
        families = set()
        for model in self:
            for family in model.families:
                families.add(family)
        yield from families

    def families_to_model(self) -> Dict[Family, Model]:
        """
        Generates a mapping of Family objects to their corresponding Model objects.

        Constructs a dictionary mapping each Family object in the `families`
        attribute to its associated Model object, derived from the `model`
        attribute of each Family.

        Returns:
            Dict[Family, Model]: A dictionary where keys are Family objects
            and values are their corresponding Model objects.
        """
        fam2model = {}
        for family in self.families:
            print(family, family.model)
            fam2model[family] = family.model
        return fam2model

    def get_model(self, name: str) -> Model:
        """
        Retrieves a model from the internal collection based on its name.

        This method attempts to fetch a model from an internal mapping using the
        provided name. If the name does not exist in the mapping, a KeyError is
        raised. On success, it returns the corresponding model.

        Args:
            name (str): The name of the model to retrieve.

        Returns:
            Model: The model corresponding to the provided name.

        Raises:
            KeyError: If the provided name does not exist in the internal mapping.
        """
        try:
            model = self._model_getter[name]
        except KeyError:
            raise KeyError("Model not present in set of value")
        else:
            return model

    def add_model(self, model: Model):
        """
        Adds a new model to the collection if it does not already exist.

        This method allows adding a model to the collection, provided that no model
        with the same name exists in the current collection. If a model with the
        same name is already present, a KeyError will be raised. Before adding
        the model, it ensures the validity of the model by invoking its `check_model`
        method.

        Args:
            model (Model): The model instance to be added to the collection.

        Raises:
            KeyError: If a model with the same name is already present in the collection.
            Any exception raised by model.check_model() if the model is invalid.
        """
        try:
            self.get_model(model.name)
        except KeyError:
            model.check_model()
            self._model_getter[model.name] = model
        else:
            raise KeyError(f"Model {model.name} already in set of value")

    def read(self, model_path: Path):
        """
        Reads a model configuration from a JSON file, processes it, and adds the parsed model.

        This function attempts to read the provided file path as a JSON file, extract model
        configuration data from it, and handle any parsing errors that may occur. Successful
        processing results in the extracted model being added to the system.

        Args:
            model_path (Path): Path to the configuration JSON file to be read.

        Raises:
            KeyError: If one or more required keys are missing in the JSON file.
            TypeError: If one or more attributes in the JSON file are not correctly structured.
            ValueError: If one or more attributes have unacceptable values in the JSON file.
            Exception: For any unexpected issues encountered while reading the JSON file.
        """
        with open(model_path.resolve().as_posix()) as json_file:
            data = json.load(json_file)
            try:
                model = Model.read_model(data)
            except KeyError:
                raise KeyError(
                    f"Problem with one or more keys in {model_path} are missing."
                )
            except TypeError:
                raise TypeError(
                    f"One or more attributes are not with the correct presence in {model_path}."
                )
            except ValueError:
                raise ValueError(
                    f"One or more attributes are not with an acceptable value in {model_path}."
                )
            except Exception:
                raise Exception(f"Unexpected problem to read JSON {model_path}")
            else:
                self.add_model(model)

class _BasicFeatures:
    """
    Handles basic features typically required for configurations.

    This class is used as a foundation for managing basic configurations, ensuring
    parameter consistency, and providing utility methods for its derived classes. It includes
    mechanisms for naming, transitivity, and contextual window definitions. The class also
    provides methods to override and validate parameters dynamically by integrating it with
    external sources or inheriting configurations from parent instances.

    Attributes:
        name (str): Name of the element.
        transitivity (int): Size of the transitive closure used to build the graph.
        window (int): Number of neighboring genes considered on each side of a gene of
        interest when searching for conserved genomic contexts.
    """

    def __init__(self, name: str = "", transitivity: int = 0, window: int = 1):
        """
        Initializes the class with specific attributes for name, transitivity, and
        window. These attributes are used to define the characteristics and behavior
        of an object of this class.

        Args:
            name (str): The name associated with the object. Defaults to an
                empty string.
            transitivity (int): An integer value representing the transitivity
                associated with the object. Defaults to 0.
            window (int): An integer that defines the window size or configuration
                parameter. Defaults to 1.
        """
        self.name = name
        self.transitivity = transitivity
        self.window = window

    def __repr__(self):
        """
        Provides a string representation of the object. This method is intended to
        provide a clear and concise human-readable representation of the object in
        the form of its class name and name attribute.

        Returns:
            str: A string containing the class name and the object's name attribute.
        """
        return f"{self.__class__.__name__} name: {self.name}"

    def __str__(self):
        """
        Converts the object representation to its string form.

        This method is used to return a string representation of the class instance,
        primarily for debugging or logging purposes. It includes the name of the
        class and the value of the `name` attribute.

        Returns:
            str: A formatted string containing the class name and the `name` attribute.
        """
        return f"{self.__class__.__name__} name: {self.name}"

    def read_parameters(
        self, parameters: Dict[str, Union[str, int, bool]], param_keys: Set[str]
    ):
        """
        Reads and assigns parameters from a provided dictionary or falls back to parent attributes
        if the parameter is not found.

        Args:
            parameters (Dict[str, Union[str, int, bool]]): A dictionary containing parameter key-value pairs.
            param_keys (List[str]): A list of keys to be retrieved from the parameters' dictionary.

        """
        for param in param_keys:
            if param in parameters:
                self.__setattr__(param, parameters[param])
            else:
                if hasattr(self, "_parent"):
                    parent = self.__getattribute__("_parent")
                    if hasattr(parent, param):
                        self.__setattr__(param, parent.__getattribute__(param))


class _FuFamFeatures:
    """Represents features related to functional family rules and their properties.

    This class is used to define and manage the attributes and rules associated
    with a functional family. It includes properties such as presence types,
    duplication count, exchangeable families, and whether a family can be present
    in multiple systems or models.

    Attributes:
        presence (str): Type of the rule, such as mandatory, accessory, forbidden,
            or neutral.
        duplicate (int): Number of duplicates for the functional family.
        exchangeable (set): Set of exchangeable families that can replace one
            another in functionality.
        multi_system (bool): Determines if the family can exist across multiple
            systems.
        multi_model (bool): Determines if the family can exist across multiple
            models.
    """

    def __init__(
        self,
        presence: str = "",
        parent: Union[FuncUnit, Model] = None,
        duplicate: int = 0,
        exchangeable: Set[str] = None,
        multi_system: bool = False,
        multi_model: bool = True,
    ):
        """
        Initializes the instance of the class with given parameters.

        Args:
            presence (str): Describes the presence attribute, default is an empty string.
            parent (Union[FuncUnit, Model]): The parent object which can be of type FuncUnit
                or Model. Defaults to None.
            duplicate (int): Specifies the count of duplicates. Defaults to 0.
            exchangeable (Set[str]): A set of values indicating exchangeable items. Defaults
                to an empty set if not provided.
            multi_system (bool): Indicates if the functionality spans across multiple systems.
                Defaults to False.
            multi_model (bool): Indicates if the functionality spans across multiple models.
                Defaults to True.
        """

        self.presence = presence
        self.duplicate = duplicate
        self.exchangeable = exchangeable if exchangeable is not None else set()
        self.multi_system = multi_system
        self.multi_model = multi_model

        self._parent = parent


class _ModFuFeatures:
    """
    Represents the features of functional units or families in a model.

    This class is responsible for managing sets of child elements categorized into
    mandatory, accessory, forbidden, and neutral groups. It ensures consistency
    by enforcing rules about the minimum number of mandatory elements, total
    elements, and whether the elements belong to the same strand. The class also
    provides methods to access, validate, and manipulate the functional units or
    families.

    Attributes:
        mandatory (Set[FuncUnit, Family]): Set of mandatory child elements for the model.
        min_mandatory (int): Minimum number of mandatory child elements required.
        accessory (Set[FuncUnit, Family]): Set of accessory child elements for the model.
        min_total (int): Minimum total number of child elements required.
        forbidden (Set[FuncUnit, Family]): Set of forbidden child elements for the model.
        neutral (Set[FuncUnit, Family]): Set of neutral child elements for the model.
        same_strand (bool): Indicates whether all child elements must reside on the same strand.
    """

    def __init__(
        self,
        mandatory: Set[FuncUnit, Family] = None,
        accessory: Set[FuncUnit, Family] = None,
        forbidden: Set[FuncUnit, Family] = None,
        neutral: Set[FuncUnit, Family] = None,
        min_mandatory: int = 1,
        min_total: int = 1,
        same_strand: bool = False,
    ):
        """
        Initializes an object with specified sets of functional units or families categorized as
        mandatory, accessory, forbidden, and neutral. It also allows configuration of constraints
        on minimum requirements and positional restrictions.

        Args:
            mandatory (Set[FuncUnit, Family], optional): A set of functional units or families
                that must be included. Defaults to an empty set if not specified.
            accessory (Set[FuncUnit, Family], optional): A set of functional units or families
                that are optional to include. Defaults to an empty set if not specified.
            forbidden (Set[FuncUnit, Family], optional): A set of functional units or families
                that must not be included. Defaults to an empty set if not specified.
            neutral (Set[FuncUnit, Family], optional): A set of functional units or families
                that are neutral in relevance. Defaults to an empty set if not specified.
            min_mandatory (int, optional): The minimum number of functional units or families
                that must be included from the mandatory set. Defaults to 1.
            min_total (int, optional): The minimum total number of functional units or families
                that must be included overall. Defaults to 1.
            same_strand (bool, optional): A flag indicating whether all included functional units
                or families must be on the same strand. Defaults to False.
        """
        self.mandatory = mandatory if mandatory is not None else set()
        self.min_mandatory = min_mandatory
        self.accessory = accessory if accessory is not None else set()
        self.min_total = min_total
        self.forbidden = forbidden if forbidden is not None else set()
        self.neutral = neutral if neutral is not None else set()
        self.same_strand = same_strand
        self._child_type = None
        self._child_getter = None

    @property
    def _children(self) -> Generator[Union[FuncUnit, Family], None, None]:
        """
        Provides an iterator for child elements, including mandatory, accessory,
        forbidden, and neutral elements.

        Yields:
            Union[FuncUnit, Family]: An instance representing a child element from
            the combined set of mandatory, accessory, forbidden, and neutral elements.
        """
        for child in self.mandatory.union(self.accessory, self.forbidden, self.neutral):
            yield child

    def _child_names(self, presence: str = None):
        """
        Retrieves the names of child objects based on their presence status.

        Args:
            presence: A string representing the presence status of child objects
                to filter (e.g., "active", "inactive"). If None, retrieves names
                of all child objects.

        Returns:
            set: A set containing the names of child objects that match the given
                presence status, or all child names if presence is None.
        """
        if presence is None:
            return {child.name for child in self._children}
        else:
            return {
                child.name for child in self._children if child.presence == presence
            }

    def _duplicate(self, filter_type: str = None):
        """
        Filters and yields child elements based on a specific filter type.

        This function allows filtering of child elements according to a specified
        attribute type, such as 'mandatory', 'accessory', 'forbidden', or 'neutral'.
        If no filter type is provided, it selects all available children. It iterates
        through the filtered children and yields those with a `duplicate` value of
        1 or higher.

        Args:
            filter_type (str, optional): The type of filter to apply to the children.
                Acceptable values are 'mandatory', 'accessory', 'forbidden', 'neutral',
                or None. Defaults to None.

        Yields:
            object: Each child element that matches the filter condition and has a
                `duplicate` value of 1 or higher.
        """
        assert filter_type in [None, "mandatory", "accessory", "forbidden", "neutral"]
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
        Validates the configuration of mandatory and total required elements for a specific type.
        Ensures that the mandatory elements meet the minimum requirements, both in count and presence,
        and that the total elements conform to the defined minimum limits.

        Raises:
            Exception: If the number of mandatory elements is lower than the required minimum.
            Exception: If the number of total elements is lower than the required minimum.
            Exception: If the minimum mandatory count exceeds the minimum total count.
            Exception: If no mandatory elements are present.

        """
        if len(self.mandatory) == 0:
            raise Exception(
                f"There are no mandatory {self.child_type}. "
                f"You should have at least one mandatory {self.child_type} with mandatory presence."
            )
        if self.min_mandatory > len(self.mandatory) + sum(
            [child.duplicate for child in self._duplicate("mandatory")]
        ):
            raise Exception(
                f"There are less mandatory {self.child_type} than the minimum mandatory"
            )
        if self.min_total > len(self.mandatory.union(self.accessory)) + sum(
            [
                child.duplicate
                for child in set(self._duplicate("mandatory")).union(
                    set(self._duplicate("accessory"))
                )
            ]
        ):
            raise Exception(f"There are less {self.child_type} than the minimum total")
        if self.min_mandatory > self.min_total:
            raise Exception(
                f"Minimum mandatory {self.child_type} value is greater than minimum total."
            )

    @property
    def child_type(self):
        """
        Determines the consistent child type for the instance.

        Raises:
            Exception: If the child type among children is inconsistent.

        Returns:
            Any: The consistent child type of the instance.
        """
        if self._child_type is not None:
            return self._child_type
        else:
            child_type = None
            for child in self._children:
                curr_child_type = type(child).__name__
                if child_type is None:
                    child_type = curr_child_type
                elif child_type != curr_child_type:
                    raise TypeError(
                        f"The child type is inconsistent. "
                        f"It contains {child_type} and {curr_child_type}"
                    )
            self._child_type = child_type
            return self._child_type

    def add(self, child: Union[FuncUnit, Family]):
        """
        Adds a child element to one of the sets in the instance based on its `presence`
        attribute. The child can either be of the type `FuncUnit` or `Family`. The `presence`
        attribute of the child determines which set (`mandatory`, `accessory`, `forbidden`,
        or `neutral`) the child will be added to.

        Args:
            child (Union[FuncUnit, Family]): The element to be added. The element must
                have a `presence` attribute which determines the specific set it belongs
                to. If it is of the type `FuncUnit`, its `check_func_unit` method will be
                invoked before proceeding.

        """
        if set(self._children) and type(child).__name__ != self.child_type:
            raise TypeError(
                f"The child type is inconsistent. Expected {self.child_type} but found {type(child).__name__}."
            )
        if child.presence == "mandatory":
            self.mandatory.add(child)
        elif child.presence == "accessory":
            self.accessory.add(child)
        elif child.presence == "forbidden":
            self.forbidden.add(child)
        elif child.presence == "neutral":
            self.neutral.add(child)
        else:
            raise ValueError(
                f"The child {child.name} does not have a valid presence attribute ({child.presence})."
            )
        child._parent = self

    def _mk_child_getter(self):
        """
        Creates a dictionary for retrieving child objects by their names.

        This method iterates over the `_children` attribute and maps each child's
        name to the child object itself, storing these key-value pairs in the
        `_child_getter` dictionary. This allows efficient access to child objects
        later based on their names.
        """
        self._child_getter = {}
        for child in self._children:
            self._child_getter[child.name] = child

    def get(self, name: str) -> Union[FuncUnit, Family]:
        """
        Retrieves a child object associated with the given name.

        Args:
            name (str): The name of the child object to retrieve.

        Returns:
            Union[FuncUnit, Family]: The child object associated with the given name.

        Raises:
            KeyError: If no child object is found for the given name.
        """
        if self._child_getter is None:
            self._mk_child_getter()
        try:
            child = self._child_getter[name]
        except KeyError:
            raise KeyError(f"No such {self._child_type} with name {name} in {type(self)}")
        else:
            return child


class Model(_BasicFeatures, _ModFuFeatures):
    """
    A Model class representing rules which describe a biological system.

    This class provides an abstraction to handle a set of biological rules
    and their associated functional units and families. It supports operations
    to query various characteristics of the biological model, including its size,
    functional units, and families. The Model can also perform checks for consistency
    and parse data from a given dictionary for initialization and configuration.

    Attributes:
            name (str): Name of the instance.
            mandatory (Set[FuncUnit, Family]): Set of mandatory functional units or families.
            accessory (Set[FuncUnit, Family]): Set of accessory functional units or families.
            forbidden (Set[FuncUnit, Family]): Set of forbidden functional units or families.
            neutral (Set[FuncUnit, Family]): Set of neutral functional units or families.
            min_mandatory (int): Minimum number of mandatory functional units required.
            min_total (int): Minimum total number of functional elements required.
            transitivity (int): Transitivity value for relationships in functional elements.
            window (int): Window size used for analysis or traversal.
            same_strand (bool): Indicates if the functional units should be on the same strand.
            canonical (list): List containing canonical representations.
    """

    def __init__(
        self,
        name: str = "",
        mandatory: Set[FuncUnit, Family] = None,
        accessory: Set[FuncUnit, Family] = None,
        forbidden: Set[FuncUnit, Family] = None,
        neutral: Set[FuncUnit, Family] = None,
        min_mandatory: int = 1,
        min_total: int = 1,
        transitivity: int = 0,
        window: int = 1,
        same_strand: bool = False,
        canonical: list = None,
    ):
        """
        Initializes attributes of the class and sets up the properties related to the functional units.

        Args:
            name (str): Name of the instance.
            mandatory (Set[FuncUnit, Family]): Set of mandatory functional units or families.
            accessory (Set[FuncUnit, Family]): Set of accessory functional units or families.
            forbidden (Set[FuncUnit, Family]): Set of forbidden functional units or families.
            neutral (Set[FuncUnit, Family]): Set of neutral functional units or families.
            min_mandatory (int): Minimum number of mandatory functional units required.
            min_total (int): Minimum total number of functional elements required.
            transitivity (int): Transitivity value for relationships in functional elements.
            window (int): Window size used for analysis or traversal.
            same_strand (bool): Indicates if the functional units should be on the same strand.
            canonical (list): List containing canonical representations.
        """

        super().__init__(name=name, transitivity=transitivity, window=window)
        super(_BasicFeatures, self).__init__(
            mandatory=mandatory,
            min_mandatory=min_mandatory,
            accessory=accessory,
            neutral=neutral,
            min_total=min_total,
            forbidden=forbidden,
            same_strand=same_strand,
        )
        self.canonical = canonical if canonical is not None else []

    @property
    def func_units(self) -> Generator[FuncUnit, None, None]:
        """
        Retrieves generator of functional units.

        Returns:
            Generator[FuncUnit, None, None]: A generator yielding functional units from
                the collection of child elements.
        """
        yield from self._children

    def func_units_names(self, presence: str):
        """
        Returns the names of child units based on their presence.

        Args:
            presence (str): The presence status used for filtering child units.

        Returns:
            Any: A list or collection of child unit names filtered by presence.
        """
        return self._child_names(presence)

    @property
    def families(self) -> Generator[Family, None, None]:
        """
        Gets a generator that yields `Family` objects from all functional units.

        Yields:
            Generator[Family, None, None]: A generator that yields `Family` objects.
        """
        for func_unit in self.func_units:
            yield from func_unit.families

    @property
    def size(self) -> Tuple[int, int]:
        """
        Returns the size of the current object in terms of functional units and families.

        Returns:
            Tuple[int, int]: A tuple where the first element is the number of
            functional units, and the second element is the number of families.
        """
        return len(list(self.func_units)), len(list(self.families))

    def duplicate_fu(self, filter_type: str = None):
        """
        Duplicates items based on the specified filter type.

        Args:
            filter_type (str, optional): The type of filter to apply for duplicates. If `None`,
                no specific filter is applied.

        Yields:
            Any: The duplicated items obtained from the filtering process.
        """
        yield from self._duplicate(filter_type)

    def add(self, func_unit: FuncUnit):
        """
        Adds a functional unit to one of the sets in the instance based on its `presence`
        attribute.

        Args:
            func_unit (FuncUnit): The functional to be added.
        """
        func_unit.check_func_unit()
        super().add(func_unit)

    def get(self, name: str) -> Union[FuncUnit]:
        """
        Retrieves a function unit by its name.

        Args:
            name (str): The name of the function unit to be retrieved.

        Returns:
            Union[FuncUnit]: The function unit matching the specified name, or
            any applicable type derived from FuncUnit.
        """
        return super().get(name)

    def check_model(self):
        """
        Validates the consistency of a model by invoking an internal check method.

        Raises:
            Exception: Indicates a failure in model consistency, including the name of
                       the model and the specific error message.
        """
        try:
            self._check()
        except Exception as err:
            raise Exception(f"Consistency not respected in {self.name}. {err}")

    def read(self, data_model: dict):
        """
        Reads a data model dictionary and initializes or updates the object's attributes.

        Args:
            data_model (dict): Dictionary containing the data model information. Must include the mandatory
                keys 'name', 'parameters', and 'func_units'. Optionally, it may include 'window' and
                'canonical' attributes.
        """
        mandatory_keys = {"name", "parameters", "func_units"}
        param_mandatory = {"transitivity", "min_mandatory", "min_total"}

        check_dict(
            data_model, mandatory_keys=mandatory_keys, param_keys=param_mandatory
        )

        self.name = data_model["name"]
        self.read_parameters(data_model["parameters"], param_keys=param_mandatory)
        self.window = (
            data_model["window"] if "window" in data_model else self.transitivity + 1
        )
        for dict_fu in data_model["func_units"]:
            func_unit = FuncUnit()
            func_unit.model = self
            func_unit.read(dict_fu)
            self.add(func_unit)
        if "canonical" in data_model and data_model["canonical"] is not None:
            self.canonical = data_model["canonical"]

    @staticmethod
    def read_model(data_model: dict) -> Model:
        """
        Reads a data model dictionary and initializes a Model object with it.

        Args:
            data_model (dict): A dictionary containing the data model to be read.

        Returns:
            Model: An instance of the Model class populated with the provided data.
        """
        model = Model()
        model.read(data_model)
        return model


class FuncUnit(_BasicFeatures, _FuFamFeatures, _ModFuFeatures):
    """
    Represents a Functional Unit class that models functional and operational parameters, structural constraints, and various
    functional units and families. This class provides utilities to manage functional units, their relationships, and interactions.

    This class is designed to encapsulate features essential for analyzing and validating relationships between
    functional sub-elements, families, and systems. It defines the functional unit's behaviors, attributes, and operations.

    Attributes:
        name (str): Name of the functional unit.
        presence (str): Defines the presence rule (mandatory, accessory, forbidden, or neutral).
        mandatory (Set[FuncUnit, Family]): Set of mandatory sub-elements.
        accessory (Set[FuncUnit, Family]): Set of accessory sub-elements.
        forbidden (Set[FuncUnit, Family]): Set of forbidden sub-elements.
        neutral (Set[FuncUnit, Family]): Set of neutral sub-elements.
        min_mandatory (int): Minimum number of mandatory sub-elements.
        min_total (int): Minimum number of total sub-elements.
        same_strand (bool): Specifies if sub-elements must be on the same strand.
        transitivity (int): Size of the transitive closure used to build the graph.
        window (int): Number of neighboring genes considered in genomic searches.
        duplicate (int): Number of duplicates allowed.
        model (Model): Model in which the functional unit is defined.
        exchangeable (Set[str]): List of exchangeable families.
        multi_system (bool): Flag indicating if the unit can span multiple systems.
        multi_model (bool): Flag indicating if the unit can span multiple models.
    """

    def __init__(
        self,
        name: str = "",
        presence: str = "",
        mandatory: Set[FuncUnit, Family] = None,
        accessory: Set[FuncUnit, Family] = None,
        forbidden: Set[FuncUnit, Family] = None,
        neutral: Set[FuncUnit, Family] = None,
        min_mandatory: int = 1,
        min_total: int = 1,
        same_strand: bool = False,
        transitivity: int = 0,
        window: int = 1,
        duplicate: int = 0,
        model: Model = None,
        exchangeable: Set[str] = None,
        multi_system: bool = False,
        multi_model: bool = False,
    ):
        """
        Initializes an instance of a class with defined attributes, including
        functional, operational parameters, and structural constraints. This
        constructor manages various functionalities and initializes a robust
        system by leveraging functional units, families, and other attributes.

        Args:
            name (str): Name of the functional unit.
            presence (str): Defines the presence rule (mandatory, accessory, forbidden, or neutral).
            mandatory (Set[FuncUnit, Family]): Set of mandatory sub-elements.
            accessory (Set[FuncUnit, Family]): Set of accessory sub-elements.
            forbidden (Set[FuncUnit, Family]): Set of forbidden sub-elements.
            neutral (Set[FuncUnit, Family]): Set of neutral sub-elements.
            min_mandatory (int): Minimum number of mandatory sub-elements.
            min_total (int): Minimum number of total sub-elements.
            same_strand (bool): Specifies if sub-elements must be on the same strand.
            transitivity (int): Size of the transitive closure used to build the graph.
            window (int): Number of neighboring genes considered in genomic searches.
            duplicate (int): Number of duplicates allowed.
            model (Model): Model in which the functional unit is defined.
            exchangeable (Set[str]): List of exchangeable families.
            multi_system (bool): Flag indicating if the unit can span multiple systems.
            multi_model (bool): Flag indicating if the unit can span multiple models.
        """
        super().__init__(name=name, transitivity=transitivity, window=window)
        super(_BasicFeatures, self).__init__(
            presence=presence,
            duplicate=duplicate,
            parent=model,
            exchangeable=exchangeable,
            multi_system=multi_system,
            multi_model=multi_model,
        )
        super(_FuFamFeatures, self).__init__(
            mandatory=mandatory,
            min_mandatory=min_mandatory,
            accessory=accessory,
            neutral=neutral,
            min_total=min_total,
            forbidden=forbidden,
            same_strand=same_strand,
        )

    @property
    def model(self) -> Model:
        """
        Retrieves the parent model associated with this instance.

        Returns:
            Model: The parent model object.
        """
        return self._parent

    @model.setter
    def model(self, model: Model):
        """
        Sets the model attribute for the instance.

        Args:
            model (Model): The model instance to be set as the parent object.
        """
        self._parent = model

    @model.deleter
    def model(self):
        """
        Deletes the model property of the object.

        Raises:
            AttributeError: If the `_parent` attribute is not set or does not exist.

        """
        self._parent = None

    @property
    def families(self) -> Generator[Family, None, None]:
        """
        Yields Family objects contained in the current object's children.

        Yields:
            Generator[Family, None, None]: A generator producing `Family` objects.
        """
        yield from self._children

    def add(self, family: Family):
        """
        Adds a family to one of the sets in the instance based on its `presence`
        attribute.

        Args:
            family (Family): The functional to be added.
        """
        super().add(family)

    def get(self, name: str) -> Union[Family]:
        """
        Retrieves an instance of the 'Family' class using the provided name.

        Args:
            name (str): The name of the 'Family' instance to retrieve.

        Returns:
            Union[Family]: The retrieved 'Family' instance.
        """
        return super().get(name)

    def families_names(self, presence: str = None):
        """
        Returns a filtered list of family names based on the provided presence criteria.

        Args:
            presence (str, optional): Filter criteria for the presence of family members. Defaults to None.

        Returns:
            list: A list of names matching the given presence criteria.
        """
        return self._child_names(presence)

    @property
    def size(self) -> int:
        """
        Gets the size of the collection of families.

        The size property calculates and returns the total number of families
        present in the collection.

        Returns:
            int: The total number of families in the collection.
        """
        return len(list(self.families))

    def duplicate_fam(self, filter_type: str = None):
        """
        Generates items based on the specified filter type, if provided. This method uses an internal
        mechanism to yield results that match the filtering criteria.

        Args:
            filter_type (str, optional): The type of filter to apply. If None, no filtering
                is applied and all items are yielded.

        Yields:
            Any: Items matching the filter criteria.
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
            raise Exception(
                f"Consistency not respected in model {self.model.name} at functional unit {self.name}"
            )

    def read(self, data_fu: dict):
        """
        Read functional unit.

        Args:
            data_fu (dict): Data JSON file of all functional units.
        """
        mandatory_keys = {"name", "families", "presence"}
        fu_params = {
            "duplicate",
            "min_total",
            "min_mandatory",
            "transitivity",
            "multi_system",
            "multi_model",
        }
        check_dict(data_fu, mandatory_keys=mandatory_keys)

        self.name = data_fu["name"]
        self.presence = data_fu["presence"]
        self.read_parameters(
            data_fu["parameters"] if "parameters" in data_fu else {},
            param_keys=fu_params,
        )
        self.window = (
            data_fu["parameters"]["window"] if "window" in data_fu["parameters"] else self.transitivity + 1
        )
        for fam_dict in data_fu["families"]:
            family = Family()
            family.func_unit = self
            family.read(fam_dict)
            self.add(family)


class Family(_BasicFeatures, _FuFamFeatures):
    """
    Represents a family entity with configurable attributes and methods to manage its properties.

    This class is designed to encapsulate the attributes and behavior of a "family" entity,
    providing methods for initialization, property management, and reading family data from
    structured inputs like JSON files. It supports various configurations related to transitivity,
    duplicate handling, and exchangeable properties while interacting with functional units.

    Attributes:
        name (str): The name associated with the family instance, representing its identifier or label.
        transitivity (int): Represents the transitivity value, indicating the structural property of the family.
        window (int): Defines the processing window or range, typically a limit or scope.
        presence (str): Indicates the presence condition or state of the family instance.
        func_unit (FuncUnit): Functional unit linked to the family, representing a unit of operation or behavior.
        duplicate (int): Parameter controlling the duplication count or flag for this family.
        exchangeable (Set[str]): Defines a set of exchangeable features or components related to the family.
        multi_system (bool): Specifies if the family supports operation across multiple systems.
        multi_model (bool): Specifies if the family is configured to handle multiple models.
    """

    def __init__(
        self,
        name: str = "",
        transitivity: int = 0,
        window: int = 1,
        presence: str = "",
        func_unit: FuncUnit = None,
        duplicate: int = 0,
        exchangeable: Set[str] = None,
        multi_system: bool = False,
        multi_model: bool = False,
    ):
        """
        Initializes an instance of the class with specified attributes.

        Args:
            name (str): The name associated with the instance.
            transitivity (int): Transitivity value, typically indicating the relationship type in a structure.
            window (int): Specifies the window parameter for processing, generally indicates a range or limit.
            presence (str): Describes the presence condition or state of the instance.
            func_unit (FuncUnit): Functional unit associated with the instance represents a unit of functionality.
            duplicate (int): Duplicate control parameter typically denotes count or flag for duplicates.
            exchangeable (Set[str]): Set of exchangeable elements or features related to instance behavior.
            multi_system (bool): Indicates if the instance supports multi-system functionality.
            multi_model (bool): Indicates if the instance supports multiple models or configurations.
        """
        super().__init__(name=name, transitivity=transitivity, window=window)
        super(_BasicFeatures, self).__init__(
            presence=presence,
            duplicate=duplicate,
            parent=func_unit,
            exchangeable=exchangeable,
            multi_system=multi_system,
            multi_model=multi_model,
        )

    @property
    def func_unit(self) -> FuncUnit:
        """
        Gets the parent `FuncUnit` associated with this instance.

        Returns:
            FuncUnit: The parent functional unit.
        """
        return self._parent

    @func_unit.setter
    def func_unit(self, model: FuncUnit):
        """
        Sets the FuncUnit model for this instance. This method associates the provided FuncUnit model with
        the object by assigning it to the internal parent attribute.

        Args:
            model (FuncUnit): The FuncUnit model to set for this instance.
        """
        self._parent = model

    @func_unit.deleter
    def func_unit(self):
        """
        Deletes the `func_unit` property and removes the reference to the `_parent` attribute.

        Raises:
            AttributeError: If the `_parent` attribute does not exist or cannot be deleted.
        """
        self._parent = None

    @property
    def model(self) -> Model:
        """
        Returns the model associated with the functional unit.

        Returns:
            Model: The model instance associated with the functional unit.
        """
        return self.func_unit.model

    def read(self, data_fam: dict):
        """
        Read family.

        Args:
            data_fam (dict): Data JSON file with families.
        """
        fam_param = {"transitivity", "duplicate", "multi_system", "multi_model"}

        check_dict(data_fam, mandatory_keys={"name", "presence"})
        self.name = data_fam["name"]
        self.presence = data_fam["presence"]
        self.read_parameters(
            data_fam["parameters"] if "parameters" in data_fam else {},
            param_keys=fam_param,
        )
        self.window = (
            data_fam["window"] if "window" in data_fam else self.transitivity + 1
        )
        if "exchangeable" in data_fam:
            self.exchangeable = set(data_fam["exchangeable"])
