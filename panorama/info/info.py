#!/usr/bin/env python3
# coding:utf-8

"""
Panorama Info Module - Pangenome Information Extraction and Visualization

This module provides functionality to extract, process, and export information
from pangenome HDF5 files, generating interactive HTML reports for analysis.
"""

# Standard library imports
import argparse
import logging
from collections import defaultdict
from copy import copy
from math import ceil, floor
from pathlib import Path
from typing import Any, Dict, List, Union

# Third-party imports
import pandas as pd
import tables
from bokeh.embed import file_html
from bokeh.io import curdoc
from bokeh.layouts import column, row
from bokeh.models import (
    Button,
    CheckboxGroup,
    ColumnDataSource,
    DataTable,
    Div,
    RadioButtonGroup,
    RangeSlider,
    TableColumn,
)
from bokeh.models.callbacks import CustomJS
from ppanggolin.formats import get_pangenome_parameters, read_info
from ppanggolin.info.info import read_metadata_status, read_status
from tqdm import tqdm

# Local imports
from panorama.utils import check_tsv_sanity

# Type aliases for better readability
PangenomeInfo = Dict[str, Union[int, str]]
PangenomePaths = Dict[str, PangenomeInfo]
InfoDict = Dict[str, Dict[str, Any]]
ContentDict = Dict[str, Dict[str, Union[int, float, Dict[str, Union[int, float]]]]]


class PangenomeInfoExtractor:
    """
    A class to extract and process information from pangenome HDF5 files.

    This class handles the extraction of status, content, parameters, and metadata
    information from multiple pangenome files and provides methods to export
    this information as interactive HTML reports.
    """

    def __init__(self, disable_bar: bool = False):
        """
        Initialize the PangenomeInfoExtractor.

        Args:
            disable_bar (bool): Whether to disable the progress bar during processing.
        """
        self.disable_bar = disable_bar
        self.logger = logging.getLogger("PANORAMA")

    def extract_info(
        self,
        pangenomes_path: PangenomePaths,
        status: bool = False,
        content: bool = False,
        parameters: bool = False,
        metadata: bool = False,
    ) -> InfoDict:
        """
        Extract information from multiple pangenome files.

        Args:
            pangenomes_path (PangenomePaths): Dictionary mapping pangenome names to their file paths.
            status (bool, optional): Whether to extract status information. Defaults to False.
            content (bool, optional): Whether to extract content information. Defaults to False.
            parameters (bool, optional): Whether to extract parameter information. Defaults to False.
            metadata (bool, optional): Whether to extract metadata information. Defaults to False.

        Returns:
            InfoDict: Dictionary containing all extracted information organized by type.

        Note:
            If no specific information type is requested, all types will be extracted.
        """
        # If no specific info type is requested, extract all
        if not any([status, content, parameters, metadata]):
            status = content = parameters = metadata = True

        info_dict = defaultdict(dict)

        self.logger.info(f"Processing {len(pangenomes_path)} pangenome files...")

        for pangenome_name, pan_info in tqdm(
            pangenomes_path.items(),
            unit="Pangenome",
            disable=self.disable_bar,
            desc="Extracting pangenome information",
        ):
            try:
                with tables.open_file(pan_info["path"], "r") as h5f:
                    if status:
                        info_dict["status"][pangenome_name] = self._extract_status(h5f)
                    if content:
                        info_dict["content"][pangenome_name] = self._extract_content(
                            h5f
                        )
                    if parameters:
                        info_dict["parameters"][pangenome_name] = (
                            self._extract_parameters(h5f)
                        )
                    if metadata:
                        info_dict["metadata"][pangenome_name] = self._extract_metadata(
                            h5f
                        )

            except Exception as e:
                self.logger.error(f"Error processing {pangenome_name}: {str(e)}")
                raise

        return dict(info_dict)  # Convert defaultdict to regular dict

    @staticmethod
    def _extract_status(h5f: tables.File) -> Dict[str, Union[bool, str]]:
        """
        Extract status information from an HDF5 file.

        Args:
            h5f (tables.File): Open HDF5 file handle.

        Returns:
            Dict[str, Union[bool, str]]: Status information.
        """
        return read_status(h5f)["Status"]

    @staticmethod
    def _extract_content(h5f: tables.File) -> Dict[str, Any]:
        """
        Extract content information from an HDF5 file.

        Args:
            h5f (tables.File): Open HDF5 file handle.

        Returns:
            Dict[str, Any]: Content information.
        """
        return read_info(h5f)["Content"]

    @staticmethod
    def _extract_parameters(h5f: tables.File) -> Dict[str, Any]:
        """
        Extract parameter information from an HDF5 file.

        Args:
            h5f (tables.File): Open HDF5 file handle.

        Returns:
            Dict[str, Any]: Parameter information.

        Note:
            This method currently prints parameters for debugging but doesn't return them.
            Implementation needs to be completed.
        """
        step_to_parameters = get_pangenome_parameters(h5f)
        # TODO: Implement proper parameter extraction and return
        return {}

    @staticmethod
    def _extract_metadata(h5f: tables.File) -> Dict[str, Any]:
        """
        Extract metadata information from an HDF5 file.

        Args:
            h5f (tables.File): Open HDF5 file handle.

        Returns:
            Dict[str, Any]: Metadata information.
        """
        return read_metadata_status(h5f)


class HTMLExporter:
    """
    A class to export pangenome information as interactive HTML reports.

    This class provides methods to create interactive HTML tables and visualizations
    for different types of pangenome information.
    """

    def __init__(self, output_dir: Path):
        """
        Initialize the HTMLExporter.

        Args:
            output_dir (Path): Directory where HTML files will be saved.
        """
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.logger = logging.getLogger("PANORAMA")

        # Ensure JavaScript files exist
        self._validate_js_files()

    def _validate_js_files(self) -> None:
        """
        Validate that required JavaScript files exist.

        Raises:
            FileNotFoundError: If required JavaScript files are not found.
        """
        js_files = ["filterData.js", "download.js", "hide_show.js"]
        js_dir = Path(__file__).parent

        for js_file in js_files:
            if not (js_dir / js_file).exists():
                self.logger.warning(f"JavaScript file {js_file} not found in {js_dir}")

    def export_status(self, status_dict: Dict[str, Union[bool, str]]) -> None:
        """
        Export status information to an interactive HTML file.

        Args:
            status_dict (Dict[str, Union[bool, str]]): Dictionary containing status information
                for each pangenome, where keys are pangenome names and values are status data.

        Note:
            Creates an HTML file with filtering capabilities for boolean columns and
            download functionality.
        """
        self.logger.debug("Exporting status information")

        # Convert dictionary to DataFrame for easier manipulation
        status_df = pd.DataFrame.from_dict(status_dict, orient="index")
        status_df = status_df.reset_index().rename(columns={"index": "Pangenome name"})

        # Create a Bokeh data source
        source = ColumnDataSource(status_df)

        # Create table columns
        columns = [
            TableColumn(field=col_name, title=col_name)
            for col_name in status_df.columns
        ]

        # Create a data table with appropriate sizing
        table_height = min(30 * status_df.shape[0], 800)  # Cap height for usability
        status_table = DataTable(
            source=source,
            columns=columns,
            index_position=None,
            width=1900,
            height=table_height,
        )

        # Create radio buttons for boolean column filtering
        radio_buttons = self._create_boolean_filters(source)

        # Create control buttons
        filter_button = self._create_filter_button(source, radio_buttons)
        download_button = self._create_download_button(source, "pangenomes_status.tsv")

        # Layout the components
        layout = column(
            status_table,
            row(
                Div(text="", width=160),
                row(*radio_buttons, spacing=29),
                Div(text="", width=20),
                filter_button,
                download_button,
            ),
        )

        # Export to HTML
        self._save_html(layout, "status_info.html", "Pangenome Status Information")

    def export_content(self, content_dict: ContentDict) -> None:
        """
        Export content information to an interactive HTML file.

        Args:
            content_dict (ContentDict): Dictionary containing content information for each
                pangenome, including statistics about gene families, modules, and other metrics.

        Note:
            Creates an HTML file with column visibility controls, range sliders for filtering,
            and download functionality.
        """
        self.logger.debug("Exporting content information")

        # Unpack nested dictionaries and create DataFrame
        unpacked_content = self._unpack_content_dict(content_dict)
        content_df = pd.DataFrame.from_dict(unpacked_content, orient="index")
        content_df = content_df.reset_index().rename(
            columns={"index": "Pangenome name"}
        )

        # Rename columns for better readability
        content_df = self._rename_content_columns(content_df)

        # Reorder columns for logical grouping
        content_df = self._reorder_content_columns(content_df)

        # Fill NaN values with 0
        content_df = content_df.fillna(value=0)

        # Create a Bokeh data source
        source = ColumnDataSource(content_df)

        # Define which columns should be visible by default
        visible_indices = [0, 1, 2, 7, 8, 9, 10, 11, 16, 17, 18]

        # Create table columns with visibility settings
        columns = [
            TableColumn(
                field=col_name, title=col_name, visible=(idx in visible_indices)
            )
            for idx, col_name in enumerate(content_df.columns)
        ]

        # Create a data table
        table_height = min(32 * content_df.shape[0], 800)
        content_table = DataTable(
            source=source,
            columns=columns,
            index_position=None,
            width=1900,
            height=table_height,
        )

        # Create controls
        checkbox_group = self._create_column_visibility_control(
            source, columns, visible_indices
        )
        sliders = self._create_range_sliders(source)
        download_button = self._create_download_button(source, "pangenomes_content.tsv")

        # Layout components
        layout = self._layout_content_components(
            content_table, checkbox_group, sliders, download_button
        )

        # Export to HTML
        self._save_html(layout, "content_info.html", "Pangenome Content Information")

    @staticmethod
    def _unpack_content_dict(content_dict: ContentDict) -> Dict[str, Dict[str, Any]]:
        """
        Unpack nested dictionaries in the content dictionary for flat DataFrame creation.

        Args:
            content_dict (ContentDict): Dictionary with potentially nested content information.

        Returns:
            Dict[str, Dict[str, Any]]: Flattened dictionary suitable for DataFrame creation.

        Note:
            This method flattens nested dictionaries using dot notation for keys and
            handles special cases for module counts and other metrics.
        """
        unpacked_dict = defaultdict(dict)

        for pangenome_name, content in content_dict.items():
            for key, value in content.items():
                if isinstance(value, dict):
                    # Handle nested dictionaries
                    for nested_key, nested_value in value.items():
                        if isinstance(nested_value, dict):
                            # Double-nested: use dot notation
                            for deep_key, deep_value in nested_value.items():
                                unpacked_dict[pangenome_name][
                                    f"{nested_key}.{deep_key}"
                                ] = deep_value
                        elif nested_key == "Number_of_modules" or nested_key.endswith(
                            "_count"
                        ):
                            # Special handling for count fields
                            unpacked_dict[pangenome_name][key] = nested_value
                        else:
                            unpacked_dict[pangenome_name][nested_key] = nested_value
                else:
                    # Simple key-value pair
                    unpacked_dict[pangenome_name][key] = value

        return dict(unpacked_dict)

    @staticmethod
    def _rename_content_columns(df: pd.DataFrame) -> pd.DataFrame:
        """
        Rename DataFrame columns for better readability.

        Args:
            df (pd.DataFrame): DataFrame with content information.

        Returns:
            pd.DataFrame: DataFrame with renamed columns.
        """
        column_renames = {
            "min_genomes_frequency": "Genomes_frequency.min",
            "max_genomes_frequency": "Genomes_frequency.max",
            "mean_genomes_frequency": "Genomes_frequency.mean",
            "sd_genomes_frequency": "Genomes_frequency.sd",
        }
        return df.rename(columns=column_renames)

    @staticmethod
    def _reorder_content_columns(df: pd.DataFrame) -> pd.DataFrame:
        """
        Reorder DataFrame columns for logical grouping and presentation.

        Args:
            df (pd.DataFrame): DataFrame with content information.

        Returns:
            pd.DataFrame: DataFrame with reordered columns.
        """
        columns = df.columns.values.tolist()

        # Define column order for logical grouping
        # This groups related metrics together for easier analysis
        column_order = (
            columns[:3]  # Basic info columns
            + columns[6:8]  # Frequency stats
            + columns[9:7:-1]  # Reversed frequency range
            + columns[3:6]  # Gene family metrics
            + columns[10:13]  # Module metrics
            + columns[17:20]  # Additional statistics
            + columns[13:17]  # Core/accessory metrics
            + columns[20:]  # Remaining columns
        )

        return df.reindex(columns=column_order)

    @staticmethod
    def _create_boolean_filters(source: ColumnDataSource) -> List[RadioButtonGroup]:
        """
        Create radio button groups for filtering boolean columns.

        Args:
            source (ColumnDataSource): Bokeh data source.

        Returns:
            List[RadioButtonGroup]: List of radio button groups for boolean columns.
        """
        radio_buttons = []

        for column_name in source.column_names:
            if (
                hasattr(source.data[column_name], "dtype")
                and source.data[column_name].dtype == "bool"
            ):
                radio_buttons.append(
                    RadioButtonGroup(
                        labels=["all", "true", "false"], active=0, name=column_name
                    )
                )

        return radio_buttons

    def _create_filter_button(
        self, source: ColumnDataSource, radio_buttons: List[RadioButtonGroup]
    ) -> Button:
        """
        Create a filter button with a JavaScript callback for boolean filtering.

        Args:
            source (ColumnDataSource): Bokeh data source.
            radio_buttons (List[RadioButtonGroup]): Radio buttons for filtering.

        Returns:
            Button: Configured filter button.
        """
        filter_button = Button(label="Filter", button_type="primary")

        js_code = (
            self._load_js_file("filterData.js")
            + "filterDataBool(source, radio_buttons);"
        )
        filter_button.js_on_event(
            "button_click",
            CustomJS(
                args={
                    "source": source,
                    "radio_buttons": radio_buttons,
                    "original_source": copy(source.data),
                },
                code=js_code,
            ),
        )

        return filter_button

    def _create_download_button(
        self, source: ColumnDataSource, filename: str
    ) -> Button:
        """
        Create a download button with a JavaScript callback for data export.

        Args:
            source (ColumnDataSource): Bokeh data source.
            filename (str): Name of the file to download.

        Returns:
            Button: Configured download button.
        """
        download_button = Button(label="Download", button_type="success")

        download_button.js_on_event(
            "button_click",
            CustomJS(
                args={"source": source, "filename": filename},
                code=self._load_js_file("download.js"),
            ),
        )

        return download_button

    def _create_column_visibility_control(
        self,
        source: ColumnDataSource,
        columns: List[TableColumn],
        visible_indices: List[int],
    ) -> CheckboxGroup:
        """
        Create a checkbox group for controlling column visibility.

        Args:
            source (ColumnDataSource): Bokeh data source.
            columns (List[TableColumn]): Table columns.
            visible_indices (List[int]): Indices of initially visible columns.

        Returns:
            CheckboxGroup: Configured a checkbox group.
        """
        checkbox_group = CheckboxGroup(
            labels=source.column_names[1:],  # Exclude first column (pangenome name)
            active=visible_indices,
        )

        checkbox_group.js_on_change(
            "active",
            CustomJS(
                args={
                    "source": source,
                    "columns": columns,
                    "checkbox_group": checkbox_group,
                },
                code=self._load_js_file("hide_show.js"),
            ),
        )

        return checkbox_group

    def _create_range_sliders(self, source: ColumnDataSource) -> List[RangeSlider]:
        """
        Create range sliders for numeric column filtering.

        Args:
            source (ColumnDataSource): Bokeh data source.

        Returns:
            List[RangeSlider]: List of configured range sliders.
        """
        sliders = []

        for column_name in source.column_names[1:]:  # Skip pangenome name column
            column_data = source.data[column_name]

            if hasattr(column_data, "dtype") and column_data.dtype in (int, float):
                # Calculate slider parameters based on data type
                if column_data.dtype == int:
                    start = min(column_data)
                    end = max(column_data)
                    step = 1
                else:  # float
                    start = floor(min(column_data))
                    end = ceil(max(column_data))
                    step = 0.1

                # Ensure start != end for slider functionality
                if start == end:
                    if start > 1:
                        start -= 1
                    end += 1

                slider = RangeSlider(
                    start=start,
                    end=end,
                    value=(start, end),
                    step=step,
                    title=column_name,
                )

                sliders.append(slider)

        # Set up slider interactions
        self._setup_slider_callbacks(sliders, source)

        return sliders

    def _setup_slider_callbacks(
        self, sliders: List[RangeSlider], source: ColumnDataSource
    ) -> None:
        """
        Set up JavaScript callbacks for slider interactions.

        Args:
            sliders (List[RangeSlider]): List of range sliders.
            source (ColumnDataSource): Bokeh data source.
        """
        js_code_base = (
            self._load_js_file("filterData.js")
            + "filterDataSliders(source, slider, other_sliders);"
        )

        for idx, slider in enumerate(sliders):
            other_sliders = sliders[:idx] + sliders[idx + 1 :]

            slider.js_on_change(
                "value_throttled",
                CustomJS(
                    args={
                        "source": source,
                        "slider": slider,
                        "other_sliders": other_sliders,
                        "original_source": source.data.copy(),
                    },
                    code=js_code_base,
                ),
            )

    @staticmethod
    def _layout_content_components(
        table: DataTable,
        checkbox_group: CheckboxGroup,
        sliders: List[RangeSlider],
        download_button: Button,
    ) -> column:
        """
        Layout content export components in an organized manner.

        Args:
            table (DataTable): Main data table.
            checkbox_group (CheckboxGroup): Column visibility controls.
            sliders (List[RangeSlider]): Range sliders for filtering.
            download_button (Button): Download button.

        Returns:
            column: Bokeh layout containing all components.
        """
        # Group sliders into columns for better organization
        slider_columns = []
        sliders_per_column = 6

        for i in range(0, len(sliders), sliders_per_column):
            slider_group = sliders[i : i + sliders_per_column]
            if i + sliders_per_column >= len(sliders):
                # Add the download button to the last column
                slider_group.append(download_button)
            slider_columns.append(column(*slider_group))

        controls_row = row(checkbox_group, *slider_columns, spacing=20)

        return column(table, controls_row, spacing=50)

    def _load_js_file(self, filename: str) -> str:
        """
        Load JavaScript code from the file.

        Args:
            filename (str): Name of the JavaScript file.

        Returns:
            str: JavaScript code content.
        """
        js_path = Path(__file__).parent / filename
        try:
            return js_path.read_text()
        except FileNotFoundError:
            self.logger.warning(
                f"JavaScript file {filename} not found, using empty string"
            )
            return ""

    def _save_html(self, layout: column, filename: str, title: str) -> None:
        """
        Save Bokeh layout to the HTML file.

        Args:
            layout (column): Bokeh layout to save.
            filename (str): Output filename.
            title (str): HTML page title.
        """
        curdoc().add_root(layout)

        html_content = file_html(layout, "cdn", title=title)
        output_path = self.output_dir / filename

        with open(output_path, "w", encoding="utf-8") as file:
            file.write(html_content)

        self.logger.info(f"Exported {title} to {output_path}")


def export_info(info_dict: InfoDict, output: Path) -> None:
    """
    Export extracted information to HTML files.

    Args:
        info_dict (InfoDict): Dictionary containing extracted pangenome information.
        output (Path): Output directory path where HTML files will be saved.

    Raises:
        NotImplementedError: If parameter or metadata export is requested (not yet implemented).
        KeyError: If an unrecognized information type is provided.

    Note:
        Currently supports exporting status and content information. Parameter and
        metadata export functionality is planned for future implementation.
    """
    exporter = HTMLExporter(output)

    for info_type, data in info_dict.items():
        if info_type == "status":
            exporter.export_status(data)
        elif info_type == "content":
            exporter.export_content(data)
        elif info_type == "parameters":
            raise NotImplementedError(
                "Parameter export is not yet implemented. "
                "This feature is planned for a future release."
            )
        elif info_type == "metadata":
            raise NotImplementedError(
                "Metadata export is not yet implemented. "
                "This feature is planned for a future release."
            )
        else:
            raise KeyError(
                f"Unrecognized information type: '{info_type}'. "
                f"Supported types are: status, content, parameters, metadata"
            )


def launch(args: argparse.Namespace) -> None:
    """
    Main command launcher for the info extraction functionality.

    Args:
        args (argparse.Namespace): Parsed command-line arguments containing all user options.

    Note:
        This function orchestrates the entire workflow: validates input files,
        extracts information, and exports results to HTML files.
    """
    logger = logging.getLogger("PANORAMA")
    logger.debug("Launching info command")

    try:
        # Validate and parse input pangenome file list
        pangenomes_to_path = check_tsv_sanity(args.pangenomes)
        logger.info(f"Validated {len(pangenomes_to_path)} pangenome files")

        # Extract information from pangenome files
        extractor = PangenomeInfoExtractor(disable_bar=args.disable_prog_bar)
        info_dict = extractor.extract_info(
            pangenomes_path=pangenomes_to_path,
            status=args.status,
            content=args.content,
            parameters=args.parameters,
            metadata=args.metadata,
        )

        # Export information to HTML files
        export_info(info_dict, args.output)

        logger.info("Info extraction and export completed successfully")

    except Exception as e:
        logger.error(f"Error during info command execution: {str(e)}")
        raise


def subparser(sub_parser: argparse._SubParsersAction) -> argparse.ArgumentParser:
    """
    Create the subparser for the info command.

    Args:
        sub_parser (argparse._SubParsersAction): Parent subparser to add the info command to.

    Returns:
        argparse.ArgumentParser: Configured argument parser for the info command.
    """
    parser = sub_parser.add_parser(
        "info",
        help="Extract and export information from pangenome files",
        description="Extract status, content, parameters, and metadata information "
        "from pangenome HDF5 files and export as interactive HTML reports.",
    )
    parser_info(parser)
    return parser


def parser_info(parser: argparse.ArgumentParser) -> None:
    """
    Configure an argument parser with info command-specific options.

    Args:
        parser (argparse.ArgumentParser): Argument parser to configure.

    Note:
        Sets up required arguments (input files, output directory) and optional
        arguments for controlling which information types to extract and export.
    """
    # Required arguments
    required = parser.add_argument_group(
        title="Required arguments",
        description="All of the following arguments are required:",
    )
    required.add_argument(
        "-p",
        "--pangenomes",
        required=True,
        type=Path,
        help="Path to a TSV file listing pangenome .h5 files with their names and paths",
    )
    required.add_argument(
        "-o",
        "--output",
        required=True,
        type=Path,
        help="Output directory where HTML reports will be saved",
    )

    # Information display options
    display = parser.add_argument_group(
        title="Information Display Options",
        description="Select which types of information to extract (default: all types)",
    )
    display.add_argument(
        "-s",
        "--status",
        action="store_true",
        help="Extract and export status information showing completion status "
        "of different analysis steps for each pangenome",
    )
    display.add_argument(
        "-c",
        "--content",
        action="store_true",
        help="Extract and export content information including gene family "
        "statistics, core/accessory genome metrics, and module information",
    )
    display.add_argument(
        "-a",
        "--parameters",
        action="store_true",
        help="Extract and export parameters used at each step of pangenome "
        "generation (currently not implemented)",
    )
    display.add_argument(
        "-m",
        "--metadata",
        action="store_true",
        help="Extract and export metadata information stored in pangenome files "
        "(currently not implemented)",
    )
