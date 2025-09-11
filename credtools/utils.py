"""Functions and decorators for common tasks in Python programming."""

import logging
import os
import shutil
import subprocess
import tempfile
from functools import wraps
from pathlib import Path
from typing import Any, Callable, Dict, List, Optional, TypeVar, Union

logger = logging.getLogger("Utils")


# Type variable for decorator
F = TypeVar("F", bound=Callable[..., Any])


def io_in_tempdir(dir: str = "./tmp") -> Callable[[F], F]:
    """
    Create a temporary directory for I/O operations during function execution.

    This decorator creates a temporary directory before executing the decorated function and
    provides the path to this directory via the `temp_dir` keyword argument. After the function
    execution, the temporary directory is removed based on the logging level:
    - If the logging level is set to `INFO` or higher, the temporary directory is deleted.
    - If the logging level is lower than `INFO` (e.g., `DEBUG`), the directory is retained for inspection.

    Parameters
    ----------
    dir : str, optional
        The parent directory where the temporary directory will be created, by default "./tmp".

    Returns
    -------
    Callable[[F], F]
        A decorator that manages a temporary directory for the decorated function.

    Raises
    ------
    OSError
        If the temporary directory cannot be created.

    Examples
    --------
    ```python
    @io_in_tempdir(dir="./temporary")
    def process_data(temp_dir: str, data: str) -> None:
        # Perform I/O operations using temp_dir
        with open(f"{temp_dir}/data.txt", "w") as file:
            file.write(data)

    process_data(data="Sample data")
    ```
    """

    def decorator(func: F) -> F:
        @wraps(func)
        def wrapper(*args: Any, **kwargs: Any) -> Any:
            if not os.path.exists(dir):
                os.makedirs(dir, exist_ok=True)
            temp_dir = tempfile.mkdtemp(dir=dir)
            logger.debug(f"Created temporary directory: {temp_dir}")

            try:
                # Inject temp_dir into the function's keyword arguments
                result = func(*args, temp_dir=temp_dir, **kwargs)
            except Exception as e:
                logger.error(f"An error occurred in function '{func.__name__}': {e}")
                raise
            else:
                # Determine whether to remove the temporary directory based on the logging level
                if logger.getEffectiveLevel() >= logging.INFO:
                    try:
                        shutil.rmtree(temp_dir)
                        logger.debug(f"Removed temporary directory: {temp_dir}")
                    except Exception as cleanup_error:
                        logger.warning(
                            f"Failed to remove temporary directory '{temp_dir}': {cleanup_error}"
                        )
                else:
                    logger.debug(
                        f"Retaining temporary directory '{temp_dir}' for inspection due to logging level."
                    )
                return result

        return wrapper  # type: ignore

    return decorator


class ExternalTool:
    """
    A class to manage and run external tools.

    This class provides a unified interface for managing external bioinformatics tools,
    handling path resolution, and executing commands with proper error checking.

    Parameters
    ----------
    name : str
        The name of the external tool.
    default_path : Optional[str], optional
        The default path to the tool if not found in the system PATH, by default None.

    Attributes
    ----------
    name : str
        The name of the external tool.
    default_path : Optional[str]
        The default path to the tool if not found in the system PATH.
    custom_path : Optional[str]
        A custom path set by the user.

    Methods
    -------
    set_custom_path(path: str) -> None
        Sets a custom path for the tool if it exists.
    get_path() -> str
        Retrieves the path to the tool, checking custom, system, and default paths.
    run(command: List[str], log_file: str, output_file_path: Optional[Union[str, List[str]]]) -> None
        Runs the tool with the given arguments.

    Examples
    --------
    >>> tool = ExternalTool("samtools", "/usr/local/bin/samtools")
    >>> tool.set_custom_path("/opt/samtools/bin/samtools")
    >>> tool.run(["view", "-h", "input.bam"], "samtools.log", "output.sam")
    """

    def __init__(self, name: str, default_path: Optional[str] = None) -> None:
        """
        Initialize the ExternalTool with a name and an optional default path.

        Parameters
        ----------
        name : str
            The name of the external tool.
        default_path : Optional[str], optional
            The default path to the tool if not found in the system PATH, by default None.
        """
        self.name = name
        self.default_path = default_path
        self.custom_path: Optional[str] = None

    def set_custom_path(self, path: str) -> None:
        """
        Set a custom path for the tool if it exists.

        Parameters
        ----------
        path : str
            The custom path to set.

        Raises
        ------
        FileNotFoundError
            If the custom path does not exist.
        """
        if os.path.exists(path):
            self.custom_path = path
        else:
            raise FileNotFoundError(
                f"Custom path for {self.name} does not exist: {path}"
            )

    def get_path(self) -> str:
        """
        Retrieve the path to the tool, checking custom, system, and default paths.

        The function checks paths in the following order:
        1. Custom path (if set via set_custom_path)
        2. System PATH (using shutil.which)
        3. Default path (relative to package directory)

        Returns
        -------
        str
            The path to the tool.

        Raises
        ------
        FileNotFoundError
            If the tool cannot be found in any of the paths.
        """
        if self.custom_path:
            return self.custom_path

        system_tool = shutil.which(self.name)
        if system_tool:
            return system_tool

        if self.default_path:
            package_dir = Path(__file__).parent
            internal_tool = package_dir / self.default_path
            if internal_tool.exists():
                return str(internal_tool)

        raise FileNotFoundError(f"Could not find {self.name} executable")

    def run(
        self,
        command: List[str],
        log_file: str,
        output_file_path: Optional[Union[str, List[str]]] = None,
    ) -> None:
        """
        Execute a command line instruction, log the output, and handle errors.

        This function runs the given command, captures stdout and stderr,
        logs them using logging.debug, and raises exceptions for command failures
        or missing output files.

        Parameters
        ----------
        command : List[str]
            The command line instruction to be executed (without the tool name).
        log_file : str
            The file to log the output.
        output_file_path : Optional[Union[str, List[str]]], optional
            The expected output file path(s). If provided, the function will check
            if these files exist after command execution, by default None.

        Raises
        ------
        subprocess.CalledProcessError
            If the command execution fails.
        FileNotFoundError
            If the specified output file is not found after command execution.

        Examples
        --------
        >>> tool = ExternalTool("finemap")
        >>> tool.run(["--help"], "finemap.log")
        >>> tool.run(["--in-files", "data.master"], "finemap.log", "output.snp")
        """
        full_command = [self.get_path()] + command
        try:
            # Run the command and capture output
            logger.debug(f"Run command: {' '.join(full_command)}")
            with open(log_file, "w") as log:
                subprocess.run(
                    full_command, shell=False, check=True, stdout=log, stderr=log
                )

            # Check for output file if path is provided
            if output_file_path:
                if isinstance(output_file_path, str):
                    output_file_path = [output_file_path]
                for path in output_file_path:
                    if not os.path.exists(path):
                        raise FileNotFoundError(
                            f"Expected output file not found: {path}"
                        )

        except Exception as e:
            logger.error(f"Command execution failed: {e}\nSee {log_file} for details.")
            raise


class ToolManager:
    """
    A class to manage multiple external tools.

    This class provides a centralized registry for managing multiple external tools,
    allowing for easy registration, configuration, and execution of bioinformatics software.

    Attributes
    ----------
    tools : Dict[str, ExternalTool]
        A dictionary to store registered tools by their names.

    Methods
    -------
    register_tool(name: str, default_path: Optional[str] = None) -> None
        Registers a new tool with an optional default path.
    set_tool_path(name: str, path: str) -> None
        Sets a custom path for a registered tool.
    get_tool(name: str) -> ExternalTool
        Retrieves a registered tool by its name.
    run_tool(name: str, args: List[str], log_file: str, output_file_path: Optional[Union[str, List[str]]]) -> None
        Runs a registered tool with the given arguments.

    Examples
    --------
    >>> manager = ToolManager()
    >>> manager.register_tool("finemap", "bin/finemap")
    >>> manager.set_tool_path("finemap", "/usr/local/bin/finemap")
    >>> manager.run_tool("finemap", ["--help"], "finemap.log")
    """

    def __init__(self) -> None:
        """Initialize the ToolManager with an empty dictionary of tools."""
        self.tools: Dict[str, ExternalTool] = {}

    def register_tool(self, name: str, default_path: Optional[str] = None) -> None:
        """
        Register a new tool with an optional default path.

        Parameters
        ----------
        name : str
            The name of the tool to register.
        default_path : Optional[str], optional
            The default path to the tool if not found in the system PATH, by default None.
        """
        self.tools[name] = ExternalTool(name, default_path)

    def set_tool_path(self, name: str, path: str) -> None:
        """
        Set a custom path for a registered tool.

        Parameters
        ----------
        name : str
            The name of the registered tool.
        path : str
            The custom path to set for the tool.

        Raises
        ------
        KeyError
            If the tool is not registered.
        """
        if name not in self.tools:
            raise KeyError(f"Tool {name} is not registered")
        self.tools[name].set_custom_path(path)

    def get_tool(self, name: str) -> ExternalTool:
        """
        Retrieve a registered tool by its name.

        Parameters
        ----------
        name : str
            The name of the registered tool.

        Returns
        -------
        ExternalTool
            The registered tool.

        Raises
        ------
        KeyError
            If the tool is not registered.
        """
        if name not in self.tools:
            raise KeyError(f"Tool {name} is not registered")
        return self.tools[name]

    def run_tool(
        self,
        name: str,
        args: List[str],
        log_file: str,
        output_file_path: Optional[Union[str, List[str]]] = None,
    ) -> None:
        """
        Run a registered tool with the given arguments.

        Parameters
        ----------
        name : str
            The name of the registered tool.
        args : List[str]
            The arguments to pass to the tool.
        log_file : str
            The file to log the output.
        output_file_path : Optional[Union[str, List[str]]], optional
            The expected output file path(s). If provided, the function will check
            if these files exist after command execution, by default None.

        Raises
        ------
        KeyError
            If the tool is not registered.
        subprocess.CalledProcessError
            If the subprocess call fails.
        FileNotFoundError
            If expected output files are not found after execution.
        """
        if name not in self.tools:
            raise KeyError(f"Tool {name} is not registered")
        return self.get_tool(name).run(args, log_file, output_file_path)


tool_manager = ToolManager()
for tool in ["finemap", "SuSiEx"]:
    tool_manager.register_tool(tool, f"bin/{tool}")
