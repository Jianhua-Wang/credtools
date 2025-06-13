"""Main web application for CREDTOOLS visualization."""

import os
from pathlib import Path
from typing import Optional

try:
    import dash
    import dash_bootstrap_components as dbc
    import pandas as pd
    from dash import html
except ImportError as e:
    raise ImportError(
        "Web dependencies not found. Please install with: "
        "pip install dash dash-bootstrap-components plotly"
    ) from e


def create_app(
    webdata_dir: str = "webdata",
    debug: bool = True,
    port: int = 8080,
    host: str = "0.0.0.0",
) -> dash.Dash:
    """Create and configure the Dash app.

    Args
    ----
    webdata_dir : str
        Directory containing the processed web data
    debug : bool
        Whether to run in debug mode
    port : int
        Port to run the server on
    host : str
        Host to bind the server to

    Returns
    -------
    dash.Dash
        Configured Dash app instance
    """
    # Initialize the Dash app with a modern theme
    app = dash.Dash(
        __name__,
        external_stylesheets=[dbc.themes.BOOTSTRAP],
        use_pages=True,
        suppress_callback_exceptions=True,
    )

    # Set the webdata directory for the app
    # Store in server config rather than as attribute
    app.server.config["WEBDATA_DIR"] = webdata_dir

    # Check if data exists
    summary_file = os.path.join(webdata_dir, "all_loci_info.txt")
    if not os.path.exists(summary_file):
        raise FileNotFoundError(
            f"Web data not found at {summary_file}. "
            "Please run data export first using the export functionality."
        )

    # Load the summary data to validate
    df_summary = pd.read_csv(summary_file, sep="\t")
    meta_methods = sorted(df_summary["meta_type"].unique())
    fine_mapping_methods = sorted(df_summary["tool"].unique())

    # Create the app layout
    app.layout = create_layout()

    return app


def create_layout():
    """Create the main layout for the web app."""
    return html.Div(
        [
            # Banner
            dbc.Navbar(
                dbc.Container(
                    [
                        dbc.Row(
                            [
                                dbc.Col(
                                    [
                                        html.H2(
                                            "Multi-ancestry GWAS Fine-Mapping",
                                            className="mb-0 text-white",
                                        ),
                                    ],
                                    width="auto",
                                ),
                            ],
                            align="center",
                        ),
                    ]
                ),
                color="primary",
                dark=True,
                className="mb-4",
            ),
            # Main content
            dbc.Container(
                [
                    # Page content will be inserted here
                    dash.page_container
                ]
            ),
        ],
        className="bg-light min-vh-100",
    )


def run_app(
    webdata_dir: str = "webdata",
    debug: bool = True,
    port: int = 8080,
    host: str = "0.0.0.0",
) -> None:
    """Run the web app.

    Args
    ----
    webdata_dir : str
        Directory containing the processed web data
    debug : bool
        Whether to run in debug mode
    port : int
        Port to run the server on
    host : str
        Host to bind the server to
    """
    app = create_app(webdata_dir=webdata_dir, debug=debug, port=port, host=host)

    # Import pages to register them
    try:
        from . import pages
    except ImportError:
        print("Warning: Web pages not found. Some functionality may be limited.")

    app.run_server(debug=debug, port=port, host=host)
