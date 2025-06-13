"""Locus page for CREDTOOLS web visualization."""

import os
from pathlib import Path

try:
    import dash
    import dash_bootstrap_components as dbc
    import dash_mantine_components as dmc
    import numpy as np
    import pandas as pd
    import plotly.express as px
    import plotly.graph_objects as go
    from dash import Input, Output, State, callback, dash_table, dcc, html
except ImportError:
    # Web dependencies not available
    pass

# Register this page
try:
    dash.register_page(
        __name__, path="/locus/<locus_id>", path_template="/locus/<locus_id>"
    )
except (ImportError, AttributeError):
    pass


def get_r2_color(r2):
    """Get color for r2 value."""
    if pd.isna(r2):
        return "#BEBEBE"
    elif r2 <= 0.2:
        return "#4A68D9"
    elif r2 <= 0.4:
        return "#6CEAED"
    elif r2 <= 0.6:
        return "#5DCA3B"
    elif r2 <= 0.8:
        return "#F2A93B"
    else:
        return "#EA4025"


def get_summary_data(webdata_dir: str = "webdata"):
    """Load summary data from webdata directory."""
    summary_file = os.path.join(webdata_dir, "all_loci_info.txt")
    if os.path.exists(summary_file):
        return pd.read_csv(summary_file, sep="\t")
    return pd.DataFrame()


def layout(locus_id=None):
    """Create layout for locus page."""
    if locus_id is None:
        return html.Div("No locus selected")

    # Load summary data for the dropdowns
    df_summary = get_summary_data()
    if df_summary.empty:
        return html.Div("No data available")

    meta_methods = sorted(df_summary["meta_type"].unique())

    return dbc.Container(
        [
            dcc.Location(id="url", refresh=False),
            dbc.Row(
                [
                    html.H3(id="locus-title", className="mb-4"),
                ]
            ),
            dbc.Row(
                [
                    dbc.Col(
                        [
                            dmc.Select(
                                id="meta-method-select-locus",
                                label="Meta-analysis Method",
                                data=[
                                    {"value": method, "label": method}
                                    for method in meta_methods
                                ],
                                value=meta_methods[0] if meta_methods else None,
                                style={"width": "100%"},
                                className="mb-3",
                            ),
                        ],
                        width=4,
                    ),
                    dbc.Col(
                        [
                            dmc.Select(
                                id="fine-mapping-method-select-locus",
                                label="Fine-mapping Method",
                                data=[],  # Will be populated by callback
                                style={"width": "100%"},
                                className="mb-3",
                            ),
                        ],
                        width=4,
                    ),
                    dbc.Col(
                        [
                            dmc.Select(
                                id="sumstat-select-locus",
                                label="Summary Statistics",
                                data=[],  # Will be populated by callback
                                style={"width": "100%"},
                                className="mb-3",
                            ),
                        ],
                        width=4,
                    ),
                ]
            ),
            dbc.Row(
                [
                    dbc.Col(
                        [
                            dcc.Graph(id="locuszoom-plot"),
                        ],
                        width=9,
                    ),
                    dbc.Col(
                        [
                            dcc.Graph(id="observed-expected-plot"),
                        ],
                        width=3,
                    ),
                ],
                className="mb-4",
            ),
            dbc.Row(
                [
                    dbc.Col(
                        [
                            html.H4("QC Metrics", className="mb-3"),
                            dash_table.DataTable(
                                id="qc-table",
                                columns=[
                                    {"name": "Population", "id": "popu"},
                                    {"name": "Cohort", "id": "cohort"},
                                    {"name": "# SNPs", "id": "nsnp"},
                                    {"name": "# SNPs (p<1e-5)", "id": "nsnp_1e-5"},
                                    {"name": "# SNPs (p<5e-8)", "id": "nsnp_5e-8"},
                                    {"name": "Lambda", "id": "lambda"},
                                    {"name": "DENTIST-S", "id": "dentist-s"},
                                    {"name": "MAF Correlation", "id": "maf_corr"},
                                ],
                                style_table={"overflowX": "auto"},
                                style_cell={"textAlign": "left", "padding": "10px"},
                            ),
                        ],
                        width=12,
                    ),
                ]
            ),
        ]
    )


# Add simplified callbacks
try:

    @callback(Output("locus-title", "children"), Input("url", "pathname"))
    def update_locus_title(pathname):
        """Update locus title based on URL."""
        if pathname and "/locus/" in pathname:
            locus_id = pathname.split("/locus/")[-1]
            return f"Locus: {locus_id}"
        return "Locus View"

    @callback(
        [
            Output("fine-mapping-method-select-locus", "data"),
            Output("fine-mapping-method-select-locus", "value"),
        ],
        Input("meta-method-select-locus", "value"),
    )
    def update_fine_mapping_methods_locus(meta_method):
        """Update fine mapping methods for locus page."""
        if not meta_method:
            return [], None
        df_summary = get_summary_data()
        if df_summary.empty:
            return [], None
        methods = sorted(
            df_summary[df_summary["meta_type"] == meta_method]["tool"].unique()
        )
        return [{"value": method, "label": method} for method in methods], (
            methods[0] if methods else None
        )

except (ImportError, AttributeError):
    # Callbacks not available without dash
    pass
