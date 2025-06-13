"""Home page for CREDTOOLS web visualization."""

import os

try:
    import dash
    import dash_bootstrap_components as dbc
    import dash_mantine_components as dmc
    import pandas as pd
    import plotly.graph_objects as go
    from dash import dash_table  # type: ignore
    from dash import Input, Output, callback, dcc, html
except ImportError:
    # Web dependencies not available
    dash_table = None  # type: ignore

# Register this page
try:
    dash.register_page(__name__, path="/")
except (ImportError, AttributeError):
    pass


def get_summary_data(webdata_dir: str = "webdata"):
    """Load summary data from webdata directory."""
    summary_file = os.path.join(webdata_dir, "all_loci_info.txt")
    if os.path.exists(summary_file):
        return pd.read_csv(summary_file, sep="\t")
    return pd.DataFrame()


# Main layout
layout = dbc.Container(
    [
        dbc.Row(
            [
                dbc.Col(
                    [
                        dmc.Select(
                            id="meta-method-select",
                            label="Meta-analysis Method",
                            data=[],
                            value=None,
                            style={"width": "100%"},
                            className="mb-3",
                        ),
                    ],
                    width=4,
                ),
                dbc.Col(
                    [
                        dmc.Select(
                            id="fine-mapping-method-select",
                            label="Fine-mapping Method",
                            data=[],
                            value=None,
                            style={"width": "100%"},
                            className="mb-3",
                        ),
                    ],
                    width=4,
                ),
                dbc.Col(
                    [
                        dmc.Select(
                            id="plot-type-select",
                            label="Plot Type",
                            data=[
                                {
                                    "value": "n_credsets",
                                    "label": "Number of Credible Sets",
                                },
                                {
                                    "value": "whole_credsize",
                                    "label": "Whole Credible Set Size",
                                },
                                {
                                    "value": "n_PIP_gt_0.1",
                                    "label": "Number of SNPs with PIP > 0.1",
                                },
                            ],
                            value="n_credsets",
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
                        dcc.Graph(id="barplot"),
                    ],
                    width=12,
                ),
            ],
            className="mb-4",
        ),
        dbc.Row(
            [
                dbc.Col(
                    [
                        dash_table.DataTable(  # type: ignore
                            id="locus-table",
                            columns=[
                                {
                                    "name": "Locus ID",
                                    "id": "locus_link",
                                    "presentation": "markdown",
                                },
                                {"name": "Chr", "id": "chr"},
                                {"name": "Start", "id": "start"},
                                {"name": "End", "id": "end"},
                                {"name": "# SNPs", "id": "nsnp"},
                                {"name": "# SNPs (PIP>1e-5)", "id": "nsnp_1e-5"},
                                {"name": "# SNPs (PIP>5e-8)", "id": "nsnp_5e-8"},
                                {"name": "# Credible Sets", "id": "n_credsets"},
                                {"name": "Credible Set Size", "id": "whole_credsize"},
                                {"name": "# SNPs (PIP>0.1)", "id": "n_PIP_gt_0.1"},
                            ],
                            sort_action="native",
                            sort_mode="multi",
                            page_size=20,
                            style_table={"overflowX": "auto"},
                            style_cell={
                                "textAlign": "left",
                                "padding": "10px",
                                "whiteSpace": "normal",
                                "height": "auto",
                            },
                            style_header={
                                "backgroundColor": "rgb(230, 230, 230)",
                                "fontWeight": "bold",
                            },
                            style_data_conditional=[
                                {
                                    "if": {"row_index": "odd"},
                                    "backgroundColor": "rgb(248, 248, 248)",
                                }
                            ],
                        ),
                    ],
                    width=12,
                ),
            ]
        ),
    ]
)

try:

    @callback(
        Output("meta-method-select", "data"),
        Output("meta-method-select", "value"),
        Input("meta-method-select", "id"),  # Dummy input to trigger initial load
    )
    def populate_meta_methods(_):
        """Populate meta method dropdown with available options."""
        df_summary = get_summary_data()
        if df_summary.empty:
            return [], None
        meta_methods = sorted(df_summary["meta_type"].unique())
        return [{"value": method, "label": method} for method in meta_methods], (
            meta_methods[0] if meta_methods else None
        )

    @callback(
        [
            Output("fine-mapping-method-select", "data"),
            Output("fine-mapping-method-select", "value"),
        ],
        Input("meta-method-select", "value"),
    )
    def update_fine_mapping_methods(meta_method):
        """Update fine mapping methods based on selected meta method."""
        if not meta_method:
            return [], None
        df_summary = get_summary_data()
        if df_summary.empty:
            return [], None
        filtered_methods = df_summary[df_summary["meta_type"] == meta_method][
            "tool"
        ].unique()
        sorted_methods = sorted(filtered_methods)
        return [{"value": method, "label": method} for method in sorted_methods], (
            sorted_methods[0] if sorted_methods else None
        )

    @callback(
        [Output("barplot", "figure"), Output("locus-table", "data")],
        [
            Input("meta-method-select", "value"),
            Input("fine-mapping-method-select", "value"),
            Input("plot-type-select", "value"),
        ],
    )
    def update_plots_and_table(meta_method, fine_mapping_method, plot_type):
        """Update plots and table based on selected filters."""
        df_summary = get_summary_data()
        if df_summary.empty or not meta_method or not fine_mapping_method:
            return go.Figure(), []

        # Filter data based on selections
        filtered_df = df_summary[
            (df_summary["meta_type"] == meta_method)
            & (df_summary["tool"] == fine_mapping_method)
        ].copy()

        if filtered_df.empty:
            return go.Figure(), []

        filtered_df.sort_values(by=["nsnp_1e-5"], ascending=False, inplace=True)
        filtered_df.drop_duplicates(subset=["locus_id"], inplace=True)

        # Sort by chr and start position
        filtered_df = filtered_df.sort_values(["chr", "start"])

        # Add clickable links
        filtered_df["locus_link"] = filtered_df["locus_id"].apply(
            lambda x: f"[**{x}**](/locus/{x})"
        )

        # Create barplot
        fig = go.Figure()
        fig.add_trace(
            go.Bar(
                x=filtered_df["locus_id"],
                y=filtered_df[plot_type],
                marker_color="rgb(55, 83, 109)",
            )
        )

        title_map = {
            "n_credsets": "Number of Credible Sets",
            "whole_credsize": "Credible Set Size",
            "n_PIP_gt_0.1": "Number of SNPs with PIP > 0.1",
        }

        fig.update_layout(
            title=f"{title_map.get(plot_type, plot_type)} per Locus",
            xaxis_title="Locus ID",
            yaxis_title="Count",
            template="plotly_white",
            showlegend=False,
            xaxis={"tickangle": 45},
            margin={"t": 50, "b": 100},
        )

        # Prepare table data
        table_columns = [
            "locus_link",
            "chr",
            "start",
            "end",
            "popu",
            "cohort",
            "n_credsets",
            "whole_credsize",
            "n_PIP_gt_0.1",
            "nsnp",
            "nsnp_1e-5",
            "nsnp_5e-8",
        ]
        available_columns = [col for col in table_columns if col in filtered_df.columns]
        table_data = filtered_df[available_columns].to_dict("records")

        return fig, table_data

except (ImportError, AttributeError):
    # Callbacks not available without dash
    pass
