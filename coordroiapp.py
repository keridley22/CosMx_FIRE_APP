import dash
from dash import dcc, html, Input, Output
import plotly.graph_objs as go
import pandas as pd
import os
import glob

# --- Load Metadata ---
metadata = pd.read_csv("/Users/katherineridley/Projects/CosMx/APP/B1_metadata.csv").head(150000)
x_min = metadata['x_slide_mm'].min()
x_max = metadata['x_slide_mm'].max()
y_min = metadata['y_slide_mm'].min()
y_max = metadata['y_slide_mm'].max()

# Create a base scatter trace for metadata
metadata_trace = go.Scatter(
    x = metadata['x_slide_mm'],
    y = metadata['y_slide_mm'],
    mode = 'markers',
    marker = dict(size=1, color='blue', opacity=0.7, line=dict(width=0)),
    text = [f"({x:.2f}, {y:.2f})" for x, y in zip(metadata['x_slide_mm'], metadata['y_slide_mm'])],
    hoverinfo = 'text',
    name = 'Metadata'
)

def get_roi_traces():
    """Search the ROI folder for _transformed.csv files, load them, and return a list of scatter traces."""
    roi_dir = "/Users/katherineridley/Projects/CosMx/APP/B1_rois/"
    pattern = os.path.join(roi_dir, "*_transformed.csv")
    roi_files = glob.glob(pattern)
    roi_traces = []
    for file in roi_files:
        df = pd.read_csv(file)
        # Ensure columns are trimmed and lowercased if needed.
        df.columns = [col.strip() for col in df.columns]
        # Expecting columns 'x_mm' and 'y_mm'
        if 'x_mm' not in df.columns or 'y_mm' not in df.columns:
            print(f"Skipping {file}: Missing 'x_mm' or 'y_mm'")
            continue
        trace = go.Scatter(
            x = df['x_mm'],
            y = df['y_mm'],
            mode = 'markers',
            marker = dict(size=3, color='red', opacity=0.7),
            name = os.path.basename(file)
        )
        roi_traces.append(trace)
    return roi_traces

# Initial figure with only metadata.
initial_fig = go.Figure(
    data=[metadata_trace],
    layout=go.Layout(
        title="Metadata Scatter Plot",
        xaxis=dict(title="x_slide_mm", range=[x_min, x_max]),
        yaxis=dict(title="y_slide_mm", range=[y_min, y_max], scaleanchor="x", scaleratio=1),
        width=1500,
        height=2000,
        hovermode="closest"
    )
)

# --- Dash App Setup ---
app = dash.Dash(__name__)
app.title = "Metadata Scatter Plot with ROI Overlay"

app.layout = html.Div([
    html.H1("Metadata Scatter Plot"),
    dcc.Checklist(
        id='toggle-roi',
        options=[{'label': 'Show ROI Coordinates', 'value': 'show'}],
        value=[],
        labelStyle={'display': 'inline-block', 'marginRight': '20px'}
    ),
    dcc.Graph(id='scatter-plot', figure=initial_fig, style={'width': '100%', 'height': '90vh'})
])

@app.callback(
    Output('scatter-plot', 'figure'),
    Input('toggle-roi', 'value')
)
def update_figure(toggle_value):
    # Start with the metadata trace
    traces = [metadata_trace]
    if 'show' in toggle_value:
        roi_traces = get_roi_traces()
        traces.extend(roi_traces)
    fig = go.Figure(
        data=traces,
        layout=go.Layout(
            title="Metadata Scatter Plot",
            xaxis=dict(title="x_slide_mm", range=[x_min, x_max]),
            yaxis=dict(title="y_slide_mm", range=[y_min, y_max], scaleanchor="x", scaleratio=1),
            width=1500,
            height=2000,
            hovermode="closest"
        )
    )
    return fig

if __name__ == '__main__':
    app.run(debug=True, port=8060)
