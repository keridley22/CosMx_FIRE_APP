import dash
from dash import dcc, html, Input, Output, State
import pandas as pd
import plotly.graph_objects as go
import numpy as np
import os
import json

# Load your data
metadata_1 = pd.read_csv("/Users/katherineridley/Projects/CosMx/APP/B1_metadata.csv")
metadata_2 = pd.read_csv("/Users/katherineridley/Projects/CosMx/APP/A2_metadata.csv")

# Add 'slide_id' column
metadata_1['slide_id'] = 'B1'
metadata_2['slide_id'] = 'A2'

# Combine data into a single DataFrame
combined_metadata = pd.concat([metadata_1, metadata_2], ignore_index=True)

# Set up the output directory for ROIs
output_directory = "/Users/katherineridley/Projects/CosMx/APP/region_rois"
os.makedirs(output_directory, exist_ok=True)

# Initialize the Dash app
app = dash.Dash(__name__)
app.title = "Interactive ROI Annotation"

# Dropdown to select slide ID
slide_ids = combined_metadata['slide_id'].unique()

app.layout = html.Div([
    html.H1("Scatterplot with ROI Drawing"),
    
    html.Label("Select Slide:"),
    dcc.Dropdown(
        id='slide-dropdown',
        options=[{'label': slide, 'value': slide} for slide in slide_ids],
        value=slide_ids[0]
    ),
    
    dcc.Graph(id='scatter-plot', config={'modeBarButtonsToAdd': ['drawclosedpath', 'eraseshape']}),
    
    html.Label("Enter ROI Name:"),
    dcc.Input(id='roi-name-input', type='text', placeholder='Enter ROI name'),
    
    html.Button('Save ROIs', id='save-button'),
    html.Div(id='save-message')
])

# Callback to update the scatterplot based on selected slide
@app.callback(
    Output('scatter-plot', 'figure'),
    Input('slide-dropdown', 'value')
)
def update_scatter_plot(selected_slide):
    # Filter data for the selected slide
    slide_data = combined_metadata[combined_metadata['slide_id'] == selected_slide]

    if len(slide_data) > 250000:
        slide_data = slide_data.sample(250000)

    # Ensure columns are numeric
    slide_data['x_slide_mm'] = pd.to_numeric(slide_data['x_slide_mm'], errors='coerce')
    slide_data['y_slide_mm'] = pd.to_numeric(slide_data['y_slide_mm'], errors='coerce')

    # Remove rows with NaN values
    slide_data = slide_data.dropna(subset=['x_slide_mm', 'y_slide_mm'])

    # If no data is found, return an empty figure
    if slide_data.empty:
        return go.Figure()

    # Create a scatter plot
    fig = go.Figure()
    fig.add_trace(go.Scattergl(
        x=slide_data['x_slide_mm'],
        y=slide_data['y_slide_mm'],
        mode='markers',
        marker=dict(size=1.2, color='blue', opacity=0.5),
        name='Cells'
    ))

    # Set layout properties using 'drawclosedpath'
    fig.update_layout(
        title=f"Slide: {selected_slide}",
        xaxis_title='X (mm)',
        yaxis_title='Y (mm)',
        dragmode='drawclosedpath',
        newshape=dict(line_color='red'),
        height=800,
        width=800
    )

    return fig

# Callback to save drawn ROIs

@app.callback(
    Output('save-message', 'children'),
    Input('save-button', 'n_clicks'),
    State('scatter-plot', 'relayoutData'),
    State('slide-dropdown', 'value'),
    State('roi-name-input', 'value')
)
def save_rois(n_clicks, relayout_data, selected_slide, roi_name):
    # Debugging: Print the input values
    print(f"Button Clicks: {n_clicks}")
    print(f"ROI Name: {roi_name}")
    print(f"Selected Slide: {selected_slide}")
    print(f"Relayout Data: {json.dumps(relayout_data, indent=2) if relayout_data else 'None'}")

    # Check if any inputs are missing
    if n_clicks is None or relayout_data is None or not roi_name:
        return "No ROIs drawn or ROI name missing."

    # Extract the last drawn shape from relayoutData
    shapes = relayout_data.get('shapes', [])
    
    if not shapes:
        return "No shapes found in the relayout data."

    # Focus only on the most recent shape
    last_shape = shapes[-1]  # Get the last drawn shape
    roi_data = []

    if last_shape['type'] == 'path':
        path_data = last_shape['path'].replace('M', '').replace('Z', '').split('L')
        for coord in path_data:
            try:
                x, y = coord.split(',')
                roi_data.append({
                    "x": float(x),
                    "y": float(y),
                    "roi_name": roi_name,
                    "slide_id": selected_slide
                })
            except ValueError:
                print(f"Skipping invalid coordinate: {coord}")
                continue

    # Check if any ROI data was collected
    if not roi_data:
        return "No valid ROI coordinates found for the most recent shape."

    # Save ROI data to CSV
    roi_csv_path = os.path.join(output_directory, f"{selected_slide}_rois.csv")
    
    # If the CSV file already exists, append to it
    if os.path.exists(roi_csv_path):
        existing_data = pd.read_csv(roi_csv_path)
        new_data = pd.DataFrame(roi_data)
        combined_data = pd.concat([existing_data, new_data], ignore_index=True)
        combined_data.to_csv(roi_csv_path, index=False)
    else:
        # Create a new CSV file if it doesn't exist
        new_data = pd.DataFrame(roi_data)
        new_data.to_csv(roi_csv_path, index=False)

    return f"Most recent ROI saved at {roi_csv_path}."




if __name__ == '__main__':
    app.run(debug=True)
