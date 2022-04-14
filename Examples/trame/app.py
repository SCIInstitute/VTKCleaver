import os
import tempfile

from trame import state
from trame.layouts import SinglePageWithDrawer
from trame.html import vtk, vuetify

from vtkmodules.vtkFiltersSources import vtkConeSource
from vtkmodules.vtkRenderingCore import (
    vtkActor,
    vtkPolyDataMapper,
    vtkRenderer,
    vtkRenderWindow,
    vtkRenderWindowInteractor,
)

# Required for interactor initialization
from vtkmodules.vtkInteractionStyle import vtkInteractorStyleSwitch  # noqa

# Required for rendering initialization, not necessary for
# local rendering, but doesn't hurt to include it
import vtkmodules.vtkRenderingOpenGL2  # noqa

# vtk
import vtk as vtkpackage

# vtkCleaver
try:
  from vtkCleaver import vtkCleaverImageToUnstructuredGridFilter
except ImportError:
  from vtk import vtkCleaverImageToUnstructuredGridFilter


# -----------------------------------------------------------------------------
# Constants
# -----------------------------------------------------------------------------

DEFAULT_SAMPLING_RATE = 1.0
DEFAULT_FEATURE_SCALING = 1.0
DEFAULT_RATE_OF_CHANGE = 0.2


# -----------------------------------------------------------------------------
# VTK pipeline
# -----------------------------------------------------------------------------

renderer = vtkRenderer()
render_window = vtkRenderWindow()
render_window.AddRenderer(renderer)

render_window_interactor = vtkRenderWindowInteractor()
render_window_interactor.SetRenderWindow(render_window)
render_window_interactor.GetInteractorStyle().SetCurrentStyleToTrackballCamera()

cleaver_mesher = vtkCleaverImageToUnstructuredGridFilter()

# remove the background mesh with label 0
threshold = vtkpackage.vtkThreshold()
threshold.SetInputConnection(cleaver_mesher.GetOutputPort())
threshold.SetThresholdFunction(threshold.THRESHOLD_UPPER)
threshold.SetUpperThreshold(0.99)

sF = vtkpackage.vtkDataSetSurfaceFilter()
sF.SetInputConnection(threshold.GetOutputPort())

mapper = vtkpackage.vtkPolyDataMapper()
mapper.SetInputConnection(sF.GetOutputPort())
mapper.ScalarVisibilityOn()
mapper.SetScalarModeToUseCellData()
mapper.SetColorModeToMapScalars()

actor = vtkActor()
actor.GetProperty().SetEdgeVisibility(True)
actor.SetMapper(mapper)

xmlReader = vtkpackage.vtkXMLImageDataReader()
nrrdReader = vtkpackage.vtkNrrdReader()

input_images = {}

state.cleaver_running = False
state.input_available = False
state.input_file_names = []

# -----------------------------------------------------------------------------
# Functions
# -----------------------------------------------------------------------------

@state.change("files")
def load_client_files(files, **kwargs):

    input_images.clear()
    
    if files and len(files):
        if not files[0].get("content"):
            return
        for file in files:
            fileName = file.get("name")
            print(f"Load {fileName}")
            ext = os.path.splitext(fileName)[1]
            reader = xmlReader
            if ext == ".nrrd":
                reader = nrrdReader
            
            with tempfile.NamedTemporaryFile(suffix=ext) as fp:
                fp.write(file.get("content"))
                fp.seek(0)

                reader.SetFileName(fp.name)
                reader.Update()

                img = vtkpackage.vtkImageData()
                img.ShallowCopy(reader.GetOutput())

                input_images[fileName] = img

    state.input_file_names = list(input_images.keys())


def cleave_inputs():
    state.cleaver_running = True

    cleaver_mesher.RemoveAllInputs()
    renderer.RemoveActor(actor)

    for fileName in state.input_file_names:
        cleaver_mesher.AddInputData(0, input_images[fileName])
    
    print(state.rate_of_change)

    cleaver_mesher.SetSamplingRate(state.sampling_rate)
    cleaver_mesher.SetFeatureScaling(state.feature_scaling)
    cleaver_mesher.SetRateOfChange(state.rate_of_change)
    cleaver_mesher.Update()

    ugrid = cleaver_mesher.GetOutput()
    mapper.SetScalarRange(ugrid.GetCellData().GetScalars().GetRange())

    renderer.AddActor(actor)
    renderer.ResetCamera()

    state.cleaver_running = False

    html_view.update()


@state.change("input_file_names")
def update_input_available(input_file_names, **kwargs):
    state.input_available = len(input_file_names) > 0
    print(f"update_input_available: {state.input_available}")


@state.change("sampling_rate")
def update_sampling_rate(sampling_rate, **kwargs):
    print(f"sampling_rate: {sampling_rate}")


def reset_sampling_rate():
    state.sampling_rate = DEFAULT_SAMPLING_RATE


@state.change("feature_scaling")
def update_feature_scaling(feature_scaling, **kwargs):
    print(f"feature_scaling: {feature_scaling}")


def reset_feature_scaling():
    state.feature_scaling = DEFAULT_FEATURE_SCALING


@state.change("rate_of_change")
def update_rate_of_change(rate_of_change, **kwargs):
    print(f"rate_of_change: {rate_of_change}")


def reset_rate_of_change():
    state.rate_of_change = DEFAULT_RATE_OF_CHANGE


# -----------------------------------------------------------------------------
# GUI Cards
# -----------------------------------------------------------------------------

compact_style = {
    "hide_details": True,
    "dense": True,
}


def ui_card(title, ui_name):
    card_style = {}
    title_style = {
        "classes": "grey lighten-1 py-1 grey--text text--darken-3",
        "style": "user-select: none; cursor: pointer",
        **compact_style,
    }
    content_style = {"classes": "py-2"}
    with vuetify.VCard(**card_style):
        vuetify.VCardTitle(title, **title_style)
        content = vuetify.VCardText(**content_style)

    return content

def input_card():
    with ui_card("Input", "input"):
        vuetify.VFileInput(
            multiple=True,
            show_size=True,
            small_chips=True,
            truncate_length=25,
            v_model=("files", None),
            dense=True,
            hide_details=True,
            style="max-width: 300px;",
            accept=".nrrd,.vti",
            __properties=["accept"],
        )
    
    # TODO Add a component allowing to re-order "state.input_file_names"
    # Possible approaches:
    # - v-simple-table + SortableJS
    #   See https://codepen.io/mykysyk/pen/qBdBRMB
    # - v-chip + draggable
    #   See https://github.com/vuetifyjs/vuetify/issues/11614 and https://jsfiddle.net/cjrqzkf0/


def sizing_field_row(**kwargs):
    with vuetify.VRow():
        vuetify.VSlider(
            **kwargs,
            thumb_label=True,
            classes="my-1",
            **compact_style,
        )
        with vuetify.VBtn(icon=True, click=reset_sampling_rate, **compact_style):
            vuetify.VIcon("mdi-restore")


def sizing_field_card():
    with ui_card("Sizing Field", "sizing-field"):
        sizing_field_row(
            v_model=("sampling_rate", DEFAULT_SAMPLING_RATE),
            min=0.01,
            max=99.9,
            step=1.0,
            label="Sampling Rate"
        )
        sizing_field_row(
            v_model=("feature_scaling", DEFAULT_FEATURE_SCALING),
            min=0.0001,
            max=9999.0,
            step=0.0001,
            label="Feature Scaling"
        )
        sizing_field_row(
            v_model=("rate_of_change", DEFAULT_RATE_OF_CHANGE),
            min=0.0001,
            max=10.0,
            step=0.0001,
            label="Rate of Change (Lipschitz)"
        )


def cleaving_tool_card():
    with ui_card("Cleaving Tool", "cleaving-tool"):
        with vuetify.VRow(justify="center", **compact_style):
            vuetify.VBtn(
                "Cleave Mesh",
                disabled=("!input_available || cleaver_running",),
                loading=("cleaver_running", True),
                click=cleave_inputs,
                **compact_style,
            )

# -----------------------------------------------------------------------------
# GUI
# -----------------------------------------------------------------------------

html_view = vtk.VtkLocalView(render_window)

layout = SinglePageWithDrawer("Cleaver", on_ready=html_view.update)
layout.title.set_text("Cleaver")

with layout.drawer as drawer:
    drawer.width = 325
    vuetify.VDivider()
    input_card()
    vuetify.VDivider()
    sizing_field_card()
    vuetify.VDivider()
    cleaving_tool_card()

with layout.content:
    vuetify.VContainer(
        fluid=True,
        classes="pa-0 fill-height",
        children=[html_view],
    )

# -----------------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------------

if __name__ == "__main__":
    layout.start()
