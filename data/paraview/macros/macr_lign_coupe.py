from __future__ import annotations
from packaging.version import Version
import enum
import math
import medcoupling as mc

# state file generated using paraview version 5.11.2
import paraview

PV_VERSION = Version(paraview.__version__)
paraview.compatibility.major = 5
paraview.compatibility.minor = 11

#### import the simple module from the paraview
from paraview.simple import (
    Glyph,
    GetSettingsProxy,
    GetActiveSource,
    GetScalarBar,
    Calculator,
    ExtractGroup,
    GetActiveView,
    GetColorTransferFunction,
    ColorBy,
    HideScalarBarIfNotNeeded,
    PythonAnnotation,
    Render,
    Show,
    Hide,
    UpdatePipeline,
    
)

#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()


PREFIX_CUT = "CUT_"


class AxisFields(enum.Enum):
    VECT_X = f"{PREFIX_CUT}AXE_X"
    VECT_Y = f"{PREFIX_CUT}AXE_Y"
    VECT_Z = f"{PREFIX_CUT}AXE_Z"


AXIS_LOCATION_MASK = f"{PREFIX_CUT}AXE_MASK"

AXIS_TO_COLOR = {
    AxisFields.VECT_X: [1.0, 0.0, 0.0],
    AxisFields.VECT_Y: [1.0, 1.0, 0.0],
    AxisFields.VECT_Z: [0.0, 1.0, 0.0],
}


def get_available_point_arrays(source) -> dict[str, "PointData"]:
    return {point_data.Name: point_data for point_data in source.PointData if point_data.Name.startswith(PREFIX_CUT)}


def get_available_components(point_fields_data, field_name: str) -> list[str]:
    ai = point_fields_data.GetArrayInformation(field_name)
    return [ai.GetComponentName(compo_num) for compo_num in range(ai.GetNumberOfComponents())]
def create_cut_labels(base_source, view):
    parent_filter = base_source

    groups_source = ExtractGroup(
        registrationName="LABELS",
        Input=parent_filter
    )

    groups_display = Show(
        groups_source,
        view,
        'UnstructuredGridRepresentation'
    )
    groups_display.Representation = 'Surface'

    file_name = base_source.FileNames[0]
    data = mc.MEDFileMesh.New(file_name)

    group_names = data.getGroupsNames()
    coords = data.getCoords()

    # Number of coordinate components: 2 for a 2D mesh, 3 for a 3D mesh
    space_dim = coords.getNumberOfComponents()

    print("Number of coordinate components:", space_dim)

    labels = []

    for g in group_names:

        if not g.startswith("CUT_"):
            continue

        nodes = data.getNodeGroupArr(g)

        # Security
        if nodes is None or len(nodes) == 0:
            continue

        group_coords = coords[nodes]

        # ------------------------------------------------------------
        # Compute the center
        # ------------------------------------------------------------
        center = [
            sum(p[i] for p in group_coords) / len(group_coords)
            for i in range(space_dim)
        ]

        # ------------------------------------------------------------
        # Look for the farthest point from the center
        # ------------------------------------------------------------
        max_dist = -1.0
        end_point = group_coords[0]

        for p in group_coords:

            dist = sum(
                (p[i] - center[i]) ** 2
                for i in range(space_dim)
            )

            if dist > max_dist:
                max_dist = dist
                end_point = p

        # ------------------------------------------------------------
        # Coordinates of the selected point
        # ------------------------------------------------------------
        cx = end_point[0]
        cy = end_point[1]

        # In 2D, z = 0
        if space_dim >= 3:
            cz = end_point[2]
        else:
            cz = 0.0

        # ------------------------------------------------------------
        # Bounding box
        # ------------------------------------------------------------
        xmin = min(p[0] for p in group_coords)
        xmax = max(p[0] for p in group_coords)

        ymin = min(p[1] for p in group_coords)
        ymax = max(p[1] for p in group_coords)

        if space_dim >= 3:
            zmin = min(p[2] for p in group_coords)
            zmax = max(p[2] for p in group_coords)
        else:
            zmin = 0.0
            zmax = 0.0

        # ------------------------------------------------------------
        # Diagonal of the bounding box
        # ------------------------------------------------------------
        group_diag = math.sqrt(
            (xmax - xmin) ** 2 +
            (ymax - ymin) ** 2 +
            (zmax - zmin) ** 2
        )

        # Avoid zero offset for degenerate groups
        if group_diag < 1e-12:
            group_diag = 1.0

        # ------------------------------------------------------------
        # Position with respect to the camera
        # ------------------------------------------------------------
        camera = view.GetActiveCamera()
        cam_pos = camera.GetPosition()

        vx = cam_pos[0] - cx
        vy = cam_pos[1] - cy
        vz = cam_pos[2] - cz

        norm = math.sqrt(
            vx * vx +
            vy * vy +
            vz * vz
        )

        if norm > 1e-12:
            vx /= norm
            vy /= norm
            vz /= norm

        # ------------------------------------------------------------
        # Small offset in front of the line
        # ------------------------------------------------------------
        eps = 0.005 * group_diag

        cx2 = cx + eps * vx
        cy2 = cy + eps * vy
        cz2 = cz + eps * vz

        # ------------------------------------------------------------
        # Add text
        # ------------------------------------------------------------
        python_annotation = PythonAnnotation(
            registrationName=g.upper(),
            Input=groups_source
        )

        python_annotation.Expression = '"%s"' % g.upper()

        python_annotation_display = Show(
            python_annotation,
            view,
            'TextSourceRepresentation'
        )

        python_annotation_display.TextPropMode = 'Billboard 3D Text'

        python_annotation_display.BillboardPosition = [
            cx2,
            cy2,
            cz2
        ]

        python_annotation_display.Color = [
            0.0,
            0.0,
            0.0
        ]

        python_annotation_display.FontSize = 15

        labels.append({
            "source": python_annotation,
            "display": python_annotation_display,
            "center": (cx2, cy2, cz2),
            "size": group_diag
        })

    # ------------------------------------------------------------
    # Manage the size of the labels
    # ------------------------------------------------------------
    camera = view.GetActiveCamera()

    reference_scale = camera.GetParallelScale()
    reference_fontsize = 15

    def update_labels(caller=None, event=None):

        camera = view.GetActiveCamera()

        scale = camera.GetParallelScale()

        if scale < 1e-12:
            return

        fontsize = int(
            reference_fontsize *
            reference_scale /
            scale
        )

        fontsize = max(8, min(fontsize, 11))

        for label in labels:

            label["display"].FontSize = fontsize
            label["display"].UpdateVTKObjects()

        view.StillRender()

    # Hide labels initially
    for label in labels:
        Hide(label["source"], view)

    # Update labels when the view is rendered
    interactor = view.GetInteractor()

    interactor.AddObserver(
        "RenderEvent",
        update_labels
    )

    update_labels()

    # Hide the group
    Hide(groups_source, view)

    return labels

def main():
    base_source = GetActiveSource()
    view = GetActiveView()
    Hide(base_source, view)

    palette = GetSettingsProxy("ColorPalette")
    edge_color =palette.Foreground
    point_color = [0.5 * val for val in edge_color]

    available_points_arrays = get_available_point_arrays(source=base_source)
    available_axis_arrays = tuple((field.value for field in AxisFields if field.value in available_points_arrays))
    available_non_axis_arrays = list((field for field in available_points_arrays if field not in available_axis_arrays))

    if AXIS_LOCATION_MASK in available_non_axis_arrays:
        available_non_axis_arrays.remove(AXIS_LOCATION_MASK)

    xmin, xmax, ymin, ymax, zmin, zmax = base_source.GetDataInformation().DataInformation.GetBounds()
    max_bound = max((xmax-xmin, ymax - ymin, zmax - zmin))

    point_fields_data = base_source.PointData.GetFieldData()

    parent_filter = base_source


    line_source = Calculator(registrationName="LINES", Input=parent_filter)
    line_display = Show(line_source, view, )
    line_display.SetRepresentationType("Wireframe")
    line_display.RenderLinesAsTubes = 1
    line_display.LineWidth = 7.0
    # change solid color
    line_display.AmbientColor = edge_color
    line_display.DiffuseColor = edge_color
    point_source = Calculator(registrationName="POINTS", Input=parent_filter)
    point_display = Show(point_source, view, )
    point_display.SetRepresentationType("Points")
    point_display.RenderPointsAsSpheres = 1
    point_display.PointSize = 10.0
    point_display.AmbientColor = point_color
    point_display.DiffuseColor = point_color

    Hide(point_source, view)

    if available_axis_arrays:
        glyph_scale = Calculator(registrationName="SET GLYPH SIZE", Input=parent_filter)

        glyph_scale.Function = f"{max_bound / 10.} * {AXIS_LOCATION_MASK}"
        glyph_scale.ResultArrayName = "GLYPH_SCALE"

        parent_filter = glyph_scale

    for field_name in available_axis_arrays:
        # create a new 'Glyph'
        glyph = Glyph(
            registrationName=field_name,
            Input=parent_filter,
            GlyphType="Arrow",
        )
        glyph.OrientationArray = ["POINTS", field_name]
        glyph.GlyphTransform = "Transform2"
        glyph.GlyphMode = "All Points"
        glyph.ScaleArray = ['POINTS', 'GLYPH_SCALE']
        glyph.ScaleFactor = 1
        # get display properties
        glyph_display = Show(glyph, view, "UnstructuredGridRepresentation")
        # disable scalar coloring
        ColorBy(glyph_display, None)
        color = AXIS_TO_COLOR[AxisFields(field_name)]
        # change solid color
        glyph_display.AmbientColor = color
        glyph_display.DiffuseColor = color
        Hide(glyph, view)

    for field_name in available_non_axis_arrays:
        components = get_available_components(point_fields_data=point_fields_data, field_name=field_name)
        # create a new 'Extract Component'
        for compo_name in components:
            output_name = f"{field_name}_{compo_name}"
            extract_compo = Calculator(registrationName=output_name, Input=parent_filter)
            extract_compo.Function = output_name
            extract_compo.ResultArrayName = output_name
            extract_compo.UpdatePipeline()
            # show data in view
            extract_compo_display = Show(extract_compo, view, "UnstructuredGridRepresentation")
            if PV_VERSION >= Version("5.13.0"):
                extract_compo_display.DisableLighting = 1
            extract_compo_display.SetRepresentationType("Wireframe")
            extract_compo_display.RenderLinesAsTubes = 1
            extract_compo_display.RenderPointsAsSpheres = 1
            extract_compo_display.LineWidth = 10.0
            extract_compo_display.PointSize = 20.0

            extract_compo_display.GaussianRadius = 0.2
            extract_compo_display.ShaderPreset = "Sphere"

            # set scalar coloring
            ColorBy(extract_compo_display, ("POINTS", output_name))

            # get separate color transfer function/color map
            extract_compo_lut = GetColorTransferFunction(output_name, extract_compo_display)

            # Apply a preset using its name.
            extract_compo_lut.ApplyPreset("Turbo", True)
            extract_compo_lut.NumberOfTableValues = 12
            extract_compo_lut.RescaleOnVisibilityChange = 1
            extract_compo_lut.AutomaticRescaleRangeMode = 'Clamp and update every timestep'

            # Scalar bar settings
            colorbar = GetScalarBar(extract_compo_lut, view)
            colorbar.ScalarBarLength = 0.7
            colorbar.DrawTickMarks = 0
            colorbar.DrawTickLabels = 0
            colorbar.AddRangeLabels = 0
            colorbar.AddRangeAnnotations = 1
            colorbar.AutomaticAnnotations = 1
            colorbar.WindowLocation = 'Lower Right Corner'

            # Hide the scalar bar for this color map if no visible data is colored by it.
            HideScalarBarIfNotNeeded(extract_compo_lut, view)

            # rescale color and/or opacity maps used to include current data range
            extract_compo_display.RescaleTransferFunctionToDataRange(True)

            # show color bar/color legend
            extract_compo_display.SetScalarBarVisibility(view, True)
            Hide(extract_compo, view)

    create_cut_labels(base_source, view)
    Render(view)
    UpdatePipeline()


main()
