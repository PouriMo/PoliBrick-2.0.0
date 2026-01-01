![Plugin_logo](Resources/PBLogo.JPG)
# PoliBrick

PoliBrick is a Grasshopper plugin that streamlines the export of complex geometries (like brick assemblies) into **Itasca 3DEC 5.2** format. 

It handles the generation of `.dat` files with correct Polyhedra syntax, ensuring watertight geometry for numerical modeling.

![Example Workflow](Resources/PBFamily.png)
*(Put a nice screenshot of the component in Grasshopper here)*

## Features
- Exports Grasshopper BRep/Mesh geometry to 3DEC Polyhedra blocks.
- **Auto-correction**: Ensures correct vertex winding (CCW) for faces to avoid "Negative Volume" errors.
- Assigns material (`mat`) and constitutive model (`con`) IDs.

## Installation
1. Go to the [Releases](../../releases) page.
2. Download the `PoliBrick.gha` file.
3. In Grasshopper, go to `File > Special Folders > Components Folder`.
4. Paste the `.gha` file there.
5. Right-click the file > **Properties** > **Unblock** (if necessary).
6. Restart Rhino/Grasshopper.

## Usage
1. **Input Geometry:** Connect your closed Breps or Meshes to the `Geo` input.
2. **Settings:** define your Material ID and Constitutive Model ID.
3. **Export:** Connect a panel or a "Text Save" component to the output to generate your `.dat` file.




## Example
Check the `examples/` folder for a sample:
- `wall_test.gh`: The Grasshopper generation script.
- `wall_test.dat`: The resulting 3DEC input file.

## Requirements
- Rhinoceros 6 or 7
- Itasca 3DEC 5.2 (for running the output)

## License
This project is licensed under the MIT License.
