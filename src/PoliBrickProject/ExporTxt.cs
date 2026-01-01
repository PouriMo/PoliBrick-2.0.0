using Grasshopper;
using Grasshopper.Kernel;
using Grasshopper.Kernel.Data;
using Grasshopper.Kernel.Types;
using Grasshopper.Kernel.Parameters;
using Rhino.Geometry;
using System;
using System.Collections.Generic;
using System.Drawing;
using System.IO;
using System.Linq;
using System.Text;
using System.Text.RegularExpressions;

namespace PoliBrick.Util
{
    // ==========================================================================================
    // COMPONENT 1: IMPORT (READ ONLY)
    // ==========================================================================================
    public class ImportTxt : GH_Component
    {
        public ImportTxt()
          : base("4.Import Txt Mesh", "InTxt",
              "Import mesh data from PB .txt or 3DEC .3ddat files",
              "PoliBrick", "Utilities")
        {
        }

        public override GH_Exposure Exposure => GH_Exposure.secondary;

        protected override void RegisterInputParams(GH_InputParamManager pManager)
        {
            var param = new Param_FilePath();
            param.Name = "Path";
            param.NickName = "P";
            param.Description = "The file path to import from";
            param.Access = GH_ParamAccess.item;
            param.FileFilter = "All Supported|*.txt;*.3ddat|Text Files (*.txt)|*.txt|3DEC Files (*.3ddat)|*.3ddat|All Files (*.*)|*.*";
            pManager.AddParameter(param);

            pManager.AddBooleanParameter("Read", "R", "Set to true to read the file", GH_ParamAccess.item, false);
        }

        protected override void RegisterOutputParams(GH_OutputParamManager pManager)
        {
            pManager.AddMeshParameter("Meshes", "M", "Imported meshes", GH_ParamAccess.tree);
            pManager.AddTextParameter("Info", "I", "Status information", GH_ParamAccess.item);
        }

        protected override void SolveInstance(IGH_DataAccess DA)
        {
            string path = string.Empty;
            bool read = false;

            if (!DA.GetData(0, ref path)) return;
            DA.GetData(1, ref read);

            if (!read)
            {
                DA.SetData(1, "Inactive. Set Read to true.");
                return;
            }

            DataTree<Mesh> outputMeshes;
            string info;

            string extension = Path.GetExtension(path).ToLowerInvariant();
            if (extension == ".3ddat")
            {
                ImportExportShared.DoImport3DDAT(path, out outputMeshes, out info);
            }
            else
            {
                ImportExportShared.DoImport(path, out outputMeshes, out info);
            }

            var outputMeshesGh = new GH_Structure<GH_Mesh>();
            foreach (var ghPath in outputMeshes.Paths)
            {
                var meshes = outputMeshes.Branch(ghPath);
                foreach (var mesh in meshes)
                {
                    outputMeshesGh.Append(new GH_Mesh(mesh), ghPath);
                }
            }

            DA.SetDataTree(0, outputMeshesGh);
            DA.SetData(1, info);
        }

        public static Bitmap ConvertByteArrayToBitmap(byte[] imageData)
        {
            using (MemoryStream ms = new MemoryStream(imageData)) { return new Bitmap(ms); }
        }

        protected override Bitmap Icon => ConvertByteArrayToBitmap(PoliBrick.Properties.Resources.ImportIcon);
        public override Guid ComponentGuid => new Guid("11111111-AAAA-BBBB-CCCC-000000000001");
    }

    // ==========================================================================================
    // COMPONENT 2: EXPORT (WRITE ONLY)
    // ==========================================================================================
    public class ExportTxt : GH_Component
    {
        public ExportTxt()
          : base("3.Export Txt Mesh", "OutTxt",
              "Export mesh data to PB .txt or 3DEC .3ddat files",
              "PoliBrick", "Utilities")
        {
        }

        public override GH_Exposure Exposure => GH_Exposure.secondary;

        protected override void RegisterInputParams(GH_InputParamManager pManager)
        {
            pManager.AddMeshParameter("Meshes", "M", "Meshes to export (tree structure)", GH_ParamAccess.tree);

            var param = new Param_FilePath();
            param.Name = "Path";
            param.NickName = "P";
            param.Description = "Location to save the file";
            param.Access = GH_ParamAccess.item;
            param.FileFilter = "Text Files (*.txt)|*.txt|3DEC Files (*.3ddat)|*.3ddat";
            pManager.AddParameter(param);

            pManager.AddIntegerParameter("Format", "F",
                "0 = Standard (with Topology)\n1 = Straus (Geometry only)\n2 = 3DEC 5.2 (.3ddat, one file per branch)",
                GH_ParamAccess.item, 0);

            pManager.AddNumberParameter("Threshold", "T", "Distance tolerance to merge nearby nodes (for formats 0,1)",
                GH_ParamAccess.item, 0.001);

            pManager.AddBooleanParameter("Write", "W", "Set to true to write the file", GH_ParamAccess.item, false);
        }

        protected override void RegisterOutputParams(GH_OutputParamManager pManager)
        {
            pManager.AddTextParameter("Info", "I", "Status information", GH_ParamAccess.item);
        }

        protected override void SolveInstance(IGH_DataAccess DA)
        {
            var inputMeshesGh = new GH_Structure<GH_Mesh>();
            if (!DA.GetDataTree(0, out inputMeshesGh)) return;

            string path = string.Empty;
            if (!DA.GetData(1, ref path)) return;

            int formatMode = 0;
            DA.GetData(2, ref formatMode);

            double threshold = 0.001;
            DA.GetData(3, ref threshold);

            bool write = false;
            DA.GetData(4, ref write);

            if (!write)
            {
                DA.SetData(0, "Inactive. Set Write to true.");
                return;
            }

            var inputMeshes = new DataTree<Mesh>();
            foreach (var ghPath in inputMeshesGh.Paths)
            {
                var ghMeshes = inputMeshesGh.get_Branch(ghPath);
                foreach (GH_Mesh ghMesh in ghMeshes)
                {
                    if (ghMesh != null && ghMesh.Value != null)
                        inputMeshes.Add(ghMesh.Value, ghPath);
                }
            }

            string info;

            if (formatMode == 2)
            {
                ImportExportShared.DoExport3DDAT(inputMeshes, path, out info);
            }
            else
            {
                bool toStraus = (formatMode == 1);
                ImportExportShared.DoExport(inputMeshes, path, toStraus, threshold, out info);
            }

            DA.SetData(0, info);
        }

        public static Bitmap ConvertByteArrayToBitmap(byte[] imageData)
        {
            using (MemoryStream ms = new MemoryStream(imageData)) { return new Bitmap(ms); }
        }

        protected override Bitmap Icon => ConvertByteArrayToBitmap(PoliBrick.Properties.Resources.ExportIcon);
        public override Guid ComponentGuid => new Guid("22222222-AAAA-BBBB-CCCC-000000000002");
    }

    // ==========================================================================================
    // SHARED LOGIC (STATIC HELPERS)
    // ==========================================================================================
    public static class ImportExportShared
    {
        public class ElementData
        {
            public string ElementType { get; set; }
            public int BranchID { get; set; }
            public List<int> NodeIDs { get; set; }
        }

        // -------------------------
        // IMPORT LOGIC (STANDARD .txt)
        // -------------------------
        public static void DoImport(string P, out DataTree<Mesh> R, out string A)
        {
            R = new DataTree<Mesh>();
            A = string.Empty;

            if (string.IsNullOrEmpty(P) || !File.Exists(P))
            {
                A = "Error: File path is invalid or file does not exist.";
                return;
            }

            var timer = System.Diagnostics.Stopwatch.StartNew();
            var allNodes = new Dictionary<int, Point3d>();
            var allElements = new List<ElementData>();
            string currentSection = "";
            var outputMeshes = new DataTree<Mesh>();

            try
            {
                using (var reader = new StreamReader(P))
                {
                    string line;
                    while ((line = reader.ReadLine()) != null)
                    {
                        line = line.Trim();
                        if (string.IsNullOrWhiteSpace(line)) continue;

                        if (line.StartsWith("/"))
                        {
                            if (line.Contains("NODE COORDINATES")) currentSection = "NODE";
                            else if (line.Contains("BRICK ELEMENTS")) currentSection = "ELEMENT";
                            else if (line.Contains("INTERFACE TOPOLOGY")) currentSection = "INTERFACE";
                            else currentSection = "";
                            continue;
                        }

                        var parts = line.Split(null as char[], StringSplitOptions.RemoveEmptyEntries);
                        if (parts.Length == 0) continue;

                        if (currentSection == "NODE" && parts[0] == "Node")
                        {
                            int nodeID = int.Parse(parts[1]);
                            double x = double.Parse(parts[2], System.Globalization.CultureInfo.InvariantCulture);
                            double y = double.Parse(parts[3], System.Globalization.CultureInfo.InvariantCulture);
                            double z = double.Parse(parts[4], System.Globalization.CultureInfo.InvariantCulture);
                            allNodes[nodeID] = new Point3d(x, y, z);
                        }
                        else if (currentSection == "ELEMENT")
                        {
                            string elemType = parts[0];
                            int branchID = int.Parse(parts[3]);
                            var nodeIDs = new List<int>();
                            for (int i = 4; i < parts.Length; i++) nodeIDs.Add(int.Parse(parts[i]));
                            allElements.Add(new ElementData { ElementType = elemType, BranchID = branchID, NodeIDs = nodeIDs });
                        }
                    }
                }

                if (allNodes.Count == 0 || allElements.Count == 0)
                    throw new Exception("No nodes or elements were read.");

                foreach (var element in allElements)
                {
                    var meshPoints = new List<Point3d>();
                    foreach (int nodeID in element.NodeIDs)
                    {
                        if (allNodes.ContainsKey(nodeID)) meshPoints.Add(allNodes[nodeID]);
                        else throw new Exception($"Missing node ID: {nodeID}");
                    }

                    Mesh newMesh = null;
                    if (element.ElementType == "Hexa8" && meshPoints.Count == 8)
                        newMesh = CreateHexa8Mesh(meshPoints);
                    else if (element.ElementType == "Wedge6" && meshPoints.Count == 6)
                        newMesh = CreateWedge6Mesh(meshPoints);

                    if (newMesh != null)
                    {
                        newMesh.UnifyNormals();
                        newMesh.RebuildNormals();
                        GH_Path path = new GH_Path(element.BranchID - 1);
                        outputMeshes.Add(newMesh, path);
                    }
                }

                timer.Stop();
                R = outputMeshes;
                A = $"Successfully imported {outputMeshes.DataCount} meshes in {timer.Elapsed.TotalSeconds:0.00}s";
            }
            catch (Exception ex)
            {
                A = $"Error: {ex.Message}";
            }
        }

        // -------------------------
        // IMPORT LOGIC (3DEC .3ddat)
        // -------------------------
        public static void DoImport3DDAT(string path, out DataTree<Mesh> meshes, out string info)
        {
            meshes = new DataTree<Mesh>();
            info = string.Empty;

            if (string.IsNullOrEmpty(path) || !File.Exists(path))
            {
                info = "Error: File path is invalid or file does not exist.";
                return;
            }

            var timer = System.Diagnostics.Stopwatch.StartNew();

            try
            {
                // Read the entire file
                string fileContent = File.ReadAllText(path);

                // Step 1: Remove comments (lines starting with ;)
                string[] rawLines = fileContent.Split(new[] { '\r', '\n' }, StringSplitOptions.RemoveEmptyEntries);
                StringBuilder cleanContent = new StringBuilder();
                foreach (string line in rawLines)
                {
                    string trimmed = line.Trim();
                    if (!trimmed.StartsWith(";") && trimmed.ToLower() != "ret")
                    {
                        cleanContent.AppendLine(line);
                    }
                }
                fileContent = cleanContent.ToString();

                // Step 2: Remove line continuations (&) - join into single lines per block
                fileContent = Regex.Replace(fileContent, @"&\s*\r?\n\s*", " ");
                fileContent = Regex.Replace(fileContent, @"&", " "); // Remove any remaining &

                // Step 3: Split into lines (each line should now be a complete poly reg command)
                string[] lines = fileContent.Split(new[] { '\r', '\n' }, StringSplitOptions.RemoveEmptyEntries);

                Dictionary<int, List<List<Point3d>>> blockFaces = new Dictionary<int, List<List<Point3d>>>();
                int blocksFound = 0;

                foreach (string rawLine in lines)
                {
                    string line = rawLine.Trim();

                    if (string.IsNullOrWhiteSpace(line))
                        continue;

                    // Check if this line contains a poly reg command
                    Match polyMatch = Regex.Match(line, @"poly\s+(?:reg|region)\s+(\d+)", RegexOptions.IgnoreCase);
                    if (!polyMatch.Success)
                        continue;

                    int blockId = int.Parse(polyMatch.Groups[1].Value);
                    blocksFound++;

                    if (!blockFaces.ContainsKey(blockId))
                        blockFaces[blockId] = new List<List<Point3d>>();

                    // Find all face definitions in this line
                    MatchCollection faceMatches = Regex.Matches(line,
                        @"face\s+ID\s+\d+\s+((?:[+-]?\d+\.?\d*[eE][+-]?\d+\s*)+)",
                        RegexOptions.IgnoreCase);

                    foreach (Match faceMatch in faceMatches)
                    {
                        string coordString = faceMatch.Groups[1].Value;
                        List<Point3d> facePoints = Parse3DDATCoordinates(coordString);

                        if (facePoints.Count >= 3)
                        {
                            blockFaces[blockId].Add(facePoints);
                        }
                    }
                }

                // Create meshes from parsed blocks
                int meshCount = 0;
                foreach (var kvp in blockFaces.OrderBy(x => x.Key))
                {
                    int blockId = kvp.Key;
                    List<List<Point3d>> faces = kvp.Value;

                    // Need at least 4 faces for a closed solid
                    if (faces.Count < 4)
                    {
                        continue;
                    }

                    Mesh blockMesh = Create3DDATMeshFromFaces(faces);
                    if (blockMesh != null && blockMesh.Vertices.Count > 0 && blockMesh.Faces.Count > 0)
                    {
                        meshes.Add(blockMesh, new GH_Path(0));
                        meshCount++;
                    }
                }

                timer.Stop();

                if (meshCount == 0)
                {
                    info = $"Warning: Found {blocksFound} block definitions but could not create valid meshes. Check file format.";
                }
                else
                {
                    info = $"Successfully imported {meshCount} blocks from 3DEC file in {timer.Elapsed.TotalSeconds:0.00}s";
                }
            }
            catch (Exception ex)
            {
                info = $"Error: {ex.Message}";
            }
        }

        private static List<Point3d> Parse3DDATCoordinates(string coordString)
        {
            List<Point3d> points = new List<Point3d>();

            // Match scientific notation numbers
            MatchCollection matches = Regex.Matches(coordString,
                @"[+-]?\d+\.?\d*[eE][+-]?\d+|[+-]?\d+\.\d+|[+-]?\d+");

            List<double> values = new List<double>();
            foreach (Match m in matches)
            {
                if (double.TryParse(m.Value, System.Globalization.NumberStyles.Float,
                    System.Globalization.CultureInfo.InvariantCulture, out double val))
                {
                    values.Add(val);
                }
            }

            // Group into xyz triplets
            for (int i = 0; i + 2 <= values.Count; i += 3)
            {
                points.Add(new Point3d(values[i], values[i + 1], values[i + 2]));
            }

            return points;
        }

        private static Mesh Create3DDATMeshFromFaces(List<List<Point3d>> faces)
        {
            Mesh mesh = new Mesh();
            Dictionary<string, int> vertexMap = new Dictionary<string, int>();

            foreach (var face in faces)
            {
                if (face.Count < 3) continue;

                List<int> faceIndices = new List<int>();

                foreach (Point3d pt in face)
                {
                    // Create a key with reasonable precision
                    string key = $"{pt.X:F8}_{pt.Y:F8}_{pt.Z:F8}";

                    if (!vertexMap.TryGetValue(key, out int vertIndex))
                    {
                        vertIndex = mesh.Vertices.Count;
                        mesh.Vertices.Add(pt);
                        vertexMap[key] = vertIndex;
                    }

                    faceIndices.Add(vertIndex);
                }

                // Add the face - note: 3DEC faces have inward normals,
                // so we reverse to get outward normals for Rhino display
                if (faceIndices.Count == 3)
                {
                    mesh.Faces.AddFace(faceIndices[2], faceIndices[1], faceIndices[0]);
                }
                else if (faceIndices.Count == 4)
                {
                    mesh.Faces.AddFace(faceIndices[3], faceIndices[2], faceIndices[1], faceIndices[0]);
                }
                else if (faceIndices.Count > 4)
                {
                    // Triangulate polygon (reversed)
                    for (int i = 1; i < faceIndices.Count - 1; i++)
                    {
                        mesh.Faces.AddFace(faceIndices[0], faceIndices[i + 1], faceIndices[i]);
                    }
                }
            }

            if (mesh.Vertices.Count > 0 && mesh.Faces.Count > 0)
            {
                mesh.Vertices.CombineIdentical(true, true);
                mesh.Compact();
                mesh.UnifyNormals();
                mesh.RebuildNormals();
            }

            return mesh;
        }

        // -------------------------
        // EXPORT LOGIC (STANDARD .txt)
        // -------------------------
        public static void DoExport(DataTree<Mesh> M, string P, bool toStraus, double threshold, out string R)
        {
            R = string.Empty;
            try
            {
                var timer = System.Diagnostics.Stopwatch.StartNew();

                List<Mesh> allFlatMeshes = new List<Mesh>();
                List<int> meshBranchIndices = new List<int>();
                int bIndex = 0;
                foreach (var branch in M.Branches)
                {
                    bIndex++;
                    foreach (var mesh in branch)
                    {
                        if (mesh != null && mesh.IsValid)
                        {
                            allFlatMeshes.Add(mesh);
                            meshBranchIndices.Add(bIndex);
                        }
                    }
                }

                List<Point3d> nodeAccumulatorSum = new List<Point3d>();
                List<int> nodeAccumulatorCount = new List<int>();
                RTree rTree = new RTree();

                int approxSize = allFlatMeshes.Sum(m => m.Vertices.Count);
                nodeAccumulatorSum.Capacity = approxSize;
                nodeAccumulatorCount.Capacity = approxSize;

                for (int m = 0; m < allFlatMeshes.Count; m++)
                {
                    Mesh msh = allFlatMeshes[m];
                    for (int v = 0; v < msh.Vertices.Count; v++)
                    {
                        Point3d pt = msh.Vertices[v];
                        int foundIndex = -1;

                        rTree.Search(new Sphere(pt, threshold), (sender, args) =>
                        {
                            foundIndex = args.Id;
                            args.Cancel = true;
                        });

                        if (foundIndex != -1)
                        {
                            nodeAccumulatorSum[foundIndex] += pt;
                            nodeAccumulatorCount[foundIndex]++;
                        }
                        else
                        {
                            int newIndex = nodeAccumulatorSum.Count;
                            nodeAccumulatorSum.Add(pt);
                            nodeAccumulatorCount.Add(1);
                            rTree.Insert(pt, newIndex);
                        }
                    }
                }

                List<Point3d> averagedNodes = new List<Point3d>(nodeAccumulatorSum.Count);
                for (int i = 0; i < nodeAccumulatorSum.Count; i++)
                {
                    averagedNodes.Add(nodeAccumulatorSum[i] / nodeAccumulatorCount[i]);
                }

                List<string> elementTypes = new List<string>();
                List<List<int>> elementNodes = new List<List<int>>();
                int processedMeshes = 0;

                for (int m = 0; m < allFlatMeshes.Count; m++)
                {
                    Mesh msh = allFlatMeshes[m];
                    bool isWedge = msh.Faces.Count == 5;
                    string eType = isWedge ? "Wedge6" : "Hexa8";

                    List<Point3d> orderedPoints = BuildTrueMesh(msh, isWedge);
                    List<int> finalNodeIDs = new List<int>();

                    foreach (Point3d pt in orderedPoints)
                    {
                        int foundID = -1;
                        rTree.Search(new Sphere(pt, threshold * 1.01 + 1e-9), (sender, args) =>
                        {
                            foundID = args.Id + 1;
                            args.Cancel = true;
                        });

                        if (foundID != -1)
                            finalNodeIDs.Add(foundID);
                        else
                            throw new Exception("Topology Error: A mesh vertex could not be mapped to a global node.");
                    }

                    elementTypes.Add(eType);
                    elementNodes.Add(finalNodeIDs);
                    processedMeshes++;
                }

                List<string> interfaceTopology = new List<string>();
                if (!toStraus)
                {
                    for (int i = 0; i < elementNodes.Count; i++)
                    {
                        var facesI = GetElementFaces(elementNodes[i], elementTypes[i]);
                        for (int j = i + 1; j < elementNodes.Count; j++)
                        {
                            var facesJ = GetElementFaces(elementNodes[j], elementTypes[j]);
                            foreach (var faceI in facesI)
                            {
                                bool matchFound = false;
                                foreach (var faceJ in facesJ)
                                {
                                    if (AreFacesEqual(faceI, faceJ))
                                    {
                                        string interfaceLine = string.Format("{0,10}{1,10}", i + 1, j + 1);
                                        foreach (int nodeIndex in faceI)
                                            interfaceLine += string.Format("{0,10}", nodeIndex);
                                        interfaceTopology.Add(interfaceLine);
                                        matchFound = true;
                                        break;
                                    }
                                }
                                if (matchFound) break;
                            }
                        }
                    }
                }

                using (StreamWriter writer = new StreamWriter(P))
                {
                    writer.WriteLine("/ ______________________________________________________________________________");
                    writer.WriteLine("/ NODE COORDINATES");
                    writer.WriteLine();

                    for (int i = 0; i < averagedNodes.Count; i++)
                    {
                        Point3d pt = averagedNodes[i];
                        writer.WriteLine(" Node {0} {1} {2} {3}", i + 1,
                            FormatCoordinate(pt.X), FormatCoordinate(pt.Y), FormatCoordinate(pt.Z));
                    }

                    writer.WriteLine();
                    writer.WriteLine("/ ______________________________________________________________________________");
                    writer.WriteLine("/ BRICK ELEMENTS");
                    writer.WriteLine();

                    for (int i = 0; i < elementNodes.Count; i++)
                    {
                        writer.Write("{0,-15} {1} 1 {2}", elementTypes[i], i + 1, meshBranchIndices[i]);
                        foreach (int nodeID in elementNodes[i])
                            writer.Write(" {0}", nodeID);
                        writer.WriteLine();
                    }

                    if (!toStraus && interfaceTopology.Count > 0)
                    {
                        writer.WriteLine();
                        writer.WriteLine("/ ______________________________________________________________________________");
                        writer.WriteLine("/ INTERFACE TOPOLOGY");
                        writer.WriteLine();
                        foreach (string line in interfaceTopology)
                            writer.WriteLine(line);
                    }
                }

                timer.Stop();
                R = $"Exported {processedMeshes} meshes. Nodes: {averagedNodes.Count}. Time: {timer.Elapsed.TotalSeconds:0.00}s";
            }
            catch (Exception ex)
            {
                R = $"Error: {ex.Message}";
            }
        }

        // -------------------------
        // EXPORT LOGIC (3DEC .3ddat)
        // -------------------------
        public static void DoExport3DDAT(DataTree<Mesh> meshes, string basePath, out string info)
        {
            info = string.Empty;

            try
            {
                var timer = System.Diagnostics.Stopwatch.StartNew();

                if (meshes == null || meshes.BranchCount == 0)
                {
                    info = "Error: No meshes to export.";
                    return;
                }

                string directory = Path.GetDirectoryName(basePath);
                string fileNameWithoutExt = Path.GetFileNameWithoutExtension(basePath);
                string extension = Path.GetExtension(basePath);

                if (string.IsNullOrEmpty(extension) || extension.ToLower() != ".3ddat")
                    extension = ".3ddat";

                int totalBlocks = 0;
                List<string> createdFiles = new List<string>();

                const double planarityTolerance = 1e-6;

                for (int branchIndex = 0; branchIndex < meshes.BranchCount; branchIndex++)
                {
                    GH_Path path = meshes.Paths[branchIndex];
                    List<Mesh> branchMeshes = meshes.Branch(path);

                    if (branchMeshes == null || branchMeshes.Count == 0)
                        continue;

                    string branchFilePath;
                    if (meshes.BranchCount == 1)
                    {
                        branchFilePath = Path.Combine(directory, $"{fileNameWithoutExt}{extension}");
                    }
                    else
                    {
                        branchFilePath = Path.Combine(directory, $"{fileNameWithoutExt}_{branchIndex}{extension}");
                    }

                    StringBuilder sb = new StringBuilder();
                    sb.AppendLine(";");
                    sb.AppendLine("; 3DEC input deck produced by PoliBrick");
                    sb.AppendLine($"; geometry built: {DateTime.Now}");
                    sb.AppendLine(";");

                    int blockId = 1;
                    int globalFaceId = 1;

                    foreach (Mesh mesh in branchMeshes)
                    {
                        if (mesh == null || !mesh.IsValid) continue;

                        Mesh workMesh = mesh.DuplicateMesh();
                        workMesh.Vertices.CombineIdentical(true, true);
                        workMesh.UnifyNormals();
                        workMesh.RebuildNormals();

                        Point3d centroid = Point3d.Origin;
                        foreach (Point3f v in workMesh.Vertices)
                            centroid += new Point3d(v);
                        centroid /= workMesh.Vertices.Count;

                        List<List<Point3d>> allFaces = new List<List<Point3d>>();

                        for (int f = 0; f < workMesh.Faces.Count; f++)
                        {
                            MeshFace face = workMesh.Faces[f];

                            Point3d pA = new Point3d(workMesh.Vertices[face.A]);
                            Point3d pB = new Point3d(workMesh.Vertices[face.B]);
                            Point3d pC = new Point3d(workMesh.Vertices[face.C]);

                            if (face.IsQuad)
                            {
                                Point3d pD = new Point3d(workMesh.Vertices[face.D]);

                                if (IsQuadPlanar(pA, pB, pC, pD, planarityTolerance))
                                {
                                    allFaces.Add(new List<Point3d> { pA, pB, pC, pD });
                                }
                                else
                                {
                                    allFaces.Add(new List<Point3d> { pA, pB, pC });
                                    allFaces.Add(new List<Point3d> { pA, pC, pD });
                                }
                            }
                            else
                            {
                                allFaces.Add(new List<Point3d> { pA, pB, pC });
                            }
                        }

                        sb.AppendLine($"; block {blockId}");
                        sb.AppendLine($"poly reg {blockId,9} mat   1 con   1 &");

                        foreach (var pts in allFaces)
                        {
                            Point3d faceCent = Point3d.Origin;
                            foreach (Point3d pt in pts)
                                faceCent += pt;
                            faceCent /= pts.Count;

                            Vector3d normal = Vector3d.CrossProduct(
                                pts[1] - pts[0],
                                pts[2] - pts[0]);

                            List<Point3d> orderedPts = new List<Point3d>(pts);

                            // 3DEC requires normals pointing INWARD
                            if (normal * (faceCent - centroid) > 0)
                                orderedPts.Reverse();

                            sb.Append($"face ID {globalFaceId,9}  ");
                            sb.AppendLine($"{orderedPts[0].X:e} {orderedPts[0].Y:e} {orderedPts[0].Z:e} &");

                            for (int i = 1; i < orderedPts.Count; i++)
                            {
                                sb.Append("                   ");
                                sb.AppendLine($"{orderedPts[i].X:e} {orderedPts[i].Y:e} {orderedPts[i].Z:e} &");
                            }

                            globalFaceId++;
                        }

                        sb.AppendLine(" ");
                        blockId++;
                        totalBlocks++;
                    }

                    sb.AppendLine(" ");
                    sb.AppendLine("ret");

                    File.WriteAllText(branchFilePath, sb.ToString());
                    createdFiles.Add(Path.GetFileName(branchFilePath));
                }

                timer.Stop();

                if (createdFiles.Count == 1)
                {
                    info = $"Exported {totalBlocks} blocks to {createdFiles[0]} in {timer.Elapsed.TotalSeconds:0.00}s";
                }
                else
                {
                    info = $"Exported {totalBlocks} blocks to {createdFiles.Count} files: {string.Join(", ", createdFiles)} in {timer.Elapsed.TotalSeconds:0.00}s";
                }
            }
            catch (Exception ex)
            {
                info = $"Error: {ex.Message}";
            }
        }

        // -------------------------
        // HELPER METHODS
        // -------------------------

        private static bool IsQuadPlanar(Point3d p1, Point3d p2, Point3d p3, Point3d p4, double tolerance)
        {
            Vector3d v1 = p2 - p1;
            Vector3d v2 = p3 - p1;
            Vector3d normal = Vector3d.CrossProduct(v1, v2);

            if (normal.IsTiny(1e-10))
                return false;

            normal.Unitize();

            Vector3d v4 = p4 - p1;
            double distance = Math.Abs(v4 * normal);

            return distance <= tolerance;
        }

        private static Mesh CreateHexa8Mesh(List<Point3d> n)
        {
            Mesh m = new Mesh();
            m.Vertices.AddVertices(n);
            m.Faces.AddFace(0, 1, 2, 3);
            m.Faces.AddFace(4, 5, 6, 7);
            m.Faces.AddFace(0, 3, 7, 4);
            m.Faces.AddFace(1, 2, 6, 5);
            m.Faces.AddFace(0, 4, 5, 1);
            m.Faces.AddFace(3, 2, 6, 7);
            return m;
        }

        private static Mesh CreateWedge6Mesh(List<Point3d> n)
        {
            Mesh m = new Mesh();
            m.Vertices.AddVertices(n);
            m.Faces.AddFace(0, 1, 2);
            m.Faces.AddFace(3, 4, 5);
            m.Faces.AddFace(0, 2, 5, 3);
            m.Faces.AddFace(1, 0, 3, 4);
            m.Faces.AddFace(2, 1, 4, 5);
            return m;
        }

        private static List<List<int>> GetElementFaces(List<int> nodes, string elementType)
        {
            if (elementType == "Hexa8" && nodes.Count == 8) return GetHexa8Faces(nodes);
            if (elementType == "Wedge6" && nodes.Count == 6) return GetWedge6Faces(nodes);
            return new List<List<int>>();
        }

        private static List<List<int>> GetWedge6Faces(List<int> n)
        {
            return new List<List<int>> {
                new List<int> { n[0], n[1], n[2] },
                new List<int> { n[3], n[4], n[5] },
                new List<int> { n[0], n[2], n[5], n[3] },
                new List<int> { n[1], n[0], n[3], n[4] },
                new List<int> { n[2], n[1], n[4], n[5] }
            };
        }

        private static List<List<int>> GetHexa8Faces(List<int> n)
        {
            return new List<List<int>> {
                new List<int> { n[0], n[1], n[2], n[3] },
                new List<int> { n[4], n[5], n[6], n[7] },
                new List<int> { n[0], n[3], n[7], n[4] },
                new List<int> { n[1], n[2], n[6], n[5] },
                new List<int> { n[0], n[4], n[5], n[1] },
                new List<int> { n[3], n[2], n[6], n[7] }
            };
        }

        private static bool AreFacesEqual(List<int> face1, List<int> face2)
        {
            if (face1.Count != face2.Count || face1.Count < 3) return false;
            var face1Set = new HashSet<int>(face1);
            return face2.All(node => face1Set.Contains(node));
        }

        private static List<Point3d> BuildTrueMesh(Mesh m, bool isWedge)
        {
            Mesh result = m.DuplicateMesh();
            result.Vertices.CombineIdentical(true, true);
            result.Vertices.CullUnused();
            result.UnifyNormals();
            result.RebuildNormals();
            FlipCheck(result);

            List<Point3d> faceAverts = new List<Point3d>(8);
            List<int> inds = new List<int>(4);
            MeshFace firstFace = isWedge ? result.Faces.First(f => !f.IsQuad) : result.Faces.First(f => f.IsQuad);

            faceAverts.Add(result.Vertices[firstFace.A]);
            if (firstFace.IsQuad) faceAverts.Add(result.Vertices[firstFace.D]);
            faceAverts.Add(result.Vertices[firstFace.C]);
            faceAverts.Add(result.Vertices[firstFace.B]);

            inds.Add(firstFace.A);
            if (firstFace.IsQuad) inds.Add(firstFace.D);
            inds.Add(firstFace.C);
            inds.Add(firstFace.B);

            HashSet<Point3d> uniquePoints = new HashSet<Point3d>(faceAverts);
            for (int j = 0; j < inds.Count; j++)
            {
                var connectedVerts = result.TopologyVertices.ConnectedTopologyVertices(inds[j], false);
                foreach (int connectedIndex in connectedVerts)
                {
                    var targetPoint = result.Vertices[connectedIndex];
                    if (uniquePoints.Add(targetPoint))
                    {
                        faceAverts.Add(targetPoint);
                        if ((isWedge && faceAverts.Count == 6) || (!isWedge && faceAverts.Count == 8)) break;
                    }
                }
                if ((isWedge && faceAverts.Count == 6) || (!isWedge && faceAverts.Count == 8)) break;
            }
            return faceAverts;
        }

        private static void FlipCheck(Mesh mesh)
        {
            var amp = AreaMassProperties.Compute(mesh);
            if (amp == null) return;
            Vector3d normal = (Vector3d)mesh.Normals[0];
            Vector3d toCentroid = amp.Centroid - (Point3d)mesh.Vertices[0];
            if (Vector3d.Multiply(normal, toCentroid) > 0) mesh.Flip(true, true, true);
        }

        private static string FormatCoordinate(double value)
        {
            return value.ToString("0.00000000000E+00", System.Globalization.CultureInfo.InvariantCulture);
        }
    }
}