using Grasshopper;
using Grasshopper.GUI;
using Grasshopper.GUI.Canvas;
using Grasshopper.Kernel;
using Grasshopper.Kernel.Attributes;
using Grasshopper.Kernel.Data;
using Grasshopper.Kernel.Types;
using Rhino.Geometry;
using System;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.Drawing;
using System.IO;
using System.Threading.Tasks;
using System.Collections;
using System.Drawing.Drawing2D;
using System.Linq;
using System.Windows.Forms;

namespace PoliBrick.Util
{
    // --- 1. Custom Attributes (UI) ---
    public class QuadMeshAttributes : GH_ComponentAttributes
    {
        private QuadConvertMeshSolver _owner;

        public QuadMeshAttributes(QuadConvertMeshSolver owner) : base(owner)
        {
            _owner = owner;
        }

        protected override void Layout()
        {
            base.Layout();
            // Add space at bottom for the button
            Bounds = new RectangleF(Bounds.X, Bounds.Y, Bounds.Width, Bounds.Height + 25);
        }

        protected override void Render(GH_Canvas canvas, Graphics graphics, GH_CanvasChannel channel)
        {
            base.Render(canvas, graphics, channel);

            if (channel == GH_CanvasChannel.Objects)
            {
                // 1. Define Geometry
                int buttonHeight = 18; // Standard height from previous style
                int gap = 3;

                // Calculate integer coordinates for sharper drawing
                Rectangle r = new Rectangle(
                    (int)Bounds.X + gap,
                    (int)Bounds.Bottom - buttonHeight - gap,
                    (int)Bounds.Width - (gap * 2),
                    buttonHeight
                );

                // 2. Define Colors (Matches Previous Style)
                bool active = _owner.FlipNormals;
                Color colorOff = Color.FromArgb(225, 225, 225); // Light Gray
                Color colorOn = Color.Orange;
                Color colorBorder = Color.Gray;
                Color textOff = Color.DimGray;
                Color textOn = Color.White;

                // 3. Draw Rounded Background
                using (GraphicsPath path = RoundedRect(r, r.Height / 4))
                {
                    using (Brush b = new SolidBrush(active ? colorOn : colorOff))
                    {
                        graphics.FillPath(b, path);
                    }
                    using (Pen p = new Pen(colorBorder))
                    {
                        graphics.DrawPath(p, path);
                    }
                }

                // 4. Draw Text
                using (Font font = new Font("Arial", 7f, FontStyle.Bold)) // Font from previous style
                using (Brush textBrush = new SolidBrush(active ? textOn : textOff))
                {
                    StringFormat fmt = new StringFormat
                    {
                        Alignment = StringAlignment.Center,
                        LineAlignment = StringAlignment.Center
                    };
                    graphics.DrawString("Flip", font, textBrush, r, fmt);
                }
            }
        }

        public override GH_ObjectResponse RespondToMouseDown(GH_Canvas sender, GH_CanvasMouseEvent e)
        {
            if (e.Button == MouseButtons.Left)
            {
                // Re-calculate the specific button rectangle to ensure accurate clicking
                int buttonHeight = 18;
                int gap = 3;
                RectangleF buttonRect = new RectangleF(
                    Bounds.X + gap,
                    Bounds.Bottom - buttonHeight - gap,
                    Bounds.Width - (gap * 2),
                    buttonHeight
                );

                if (buttonRect.Contains(e.CanvasLocation))
                {
                    _owner.FlipNormals = !_owner.FlipNormals;
                    _owner.ExpireSolution(true);
                    return GH_ObjectResponse.Handled;
                }
            }
            return base.RespondToMouseDown(sender, e);
        }

        // Helper method for the specific rounded look
        
        private GraphicsPath RoundedRect(Rectangle r, int rad)
        {
            GraphicsPath p = new GraphicsPath();
            int d = rad * 2;
            p.AddArc(r.X, r.Y, d, d, 180, 90);
            p.AddArc(r.Right - d, r.Y, d, d, 270, 90);
            p.AddArc(r.Right - d, r.Bottom - d, d, d, 0, 90);
            p.AddArc(r.X, r.Bottom - d, d, d, 90, 90);
            p.CloseFigure();
            return p;
        }
    }

    // --- 2. Main Solver Class ---
    public class QuadConvertMeshSolver : GH_Component
    {
        public bool FlipNormals { get; set; } = false;
        private class FaceMetadata
        {
            public int FaceIndex;
            public int StripId;
            public int StripType;         // 0=main, 1=helper, 2=cap
            public int CreationOrder;     // Global insertion order
            public double ArcParameter;   // 0-1 position along strip
            public Point3d Center;
            public int LocalIndex;        // Position within strip (after sorting)

            // Assigned during pattern logic
            public GH_Path FinalPath;
            public char PatternChar;
            // ===== NEW: Store boundary curve start points =====
            public Point3d BoundaryStartA;
            public Point3d BoundaryStartB;
        }

        private List<FaceMetadata> faceMetadataList = new List<FaceMetadata>();
        public QuadConvertMeshSolver()
         : base("2.Propagating Quad Mesh", "PropQuad",
           "Creates a quad-priority mesh with Pattern Cycling, Global Colors, and Consistent Ordering",
           "PoliBrick", "Utilities")
        {
        }
        public override GH_Exposure Exposure
        {
            get { return GH_Exposure.secondary; }
        }
        public override Guid ComponentGuid => new Guid("9A3E8D5F-1165-4F8A-A0D3-97B2A91D3C01");

        public override void CreateAttributes()
        {
            m_attributes = new QuadMeshAttributes(this);
        }

        public override bool Write(GH_IO.Serialization.GH_IWriter writer)
        {
            writer.SetBoolean("FlipNormals", FlipNormals);
            return base.Write(writer);
        }

        public override bool Read(GH_IO.Serialization.GH_IReader reader)
        {
            if (reader.ItemExists("FlipNormals"))
                FlipNormals = reader.GetBoolean("FlipNormals");
            return base.Read(reader);
        }

        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddCurveParameter("Meshing Curves", "C", "Curves used for meshing (already refined/ordered). Tree must be in strip order.", GH_ParamAccess.tree);
            pManager.AddTextParameter("Patterns", "P", "Pattern for block sub-branching (e.g., 'AAB'). Will cycle through rows.", GH_ParamAccess.tree);
            pManager.AddNumberParameter("Length", "L", "Target division length", GH_ParamAccess.item);
            pManager.AddNumberParameter("Thickness", "T", "Extrusion thickness for discrete blocks", GH_ParamAccess.item, 0.0);
            pManager.AddNumberParameter("Angle Tolerance", "A", "Angle tolerance for quad conversion (degrees)", GH_ParamAccess.item, 3.0);
            pManager.AddNumberParameter("Aspect Ratio", "AR", "Aspect ratio tolerance for quad conversion", GH_ParamAccess.item, 0.075);
            pManager.AddNumberParameter("Minimum Corner Angle", "MinA", "Minimum corner angle tolerance for quad conversion (degrees)", GH_ParamAccess.item, 10.0);
            pManager.AddSurfaceParameter("Surface", "S", "Base surface for projections", GH_ParamAccess.item);
        }

        protected override void RegisterOutputParams(GH_OutputParamManager pManager)
        {
            pManager.AddMeshParameter("Homogenized Mesh", "M", "The unified mesh, built from all branches", GH_ParamAccess.item);
            pManager.AddMeshParameter("Discrete Mesh Blocks", "B", "DataTree of discrete mesh blocks, branched by {Row;Sub}", GH_ParamAccess.tree);
        }

        protected override void SolveInstance(IGH_DataAccess DA)
        {
            GH_Structure<GH_Curve> curvesTree;
            GH_Structure<GH_String> patternsTree;
            double length = 0.0;
            double thickness = 0.0;
            double angleTolDegrees = 3.0;
            double aspectRatio = 0.075;
            double minAngDegrees = 10.0;
            Surface surface = null;

            if (!DA.GetDataTree(0, out curvesTree)) return;
            if (!DA.GetDataTree(1, out patternsTree)) return;
            if (!DA.GetData(2, ref length)) return;
            if (!DA.GetData(3, ref thickness)) return;
            DA.GetData(4, ref angleTolDegrees);
            DA.GetData(5, ref aspectRatio);
            DA.GetData(6, ref minAngDegrees);
            if (!DA.GetData(7, ref surface)) return;

            if (length <= 0.0)
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Warning, "Length must be positive.");
                return;
            }

            // --- 2. MESHING LOGIC ---
            Brep brep = Brep.CreateFromSurface(surface);
            double pullTol = 0.001;
            Mesh masterMesh = new Mesh();
            DataTree<Mesh> blocksTree = new DataTree<Mesh>();
            var vertexIndexToRow = new Dictionary<int, int>();
            var pointToVertexIndex = new Dictionary<Point3d, int>(new Point3dComparer());
            var vertexToArcLength = new Dictionary<int, double>();
            var vertexToUV = new Dictionary<int, Point2d>();
            var branchToCurves = new Dictionary<int, List<Curve>>();

            faceMetadataList.Clear(); // Clear from previous runs

            IList<GH_Path> paths = curvesTree.Paths;
            if (paths.Count < 2) return;

            // ===== PASS 1: Sample branches intelligently =====
            var allBranchVertexIndices = new Dictionary<int, List<List<int>>>();
            var sub0VerticesByBase = new Dictionary<int, List<List<int>>>();

            // First pass: identify which branches have {n;0}
            var basesWithSub0 = new HashSet<int>();
            for (int i = 0; i < paths.Count; i++)
            {
                GH_Path pth = paths[i];
                int baseIndex = pth.Indices[0];
                int subIndex = (pth.Length > 1) ? pth.Indices[1] : -1;
                if (subIndex == 0) basesWithSub0.Add(baseIndex);
            }

            // Second pass: sample branches (skip {n+1} if {n;0} exists)
            for (int i = 0; i < paths.Count; i++)
            {
                IList branchCurves = curvesTree.get_Branch(paths[i]);
                if (branchCurves == null) continue;

                GH_Path pth = paths[i];
                int baseIndex = pth.Indices[0];
                int subIndex = (pth.Length > 1) ? pth.Indices[1] : -1;
                bool isSub0 = (subIndex == 0);
                bool isMainWithSub0 = (subIndex == -1 && basesWithSub0.Contains(baseIndex));

                // SKIP SAMPLING if this is {n+1} and {n;0} exists - we'll build it from {n;0}
                bool shouldSkipSampling = false;
                if (subIndex == -1 && baseIndex > 0)
                {
                    // Check if previous base has {n;0}
                    if (basesWithSub0.Contains(baseIndex - 1))
                    {
                        shouldSkipSampling = true;
                    }
                }

                List<Curve> thisBranchCurves = new List<Curve>();
                var thisBranchIndices = new List<List<int>>();

                for (int jj = 0; jj < branchCurves.Count; jj++)
                {
                    GH_Curve ghCurve = branchCurves[jj] as GH_Curve;
                    if (ghCurve == null) continue;
                    Curve curve = ghCurve.Value;
                    if (curve == null) continue;
                    thisBranchCurves.Add(curve);

                    if (shouldSkipSampling)
                    {
                        
                        thisBranchIndices.Add(new List<int>());
                        continue;
                    }

                    // Normal sampling for all other branches
                    List<Point3d> points = GetProjectedPointsFromCurve(curve, length, surface, pullTol);
                    if (points == null || points.Count < 2) continue;

                    double curveCumArc = 0.0;
                    List<double> arcLengths = new List<double>();
                    arcLengths.Add(curveCumArc);
                    for (int k = 1; k < points.Count; k++)
                    {
                        curveCumArc += points[k - 1].DistanceTo(points[k]);
                        arcLengths.Add(curveCumArc);
                    }

                    List<int> curveIndices = new List<int>();
                    for (int k = 0; k < points.Count; k++)
                    {
                        Point3d pt = points[k];
                        int vIndex = AddVertex(masterMesh, pt, pointToVertexIndex);
                        curveIndices.Add(vIndex);
                        if (!vertexIndexToRow.ContainsKey(vIndex)) vertexIndexToRow.Add(vIndex, i);
                        if (!vertexToArcLength.ContainsKey(vIndex)) vertexToArcLength.Add(vIndex, arcLengths[k]);

                        if (!vertexToUV.ContainsKey(vIndex))
                        {
                            surface.ClosestPoint(pt, out double u, out double v);
                            vertexToUV[vIndex] = new Point2d(u, v);
                        }
                    }
                    thisBranchIndices.Add(curveIndices);

                    if (isSub0)
                    {
                        if (!sub0VerticesByBase.ContainsKey(baseIndex))
                            sub0VerticesByBase[baseIndex] = new List<List<int>>();
                        sub0VerticesByBase[baseIndex].Add(new List<int>(curveIndices));
                    }
                }

                branchToCurves.Add(i, thisBranchCurves);
                allBranchVertexIndices.Add(i, thisBranchIndices);
            }

            // ===== PASS 2: Build {n+1} branches from {n;0} vertices =====
            for (int i = 0; i < paths.Count; i++)
            {
                GH_Path pth = paths[i];
                int baseIndex = pth.Indices[0];
                int subIndex = (pth.Length > 1) ? pth.Indices[1] : -1;
                bool isMain = (subIndex == -1);

                if (!isMain || baseIndex == 0) continue;
                if (!sub0VerticesByBase.TryGetValue(baseIndex - 1, out List<List<int>> sub0CurvesList)) continue;
                if (!allBranchVertexIndices.TryGetValue(i, out var thisBranchIndices)) continue;
                if (!branchToCurves.TryGetValue(i, out var thisBranchCurves)) continue;

                int numSub0Curves = sub0CurvesList.Count;
                int numN1Curves = thisBranchIndices.Count;
                double maxDist = Math.Max(length * 0.75, (Rhino.RhinoDoc.ActiveDoc?.ModelAbsoluteTolerance ?? 1e-3) * 20.0);

                if (numSub0Curves == numN1Curves)
                {
                    // TOPOLOGY 1: Corner Connection - Direct replacement
                    for (int jj = 0; jj < thisBranchIndices.Count; jj++)
                    {
                        List<int> currentIndices = thisBranchIndices[jj];
                        if (jj >= thisBranchCurves.Count || jj >= sub0CurvesList.Count) continue;

                        List<int> sub0Vertices = sub0CurvesList[jj];
                        if (sub0Vertices == null || sub0Vertices.Count < 2) continue;

                        Curve curve = thisBranchCurves[jj];

                        // Filter {n;0} vertices to only those ON this {n+1} curve
                        var verticesOnCurve = new List<(int vIdx, double t)>();
                        foreach (int vIdx in sub0Vertices)
                        {
                            Point3d pt = (Point3d)masterMesh.Vertices[vIdx];
                            if (curve.ClosestPoint(pt, out double t))
                            {
                                Point3d closestPt = curve.PointAt(t);
                                if (pt.DistanceTo(closestPt) <= maxDist)
                                {
                                    verticesOnCurve.Add((vIdx, t));
                                }
                            }
                        }

                        // Sort by parameter and use as the vertex list
                        var sorted = verticesOnCurve.OrderBy(x => x.t).Select(x => x.vIdx).ToList();

                        currentIndices.Clear();
                        currentIndices.AddRange(sorted);

                        foreach (int vIdx in sorted)
                        {
                            if (!vertexIndexToRow.ContainsKey(vIdx))
                                vertexIndexToRow[vIdx] = i;
                        }
                    }
                }
                else if (numSub0Curves < numN1Curves)
                {
                    // TOPOLOGY 2: Border Connection
                    var n1ToSub0Mapping = new Dictionary<int, int>();

                    // Map each {n+1} curve to best matching {n;0} curve
                    for (int jj = 0; jj < thisBranchCurves.Count; jj++)
                    {
                        Curve n1Curve = thisBranchCurves[jj];
                        int bestSub0Index = -1;
                        int maxMatches = 0;

                        for (int s0 = 0; s0 < sub0CurvesList.Count; s0++)
                        {
                            List<int> sub0Verts = sub0CurvesList[s0];
                            if (sub0Verts == null || sub0Verts.Count == 0) continue;
                            int matchCount = 0;

                            foreach (int vIdx in sub0Verts)
                            {
                                Point3d pt = (Point3d)masterMesh.Vertices[vIdx];
                                if (n1Curve.ClosestPoint(pt, out double t))
                                {
                                    if (pt.DistanceTo(n1Curve.PointAt(t)) <= maxDist)
                                        matchCount++;
                                }
                            }
                            if (matchCount > maxMatches)
                            {
                                maxMatches = matchCount;
                                bestSub0Index = s0;
                            }
                        }
                        n1ToSub0Mapping[jj] = bestSub0Index;
                    }

                    // Build vertex lists from {n;0}
                    for (int jj = 0; jj < thisBranchIndices.Count; jj++)
                    {
                        List<int> currentIndices = thisBranchIndices[jj];
                        if (!n1ToSub0Mapping.TryGetValue(jj, out int sub0Index) || sub0Index < 0) continue;
                        if (sub0Index >= sub0CurvesList.Count) continue;

                        List<int> sub0Vertices = sub0CurvesList[sub0Index];
                        if (sub0Vertices == null || sub0Vertices.Count < 2) continue;
                        Curve curve = thisBranchCurves[jj];

                        // Filter to vertices ON this curve
                        var verticesOnCurve = new List<(int vIdx, double t)>();
                        foreach (int vIdx in sub0Vertices)
                        {
                            Point3d pt = (Point3d)masterMesh.Vertices[vIdx];
                            if (curve.ClosestPoint(pt, out double t))
                            {
                                Point3d closestPt = curve.PointAt(t);
                                if (pt.DistanceTo(closestPt) <= maxDist)
                                {
                                    verticesOnCurve.Add((vIdx, t));
                                }
                            }
                        }

                        // Sort and assign
                        var sorted = verticesOnCurve.OrderBy(x => x.t).Select(x => x.vIdx).ToList();

                        currentIndices.Clear();
                        currentIndices.AddRange(sorted);

                        foreach (int vIdx in sorted)
                        {
                            if (!vertexIndexToRow.ContainsKey(vIdx))
                                vertexIndexToRow[vIdx] = i;
                        }
                    }
                }
            }

            // --- EXECUTE MESHING AND CAPS ---
            var pairs = BuildMeshingPairs(paths);
            // create strip ids in the same order
            var pairToStripId = new Dictionary<(int, int), int>();
            for (int s = 0; s < pairs.Count; s++)
            {
                var a = pairs[s].idxA;
                var b = pairs[s].idxB;
                var key = a < b ? (a, b) : (b, a);
                pairToStripId[key] = s;
            }

            // Pass faceRowIds to track created faces
            ExecuteMeshingLogic(pairs, allBranchVertexIndices, branchToCurves, masterMesh,
            vertexIndexToRow, vertexToArcLength, pointToVertexIndex, vertexToUV,
            surface, pullTol, out Dictionary<int, List<Curve>> stripCenters);

            var capPairs = BuildCapPairs_Sub1(paths);

            // Pass faceRowIds to track cap faces
            BuildEndpointCaps_Sub1(capPairs, allBranchVertexIndices, masterMesh, stripCenters);

            // Pass faceRowIds to track loop closure faces
            CloseSmallNakedLoops(masterMesh, surface, length * 6.0);

            double angleTolRad = Math.PI / 180.0 * angleTolDegrees;
            double minCornerRad = Math.PI / 180.0 * minAngDegrees;
            FixBoundaryVertices(masterMesh, surface, vertexToUV);

            // Pass faceRowIds to maintain sync during quad conversion
            OptimizedConvertTrianglesToQuads(masterMesh, angleTolRad, aspectRatio, minCornerRad);

            masterMesh.Vertices.CombineIdentical(true, true);
            masterMesh.Faces.CullDegenerateFaces(); // This might shift indices slightly, but usually safe if limited to degenerate
            masterMesh.UnifyNormals();
            masterMesh.Normals.ComputeNormals();
            masterMesh.Compact();

            EnsureConsistentNormals(masterMesh, vertexIndexToRow, branchToCurves);

            if (thickness != 0.0)
            {
                BuildDiscreteBlocksWithCorrectOrder(masterMesh, blocksTree, thickness,
                    FlipNormals, patternsTree, stripCenters);
            }

            DA.SetData(0, masterMesh);
            DA.SetDataTree(1, blocksTree);

        }

        // --- 3. GLOBAL CORNER DETECTION & PROCESSING ---

        class CurveId
        {
            public GH_Path Path;
            public int Index;
            public Curve CurveGeometry;
            public CurveId(GH_Path p, int i, Curve c) { Path = p; Index = i; CurveGeometry = c; }
        }

        // --- HELPER METHODS ---
        // Helper to find neighbors using standard RhinoCommon Topology
        private List<FaceMetadata> SortFacesAlongStrip(
    Mesh mesh,
    List<FaceMetadata> facesInStrip,
    Curve centerCurve)
        {
            if (facesInStrip.Count <= 1) return facesInStrip;

            // Build adjacency graph
            var faceSet = new HashSet<int>(facesInStrip.Select(f => f.FaceIndex));
            var adjacency = new Dictionary<int, List<int>>();
            var topologyEdges = mesh.TopologyEdges;

            foreach (var faceMeta in facesInStrip)
            {
                int fIdx = faceMeta.FaceIndex;
                adjacency[fIdx] = new List<int>();

                int[] faceEdges = topologyEdges.GetEdgesForFace(fIdx);
                foreach (int edgeIdx in faceEdges)
                {
                    int[] connectedFaces = topologyEdges.GetConnectedFaces(edgeIdx);
                    foreach (int neighbor in connectedFaces)
                    {
                        if (neighbor != fIdx && faceSet.Contains(neighbor))
                        {
                            if (!adjacency[fIdx].Contains(neighbor))
                                adjacency[fIdx].Add(neighbor);
                        }
                    }
                }
            }

            // ===== FIND START POINT USING STORED BOUNDARY INFO =====
            int startFace = -1;

            // Get boundary start points from first face (all faces in strip share same boundaries)
            Point3d boundaryStartA = facesInStrip.First().BoundaryStartA;
            Point3d boundaryStartB = facesInStrip.First().BoundaryStartB;

            // Check if boundary info is valid
            bool hasBoundaryInfo = (boundaryStartA != Point3d.Origin || boundaryStartB != Point3d.Origin);

            if (hasBoundaryInfo)
            {
                // Calculate strip start as midpoint of both boundary starts
                Point3d stripStart = (boundaryStartA + boundaryStartB) * 0.5;

                // Find endpoint faces (faces with 1-2 neighbors = strip ends)
                var endpointFaces = adjacency
                    .Where(kvp => kvp.Value.Count <= 2)
                    .Select(kvp => facesInStrip.First(f => f.FaceIndex == kvp.Key))
                    .ToList();

                if (endpointFaces.Count >= 2)
                {
                    // Analyze each endpoint
                    var endpointAnalysis = endpointFaces
                        .Select(f => new
                        {
                            Face = f,
                            DistToStripStart = f.Center.DistanceTo(stripStart),
                            DistToStartA = f.Center.DistanceTo(boundaryStartA),
                            DistToStartB = f.Center.DistanceTo(boundaryStartB),
                            Arc = f.ArcParameter
                        })
                        .OrderBy(x => x.DistToStripStart)
                        .ToList();

                    var closestToStart = endpointAnalysis.First();
                    var farthestFromStart = endpointAnalysis.Last();

                    // Calculate strip length for threshold
                    Point3d boundaryEndA = Point3d.Origin;
                    Point3d boundaryEndB = Point3d.Origin;

                    // Estimate strip end (use face with highest arc parameter)
                    var endFace = facesInStrip.OrderByDescending(f => f.ArcParameter).First();
                    double stripLength = stripStart.DistanceTo(endFace.Center);
                    double overhangThreshold = stripLength * 0.4;

                    // Check if farthest endpoint is actually near the start (overhang case)
                    bool farthestIsNearStart = (farthestFromStart.DistToStripStart < overhangThreshold);

                    if (farthestIsNearStart && farthestFromStart.Arc < 0.35)
                    {
                        // FARTHEST endpoint is near strip start AND has low arc → START OVERHANG
                        startFace = farthestFromStart.Face.FaceIndex;
                    }
                    else
                    {
                        // Normal case - start from closest to strip start
                        startFace = closestToStart.Face.FaceIndex;
                    }
                }
                else if (endpointFaces.Count == 1)
                {
                    startFace = endpointFaces.First().FaceIndex;
                }
                else
                {
                    // No clear endpoints - use closest to strip start
                    startFace = facesInStrip
                        .OrderBy(f => f.Center.DistanceTo(stripStart))
                        .First()
                        .FaceIndex;
                }
            }
            else if (centerCurve != null)
            {
                // Fallback: use center curve only
                Point3d curveStart = centerCurve.PointAtStart;

                var endpointFaces = adjacency
                    .Where(kvp => kvp.Value.Count <= 2)
                    .Select(kvp => facesInStrip.First(f => f.FaceIndex == kvp.Key))
                    .ToList();

                if (endpointFaces.Count >= 2)
                {
                    var endpointsSorted = endpointFaces
                        .Select(f => new { Face = f, Dist = f.Center.DistanceTo(curveStart), Arc = f.ArcParameter })
                        .OrderBy(x => x.Dist)
                        .ToList();

                    var closest = endpointsSorted.First();
                    var farthest = endpointsSorted.Last();

                    if (farthest.Arc < 0.3)
                        startFace = farthest.Face.FaceIndex;
                    else
                        startFace = closest.Face.FaceIndex;
                }
                else
                {
                    startFace = facesInStrip
                        .OrderBy(f => f.Center.DistanceTo(curveStart))
                        .First()
                        .FaceIndex;
                }
            }
            else
            {
                // No curve info - topology fallback
                if (adjacency.Count > 0)
                    startFace = adjacency.OrderBy(kvp => kvp.Value.Count).First().Key;
                else
                    startFace = facesInStrip.First().FaceIndex;
            }

            // ===== BFS TRAVERSAL =====
            var sorted = new List<FaceMetadata>();
            var visited = new HashSet<int>();
            int current = startFace;
            visited.Add(current);

            while (sorted.Count < facesInStrip.Count)
            {
                var currentMeta = facesInStrip.First(f => f.FaceIndex == current);
                sorted.Add(currentMeta);

                int next = -1;
                double minDist = double.MaxValue;

                if (adjacency.TryGetValue(current, out var neighbors))
                {
                    Point3d currentCenter = currentMeta.Center;

                    foreach (int neighborIdx in neighbors)
                    {
                        if (visited.Contains(neighborIdx)) continue;

                        var neighborMeta = facesInStrip.First(f => f.FaceIndex == neighborIdx);
                        double dist = currentCenter.DistanceTo(neighborMeta.Center);

                        if (dist < minDist)
                        {
                            minDist = dist;
                            next = neighborIdx;
                        }
                    }
                }

                if (next == -1)
                {
                    if (sorted.Count < facesInStrip.Count)
                    {
                        Point3d lastCenter = sorted.Last().Center;
                        var unvisited = facesInStrip.Where(f => !visited.Contains(f.FaceIndex)).ToList();

                        if (unvisited.Count > 0)
                        {
                            next = unvisited
                                .OrderBy(f => f.Center.DistanceTo(lastCenter))
                                .First()
                                .FaceIndex;
                        }
                    }

                    if (next == -1) break;
                }

                visited.Add(next);
                current = next;
            }

            return sorted;
        }
        private Point3d GetFaceCenter(Mesh mesh, MeshFace face)
        {
            Point3d center = Point3d.Origin;
            int count = 0;

            center += (Point3d)mesh.Vertices[face.A];
            center += (Point3d)mesh.Vertices[face.B];
            center += (Point3d)mesh.Vertices[face.C];
            count = 3;

            if (face.IsQuad)
            {
                center += (Point3d)mesh.Vertices[face.D];
                count = 4;
            }

            center /= count;
            return center;
        }

        private double GetParameterOnCurve(Point3d point, Curve curve)
        {
            if (curve == null) return 0.0;

            if (curve.ClosestPoint(point, out double t))
                return t;

            return 0.0;
        }

        private void CloseSmallNakedLoops(Mesh mesh, Surface surface, double maxLoopPerimeter)
        {
            var loops = mesh.GetNakedEdges();
            if (loops == null) return;

            foreach (var pl in loops)
            {
                if (pl == null || pl.Count < 4) continue;
                double per = 0.0; for (int i = 1; i < pl.Count; i++) per += pl[i - 1].DistanceTo(pl[i]);
                if (per > maxLoopPerimeter) continue;

                Point3d c = Point3d.Origin; int n = pl.Count - 1; for (int i = 0; i < n; i++) c += pl[i]; c /= n;
                surface.ClosestPoint(c, out double u, out double v); Point3d cp = surface.PointAt(u, v);
                int cIdx = mesh.Vertices.Add(cp);

                int facesBefore = mesh.Faces.Count;
                for (int i = 0; i < n; i++)
                {
                    int a = FindOrAddVertex(mesh, pl[i]);
                    int b = FindOrAddVertex(mesh, pl[(i + 1) % n]);
                    if (a == b || a == cIdx || b == cIdx) continue;
                    mesh.Faces.AddFace(cIdx, a, b);
                }
                int facesAfter = mesh.Faces.Count;

                // ===== NEW: Track loop closure faces =====
                for (int k = facesBefore; k < facesAfter; k++)
                {
                    Point3d faceCenter = GetFaceCenter(mesh, mesh.Faces[k]);

                    faceMetadataList.Add(new FaceMetadata
                    {
                        FaceIndex = k,
                        StripId = -1, // No strip assignment
                        StripType = -1,
                        CreationOrder = faceMetadataList.Count,
                        ArcParameter = 0.0,
                        Center = faceCenter
                    });
                }
            }
        }

        private int FindOrAddVertex(Mesh mesh, Point3d p)
        {
            // brute force for tiny loops; ok performance
            double tol = (Rhino.RhinoDoc.ActiveDoc?.ModelAbsoluteTolerance ?? 1e-3) * 2.0;
            for (int i = 0; i < mesh.Vertices.Count; i++)
                if (((Point3d)mesh.Vertices[i]).DistanceTo(p) <= tol)
                    return i;

            return mesh.Vertices.Add(p);
        }
        private void BuildEndpointCaps_Sub1(
    List<(int mainIdx, int sub1Idx)> capPairs,
    Dictionary<int, List<List<int>>> allBranchVertexIndices,
    Mesh mesh,
    Dictionary<int, List<Curve>> stripCenterCurves) 
        {
            if (capPairs == null || mesh == null) return;

            foreach (var (mainIdx, sub1Idx) in capPairs)
            {
                if (!allBranchVertexIndices.TryGetValue(mainIdx, out var mainCurves)) continue;
                if (!allBranchVertexIndices.TryGetValue(sub1Idx, out var subCurves)) continue;

                foreach (var sub in subCurves)
                {
                    // [Find bestMain/Start logic omitted for brevity, same as previous]
                    int bestMain = -1;
                    bool mainAtStart = true; bool subAtStart = true; double bestDist = double.MaxValue;
                    Point3d subS = (Point3d)mesh.Vertices[sub[0]]; Point3d subE = (Point3d)mesh.Vertices[sub[sub.Count - 1]];
                    for (int i = 0; i < mainCurves.Count; i++)
                    {
                        var main = mainCurves[i]; Point3d mainS = (Point3d)mesh.Vertices[main[0]]; Point3d mainE = (Point3d)mesh.Vertices[main[main.Count - 1]];
                        void Test(Point3d mp, bool mStart, Point3d sp, bool sStart) { double d = mp.DistanceTo(sp); if (d < bestDist) { bestDist = d; bestMain = i; mainAtStart = mStart; subAtStart = sStart; } }
                        Test(mainS, true, subS, true); Test(mainS, true, subE, false); Test(mainE, false, subS, true); Test(mainE, false, subE, false);
                    }

                    if (bestMain < 0) continue;
                    var mainPts = mainCurves[bestMain];

                    double tol = (Rhino.RhinoDoc.ActiveDoc?.ModelAbsoluteTolerance ?? 1e-3) * 10.0;
                    if (bestDist > tol) continue;

                    // Build Patch and Track Faces
                    int facesBefore = mesh.Faces.Count;
                    AddEndpointStitchPatch(mesh, mainPts, sub, mainAtStart, subAtStart, steps: 8);
                    int facesAfter = mesh.Faces.Count;

                    // ===== NEW: Track cap faces =====
                    for (int k = facesBefore; k < facesAfter; k++)
                    {
                        Point3d faceCenter = GetFaceCenter(mesh, mesh.Faces[k]);

                        faceMetadataList.Add(new FaceMetadata
                        {
                            FaceIndex = k,
                            StripId = mainIdx,
                            StripType = 2, // Cap type
                            CreationOrder = faceMetadataList.Count,
                            ArcParameter = 0.0, // Caps don't need ordering
                            Center = faceCenter
                        });
                    }
                }
            }
        }

        private void AddEndpointStitchPatch(
    Mesh mesh,
    List<int> main,
    List<int> sub,
    bool mainAtStart,
    bool subAtStart,
    int steps)
        {
            if (mesh == null || main == null || sub == null) return;
            if (main.Count < 2 || sub.Count < 2) return;
            if (steps < 1) steps = 1;

            List<int> A = BuildChainFromEndpoint(main, mainAtStart, steps + 1);
            List<int> B = BuildChainFromEndpoint(sub, subAtStart, steps + 1);

            if (A.Count < 2 || B.Count < 2) return;

            // Always force shared apex index for the patch to ensure topology closes
            B[0] = A[0];

            // Start pointers
            int i = 0, j = 0;

            // Stitch like CreateAdaptiveMesh_NoFans but constrained to this local patch
            while (i < A.Count - 1 || j < B.Count - 1)
            {
                int iNext = (i < A.Count - 1) ? i + 1 : i;
                int jNext = (j < B.Count - 1) ? j + 1 : j;

                int Ai = A[i];
                int Bj = B[j];
                int Ai1 = A[iNext];
                int Bj1 = B[jNext];

                // If both can advance, pick diagonal like your strip code
                if (i < A.Count - 1 && j < B.Count - 1)
                {
                    Point3d pAi = (Point3d)mesh.Vertices[Ai];
                    Point3d pAi1 = (Point3d)mesh.Vertices[Ai1];
                    Point3d pBj = (Point3d)mesh.Vertices[Bj];
                    Point3d pBj1 = (Point3d)mesh.Vertices[Bj1];

                    double diagA = pAi1.DistanceTo(pBj);  // advance A
                    double diagB = pAi.DistanceTo(pBj1);  // advance B

                    if (diagA <= diagB)
                    {
                        // Triangle (Ai, Bj, Ai1). Check for degeneracy first!
                        if (Ai != Bj && Ai != Ai1 && Bj != Ai1)
                        {
                            if (!TriangleExists(mesh, Ai, Bj, Ai1))
                                mesh.Faces.AddFace(Ai, Bj, Ai1);
                        }
                        i++;
                    }
                    else
                    {
                        // Triangle (Ai, Bj, Bj1). Check for degeneracy first!
                        if (Ai != Bj && Ai != Bj1 && Bj != Bj1)
                        {
                            if (!TriangleExists(mesh, Ai, Bj, Bj1))
                                mesh.Faces.AddFace(Ai, Bj, Bj1);
                        }
                        j++;
                    }
                }
                else if (i < A.Count - 1)
                {
                    // Only A can advance: triangle (Ai, Bj, Ai1)
                    if (Ai != Bj && Ai != Ai1 && Bj != Ai1)
                    {
                        if (!TriangleExists(mesh, Ai, Bj, Ai1))
                            mesh.Faces.AddFace(Ai, Bj, Ai1);
                    }
                    i++;
                }
                else if (j < B.Count - 1)
                {
                    // Only B can advance: triangle (Ai, Bj, Bj1)
                    if (Ai != Bj && Ai != Bj1 && Bj != Bj1)
                    {
                        if (!TriangleExists(mesh, Ai, Bj, Bj1))
                            mesh.Faces.AddFace(Ai, Bj, Bj1);
                    }
                    j++;
                }
                else break;
            }
        }

        private List<int> BuildChainFromEndpoint(List<int> poly, bool atStart, int count)
        {
            var chain = new List<int>(count);
            if (poly == null || poly.Count == 0) return chain;

            if (atStart)
            {
                for (int k = 0; k < poly.Count && chain.Count < count; k++)
                    chain.Add(poly[k]);
            }
            else
            {
                for (int k = poly.Count - 1; k >= 0 && chain.Count < count; k--)
                    chain.Add(poly[k]);
            }

            return chain;
        }

        private bool TriangleExists(Mesh mesh, int a, int b, int c)
        {
            // Checks if a triangle with the same vertex set already exists (any winding)
            // This is O(faceCount) but only used for a few corner caps so it's fine.
            var set = new HashSet<int> { a, b, c };

            for (int i = 0; i < mesh.Faces.Count; i++)
            {
                var f = mesh.Faces[i];
                if (!f.IsTriangle) continue;
                if (set.Contains(f.A) && set.Contains(f.B) && set.Contains(f.C))
                    return true;
            }
            return false;
        }
        private static List<(int mainIdx, int sub1Idx)> BuildCapPairs_Sub1(IList<GH_Path> paths)
        {
            var keyToIdx = new Dictionary<PathKey, int>();
            var bases = new SortedSet<int>();

            for (int i = 0; i < paths.Count; i++)
            {
                var k = KeyOf(paths[i]);
                keyToIdx[k] = i;
                bases.Add(k.Base);
            }

            var caps = new List<(int, int)>();

            foreach (int n in bases)
            {
                if (!keyToIdx.TryGetValue(new PathKey(n, -1), out int mainIdx))
                    continue;
                if (!keyToIdx.TryGetValue(new PathKey(n, 1), out int sub1Idx))
                    continue;

                caps.Add((mainIdx, sub1Idx));
            }

            return caps;
        }
        private (double a0, double a1) ArcRangeFromEndpointsOnPolyline(
    Mesh mesh,
    List<int> polylineVertexIndices,
    Dictionary<int, double> vertexToArcLength,
    Point3d end0,
    Point3d end1)
        {
            // find closest vertex (by 3D) on the sampled polyline to each endpoint
            int best0 = -1, best1 = -1;
            double d0 = double.MaxValue, d1 = double.MaxValue;

            for (int i = 0; i < polylineVertexIndices.Count; i++)
            {
                int vi = polylineVertexIndices[i];
                Point3d p = (Point3d)mesh.Vertices[vi];

                double dd0 = p.DistanceTo(end0);
                if (dd0 < d0) { d0 = dd0; best0 = vi; }

                double dd1 = p.DistanceTo(end1);
                if (dd1 < d1) { d1 = dd1; best1 = vi; }
            }

            if (best0 == -1 || best1 == -1) return (0, 0);

            double a0 = vertexToArcLength.TryGetValue(best0, out double aa0) ? aa0 : 0.0;
            double a1 = vertexToArcLength.TryGetValue(best1, out double aa1) ? aa1 : 0.0;

            if (a1 < a0) { double tmp = a0; a0 = a1; a1 = tmp; }
            return (a0, a1);
        }
        
        private readonly struct PathKey : IEquatable<PathKey>
        {
            public readonly int Base;
            public readonly int Sub; // -1 for {n}, otherwise {n;sub}
            public PathKey(int b, int s) { Base = b; Sub = s; }
            public bool Equals(PathKey other) => Base == other.Base && Sub == other.Sub;
            public override bool Equals(object obj) => obj is PathKey k && Equals(k);
            public override int GetHashCode() => (Base * 397) ^ Sub;
        }

        private static PathKey KeyOf(GH_Path p)
        {
            int b = p.Indices[0];
            int s = (p.Length > 1) ? p.Indices[1] : -1;
            return new PathKey(b, s);
        }

        private static List<(int idxA, int idxB, int stripType)> BuildMeshingPairs(IList<GH_Path> paths)
        {
            var keyToIdx = new Dictionary<PathKey, int>();
            var bases = new SortedSet<int>();

            for (int i = 0; i < paths.Count; i++)
            {
                var k = KeyOf(paths[i]);
                keyToIdx[k] = i;
                bases.Add(k.Base);
            }

            var pairs = new List<(int, int, int)>(); // Added stripType
            var processedAsNextStrip = new HashSet<int>();

            foreach (int n in bases)
            {
                if (processedAsNextStrip.Contains(n))
                    continue;

                if (!keyToIdx.TryGetValue(new PathKey(n, -1), out int mainIdx))
                    continue;

                bool hasSub0 = keyToIdx.TryGetValue(new PathKey(n, 0), out int sub0Idx);

                if (hasSub0)
                {
                    // Helper strip: type = 1
                    pairs.Add((mainIdx, sub0Idx, 1));

                    if (keyToIdx.TryGetValue(new PathKey(n + 1, -1), out int nextMainIdx))
                    {
                        bool nextHasSub0 = keyToIdx.TryGetValue(new PathKey(n + 1, 0), out int nextSub0Idx);

                        if (nextHasSub0)
                        {
                            pairs.Add((nextMainIdx, nextSub0Idx, 1));
                            processedAsNextStrip.Add(n + 1);
                        }
                        else if (keyToIdx.TryGetValue(new PathKey(n + 2, -1), out int n2MainIdx))
                        {
                            // Main strip: type = 0
                            pairs.Add((nextMainIdx, n2MainIdx, 0));
                            processedAsNextStrip.Add(n + 1);
                        }
                    }
                }
                else
                {
                    if (keyToIdx.TryGetValue(new PathKey(n + 1, -1), out int nextMainIdx))
                    {
                        // Main strip: type = 0
                        pairs.Add((mainIdx, nextMainIdx, 0));
                    }
                }
            }

            return pairs;
        }
        private void ExecuteMeshingLogic(
    List<(int idxA, int idxB, int stripType)> pairs,
    Dictionary<int, List<List<int>>> allBranchVertexIndices,
    Dictionary<int, List<Curve>> branchToCurves,
    Mesh masterMesh,
    Dictionary<int, int> vertexIndexToRow,
    Dictionary<int, double> vertexToArcLength,
    Dictionary<Point3d, int> pointToVertexIndex,
    Dictionary<int, Point2d> vertexToUV,
    Surface surface,
    double pullTol,
    out Dictionary<int, List<Curve>> stripCenterCurves)
        {

            var localStripCurves = new Dictionary<int, List<Curve>>();
            int globalCreationOrder = 0;

            foreach (var (idxA, idxB, stripType) in pairs)
            {
                if (!allBranchVertexIndices.ContainsKey(idxA) || !allBranchVertexIndices.ContainsKey(idxB))
                    continue;

                List<List<int>> branchA = allBranchVertexIndices[idxA];
                List<List<int>> branchB = allBranchVertexIndices[idxB];
                List<Curve> curvesA = branchToCurves[idxA];
                List<Curve> curvesB = branchToCurves[idxB];

                int numA = branchA.Count;
                int numB = branchB.Count;

                void ProcessStrip(List<int> indsA, List<int> indsB, Curve cA, Curve cB)
                {
                    if (indsA.Count < 2 || indsB.Count < 2) return;

                    // Capture face count BEFORE adding new faces
                    int facesBefore = masterMesh.Faces.Count;

                    List<Point3d> ptsA = indsA.Select(x => (Point3d)masterMesh.Vertices[x]).ToList();
                    List<Point3d> ptsB = indsB.Select(x => (Point3d)masterMesh.Vertices[x]).ToList();

                    int startA, startB;
                    FindClosestStartIndices(ptsA, ptsB, out startA, out startB);

                    bool isClosed = cA.IsClosed && cB.IsClosed;

                    int minA, minB, maxA, maxB;
                    if (isClosed)
                    {
                        CreateAdaptiveMesh_Closed(masterMesh, indsA, indsB, startA, startB);
                        minA = 0; minB = 0; maxA = ptsA.Count - 1; maxB = ptsB.Count - 1;
                    }
                    else
                    {
                        CreateAdaptiveMesh_NoFans(masterMesh, indsA, indsB, startA, startB,
                            out minA, out minB, out maxA, out maxB);
                    }

                    List<Point3d> startOverhangPts, endOverhangPts;
                    ProcessOverhangs(masterMesh, idxA, indsA, indsB, ptsA, ptsB, minA, minB, maxA, maxB,
                        vertexIndexToRow, vertexToArcLength, pointToVertexIndex, vertexToUV,
                        surface, pullTol, out startOverhangPts, out endOverhangPts);

                    int facesAfter = masterMesh.Faces.Count;

                    // --- Construct Debug Curve FIRST (we need it for metadata) ---
                    List<Point3d> debugPoints = new List<Point3d>();
                    if (startOverhangPts != null && startOverhangPts.Count > 0)
                        debugPoints.AddRange(startOverhangPts);

                    List<Point3d> bodyPts = GetMainBodyCenters(ptsA, ptsB, startA, startB, minA, minB, maxA, maxB, isClosed);
                    debugPoints.AddRange(bodyPts);

                    if (endOverhangPts != null && endOverhangPts.Count > 0)
                    {
                        endOverhangPts.Reverse();
                        debugPoints.AddRange(endOverhangPts);
                    }

                    Curve centerCurve = null;
                    if (debugPoints.Count >= 2)
                    {
                        List<Point3d> cleanPts = new List<Point3d> { debugPoints[0] };
                        for (int k = 1; k < debugPoints.Count; k++)
                        {
                            if (cleanPts.Last().DistanceTo(debugPoints[k]) > 1e-6)
                                cleanPts.Add(debugPoints[k]);
                        }

                        if (cleanPts.Count >= 2)
                        {
                            if (isClosed && cleanPts.Count > 2 && cleanPts[0].DistanceTo(cleanPts.Last()) < 1e-6)
                                cleanPts.RemoveAt(cleanPts.Count - 1);

                            centerCurve = isClosed
                                ? Curve.CreateInterpolatedCurve(cleanPts, 3, CurveKnotStyle.ChordPeriodic)
                                : Curve.CreateInterpolatedCurve(cleanPts, 3);

                            if (centerCurve != null)
                            {
                                if (!localStripCurves.ContainsKey(idxA))
                                    localStripCurves[idxA] = new List<Curve>();
                                localStripCurves[idxA].Add(centerCurve);
                            }
                        }
                    }

                    // ===== NEW METADATA TRACKING =====
                    for (int k = facesBefore; k < facesAfter; k++)
                    {
                        Point3d faceCenter = GetFaceCenter(masterMesh, masterMesh.Faces[k]);

                        // Compute arc parameter (0-1 along strip)
                        double arcParam = 0.0;
                        if (centerCurve != null && centerCurve.ClosestPoint(faceCenter, out double t))
                        {
                            arcParam = centerCurve.Domain.NormalizedParameterAt(t);
                        }

                        faceMetadataList.Add(new FaceMetadata
                        {
                            FaceIndex = k,
                            StripId = idxA,
                            StripType = stripType,
                            CreationOrder = globalCreationOrder++,
                            ArcParameter = arcParam,
                            Center = faceCenter,
                            // ===== NEW: Store boundary curve start points =====
                            BoundaryStartA = cA != null ? cA.PointAtStart : Point3d.Origin,
                            BoundaryStartB = cB != null ? cB.PointAtStart : Point3d.Origin
                        });
                    }
                }

                if (numA == numB)
                {
                    for (int j = 0; j < numA; j++) ProcessStrip(branchA[j], branchB[j], curvesA[j], curvesB[j]);
                }
                else if (numA == 1 && numB > 1)
                {
                    Curve cA = curvesA[0]; List<int> iA = branchA[0];
                    for (int j = 0; j < numB; j++)
                    {
                        Curve cB = curvesB[j]; List<int> iB = branchB[j];
                        var (a0, a1) = ArcRangeFromEndpointsOnPolyline(masterMesh, iA, vertexToArcLength, cB.PointAtStart, cB.PointAtEnd);
                        List<int> subA = iA.Where(v => { if (!vertexToArcLength.TryGetValue(v, out double a)) return false; return a >= a0 && a <= a1; }).ToList();
                        ProcessStrip(subA, iB, cA, cB);
                    }
                }
                else if (numB == 1 && numA > 1)
                {
                    Curve cB = curvesB[0]; List<int> iB = branchB[0];
                    for (int j = 0; j < numA; j++)
                    {
                        Curve cA = curvesA[j]; List<int> iA = branchA[j];
                        var (a0, a1) = ArcRangeFromEndpointsOnPolyline(masterMesh, iB, vertexToArcLength, cA.PointAtStart, cA.PointAtEnd);
                        List<int> subB = iB.Where(v => { if (!vertexToArcLength.TryGetValue(v, out double a)) return false; return a >= a0 && a <= a1; }).ToList();
                        ProcessStrip(iA, subB, cA, cB);
                    }
                }
            }
            stripCenterCurves = localStripCurves;
        }

        // --- NEW HELPER: Replicates the meshing logic to get accurate center points ---
        private List<Point3d> GetMainBodyCenters(List<Point3d> pA, List<Point3d> pB,
            int startA, int startB, int minA, int minB, int maxA, int maxB, bool isClosed)
        {
            List<Point3d> centers = new List<Point3d>();

            // 1. Backward loop (Reverse results later)
            List<Point3d> backward = new List<Point3d>();
            int i = startA; int j = startB;

            // Fix for Intersection Case: Check if start points are effectively the same
            bool startsIntersect = pA[startA].DistanceTo(pB[startB]) < 0.001;

            if (!isClosed)
            {
                while (i > minA && j > minB)
                {
                    double d1 = pA[i - 1].DistanceTo(pB[j]);
                    double d2 = pA[i].DistanceTo(pB[j - 1]);
                    Point3d c;
                    if (d1 < d2)
                    {
                        c = (pA[i] + pB[j] + pA[i - 1]) / 3.0; // Triangle center
                        i--;
                    }
                    else
                    {
                        c = (pA[i] + pB[j] + pB[j - 1]) / 3.0; // Triangle center
                        j--;
                    }
                    backward.Add(c);
                }
            }
            backward.Reverse();
            centers.AddRange(backward);

            // 2. Intersection Anchor
            // If the curves intersect at the start, explicitly anchor the curve there.
            if (startsIntersect)
            {
                centers.Add((pA[startA] + pB[startB]) * 0.5);
            }

            // 3. Forward loop
            i = startA; j = startB;
            int limitA = isClosed ? pA.Count : maxA;
            int limitB = isClosed ? pB.Count : maxB;

            // Simple loop guard for closed curves to prevent infinite loops (naive approach)
            int steps = 0;
            int maxSteps = (pA.Count + pB.Count) * 2;

            while ((isClosed && steps < maxSteps) || (!isClosed && i < limitA && j < limitB))
            {
                // Break for closed curves if we wrapped around
                if (isClosed && steps > 0 && i == startA && j == startB) break;

                int nextI = (i + 1) % pA.Count;
                int nextJ = (j + 1) % pB.Count;

                // Stop at ends for open curves
                if (!isClosed && (i == limitA || j == limitB)) break;
                if (!isClosed && (nextI > limitA || nextJ > limitB)) break; // Index check

                double d1 = pA[nextI].DistanceTo(pB[j]);
                double d2 = pA[i].DistanceTo(pB[nextJ]);

                Point3d c;
                if (d1 < d2)
                {
                    c = (pA[i] + pB[j] + pA[nextI]) / 3.0;
                    i = nextI;
                }
                else
                {
                    c = (pA[i] + pB[j] + pB[nextJ]) / 3.0;
                    j = nextJ;
                }
                centers.Add(c);
                steps++;
            }

            return centers;
        }

        // --- UPDATED OVERHANG LOGIC: Captures Center Points ---

        private void ProcessOverhangs(Mesh masterMesh, int rowIndex,
            List<int> indicesA, List<int> indicesB,
            List<Point3d> pointsA, List<Point3d> pointsB,
            int minA, int minB, int maxA, int maxB,
            Dictionary<int, int> vertexIndexToRow, Dictionary<int, double> vertexToArcLength,
            Dictionary<Point3d, int> pointToVertexIndex, Dictionary<int, Point2d> vertexToUV,
            Surface surface, double pullTol,
            out List<Point3d> startOverhangPts, out List<Point3d> endOverhangPts)
        {
            startOverhangPts = new List<Point3d>();
            endOverhangPts = new List<Point3d>();

            int lastA = pointsA.Count - 1;
            int lastB = pointsB.Count - 1;

            // End Overhangs
            if (maxA < lastA)
            {
                int k = lastA - maxA;
                if (k > 0) FillOverhang_End(masterMesh, rowIndex, indicesA, maxA, lastA, indicesB[lastB], pointsA, pointsB[lastB], k, vertexIndexToRow, vertexToArcLength, pointToVertexIndex, vertexToUV, surface, true, pullTol, out endOverhangPts);
            }
            else if (maxB < lastB)
            {
                int k = lastB - maxB;
                if (k > 0) FillOverhang_End(masterMesh, rowIndex, indicesB, maxB, lastB, indicesA[lastA], pointsB, pointsA[lastA], k, vertexIndexToRow, vertexToArcLength, pointToVertexIndex, vertexToUV, surface, false, pullTol, out endOverhangPts);
            }

            // Start Overhangs
            if (minA > 0)
            {
                int k = minA;
                if (k > 0) FillOverhang_Start(masterMesh, rowIndex, indicesA, minA, indicesB[0], pointsA, pointsB[0], k, vertexIndexToRow, vertexToArcLength, pointToVertexIndex, vertexToUV, surface, true, pullTol, out startOverhangPts);
            }
            else if (minB > 0)
            {
                int k = minB;
                if (k > 0) FillOverhang_Start(masterMesh, rowIndex, indicesB, minB, indicesA[0], pointsB, pointsA[0], k, vertexIndexToRow, vertexToArcLength, pointToVertexIndex, vertexToUV, surface, false, pullTol, out startOverhangPts);
            }
        }

        private void FillOverhang_Start(Mesh mesh, int row, List<int> baseIndices, int minIdx, int apexV,
          List<Point3d> basePoints, Point3d apexP, int k,
          Dictionary<int, int> rowMap, Dictionary<int, double> arcMap, Dictionary<Point3d, int> ptMap,
          Dictionary<int, Point2d> uvMap, Surface surface, bool isA, double pullTol,
          out List<Point3d> createdCenters)
        {
            createdCenters = new List<Point3d>();
            int baseEndV = baseIndices[0];
            Point2d uvApex = uvMap[apexV];
            Point2d uvEnd = uvMap[baseEndV];

            // Build the "ladder" vertices
            List<double> segs = new List<double>();
            double total = 0.0;
            for (int s = 0; s < k; s++)
            {
                double l = basePoints[minIdx - s].DistanceTo(basePoints[minIdx - s - 1]);
                segs.Add(l);
                total += l;
            }
            List<int> ladderIndices = new List<int> { apexV };
            List<Point3d> ladderPoints = new List<Point3d> { apexP };

            double cum = 0.0;
            for (int s = 0; s < k - 1; s++)
            {
                cum += segs[s] / total;
                double t = cum;
                double u = uvApex.X + (uvEnd.X - uvApex.X) * t;
                double v = uvApex.Y + (uvEnd.Y - uvApex.Y) * t;
                Point3d p = surface.PointAt(u, v);
                int newIndex = AddVertex(mesh, p, ptMap);
                if (!rowMap.ContainsKey(newIndex)) rowMap.Add(newIndex, row);
                if (!uvMap.ContainsKey(newIndex)) uvMap.Add(newIndex, new Point2d(u, v));
                double arc1 = arcMap[baseIndices[minIdx - s]];
                double arc2 = arcMap[baseIndices[minIdx - s - 1]];
                if (!arcMap.ContainsKey(newIndex)) arcMap.Add(newIndex, (arc1 + arc2) * 0.5);

                ladderIndices.Add(newIndex);
                ladderPoints.Add(p);
            }
            ladderIndices.Add(baseEndV);
            ladderPoints.Add(mesh.Vertices[baseEndV]);

            // Create Faces and Capture Centers
            for (int i = 0; i < k; i++)
            {
                int L_curr = ladderIndices[i];
                int L_next = ladderIndices[i + 1];
                int B_curr = baseIndices[minIdx - i];
                int B_next = baseIndices[minIdx - i - 1];

                Point3d pLc = ladderPoints[i];
                Point3d pLn = ladderPoints[i + 1];
                Point3d pBc = basePoints[minIdx - i];
                Point3d pBn = basePoints[minIdx - i - 1];

                bool isTip = (i == k - 1);
                Point3d center;

                if (isTip)
                {
                    if (isA) mesh.Faces.AddFace(L_curr, B_curr, B_next);
                    else mesh.Faces.AddFace(L_curr, B_next, B_curr);
                    // Centroid of Tip Triangle
                    center = (pLc + pBc + pBn) / 3.0;
                }
                else
                {
                    if (isA) mesh.Faces.AddFace(L_curr, B_curr, B_next, L_next);
                    else mesh.Faces.AddFace(L_curr, L_next, B_next, B_curr);
                    // Centroid of Quad
                    center = (pLc + pBc + pBn + pLn) / 4.0;
                }
                createdCenters.Add(center);
            }
        }

        private void FillOverhang_End(Mesh mesh, int row, List<int> baseIndices, int maxIdx, int lastIdx, int apexV,
          List<Point3d> basePoints, Point3d apexP, int k,
          Dictionary<int, int> rowMap, Dictionary<int, double> arcMap, Dictionary<Point3d, int> ptMap,
          Dictionary<int, Point2d> uvMap, Surface surface, bool isA, double pullTol,
          out List<Point3d> createdCenters)
        {
            createdCenters = new List<Point3d>();
            int baseEndV = baseIndices[lastIdx];
            Point2d uvApex = uvMap[apexV];
            Point2d uvEnd = uvMap[baseEndV];

            List<double> segs = new List<double>();
            double total = 0.0;
            for (int s = 0; s < k; s++)
            {
                double l = basePoints[maxIdx + s].DistanceTo(basePoints[maxIdx + s + 1]);
                segs.Add(l);
                total += l;
            }
            List<int> ladderIndices = new List<int> { apexV };
            List<Point3d> ladderPoints = new List<Point3d> { apexP };

            double cum = 0.0;
            for (int s = 0; s < k - 1; s++)
            {
                cum += segs[s] / total;
                double t = cum;
                double u = uvApex.X + (uvEnd.X - uvApex.X) * t;
                double v = uvApex.Y + (uvEnd.Y - uvApex.Y) * t;
                Point3d p = surface.PointAt(u, v);
                int newIndex = AddVertex(mesh, p, ptMap);
                if (!rowMap.ContainsKey(newIndex)) rowMap.Add(newIndex, row);
                if (!uvMap.ContainsKey(newIndex)) uvMap.Add(newIndex, new Point2d(u, v));
                double arc1 = arcMap[baseIndices[maxIdx + s]];
                double arc2 = arcMap[baseIndices[maxIdx + s + 1]];
                if (!arcMap.ContainsKey(newIndex)) arcMap.Add(newIndex, (arc1 + arc2) * 0.5);

                ladderIndices.Add(newIndex);
                ladderPoints.Add(p);
            }
            ladderIndices.Add(baseEndV);
            ladderPoints.Add(mesh.Vertices[baseEndV]);

            for (int i = 0; i < k; i++)
            {
                int L_curr = ladderIndices[i];
                int L_next = ladderIndices[i + 1];
                int B_curr = baseIndices[maxIdx + i];
                int B_next = baseIndices[maxIdx + i + 1];

                Point3d pLc = ladderPoints[i];
                Point3d pLn = ladderPoints[i + 1];
                Point3d pBc = basePoints[maxIdx + i];
                Point3d pBn = basePoints[maxIdx + i + 1];

                bool isTip = (i == k - 1);
                Point3d center;

                if (isTip)
                {
                    if (isA) mesh.Faces.AddFace(L_curr, B_curr, B_next);
                    else mesh.Faces.AddFace(L_curr, B_next, B_curr);
                    // Centroid of Tip Triangle
                    center = (pLc + pBc + pBn) / 3.0;
                }
                else
                {
                    if (isA) mesh.Faces.AddFace(L_curr, B_curr, B_next, L_next);
                    else mesh.Faces.AddFace(L_curr, L_next, B_next, B_curr);
                    // Centroid of Quad
                    center = (pLc + pBc + pBn + pLn) / 4.0;
                }
                createdCenters.Add(center);
            }
        }

        private List<Point3d> GetProjectedPointsFromCurve(Curve curve, double length, Surface surface, double pullTol)
        {
            List<Point3d> points = new List<Point3d>();
            if (curve == null) return points;
            if (pullTol > 0)
            {
                var brep = Brep.CreateFromSurface(surface);
                var pulled = curve.PullToBrepFace(brep.Faces[0], pullTol);
                if (pulled != null && pulled.Length > 0) curve = pulled[0];
            }
            double curveLength = curve.GetLength();
            int optimalDivisions = Math.Max(2, (int)Math.Round(curveLength / length));
            double[] t_params = curve.DivideByCount(optimalDivisions, true);
            if (t_params == null || t_params.Length < 2) t_params = curve.DivideByLength(length, true);
            if (t_params == null) return points;
            foreach (double t in t_params)
            {
                Point3d pt = curve.PointAt(t);
                surface.ClosestPoint(pt, out double u, out double v);
                points.Add(surface.PointAt(u, v));
            }
            if (curve.IsClosed && points.Count > 1 && points[0].DistanceTo(points[points.Count - 1]) < 1e-9)
                points.RemoveAt(points.Count - 1);
            return points;
        }

        private void FixBoundaryVertices(Mesh mesh, Surface surface, Dictionary<int, Point2d> vertexToUV)
        {
            var boundaryEdges = mesh.GetNakedEdges();
            foreach (Polyline boundary in boundaryEdges)
            {
                for (int i = 0; i < boundary.Count - 1; i++)
                {
                    for (int q = 0; q < mesh.Vertices.Count; q++)
                    {
                        Point3d meshPt = (Point3d)mesh.Vertices[q];
                        if (meshPt.DistanceTo(boundary[i]) < 0.001 || meshPt.DistanceTo(boundary[i + 1]) < 0.001)
                        {
                            if (vertexToUV.ContainsKey(q))
                            {
                                Point2d uv = vertexToUV[q];
                                mesh.Vertices[q] = (Point3f)surface.PointAt(uv.X, uv.Y);
                            }
                            else
                            {
                                surface.ClosestPoint(meshPt, out double u, out double v);
                                mesh.Vertices[q] = (Point3f)surface.PointAt(u, v);
                                vertexToUV[q] = new Point2d(u, v);
                            }
                        }
                    }
                }
            }
        }

        private int AddVertex(Mesh mesh, Point3d pt, Dictionary<Point3d, int> lookup)
        {
            if (lookup.ContainsKey(pt)) return lookup[pt];
            int vIndex = mesh.Vertices.Add(pt);
            lookup.Add(pt, vIndex);
            return vIndex;
        }

        private void FindClosestStartIndices(List<Point3d> pointsA, List<Point3d> pointsB, out int bestIdxA, out int bestIdxB)
        {
            bestIdxA = 0; bestIdxB = 0;
            double minDist = double.MaxValue;
            for (int ii = 0; ii < pointsA.Count; ii++)
                for (int jj = 0; jj < pointsB.Count; jj++)
                {
                    double dist = pointsA[ii].DistanceTo(pointsB[jj]);
                    if (dist < minDist) { minDist = dist; bestIdxA = ii; bestIdxB = jj; }
                }
        }

        private void CreateAdaptiveMesh_NoFans(Mesh mesh, List<int> indicesA, List<int> indicesB, int startA, int startB, out int minA, out int minB, out int maxA, out int maxB)
        {
            int i_fwd = startA; int j_fwd = startB;
            while (i_fwd < indicesA.Count - 1 && j_fwd < indicesB.Count - 1)
            {
                Point3d pA0 = mesh.Vertices[indicesA[i_fwd]]; Point3d pA1 = mesh.Vertices[indicesA[i_fwd + 1]];
                Point3d pB0 = mesh.Vertices[indicesB[j_fwd]]; Point3d pB1 = mesh.Vertices[indicesB[j_fwd + 1]];
                double diag1 = pA1.DistanceTo(pB0); double diag2 = pA0.DistanceTo(pB1);
                if (diag1 < diag2) { mesh.Faces.AddFace(indicesA[i_fwd], indicesB[j_fwd], indicesA[i_fwd + 1]); i_fwd++; }
                else { mesh.Faces.AddFace(indicesA[i_fwd], indicesB[j_fwd], indicesB[j_fwd + 1]); j_fwd++; }
            }
            int i_bwd = startA; int j_bwd = startB;
            while (i_bwd > 0 && j_bwd > 0)
            {
                Point3d pA0 = mesh.Vertices[indicesA[i_bwd]]; Point3d pA1 = mesh.Vertices[indicesA[i_bwd - 1]];
                Point3d pB0 = mesh.Vertices[indicesB[j_bwd]]; Point3d pB1 = mesh.Vertices[indicesB[j_bwd - 1]];
                double diag1 = pA1.DistanceTo(pB0); double diag2 = pA0.DistanceTo(pB1);
                if (diag1 < diag2) { mesh.Faces.AddFace(indicesA[i_bwd], indicesB[j_bwd], indicesA[i_bwd - 1]); i_bwd--; }
                else { mesh.Faces.AddFace(indicesA[i_bwd], indicesB[j_bwd], indicesB[j_bwd - 1]); j_bwd--; }
            }
            minA = i_bwd; minB = j_bwd; maxA = i_fwd; maxB = j_fwd;
        }

        private void CreateAdaptiveMesh_Closed(Mesh mesh, List<int> indicesA, List<int> indicesB, int startA, int startB)
        {
            int nA = indicesA.Count;
            int nB = indicesB.Count;
            int i = startA;
            int j = startB;
            int stepsA = 0;
            int stepsB = 0;
            while (stepsA < nA || stepsB < nB)
            {
                int iNext = (i + 1) % nA;
                int jNext = (j + 1) % nB;
                bool moveA = false;
                if (stepsA == nA) { moveA = false; }
                else if (stepsB == nB) { moveA = true; }
                else
                {
                    Point3d pA_curr = mesh.Vertices[indicesA[i]];
                    Point3d pA_next = mesh.Vertices[indicesA[iNext]];
                    Point3d pB_curr = mesh.Vertices[indicesB[j]];
                    Point3d pB_next = mesh.Vertices[indicesB[jNext]];
                    double d1 = pA_next.DistanceTo(pB_curr);
                    double d2 = pA_curr.DistanceTo(pB_next);
                    moveA = (d1 < d2);
                }
                if (moveA)
                {
                    mesh.Faces.AddFace(indicesA[i], indicesB[j], indicesA[iNext]);
                    i = iNext;
                    stepsA++;
                }
                else
                {
                    mesh.Faces.AddFace(indicesA[i], indicesB[j], indicesB[jNext]);
                    j = jNext;
                    stepsB++;
                }
            }
        }

        private void EnsureConsistentNormals(Mesh masterMesh, Dictionary<int, int> vertexIndexToRow, Dictionary<int, List<Curve>> branchToCurves)
        {
            if (vertexIndexToRow.Count > 0 && masterMesh.Faces.Count > 0)
            {
                MeshFace f = masterMesh.Faces[0];
                List<int> vs = new List<int> { f.A, f.B, f.C }; if (f.IsQuad) vs.Add(f.D);
                var rows = vs.Select(v => vertexIndexToRow.ContainsKey(v) ? vertexIndexToRow[v] : -1).Distinct().Where(r => r != -1).ToList();
                if (rows.Count == 2)
                {
                    int row1 = rows[0], row2 = rows[1];
                    if (row1 > row2) { int temp = row1; row1 = row2; row2 = temp; }
                    Point3d avg1 = Point3d.Origin; int count1 = 0; Point3d avg2 = Point3d.Origin; int count2 = 0;
                    foreach (int v in vs)
                    {
                        if (!vertexIndexToRow.ContainsKey(v)) continue;
                        if (vertexIndexToRow[v] == row1) { avg1 += (Point3d)masterMesh.Vertices[v]; count1++; }
                        else { avg2 += (Point3d)masterMesh.Vertices[v]; count2++; }
                    }
                    if (count1 > 0) avg1 /= count1; if (count2 > 0) avg2 /= count2;
                    Vector3d dir = avg2 - avg1;
                    Vector3d fallbackNormal = masterMesh.FaceNormals[0];
                    if (Vector3d.Multiply(fallbackNormal, dir) < 0) masterMesh.Flip(true, true, true);
                }
            }
        }

        private void OptimizedConvertTrianglesToQuads(Mesh mesh, double angleTolRad, double aspectTol, double minCornerRad)
        {
            double maxCornerRad = Math.PI - minCornerRad;
            bool changed = true;
            while (changed)
            {
                changed = false;
                var edgeToFaces = new Dictionary<EdgeKey, List<int>>();
                int faceCount = mesh.Faces.Count;
                for (int i = 0; i < faceCount; i++)
                {
                    if (!mesh.Faces[i].IsTriangle) continue;
                    MeshFace f = mesh.Faces[i];
                    AddEdge(f.A, f.B, i, edgeToFaces); AddEdge(f.B, f.C, i, edgeToFaces); AddEdge(f.C, f.A, i, edgeToFaces);
                }

                var candidates = edgeToFaces.Where(kv => kv.Value.Count == 2).ToList();
                var facesToDelete = new HashSet<int>();
                var newFaces = new List<MeshFace>();
                var newFaceMetadata = new List<FaceMetadata>(); // ===== NEW =====

                foreach (var kv in candidates)
                {
                    int f1 = kv.Value[0];
                    int f2 = kv.Value[1];
                    if (facesToDelete.Contains(f1) || facesToDelete.Contains(f2)) continue;

                    if (CanMerge(mesh, f1, f2, angleTolRad, aspectTol, minCornerRad, maxCornerRad))
                    {
                        MeshFace quad = CreateMergedFace(mesh, f1, f2);
                        newFaces.Add(quad);

                        // ===== NEW: Inherit metadata from first parent =====
                        var parent1Meta = faceMetadataList.FirstOrDefault(m => m.FaceIndex == f1);
                        if (parent1Meta != null)
                        {
                            newFaceMetadata.Add(new FaceMetadata
                            {
                                FaceIndex = -1, // Will be updated after adding to mesh
                                StripId = parent1Meta.StripId,
                                StripType = parent1Meta.StripType,
                                CreationOrder = parent1Meta.CreationOrder,
                                ArcParameter = parent1Meta.ArcParameter,
                                Center = parent1Meta.Center
                            });
                        }

                        facesToDelete.Add(f1);
                        facesToDelete.Add(f2);
                        changed = true;
                    }
                }

                if (changed)
                {
                    int firstNewFaceIdx = mesh.Faces.Count;
                    mesh.Faces.AddFaces(newFaces);

                    // ===== NEW: Update face indices and add metadata =====
                    for (int i = 0; i < newFaceMetadata.Count; i++)
                    {
                        newFaceMetadata[i].FaceIndex = firstNewFaceIdx + i;
                        faceMetadataList.Add(newFaceMetadata[i]);
                    }

                    mesh.Faces.DeleteFaces(facesToDelete);

                    // ===== NEW: Remove deleted face metadata =====
                    faceMetadataList.RemoveAll(m => facesToDelete.Contains(m.FaceIndex));

                    // ===== NEW: Reindex remaining faces =====
                    var indexMapping = new Dictionary<int, int>();
                    int newIdx = 0;
                    for (int oldIdx = 0; oldIdx < mesh.Faces.Count + facesToDelete.Count; oldIdx++)
                    {
                        if (!facesToDelete.Contains(oldIdx))
                        {
                            indexMapping[oldIdx] = newIdx;
                            newIdx++;
                        }
                    }

                    foreach (var meta in faceMetadataList)
                    {
                        if (indexMapping.TryGetValue(meta.FaceIndex, out int mappedIdx))
                            meta.FaceIndex = mappedIdx;
                    }
                }
            }
        }
        private MeshFace CreateMergedFace(Mesh mesh, int f1, int f2)
        {
            MeshFace face1 = mesh.Faces[f1]; int a = face1.A, b = face1.B, c = face1.C; List<int> tri1 = new List<int> { a, b, c };
            MeshFace face2 = mesh.Faces[f2]; int d = face2.A, e = face2.B, ff = face2.C; List<int> tri2 = new List<int> { d, e, ff };
            var shared = tri1.Intersect(tri2).ToList();
            int sh1 = shared[0]; int sh2 = shared[1];
            var uniq1 = tri1.Except(shared).First(); var uniq2 = tri2.Except(shared).First();
            return new MeshFace(sh1, uniq1, sh2, uniq2);
        }

        private void AddEdge(int a, int b, int fi, Dictionary<EdgeKey, List<int>> dict)
        {
            EdgeKey key = new EdgeKey(a, b);
            if (!dict.ContainsKey(key)) dict[key] = new List<int>();
            dict[key].Add(fi);
        }

        class EdgeKey : IEquatable<EdgeKey>
        {
            public int V1, V2;
            public EdgeKey(int v1, int v2) { V1 = Math.Min(v1, v2); V2 = Math.Max(v1, v2); }
            public bool Equals(EdgeKey other) => other != null && V1 == other.V1 && V2 == other.V2;
            public override bool Equals(object obj) => Equals(obj as EdgeKey);
            public override int GetHashCode() => V1 ^ V2;
        }

        private bool CanMerge(Mesh mesh, int f1, int f2, double angleTolRad, double aspectTol, double minCornerRad, double maxCornerRad)
        {
            MeshFace face1 = mesh.Faces[f1]; int a = face1.A, b = face1.B, c = face1.C; List<int> tri1 = new List<int> { a, b, c };
            MeshFace face2 = mesh.Faces[f2]; int d = face2.A, e = face2.B, ff = face2.C; List<int> tri2 = new List<int> { d, e, ff };
            var shared = tri1.Intersect(tri2).ToList();
            if (shared.Count != 2) return false;
            var uniq1 = tri1.Except(shared).ToList(); var uniq2 = tri2.Except(shared).ToList();
            if (uniq1.Count != 1 || uniq2.Count != 1) return false;
            int sh1 = shared[0]; int sh2 = shared[1]; int u1 = uniq1[0]; int u2 = uniq2[0];
            Point3d pa = mesh.Vertices[a]; Point3d pb = mesh.Vertices[b]; Point3d pc = mesh.Vertices[c];
            Vector3d norm1 = Vector3d.CrossProduct(pb - pa, pc - pa); if (norm1.Length == 0) return false; norm1.Unitize();
            Point3d pd = mesh.Vertices[d]; Point3d pe = mesh.Vertices[e]; Point3d pff = mesh.Vertices[ff];
            Vector3d norm2 = Vector3d.CrossProduct(pe - pd, pff - pd); if (norm2.Length == 0) return false; norm2.Unitize();
            if (Vector3d.VectorAngle(norm1, norm2) > angleTolRad) return false;
            double diagShared = mesh.Vertices[sh1].DistanceTo(mesh.Vertices[sh2]);
            double diagNew = mesh.Vertices[u1].DistanceTo(mesh.Vertices[u2]);
            double maxD = Math.Max(diagShared, diagNew); double minD = Math.Min(diagShared, diagNew);
            if (minD <= 1e-9 || minD / maxD < aspectTol) return false;
            double angSh1 = Vector3d.VectorAngle(mesh.Vertices[u1] - mesh.Vertices[sh1], mesh.Vertices[u2] - mesh.Vertices[sh1]);
            double angSh2 = Vector3d.VectorAngle(mesh.Vertices[u1] - mesh.Vertices[sh2], mesh.Vertices[u2] - mesh.Vertices[sh2]);
            double angU1 = Vector3d.VectorAngle(mesh.Vertices[sh1] - mesh.Vertices[u1], mesh.Vertices[sh2] - mesh.Vertices[u1]);
            double angU2 = Vector3d.VectorAngle(mesh.Vertices[sh1] - mesh.Vertices[u2], mesh.Vertices[sh2] - mesh.Vertices[u2]);
            if (angSh1 < minCornerRad || angSh1 > maxCornerRad) return false;
            if (angSh2 < minCornerRad || angSh2 > maxCornerRad) return false;
            if (angU1 < minCornerRad || angU1 > maxCornerRad) return false;
            if (angU2 < minCornerRad || angU2 > maxCornerRad) return false;
            return true;
        }
        private void BuildDiscreteBlocksWithCorrectOrder(
     Mesh mesh,
     DataTree<Mesh> blocksTree,
     double thickness,
     bool flip,
     GH_Structure<GH_String> patternsTree,
     Dictionary<int, List<Curve>> stripCenterCurves) // ← Parameter (already declared)
        {
            if (flip) thickness = -thickness;

            // --- PATTERN SETUP ---
            List<string> patternList = new List<string>();
            HashSet<char> uniqueChars = new HashSet<char>();

            if (patternsTree != null && patternsTree.Branches.Count > 0)
            {
                foreach (var branch in patternsTree.Branches)
                {
                    string pVal = (branch.Count > 0 && branch[0] != null) ? branch[0].Value : "A";
                    patternList.Add(pVal);
                    string cleanP = GetBasePattern(pVal);
                    foreach (char c in cleanP) uniqueChars.Add(c);
                }
            }
            else
            {
                patternList.Add("A");
                uniqueChars.Add('A');
            }

            // --- COLOR MAP ---
            Dictionary<char, Color> globalColorMap = new Dictionary<char, Color>();
            Random globalRand = new Random(42);
            HashSet<int> usedHues = new HashSet<int>();

            foreach (char c in uniqueChars)
            {
                if (c == 'S' || c == '_')
                {
                    globalColorMap[c] = Color.Gray;
                    continue;
                }

                int hue;
                int attempts = 0;
                do
                {
                    hue = globalRand.Next(0, 360);
                    attempts++;
                }
                while (usedHues.Any(h => Math.Abs(h - hue) < 30) && attempts < 100);

                usedHues.Add(hue);
                globalColorMap[c] = FromHsv(hue, 0.85, 0.9);
            }

            int totalPatterns = patternList.Count;

            // ========== CORRECTED SORTING SECTION ==========

            // --- SORT FACES WITHIN EACH STRIP ---
            var stripGroups = faceMetadataList
                .Where(f => f.StripId >= 0) // Exclude orphan faces
                .GroupBy(f => f.StripId)
                .ToList();

            // Sort each strip using hybrid method
            foreach (var group in stripGroups)
            {
                int stripId = group.Key;
                var facesInStrip = group.ToList();

                // Get center curve for this strip (CORRECTED)
                Curve centerCurve = null;
                if (stripCenterCurves != null &&
                    stripCenterCurves.TryGetValue(stripId, out var curveList) &&
                    curveList != null &&
                    curveList.Count > 0)
                {
                    centerCurve = curveList[0]; // Get first curve from the List<Curve>
                }

                // Use hybrid sorting
                var sortedFaces = SortFacesAlongStrip(mesh, facesInStrip, centerCurve);

                // Update local indices
                for (int i = 0; i < sortedFaces.Count; i++)
                {
                    sortedFaces[i].LocalIndex = i;
                }
              }

            // --- ORDER STRIPS BY CREATION TIME ---
            var orderedStripGroups = stripGroups
                .OrderBy(g => g.Min(f => f.CreationOrder))
                .ToList();

            // --- SEQUENTIAL PATTERN ASSIGNMENT ---
            int globalStripCounter = 0;

            foreach (var stripGroup in orderedStripGroups)
            {
                int stripId = stripGroup.Key;
                var facesInStrip = stripGroup.OrderBy(f => f.LocalIndex).ToList();

                // --- ASSIGN PATTERN ---
                int patternIndex = globalStripCounter % totalPatterns;
                globalStripCounter++;

                string rawPattern = patternList[patternIndex];
                bool repeatPattern = rawPattern.StartsWith("__") && rawPattern.EndsWith("__");
                string basePattern = GetBasePattern(rawPattern);

                int faceCount = facesInStrip.Count;
                string finalPattern = GetFinalPattern(basePattern, faceCount, repeatPattern);

                // --- MAP PATTERN TO SUB-BRANCHES ---
                var patternMap = new Dictionary<char, int>();
                int subBranchCounter = 0;

                for (int i = 0; i < faceCount; i++)
                {
                    var faceMeta = facesInStrip[i];
                    char p = (i < finalPattern.Length) ? finalPattern[i] : 'A';

                    if (!patternMap.ContainsKey(p))
                    {
                        patternMap[p] = subBranchCounter;
                        subBranchCounter++;
                    }

                    faceMeta.FinalPath = new GH_Path(stripId, patternMap[p]);
                    faceMeta.PatternChar = p;
                }
            }

            // --- PARALLEL BLOCK CREATION ---
            var results = new ConcurrentBag<KeyValuePair<GH_Path, Mesh>>();

            Parallel.ForEach(faceMetadataList.Where(f => f.FinalPath != null), faceMeta =>
            {
                MeshFace face = mesh.Faces[faceMeta.FaceIndex];
                Mesh block = new Mesh();

                if (face.IsQuad)
                    BuildQuadBlock(mesh, face, block, thickness);
                else if (face.IsTriangle)
                    BuildTriangleBlock(mesh, face, block, thickness);

                if (globalColorMap.TryGetValue(faceMeta.PatternChar, out Color color))
                    block.VertexColors.CreateMonotoneMesh(color);

                block.Normals.ComputeNormals();
                block.Compact();

                results.Add(new KeyValuePair<GH_Path, Mesh>(faceMeta.FinalPath, block));
            });

            // --- OUTPUT ---
            var sortedResults = results.OrderBy(x => x.Key);
            foreach (var kvp in sortedResults)
                blocksTree.Add(kvp.Value, kvp.Key);
        }


        private string GetBasePattern(string rawPattern)
        {
            if (string.IsNullOrEmpty(rawPattern)) return "A";
            if (rawPattern.Length > 4 && rawPattern.StartsWith("__") && rawPattern.EndsWith("__"))
                return rawPattern.Substring(2, rawPattern.Length - 4);
            return rawPattern;
        }

        private string GetFinalPattern(string basePattern, int count, bool repeat)
        {
            System.Text.StringBuilder sb = new System.Text.StringBuilder(count);
            for (int i = 0; i < count; i++)
            {
                if (repeat) sb.Append(basePattern[i % basePattern.Length]);
                else sb.Append(i < basePattern.Length ? basePattern[i] : basePattern[basePattern.Length - 1]);
            }
            return sb.ToString();
        }

        private Color FromHsv(double h, double s, double v)
        {
            double c = v * s; double x = c * (1 - Math.Abs((h / 60) % 2 - 1)); double m = v - c; double r = 0, g = 0, b = 0;
            if (h < 60) { r = c; g = x; } else if (h < 120) { r = x; g = c; } else if (h < 180) { g = c; b = x; } else if (h < 240) { g = x; b = c; } else if (h < 300) { r = x; b = c; } else { r = c; b = x; }
            return Color.FromArgb(255, (int)((r + m) * 255), (int)((g + m) * 255), (int)((b + m) * 255));
        }

        private void BuildQuadBlock(Mesh mesh, MeshFace face, Mesh block, double thickness)
        {
            int vA = face.A; int vB = face.B; int vC = face.C; int vD = face.D;
            Point3d pA = (Point3d)mesh.Vertices[vA]; Point3d pB = (Point3d)mesh.Vertices[vB];
            Point3d pC = (Point3d)mesh.Vertices[vC]; Point3d pD = (Point3d)mesh.Vertices[vD];
            Vector3d nA = new Vector3d(mesh.Normals[vA]) * thickness; Vector3d nB = new Vector3d(mesh.Normals[vB]) * thickness;
            Vector3d nC = new Vector3d(mesh.Normals[vC]) * thickness; Vector3d nD = new Vector3d(mesh.Normals[vD]) * thickness;
            block.Vertices.Add(pA); block.Vertices.Add(pB); block.Vertices.Add(pC); block.Vertices.Add(pD);
            block.Vertices.Add(pA + nA); block.Vertices.Add(pB + nB); block.Vertices.Add(pC + nC); block.Vertices.Add(pD + nD);
            block.Faces.AddFace(0, 1, 2, 3); block.Faces.AddFace(4, 7, 6, 5);
            block.Faces.AddFace(0, 1, 5, 4); block.Faces.AddFace(1, 2, 6, 5);
            block.Faces.AddFace(2, 3, 7, 6); block.Faces.AddFace(3, 0, 4, 7);
        }

        private void BuildTriangleBlock(Mesh mesh, MeshFace face, Mesh block, double thickness)
        {
            int vA = face.A; int vB = face.B; int vC = face.C;
            Point3d pA = (Point3d)mesh.Vertices[vA]; Point3d pB = (Point3d)mesh.Vertices[vB]; Point3d pC = (Point3d)mesh.Vertices[vC];
            Vector3d nA = new Vector3d(mesh.Normals[vA]) * thickness; Vector3d nB = new Vector3d(mesh.Normals[vB]) * thickness;
            Vector3d nC = new Vector3d(mesh.Normals[vC]) * thickness;
            block.Vertices.Add(pA); block.Vertices.Add(pB); block.Vertices.Add(pC);
            block.Vertices.Add(pA + nA); block.Vertices.Add(pB + nB); block.Vertices.Add(pC + nC);
            block.Faces.AddFace(0, 1, 2); block.Faces.AddFace(3, 5, 4);
            block.Faces.AddFace(0, 1, 4, 3); block.Faces.AddFace(1, 2, 5, 4);
            block.Faces.AddFace(2, 0, 3, 5);
        }

        public static Bitmap ConvertByteArrayToBitmap(byte[] imageData) { using (MemoryStream ms = new MemoryStream(imageData)) { return new Bitmap(ms); } }
        protected override Bitmap Icon => ConvertByteArrayToBitmap(PoliBrick.Properties.Resources.Solver);
    }

    class Point3dComparer : IEqualityComparer<Point3d>
    {
        public bool Equals(Point3d x, Point3d y)
        {
            return x.DistanceTo(y) < 0.001;
        }

        public int GetHashCode(Point3d obj)
        {
            // Round coordinates to ensure nearby points have the same hash
            int x = (int)Math.Round(obj.X * 1000);
            int y = (int)Math.Round(obj.Y * 1000);
            int z = (int)Math.Round(obj.Z * 1000);
            return x ^ y ^ z;
        }
    }
}