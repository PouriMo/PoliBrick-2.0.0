using Grasshopper;
using Grasshopper.Kernel;
using Grasshopper.Kernel.Data;
using Grasshopper.Kernel.Types;
using Rhino;
using Rhino.Geometry;
using Rhino.Geometry.Intersect;
using System;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.Drawing;
using System.Linq;
using System.Threading.Tasks;

namespace PoliBrick.Util
{
    public class CullIntersection : GH_Component
    {
        public CullIntersection()
          : base("6.CullIntersection", "Cull",
              "Refine brick intersections using pre-calculated planes from Solver.",
              "PoliBrick", "Utilities")
        {
        }

        public override GH_Exposure Exposure => GH_Exposure.primary;

        protected override void RegisterInputParams(GH_InputParamManager pManager)
        {
            pManager.AddMeshParameter("Bricks", "B", "Input Bricks (Ordered Tree from Solver)", GH_ParamAccess.tree);
            pManager.AddPlaneParameter("FacePlanes", "P", "Input Face Planes (Output FP from Solver)", GH_ParamAccess.tree);
            pManager.AddTextParameter("FrontPattern", "FP", "Horizontal Pattern (< = Cut Current, > = Cut Next, % = Average)", GH_ParamAccess.list);
            pManager.AddTextParameter("TopPattern", "TP", "Vertical Pattern", GH_ParamAccess.list);
            pManager.AddNumberParameter("Offset", "O", "Mortar Gap", GH_ParamAccess.item, 0.0);
        }

        protected override void RegisterOutputParams(GH_OutputParamManager pManager)
        {
            pManager.AddMeshParameter("RefinedBricks", "R", "Refined meshes", GH_ParamAccess.tree);
        }

        class BrickData
        {
            public Mesh Geometry;
            public Color BaseColor;
            public int CurveIndex;
            public int BlockIndex;
            public int GlobalID;
            public GH_Path OriginalPath;
            public GH_Path PlanePath;
            public Point3d BBoxCenter;

            public BrickData(Mesh m, int cIdx, int bIdx, GH_Path origPath, int id)
            {
                Geometry = m;
                BaseColor = (m.VertexColors.Count > 0) ? m.VertexColors[0] : Color.Gray;
                CurveIndex = cIdx;
                BlockIndex = bIdx;
                OriginalPath = origPath;
                GlobalID = id;
                PlanePath = new GH_Path(cIdx).AppendElement(bIdx);
                BBoxCenter = m.GetBoundingBox(true).Center;
            }
        }

        struct CutTask
        {
            public int TargetBrickID;
            public Plane CutPlane;
            public Point3d KeepPoint;
        }

        protected override void SolveInstance(IGH_DataAccess DA)
        {
            if (!DA.GetDataTree(0, out GH_Structure<GH_Mesh> bricksTree)) return;
            if (!DA.GetDataTree(1, out GH_Structure<GH_Plane> planesTree)) return;

            List<string> frontPat = new List<string>();
            List<string> topPat = new List<string>();
            if (!DA.GetDataList(2, frontPat)) return;
            if (!DA.GetDataList(3, topPat)) return;

            double offset = 0;
            if (!DA.GetData(4, ref offset)) return;

            List<BrickData> database = new List<BrickData>();
            int idCounter = 0;

            // Track block counts per curve for closed loop detection
            Dictionary<int, int> blocksPerCurve = new Dictionary<int, int>();

            foreach (GH_Path path in bricksTree.Paths)
            {
                int cIdx = path.Indices[0];
                var branch = bricksTree.get_Branch(path);

                for (int bIdx = 0; bIdx < branch.Count; bIdx++)
                {
                    GH_Mesh ghMesh = branch[bIdx] as GH_Mesh;
                    if (ghMesh?.Value == null) continue;

                    Mesh m = ghMesh.Value.DuplicateMesh();
                    database.Add(new BrickData(m, cIdx, bIdx, path, idCounter++));

                    // Track max block index per curve
                    if (!blocksPerCurve.ContainsKey(cIdx))
                        blocksPerCurve[cIdx] = 0;
                    blocksPerCurve[cIdx] = Math.Max(blocksPerCurve[cIdx], bIdx + 1);
                }
            }

            RTree rTree = new RTree();
            for (int i = 0; i < database.Count; i++)
                rTree.Insert(database[i].Geometry.GetBoundingBox(true), i);

            ConcurrentBag<CutTask> tasks = new ConcurrentBag<CutTask>();

            Parallel.For(0, database.Count, i =>
            {
                BrickData bA = database[i];
                BoundingBox searchBox = bA.Geometry.GetBoundingBox(true);
                searchBox.Inflate(0.01);

                List<int> candidates = new List<int>();
                rTree.Search(searchBox, (s, e) => { if (e.Id > i) candidates.Add(e.Id); });

                foreach (int j in candidates)
                {
                    BrickData bB = database[j];
                    if (CheckCollision(bA.Geometry, bB.Geometry))
                    {
                        ResolveCollision(bA, bB, planesTree, frontPat, topPat, offset, tasks, blocksPerCurve);
                    }
                }
            });

            var groupedTasks = tasks.GroupBy(t => t.TargetBrickID);
            Parallel.ForEach(groupedTasks, group =>
            {
                BrickData targetBrickData = database[group.Key];
                Mesh targetMesh = targetBrickData.Geometry.DuplicateMesh();
                Color color = targetBrickData.BaseColor;

                foreach (var task in group)
                {
                    targetMesh = CutMeshClean(targetMesh, task.CutPlane, color, task.KeepPoint);
                    if (targetMesh == null) break;
                }
                database[group.Key].Geometry = targetMesh;
            });

            DataTree<GH_Mesh> outTree = new DataTree<GH_Mesh>();
            foreach (var b in database)
            {
                if (b.Geometry != null && b.Geometry.Vertices.Count > 0 && b.Geometry.Faces.Count > 0)
                {
                    outTree.Add(new GH_Mesh(b.Geometry), b.OriginalPath);
                }
            }
            DA.SetDataTree(0, outTree);
        }

        private bool CheckCollision(Mesh a, Mesh b)
        {
            if (!BoundingBox.Intersection(a.GetBoundingBox(true), b.GetBoundingBox(true)).IsValid)
                return false;
            try
            {
                Mesh[] overlaps = Mesh.CreateBooleanIntersection(new List<Mesh> { a }, new List<Mesh> { b });
                return overlaps != null && overlaps.Length > 0;
            }
            catch { return false; }
        }

        private void ResolveCollision(BrickData bA, BrickData bB, GH_Structure<GH_Plane> planesTree,
            List<string> frontPat, List<string> topPat, double offset, ConcurrentBag<CutTask> tasks,
            Dictionary<int, int> blocksPerCurve)
        {
            bool sameCurve = bA.CurveIndex == bB.CurveIndex;

            BrickData current, next;
            bool isClosedLoopConnection = false;

            if (sameCurve)
            {
                int blockCount = blocksPerCurve.ContainsKey(bA.CurveIndex) ? blocksPerCurve[bA.CurveIndex] : 0;

                // Detect closed loop connection: first block (0) meets last block (blockCount-1)
                int minIdx = Math.Min(bA.BlockIndex, bB.BlockIndex);
                int maxIdx = Math.Max(bA.BlockIndex, bB.BlockIndex);

                // It's a closed loop connection if one is index 0 and the other is the last index
                // AND there are at least 3 blocks (otherwise it's just two adjacent blocks)
                isClosedLoopConnection = (minIdx == 0 && maxIdx == blockCount - 1 && blockCount >= 3);

                if (isClosedLoopConnection)
                {
                    // In closed loop: the LAST block's front face meets the FIRST block's back face
                    // So the last block is "current" (contributing front face)
                    // And the first block is "next" (contributing back face)
                    current = (bA.BlockIndex > bB.BlockIndex) ? bA : bB;
                    next = (bA.BlockIndex > bB.BlockIndex) ? bB : bA;
                }
                else
                {
                    // Normal sequential case: lower index is current
                    current = (bA.BlockIndex < bB.BlockIndex) ? bA : bB;
                    next = (bA.BlockIndex < bB.BlockIndex) ? bB : bA;
                }
            }
            else
            {
                current = (bA.CurveIndex < bB.CurveIndex) ? bA : bB;
                next = (bA.CurveIndex < bB.CurveIndex) ? bB : bA;
            }

            string symbol = "";
            if (sameCurve)
            {
                if (frontPat.Count == 0) return;
                symbol = frontPat[current.BlockIndex % frontPat.Count];
            }
            else
            {
                if (topPat.Count == 0) return;
                symbol = topPat[current.BlockIndex % topPat.Count];
            }

            var currentPlanes = planesTree.get_Branch(current.PlanePath) as List<GH_Plane>;
            var nextPlanes = planesTree.get_Branch(next.PlanePath) as List<GH_Plane>;

            if (currentPlanes == null || nextPlanes == null || currentPlanes.Count < 4 || nextPlanes.Count < 4)
                return;

            // Plane indices: 0=Front, 1=Back, 2=Top, 3=Bottom

            if (sameCurve)
            {
                // Horizontal intersection (along curve)
                // Current's FRONT face meets Next's BACK face
                Plane currentFront = currentPlanes[0].Value;  // Front of current
                Plane nextBack = nextPlanes[1].Value;          // Back of next

                // Anchors: points we want to KEEP
                // For current: keep the back portion
                // For next: keep the front portion
                Point3d currentAnchor = currentPlanes[1].Value.Origin;  // Back of current
                Point3d nextAnchor = nextPlanes[0].Value.Origin;        // Front of next

                if (symbol == ">")
                {
                    // Cut next block using current's front plane
                    Plane knife = currentFront;
                    // Move plane away from next (towards current) to create gap
                    knife.Translate(-knife.Normal * offset * 0.5);
                    tasks.Add(new CutTask { TargetBrickID = next.GlobalID, CutPlane = knife, KeepPoint = nextAnchor });
                }
                else if (symbol == "<")
                {
                    // Cut current block using next's back plane
                    Plane knife = nextBack;
                    // Move plane away from current (towards next) to create gap
                    knife.Translate(-knife.Normal * offset * 0.5);
                    tasks.Add(new CutTask { TargetBrickID = current.GlobalID, CutPlane = knife, KeepPoint = currentAnchor });
                }
                else if (symbol == "%")
                {
                    // Cut both using average plane
                    Plane knife = CalculateAveragePlane(currentFront, nextBack);

                    // Create two cutting planes offset from average
                    // knifeA cuts current: move TOWARDS current (away from intersection) to create gap
                    // knifeB cuts next: move TOWARDS next (away from intersection) to create gap
                    Plane knifeA = knife;
                    Plane knifeB = knife;

                    // Knife normal points from current towards next
                    // For current: move plane in NEGATIVE normal direction (towards current's back)
                    // For next: move plane in POSITIVE normal direction (towards next's front)
                    knifeA.Translate(-knife.Normal * offset * 0.5);
                    knifeB.Translate(knife.Normal * offset * 0.5);

                    tasks.Add(new CutTask { TargetBrickID = current.GlobalID, CutPlane = knifeA, KeepPoint = currentAnchor });
                    tasks.Add(new CutTask { TargetBrickID = next.GlobalID, CutPlane = knifeB, KeepPoint = nextAnchor });
                }
            }
            else
            {
                // Vertical intersection (between curves/rows)
                // Current's TOP face meets Next's BOTTOM face
                Plane currentTop = currentPlanes[2].Value;     // Top of current
                Plane nextBottom = nextPlanes[3].Value;        // Bottom of next

                // Anchors: points we want to KEEP
                // For current: keep the bottom portion
                // For next: keep the top portion
                Point3d currentAnchor = currentPlanes[3].Value.Origin;  // Bottom of current
                Point3d nextAnchor = nextPlanes[2].Value.Origin;        // Top of next

                if (symbol == ">")
                {
                    // Cut next block using current's top plane
                    Plane knife = currentTop;
                    knife.Translate(-knife.Normal * offset * 0.5);
                    tasks.Add(new CutTask { TargetBrickID = next.GlobalID, CutPlane = knife, KeepPoint = nextAnchor });
                }
                else if (symbol == "<")
                {
                    // Cut current block using next's bottom plane
                    Plane knife = nextBottom;
                    knife.Translate(-knife.Normal * offset * 0.5);
                    tasks.Add(new CutTask { TargetBrickID = current.GlobalID, CutPlane = knife, KeepPoint = currentAnchor });
                }
                else if (symbol == "%")
                {
                    // Cut both using average plane
                    Plane knife = CalculateAveragePlane(currentTop, nextBottom);

                    Plane knifeA = knife;
                    Plane knifeB = knife;

                    // Knife normal points from current towards next (upward typically)
                    knifeA.Translate(-knife.Normal * offset * 0.5);
                    knifeB.Translate(knife.Normal * offset * 0.5);

                    tasks.Add(new CutTask { TargetBrickID = current.GlobalID, CutPlane = knifeA, KeepPoint = currentAnchor });
                    tasks.Add(new CutTask { TargetBrickID = next.GlobalID, CutPlane = knifeB, KeepPoint = nextAnchor });
                }
            }
        }

        /// <summary>
        /// Calculates average plane between two faces with consistent normal direction.
        /// Normal points from p1 (current's outward face) towards p2 (next's outward face).
        /// </summary>
        private Plane CalculateAveragePlane(Plane p1, Plane p2)
        {
            // p1: outward face of current block (e.g., front face - normal points forward)
            // p2: outward face of next block (e.g., back face - normal points backward)

            // Calculate origin at midpoint
            Point3d newOrigin;

            // Check if planes are nearly parallel
            double normalDot = p1.Normal * p2.Normal;

            if (Math.Abs(normalDot) > 0.99)
            {
                // Nearly parallel planes - just average the origins
                newOrigin = (p1.Origin + p2.Origin) / 2.0;
            }
            else
            {
                // Non-parallel planes - find intersection line and place origin on it
                if (Intersection.PlanePlane(p1, p2, out Line intersectLine))
                {
                    Point3d midOrigins = (p1.Origin + p2.Origin) / 2.0;
                    double t = intersectLine.ClosestParameter(midOrigins);
                    newOrigin = intersectLine.PointAt(t);
                }
                else
                {
                    newOrigin = (p1.Origin + p2.Origin) / 2.0;
                }
            }

            // Calculate normal direction
            // p1.Normal points outward from current (towards next)
            // p2.Normal points outward from next (towards current, i.e., opposite direction)
            // We want the average plane normal to point from current towards next

            // Since p2.Normal points opposite, we subtract it (or equivalently, add -p2.Normal)
            Vector3d avgNormal = p1.Normal - p2.Normal;

            if (avgNormal.IsTiny(0.001))
            {
                // Normals are nearly opposite and cancel out - use p1's normal
                avgNormal = p1.Normal;
            }

            avgNormal.Unitize();

            // Create plane with this normal
            // Find a reasonable X axis by averaging
            Vector3d avgX = (p1.XAxis + p2.XAxis) / 2.0;
            if (avgX.IsTiny(0.001))
            {
                avgX = p1.XAxis;
            }

            // Make X perpendicular to normal
            avgX = avgX - (avgX * avgNormal) * avgNormal;
            avgX.Unitize();

            // Y axis is cross product
            Vector3d avgY = Vector3d.CrossProduct(avgNormal, avgX);
            avgY.Unitize();

            // Construct plane (XAxis, YAxis) where Normal = X × Y
            // Since we computed Y = Normal × X, we have Normal = X × Y by right-hand rule
            return new Plane(newOrigin, avgX, avgY);
        }

        // ==================== CLEAN MESH CUTTING - CORNERS ONLY ====================

        private Mesh CutMeshClean(Mesh mesh, Plane cutter, Color color, Point3d keepPoint)
        {
            if (mesh == null || !mesh.IsValid) return mesh;

            const double tolerance = 0.0001;

            mesh = mesh.DuplicateMesh();
            mesh.Vertices.CombineIdentical(true, true);
            mesh.Weld(Math.PI);

            bool keepPositive = cutter.DistanceTo(keepPoint) >= 0;
            int keepSign = keepPositive ? 1 : -1;

            int vertexCount = mesh.Vertices.Count;
            int[] vertexSide = new int[vertexCount];

            for (int i = 0; i < vertexCount; i++)
            {
                double dist = cutter.DistanceTo(mesh.Vertices[i]);
                if (Math.Abs(dist) < tolerance)
                    vertexSide[i] = 0;
                else
                    vertexSide[i] = dist > 0 ? 1 : -1;
            }

            bool hasKept = vertexSide.Any(s => s == keepSign || s == 0);
            bool hasCut = vertexSide.Any(s => s == -keepSign);

            if (!hasCut) return mesh;
            if (!hasKept) return null;

            Mesh result = new Mesh();
            Dictionary<int, int> oldToNewVertex = new Dictionary<int, int>();

            for (int i = 0; i < vertexCount; i++)
            {
                if (vertexSide[i] == keepSign || vertexSide[i] == 0)
                {
                    oldToNewVertex[i] = result.Vertices.Count;
                    result.Vertices.Add(mesh.Vertices[i]);
                }
            }

            Dictionary<long, int> edgeIntersectionVertex = new Dictionary<long, int>();
            List<int> capVertices = new List<int>();

            foreach (MeshFace face in mesh.Faces)
            {
                int[] faceVerts = face.IsQuad
                    ? new[] { face.A, face.B, face.C, face.D }
                    : new[] { face.A, face.B, face.C };

                int[] sides = faceVerts.Select(v => vertexSide[v]).ToArray();

                if (sides.All(s => s == -keepSign))
                    continue;

                if (sides.All(s => s == keepSign || s == 0))
                {
                    int[] newVerts = faceVerts.Select(v => oldToNewVertex[v]).ToArray();
                    AddFaceToMesh(result, newVerts);

                    for (int i = 0; i < faceVerts.Length; i++)
                    {
                        if (sides[i] == 0)
                        {
                            int mappedVert = oldToNewVertex[faceVerts[i]];
                            if (!capVertices.Contains(mappedVert))
                                capVertices.Add(mappedVert);
                        }
                    }
                    continue;
                }

                List<int> newFaceVerts = new List<int>();

                for (int i = 0; i < faceVerts.Length; i++)
                {
                    int currVert = faceVerts[i];
                    int nextVert = faceVerts[(i + 1) % faceVerts.Length];
                    int currSide = sides[i];
                    int nextSide = sides[(i + 1) % faceVerts.Length];

                    if (currSide == keepSign || currSide == 0)
                    {
                        int mappedVert = oldToNewVertex[currVert];
                        if (!newFaceVerts.Contains(mappedVert))
                            newFaceVerts.Add(mappedVert);

                        if (currSide == 0 && !capVertices.Contains(mappedVert))
                            capVertices.Add(mappedVert);
                    }

                    if ((currSide == keepSign && nextSide == -keepSign) ||
                        (currSide == -keepSign && nextSide == keepSign))
                    {
                        long edgeKey = MakeEdgeKey(currVert, nextVert);

                        if (!edgeIntersectionVertex.ContainsKey(edgeKey))
                        {
                            Line edge = new Line(
                                new Point3d(mesh.Vertices[currVert]),
                                new Point3d(mesh.Vertices[nextVert]));

                            if (Intersection.LinePlane(edge, cutter, out double t))
                            {
                                int newVertIdx = result.Vertices.Count;
                                result.Vertices.Add(edge.PointAt(t));
                                edgeIntersectionVertex[edgeKey] = newVertIdx;
                            }
                        }

                        if (edgeIntersectionVertex.TryGetValue(edgeKey, out int intersectVert))
                        {
                            if (!newFaceVerts.Contains(intersectVert))
                                newFaceVerts.Add(intersectVert);

                            if (!capVertices.Contains(intersectVert))
                                capVertices.Add(intersectVert);
                        }
                    }
                }

                if (newFaceVerts.Count >= 3)
                {
                    CreateFacesFromPolygon(result, newFaceVerts);
                }
            }

            if (capVertices.Count >= 3)
            {
                CreateCapFromVertices(result, capVertices, cutter, keepPositive);
            }

            result.Vertices.CombineIdentical(true, true);
            result.Faces.CullDegenerateFaces();
            result.Compact();
            result.Normals.ComputeNormals();
            result.UnifyNormals();

            if (!result.IsClosed)
            {
                result = TryCloseMesh(result, cutter, keepPositive);
            }

            result.VertexColors.CreateMonotoneMesh(color);

            return result.Faces.Count > 0 ? result : null;
        }

        private long MakeEdgeKey(int v1, int v2)
        {
            int min = Math.Min(v1, v2);
            int max = Math.Max(v1, v2);
            return ((long)min << 32) | (uint)max;
        }

        private void AddFaceToMesh(Mesh mesh, int[] verts)
        {
            if (verts.Length == 3)
                mesh.Faces.AddFace(verts[0], verts[1], verts[2]);
            else if (verts.Length == 4)
                mesh.Faces.AddFace(verts[0], verts[1], verts[2], verts[3]);
        }

        private void CreateFacesFromPolygon(Mesh mesh, List<int> vertIndices)
        {
            int n = vertIndices.Count;
            if (n < 3) return;

            if (n == 3)
            {
                mesh.Faces.AddFace(vertIndices[0], vertIndices[1], vertIndices[2]);
                return;
            }

            if (n == 4)
            {
                Point3d[] pts = vertIndices.Select(i => new Point3d(mesh.Vertices[i])).ToArray();

                if (IsValidConvexQuad(pts))
                {
                    mesh.Faces.AddFace(vertIndices[0], vertIndices[1], vertIndices[2], vertIndices[3]);
                }
                else
                {
                    int[] bestSplit = GetBestQuadDiagonal(pts);
                    mesh.Faces.AddFace(vertIndices[bestSplit[0]], vertIndices[bestSplit[1]], vertIndices[bestSplit[2]]);
                    mesh.Faces.AddFace(vertIndices[bestSplit[0]], vertIndices[bestSplit[2]], vertIndices[bestSplit[3]]);
                }
                return;
            }

            Point3d[] points = vertIndices.Select(i => new Point3d(mesh.Vertices[i])).ToArray();
            Point3d centroid = new Point3d(
                points.Average(p => p.X),
                points.Average(p => p.Y),
                points.Average(p => p.Z));

            int pivotIdx = 0;
            double minDist = double.MaxValue;
            for (int i = 0; i < n; i++)
            {
                double d = points[i].DistanceTo(centroid);
                if (d < minDist)
                {
                    minDist = d;
                    pivotIdx = i;
                }
            }

            List<int> reordered = new List<int>();
            for (int i = 0; i < n; i++)
                reordered.Add(vertIndices[(pivotIdx + i) % n]);

            for (int i = 1; i < n - 1; i++)
            {
                mesh.Faces.AddFace(reordered[0], reordered[i], reordered[i + 1]);
            }
        }

        private bool IsValidConvexQuad(Point3d[] pts)
        {
            if (pts.Length != 4) return false;

            if (Plane.FitPlaneToPoints(pts.ToList(), out Plane plane) != PlaneFitResult.Success)
                return false;

            const double planeTol = 0.01;
            foreach (var pt in pts)
                if (Math.Abs(plane.DistanceTo(pt)) > planeTol)
                    return false;

            Transform toPlane = Transform.PlaneToPlane(plane, Plane.WorldXY);
            Point2d[] pts2d = new Point2d[4];
            for (int i = 0; i < 4; i++)
            {
                Point3d p = new Point3d(pts[i]);
                p.Transform(toPlane);
                pts2d[i] = new Point2d(p.X, p.Y);
            }

            double? expectedSign = null;
            for (int i = 0; i < 4; i++)
            {
                Point2d a = pts2d[i];
                Point2d b = pts2d[(i + 1) % 4];
                Point2d c = pts2d[(i + 2) % 4];

                double cross = (b.X - a.X) * (c.Y - b.Y) - (b.Y - a.Y) * (c.X - b.X);

                if (Math.Abs(cross) < 1e-10) continue;

                if (expectedSign == null)
                    expectedSign = Math.Sign(cross);
                else if (Math.Sign(cross) != expectedSign)
                    return false;
            }

            return true;
        }

        private int[] GetBestQuadDiagonal(Point3d[] pts)
        {
            double area1_02 = TriangleArea(pts[0], pts[1], pts[2]);
            double area2_02 = TriangleArea(pts[0], pts[2], pts[3]);
            double ratio_02 = Math.Min(area1_02, area2_02) / (Math.Max(area1_02, area2_02) + 1e-10);

            double area1_13 = TriangleArea(pts[1], pts[2], pts[3]);
            double area2_13 = TriangleArea(pts[1], pts[3], pts[0]);
            double ratio_13 = Math.Min(area1_13, area2_13) / (Math.Max(area1_13, area2_13) + 1e-10);

            double aspect_02 = Math.Min(GetTriangleAspectRatio(pts[0], pts[1], pts[2]),
                                        GetTriangleAspectRatio(pts[0], pts[2], pts[3]));
            double aspect_13 = Math.Min(GetTriangleAspectRatio(pts[1], pts[2], pts[3]),
                                        GetTriangleAspectRatio(pts[1], pts[3], pts[0]));

            double score_02 = ratio_02 * aspect_02;
            double score_13 = ratio_13 * aspect_13;

            if (score_02 >= score_13)
                return new[] { 0, 1, 2, 3 };
            else
                return new[] { 1, 2, 3, 0 };
        }

        private double TriangleArea(Point3d a, Point3d b, Point3d c)
        {
            Vector3d v1 = b - a;
            Vector3d v2 = c - a;
            return Vector3d.CrossProduct(v1, v2).Length * 0.5;
        }

        private double GetTriangleAspectRatio(Point3d a, Point3d b, Point3d c)
        {
            double ab = a.DistanceTo(b);
            double bc = b.DistanceTo(c);
            double ca = c.DistanceTo(a);

            double s = (ab + bc + ca) / 2;
            double area = Math.Sqrt(Math.Max(0, s * (s - ab) * (s - bc) * (s - ca)));

            if (area < 1e-10) return 0;

            double longest = Math.Max(ab, Math.Max(bc, ca));
            double inradius = area / s;

            return inradius / longest * 2;
        }

        private void CreateCapFromVertices(Mesh mesh, List<int> capVerts, Plane cutter, bool keepPositive)
        {
            List<int> uniqueVerts = capVerts.Distinct().ToList();

            List<int> onPlaneVerts = new List<int>();
            foreach (int idx in uniqueVerts)
            {
                if (Math.Abs(cutter.DistanceTo(mesh.Vertices[idx])) < 0.01)
                    onPlaneVerts.Add(idx);
            }

            if (onPlaneVerts.Count < 3) return;

            List<Point3d> points = onPlaneVerts.Select(i => new Point3d(mesh.Vertices[i])).ToList();
            Point3d centroid = new Point3d(
                points.Average(p => p.X),
                points.Average(p => p.Y),
                points.Average(p => p.Z));

            var sorted = onPlaneVerts
                .Select((idx, i) => new { Index = idx, Point = points[i] })
                .OrderBy(x =>
                {
                    Vector3d v = x.Point - centroid;
                    return Math.Atan2(
                        Vector3d.CrossProduct(cutter.XAxis, v) * cutter.Normal,
                        v * cutter.XAxis);
                })
                .ToList();

            List<int> orderedVerts = sorted.Select(x => x.Index).ToList();
            int n = orderedVerts.Count;

            // Normal should point AWAY from the kept side (into the cut region)
            // If we kept positive side, cap normal should point negative (flip)
            // If we kept negative side, cap normal should point positive (no flip)
            bool flipNormal = keepPositive;

            if (n == 3)
            {
                if (flipNormal)
                    mesh.Faces.AddFace(orderedVerts[0], orderedVerts[2], orderedVerts[1]);
                else
                    mesh.Faces.AddFace(orderedVerts[0], orderedVerts[1], orderedVerts[2]);
            }
            else if (n == 4)
            {
                Point3d[] capPts = orderedVerts.Select(i => new Point3d(mesh.Vertices[i])).ToArray();

                if (IsValidConvexQuad(capPts))
                {
                    if (flipNormal)
                        mesh.Faces.AddFace(orderedVerts[0], orderedVerts[3], orderedVerts[2], orderedVerts[1]);
                    else
                        mesh.Faces.AddFace(orderedVerts[0], orderedVerts[1], orderedVerts[2], orderedVerts[3]);
                }
                else
                {
                    int[] split = GetBestQuadDiagonal(capPts);
                    if (flipNormal)
                    {
                        mesh.Faces.AddFace(orderedVerts[split[2]], orderedVerts[split[1]], orderedVerts[split[0]]);
                        mesh.Faces.AddFace(orderedVerts[split[0]], orderedVerts[split[3]], orderedVerts[split[2]]);
                    }
                    else
                    {
                        mesh.Faces.AddFace(orderedVerts[split[0]], orderedVerts[split[1]], orderedVerts[split[2]]);
                        mesh.Faces.AddFace(orderedVerts[split[0]], orderedVerts[split[2]], orderedVerts[split[3]]);
                    }
                }
            }
            else
            {
                for (int i = 1; i < n - 1; i++)
                {
                    if (flipNormal)
                        mesh.Faces.AddFace(orderedVerts[0], orderedVerts[i + 1], orderedVerts[i]);
                    else
                        mesh.Faces.AddFace(orderedVerts[0], orderedVerts[i], orderedVerts[i + 1]);
                }
            }
        }

        private Mesh TryCloseMesh(Mesh mesh, Plane cutter, bool keepPositive)
        {
            var nakedEdges = mesh.GetNakedEdges();
            if (nakedEdges == null || nakedEdges.Length == 0)
                return mesh;

            foreach (var polyline in nakedEdges)
            {
                if (polyline == null || polyline.Count < 3)
                    continue;

                List<int> boundaryVerts = new List<int>();
                foreach (Point3d pt in polyline)
                {
                    int closestIdx = -1;
                    double minDist = double.MaxValue;
                    for (int i = 0; i < mesh.Vertices.Count; i++)
                    {
                        double d = pt.DistanceTo(mesh.Vertices[i]);
                        if (d < minDist)
                        {
                            minDist = d;
                            closestIdx = i;
                        }
                    }
                    if (closestIdx >= 0 && minDist < 0.01 && !boundaryVerts.Contains(closestIdx))
                    {
                        boundaryVerts.Add(closestIdx);
                    }
                }

                if (boundaryVerts.Count >= 3)
                {
                    bool onCuttingPlane = boundaryVerts.All(i =>
                        Math.Abs(cutter.DistanceTo(mesh.Vertices[i])) < 0.01);

                    if (onCuttingPlane)
                    {
                        CreateCapFromVertices(mesh, boundaryVerts, cutter, keepPositive);
                    }
                    else
                    {
                        CreateFacesFromPolygon(mesh, boundaryVerts);
                    }
                }
            }

            mesh.Vertices.CombineIdentical(true, true);
            mesh.Normals.ComputeNormals();
            mesh.UnifyNormals();

            return mesh;
        }

        protected override Bitmap Icon => ConvertByteArrayToBitmap(PoliBrick.Properties.Resources.Cullicon);
        public override Guid ComponentGuid => new Guid("01B09AB5-931F-4447-BF8F-E174519B6008");

        private static Bitmap ConvertByteArrayToBitmap(byte[] imageData)
        {
            using (var ms = new System.IO.MemoryStream(imageData))
                return new Bitmap(ms);
        }
    }
}