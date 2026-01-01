using System;
using System.Collections.Generic;
using System.Collections.Concurrent;
using System.Drawing;
using System.IO;
using System.Linq;
using System.Threading.Tasks;
using Grasshopper.Kernel;
using Rhino.Geometry;
using Rhino.Geometry.Intersect;

namespace PoliBrick.Util
{
    public class Cutter : GH_Component
    {
        private List<Curve> debugCurves = new List<Curve>();

        public Cutter()
          : base("5.Cutter", "Cutter",
              "Slices meshes at boundaries and filters them by checking if they sit on the valid surface interior.",
              "PoliBrick", "Utilities")
        {
        }

        public override GH_Exposure Exposure => GH_Exposure.primary;

        protected override void RegisterInputParams(GH_InputParamManager pManager)
        {
            pManager.AddBrepParameter("Surface", "S", "Brep Surface to use", GH_ParamAccess.item);
            pManager.AddMeshParameter("Bricks", "B", "Meshes to refine", GH_ParamAccess.list);
            pManager.AddIntegerParameter("Edge Indices", "I", "Edge indices", GH_ParamAccess.list);
            pManager.AddIntegerParameter("Search Spheres", "N", "Number of search spheres", GH_ParamAccess.item, 10);
        }

        protected override void RegisterOutputParams(GH_OutputParamManager pManager)
        {
            pManager.AddMeshParameter("Refined Blocks", "R", "Refined meshes", GH_ParamAccess.list);
        }

        protected override void SolveInstance(IGH_DataAccess DA)
        {
            Brep brep = null;
            List<Mesh> blocks = new List<Mesh>();
            List<int> indices = new List<int>();
            int searchSpheres = 10;

            if (!DA.GetData(0, ref brep)) return;
            if (!DA.GetDataList(1, blocks)) return;
            if (!DA.GetDataList(2, indices)) return;
            if (!DA.GetData(3, ref searchSpheres)) return;

            object refineBlocks = null;
            debugCurves.Clear();

            RunScript(brep, blocks, indices, searchSpheres, ref refineBlocks);

            DA.SetDataList(0, (List<Mesh>)refineBlocks);
        }

        private void RunScript(Brep brep, List<Mesh> blocks, List<int> indices, int searchSpheres, ref object refineBlocks)
        {
            if (brep == null || blocks == null || indices == null || indices.Count == 0) return;

            Point3d[] searchPts = GetCentroidPoints(blocks).ToArray();
            var rTree = new RTree();
            for (int i = 0; i < searchPts.Length; i++) rTree.Insert(searchPts[i], i);

            double[] edgeLengths = new double[brep.Edges.Count];
            for (int i = 0; i < brep.Edges.Count; i++) edgeLengths[i] = brep.Edges[i].GetLength();

            const int batchSize = 50;
            var batches = Enumerable.Range(0, (indices.Count + batchSize - 1) / batchSize)
                .Select(b => indices.Skip(b * batchSize).Take(batchSize).ToList()).ToList();

            Surface surfaceForCutting = brep.Faces[0];

            Parallel.ForEach(batches, currentBatch =>
            {
                foreach (int i in currentBatch)
                {
                    ProcessEdge(brep, blocks, searchPts, surfaceForCutting, rTree, edgeLengths, i, searchSpheres);
                }
            });

            List<Mesh> finalMeshes = FilterMeshesByTopology(blocks, brep);
            refineBlocks = finalMeshes;
        }

        private List<Mesh> FilterMeshesByTopology(List<Mesh> blocks, Brep brep)
        {
            var keepList = new ConcurrentBag<Mesh>();

            Parallel.ForEach(blocks, (mesh) =>
            {
                if (mesh == null || mesh.Vertices.Count == 0) return;

                var amp = AreaMassProperties.Compute(mesh);
                if (amp == null) return;

                Point3d centroid = amp.Centroid;
                Point3d closestPoint;
                ComponentIndex ci;
                double s, t;
                Vector3d normal;

                bool success = brep.ClosestPoint(centroid, out closestPoint, out ci, out s, out t, 0.0, out normal);

                if (success)
                {
                    bool keep = false;

                    if (ci.ComponentIndexType == ComponentIndexType.BrepFace)
                    {
                        BrepFace face = brep.Faces[ci.Index];
                        PointFaceRelation relation = face.IsPointOnFace(s, t);

                        if (relation == PointFaceRelation.Interior)
                        {
                            keep = true;
                        }
                        else if (relation == PointFaceRelation.Boundary)
                        {
                            if (centroid.DistanceTo(closestPoint) < 0.01) keep = true;
                        }
                    }

                    if (keep)
                    {
                        keepList.Add(mesh);
                    }
                }
            });

            return keepList.ToList();
        }

        private void ProcessEdge(Brep brep, List<Mesh> blocks, Point3d[] searchPts, Surface surface,
            RTree rTree, double[] edgeLengths, int index, int searchSpheres)
        {
            BrepEdge edge = brep.Edges[index % brep.Edges.Count];
            lock (debugCurves) { debugCurves.Add(edge); }

            List<Point3d> needles = GetEdgeNeedlePoints(edge, searchSpheres);
            List<int> nearbyIndices = new List<int>();
            double searchRadius = edgeLengths[index % edgeLengths.Length] / (2 * searchSpheres);

            foreach (var needle in needles)
                rTree.Search(new Sphere(needle, searchRadius), (sender, args) => nearbyIndices.Add(args.Id));

            nearbyIndices = nearbyIndices.Distinct().ToList();

            foreach (int j in nearbyIndices)
            {
                Mesh mesh = blocks[j];
                if (mesh == null) continue;

                double t, u, v;
                Point3d searchPt = searchPts[j];

                edge.ClosestPoint(searchPt, out t);
                Point3d closestPt = edge.PointAt(t);

                surface.ClosestPoint(closestPt, out u, out v);
                Vector3d normal = surface.NormalAt(u, v);

                Plane cuttingPlane = new Plane(closestPt, normal, edge.TangentAt(t));

                lock (blocks)
                {
                    CutMeshClean(blocks, mesh, cuttingPlane, j);
                }
            }
        }

        private List<Point3d> GetCentroidPoints(List<Mesh> meshes)
        {
            var centroids = new Point3d[meshes.Count];
            Parallel.For(0, meshes.Count, i =>
            {
                var amp = AreaMassProperties.Compute(meshes[i]);
                centroids[i] = amp != null ? amp.Centroid : meshes[i].GetBoundingBox(true).Center;
            });
            return centroids.ToList();
        }

        private List<Point3d> GetEdgeNeedlePoints(BrepEdge edge, int searchSpheres)
        {
            var needles = new Point3d[searchSpheres + 2];
            double min = edge.Domain.Min;
            double length = edge.Domain.Length;
            Parallel.For(0, searchSpheres, k =>
            {
                needles[k] = edge.PointAt(min + ((double)k / searchSpheres) * length);
            });
            needles[searchSpheres] = edge.PointAtStart;
            needles[searchSpheres + 1] = edge.PointAtEnd;
            return needles.ToList();
        }

        // ==================== CLEAN MESH CUTTING - CORNERS ONLY ====================

        /// <summary>
        /// Cuts mesh using only corner vertices - no mid-edge points created
        /// Replaces the old SplitAndModifyMesh method
        /// </summary>
        private void CutMeshClean(List<Mesh> meshes, Mesh mesh, Plane cutter, int index)
        {
            if (mesh == null || !mesh.IsValid) return;

            // Check if plane actually intersects the mesh
            if (!MeshPlaneIntersection(mesh, cutter)) return;

            // Get colors before cutting
            var originalColors = mesh.VertexColors.ToArray();
            bool hasColors = originalColors != null && originalColors.Length > 0;
            Color mainColor = hasColors ? GetAverageColor(originalColors) : Color.White;

            const double tolerance = 0.0001;

            // Prepare mesh
            mesh = mesh.DuplicateMesh();
            mesh.Vertices.CombineIdentical(true, true);
            mesh.Weld(Math.PI);

            // Determine which side to keep (positive side of plane - away from surface)
            // The normal points outward from surface, so we keep the positive side
            int keepSign = 1;

            // Classify all vertices
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

            // Check if cutting is needed
            bool hasKept = vertexSide.Any(s => s == keepSign || s == 0);
            bool hasCut = vertexSide.Any(s => s == -keepSign);

            if (!hasCut) return; // Nothing to cut
            if (!hasKept)
            {
                meshes[index] = null; // Everything would be cut
                return;
            }

            // Build new mesh
            Mesh result = new Mesh();

            // Map old vertex indices to new
            Dictionary<int, int> oldToNewVertex = new Dictionary<int, int>();

            // Add kept original vertices
            for (int i = 0; i < vertexCount; i++)
            {
                if (vertexSide[i] == keepSign || vertexSide[i] == 0)
                {
                    oldToNewVertex[i] = result.Vertices.Count;
                    result.Vertices.Add(mesh.Vertices[i]);
                }
            }

            // Track edge intersections
            Dictionary<long, int> edgeIntersectionVertex = new Dictionary<long, int>();

            // Track boundary vertices for cap
            List<int> capVertices = new List<int>();

            // Process each face
            foreach (MeshFace face in mesh.Faces)
            {
                int[] faceVerts = face.IsQuad
                    ? new[] { face.A, face.B, face.C, face.D }
                    : new[] { face.A, face.B, face.C };

                int[] sides = faceVerts.Select(v => vertexSide[v]).ToArray();

                // All on cut side - skip
                if (sides.All(s => s == -keepSign))
                    continue;

                // All on keep side or on plane - copy directly
                if (sides.All(s => s == keepSign || s == 0))
                {
                    int[] newVerts = faceVerts.Select(v => oldToNewVertex[v]).ToArray();
                    AddFaceToMesh(result, newVerts);

                    // Track vertices on the cutting plane for cap
                    for (int i = 0; i < faceVerts.Length; i++)
                    {
                        if (sides[i] == 0 && !capVertices.Contains(oldToNewVertex[faceVerts[i]]))
                            capVertices.Add(oldToNewVertex[faceVerts[i]]);
                    }
                    continue;
                }

                // Face is split
                List<int> newFaceVerts = new List<int>();

                for (int i = 0; i < faceVerts.Length; i++)
                {
                    int currVert = faceVerts[i];
                    int nextVert = faceVerts[(i + 1) % faceVerts.Length];
                    int currSide = sides[i];
                    int nextSide = sides[(i + 1) % faceVerts.Length];

                    // Add current vertex if on keep side or on plane
                    if (currSide == keepSign || currSide == 0)
                    {
                        int mappedVert = oldToNewVertex[currVert];
                        if (!newFaceVerts.Contains(mappedVert))
                            newFaceVerts.Add(mappedVert);

                        if (currSide == 0 && !capVertices.Contains(mappedVert))
                            capVertices.Add(mappedVert);
                    }

                    // If edge crosses plane, add intersection
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

                // Create faces from collected vertices
                if (newFaceVerts.Count >= 3)
                {
                    CreateFacesFromPolygon(result, newFaceVerts);
                }
            }

            // Create cap to close mesh
            if (capVertices.Count >= 3)
            {
                CreateCapFromVertices(result, capVertices, cutter, true);
            }

            // Final cleanup
            result.Vertices.CombineIdentical(true, true);
            result.Faces.CullDegenerateFaces();
            result.Compact();

            // Try to close any remaining holes
            if (!result.IsClosed)
            {
                result = TryCloseMesh(result, cutter);
            }

            result.Normals.ComputeNormals();
            result.UnifyNormals();

            // Apply colors
            if (hasColors)
            {
                result.VertexColors.Clear();
                result.VertexColors.CreateMonotoneMesh(mainColor);
            }

            if (result.Faces.Count > 0)
                meshes[index] = result;
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

        /// <summary>
        /// Creates faces from polygon using only corner vertices
        /// Prefers quads for structural analysis
        /// </summary>
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
                    int[] split = GetBestQuadDiagonal(pts);
                    mesh.Faces.AddFace(vertIndices[split[0]], vertIndices[split[1]], vertIndices[split[2]]);
                    mesh.Faces.AddFace(vertIndices[split[0]], vertIndices[split[2]], vertIndices[split[3]]);
                }
                return;
            }

            // For n > 4: Fan triangulation from best pivot
            Point3d[] points = vertIndices.Select(i => new Point3d(mesh.Vertices[i])).ToArray();
            Point3d centroid = new Point3d(
                points.Average(p => p.X),
                points.Average(p => p.Y),
                points.Average(p => p.Z));

            // Find vertex closest to centroid for pivot
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

            // Reorder so pivot is first
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

            // Check coplanarity
            if (Plane.FitPlaneToPoints(pts.ToList(), out Plane plane) != PlaneFitResult.Success)
                return false;

            const double planeTol = 0.01;
            foreach (var pt in pts)
                if (Math.Abs(plane.DistanceTo(pt)) > planeTol)
                    return false;

            // Project to 2D
            Transform toPlane = Transform.PlaneToPlane(plane, Plane.WorldXY);
            Point2d[] pts2d = new Point2d[4];
            for (int i = 0; i < 4; i++)
            {
                Point3d p = new Point3d(pts[i]);
                p.Transform(toPlane);
                pts2d[i] = new Point2d(p.X, p.Y);
            }

            // Check convexity
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
            // Compare aspect ratios for both diagonal options
            double aspect_02 = Math.Min(
                GetTriangleAspectRatio(pts[0], pts[1], pts[2]),
                GetTriangleAspectRatio(pts[0], pts[2], pts[3]));

            double aspect_13 = Math.Min(
                GetTriangleAspectRatio(pts[1], pts[2], pts[3]),
                GetTriangleAspectRatio(pts[1], pts[3], pts[0]));

            if (aspect_02 >= aspect_13)
                return new[] { 0, 1, 2, 3 };
            else
                return new[] { 1, 2, 3, 0 };
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

        /// <summary>
        /// Creates cap face from boundary vertices
        /// </summary>
        private void CreateCapFromVertices(Mesh mesh, List<int> capVerts, Plane cutter, bool flipNormal)
        {
            if (capVerts.Count < 3) return;

            // Get unique vertices on cutting plane
            List<int> uniqueVerts = capVerts.Distinct().ToList();

            // Filter to only vertices actually on the plane
            List<int> onPlaneVerts = new List<int>();
            foreach (int idx in uniqueVerts)
            {
                if (Math.Abs(cutter.DistanceTo(mesh.Vertices[idx])) < 0.01)
                    onPlaneVerts.Add(idx);
            }

            if (onPlaneVerts.Count < 3) return;

            // Get points and compute centroid
            List<Point3d> points = onPlaneVerts.Select(i => new Point3d(mesh.Vertices[i])).ToList();
            Point3d centroid = new Point3d(
                points.Average(p => p.X),
                points.Average(p => p.Y),
                points.Average(p => p.Z));

            // Sort by angle around centroid
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

            // Create cap with correct orientation
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
                // Fan triangulation for n > 4
                for (int i = 1; i < n - 1; i++)
                {
                    if (flipNormal)
                        mesh.Faces.AddFace(orderedVerts[0], orderedVerts[i + 1], orderedVerts[i]);
                    else
                        mesh.Faces.AddFace(orderedVerts[0], orderedVerts[i], orderedVerts[i + 1]);
                }
            }
        }

        /// <summary>
        /// Attempts to close remaining holes
        /// </summary>
        private Mesh TryCloseMesh(Mesh mesh, Plane cutter)
        {
            var nakedEdges = mesh.GetNakedEdges();
            if (nakedEdges == null || nakedEdges.Length == 0)
                return mesh;

            foreach (var polyline in nakedEdges)
            {
                if (polyline == null || polyline.Count < 3)
                    continue;

                // Find vertex indices for naked edge
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
                        CreateCapFromVertices(mesh, boundaryVerts, cutter, true);
                    }
                    else
                    {
                        // Fit a plane to these vertices
                        List<Point3d> pts = boundaryVerts.Select(i => new Point3d(mesh.Vertices[i])).ToList();
                        if (Plane.FitPlaneToPoints(pts, out Plane holePlane) == PlaneFitResult.Success)
                        {
                            CreateCapFromVertices(mesh, boundaryVerts, holePlane, false);
                        }
                    }
                }
            }

            mesh.Vertices.CombineIdentical(true, true);
            mesh.Normals.ComputeNormals();
            mesh.UnifyNormals();

            return mesh;
        }

        // ==================== UTILITY METHODS ====================

        private bool MeshPlaneIntersection(Mesh mesh, Plane plane)
        {
            var bbox = mesh.GetBoundingBox(false);
            Point3d[] corners = bbox.GetCorners();
            bool pos = false, neg = false;
            foreach (var c in corners)
            {
                double d = plane.DistanceTo(c);
                if (d > 0.001) pos = true;
                if (d < -0.001) neg = true;
                if (pos && neg) return true;
            }
            return false;
        }

        private Color GetAverageColor(Color[] colors)
        {
            int r = 0, g = 0, b = 0, count = 0;
            foreach (var c in colors)
            {
                if (c.A > 0)
                {
                    r += c.R;
                    g += c.G;
                    b += c.B;
                    count++;
                }
            }
            return count > 0 ? Color.FromArgb(r / count, g / count, b / count) : Color.White;
        }

        public override BoundingBox ClippingBox
        {
            get
            {
                var bbox = BoundingBox.Empty;
                foreach (var c in debugCurves) bbox.Union(c.GetBoundingBox(true));
                return bbox;
            }
        }

        public override void DrawViewportWires(IGH_PreviewArgs args)
        {
            base.DrawViewportWires(args);
        }

        public static Bitmap ConvertByteArrayToBitmap(byte[] imageData)
        {
            using (MemoryStream ms = new MemoryStream(imageData)) return new Bitmap(ms);
        }

        protected override Bitmap Icon => ConvertByteArrayToBitmap(PoliBrick.Properties.Resources.CutterIcon);
        public override Guid ComponentGuid => new Guid("5A3F6A0E-CB24-4E98-99F3-61833CA0F8A8");
    }
}