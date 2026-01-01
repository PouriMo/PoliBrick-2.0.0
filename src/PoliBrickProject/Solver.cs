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
using System.Drawing.Drawing2D;
using System.IO;
using System.Linq;
using System.Threading.Tasks;
using System.Windows.Forms;

namespace PoliBrick.Util
{
    public class Solver : GH_Component
    {
        public bool FlipNormals { get; set; } = false;

        public Solver() : base("1.Solver", "solver", "Generate bricks using a Data Library", "PoliBrick", "Utilities") { }
        public override GH_Exposure Exposure => GH_Exposure.secondary;
        public override void CreateAttributes() { m_attributes = new SolverAttributes(this); }

        public override bool Write(GH_IO.Serialization.GH_IWriter writer)
        {
            writer.SetBoolean("FlipNormals", this.FlipNormals);
            return base.Write(writer);
        }

        public override bool Read(GH_IO.Serialization.GH_IReader reader)
        {
            if (reader.ItemExists("FlipNormals"))
                this.FlipNormals = reader.GetBoolean("FlipNormals");
            return base.Read(reader);
        }

        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddCurveParameter("Curves", "C", "Input curves (List ordered from bottom to top)", GH_ParamAccess.tree);
            pManager.AddTextParameter("Pattern", "P", "Brick Type Pattern string (Tree/List). Cycles through rows.", GH_ParamAccess.tree);
            pManager.AddGenericParameter("Data", "D", "Data Library (Tree) from BrData", GH_ParamAccess.tree);
            pManager.AddSurfaceParameter("Surface", "S", "Target surface", GH_ParamAccess.item);
        }

        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddGenericParameter("Bricks (Ordered)", "B", "Bricks ordered by curve sequence", GH_ParamAccess.tree);
            pManager.AddGenericParameter("Bricks (Type)", "BT", "Bricks separated by Type", GH_ParamAccess.tree);
            int planeIdx = pManager.AddGenericParameter("Face Planes", "FP", "Ordered planes per block", GH_ParamAccess.tree);
            if (pManager[planeIdx] is IGH_PreviewObject previewParam) previewParam.Hidden = true;
        }

        class CurveJob
        {
            public Curve Curve;
            public GH_Path OriginalPath;
            public string SpecificPattern; // The Brick Type Pattern (e.g. "ABAB")
            public int JobIndex;
            public Vector3d GuideVector;
        }

        struct BrickResult
        {
            public Mesh Mesh;
            public int TypeIndex;
            public List<Plane> Planes;
        }

        protected override void SolveInstance(IGH_DataAccess DA)
        {
            // 1. INPUTS
            if (!DA.GetDataTree(0, out GH_Structure<GH_Curve> curvesTree)) return;
            if (!DA.GetDataTree(1, out GH_Structure<GH_String> patternTree)) return;
            if (!DA.GetDataTree(2, out GH_Structure<IGH_Goo> rawDataTree)) return;

            Surface surface = null;
            if (!DA.GetData(3, ref surface)) return;

            // 2. BUILD DATA LIBRARY
            Dictionary<string, List<string>> dataLibrary = new Dictionary<string, List<string>>();
            foreach (GH_Path path in rawDataTree.Paths)
            {
                var branch = rawDataTree.get_Branch(path);
                List<string> strList = new List<string>();
                foreach (var item in branch)
                {
                    if (item is GH_String s) strList.Add(s.Value);
                    else if (item != null) strList.Add(item.ToString());
                }

                if (strList.Count > 0 && !string.IsNullOrEmpty(strList[0]))
                {
                    string key = strList[0];
                    if (!dataLibrary.ContainsKey(key)) dataLibrary.Add(key, strList);
                }
            }

            if (dataLibrary.Count == 0) return;

            // 3. SORT KEYS
            List<string> sortedKeys = dataLibrary.Keys.OrderBy(k => k).ToList();
            Dictionary<string, int> keyMap = new Dictionary<string, int>();
            for (int i = 0; i < sortedKeys.Count; i++) keyMap[sortedKeys[i]] = i;

            // 4. PREPARE CURVES & VECTORS
            List<GH_Curve> flatCurves = new List<GH_Curve>();
            foreach (var path in curvesTree.Paths)
                flatCurves.AddRange(curvesTree.get_Branch(path).OfType<GH_Curve>());

            if (flatCurves.Count == 0) return;

            List<Vector3d> guideVectors = new List<Vector3d>();
            for (int i = 0; i < flatCurves.Count; i++)
            {
                Curve current = flatCurves[i].Value;
                if (current == null) { guideVectors.Add(Vector3d.ZAxis); continue; }

                Point3d midCurr = current.PointAtNormalizedLength(0.5);
                Vector3d up;

                if (flatCurves.Count > 1)
                {
                    if (i < flatCurves.Count - 1)
                    {
                        Curve next = flatCurves[i + 1].Value;
                        Point3d midNext = next.PointAtNormalizedLength(0.5);
                        up = midNext - midCurr;
                    }
                    else
                    {
                        Curve prev = flatCurves[i - 1].Value;
                        Point3d midPrev = prev.PointAtNormalizedLength(0.5);
                        up = midCurr - midPrev;
                    }
                }
                else
                {
                    surface.ClosestPoint(midCurr, out double u, out double v);
                    Vector3d norm = surface.NormalAt(u, v);
                    Vector3d tan = current.TangentAt(current.Domain.Mid);
                    up = Vector3d.CrossProduct(tan, norm);
                }

                if (up.IsTiny(1e-6)) up = Vector3d.ZAxis;
                guideVectors.Add(up);
            }

            // 5. CREATE JOBS
            List<CurveJob> allJobs = new List<CurveJob>();

            // Flatten Pattern Tree to allow cycling through a simple list (Row 1, Row 2...)
            List<string> flatPatterns = new List<string>();
            foreach (var item in patternTree.AllData(true))
            {
                if (item is GH_String s) flatPatterns.Add(s.Value);
            }

            if (flatPatterns.Count == 0) return;

            int branchCounter = 0;
            int globalCurveIndex = 0;

            foreach (GH_Path path in curvesTree.Paths)
            {
                var curveBranch = curvesTree.get_Branch(path);

                // Cycle through flattened patterns using Branch Index
                string assignedPattern = flatPatterns[branchCounter % flatPatterns.Count];

                foreach (var item in curveBranch)
                {
                    if (item is GH_Curve ghCrv && ghCrv.Value != null)
                    {
                        allJobs.Add(new CurveJob
                        {
                            Curve = ghCrv.Value,
                            OriginalPath = path,
                            SpecificPattern = assignedPattern,
                            JobIndex = globalCurveIndex,
                            GuideVector = guideVectors[globalCurveIndex]
                        });
                        globalCurveIndex++;
                    }
                }
                branchCounter++;
            }

            // 6. GENERATED FALLBACK COLORS
            Dictionary<char, Color> fallbackColors = new Dictionary<char, Color>();
            Random rand = new Random(45);
            HashSet<int> usedHues = new HashSet<int>();

            foreach (string key in sortedKeys)
            {
                if (key.Length == 0) continue;
                char c = key[0];
                if (c == 'S') { fallbackColors[c] = Color.FromArgb(0, 0, 0, 0); continue; }
                int hue; int attempts = 0;
                do { hue = rand.Next(0, 360); attempts++; } while (usedHues.Any(h => Math.Abs(h - hue) < 30) && attempts < 100);
                usedHues.Add(hue); fallbackColors[c] = FromHsv(hue, 1.0, 1.0);
            }

            // 7. EXECUTION
            ConcurrentDictionary<int, List<BrickResult>> results = new ConcurrentDictionary<int, List<BrickResult>>();

            Parallel.For(0, allJobs.Count, j =>
            {
                CurveJob job = allJobs[j];
                Curve curve = job.Curve;
                Vector3d guideVec = job.GuideVector;

                if (string.IsNullOrEmpty(job.SpecificPattern)) return;

                bool isClosed = curve.IsClosed;
                double curveTotalLength = curve.GetLength();
                double accumulatedLength = 0.0;

                Dictionary<string, int> brickCounters = new Dictionary<string, int>();

                Point3d start = curve.PointAtStart;
                Point3d originalStart = start;

                List<BrickResult> localBricks = new List<BrickResult>();

                string basePattern = job.SpecificPattern;
                bool repeatPatternDynamically = false;
                if (basePattern.StartsWith("__") && basePattern.EndsWith("__"))
                {
                    basePattern = basePattern.Substring(2, basePattern.Length - 4);
                    repeatPatternDynamically = true;
                }
                string dynamicPattern = basePattern;

                int k = 0; int maxIterations = 10000;

                while (k < maxIterations)
                {
                    if (k >= dynamicPattern.Length) { if (repeatPatternDynamically) dynamicPattern += basePattern; else break; }

                    char patternChar = dynamicPattern[k];
                    string charStr = patternChar.ToString();

                    int currentTypeIndex = -1;
                    if (keyMap.ContainsKey(charStr)) currentTypeIndex = keyMap[charStr];

                    if (dataLibrary.TryGetValue(charStr, out List<string> dataBranch))
                    {
                        if (!brickCounters.ContainsKey(charStr)) brickCounters[charStr] = -1;
                        brickCounters[charStr]++;
                        int currentCount = brickCounters[charStr];

                        if (dataBranch.Count >= 2)
                        {
                            var dimData = dataBranch[1].Split(',');
                            if (dimData.Length >= 4)
                            {
                                double u = Convert.ToDouble(dimData[0]);
                                double v = Convert.ToDouble(dimData[1]);
                                double w = Convert.ToDouble(dimData[2]);
                                double m = Convert.ToDouble(dimData[3]);

                                // LOGIC: Rotation Pattern comes from Data Library (Index 3)
                                string actualRotationPattern = (dataBranch.Count > 3) ? dataBranch[3] : "";

                                // LOGIC: Color Selection (Index 4)
                                Color finalColor = Color.Empty;
                                if (dataBranch.Count > 4)
                                {
                                    string userColorStr = dataBranch[4];
                                    if (!string.IsNullOrEmpty(userColorStr) && userColorStr != "Empty")
                                    {
                                        var argb = userColorStr.Split(',');
                                        if (argb.Length == 4)
                                        {
                                            finalColor = Color.FromArgb(
                                                int.Parse(argb[0]), int.Parse(argb[1]), int.Parse(argb[2]), int.Parse(argb[3]));
                                        }
                                    }
                                }
                                // Fallback
                                if (finalColor.IsEmpty)
                                {
                                    if (fallbackColors.TryGetValue(patternChar, out Color fb)) finalColor = fb;
                                    else finalColor = Color.Gray;
                                }

                                // Check Closure
                                if (isClosed && k > 0 && accumulatedLength > curveTotalLength * 0.75)
                                {
                                    double distToOrigin = start.DistanceTo(originalStart);
                                    if (distToOrigin < u + m)
                                    {
                                        if (patternChar != 'S' && distToOrigin > 1e-6)
                                        {
                                            PlaceFinalBrick(start, originalStart, v, w, surface, dataBranch, actualRotationPattern, patternChar, currentTypeIndex, currentCount, ref localBricks, finalColor, FlipNormals, guideVec);
                                        }
                                        goto EndOfCurveProcessing;
                                    }
                                }

                                Point3d end = FindPointAtChordLength(curve, start, u);
                                if (!end.IsValid)
                                {
                                    if (!isClosed)
                                    {
                                        end = curve.PointAtEnd;
                                        if (patternChar != 'S' && start.DistanceTo(end) > 1e-6)
                                        {
                                            PlaceFinalBrick(start, end, v, w, surface, dataBranch, actualRotationPattern, patternChar, currentTypeIndex, currentCount, ref localBricks, finalColor, FlipNormals, guideVec);
                                        }
                                    }
                                    goto EndOfCurveProcessing;
                                }

                                if (patternChar != 'S')
                                {
                                    PlaceFinalBrick(start, end, v, w, surface, dataBranch, actualRotationPattern, patternChar, currentTypeIndex, currentCount, ref localBricks, finalColor, FlipNormals, guideVec);
                                }

                                accumulatedLength += u + m;
                                start = end;
                                if (m > 0)
                                {
                                    Point3d nextStartAfterGap = FindPointAtChordLength(curve, start, m);
                                    if (nextStartAfterGap.IsValid) start = nextStartAfterGap;
                                    else goto EndOfCurveProcessing;
                                }
                            }
                        }
                    }
                    k++;
                }
            EndOfCurveProcessing:;
                results[j] = localBricks;
            });

            // 8. OUTPUT
            DataTree<GH_Mesh> bricksOrderedTree = new DataTree<GH_Mesh>();
            DataTree<GH_Mesh> bricksTypeTree = new DataTree<GH_Mesh>();
            DataTree<GH_Plane> planesTree = new DataTree<GH_Plane>();

            foreach (var kvp in results.OrderBy(x => x.Key))
            {
                int jobIdx = kvp.Key;
                GH_Path basePath = new GH_Path(jobIdx);
                int blockCounter = 0;

                foreach (var item in kvp.Value)
                {
                    Mesh mesh = item.Mesh;
                    List<Plane> faces = item.Planes;
                    int typeIdx = item.TypeIndex;

                    if (mesh != null)
                    {
                        bricksOrderedTree.Add(new GH_Mesh(mesh), basePath);

                        GH_Path typePath = new GH_Path(basePath);
                        typePath = typePath.AppendElement(typeIdx);
                        bricksTypeTree.Add(new GH_Mesh(mesh), typePath);

                        if (faces != null)
                        {
                            GH_Path blockPath = new GH_Path(jobIdx, blockCounter);
                            foreach (Plane p in faces) planesTree.Add(new GH_Plane(p), blockPath);
                        }
                        blockCounter++;
                    }
                }
            }

            DA.SetDataTree(0, bricksOrderedTree);
            DA.SetDataTree(1, bricksTypeTree);
            DA.SetDataTree(2, planesTree);
        }

        private void PlaceFinalBrick(Point3d start, Point3d end, double v, double w, Surface surface, List<string> dataBranch, string rotationPattern, char patternChar, int typeIndex, int counter, ref List<BrickResult> localBricks, Color brickColor, bool flip, Vector3d guideVec)
        {
            double actualChord = start.DistanceTo(end);
            Mesh meshBlock = CreateMesh(actualChord, v, w);
            meshBlock.VertexColors.CreateMonotoneMesh(brickColor);

            ProcessMesh(meshBlock, start, end, actualChord, v, w, surface, dataBranch, rotationPattern, patternChar, typeIndex, counter, ref localBricks, flip, guideVec);
        }

        private void ProcessMesh(Mesh mesh, Point3d start, Point3d end, double u, double v, double w, Surface surface, List<string> dataBranch, string rotationPattern, char patternChar, int typeIndex, int counter, ref List<BrickResult> results, bool manualFlip, Vector3d guideVector)
        {
            Point3d mid = 0.5 * start + 0.5 * end;
            surface.ClosestPoint(mid, out double surU, out double surV);
            Vector3d zAxis = surface.NormalAt(surU, surV);
            if (manualFlip) { zAxis.Reverse(); }

            Vector3d yAxis = end - start;
            Vector3d xAxis = Vector3d.CrossProduct(yAxis, zAxis);
            yAxis = Vector3d.CrossProduct(zAxis, xAxis);

            Plane newPlane = new Plane(mid, xAxis, yAxis);

            bool swapTopBottom = (Vector3d.Multiply(xAxis, guideVector) < 0);

            double zMid = w / 2.0;
            List<Plane> facePlanes = new List<Plane>();
            if (patternChar != 'S')
            {
                facePlanes.Add(new Plane(new Point3d(0, u / 2.0, zMid), new Vector3d(0, 1, 0)));
                facePlanes.Add(new Plane(new Point3d(0, -u / 2.0, zMid), new Vector3d(0, -1, 0)));

                if (!swapTopBottom)
                {
                    facePlanes.Add(new Plane(new Point3d(v / 2.0, 0, zMid), new Vector3d(1, 0, 0)));
                    facePlanes.Add(new Plane(new Point3d(-v / 2.0, 0, zMid), new Vector3d(-1, 0, 0)));
                }
                else
                {
                    facePlanes.Add(new Plane(new Point3d(-v / 2.0, 0, zMid), new Vector3d(-1, 0, 0)));
                    facePlanes.Add(new Plane(new Point3d(v / 2.0, 0, zMid), new Vector3d(1, 0, 0)));
                }
                facePlanes.Add(new Plane(new Point3d(0, 0, 0), new Vector3d(0, 0, -1)));
                facePlanes.Add(new Plane(new Point3d(0, 0, w), new Vector3d(0, 0, 1)));
            }

            Transform finalTransform = Transform.Identity;

            if (dataBranch.Count < 3)
            {
                finalTransform = Transform.PlaneToPlane(Plane.WorldXY, newPlane);
            }
            else
            {
                var move = dataBranch[2].Split(',');
                double alfa = Convert.ToDouble(move[0]) * Math.PI / 180;
                double beta = Convert.ToDouble(move[1]) * Math.PI / 180;
                double theta = Convert.ToDouble(move[2]) * Math.PI / 180;

                // 3-Axis Rotation applied only if pattern char is 'T'
                if (rotationPattern.Length > 0 && rotationPattern[counter % rotationPattern.Length] == 'T')
                {
                    newPlane.Rotate(alfa, newPlane.YAxis);
                    newPlane.Rotate(beta, newPlane.ZAxis);
                    newPlane.Rotate(theta, newPlane.XAxis);
                }

                finalTransform = Transform.PlaneToPlane(Plane.WorldXY, newPlane);
            }

            mesh.Transform(finalTransform);
            for (int i = 0; i < facePlanes.Count; i++)
            {
                Plane p = facePlanes[i];
                p.Transform(finalTransform);
                facePlanes[i] = p;
            }

            if (patternChar == 'S') results.Add(new BrickResult { Mesh = new Mesh(), TypeIndex = typeIndex, Planes = new List<Plane>() });
            else results.Add(new BrickResult { Mesh = mesh, TypeIndex = typeIndex, Planes = facePlanes });
        }

        private Point3d FindPointAtChordLength(Curve curve, Point3d startPoint, double chordLength)
        {
            if (curve == null || chordLength <= 1e-6) return Point3d.Unset;
            double startParam; curve.ClosestPoint(startPoint, out startParam);
            double chordLengthSquared = chordLength * chordLength;
            if (startPoint.DistanceToSquared(curve.PointAtEnd) < chordLengthSquared - 1e-9 && !curve.IsClosed) { return Point3d.Unset; }
            double minParam = startParam;
            double maxParam = curve.Domain.Max;
            if (curve.IsClosed) { maxParam = startParam + curve.Domain.Length; }
            for (int i = 0; i < 20; i++)
            {
                double midParam = (minParam + maxParam) / 2.0;
                double distSquared = startPoint.DistanceToSquared(curve.PointAt(midParam));
                if (distSquared < chordLengthSquared) minParam = midParam; else maxParam = midParam;
            }
            double finalParam = (minParam + maxParam) / 2.0;
            return curve.PointAt(finalParam);
        }

        private Mesh CreateMesh(double u, double v, double w)
        {
            Mesh mesh = new Mesh();
            mesh.Vertices.AddVertices(new List<Point3d> { new Point3d(-v / 2, -u / 2, 0), new Point3d(v / 2, -u / 2, 0), new Point3d(v / 2, u / 2, 0), new Point3d(-v / 2, u / 2, 0), new Point3d(-v / 2, -u / 2, w), new Point3d(v / 2, -u / 2, w), new Point3d(v / 2, u / 2, w), new Point3d(-v / 2, u / 2, w) });
            mesh.Faces.AddFace(0, 1, 2, 3); mesh.Faces.AddFace(4, 5, 6, 7); mesh.Faces.AddFace(0, 3, 7, 4); mesh.Faces.AddFace(1, 2, 6, 5); mesh.Faces.AddFace(0, 1, 5, 4); mesh.Faces.AddFace(2, 3, 7, 6);
            mesh.Vertices.CombineIdentical(true, true); mesh.Normals.ComputeNormals();
            mesh.UnifyNormals(); return mesh;
        }

        private Color FromHsv(double h, double s, double v)
        {
            double c = v * s; double x = c * (1 - Math.Abs((h / 60) % 2 - 1)); double m = v - c; double r = 0, g = 0, b = 0;
            if (h < 60) { r = c; g = x; } else if (h < 120) { r = x; g = c; } else if (h < 180) { g = c; b = x; } else if (h < 240) { g = x; b = c; } else if (h < 300) { r = x; b = c; } else { r = c; b = x; }
            return Color.FromArgb(255, (int)((r + m) * 255), (int)((g + m) * 255), (int)((b + m) * 255));
        }

        public static Bitmap ConvertByteArrayToBitmap(byte[] imageData) { using (MemoryStream ms = new MemoryStream(imageData)) { return new Bitmap(ms); } }
        protected override Bitmap Icon => ConvertByteArrayToBitmap(PoliBrick.Properties.Resources.SolverIcon);
        public override Guid ComponentGuid => new Guid("9C470AF1-CEE1-416A-A564-C9E91363302E");
    }

    public class SolverAttributes : GH_ComponentAttributes
    {
        private Solver _owner;
        public SolverAttributes(Solver owner) : base(owner) { _owner = owner; }
        protected override void Layout() { base.Layout(); Bounds = new RectangleF(Bounds.X, Bounds.Y, Bounds.Width, Bounds.Height + 25); }
        protected override void Render(GH_Canvas canvas, Graphics graphics, GH_CanvasChannel channel)
        {
            base.Render(canvas, graphics, channel);
            if (channel == GH_CanvasChannel.Objects)
            {
                float buttonHeight = 18; int gap = 3;
                Rectangle r = new Rectangle((int)Bounds.X + gap, (int)Bounds.Bottom - (int)buttonHeight - gap, (int)Bounds.Width - (gap * 2), (int)buttonHeight);
                bool active = _owner.FlipNormals;
                Color colorOff = Color.FromArgb(225, 225, 225); Color colorOn = Color.Orange; Color colorBorder = Color.Gray; Color textOff = Color.DimGray; Color textOn = Color.White;
                using (GraphicsPath path = RoundedRect(r, r.Height / 4)) { using (Brush b = new SolidBrush(active ? colorOn : colorOff)) { graphics.FillPath(b, path); } using (Pen p = new Pen(colorBorder)) { graphics.DrawPath(p, path); } }
                using (Font font = new Font("Arial", 7f, FontStyle.Bold)) using (Brush textBrush = new SolidBrush(active ? textOn : textOff)) { StringFormat fmt = new StringFormat { Alignment = StringAlignment.Center, LineAlignment = StringAlignment.Center }; graphics.DrawString("Flip", font, textBrush, r, fmt); }
            }
        }
        public override GH_ObjectResponse RespondToMouseDown(GH_Canvas sender, GH_CanvasMouseEvent e)
        {
            if (e.Button == MouseButtons.Left)
            {
                int buttonHeight = 18; int gap = 3;
                RectangleF buttonRect = new RectangleF(Bounds.X + gap, Bounds.Bottom - buttonHeight - gap, Bounds.Width - (gap * 2), buttonHeight);
                if (buttonRect.Contains(e.CanvasLocation)) { _owner.FlipNormals = !_owner.FlipNormals; _owner.ExpireSolution(true); return GH_ObjectResponse.Handled; }
            }
            return base.RespondToMouseDown(sender, e);
        }
        private GraphicsPath RoundedRect(Rectangle r, int rad) { GraphicsPath p = new GraphicsPath(); int d = rad * 2; p.AddArc(r.X, r.Y, d, d, 180, 90); p.AddArc(r.Right - d, r.Y, d, d, 270, 90); p.AddArc(r.Right - d, r.Bottom - d, d, d, 0, 90); p.AddArc(r.X, r.Bottom - d, d, d, 90, 90); p.CloseFigure(); return p; }
    }
}