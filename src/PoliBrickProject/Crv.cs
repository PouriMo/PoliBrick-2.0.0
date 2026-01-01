using GH_IO.Serialization;
using Grasshopper;
using Grasshopper.GUI;
using Grasshopper.GUI.Canvas;
using Grasshopper.Kernel;
using Grasshopper.Kernel.Attributes;
using Grasshopper.Kernel.Data;
using Grasshopper.Kernel.Types;
using Rhino.Geometry;
using System;
using System.Collections.Generic;
using System.Drawing;
using System.Drawing.Drawing2D;
using System.IO;
using System.Linq;
using System.Threading.Tasks;
using System.Windows.Forms;

public class Crv : GH_Component
{
    public bool RunClicked { get; set; } = false;
    public bool ResetClicked { get; set; } = false;
    private bool _runPressed = false;
    private bool _resetPressed = false;
    private bool _isProcessing = false;

    // Toggles
    public bool ExtendCurves { get; set; } = true;
    public bool FilterIntersections { get; set; } = true;

    // State
    private NurbsCurve initialCurve;
    private NurbsCurve fixedInitialCurve;
    public IterationState iterationState = null;

    public Crv() : base("2.ParallelCurve", "PC", "Generate parallel curves on mesh", "PoliBrick", "Utilities") { }
    public override GH_Exposure Exposure
    {
        get { return GH_Exposure.primary; }
    }
    public override Guid ComponentGuid => new Guid("d6fee1e8-fde9-46c9-8667-589a0a7b8b47");

    // ==========================================================
    // 1. PERSISTENCE (Saving Curves & Settings)
    // ==========================================================
    public override bool Write(GH_IWriter writer)
    {
        // 1. Save Settings
        writer.SetBoolean("ExtendCurves", ExtendCurves);
        writer.SetBoolean("FilterIntersections", FilterIntersections);

        // 2. Save Curve Data
        bool hasData = (iterationState != null &&
                       (iterationState.CurvesDir1.Count > 0 || iterationState.CurvesDir2.Count > 0));

        writer.SetBoolean("HasData", hasData);

        if (hasData)
        {
            // Serialize Dir1
            writer.SetInt32("Dir1_BranchCount", iterationState.CurvesDir1.Count);
            for (int i = 0; i < iterationState.CurvesDir1.Count; i++)
            {
                var branch = iterationState.CurvesDir1[i];
                writer.SetInt32($"D1_B{i}_Count", branch.Count);
                for (int j = 0; j < branch.Count; j++)
                {
                    byte[] bytes = GH_Convert.CommonObjectToByteArray(branch[j]);
                    writer.SetByteArray($"D1_B{i}_C{j}", bytes);
                }
            }

            // Serialize Dir2
            writer.SetInt32("Dir2_BranchCount", iterationState.CurvesDir2.Count);
            for (int i = 0; i < iterationState.CurvesDir2.Count; i++)
            {
                var branch = iterationState.CurvesDir2[i];
                writer.SetInt32($"D2_B{i}_Count", branch.Count);
                for (int j = 0; j < branch.Count; j++)
                {
                    byte[] bytes = GH_Convert.CommonObjectToByteArray(branch[j]);
                    writer.SetByteArray($"D2_B{i}_C{j}", bytes);
                }
            }

            // Save state meta-data
            writer.SetInt32("State_MaxIter", iterationState.MaxIterations);
            writer.SetInt32("State_Mode", iterationState.Mode);
            writer.SetInt32("State_K", iterationState.CurrentK);
            writer.SetBoolean("State_Complete", iterationState.IsComplete);
        }

        return base.Write(writer);
    }

    public override bool Read(GH_IReader reader)
    {
        // 1. Load Settings
        if (reader.ItemExists("ExtendCurves")) ExtendCurves = reader.GetBoolean("ExtendCurves");
        if (reader.ItemExists("FilterIntersections")) FilterIntersections = reader.GetBoolean("FilterIntersections");

        // 2. Load Curve Data
        if (reader.ItemExists("HasData") && reader.GetBoolean("HasData"))
        {
            int maxIter = reader.GetInt32("State_MaxIter");
            int mode = reader.GetInt32("State_Mode");

            iterationState = new IterationState(maxIter, mode);
            iterationState.CurrentK = reader.GetInt32("State_K");
            iterationState.IsComplete = reader.GetBoolean("State_Complete");

            // Load Dir1
            int d1Count = reader.GetInt32("Dir1_BranchCount");
            for (int i = 0; i < d1Count; i++)
            {
                int bCount = reader.GetInt32($"D1_B{i}_Count");
                List<Curve> branch = new List<Curve>();
                for (int j = 0; j < bCount; j++)
                {
                    byte[] bytes = reader.GetByteArray($"D1_B{i}_C{j}");
                    var crv = GH_Convert.ByteArrayToCommonObject<Curve>(bytes);
                    if (crv != null) branch.Add(crv);
                }
                iterationState.CurvesDir1.Add(branch);
            }

            // Load Dir2
            int d2Count = reader.GetInt32("Dir2_BranchCount");
            for (int i = 0; i < d2Count; i++)
            {
                int bCount = reader.GetInt32($"D2_B{i}_Count");
                List<Curve> branch = new List<Curve>();
                for (int j = 0; j < bCount; j++)
                {
                    byte[] bytes = reader.GetByteArray($"D2_B{i}_C{j}");
                    var crv = GH_Convert.ByteArrayToCommonObject<Curve>(bytes);
                    if (crv != null) branch.Add(crv);
                }
                iterationState.CurvesDir2.Add(branch);
            }
        }

        return base.Read(reader);
    }
    // ==========================================================

    public static Bitmap ConvertByteArrayToBitmap(byte[] imageData)
    {
        using (MemoryStream ms = new MemoryStream(imageData)) { return new Bitmap(ms); }
    }
    protected override Bitmap Icon => ConvertByteArrayToBitmap(PoliBrick.Properties.Resources.ParaIcon);

    protected override void RegisterInputParams(GH_InputParamManager pManager)
    {
        pManager.AddSurfaceParameter("Surface", "srf", "The base surface", GH_ParamAccess.item);
        pManager.AddCurveParameter("Curve", "crv", "The initial curve", GH_ParamAccess.item);
        pManager.AddNumberParameter("Distance", "dist", "The offset distance", GH_ParamAccess.item);
        pManager.AddIntegerParameter("Count", "c", "Number of parallel curves", GH_ParamAccess.item);
        pManager.AddNumberParameter("Accuracy", "acc", "Accuracy (0-1)", GH_ParamAccess.item, 0.5);
        pManager.AddIntegerParameter("Mode", "mode", "0=normal, 1=reverse, 2=dual-limited, 3=dual-unlimited", GH_ParamAccess.item, 0);
    }

    protected override void RegisterOutputParams(GH_OutputParamManager pManager)
    {
        pManager.AddCurveParameter("Curves", "crvs", "Generated curves", GH_ParamAccess.tree);
        pManager.AddCurveParameter("Ordered Curves", "OrderedCrvs", "Curves re-ordered spatially", GH_ParamAccess.tree);
        pManager.AddCurveParameter("Refined Curves", "RefCrv","Refined curve rows for meshing. Paths: {n}=base row, {n;0}=stitched/refined, {n;1}=unjoined boundary pieces.",
        GH_ParamAccess.tree);
    }

    protected override void SolveInstance(IGH_DataAccess DA)
    {
        Surface surface = null;
        Curve curve = null;
        double distance = 0;
        int count = 0;
        double accuracy = 0.5;
        int flip = 0;

        if (!DA.GetData(0, ref surface)) return;
        if (!DA.GetData(1, ref curve)) return;
        if (!DA.GetData(2, ref distance)) return;
        if (!DA.GetData(3, ref count)) return;
        DA.GetData(4, ref accuracy);
        DA.GetData(5, ref flip);

        // Reset Logic
        if (_resetPressed)
        {
            _isProcessing = false;
            iterationState = new IterationState((flip == 3) ? 0 : count, flip);
            _resetPressed = false;
            DA.SetDataTree(0, new GH_Structure<GH_Curve>());
            DA.SetDataTree(1, new GH_Structure<GH_Curve>());
            DA.SetDataTree(2, new GH_Structure<GH_Curve>());
            return;
        }

        // Initialize State if null (First Run)
        if (iterationState == null)
        {
            iterationState = new IterationState((flip == 3) ? 0 : count, flip);
        }
        
        if (!_isProcessing && !_runPressed && iterationState.CurvesDir1.Count == 0 && iterationState.CurvesDir2.Count == 0)
        {
            // Only update config if empty, otherwise we preserve loaded data
            if (iterationState.Mode != flip) iterationState = new IterationState((flip == 3) ? 0 : count, flip);
        }

        // Run Logic
        if (_runPressed)
        {
            // If parameters changed, re-init state
            if (iterationState.Mode != flip || (iterationState.MaxIterations != count && flip != 3))
            {
                iterationState = new IterationState((flip == 3) ? 0 : count, flip);
            }
            StartAutomaticProcessing(surface, curve, accuracy, distance, flip);
            _runPressed = false;
        }

        // --- OUTPUT GENERATION (Run from State) ---
        // This runs both during active processing AND after loading from file
        var tree = new GH_Structure<GH_Curve>();
        var orderedTree = new GH_Structure<GH_Curve>();
        int branchCounter = 0;
        List<List<Curve>> flatDir1 = new List<List<Curve>>();
        List<List<Curve>> flatDir2 = new List<List<Curve>>();

        if (iterationState != null)
        {
            if (iterationState.CurvesDir1 != null)
            {
                for (int i = 0; i < iterationState.CurvesDir1.Count; i++)
                {
                    flatDir1.Add(iterationState.CurvesDir1[i]);
                    foreach (var cr in iterationState.CurvesDir1[i])
                        tree.Append(new GH_Curve(cr), new GH_Path(branchCounter));
                    branchCounter++;
                }
            }
            if (iterationState.CurvesDir2 != null)
            {
                for (int i = 0; i < iterationState.CurvesDir2.Count; i++)
                {
                    flatDir2.Add(iterationState.CurvesDir2[i]);
                    foreach (var cr in iterationState.CurvesDir2[i])
                        tree.Append(new GH_Curve(cr), new GH_Path(branchCounter));
                    branchCounter++;
                }
            }
        }

        // Order Logic
        if (flip == 2 || flip == 3)
        {
            int orderedBranchCounter = 0;
            for (int i = flatDir1.Count - 1; i > 0; i--)
            {
                foreach (var cr in flatDir1[i]) orderedTree.Append(new GH_Curve(cr), new GH_Path(orderedBranchCounter));
                orderedBranchCounter++;
            }
            if (flatDir1.Count > 0)
            {
                foreach (var cr in flatDir1[0]) orderedTree.Append(new GH_Curve(cr), new GH_Path(orderedBranchCounter));
                orderedBranchCounter++;
            }
            else if (flatDir2.Count > 0)
            {
                foreach (var cr in flatDir2[0]) orderedTree.Append(new GH_Curve(cr), new GH_Path(orderedBranchCounter));
                orderedBranchCounter++;
            }
            for (int i = 0; i < flatDir2.Count; i++)
            {
                foreach (var cr in flatDir2[i]) orderedTree.Append(new GH_Curve(cr), new GH_Path(orderedBranchCounter));
                orderedBranchCounter++;
            }
        }
        else
        {
            orderedTree = tree.Duplicate();
        }

        // --- Build Meshing Curves (with deeper branches {n;0} and {n;1}) ---
        GH_Structure<GH_Curve> meshTree;

        if (surface != null && orderedTree != null && orderedTree.Paths.Count > 0)
        {
            DataTree<Curve> processed = GenerateGlobalCornerCurves(orderedTree, surface);
            processed = InsertSplitBridgeCurves(orderedTree, processed, surface);
            meshTree = ToGhCurveTree(processed);
        }
        else
        {
            meshTree = orderedTree.Duplicate();
        }

        // Outputs
        DA.SetDataTree(0, tree);
        DA.SetDataTree(1, orderedTree);
        DA.SetDataTree(2, meshTree);
    }

    public void ToggleRun() => _runPressed = !_runPressed;
    public void ToggleReset() => _resetPressed = !_resetPressed;
    public override void CreateAttributes() => m_attributes = new CrvAttributes(this);

    // ==========================================================
    // PROCESSING LOGIC 
    // ==========================================================
    private void StartAutomaticProcessing(Surface surface, Curve curve, double accuracy, double distance, int flip)
    {
        if (_isProcessing) return;
        _isProcessing = true;
        Task.Run(() =>
        {
            while (!iterationState.IsComplete && _isProcessing)
            {
                List<List<Curve>> allCurvesHistory = new List<List<Curve>>();
                if (iterationState.CurvesDir1 != null) allCurvesHistory.AddRange(iterationState.CurvesDir1);
                if (iterationState.CurvesDir2 != null) allCurvesHistory.AddRange(iterationState.CurvesDir2);

                List<Curve> newCurves = ProcessNextIteration(surface, curve, accuracy, distance, flip, allCurvesHistory);

                if (newCurves.Count > 0)
                {
                    if (iterationState.CurrentK == 0)
                        iterationState.CurvesDir1.Add(new List<Curve>(newCurves));
                    else
                        iterationState.CurvesDir2.Add(new List<Curve>(newCurves));
                }
                else
                {
                    if (flip == 2 || flip == 3)
                    {
                        if (iterationState.CurrentK == 0 && iterationState.CurrentIteration > 0)
                        {
                            iterationState.CurrentK++;
                            iterationState.CurrentIteration = 0;
                            iterationState.BaseCurves = new List<NurbsCurve> { fixedInitialCurve.DuplicateCurve() as NurbsCurve };
                        }
                        else
                        {
                            iterationState.IsComplete = true;
                        }
                    }
                    else
                    {
                        iterationState.IsComplete = true;
                    }
                }
                Rhino.RhinoApp.InvokeOnUiThread((Action)delegate { this.ExpireSolution(true); });
                System.Threading.Thread.Sleep(50);
            }
            _isProcessing = false;
        });
    }

    private List<Curve> ProcessNextIteration(Surface surface, Curve curve, double accuracy, double distance, int flip, List<List<Curve>> allCurvesHistory)
    {
        switch (flip)
        {
            case 0: case 1: return ProcessSingleDirectionIteration(surface, curve, accuracy, distance, flip, allCurvesHistory);
            case 2: return ProcessLimitedDualDirectionIteration(surface, curve, accuracy, distance, allCurvesHistory);
            case 3: return ProcessUnlimitedDualDirectionIteration(surface, curve, accuracy, distance, allCurvesHistory);
            default: return new List<Curve>();
        }
    }

    private List<Curve> ProcessSingleDirectionIteration(Surface surface, Curve curve, double accuracy, double distance, int flip, List<List<Curve>> allCurvesHistory)
    {
        if (iterationState.CurrentIteration == 0)
        {
            var initCr = curve.ToNurbsCurve();
            initialCurve = ProcessInitialCurve(initCr, surface);
            if (fixedInitialCurve == null) fixedInitialCurve = initialCurve.DuplicateCurve() as NurbsCurve;
            iterationState.BaseCurves = new List<NurbsCurve> { initialCurve };
            iterationState.CurrentIteration++;
            return new List<Curve> { initialCurve };
        }
        if (iterationState.MaxIterations > 0 && iterationState.CurrentIteration > iterationState.MaxIterations)
        {
            iterationState.IsComplete = true;
            return new List<Curve>();
        }
        return GenerateCurvesFromBase(surface, distance, flip, accuracy, allCurvesHistory);
    }

    private List<Curve> ProcessLimitedDualDirectionIteration(Surface surface, Curve curve, double accuracy, double distance, List<List<Curve>> allCurvesHistory)
    {
        if (iterationState.CurrentK > 1) { iterationState.IsComplete = true; return new List<Curve>(); }
        bool stop = false;
        if (iterationState.CurrentK == 0 && iterationState.MaxIterations > 0 && iterationState.CurrentIteration > iterationState.MaxIterations) stop = true;
        if (iterationState.CurrentK == 1 && iterationState.MaxIterations > 0 && iterationState.CurrentIteration >= iterationState.MaxIterations) stop = true;

        if (stop) return new List<Curve>();
        List<Curve> result = ProcessDirection(surface, curve, accuracy, distance, iterationState.CurrentK, allCurvesHistory);
        if (result.Count > 0) iterationState.CurrentIteration++;
        return result;
    }

    private List<Curve> ProcessUnlimitedDualDirectionIteration(Surface surface, Curve curve, double accuracy, double distance, List<List<Curve>> allCurvesHistory)
    {
        if (iterationState.CurrentK > 1) { iterationState.IsComplete = true; return new List<Curve>(); }
        List<Curve> result = ProcessDirection(surface, curve, accuracy, distance, iterationState.CurrentK, allCurvesHistory);
        if (result.Count > 0) iterationState.CurrentIteration++;
        return result;
    }

    private List<Curve> ProcessDirection(Surface surface, Curve curve, double accuracy, double distance, int direction, List<List<Curve>> allCurvesHistory)
    {
        if (iterationState.CurrentIteration == 0 && direction == 0)
        {
            var initCr = curve.ToNurbsCurve();
            var baseC = ProcessInitialCurve(initCr, surface);
            iterationState.BaseCurves = new List<NurbsCurve> { baseC };
            fixedInitialCurve = baseC.DuplicateCurve() as NurbsCurve;
            return new List<Curve>(iterationState.BaseCurves);
        }
        return GenerateCurvesFromBase(surface, distance, direction, accuracy, allCurvesHistory);
    }

    private NurbsCurve ProcessInitialCurve(NurbsCurve initCr, Surface surface)
    {
        NurbsCurve result = initCr;
        if (ExtendCurves && !initCr.IsClosed)
        {
            double gap = initCr.PointAtStart.DistanceTo(initCr.PointAtEnd);
            if (gap > initCr.GetLength() * 0.1)
            {
                var ext = initCr.ExtendOnSurface(CurveEnd.Both, surface);
                if (ext != null) result = ext.ToNurbsCurve();
            }
        }
        var pulled = result.PullToBrepFace(Brep.CreateFromSurface(surface).Faces[0], 0.001);
        if (pulled != null && pulled.Length > 0) return pulled[0].ToNurbsCurve();
        return result;
    }

    private List<Curve> GenerateCurvesFromBase(Surface surface, double distance, int direction, double accuracy, List<List<Curve>> allCurvesHistory)
    {
        if (iterationState.BaseCurves == null || iterationState.BaseCurves.Count == 0) return new List<Curve>();
        List<(Curve nonExtended, Curve extended)> curvePairs = new List<(Curve, Curve)>();
        BrepFace face = Brep.CreateFromSurface(surface).Faces[0];
        double pullTol = 0.001;

        foreach (NurbsCurve baseCrv in iterationState.BaseCurves)
        {
            if (baseCrv == null) continue;
            bool isGeometricallyClosed = baseCrv.IsClosed || baseCrv.PointAtStart.DistanceTo(baseCrv.PointAtEnd) < distance * 0.5;

            if (isGeometricallyClosed)
            {
                double minLoopLength = 2 * Math.PI * distance * 1.1;
                if (baseCrv.GetLength() < minLoopLength) continue;
            }

            double tol = CalculateTolerance(baseCrv.GetLength(), accuracy);
            double len = baseCrv.GetLength();
            int cn = (int)Math.Ceiling(len / tol);
            if (cn < 5) cn = 10;

            List<List<Point3d>> allSegments = new List<List<Point3d>>();
            List<Point3d> currentSegment = new List<Point3d>();

            for (int i = 0; i < cn + 1; i++)
            {
                Point3d? pt = TryGetOffsetPoint(surface, distance, direction, baseCrv, len, cn, i);
                if (pt.HasValue) currentSegment.Add(pt.Value);
                else
                {
                    if (currentSegment.Count > 1) allSegments.Add(new List<Point3d>(currentSegment));
                    currentSegment.Clear();
                }
            }
            if (currentSegment.Count > 1) allSegments.Add(currentSegment);

            foreach (List<Point3d> rawPts in allSegments)
            {
                if (rawPts.Count < 3) continue;
                List<Point3d> segPts = RemoveKinks(rawPts, distance);
                if (segPts.Count < 3) continue;

                Curve smoothCrv = null;
                try { smoothCrv = Curve.CreateInterpolatedCurve(segPts, 3, CurveKnotStyle.Chord); } catch { continue; }
                if (smoothCrv == null || !smoothCrv.IsValid) continue;

                bool isClosedState = false;
                double gap = smoothCrv.PointAtStart.DistanceTo(smoothCrv.PointAtEnd);
                if (baseCrv.IsClosed || gap < distance * 1.5)
                {
                    if (smoothCrv.MakeClosed(distance * 1.5)) isClosedState = true;
                }

                var p1 = smoothCrv.PullToBrepFace(face, pullTol);
                Curve onSrfCrv = (p1 != null && p1.Length > 0) ? p1[0] : smoothCrv;
                Curve nonExtended = onSrfCrv.DuplicateCurve();

                Curve extendedCandidate = null;
                if (ExtendCurves && !isClosedState)
                {
                    double surfaceGap = onSrfCrv.PointAtStart.DistanceTo(onSrfCrv.PointAtEnd);
                    if (surfaceGap > distance * 1.0)
                    {
                        Curve extended = onSrfCrv.DuplicateCurve();
                        for (int extIter = 0; extIter < 10; extIter++)
                        {
                            var tempExt = extended.ExtendOnSurface(CurveEnd.Both, surface);
                            if (tempExt == null || tempExt.GetLength() <= extended.GetLength() + pullTol) break;
                            extended = tempExt;
                        }
                        var selfInt = Rhino.Geometry.Intersect.Intersection.CurveSelf(extended, pullTol);
                        if (selfInt == null || selfInt.Count == 0)
                        {
                            var p2 = extended.PullToBrepFace(face, pullTol);
                            extendedCandidate = (p2 != null && p2.Length > 0) ? p2[0] : extended;
                        }
                    }
                }
                curvePairs.Add((nonExtended, extendedCandidate));
            }
        }

        List<Curve> proposedExtended = curvePairs.Select(p => p.extended ?? p.nonExtended).ToList();
        bool hasSignificantIntraIntersection = false;
        for (int i = 0; i < proposedExtended.Count; i++)
        {
            for (int j = i + 1; j < proposedExtended.Count; j++)
            {
                var ccx = Rhino.Geometry.Intersect.Intersection.CurveCurve(proposedExtended[i], proposedExtended[j], pullTol, pullTol);
                if (ccx != null && ccx.Count > 0) { hasSignificantIntraIntersection = true; break; }
            }
            if (hasSignificantIntraIntersection) break;
        }

        List<Curve> selectedCurves = hasSignificantIntraIntersection ? curvePairs.Select(p => p.nonExtended).ToList() : proposedExtended;
        double joinTol = distance * 0.1;
        List<Curve> processedCurves = Curve.JoinCurves(selectedCurves, joinTol, false).ToList();

        List<Curve> allNewCurves = new List<Curve>();
        List<NurbsCurve> nextBaseCurves = new List<NurbsCurve>();

        foreach (Curve crv in processedCurves)
        {
            if (crv.GetLength() < distance * 0.1) continue;
            Curve finalCrv = crv;
            if (FilterIntersections && IsSignificantIntersection(finalCrv, allCurvesHistory)) continue;
            var selfInt = Rhino.Geometry.Intersect.Intersection.CurveSelf(finalCrv, pullTol);
            if (selfInt != null && selfInt.Count > 0) continue;

            allNewCurves.Add(finalCrv);
            nextBaseCurves.Add(finalCrv.ToNurbsCurve());
        }

        if (allCurvesHistory.Count > 0)
        {
            var prevList = allCurvesHistory.Last();
            if (prevList.Count > 0)
            {
                double avgPrevLength = prevList.Average(c => c.GetLength());
                List<Curve> filteredNew = new List<Curve>();
                foreach (var nc in allNewCurves)
                {
                    if (nc.GetLength() >= avgPrevLength * 0.05) filteredNew.Add(nc);
                }
                allNewCurves = filteredNew;
                nextBaseCurves = allNewCurves.Select(c => c.ToNurbsCurve()).ToList();
            }
        }

        List<Curve> uniqueNew = new List<Curve>();
        foreach (var nc in allNewCurves)
        {
            bool isDup = false;
            foreach (var uc in uniqueNew)
            {
                var ccx = Rhino.Geometry.Intersect.Intersection.CurveCurve(nc, uc, pullTol, pullTol);
                if (ccx != null && ccx.Any(e => e.IsOverlap)) { isDup = true; break; }
            }
            if (!isDup) uniqueNew.Add(nc);
        }
        allNewCurves = uniqueNew;
        nextBaseCurves = allNewCurves.Select(c => c.ToNurbsCurve()).ToList();

        iterationState.BaseCurves = nextBaseCurves;
        return allNewCurves;
    }

    private List<Point3d> RemoveKinks(List<Point3d> pts, double dist)
    {
        if (pts.Count < 3) return pts;
        List<Point3d> clean = new List<Point3d>();
        clean.Add(pts[0]);
        for (int i = 1; i < pts.Count - 1; i++)
        {
            Vector3d v1 = pts[i] - pts[i - 1];
            Vector3d v2 = pts[i + 1] - pts[i];
            if (v1.Length < 0.001 || v2.Length < 0.001) continue;
            double angle = Vector3d.VectorAngle(v1, v2);
            if (angle < (Math.PI * 0.8)) clean.Add(pts[i]);
        }
        clean.Add(pts[pts.Count - 1]);
        return clean;
    }

    private Point3d? TryGetOffsetPoint(Surface surface, double distance, int k, NurbsCurve initCr, double len, int cn, int i)
    {
        Point3d pt0 = (i == cn) ? initCr.PointAtEnd : initCr.PointAtLength(len / cn * i);
        double t; initCr.ClosestPoint(pt0, out t);
        double u, v; surface.ClosestPoint(pt0, out u, out v);
        Point3d pt1 = surface.PointAt(u, v);
        if (pt0.DistanceToSquared(pt1) > 0.01) return null;
        Vector3d norm = surface.NormalAt(u, v); norm.Unitize();
        Vector3d vS = Vector3d.CrossProduct(initCr.TangentAt(t), norm); vS.Unitize();
        if (k == 1) { vS = -vS; norm = -norm; }
        Curve arc = ArcCurve.CreateArcBlend(pt1 - norm * distance, -vS, pt1 + norm * distance, vS, 1.0);
        double[] tt;
        bool xx = Rhino.Geometry.Intersect.Intersection.CurveBrep(arc, Brep.CreateFromSurface(surface), 0.01, 0.01, out tt);
        if (xx && tt != null && tt.Length > 0) return arc.PointAt(tt[0]);
        return null;
    }

    private bool IsSignificantIntersection(Curve newCurve, List<List<Curve>> history)
    {
        if (history == null || history.Count == 0) return false;
        foreach (var list in history)
        {
            foreach (var oldCurve in list)
            {
                var ccx = Rhino.Geometry.Intersect.Intersection.CurveCurve(newCurve, oldCurve, 0.01, 0.01);
                if (ccx != null && ccx.Count > 0)
                {
                    bool allGrazing = true;
                    foreach (var evt in ccx)
                    {
                        if (evt.IsOverlap) return true;
                        bool atEndA = Math.Abs(evt.ParameterA - newCurve.Domain.Min) < 0.01 || Math.Abs(evt.ParameterA - newCurve.Domain.Max) < 0.01;
                        if (!atEndA) { allGrazing = false; break; }
                    }
                    if (!allGrazing) return true;
                }
            }
        }
        return false;
    }

    private double CalculateTolerance(double length, double accuracy)
    {
        double tol = 5.0;
        if (accuracy < 0) tol = 5;
        else if (accuracy > 1) tol = 50;
        else tol = accuracy * 45 + 5;
        return Math.Max(length / (tol * 1.5), 0.01);
    }
    public void ScheduleSolutionCallback(GH_Document doc) => this.ExpireSolution(false);
    class CurveId
    {
        public GH_Path Path;
        public int Index;
        public Curve CurveGeometry;
        public CurveId(GH_Path p, int i, Curve c) { Path = p; Index = i; CurveGeometry = c; }
    }

    class CornerMod
    {
        public List<Curve> BoundarySegments;
        public bool IsStartOfCurve;
    }
    class BoundaryEdge
    {
        // 0 = U-iso (constant V)
        // 1 = V-iso (constant U)
        public int Iso;
        public double IsoValue;
    }

    private DataTree<Curve> GenerateGlobalCornerCurves(GH_Structure<GH_Curve> curvesTree, Surface surface)
    {
        DataTree<Curve> resultTree = new DataTree<Curve>();

        Interval uDom = surface.Domain(0);
        Interval vDom = surface.Domain(1);

        // Collect all curves with ids
        List<CurveId> allCurves = new List<CurveId>();
        foreach (GH_Path path in curvesTree.Paths)
        {
            var branch = curvesTree.get_Branch(path);
            for (int i = 0; i < branch.Count; i++)
            {
                GH_Curve gh = branch[i] as GH_Curve;
                if (gh?.Value != null)
                    allCurves.Add(new CurveId(path, i, gh.Value));
            }
        }

        if (allCurves.Count == 0)
            return resultTree;

        // Determine propagation direction 
        Curve refCrv = allCurves[0].CurveGeometry;
        surface.ClosestPoint(refCrv.PointAtStart, out double u1, out double v1);
        surface.ClosestPoint(refCrv.PointAtEnd, out double u2, out double v2);

        bool curveRunsAlongU = Math.Abs(u1 - u2) > Math.Abs(v1 - v2);
        bool walkAlongU = !curveRunsAlongU; // propagation direction

        // Corner UVs
        Point2d[] cornersUV =
        {
        new Point2d(uDom.Min, vDom.Min),
        new Point2d(uDom.Max, vDom.Min),
        new Point2d(uDom.Max, vDom.Max),
        new Point2d(uDom.Min, vDom.Max)
    };

        const double tol = 1e-4;
        double joinTol = Rhino.RhinoDoc.ActiveDoc?.ModelAbsoluteTolerance ?? 1e-3;

        // key: "path|index" => list of corner mods
        var mods = new Dictionary<string, List<CornerMod>>();

        // =============================
        // PROCESS EACH CORNER
        // =============================
        for (int k = 0; k < 4; k++)
        {
            Point2d cornerUV = cornersUV[k];
            Point3d cornerPt = surface.PointAt(cornerUV.X, cornerUV.Y);

            // Define the two boundary edges at this corner
            BoundaryEdge edgeU, edgeV;

            if (k == 0)
            {
                edgeU = new BoundaryEdge { Iso = 0, IsoValue = vDom.Min };
                edgeV = new BoundaryEdge { Iso = 1, IsoValue = uDom.Min };
            }
            else if (k == 1)
            {
                edgeU = new BoundaryEdge { Iso = 0, IsoValue = vDom.Min };
                edgeV = new BoundaryEdge { Iso = 1, IsoValue = uDom.Max };
            }
            else if (k == 2)
            {
                edgeU = new BoundaryEdge { Iso = 0, IsoValue = vDom.Max };
                edgeV = new BoundaryEdge { Iso = 1, IsoValue = uDom.Max };
            }
            else
            {
                edgeU = new BoundaryEdge { Iso = 0, IsoValue = vDom.Max };
                edgeV = new BoundaryEdge { Iso = 1, IsoValue = uDom.Min };
            }

            int countU = CountEndpointsOnEdge(edgeU);
            int countV = CountEndpointsOnEdge(edgeV);

            BoundaryEdge chosenEdge;
            if (countU < countV) chosenEdge = edgeU;
            else if (countV < countU) chosenEdge = edgeV;
            else chosenEdge = walkAlongU ? edgeU : edgeV;

            // Find closest endpoint ON chosen edge
            CurveId bestCurve = null;
            bool isStart = true;
            double bestDist = double.MaxValue;

            foreach (var cid in allCurves)
            {
                foreach (var pair in new[]
                {
                (cid.CurveGeometry.PointAtStart, true),
                (cid.CurveGeometry.PointAtEnd, false)
            })
                {
                    surface.ClosestPoint(pair.Item1, out double u, out double v);

                    bool onEdge =
                        (chosenEdge.Iso == 0 && Math.Abs(v - chosenEdge.IsoValue) < tol) ||
                        (chosenEdge.Iso == 1 && Math.Abs(u - chosenEdge.IsoValue) < tol);

                    if (!onEdge) continue;

                    double d = pair.Item1.DistanceTo(cornerPt);
                    if (d < bestDist)
                    {
                        bestDist = d;
                        bestCurve = cid;
                        isStart = pair.Item2;
                    }
                }
            }

            if (bestCurve == null)
                continue;

            // Create boundary segment ONLY along chosen edge
            Curve boundary = surface.IsoCurve(chosenEdge.Iso, chosenEdge.IsoValue);

            surface.ClosestPoint(
                isStart ? bestCurve.CurveGeometry.PointAtStart
                        : bestCurve.CurveGeometry.PointAtEnd,
                out double uH, out double vH);

            double tCorner = (chosenEdge.Iso == 0) ? cornerUV.X : cornerUV.Y;
            double tHit = (chosenEdge.Iso == 0) ? uH : vH;

            Curve seg = boundary.Trim(new Interval(
                Math.Min(tCorner, tHit),
                Math.Max(tCorner, tHit)));

            if (seg == null) continue;

            if (seg.PointAtStart.DistanceTo(cornerPt) >
                seg.PointAtEnd.DistanceTo(cornerPt))
                seg.Reverse();

            Curve crv = bestCurve.CurveGeometry;

            Point3d P = isStart ? crv.PointAtStart : crv.PointAtEnd;
            Point3d B = isStart ? seg.PointAtStart : seg.PointAtEnd;

            double tInner = isStart
                ? crv.Domain.ParameterAt(0.1)
                : crv.Domain.ParameterAt(0.9);

            Point3d C = crv.PointAt(tInner);

            Vector3d vBoundary = B - P;
            Vector3d vCurve = C - P;

            bool shouldJoin = false;
            if (vBoundary.Unitize() && vCurve.Unitize())
            {
                double ang = Vector3d.VectorAngle(vBoundary, vCurve);
                shouldJoin = (ang > Math.PI / 2.0);
            }

            // Store join intent inside the segment 
            seg.UserDictionary.Set("PB_ShouldJoin", shouldJoin);

            CornerMod mod = new CornerMod
            {
                BoundarySegments = new List<Curve> { seg },
                IsStartOfCurve = isStart
            };

            string key = bestCurve.Path.ToString() + "|" + bestCurve.Index;
            if (!mods.ContainsKey(key))
                mods[key] = new List<CornerMod>();

            mods[key].Add(mod);
        }

        // =============================
        // RECONSTRUCT TREE
        // Keep original at {path}
        // Put joined refinement at {path;0}
        // Put unjoined border segments at {path;1}
        // =============================
        foreach (GH_Path path in curvesTree.Paths)
        {
            var branch = curvesTree.get_Branch(path);
            if (branch == null) continue;

            for (int i = 0; i < branch.Count; i++)
            {
                Curve original = (branch[i] as GH_Curve)?.Value;
                if (original == null) continue;

                // Always keep original
                resultTree.Add(original, path);

                string key = path.ToString() + "|" + i;
                if (!mods.ContainsKey(key))
                    continue;

                Curve refined = original.DuplicateCurve();
                var extras = new List<Curve>();

                foreach (var mod in mods[key])
                {
                    Curve seg = mod.BoundarySegments[0];

                    bool shouldJoin = false;
                    if (seg.UserDictionary.ContainsKey("PB_ShouldJoin"))
                        shouldJoin = seg.UserDictionary.GetBool("PB_ShouldJoin");

                    if (shouldJoin)
                    {
                        Curve[] joined = mod.IsStartOfCurve
                            ? Curve.JoinCurves(new[] { seg, refined }, joinTol)
                            : Curve.JoinCurves(new[] { refined, seg }, joinTol);

                        if (joined != null && joined.Length > 0)
                            refined = joined[0];
                        else
                            extras.Add(seg);
                    }
                    else
                    {
                        extras.Add(seg);
                    }
                }

                // If refinement succeeded (heuristic: length changed), add to {path;0}
                
                if (refined != null)
                {
                    double L0 = original.GetLength();
                    double L1 = refined.GetLength();

                    if (Math.Abs(L1 - L0) > 1e-6)
                    {
                        GH_Path refinedPath;

                        if (path.Length == 1)
                        {
                            int n = path[0];
                            int targetN = Math.Max(0, n - 1);
                            refinedPath = new GH_Path(targetN, 0);   // ✅ corner: {n-1;0}
                        }
                        else
                        {
                            refinedPath = new GH_Path(path).AppendElement(0);
                        }

                        resultTree.Add(refined, refinedPath);
                    }
                }
                if (extras.Count > 0)
                {
                    GH_Path extraPath = new GH_Path(path).AppendElement(1); // {n;1}
                    resultTree.AddRange(extras, extraPath);
                }
            }
        }

        return resultTree;

        // =============================
        // LOCAL HELPERS
        // =============================
        int CountEndpointsOnEdge(BoundaryEdge edge)
        {
            int count = 0;
            foreach (var cid in allCurves)
            {
                foreach (var pt in new[]
                {
                cid.CurveGeometry.PointAtStart,
                cid.CurveGeometry.PointAtEnd
            })
                {
                    surface.ClosestPoint(pt, out double u, out double v);
                    if (edge.Iso == 0 && Math.Abs(v - edge.IsoValue) < tol) count++;
                    if (edge.Iso == 1 && Math.Abs(u - edge.IsoValue) < tol) count++;
                }
            }
            return count;
        }
    }


    // --- HELPER METHODS ---
    private class BoundaryHit
    {
        public int CurveIndex;        // index in the black curve list
        public bool IsStart;          // endpoint type
        public Point3d Pt;
        public int Iso;               // 0 = U-iso (constant V), 1 = V-iso (constant U)
        public double IsoValue;       // vMin/vMax or uMin/uMax
        public double T;              // parameter along that iso curve (u or v)
    }

    private DataTree<Curve> InsertSplitBridgeCurves(
    GH_Structure<GH_Curve> inputTree,
    DataTree<Curve> processedTree,
    Surface surface)
    {
        // Copy processedTree into outTree
        DataTree<Curve> outTree = new DataTree<Curve>();
        foreach (var p in processedTree.Paths)
            outTree.AddRange(processedTree.Branch(p), p);

        Interval uDom = surface.Domain(0);
        Interval vDom = surface.Domain(1);
        const double tol = 1e-4;
        double joinTol = Rhino.RhinoDoc.ActiveDoc?.ModelAbsoluteTolerance ?? 1e-3;

        IList<GH_Path> paths = inputTree.Paths;

        for (int pi = 0; pi < paths.Count - 1; pi++)
        {
            GH_Path pGreen = paths[pi];
            GH_Path pBlack = paths[pi + 1];

            var greenBranch = inputTree.get_Branch(pGreen);
            var blackBranch = inputTree.get_Branch(pBlack);

            if (greenBranch == null || blackBranch == null) continue;
            if (greenBranch.Count != 1) continue;
            if (blackBranch.Count <= 1) continue;

            
            Curve greenCurve = (greenBranch[0] as GH_Curve)?.Value;
            if (greenCurve == null) continue;

            // Check if the parent curve was originally open
            bool parentWasOpen = !greenCurve.IsClosed &&
                                greenCurve.PointAtStart.DistanceTo(greenCurve.PointAtEnd) > joinTol * 2;

            var blacks = new List<Curve>();
            for (int i = 0; i < blackBranch.Count; i++)
            {
                var c = (blackBranch[i] as GH_Curve)?.Value;
                if (c != null) blacks.Add(c);
            }
            if (blacks.Count <= 1) continue;

            // Find boundary endpoints for black curves
            var hits = new List<BoundaryHit>();

            for (int i = 0; i < blacks.Count; i++)
            {
                var c = blacks[i];
                foreach (var ep in new[] { (c.PointAtStart, true), (c.PointAtEnd, false) })
                {
                    surface.ClosestPoint(ep.Item1, out double u, out double v);

                    bool onUMin = Math.Abs(u - uDom.Min) < tol;
                    bool onUMax = Math.Abs(u - uDom.Max) < tol;
                    bool onVMin = Math.Abs(v - vDom.Min) < tol;
                    bool onVMax = Math.Abs(v - vDom.Max) < tol;

                    if (!(onUMin || onUMax || onVMin || onVMax))
                        continue;

                    // classify by closest boundary
                    double duMin = Math.Abs(u - uDom.Min);
                    double duMax = Math.Abs(u - uDom.Max);
                    double dvMin = Math.Abs(v - vDom.Min);
                    double dvMax = Math.Abs(v - vDom.Max);

                    double best = Math.Min(Math.Min(duMin, duMax), Math.Min(dvMin, dvMax));

                    int iso;
                    double isoValue;
                    double t;

                    if (best == dvMin)
                    {
                        iso = 0; isoValue = vDom.Min; t = u;
                    }
                    else if (best == dvMax)
                    {
                        iso = 0; isoValue = vDom.Max; t = u;
                    }
                    else if (best == duMin)
                    {
                        iso = 1; isoValue = uDom.Min; t = v;
                    }
                    else
                    {
                        iso = 1; isoValue = uDom.Max; t = v;
                    }

                    hits.Add(new BoundaryHit
                    {
                        CurveIndex = i,
                        IsStart = ep.Item2,
                        Pt = surface.PointAt(u, v),
                        Iso = iso,
                        IsoValue = isoValue,
                        T = t
                    });
                }
            }

            
            HashSet<int> originalOpenEndIndices = new HashSet<int>();

            if (parentWasOpen && greenCurve != null)
            {
                Point3d parentStart = greenCurve.PointAtStart;
                Point3d parentEnd = greenCurve.PointAtEnd;

                
                double closestStartDist = double.MaxValue;
                double closestEndDist = double.MaxValue;
                int startHitIdx = -1;
                int endHitIdx = -1;

                for (int i = 0; i < hits.Count; i++)
                {
                    double distToStart = hits[i].Pt.DistanceTo(parentStart);
                    double distToEnd = hits[i].Pt.DistanceTo(parentEnd);

                    if (distToStart < closestStartDist)
                    {
                        closestStartDist = distToStart;
                        startHitIdx = i;
                    }

                    if (distToEnd < closestEndDist)
                    {
                        closestEndDist = distToEnd;
                        endHitIdx = i;
                    }
                }

                if (startHitIdx >= 0) originalOpenEndIndices.Add(startHitIdx);
                if (endHitIdx >= 0) originalOpenEndIndices.Add(endHitIdx);
            }

            // Group hits by boundary edge (iso + isoValue)
            var groups = hits
                .GroupBy(h => (h.Iso, h.IsoValue))
                .ToList();

            var bridgeSegments = new List<Curve>();

            foreach (var g in groups)
            {
                var list = g.OrderBy(h => h.T).ToList();
                if (list.Count < 2) continue;

                // Pair consecutive hits into boundary segments: [0-1], [2-3], ...
                for (int i = 0; i + 1 < list.Count; i += 2)
                {
                    int globalIdx0 = hits.IndexOf(list[i]);
                    int globalIdx1 = hits.IndexOf(list[i + 1]);
                    if (originalOpenEndIndices.Contains(globalIdx0) ||
                        originalOpenEndIndices.Contains(globalIdx1))
                        continue;

                    var h0 = list[i];
                    var h1 = list[i + 1];

                    if (h0.CurveIndex == h1.CurveIndex) continue;

                    Curve boundary = surface.IsoCurve(h0.Iso, h0.IsoValue);
                    Curve seg = boundary.Trim(new Interval(Math.Min(h0.T, h1.T), Math.Max(h0.T, h1.T)));
                    if (seg == null) continue;

                    bridgeSegments.Add(seg);
                }
            }

            if (bridgeSegments.Count == 0)
                continue;

            // Join blacks + bridges
            var pieces = new List<Curve>();
            pieces.AddRange(blacks);
            pieces.AddRange(bridgeSegments);

            var joined = Curve.JoinCurves(pieces, joinTol);
            if (joined == null || joined.Length == 0)
                continue;

            int n = pGreen[0];
            int targetN = Math.Max(0, n - 1);
            GH_Path newPath = new GH_Path(pGreen).AppendElement(0);
            outTree.Add(joined[0], newPath);
        }

        return outTree;
    }
    private GH_Structure<GH_Curve> ToGhCurveTree(DataTree<Curve> dt)
    {
        var t = new GH_Structure<GH_Curve>();
        foreach (GH_Path p in dt.Paths)
        {
            var br = dt.Branch(p);
            if (br == null) continue;
            foreach (var c in br)
            {
                if (c != null) t.Append(new GH_Curve(c), p);
            }
        }
        return t;
    }

    // ================== ATTRIBUTES & STATE ==================
    public class IterationState
    {
        public int CurrentIteration { get; set; } = 0;
        public int MaxIterations { get; set; }
        public int CurrentK { get; set; } = 0;
        public int Mode { get; set; }
        public List<NurbsCurve> BaseCurves { get; set; } = new List<NurbsCurve>();
        public bool IsComplete { get; set; } = false;
        public List<List<Curve>> CurvesDir1 { get; set; } = new List<List<Curve>>();
        public List<List<Curve>> CurvesDir2 { get; set; } = new List<List<Curve>>();
        public IterationState(int maxIter, int mode) { MaxIterations = maxIter; Mode = mode; }
    }

    private class CrvAttributes : GH_ComponentAttributes
    {
        private readonly Crv owner;
        public CrvAttributes(Crv owner) : base(owner) { this.owner = owner; }

        protected override void Layout()
        {
            base.Layout();
            Rectangle bounds = GH_Convert.ToRectangle(Bounds);
            bounds.Height += 18 + 10 + 10 + 4 + 2 + 3;
            Bounds = bounds;
        }

        public override GH_ObjectResponse RespondToMouseDown(GH_Canvas sender, GH_CanvasMouseEvent e)
        {
            if (e.Button != MouseButtons.Left) return base.RespondToMouseDown(sender, e);
            System.Drawing.Point pt = new System.Drawing.Point((int)e.CanvasLocation.X, (int)e.CanvasLocation.Y);

            int buttonHeight = 18; int toggleHeight = 10; int gap = 3; int rowGap1 = 4; int rowGap2 = 2;
            int row3Y = (int)Bounds.Bottom - toggleHeight - gap;
            int row2Y = row3Y - toggleHeight - rowGap2;
            int row1Y = row2Y - buttonHeight - rowGap1;

            Rectangle runRect = new Rectangle((int)Bounds.X + gap, row1Y, buttonHeight, buttonHeight);
            Rectangle resetRect = new Rectangle(runRect.Right + 2, row1Y, (int)Bounds.Width - gap * 2 - buttonHeight - 2, buttonHeight);
            Rectangle extRect = new Rectangle((int)Bounds.X + gap, row2Y, buttonHeight, toggleHeight);
            Rectangle filtRect = new Rectangle((int)Bounds.X + gap, row3Y, buttonHeight, toggleHeight);

            if (runRect.Contains(pt))
            {
                owner.RunClicked = true;
                owner.ToggleRun();
                owner.ExpireSolution(true);
                Timer t = new Timer { Interval = 200 };
                t.Tick += (s, args) => { owner.RunClicked = false; t.Stop(); sender.Refresh(); };
                t.Start();
                return GH_ObjectResponse.Handled;
            }
            if (resetRect.Contains(pt))
            {
                owner.ResetClicked = true;
                owner.ToggleReset();
                owner.ExpireSolution(true);
                Timer t = new Timer { Interval = 200 };
                t.Tick += (s, args) => { owner.ResetClicked = false; t.Stop(); sender.Refresh(); };
                t.Start();
                return GH_ObjectResponse.Handled;
            }
            if (extRect.Contains(pt))
            {
                owner.ExtendCurves = !owner.ExtendCurves;
                owner.ExpireSolution(true);
                sender.Refresh();
                return GH_ObjectResponse.Handled;
            }
            if (filtRect.Contains(pt))
            {
                owner.FilterIntersections = !owner.FilterIntersections;
                owner.ExpireSolution(true);
                sender.Refresh();
                return GH_ObjectResponse.Handled;
            }
            return base.RespondToMouseDown(sender, e);
        }

        protected override void Render(GH_Canvas canvas, Graphics g, GH_CanvasChannel channel)
        {
            base.Render(canvas, g, channel);
            if (channel != GH_CanvasChannel.Objects) return;

            int buttonHeight = 18; int toggleHeight = 10; int gap = 3; int rowGap1 = 4; int rowGap2 = 2;
            int row3Y = (int)Bounds.Bottom - toggleHeight - gap;
            int row2Y = row3Y - toggleHeight - rowGap2;
            int row1Y = row2Y - buttonHeight - rowGap1;

            Rectangle runRect = new Rectangle((int)Bounds.X + gap, row1Y, buttonHeight, buttonHeight);
            Rectangle resetRect = new Rectangle(runRect.Right + 2, row1Y, (int)Bounds.Width - gap * 2 - buttonHeight - 2, buttonHeight);
            Rectangle extRect = new Rectangle((int)Bounds.X + gap, row2Y, buttonHeight, toggleHeight);
            Rectangle filtRect = new Rectangle((int)Bounds.X + gap, row3Y, buttonHeight, toggleHeight);

            
            DrawFlatButton(g, runRect, owner.RunClicked, true); // Play Icon inside
            DrawFlatButton(g, resetRect, owner.ResetClicked, false, "Reset");

            // Toggles
            DrawFlatToggle(g, extRect, owner.ExtendCurves);
            DrawFlatToggle(g, filtRect, owner.FilterIntersections);

            // Labels for Toggles
            DrawLabel(g, new Rectangle(extRect.Right + gap, row2Y, (int)Bounds.Right - extRect.Right - gap, toggleHeight), "Extend Curves");
            DrawLabel(g, new Rectangle(filtRect.Right + gap, row3Y, (int)Bounds.Right - filtRect.Right - gap, toggleHeight), "Filter Intersections");
        }

        private void DrawLabel(Graphics g, Rectangle bounds, string text)
        {
            using (Font f = new Font("Arial", 4f, FontStyle.Bold))
            using (Brush b = new SolidBrush(Color.FromArgb(50, 50, 50))) // Dark Gray Text
                g.DrawString(text, f, b, bounds, new StringFormat { Alignment = StringAlignment.Near, LineAlignment = StringAlignment.Center });
        }

        private void DrawFlatButton(Graphics g, Rectangle r, bool clicked, bool isPlayIcon = false, string text = "")
        {
            // Color Scheme
            Color normalFill = Color.FromArgb(225, 225, 225); // Light Gray
            Color activeFill = Color.Orange;
            Color borderColor = Color.Gray;

            using (Brush b = new SolidBrush(clicked ? activeFill : normalFill))
            using (GraphicsPath path = RoundedRect(r, r.Height / 4))
            {
                g.FillPath(b, path);
                using (Pen p = new Pen(borderColor)) g.DrawPath(p, path);
            }

            if (isPlayIcon)
            {
                // Draw Play Triangle
                System.Drawing.Point[] tri = { new System.Drawing.Point(r.X + 5, r.Y + 4), new System.Drawing.Point(r.X + 5, r.Bottom - 4), new System.Drawing.Point(r.Right - 5, r.Y + r.Height / 2) };
                using (Brush b = new SolidBrush(clicked ? Color.White : Color.DarkGray)) g.FillPolygon(b, tri);
            }
            else
            {
                using (Font f = new Font("Arial", 7f, FontStyle.Bold))
                using (Brush b = new SolidBrush(clicked ? Color.White : Color.DimGray))
                {
                    g.DrawString(text, f, b, r, new StringFormat { Alignment = StringAlignment.Center, LineAlignment = StringAlignment.Center });
                }
            }
        }

        private void DrawFlatToggle(Graphics g, Rectangle r, bool on)
        {
            // Toggle: Gray when Off, Orange when On
            Color onColor = Color.Orange;
            Color offColor = Color.FromArgb(225, 225, 225);
            Color borderColor = Color.Gray;

            using (GraphicsPath path = RoundedRect(r, r.Height / 2))
            using (Brush b = new SolidBrush(on ? onColor : offColor))
            {
                g.FillPath(b, path);
                using (Pen p = new Pen(borderColor)) g.DrawPath(p, path);
            }

            // Knob
            int kS = r.Height - 4;
            Rectangle kR = new Rectangle(on ? r.Right - kS - 2 : r.X + 2, r.Y + 2, kS, kS);

            using (GraphicsPath kP = RoundedRect(kR, kS / 2))
            {
                g.FillPath(Brushes.White, kP);
            }
        }

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
}