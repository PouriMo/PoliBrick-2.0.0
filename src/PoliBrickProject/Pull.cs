using System;
using System.Collections.Generic;
using System.Drawing;
using System.IO;
using Grasshopper.Kernel;
using Rhino.Geometry;

namespace PoliBrick.Util
{
    public class Pull : GH_Component
    {
        public Pull()
          : base("1.Pull", "Pull",
              "Pull a curve to a surface",
              "PoliBrick", "Utilities")
        {
        }
        public override GH_Exposure Exposure
        {
            get { return GH_Exposure.primary; }
        }
        protected override void RegisterInputParams(GH_InputParamManager pManager)
        {
            pManager.AddSurfaceParameter("Surface", "S", "Surface to pull curve onto", GH_ParamAccess.item);
            pManager.AddCurveParameter("Curve", "C", "Curve to pull", GH_ParamAccess.item);
        }

        protected override void RegisterOutputParams(GH_OutputParamManager pManager)
        {
            pManager.AddCurveParameter("Curves", "C", "Pulled curves", GH_ParamAccess.list);
        }

        protected override void SolveInstance(IGH_DataAccess DA)
        {
            // Get input data with proper type checking
            Surface surface = null;
            if (!DA.GetData(0, ref surface) || surface == null)
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "Invalid surface input");
                return;
            }

            Curve curve = null;
            if (!DA.GetData(1, ref curve) || curve == null)
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "Invalid curve input");
                return;
            }

            // Process the data
            List<Curve> result = ProcessCurvePull(surface, curve);

            // Set output
            if (result != null && result.Count > 0)
            {
                DA.SetDataList(0, result);
            }
            else
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Warning, "No curves were generated");
            }
        }

        private List<Curve> ProcessCurvePull(Surface surface, Curve curve)
        {
            var result = new List<Curve>();
            try
            {
                Brep b = Brep.CreateFromSurface(surface);

                if (!IsCurveOnSurface(curve, surface, 0.01, 10))
                {
                    Point3d pt = curve.PointAtNormalizedLength(0.5);
                    double u, v;
                    if (surface.ClosestPoint(pt, out u, out v))
                    {
                        Vector3d vec = pt - surface.PointAt(u, v);
                        Curve[] projectedCurves = Curve.ProjectToBrep(curve, b, vec, 0.01);
                        if (projectedCurves != null)
                        {
                            result.AddRange(projectedCurves);
                        }
                    }
                }
                else
                {
                    Curve[] pulledCurves = curve.PullToBrepFace(b.Faces[0], 0.01);
                    if (pulledCurves != null)
                    {
                        result.AddRange(pulledCurves);
                    }
                }
            }
            catch (Exception ex)
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Error, ex.Message);
            }

            return result;
        }

        private bool IsCurveOnSurface(Curve crv, Surface srf, double tolerance, int samples)
        {
            if (crv == null || srf == null) return false;

            for (int i = 0; i < samples; i++)
            {
                double t = crv.Domain.ParameterAt((double)i / (samples - 1));
                Point3d pt = crv.PointAt(t);
                double u, v;
                if (!srf.ClosestPoint(pt, out u, out v))
                    return false;

                Point3d srfPt = srf.PointAt(u, v);
                if (pt.DistanceTo(srfPt) > tolerance)
                    return false;
            }
            return true;
        }

        protected override Bitmap Icon =>
            ConvertByteArrayToBitmap(PoliBrick.Properties.Resources.pull);

        public static Bitmap ConvertByteArrayToBitmap(byte[] imageData)
        {
            using (var ms = new MemoryStream(imageData))
            {
                return new Bitmap(ms);
            }
        }

        public override Guid ComponentGuid =>
            new Guid("CB7F1D28-8641-4C88-8F1F-201B9584DD57");
    }
}
    