using System;
using System.Drawing;
using System.Collections.Generic;
using Grasshopper.Kernel;
using Grasshopper.Kernel.Types;
using System.IO;

namespace PpliBrick.Util
{
    public class BrData : GH_Component
    {
        public BrData()
          : base("3.BrickData", "ID",
              "The component assigns dimensional properties, transformation values (α,β,θ), and color to specific brick identifiers.",
              "PoliBrick", "Utilities")
        {
        }

        public override GH_Exposure Exposure
        {
            get { return GH_Exposure.primary; }
        }

        protected override Bitmap Icon => ConvertByteArrayToBitmap(PoliBrick.Properties.Resources.BrDataIcon);

        public override Guid ComponentGuid => new Guid("5B187761-58B0-454C-A215-C0215E97E86C");

        protected override void RegisterInputParams(GH_InputParamManager pManager)
        {
            pManager.AddTextParameter("Name", "n", "The assigned name of the brick", GH_ParamAccess.item);
            pManager.AddTextParameter("u,v,w,m", "uvwm", "Comma-separated brick(length,width,height) and mortar (thickness) dimensions (u,v,w,m)", GH_ParamAccess.item);
            pManager.AddTextParameter("α,β,θ", "αβθ", "Comma-separated Alpha, Beta, Theta values. Three rotation values.", GH_ParamAccess.item);
            pManager.AddTextParameter("Rotation Pattern", "RotPatt", "Pattern for rotations of bricks (T/F)", GH_ParamAccess.item);
            pManager.AddColourParameter("Color", "Col", "Optional color for the brick", GH_ParamAccess.item);

            pManager[2].Optional = true;
            pManager[3].Optional = true;
            pManager[4].Optional = true; // Color is optional
        }

        protected override void RegisterOutputParams(GH_OutputParamManager pManager)
        {
            pManager.AddTextParameter("Data", "D", "Brick data output as a list", GH_ParamAccess.list);
        }

        protected override void SolveInstance(IGH_DataAccess DA)
        {
            string name = string.Empty;
            string uvw_Mortar = string.Empty;
            string alfaBetaTheta = null;
            string rotPattern = null;
            Color brickColor = Color.Empty; 

            if (!DA.GetData(0, ref name)) return;
            if (!DA.GetData(1, ref uvw_Mortar)) return;
            DA.GetData(2, ref alfaBetaTheta);
            DA.GetData(3, ref rotPattern);
            DA.GetData(4, ref brickColor);

            try
            {
                var dataList = ProcessPatternData(name, uvw_Mortar, alfaBetaTheta, rotPattern, brickColor);
                DA.SetDataList(0, dataList);
            }
            catch (ArgumentException ex)
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Warning, ex.Message);
            }
            catch (Exception ex)
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Error, ex.Message);
            }
        }

        private List<string> ProcessPatternData(string name, string uvw_Mortar, string alfaBetaTheta, string rotPattern, Color brickColor)
        {
            // --- 1. Validate Dimensions ---
            string[] coordinates = uvw_Mortar.Split(new[] { ',' }, StringSplitOptions.RemoveEmptyEntries);
            if (coordinates.Length != 4)
                throw new ArgumentException("Invalid 'u,v,w,m' format. Expected 4 comma-separated values.");

            if (!double.TryParse(coordinates[0], out _) ||
                !double.TryParse(coordinates[1], out _) ||
                !double.TryParse(coordinates[2], out _) ||
                !double.TryParse(coordinates[3], out _))
                throw new ArgumentException("Invalid numeric values in 'u,v,w,m'.");

            // --- 2. Validate Rotations (α, β, θ) ---
            // Default to 0,0,0 if input is missing
            alfaBetaTheta = alfaBetaTheta ?? "0,0,0";
            string[] rotationValues = alfaBetaTheta.Split(new[] { ',' }, StringSplitOptions.RemoveEmptyEntries);

            if (rotationValues.Length != 3)
                throw new ArgumentException("Invalid 'α,β,θ' format. Expected 3 comma-separated values.");

            if (!double.TryParse(rotationValues[0], out _) ||
                !double.TryParse(rotationValues[1], out _) ||
                !double.TryParse(rotationValues[2], out _))
                throw new ArgumentException("Invalid numeric values in 'α,β,θ'.");

            // --- 3. Validate Pattern ---
            
            rotPattern = rotPattern ?? "T";

            if (!IsValidPattern(rotPattern))
                throw new ArgumentException("Invalid Rotation Pattern format. Expected 'T' or 'F' only.");

            // --- 4. Process Color ---
            
            string colorStr = (brickColor == Color.Empty) ? "Empty" : $"{brickColor.A},{brickColor.R},{brickColor.G},{brickColor.B}";

            // --- 5. Compile Output List ---
            return new List<string>
            {
                name,           // Index 0
                uvw_Mortar,     // Index 1
                alfaBetaTheta,  // Index 2 (α, β, θ)
                rotPattern,     // Index 3 (Rotation Pattern)
                colorStr        // Index 4 (Color)
            };
        }

        public static Bitmap ConvertByteArrayToBitmap(byte[] imageData)
        {
            using (var ms = new MemoryStream(imageData))
            {
                return new Bitmap(ms);
            }
        }

        private bool IsValidPattern(string pattern)
        {
            foreach (char c in pattern)
            {
                if (c != 'T' && c != 'F')
                    return false;
            }
            return true;
        }
    }
}