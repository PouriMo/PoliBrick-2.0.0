using System;
using System.Drawing;
using System.Collections.Generic;
using System.Text;
using Grasshopper.Kernel;
using Grasshopper.GUI;
using Grasshopper.GUI.Canvas;
using Grasshopper.Kernel.Attributes;
using System.IO;
using System.Drawing.Drawing2D;
using System.Windows.Forms;

namespace PB.Util
{
    public class PattComponent : GH_Component
    {
        private bool infiniteMode = false;

        public PattComponent()
          : base("4.PatternGenerator", "PattGen",
            "This component generates an alphabetical pattern as a string.",
            "PoliBrick", "Utilities")
        {
        }
        public override GH_Exposure Exposure
        {
            get { return GH_Exposure.primary; }
        }
        public override bool Write(GH_IO.Serialization.GH_IWriter writer)
        {
            writer.SetBoolean("InfiniteMode", infiniteMode);
            return base.Write(writer);
        }

        public override bool Read(GH_IO.Serialization.GH_IReader reader)
        {
            if (reader.ItemExists("InfiniteMode"))
                infiniteMode = reader.GetBoolean("InfiniteMode");
            return base.Read(reader);
        }

        protected override void RegisterInputParams(GH_InputParamManager pManager)
        {
            pManager.AddTextParameter("Pattern", "P", "String pattern to repeat", GH_ParamAccess.item);
            pManager.AddIntegerParameter("Repeat", "R", "Number of times to repeat", GH_ParamAccess.item, 1);
            pManager.AddTextParameter("Plus", "+", "Additional string to append", GH_ParamAccess.item, "Empty!");
        }

        protected override void RegisterOutputParams(GH_OutputParamManager pManager)
        {
            pManager.AddTextParameter("Pattern", "P", "Generated pattern output", GH_ParamAccess.item);
        }

        protected override void SolveInstance(IGH_DataAccess DA)
        {
            string iPattern = "";
            int repeat = 1;
            string plusString = "Empty!";

            if (!DA.GetData(0, ref iPattern)) return;
            if (!DA.GetData(1, ref repeat)) return;
            if (!DA.GetData(2, ref plusString)) return;

            var builder = new StringBuilder();

            if (!string.IsNullOrEmpty(iPattern))
            {
                for (int i = 0; i < repeat; i++)
                {
                    builder.Append(iPattern);
                }
            }

            if (plusString != "Empty!")
            {
                builder.Append(plusString);
            }

            string output = builder.ToString();

            if (!string.IsNullOrEmpty(output))
            {
                DA.SetData(0, infiniteMode ? $"__{output}__" : output);
            }
            else
            {
                DA.SetData(0, "");
            }
        }

        public bool InfiniteMode
        {
            get { return infiniteMode; }
            set { infiniteMode = value; }
        }

        public override void CreateAttributes()
        {
            m_attributes = new PattComponentAttributes(this);
        }

        public static Bitmap ConvertByteArrayToBitmap(byte[] imageData)
        {
            if (imageData == null || imageData.Length == 0) return null;
            using (MemoryStream ms = new MemoryStream(imageData))
            {
                return new Bitmap(ms);
            }
        }

        protected override Bitmap Icon => ConvertByteArrayToBitmap(PoliBrick.Properties.Resources.PattIcon);

        public override Guid ComponentGuid => new Guid("b28dbf73-823a-4873-8c61-cfda278f8a83");
    }

    public class PattComponentAttributes : GH_ComponentAttributes
    {
        private PattComponent _owner;

        public PattComponentAttributes(PattComponent owner) : base(owner)
        {
            _owner = owner;
        }

        protected override void Layout()
        {
            base.Layout();
            Bounds = new RectangleF(Bounds.X, Bounds.Y, Bounds.Width, Bounds.Height + 25);
        }

        protected override void Render(GH_Canvas canvas, Graphics graphics, GH_CanvasChannel channel)
        {
            base.Render(canvas, graphics, channel);

            if (channel == GH_CanvasChannel.Objects)
            {
                float buttonHeight = 20;
                float buttonWidth = Bounds.Width - 10;
                float buttonX = Bounds.X + 5;
                float buttonY = Bounds.Bottom - buttonHeight - 5;

                RectangleF buttonBounds = new RectangleF(buttonX, buttonY, buttonWidth, buttonHeight);

                Brush buttonColor = _owner.InfiniteMode ? Brushes.Orange : Brushes.Gray;

                using (GraphicsPath path = new GraphicsPath())
                {
                    float cornerRadius = 5f;
                    path.AddArc(buttonBounds.X, buttonBounds.Y, cornerRadius * 2, cornerRadius * 2, 180, 90);
                    path.AddArc(buttonBounds.Right - cornerRadius * 2, buttonBounds.Y, cornerRadius * 2, cornerRadius * 2, 270, 90);
                    path.AddArc(buttonBounds.Right - cornerRadius * 2, buttonBounds.Bottom - cornerRadius * 2, cornerRadius * 2, cornerRadius * 2, 0, 90);
                    path.AddArc(buttonBounds.X, buttonBounds.Bottom - cornerRadius * 2, cornerRadius * 2, cornerRadius * 2, 90, 90);
                    path.CloseFigure();

                    graphics.FillPath(buttonColor, path);
                    graphics.DrawPath(Pens.Black, path);
                }

                using (Font font = new Font(GH_FontServer.Standard.FontFamily, 10f, FontStyle.Bold))
                {
                    StringFormat fmt = new StringFormat();
                    fmt.Alignment = StringAlignment.Center;
                    fmt.LineAlignment = StringAlignment.Center;
                    graphics.DrawString("∞", font, Brushes.White, buttonBounds, fmt);
                }
            }
        }

        public override GH_ObjectResponse RespondToMouseDown(GH_Canvas sender, GH_CanvasMouseEvent e)
        {
            float buttonHeight = 20;
            float buttonWidth = Bounds.Width - 10;
            float buttonX = Bounds.X + 5;
            float buttonY = Bounds.Bottom - buttonHeight - 5;

            RectangleF buttonBounds = new RectangleF(buttonX, buttonY, buttonWidth, buttonHeight);

            if (e.Button == MouseButtons.Left && buttonBounds.Contains(e.CanvasLocation))
            {
                _owner.InfiniteMode = !_owner.InfiniteMode;
                _owner.ExpireSolution(true);
                return GH_ObjectResponse.Handled;
            }

            return base.RespondToMouseDown(sender, e);
        }
    }
}