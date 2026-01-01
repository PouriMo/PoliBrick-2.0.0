using System;
using System.Collections.Generic;
using System.Drawing;
using System.Drawing.Drawing2D;
using System.Linq;
using System.Windows.Forms;
using Grasshopper.GUI;
using Grasshopper.GUI.Canvas;
using Grasshopper.Kernel;
using Grasshopper.Kernel.Attributes;

namespace PoliBrickUtils
{
    // ==============================================================================
    // 1. Data Structure
    // ==============================================================================
    public class SliderRow
    {
        public double Min { get; set; } = 0.0;
        public double Max { get; set; } = 1.0;
        public double Current { get; set; } = 0.25;

        public void ClampCurrent()
        {
            Current = Math.Max(Min, Math.Min(Max, Current));
        }

        public float NormalizedT
        {
            get
            {
                if (Math.Abs(Max - Min) < 1e-9) return 0f;
                return (float)((Current - Min) / (Max - Min));
            }
        }
    }

    // ==============================================================================
    // 2. The Component
    // ==============================================================================
    public class MultiSliderComponent : GH_Component
    {
        public List<SliderRow> Rows { get; set; } = new List<SliderRow>();

        public MultiSliderComponent()
          : base("Multi-Slider Panel", "MSlider",
              "Dynamic vertical sliders list.",
              "PoliBrick", "Utilities")
        {
            Rows.Add(new SliderRow() { Current = 0.25 });
            Rows.Add(new SliderRow() { Current = 0.50 });
            Rows.Add(new SliderRow() { Current = 0.75 });
        }

        protected override void RegisterInputParams(GH_InputParamManager pManager) { }

        protected override void RegisterOutputParams(GH_OutputParamManager pManager)
        {
            pManager.AddTextParameter("Values", "V", "Concatenated values", GH_ParamAccess.item);
        }

        protected override void SolveInstance(IGH_DataAccess DA)
        {
            string result = string.Join(",", Rows.Select(r => r.Current.ToString("0.###")));
            DA.SetData(0, result);
        }

        public override void CreateAttributes()
        {
            m_attributes = new MultiSliderAttributes(this);
        }

        public override bool Write(GH_IO.Serialization.GH_IWriter writer)
        {
            writer.SetInt32("RowCount", Rows.Count);
            for (int i = 0; i < Rows.Count; i++)
            {
                writer.SetDouble($"R{i}_Min", Rows[i].Min);
                writer.SetDouble($"R{i}_Max", Rows[i].Max);
                writer.SetDouble($"R{i}_Cur", Rows[i].Current);
            }
            return base.Write(writer);
        }

        public override bool Read(GH_IO.Serialization.GH_IReader reader)
        {
            Rows.Clear();
            if (reader.ItemExists("RowCount"))
            {
                int count = reader.GetInt32("RowCount");
                for (int i = 0; i < count; i++)
                {
                    SliderRow row = new SliderRow();
                    if (reader.ItemExists($"R{i}_Min")) row.Min = reader.GetDouble($"R{i}_Min");
                    if (reader.ItemExists($"R{i}_Max")) row.Max = reader.GetDouble($"R{i}_Max");
                    if (reader.ItemExists($"R{i}_Cur")) row.Current = reader.GetDouble($"R{i}_Cur");
                    Rows.Add(row);
                }
            }
            return base.Read(reader);
        }

        protected override System.Drawing.Bitmap Icon => PoliBrick.Properties.Resources.sliders;
        public override Guid ComponentGuid => new Guid("B8A9D2F1-4E6C-4A2D-9B1F-5C8E3A7D6F09");
    }

    // ==============================================================================
    // 3. Custom Attributes (The UI)
    // ==============================================================================
    public class MultiSliderAttributes : GH_ComponentAttributes
    {
        private MultiSliderComponent _owner;

        // Layout Constants
        private const float TopMargin = 12f;
        private const float RowHeight = 24f;
        private const float DelBtnWidth = 15f;
        private const float TextBoxWidth = 40f;
        private const float AddBarHeight = 22f;
        private const float LeftMargin = 6f;
        private const float RightMargin = 8f;
        private const float SliderMinWidth = 100f;

        // Grip Constants
        private const float GripWidth = 16f;  // Wider for pill shape
        private const float GripHeight = 10f;
        private const float GripCornerRad = 4f;

        private int _dragRowIndex = -1;
        private bool _isDragging = false;

        // Hit Caches
        private List<RectangleF> _trackRects = new List<RectangleF>();
        private List<RectangleF> _minBoxRects = new List<RectangleF>();
        private List<RectangleF> _maxBoxRects = new List<RectangleF>();
        private List<RectangleF> _delBtnRects = new List<RectangleF>();
        private RectangleF _addBarRect;

        public MultiSliderAttributes(MultiSliderComponent owner) : base(owner)
        {
            _owner = owner;
        }

        protected override void Layout()
        {
            float width = LeftMargin + DelBtnWidth + TextBoxWidth + SliderMinWidth + TextBoxWidth + RightMargin;
            float height = TopMargin + (_owner.Rows.Count * RowHeight) + AddBarHeight + 8f;
            Bounds = new RectangleF(Pivot.X, Pivot.Y, width, height);

            if (_owner.Params.Output.Count > 0)
            {
                var outputParam = _owner.Params.Output[0];
                PointF outPivot = new PointF(Bounds.Right, Bounds.Y + TopMargin + (RowHeight / 2));
                outputParam.Attributes.Pivot = outPivot;
                outputParam.Attributes.Bounds = new RectangleF(outPivot.X - 5, outPivot.Y - 5, 10, 10);
            }
        }

        protected override void Render(GH_Canvas canvas, Graphics graphics, GH_CanvasChannel channel)
        {
            if (channel == GH_CanvasChannel.Objects)
            {
                // A. Output Grip
                if (_owner.Params.Output.Count > 0)
                {
                    PointF outPivot = _owner.Params.Output[0].Attributes.Pivot;
                    float size = 9f;
                    RectangleF gripRect = new RectangleF(Bounds.Right - (size / 2), outPivot.Y - (size / 2), size, size);

                    graphics.FillEllipse(Brushes.White, gripRect);
                    using (Pen thickPen = new Pen(Color.Black, 2f))
                    {
                        graphics.DrawEllipse(thickPen, gripRect);
                    }
                }

                // B. Capsule
                GH_Capsule capsule = GH_Capsule.CreateCapsule(Bounds, GH_Palette.Normal);
                capsule.Render(graphics, Selected, Owner.Locked, false);
                capsule.Dispose();

                // C. Content
                _trackRects.Clear(); _minBoxRects.Clear(); _maxBoxRects.Clear(); _delBtnRects.Clear();

                Font fontText = new Font("Arial", 8f, FontStyle.Regular);
                Font fontVal = new Font("Arial", 6f, FontStyle.Regular);

                Pen borderPen = new Pen(Color.FromArgb(100, 100, 100));
                Brush boxBgBrush = Brushes.White;
                Brush trackBgBrush = new SolidBrush(Color.FromArgb(220, 220, 220));

                Brush gripBrush = _isDragging ? Brushes.Orange : Brushes.White;
                Pen gripPen = new Pen(Color.LightGray); // Light Gray stroke

                float currentY = Bounds.Y + TopMargin;

                for (int i = 0; i < _owner.Rows.Count; i++)
                {
                    SliderRow row = _owner.Rows[i];
                    float rowTop = currentY + (i * RowHeight);
                    float rowCenterY = rowTop + (RowHeight / 2f);
                    float x = Bounds.X + LeftMargin;

                    // 1. Delete (x)
                    RectangleF delRect = new RectangleF(x, rowTop, DelBtnWidth, RowHeight);
                    _delBtnRects.Add(delRect);
                    graphics.DrawString("×", fontText, Brushes.Gray, delRect, GH_TextRenderingConstants.CenterCenter);
                    x += DelBtnWidth;

                    // 2. Min Box (Rounded)
                    RectangleF minRect = new RectangleF(x, rowTop + 4, TextBoxWidth, RowHeight - 8);
                    _minBoxRects.Add(minRect);
                    using (GraphicsPath minPath = RoundedRect(minRect, 3f))
                    {
                        graphics.FillPath(boxBgBrush, minPath);
                        graphics.DrawPath(borderPen, minPath);
                    }
                    graphics.DrawString(row.Min.ToString("0.##"), fontText, Brushes.Black, minRect, GH_TextRenderingConstants.CenterCenter);
                    x += TextBoxWidth + 6f;

                    // 3. Slider Track
                    float remainingW = Bounds.Right - RightMargin - TextBoxWidth - 6f - x;

                    RectangleF fullTrackArea = new RectangleF(x, rowTop, remainingW, RowHeight);
                    _trackRects.Add(fullTrackArea);

                    // Draw Visual Track (Rounded)
                    float trackH = 4f;
                    RectangleF visualTrack = new RectangleF(x, rowCenterY - trackH / 2, remainingW, trackH);
                    using (GraphicsPath trackPath = RoundedRect(visualTrack, 2f))
                    {
                        graphics.FillPath(trackBgBrush, trackPath);
                    }

                    // Ticks
                    int tickCount = 5;
                    for (int t = 0; t <= tickCount; t++)
                    {
                        float tx = x + (remainingW * (t / (float)tickCount));
                        graphics.DrawLine(Pens.LightGray, tx, rowCenterY - 3, tx, rowCenterY + 3);
                    }

                    // 4. Max Box (Rounded)
                    float maxBoxX = x + remainingW + 6f;
                    RectangleF maxRect = new RectangleF(maxBoxX, rowTop + 4, TextBoxWidth, RowHeight - 8);
                    _maxBoxRects.Add(maxRect);
                    using (GraphicsPath maxPath = RoundedRect(maxRect, 3f))
                    {
                        graphics.FillPath(boxBgBrush, maxPath);
                        graphics.DrawPath(borderPen, maxPath);
                    }
                    graphics.DrawString(row.Max.ToString("0.##"), fontText, Brushes.Black, maxRect, GH_TextRenderingConstants.CenterCenter);

                    // 5. Grip & Value (No Overlap Logic)
                    // The "Travel Distance" is the full width MINUS the width of the grip.
                    // This keeps the grip entirely inside the track area.
                    float travelW = remainingW - GripWidth;

                    // The left edge of the grip
                    float gripLeft = x + (travelW * row.NormalizedT);

                    // The center of the grip (for text alignment)
                    float gripCenterX = gripLeft + (GripWidth / 2f);

                    // Draw Grip
                    RectangleF gripRect = new RectangleF(gripLeft, rowCenterY - GripHeight / 2, GripWidth, GripHeight);
                    using (GraphicsPath gripPath = RoundedRect(gripRect, GripCornerRad))
                    {
                        graphics.FillPath(gripBrush, gripPath);
                        graphics.DrawPath(gripPen, gripPath);
                    }

                    // Floating Value
                    string valStr = row.Current.ToString("0.00");
                    SizeF valSize = graphics.MeasureString(valStr, fontVal);

                    // Align text to grip center
                    float textX = gripCenterX - (valSize.Width / 2);

                    // Clamp text so it doesn't bleed into the text boxes
                    // Min X allowed: start of track
                    // Max X allowed: end of track - text width
                    textX = Math.Max(x, Math.Min(x + remainingW - valSize.Width, textX));

                    graphics.DrawString(valStr, fontVal, Brushes.Black, textX, rowCenterY - 18);
                }

                // Add Bar
                float addY = currentY + (_owner.Rows.Count * RowHeight) + 2f;
                _addBarRect = new RectangleF(Bounds.X + LeftMargin, addY, Bounds.Width - LeftMargin - RightMargin, AddBarHeight);

                using (GraphicsPath path = RoundedRect(_addBarRect, 3f))
                {
                    using (SolidBrush fillB = new SolidBrush(Color.FromArgb(245, 245, 245)))
                        graphics.FillPath(fillB, path);
                    using (Pen borderP = new Pen(Color.DarkGray, 1f))
                        graphics.DrawPath(borderP, path);
                }

                PointF center = new PointF(_addBarRect.X + _addBarRect.Width / 2, _addBarRect.Y + _addBarRect.Height / 2);
                graphics.DrawLine(Pens.Gray, center.X - 4, center.Y, center.X + 4, center.Y);
                graphics.DrawLine(Pens.Gray, center.X, center.Y - 4, center.X, center.Y + 4);

                fontText.Dispose();
                fontVal.Dispose();
                borderPen.Dispose();
                trackBgBrush.Dispose();
                gripPen.Dispose();
            }
        }

        private GraphicsPath RoundedRect(RectangleF r, float rad)
        {
            GraphicsPath path = new GraphicsPath();
            float d = rad * 2;
            d = Math.Min(d, Math.Min(r.Width, r.Height));
            path.AddArc(r.X, r.Y, d, d, 180, 90);
            path.AddArc(r.Right - d, r.Y, d, d, 270, 90);
            path.AddArc(r.Right - d, r.Bottom - d, d, d, 0, 90);
            path.AddArc(r.X, r.Bottom - d, d, d, 90, 90);
            path.CloseFigure();
            return path;
        }

        public override GH_ObjectResponse RespondToMouseDown(GH_Canvas sender, GH_CanvasMouseEvent e)
        {
            if (e.Button == MouseButtons.Left)
            {
                if (_addBarRect.Contains(e.CanvasLocation))
                {
                    double lastMax = _owner.Rows.Count > 0 ? _owner.Rows.Last().Max : 1.0;
                    _owner.Rows.Add(new SliderRow() { Min = 0, Max = lastMax, Current = lastMax / 2.0 });
                    _owner.ExpireSolution(true);
                    return GH_ObjectResponse.Handled;
                }

                for (int i = 0; i < _owner.Rows.Count; i++)
                {
                    if (_delBtnRects[i].Contains(e.CanvasLocation))
                    {
                        _owner.Rows.RemoveAt(i);
                        _owner.ExpireSolution(true);
                        return GH_ObjectResponse.Handled;
                    }

                    // Hit Test Logic
                    SliderRow row = _owner.Rows[i];
                    RectangleF track = _trackRects[i];

                    float travelW = track.Width - GripWidth;
                    float gripLeft = track.X + (travelW * row.NormalizedT);
                    float rowCenterY = track.Top + (RowHeight / 2f);

                    float padding = 2f;
                    RectangleF gripHitRect = new RectangleF(
                        gripLeft - padding,
                        rowCenterY - GripHeight / 2 - padding,
                        GripWidth + padding * 2,
                        GripHeight + padding * 2);

                    if (gripHitRect.Contains(e.CanvasLocation))
                    {
                        _dragRowIndex = i;
                        _isDragging = true;
                        sender.Refresh();
                        return GH_ObjectResponse.Capture;
                    }
                }
            }
            return base.RespondToMouseDown(sender, e);
        }

        public override GH_ObjectResponse RespondToMouseMove(GH_Canvas sender, GH_CanvasMouseEvent e)
        {
            if (_isDragging && _dragRowIndex >= 0 && e.Button == MouseButtons.Left)
            {
                UpdateSliderValue(_dragRowIndex, e.CanvasX);
                _owner.ExpireSolution(true);
                return GH_ObjectResponse.Handled;
            }
            return base.RespondToMouseMove(sender, e);
        }

        public override GH_ObjectResponse RespondToMouseUp(GH_Canvas sender, GH_CanvasMouseEvent e)
        {
            if (_isDragging)
            {
                _isDragging = false;
                _dragRowIndex = -1;
                sender.Refresh();
                _owner.ExpireSolution(true);
                return GH_ObjectResponse.Release;
            }
            return base.RespondToMouseUp(sender, e);
        }

        public override GH_ObjectResponse RespondToMouseDoubleClick(GH_Canvas sender, GH_CanvasMouseEvent e)
        {
            if (e.Button == MouseButtons.Left)
            {
                for (int i = 0; i < _owner.Rows.Count; i++)
                {
                    SliderRow row = _owner.Rows[i];
                    bool changed = false;

                    if (_minBoxRects[i].Contains(e.CanvasLocation))
                    {
                        if (BetterNumberInputDialog.Show(row.Min, Cursor.Position, out double newVal) == DialogResult.OK)
                        {
                            row.Min = newVal; row.ClampCurrent(); changed = true;
                        }
                    }
                    else if (_maxBoxRects[i].Contains(e.CanvasLocation))
                    {
                        if (BetterNumberInputDialog.Show(row.Max, Cursor.Position, out double newVal) == DialogResult.OK)
                        {
                            row.Max = newVal; row.ClampCurrent(); changed = true;
                        }
                    }
                    // Allow double-clicking on the general track area to type value
                    else if (_trackRects[i].Contains(e.CanvasLocation))
                    {
                        if (BetterNumberInputDialog.Show(row.Current, Cursor.Position, out double newVal) == DialogResult.OK)
                        {
                            row.Current = newVal; row.ClampCurrent(); changed = true;
                        }
                    }

                    if (changed) { _owner.ExpireSolution(true); return GH_ObjectResponse.Handled; }
                }
            }
            return base.RespondToMouseDoubleClick(sender, e);
        }

        private void UpdateSliderValue(int rowIndex, float mouseX)
        {
            if (rowIndex < 0 || rowIndex >= _owner.Rows.Count) return;
            SliderRow row = _owner.Rows[rowIndex];
            RectangleF track = _trackRects[rowIndex];

            // MATH UPDATE:
            // The value is 0.0 when mouse is at (track.X + GripWidth/2)
            // The value is 1.0 when mouse is at (track.Right - GripWidth/2)

            float startX = track.X + (GripWidth / 2f);
            float endX = track.Right - (GripWidth / 2f);
            float usableWidth = endX - startX;

            float normalizedT = (mouseX - startX) / usableWidth;
            normalizedT = Math.Max(0f, Math.Min(1f, normalizedT));

            row.Current = row.Min + (normalizedT * (row.Max - row.Min));
        }
    }

    public class BetterNumberInputDialog : Form
    {
        private TextBox _textBox;
        private Button _btnOk;
        public double Value { get; private set; }

        private BetterNumberInputDialog(double initialValue, Point cursorLocation)
        {
            this.Text = "Set Value";
            this.Size = new Size(180, 90);
            this.FormBorderStyle = FormBorderStyle.None;
            this.StartPosition = FormStartPosition.Manual;
            this.Location = new Point(cursorLocation.X - 90, cursorLocation.Y - 45);
            this.Paint += (s, e) => { e.Graphics.DrawRectangle(Pens.Gray, 0, 0, Width - 1, Height - 1); };
            this.BackColor = Color.WhiteSmoke;
            _textBox = new TextBox() { Left = 10, Top = 10, Width = 160, Text = initialValue.ToString() };
            _btnOk = new Button() { Text = "Apply", Left = 10, Width = 160, Top = 40, Height = 25, DialogResult = DialogResult.OK, BackColor = Color.White };
            this.Controls.Add(_textBox);
            this.Controls.Add(_btnOk);
            this.AcceptButton = _btnOk;
            _btnOk.Click += (s, e) => {
                if (double.TryParse(_textBox.Text, out double val)) { Value = val; Close(); }
                else { MessageBox.Show("Invalid number"); this.DialogResult = DialogResult.None; }
            };
        }

        public static DialogResult Show(double initialValue, Point cursorLoc, out double result)
        {
            using (var form = new BetterNumberInputDialog(initialValue, cursorLoc))
            {
                var dr = form.ShowDialog();
                result = form.Value;
                return dr;
            }
        }
    }
}