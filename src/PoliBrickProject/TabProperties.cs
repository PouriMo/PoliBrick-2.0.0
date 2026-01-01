using System.Drawing;
using Grasshopper.Kernel;
using System.IO;

namespace PoliBrick.Util
{
    public class TabProperties : GH_AssemblyPriority
    {
        public static Bitmap ConvertByteArrayToBitmap(byte[] imageData)
        {
            using (MemoryStream ms = new MemoryStream(imageData))
            {
                return new Bitmap(ms);
            }
        }

        public override GH_LoadingInstruction PriorityLoad()
        {
            var server = Grasshopper.Instances.ComponentServer;
            server.AddCategoryShortName("PoliBrick","PB");
            server.AddCategorySymbolName("PoliBrick", 'P');
            server.AddCategoryIcon("PoliBrick", ConvertByteArrayToBitmap(PoliBrick.Properties.Resources.PBIcon));

            return GH_LoadingInstruction.Proceed;
        }
    }
}

