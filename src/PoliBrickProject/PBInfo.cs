using System;
using System.Drawing;
using System.IO;
using Grasshopper.Kernel;

namespace PoliBrick.Util
{
  public class PBInfo : GH_AssemblyInfo
  {
    public override string Name => "PoliBrick";

        //Return a 24x24 pixel bitmap to represent this GHA library.
        public override Bitmap Icon => null;
        public static Bitmap ConvertByteArrayToBitmap(byte[] imageData)
        {
            using (MemoryStream ms = new MemoryStream(imageData))
            {
                return new Bitmap(ms);
            }
        }
        

        //Return a short string describing the purpose of this GHA library.
        public override string Description => "";

    public override Guid Id => new Guid("0cb54b33-ccfe-496a-a71a-d0cbd2d6beb5");

    //Return a string identifying you or your company.
    public override string AuthorName => "Mohammad_Pourfouladi";

    //Return a string representing your preferred contact details.
    public override string AuthorContact => "mohammad.pourfouladi@polimi.it";

    //Return a string representing the version.  This returns the same version as the assembly.
    public override string AssemblyVersion => "2.0.0";
  }
}