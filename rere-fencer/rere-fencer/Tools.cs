using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using AppUtils;
using AppUtils.Misc;

namespace rere_fencer.Processors
{
    public static class Tools
    {
        private static string toolspath = Path.Combine(AppDomain.CurrentDomain.BaseDirectory, "Tools");

        public static string WrapDoubleQuotes(this string theStringToWrap)
        { return string.Format("{0}{1}{0}", '"', theStringToWrap); }

        private static string GetLinuxToolPath(string folderpath, string toolname, 
            string linuxextension = null, string windowsextension = null)
        {
            string appname;
            if (MiscUtils.RunningMono)
            {
                appname = linuxextension == null ? toolname : toolname + linuxextension;
            }
            else
            {
                appname = windowsextension == null
                           ? Path.Combine(folderpath, toolname)
                           : Path.Combine(folderpath, toolname + windowsextension);
            }
            return appname.WrapDoubleQuotes();
        }

        public static string Samtools
        {
            get { return GetLinuxToolPath(Path.Combine(toolspath, "samtools-windows"), "samtools", windowsextension: ".exe"); }
        }

        public static string Tabix
        {
            get { return GetLinuxToolPath(Path.Combine(toolspath, "samtools-windows"), "tabix", windowsextension: ".exe"); }
        }

        public static string Bgzip
        {
            get { return GetLinuxToolPath(Path.Combine(toolspath, "samtools-windows"), "bgzip", windowsextension: ".exe"); }
        }

        public static string Bcftools
        {
            get { return GetLinuxToolPath(Path.Combine(toolspath, "samtools-windows"), "bcftools", windowsextension: ".exe"); }
        }

        public static string Razip
        {
            get { return GetLinuxToolPath(Path.Combine(toolspath, "samtools-windows"), "razip", windowsextension: ".exe"); }
        }
    }
}
