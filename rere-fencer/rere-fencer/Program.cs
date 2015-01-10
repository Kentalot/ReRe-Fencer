using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using NDesk.Options;
using rere_fencer.Input;
using rere_fencer.Processors;

namespace rere_fencer
{
    class Program
    {
        private static readonly OptionSet _optionSet = new OptionSet();
        public static FileInfo GenomeFilePath { get; private set; }

        static void Main(string[] args)
        {
            InitializeOptionSet();
            new Rere_fencerLauncher(new GenomeReader(), new VCFReader(), new RRFProcessor(), new RRFResolver(), new GenomeWriter(), args).Launch();
        }

        private static void InitializeOptionSet()
        {
            _optionSet.Add("g|genomePath=", "Path to the Genome input file. Must be in one fasta file",
                s => GenomeFilePath = UpdateFileInfo(s));
        }

        private static FileInfo UpdateFileInfo(string filePath)
        {
            if (string.IsNullOrWhiteSpace(filePath)) throw new ArgumentException("Invalid Blank filePath input!");
            var file = new FileInfo(filePath);
            var invalid = !file.Exists || file.Length < 100;
            if (invalid) throw new ArgumentException("File not found or file too small to be input!");
            return file;
        }
    }
}
