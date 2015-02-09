using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using BioinformaticUtils.DataStructures.Contigs;
using BioinformaticUtils.GenomeTools;
using BioinformaticUtils.Vcf;
using NDesk.Options;
using rere_fencer.Processors;

namespace rere_fencer
{
    class Program
    {
        private static readonly OptionSet _optionSet = new OptionSet();
        internal static FileInfo GenomeFilePath { get; private set; }
        internal static FileInfo VcfFilePath { get; private set; }
        internal static FileInfo SvVcfFilePath { get; private set; }
        internal static IContigInterval GenomeRange { get; private set; }

        static void Main(string[] args)
        {
            InitializeOptionSet();
            _optionSet.Parse(args);
            ValidateOptions();
            var sw = new Stopwatch();
            sw.Start();
            using (var genomeReader = new TwoBitGenomeReader(GenomeFilePath.FullName))
            {
                sw.Stop();
                Console.WriteLine(sw.Elapsed);
                sw.Reset();
                sw.Start();
                new ReRe_fencerLauncher(genomeReader, new VcfReader(), new RRFProcessor(), new RRFResolver(), new GenomeWriter()).Launch();
                sw.Stop();
                Console.WriteLine(sw.Elapsed);
            }
        }

        private static void ValidateOptions()
        {
            if (GenomeFilePath == null) ShowHelp();
        }

        private static void ShowHelp()
        {
            _optionSet.WriteOptionDescriptions(Console.Out);
        }

        private static void InitializeOptionSet()
        {
            _optionSet.Add("g|genomePath=", "Path to the Genome input file. Must be in one fasta file or 2bit file",
                s => GenomeFilePath = UpdateFileInfo(s));
            _optionSet.Add("v|vcfPath=", "Path to the VCF file that contains the variants for the genome specified.",
                s => VcfFilePath = UpdateFileInfo(s));
            _optionSet.Add("s|svVCFPath=",
                "Path to the SV.VCF file that contains the Structural Variants (and possibly Copy Number Variants) for the genome specified.",
                s => SvVcfFilePath = UpdateFileInfo(s));
            _optionSet.Add("q|query=",
                "range in 1-based of coordinates you want to print out. Separate with a dash aka -, like this: 28-724",
                s =>
                {
                    var split = s.Split('-').Take(2).Select(uint.Parse).ToArray();
                    GenomeRange = new ContigInterval(split[0], split[1]);
                    Console.WriteLine("Range = {0}-{1}", GenomeRange.Start, GenomeRange.End);
                });
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
