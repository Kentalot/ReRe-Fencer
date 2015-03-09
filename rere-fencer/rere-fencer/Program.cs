using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using BioinformaticUtils.DataStructures.Contigs;
using BioinformaticUtils.DataStructures.Genomes;
using BioinformaticUtils.DataStructures.Nucleotides;
using BioinformaticUtils.Vcf;
using NDesk.Options;
using rere_fencer.Processors;

namespace rere_fencer
{
    class Program
    {
        private static readonly OptionSet OptionSet = new OptionSet();
        internal static FileInfo GenomeFilePath { get; private set; }
        internal static FileInfo VcfFilePath { get; private set; }
        internal static FileInfo SvVcfFilePath { get; private set; }
        internal static FileInfo OutputFilePath { get; private set; }
        internal static IContigInterval GenomeRange { get; private set; }
        internal static string QueryContig { get; private set; }

        static void Main(string[] args)
        {
            InitializeOptionSet();
            OptionSet.Parse(args);
            ValidateOptions();
            var sw = new Stopwatch();
            sw.Start();
            if (GenomeRange == null)
            {
                using (var genomeReader = GenomeFilePath == null 
                    ? (TwoBitGenomeReader) TwoBitGenomeReader.Hg19TwoBitReader : new TwoBitGenomeReader(GenomeFilePath.FullName))
                {
                    using (var vcfReader = new TabixVcfReader(VcfFilePath))
                    {
                        sw.Stop();
                        Console.WriteLine(sw.Elapsed);
                        sw.Reset();
                        sw.Start();
                        new ReRe_fencerLauncher(genomeReader, vcfReader, new SmallVariantsRRFProcessor(),
                            new RRFResolver(),
                            OutputFilePath).Launch();
                    }
                    sw.Stop();
                    Console.WriteLine(sw.Elapsed);
                }
            }
            else
            {
                using (var genomeReader = new TwoBitGenomeReader(GenomeFilePath.FullName))
                {
                    sw.Stop();
                    Console.WriteLine(sw.Elapsed);
                    sw.Reset();
                    sw.Start();
                    foreach (var contig in QueryContig == null ? genomeReader.Contigs.Values : new [] {genomeReader.Contigs[QueryContig]})
                        try
                        {
                            Console.WriteLine("Chr " + contig.Name + " from " + GenomeRange.Start + " to " +
                                              GenomeRange.End + "=" +
                                              new NucleotideString(contig.GetSequence(GenomeRange.Start,
                                                  GenomeRange.End)));
                        }
                        catch (Exception E)
                        {
                            Console.Error.WriteLine(E.Message + "\n" + E.StackTrace);
                        }
                    if (QueryContig != null)
                    {
                        var contig2 = "chr17";
                        Console.WriteLine("Chr " + contig2 + " from " + GenomeRange.Start + " to " +
                                          GenomeRange.End + "=" +
                                          new NucleotideString(
                                              genomeReader.Contigs[contig2].GetSequence(GenomeRange.Start,
                                                  GenomeRange.End)));
                    }
                }
                sw.Stop();
                Console.WriteLine(sw.Elapsed);
            }
        
        }

        private static void ValidateOptions()
        {   
            if (GenomeRange != null) return;
            if (OutputFilePath == null || VcfFilePath == null)
            {
                ShowHelp();
                Environment.Exit(0);
            }
        }

        private static void ShowHelp()
        {
            OptionSet.WriteOptionDescriptions(Console.Out);
        }

        private static void InitializeOptionSet()
        {
            OptionSet.Add("g|genomePath=", "Path to the Genome input file. Must be in 2bit file format for now. (Required!)",
                s => GenomeFilePath = UpdateFileInfo(s));
            OptionSet.Add("v|vcfPath=", "Path to the VCF file that contains the variants for the genome specified.",
                s => VcfFilePath = UpdateFileInfo(s));
            OptionSet.Add("s|svVCFPath=",
                "Path to the SV.VCF file that contains the Structural Variants (and possibly Copy Number Variants) for the genome specified.",
                s => SvVcfFilePath = UpdateFileInfo(s));
            OptionSet.Add("o|outputFilePath=", "Path to the output file you want. Only outputs in fasta format for now. e.g. genome.fa (Required!)",
                s => OutputFilePath = UpdateFileInfo(s, false));
            OptionSet.Add("q|query=",
                "range you want to print out. Separate with a dash aka -. Optional would be to include the contig name separated by colon aka :, like this: chr1:28-724",
                s =>
                {
                    var split = s.Split(':').Take(2).ToArray();
                    if (split.Length > 1)
                    {
                        QueryContig = split[0];
                        split = split[1].Split('-').Take(2).ToArray();
                    }
                    GenomeRange = new ContigInterval(uint.Parse(split[0]), uint.Parse(split[1]));
                    Console.WriteLine("Range = {0}-{1}", GenomeRange.Start, GenomeRange.End);
                });
        }

        private static FileInfo UpdateFileInfo(string filePath, bool testForExistence = true)
        {
            if (string.IsNullOrWhiteSpace(filePath)) throw new ArgumentException("Invalid Blank filePath input!");
            var file = new FileInfo(filePath);
            if (!testForExistence) return file;
            var invalid = !file.Exists || file.Length < 100;
            if (invalid) throw new ArgumentException("File not found or file too small to be input!");
            return file;
        }
    }
}
