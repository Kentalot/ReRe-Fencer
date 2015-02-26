using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using BioinformaticUtils.DataStructures.Contigs;
using BioinformaticUtils.DataStructures.Nucleotides;
using BioinformaticUtils.GenomeTools;
using BioinformaticUtils.Vcf;
using NDesk.Options;
using rere_fencer.Processors;

namespace rere_fencer
{
    public class ReRe_fencerLauncher
    {
        private IGenomeReader _genomeReader;
        private TabixVcfReader _vcfReader;
        private IRRFProcessor _rrfProcessor;
        private IRRFResolver _rrfResolver;
        private readonly FileInfo _outputFile;

        public ReRe_fencerLauncher(IGenomeReader genomeReader, TabixVcfReader vcfReader, IRRFProcessor rrfProcessor,
            IRRFResolver rrfResolver, FileInfo outputFile)
        {
            _genomeReader = genomeReader;
            _vcfReader = vcfReader;
            _rrfProcessor = rrfProcessor;
            _rrfResolver = rrfResolver;
            _outputFile = outputFile;
        }

        public void Launch()
        {
            /*foreach (var contig in _genomeReader.Contigs.Values)
                try
                {
                    Console.WriteLine("Chr " + contig.Name + " from " + Program.GenomeRange.Start + " to " + Program.GenomeRange.End + "=" +
                                      new NucleotideString(contig.GetSequence(Program.GenomeRange.Start, Program.GenomeRange.End)));
                }
            catch (Exception E){ Console.Error.WriteLine(E.Message + "\n" + E.StackTrace);}
            var contig2 = new Contig("chr17");
            Console.WriteLine("Chr " + contig2.Name + " from " + Program.GenomeRange.Start + " to " + Program.GenomeRange.End + "=" +
                                      new NucleotideString(_genomeReader.Contigs[contig2].GetSequence(Program.GenomeRange.Start, Program.GenomeRange.End)));*/
            
            new FastaWriter(GenerateReReFencedContigs(), _outputFile).WriteAllContigs();
        }

        private IEnumerable<IGenomeContig> GenerateReReFencedContigs()
        {
            return _genomeReader.Contigs.Select(s => s.Value).Select(genomeContig =>
                    new GenomeContig(genomeContig, _rrfProcessor.Process(genomeContig, _vcfReader.GetVariantsForChromosome(genomeContig))));
        }
    }
}
