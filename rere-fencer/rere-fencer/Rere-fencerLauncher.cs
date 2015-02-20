using System;
using System.Collections.Generic;
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
        private IVcfReader _vcfReader;
        private IRRFProcessor _rrfProcessor;
        private IRRFResolver _rrfResolver;
        private IGenomeWriter _genomeWriter;

        public ReRe_fencerLauncher(IGenomeReader genomeReader, IVcfReader vcfReader, IRRFProcessor rrfProcessor,
            IRRFResolver rrfResolver, IGenomeWriter genomeWriter)
        {
            _genomeReader = genomeReader;
            _vcfReader = vcfReader;
            _rrfProcessor = rrfProcessor;
            _rrfResolver = rrfResolver;
            _genomeWriter = genomeWriter;
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
            var position = _genomeReader.IsZeroBasedCoordinates ? 0U : 1U;
            foreach (var variantInfo in _vcfReader.ReadVcfVariants().Select(v => v.CreateVariantInfo())
                .Where(v => v.SampleInfo.IsHom() || v.SampleInfo.IsHemi()))
            {
                var end = variantInfo.Variant.Position;
            }
            
        }
    }
}
