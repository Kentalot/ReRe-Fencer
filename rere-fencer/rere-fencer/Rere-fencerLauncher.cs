using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using NDesk.Options;
using rere_fencer.Input;
using rere_fencer.Processors;

namespace rere_fencer
{
    public class Rere_fencerLauncher
    {
        private IGenomeReader _genomeReader;
        private IVCFReader _vcfReader;
        private IRRFProcessor _rrfProcessor;
        private IRRFResolver _rrfResolver;
        private IGenomeWriter _genomeWriter;

        public Rere_fencerLauncher(IGenomeReader genomeReader, IVCFReader vcfReader, IRRFProcessor rrfProcessor,
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
            foreach (var contig in _genomeReader.Contigs.Values)
                try
                {
                    Console.WriteLine("Chr " + contig.Name + " from " + Program.GenomeRange.Start + " to " + Program.GenomeRange.End + "=" +
                                      new string(contig.GetSequence(Program.GenomeRange.Start, Program.GenomeRange.End).ToArray()));
                }
            catch (Exception E){ Console.Error.WriteLine(E.Message + "\n" + E.StackTrace);}
            var contig2 = new Contig("chr17");
            Console.WriteLine("Chr " + contig2.Name + " from " + Program.GenomeRange.Start + " to " + Program.GenomeRange.End + "=" +
                                      new string(_genomeReader.Contigs[contig2].GetSequence(Program.GenomeRange.Start, Program.GenomeRange.End).ToArray()));
        }
    }
}
