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
            IRRFResolver rrfResolver, IGenomeWriter genomeWriter, string[] args)
        {
            _genomeReader = genomeReader;
            _vcfReader = vcfReader;
            _rrfProcessor = rrfProcessor;
            _rrfResolver = rrfResolver;
            _genomeWriter = genomeWriter;
        }

        public void Launch()
        {
            
        }
    }
}
