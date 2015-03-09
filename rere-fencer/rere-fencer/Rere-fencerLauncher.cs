using System;
using System.Collections.Generic;
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
            //new FastaWriter(GenerateReReFencedContigs(), _outputFile).WriteAllContigs();
            new ParallelFastaWriter(GenerateReReFencedContigs(), _outputFile).WriteAllContigs();
        }

        private IEnumerable<IGenomeContig> GenerateReReFencedContigs()
        {
            //var chromosomeOrder = _vcfReader.Chromosomes.ToList();
            //chromosomeOrder.AddRange(_genomeReader.Contigs.Select(v => (IContig) v.Value).Except(_vcfReader.Chromosomes));
            return _genomeReader.Contigs.Select(contig => (IContig) contig.Value).Select(contig => _vcfReader.Chromosomes.Contains(contig)
                ? new GenomeContig(contig, _rrfProcessor.Process(_genomeReader.Contigs[contig.Name], _vcfReader.GetVariantsForChromosome(contig, true)))
                : _genomeReader.Contigs[contig.Name]);
        }

    }
}
