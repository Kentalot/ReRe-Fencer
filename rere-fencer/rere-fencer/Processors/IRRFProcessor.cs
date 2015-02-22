using System;
using System.Collections.Generic;
using BioinformaticUtils.DataStructures.Nucleotides;
using BioinformaticUtils.GenomeTools;
using BioinformaticUtils.Vcf;

namespace rere_fencer.Processors
{
    public interface IRRFProcessor
    {
        IEnumerable<DnaNucleotide> Process(IGenomeReader genomeReader, IVcfReader vcfReader);
    }
}
